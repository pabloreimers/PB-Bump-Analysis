function out = pb_read_raw_tif_planes(tifPath, framesPerVol, groups, varargin)
% PB_READ_RAW_TIF_PLANES  Read a ScanImage fastZ BigTIFF and return one
% z-summed movie per group of planes (e.g. PB planes and EB planes from the
% same volume), plus optionally the per-plane mean image.
%
%   out = pb_read_raw_tif_planes(tifPath, framesPerVol, groups)
%   out = pb_read_raw_tif_planes(..., 'maxVolumes', 500, 'planeMeans', true)
%
% tifPath      : the ScanImage .tif (little-endian BigTIFF, one strip/frame,
%                int16, one channel)
% framesPerVol : frames per volume INCLUDING flyback frames
%                (SI.hStackManager.numFramesPerVolumeWithFlyback)
% groups       : cell array; each cell is a vector of 1-based plane indices
%                within the volume to sum (e.g. {1:8, 11:18}); {} = none
% 'maxVolumes' : stop after this many complete volumes (default Inf)
% 'planeMeans' : also accumulate the mean image of every plane (default false)
% 'verbose'    : progress every 1000 volumes (default true)
% 'photonThreshold' : if set, binarise each frame (raw > thr = one photon)
%                before summing, so stacks are photon counts (default [])
% 'splitHalves' : with photonThreshold, also return out.halves{g} = {A, B},
%                two independent photon halves (A + B = stacks{g}). Each run
%                of above-threshold pixels along a scan line is ONE photon
%                pulse (~6% of pulses spill into the next pixel) and goes to
%                A or B whole; thinning the summed counts pixel by pixel
%                instead would split a pulse across halves and make them
%                share noise. Same z for both halves, unlike an odd/even
%                plane split. For pb_lowrank_denoise's cross-validation.
% 'seed'       : rng seed for splitHalves (default 1)
%
% out.stacks{g}   H x W x nVolumes single, sum over groups{g}
% out.planeMean   H x W x framesPerVol single (if 'planeMeans')
% out.nVolumes, out.nFramesRead, out.chainOk (IFD chain ended at a 0 pointer),
% out.H, out.W
%
% Why not MATLAB's Tiff class: libtiff indexes directories with a 16-bit
% counter, so Tiff/nextDirectory fails ("Unable to read the next directory")
% after exactly 65535 frames and a long recording is silently truncated. This
% walks the BigTIFF IFD chain itself with fread (header: 'II', 43, offset
% size 8; each IFD: uint64 entry count, 20-byte entries, uint64 next-IFD;
% tag 273 = strip offset, 256/257 = width/height, 258 = bits, 339 = sample
% format) and reads each frame's pixels by offset, summing volumes as it goes
% so the raw frames are never all in memory. Any trailing partial volume is
% dropped. Same logic as obs_read_raw_tif_summed in
% data scripts\overshoot_bar_stim_script.m (validated against Tiff.read).

p = inputParser;
p.addParameter('maxVolumes', Inf, @isnumeric);
p.addParameter('planeMeans', false, @islogical);
p.addParameter('verbose', true, @islogical);
p.addParameter('separate', [], @(v) islogical(v) || isnumeric(v));   % per group: true = keep planes separate (H x W x nPlanes x nVol) instead of summing
p.addParameter('crop', [], @(v) isempty(v) || (iscell(v) && numel(v) == numel(groups)));  % per group: [] or [row1 row2 col1 col2] applied to that group's output
p.addParameter('photonThreshold', [], @(v) isempty(v) || isscalar(v));  % if set, each frame becomes (raw > thr) BEFORE summing/planeMeans -- see pb_photon_movie
p.addParameter('splitHalves', false, @islogical);  % with photonThreshold: also return two independent random photon halves
p.addParameter('seed', 1, @isnumeric);
p.parse(varargin{:});
maxVol = p.Results.maxVolumes; doPlane = p.Results.planeMeans; verbose = p.Results.verbose;
thr = p.Results.photonThreshold; doSplit = p.Results.splitHalves;
assert(~doSplit || ~isempty(thr), 'pb_read_raw_tif_planes: splitHalves needs photonThreshold');
if doSplit, rng(p.Results.seed); end
nGroups = numel(groups);
sep = logical(p.Results.separate); if isempty(sep), sep = false(1, nGroups); end
crop = p.Results.crop; if isempty(crop), crop = repmat({[]}, 1, nGroups); end
for g = 1:nGroups
    assert(all(groups{g} >= 1 & groups{g} <= framesPerVol), 'pb_read_raw_tif_planes: plane indices in group %d out of 1..%d', g, framesPerVol);
end

d = dir(tifPath); assert(~isempty(d), 'pb_read_raw_tif_planes: %s not found', tifPath); fileSize = d.bytes;
fid = fopen(tifPath, 'r', 'ieee-le'); cleanupFid = onCleanup(@() fclose(fid)); %#ok<NASGU>
magic = fread(fid, 2, 'uint8=>char')'; version = fread(fid, 1, 'uint16');
assert(strcmp(magic, 'II') && version == 43, 'pb_read_raw_tif_planes: %s is not a little-endian BigTIFF', tifPath);
offsetSize = fread(fid, 1, 'uint16'); fread(fid, 1, 'uint16');
assert(offsetSize == 8, 'pb_read_raw_tif_planes: BigTIFF offset size %d, expected 8', offsetSize);
ptr = fread(fid, 1, 'uint64');

H = []; W = []; stacks = {}; planeMean = []; vol = [];
nVolAlloc = 0; k = 0; v = 0; f = 0; chainOk = true;
tic
while ptr > 0 && v < maxVol
    if ptr + 8 > fileSize, chainOk = false; if verbose, fprintf('  note: next-IFD pointer past end of file after %d frames.\n', k); end, break, end
    fseek(fid, ptr, 'bof');
    nEntries = fread(fid, 1, 'uint64');
    if isempty(nEntries) || nEntries < 1 || nEntries > 500, chainOk = false; if verbose, fprintf('  note: implausible IFD after %d frames.\n', k); end, break, end
    ent = fread(fid, [10 nEntries], 'uint16=>double');
    tags = ent(1,:); counts = ent(3,:) + ent(4,:)*2^16 + ent(5,:)*2^32 + ent(6,:)*2^48; vals = ent(7,:) + ent(8,:)*2^16 + ent(9,:)*2^32 + ent(10,:)*2^48;
    nextIFD = fread(fid, 1, 'uint64');
    dataOff = vals(tags == 273); nStrips = counts(tags == 273);
    assert(~isempty(dataOff) && nStrips == 1, 'pb_read_raw_tif_planes: frame %d has %d strips, expected 1', k+1, nStrips);
    if isempty(H)
        W = vals(tags == 256); H = vals(tags == 257); bps = vals(tags == 258); fmt = vals(tags == 339);
        assert(bps == 16 && (isempty(fmt) || fmt == 2), 'pb_read_raw_tif_planes: expected int16 frames (bits %d, format %s)', bps, mat2str(fmt));
        nVolAlloc = min(ceil(fileSize / (H * W * 2) / framesPerVol) + 10, maxVol);
        stacks = cell(1, nGroups); rows = cell(1, nGroups); cols = cell(1, nGroups);
        for g = 1:nGroups
            if isempty(crop{g}), rows{g} = 1:H; cols{g} = 1:W; else, rows{g} = crop{g}(1):crop{g}(2); cols{g} = crop{g}(3):crop{g}(4); end
            if sep(g), stacks{g} = zeros(numel(rows{g}), numel(cols{g}), numel(groups{g}), nVolAlloc, 'single');
            else, stacks{g} = zeros(numel(rows{g}), numel(cols{g}), nVolAlloc, 'single'); end
        end
        if doSplit, stacksA = stacks; volA = zeros(H, W, framesPerVol, 'single'); end
        if doPlane, planeMean = zeros(H, W, framesPerVol, 'double'); end
        vol = zeros(H, W, framesPerVol, 'single');
    end
    fseek(fid, dataOff, 'bof');
    frame = fread(fid, [W H], 'int16=>single');
    if numel(frame) ~= W*H, chainOk = false; if verbose, fprintf('  note: frame %d truncated.\n', k+1); end, break, end
    k = k + 1; f = f + 1;
    if ~isempty(thr)
        ev = frame > thr;                     % W x H: dim 1 runs along each scan line
        if doSplit
            % one photon pulse can span adjacent pixels along the line, so
            % assign whole runs (not pixels) to half A or B
            starts = ev & ~[false(1, H); ev(1:end-1, :)];
            id = cumsum(starts(:)); coin = rand(id(end), 1) < 0.5;
            inA = false(size(ev)); inA(ev) = coin(id(ev));
            volA(:,:,f) = single(inA)';
        end
        frame = single(ev);
    end
    vol(:,:,f) = frame';
    if f == framesPerVol
        v = v + 1;
        for g = 1:nGroups
            if sep(g), stacks{g}(:,:,:,v) = vol(rows{g}, cols{g}, groups{g});
            else, stacks{g}(:,:,v) = sum(vol(rows{g}, cols{g}, groups{g}), 3); end
            if doSplit
                if sep(g), stacksA{g}(:,:,:,v) = volA(rows{g}, cols{g}, groups{g});
                else, stacksA{g}(:,:,v) = sum(volA(rows{g}, cols{g}, groups{g}), 3); end
            end
        end
        if doPlane, planeMean = planeMean + double(vol); end
        f = 0;
        if verbose && mod(v, 1000) == 0, fprintf('  read %d volumes so far (%.0f frames/s)\n', v, k/toc); end
    end
    ptr = nextIFD;
end
assert(v > 0, 'pb_read_raw_tif_planes: no complete volumes read from %s', tifPath);
for g = 1:nGroups
    if sep(g), stacks{g} = stacks{g}(:,:,:,1:v); else, stacks{g} = stacks{g}(:,:,1:v); end
    if doSplit
        if sep(g), stacksA{g} = stacksA{g}(:,:,:,1:v); else, stacksA{g} = stacksA{g}(:,:,1:v); end
        out.halves{g} = {stacksA{g}, stacks{g} - stacksA{g}};
    end
end
out.stacks = stacks; out.rows = rows; out.cols = cols; out.separate = sep;
out.nVolumes = v; out.nFramesRead = k; out.chainOk = chainOk && (ptr == 0 || v >= maxVol); out.H = H; out.W = W;
if doPlane, out.planeMean = single(planeMean / v); end
if verbose
    if ptr == 0, chainStr = 'ended cleanly'; elseif v >= maxVol, chainStr = 'stopped at maxVolumes'; else, chainStr = 'ended ABNORMALLY'; end
    fprintf('  %d frames read -> %d complete volumes (%d leftover frame(s) dropped); IFD chain %s\n', k, v, f, chainStr);
end
end
