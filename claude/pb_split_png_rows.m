function files = pb_split_png_rows(png_file, cut_fracs, names, out_dir, options)
%PB_SPLIT_PNG_ROWS Cut an exported full-figure PNG into horizontal panel bands.
%
%   files = pb_split_png_rows(png_file, cut_fracs, names, out_dir) cuts the
%   image at (approximately) the given fractions of its height -- e.g. [0.5
%   0.75] for three bands -- and writes one PNG per band to out_dir, named
%   names{k}.png. Each cut is snapped to the centre of the run of blank pixel
%   rows (>= 0.05 in tall, within +/- SearchFrac) nearest the requested
%   fraction, so it falls in the whitespace between figure rows rather than
%   through a label -- pass fractions close to the intended gap. Each band is
%   then trimmed to its content plus a small margin.
%
%   Why work on the PNG instead of the live figure (cf. PB_EXPORT_PANELS):
%     a figure taller than the screen (Figure_3_claude.m: 1700 px on a
%     1080 px display) gets its canvas clamped to the screen shortly after it
%     is created, even when invisible, so any later render or position query
%     sees a squashed layout. The script's own exportgraphics call, made
%     right after building the figure, still has the full layout -- so the
%     panels are cut from that file.
%
%   Name-value options:
%     DropTitle  (1,1) logical = true   drop the first content strip above the
%                                       first long blank run near the top
%                                       (the sgtitle)
%     SearchFrac (1,1) double  = 0.06   half-width of the cut search window,
%                                       as a fraction of the image height
%     PadIn      (1,1) double  = 0.1    margin kept around each band (inches,
%                                       using the PNG's stored resolution;
%                                       300 dpi assumed if none)

arguments
    png_file (1,:) char
    cut_fracs (1,:) double
    names cell
    out_dir (1,:) char
    options.DropTitle (1,1) logical = true
    options.SearchFrac (1,1) double = 0.06
    options.PadIn (1,1) double = 0.1
end

img = imread(png_file);
if size(img,3) == 1; img = repmat(img,1,1,3); end
info = imfinfo(png_file);
dpi = 300;
if isfield(info,'XResolution') && ~isempty(info.XResolution) && info.XResolution > 0
    dpi = info.XResolution;
    if isfield(info,'ResolutionUnit') && strcmpi(info.ResolutionUnit,'meter'); dpi = dpi*0.0254; end
end
pad = round(options.PadIn*dpi);
H = size(img,1);
row_blank = all(all(img >= 250, 3), 2)';  % near-white rows (antialiasing-safe)

% band edges: optional title strip, then the snapped cuts
top = 1;
if options.DropTitle
    first = find(~row_blank, 1);
    % end of the title = start of the first blank run after it that is at
    % least ~0.05 in tall
    run_min = round(0.05*dpi);
    d = diff([0 row_blank 0]); s0 = find(d==1); s1 = find(d==-1) - 1;
    k = find(s0 > first & (s1 - s0 + 1) >= run_min, 1);
    if ~isempty(k) && s0(k) < 0.08*H; top = s1(k); end
end
cuts = zeros(1,numel(cut_fracs));
for i = 1:numel(cut_fracs)
    lo = max(1, round((cut_fracs(i)-options.SearchFrac)*H));
    hi = min(H, round((cut_fracs(i)+options.SearchFrac)*H));
    b = row_blank(lo:hi);
    d = diff([0 b 0]); s0 = find(d==1); s1 = find(d==-1) - 1;
    ok = (s1 - s0 + 1) >= round(0.05*dpi); s0 = s0(ok); s1 = s1(ok);
    if isempty(s0)
        error('pb_split_png_rows:noGap','no blank band near %.2f of the image height', cut_fracs(i));
    end
    % the band nearest the requested fraction, not the widest: a panel title
    % usually sits closer to its own axes than to the row above, so the
    % widest gap can be the one between a title and its axes
    ctr = lo - 1 + (s0 + s1)/2;
    [~,m] = min(abs(ctr - cut_fracs(i)*H));
    cuts(i) = round(ctr(m));
end
edges = [top cuts H];
assert(numel(names) == numel(edges)-1, 'need one name per band (%d)', numel(edges)-1);

if ~isfolder(out_dir); mkdir(out_dir); end
files = cell(numel(names),1);
for k = 1:numel(names)
    band = img(edges(k):edges(k+1), :, :);
    nb = ~all(band >= 250, 3);
    rr = find(any(nb,2)); cc = find(any(nb,1));
    band = band(max(rr(1)-pad,1):min(rr(end)+pad,size(band,1)), max(cc(1)-pad,1):min(cc(end)+pad,size(band,2)), :);
    files{k} = fullfile(out_dir, [names{k} '.png']);
    imwrite(band, files{k}, 'ResolutionUnit','meter', 'XResolution',dpi/0.0254, 'YResolution',dpi/0.0254)
end
end
