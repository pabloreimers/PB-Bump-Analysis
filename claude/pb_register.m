function [movieReg, shifts, diagInfo] = pb_register(movie, options)
%PB_REGISTER Rigid motion-correct a [Y,X,T] PB movie and report how stable it is.
%
%   movieReg = pb_register(movie) runs NoRMCorre rigid (translation-only)
%   registration via the lab's existing normcorre_regProduct (same numerics as
%   registration/normcorre_regProduct.m used by register_batch.m), after an
%   optional temporal pre-smoothing to help the shift estimator lock onto the
%   bridge when single-frame SNR is poor. This is deliberately the same
%   algorithm already used elsewhere in this repo -- the goal here is
%   observability (does the PB actually stop moving?), not a new registration
%   method.
%
%   Name-value options:
%     smoothWindow (1,1) double = 5     frames for a movmean temporal smooth
%                                        applied before estimating shifts (0 to
%                                        disable; matches register_batch.m's
%                                        default of 5)
%     maxShift     (1,2) double = [25 25]  max allowed [row,col] shift in pixels
%     verbose      (1,1) logical = true
%
%   [movieReg, shifts, diagInfo] = ... also returns the per-frame NoRMCorre shift
%   struct array and a diagInfo struct with:
%     shiftTrace      [T x 2] estimated (row, col) shift per frame
%     colProfileBefore/After  [X x T] column-sum kymograph before/after
%                              registration -- for a stable PB this should show
%                              flat, non-drifting vertical stripes (the bridge
%                              glomeruli columns), not diagonal drift
%     figHandles      handles to the two diagnostic figures created
%
%   movie is registered as-is (whatever smoothing/denoising you already applied,
%   e.g. PB_REMOVE_LINE_NOISE, should happen before calling this).
%
%   See also PB_REMOVE_LINE_NOISE.

arguments
    movie {mustBeNumeric}
    options.smoothWindow (1,1) double {mustBeNonnegative} = 5
    options.maxShift (1,2) double = [25 25]
    options.verbose (1,1) logical = true
    options.figureVisible (1,1) logical = true
    options.outlierJumpPx (1,1) double = 3   % after estimation, frames whose shift deviates
                               % from a 31-frame median-filtered shift trace by more than
                               % this are treated as estimator failures and replaced by
                               % the filtered value (0 disables). NoRMCorre occasionally
                               % locks onto a spurious peak for 1-2 frames (20261012-2:
                               % -15..-20 px jumps every ~45 volumes with true motion
                               % < 2 px), which paints a shifted ghost into the mean.
    options.applyShifts = []   % NoRMCorre shifts struct array (from a previous call on
                               % another movie of the same xy size and frame count, e.g.
                               % the PB planes of the same volumes): skip estimation and
                               % just apply these. Used for the EB planes, whose own
                               % registration is unstable because the structure is small.
end

if ndims(movie) ~= 3
    error('pb_register:badInput', 'movie must be [Y, X, T] (sum/select z-planes before registering); got %d dims.', ndims(movie));
end

% NoRMCorre lives in the lab's flyg repo (flyg\dependencies\NoRMCorre, plus
% flyg\registration for normcorre_regProduct). It used to be on the saved
% MATLAB path, but a MATLAB upgrade (R2026b, Oct 2026) reset pathdef, so
% find it relative to this repo (GitHub\LPsP_2p\...\claude -> GitHub\flyg) or
% the flyg-master copy under Documents, and add it if it is missing.
if isempty(which('NoRMCorreSetParms'))
    here = fileparts(mfilename('fullpath'));
    cands = {fullfile(here, '..', '..', '..', '..', 'flyg'), fullfile(here, '..', '..', '..', '..', '..', 'flyg-master')};
    for c = 1:numel(cands)
        if isfolder(fullfile(cands{c}, 'dependencies', 'NoRMCorre'))
            addpath(genpath(fullfile(cands{c}, 'dependencies', 'NoRMCorre')));
            if isfolder(fullfile(cands{c}, 'registration')), addpath(fullfile(cands{c}, 'registration')); end
            break
        end
    end
    if isempty(which('NoRMCorreSetParms'))
        error('pb_register:noNoRMCorre', 'NoRMCorre not found; add <flyg>\dependencies\NoRMCorre to the path.');
    end
end

movieSingle = single(movie);
[nY, nX, nT] = size(movieSingle); %#ok<ASGLU>

if options.smoothWindow > 0
    if options.verbose
        fprintf('pb_register: smoothing (movmean, window=%d frames) before shift estimation\n', options.smoothWindow);
    end
    regInput = smoothdata(movieSingle, 3, 'movmean', options.smoothWindow);
else
    regInput = movieSingle;
end

if options.verbose
    fprintf('pb_register: running NoRMCorre rigid registration on %d frames (max_shift=[%d %d])\n', ...
        nT, options.maxShift(1), options.maxShift(2));
end

options_rigid = NoRMCorreSetParms('d1', nY, 'd2', nX, ...
    'max_shift', [options.maxShift, 2], 'init_batch', min(100, nT), 'us_fac', 50);
tic
if ~isempty(options.applyShifts)
    shifts = options.applyShifts(:);
    assert(numel(shifts) == nT, 'pb_register:applyShifts', 'applyShifts has %d entries but the movie has %d frames.', numel(shifts), nT);
    if options.verbose, fprintf('pb_register: applying %d precomputed shifts (no estimation)\n', nT); end
    movieReg = apply_shifts(movieSingle, shifts, options_rigid);
else
    [~, shifts, ~, options_out] = normcorre(regInput, options_rigid);
    if options.outlierJumpPx > 0 && nT >= 31
        tr = cell2mat(arrayfun(@(s) reshape(s.shifts, 1, []), shifts(:), 'UniformOutput', false));   % nT x 2 (row, col)
        trF = movmedian(tr, 31, 1);   % 31 frames: the spurious locks last up to ~7 frames, real drift is much slower
        bad = any(abs(tr - trF) > options.outlierJumpPx, 2);
        if any(bad)
            if true
                fprintf('pb_register: %d/%d frames have a shift jumping > %g px from the median-filtered trace -> replaced by the filtered shift\n', nnz(bad), nT, options.outlierJumpPx);
            end
            for f = find(bad)'
                shifts(f).shifts(1, 1, 1, :) = trF(f, :);
                if isfield(shifts, 'shifts_up') && ~isempty(shifts(f).shifts_up), shifts(f).shifts_up(1, 1, 1, :) = trF(f, :); end
            end
        end
    end
    movieReg = apply_shifts(movieSingle, shifts, options_out);
end
if options.verbose
    fprintf('pb_register: done in %.1f s\n', toc);
end

shiftTrace = nan(nT, 2);
for f = 1:nT
    shiftTrace(f, :) = shifts(f).shifts(1, 1, 1, :);
end

colProfileBefore = squeeze(sum(movieSingle, 1)); % X x T
colProfileAfter  = squeeze(sum(movieReg, 1));    % X x T

diagInfo.shiftTrace = shiftTrace;
diagInfo.colProfileBefore = colProfileBefore;
diagInfo.colProfileAfter = colProfileAfter;

figVisArg = {'Visible', 'on'};
if ~options.figureVisible
    figVisArg = {'Visible', 'off'};
end

fig1 = figure(figVisArg{:}, 'Name', 'pb_register: estimated shifts');
plot(shiftTrace(:,1), 'DisplayName', 'row (Y) shift'); hold on
plot(shiftTrace(:,2), 'DisplayName', 'col (X) shift');
xlabel('frame'); ylabel('shift (px)'); legend; grid on
title('pb_register: per-frame rigid shift', 'Interpreter', 'none')

fig2 = figure(figVisArg{:}, 'Name', 'pb_register: PB stability kymograph', 'Position', [0 0 900 700]);
subplot(2,1,1)
imagesc(colProfileBefore'); colormap(gca, 'parula'); colorbar
xlabel('X (px, across the bridge)'); ylabel('frame'); title('before registration')
subplot(2,1,2)
imagesc(colProfileAfter'); colormap(gca, 'parula'); colorbar
xlabel('X (px, across the bridge)'); ylabel('frame'); title('after registration')
sgtitle('column-sum kymograph -- PB glomeruli columns should be vertical, not drifting')

diagInfo.figHandles = [fig1, fig2];

fprintf('pb_register: shift range row=[%.1f %.1f] px, col=[%.1f %.1f] px\n', ...
    min(shiftTrace(:,1)), max(shiftTrace(:,1)), min(shiftTrace(:,2)), max(shiftTrace(:,2)));

end
