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
end

if ndims(movie) ~= 3
    error('pb_register:badInput', 'movie must be [Y, X, T] (sum/select z-planes before registering); got %d dims.', ndims(movie));
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
[~, shifts, ~, options_out] = normcorre(regInput, options_rigid);
movieReg = apply_shifts(movieSingle, shifts, options_out);
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
