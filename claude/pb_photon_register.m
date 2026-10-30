function [regMovie, shifts, info] = pb_photon_register(counts, options)
%PB_PHOTON_REGISTER Rigid sub-pixel motion correction for a sparse photon-count movie.
%
%   [regMovie, shifts] = pb_photon_register(counts) estimates one rigid
%   [dy dx] shift per volume for counts (Y x X x T, photons per pixel per
%   volume, e.g. from pb_read_raw_tif_planes(..., 'photonThreshold', thr)) and
%   returns the shifted movie (single) and shifts (T x 2, px, the amount the
%   content was displaced from the template; regMovie undoes it).
%
%   [regMovie, shifts] = pb_photon_register(counts, 'shifts', S) skips
%   estimation and just applies S -- used to give two independent thinned
%   halves of the same movie identical registration (pb_lowrank_denoise).
%
%   WHY NOT PB_REGISTER / NORMCORRE: at dLight photon rates (~0.03 photons per
%   pixel per frame, ~6500 photons per 10-plane volume) one volume carries
%   essentially no image, so per-volume template matching just locks onto
%   noise. Shifts are therefore estimated on bins of 'binVols' volumes,
%   Gaussian-smoothed, against the movie's own mean, with plain (not phase-
%   normalised -- that whitens the spectrum and amplifies the Poisson noise)
%   windowed cross-correlation and a parabolic sub-pixel peak fit. The
%   template is rebuilt from the registered movie and estimation repeated.
%
%   On 20260216-1 trial 1 (epg>dLight): 20-vol (2 s) bins showed only +/-1 px
%   slow x drift, but 5-vol (0.5 s) bins showed real 2-3 px dy excursions
%   lasting ~1 s at fast turn onsets, on top of ~0.7 px estimation noise. So
%   the default is 5-vol bins followed by a 3-bin moving median (kills single
%   bin outliers) and a 3-bin moving mean (halves the estimation noise), i.e.
%   ~1.5 s effective resolution. Holdout photon log-likelihood of the PB was
%   marginally better than with 2-s bins. Shifts are applied with bilinear
%   interpolation (zero fill), so regMovie is no longer integer counts -- do
%   any binomial thinning on 'counts' BEFORE calling this, then pass 'shifts'.
%
%   Name-value options:
%     shifts     (:,2) double = []   apply these instead of estimating
%     binVols    (1,1) double = 5    volumes per shift-estimation bin
%     sigma      (1,1) double = 2    px, Gaussian smoothing of each bin
%     maxShift   (1,2) double = [8 12]  [dy dx] search bound, px
%     nIter      (1,1) double = 2    template refinement passes
%     medWin     (1,1) double = 3    bins, moving median on the shift trace
%     meanWin    (1,1) double = 3    bins, moving mean after the median
%     verbose    (1,1) logical = true
%
%   info: binShifts (raw per-bin estimates), binCentres (volume index),
%   template (final, smoothed).
%
%   See also PB_READ_RAW_TIF_PLANES, PB_LOWRANK_DENOISE, PB_REGISTER.

arguments
    counts (:,:,:) {mustBeNumeric}
    options.shifts (:,2) double = zeros(0,2)
    options.binVols (1,1) double {mustBePositive, mustBeInteger} = 5
    options.sigma (1,1) double {mustBeNonnegative} = 2
    options.maxShift (1,2) double = [8 12]
    options.nIter (1,1) double {mustBePositive, mustBeInteger} = 2
    options.medWin (1,1) double = 3
    options.meanWin (1,1) double = 3
    options.verbose (1,1) logical = true
end

[H, W, T] = size(counts);
info = struct();

if isempty(options.shifts)
    nb = options.binVols; nB = floor(T/nb);
    bins = squeeze(sum(reshape(single(counts(:,:,1:nB*nb)), H, W, nb, nB), 3));
    if options.sigma > 0, bins = imgaussfilt(bins, options.sigma); end
    hn = @(n) 0.5 - 0.5*cos(2*pi*(0:n-1)'/(n-1));
    win = hn(H) * hn(W)';
    cy = floor(H/2) + 1; cx = floor(W/2) + 1;
    [ry, rx] = ndgrid((1:H) - cy, (1:W) - cx);
    allowed = abs(ry) <= options.maxShift(1) & abs(rx) <= options.maxShift(2);

    tmpl = mean(bins, 3);
    binShifts = zeros(nB, 2);
    for it = 1:options.nIter
        Ft = conj(fft2((tmpl - mean(tmpl(:))) .* win));
        for b = 1:nB
            I = bins(:,:,b);
            X = fftshift(real(ifft2(fft2((I - mean(I(:))) .* win) .* Ft)));
            X(~allowed) = -Inf;
            [~, ix] = max(X(:)); [r, c] = ind2sub([H W], ix);
            r = min(max(r, 2), H-1); c = min(max(c, 2), W-1);
            binShifts(b,:) = [r - cy + parab(X(r-1,c), X(r,c), X(r+1,c)), ...
                              c - cx + parab(X(r,c-1), X(r,c), X(r,c+1))];
        end
        sm = movmean(movmedian(binShifts, options.medWin, 1), options.meanWin, 1);
        centres = ((1:nB)' - 0.5) * nb + 0.5;
        shifts = interp1(centres, sm, (1:T)', 'linear', 'extrap');
        % rebuild the template from the bins shifted by their smoothed shift
        acc = zeros(H, W, 'single');
        for b = 1:nB, acc = acc + imtranslate(bins(:,:,b), -fliplr(sm(b,:)), 'linear', 'FillValues', 0); end
        tmpl = acc / nB;
        if options.verbose
            fprintf('pb_photon_register: pass %d, %d bins of %d vols: dy %.1f..%.1f, dx %.1f..%.1f px (smoothed)\n', ...
                it, nB, nb, min(sm(:,1)), max(sm(:,1)), min(sm(:,2)), max(sm(:,2)));
        end
    end
    info.binShifts = binShifts; info.binCentres = centres; info.template = tmpl;
else
    shifts = options.shifts;
    assert(size(shifts,1) == T, 'pb_photon_register: %d shifts for %d volumes', size(shifts,1), T);
end

regMovie = zeros(H, W, T, 'single');
for t = 1:T
    regMovie(:,:,t) = imtranslate(single(counts(:,:,t)), -fliplr(shifts(t,:)), 'linear', 'FillValues', 0);
end
end

function d = parab(a, b, c)
% sub-pixel offset of the vertex of the parabola through (-1,a), (0,b), (1,c)
den = a - 2*b + c;
if den >= 0, d = 0; else, d = 0.5 * (a - c) / den; end
end
