function movieOut = pb_remove_scan_noise(movieIn, window, options)
%PB_REMOVE_SCAN_NOISE Suppress scan-line noise with a short temporal mean, not a spatial filter.
%
%   movieOut = pb_remove_scan_noise(movieIn) applies a short moving-average
%   filter along the TIME dimension (the last dimension of movieIn) with a
%   10-frame window. Nothing is done spatially -- see WHY THIS EXISTS below.
%
%   THIS REPLACES PB_REMOVE_LINE_NOISE AS THE RECOMMENDED DEFAULT. That function
%   notches spatial frequencies 7-23 cyc/line out of the fast-scan axis. On
%   real data that turned out to be wrong: a slice through a densely-packed
%   glomerulus cluster has continuous real spatial-frequency energy from ~5
%   cyc/line out past 35 cyc/line -- there is no clean gap between "noise" and
%   "signal" in that domain, so any spatial notch removes real structure along
%   with the noise. On this dataset it visibly smeared/blurred the glomeruli in
%   the bridge into blobby streaks (confirmed directly: see pb_remove_line_noise.m
%   for the full account and conversation/session history around 2026-07-10).
%
%   WHY THIS WORKS INSTEAD: the noise's spatial phase is essentially
%   uncorrelated frame-to-frame (checked directly: cross-correlating the
%   isolated ~13 cyc/line component between different frames gave ~0.2, close
%   to uncorrelated), while real image content (the PB and its glomeruli) is
%   stable over many consecutive frames at this ~10.5 Hz frame rate. Averaging
%   a short window of frames cancels differently-phased copies of the noise
%   while leaving spatial structure -- at ANY spatial frequency -- completely
%   untouched, because this operation never looks across pixels, only across
%   time. There is no equivalent of the notch-filter tradeoff here: a temporal
%   mean cannot smear something sideways it never looks at.
%
%   Checked with PB_DIAGNOSE_LINE_NOISE on the 20260708-1 dataset: a plain MEAN
%   over ~7-10 frames is enough that the 13 cyc/line noise peak is no longer
%   even the dominant non-DC spectral feature (real low-frequency PB structure
%   takes over instead). A temporal MEDIAN of the same window is NOT as
%   effective (still dominated by the 13 cyc/line peak at window=10) -- median
%   is the right tool for rejecting rare outliers (see PB_REMOVE_SHOT_NOISE),
%   not for cancelling a phase-random additive oscillation, which is exactly
%   what a mean is good at.
%
%   movieOut = pb_remove_scan_noise(movieIn, window) sets the averaging window
%   in frames (default 10). This trades noise suppression for temporal
%   resolution -- it is NOT the same 30-frame smoothing already used downstream
%   in show_single_PB_trial.m; the point of this function is to get away with a
%   much shorter window than that by targeting the noise specifically. Re-check
%   with PB_DIAGNOSE_LINE_NOISE on your own data rather than trusting the
%   default blindly -- how short a window works depends on how phase-randomized
%   the noise really is on a given day/rig.
%
%   Name-value options:
%     verbose (1,1) logical = true
%
%   See also PB_DIAGNOSE_LINE_NOISE, PB_REMOVE_SHOT_NOISE, PB_REMOVE_LINE_NOISE.

arguments
    movieIn {mustBeNumeric}
    window (1,1) double {mustBeInteger, mustBePositive} = 10
    options.verbose (1,1) logical = true
end

sz = size(movieIn);
if numel(sz) < 3
    error('pb_remove_scan_noise:badInput', ...
        'movieIn must have a time dimension (e.g. [Y,X,T]); got %d dims.', numel(sz));
end
timeDim = numel(sz);
originalClass = class(movieIn);

if options.verbose
    fprintf('pb_remove_scan_noise: %d-frame temporal moving mean (dim %d) on a %s array\n', ...
        window, timeDim, mat2str(sz));
end

tic
movieOut = movmean(double(movieIn), window, timeDim);
if options.verbose
    fprintf('pb_remove_scan_noise: done in %.2f s\n', toc);
end

if any(strcmp(originalClass, {'int8','int16','int32','int64','uint8','uint16','uint32','uint64'}))
    warning('pb_remove_scan_noise:castToSingle', ...
        'Input was %s; output is single because a moving average is not integer-valued. Cast explicitly if you really want integers.', ...
        originalClass);
    movieOut = single(movieOut);
else
    movieOut = cast(movieOut, originalClass);
end

end
