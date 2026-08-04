function [movieOut, filt] = pb_remove_line_noise(movieIn, stopband, options)
%PB_REMOVE_LINE_NOISE Exactly zero a narrow spatial-frequency band along the fast-scan axis.
%
%   *** DO NOT USE THIS FOR DENOISING -- see PB_REMOVE_SCAN_NOISE INSTEAD. ***
%   Kept only so the reasoning trail (and the reason this approach was
%   abandoned) stays in the repo instead of silently vanishing.
%
%   On real data, this visibly damaged the actual PB signal: a slice through a
%   densely-packed glomerulus cluster has continuous real spatial-frequency
%   energy from ~5 cyc/line out past 35 cyc/line (checked directly on the
%   20260708-1 dataset), so there is no clean gap between "noise" (~13
%   cyc/line) and "signal" in this domain -- any spatial notch here removes
%   real glomerulus structure along with the noise. Concretely: running this on
%   a 300-volume chunk and comparing the time-averaged bridge region before/
%   after showed the sharp, distinct glomeruli blurring into wide, blobby,
%   horizontally-smeared bands (reported by the user as looking like "3 copies
%   of the image overlaid" -- an apt description of what removing a chunk of
%   the glomeruli's own real spatial spectrum does to their reconstruction).
%   This was true regardless of stopband width/taper or fit robustness (all
%   tried -- FIR+filtfilt, exact FFT projection, tapered mask, outlier-robust
%   fit -- see conversation/session history around 2026-07-10): the problem is
%   that the noise and the signal are not spectrally separable on the
%   fast-scan axis for this data, not a tuning mistake.
%
%   PB_REMOVE_SCAN_NOISE instead exploits the fact that the noise is
%   temporally incoherent while the PB signal is temporally stable, using a
%   short moving average across frames. That has no equivalent failure mode:
%   it cannot smear something sideways it never looks at.
%
%   ----------------------------------------------------------------------
%   Original docstring below, describing what this function does (not
%   whether you should use it):
%
%   movieOut = pb_remove_line_noise(movieIn) removes cycles 7-23 per line from
%   movieIn along dimension 2 (X, the fast-scan/resonant-galvo axis -- the
%   direction within a single scan line, matching how ScanImage and every other
%   loader in this repo lay out frames as [Y, X, ...]). PB_DIAGNOSE_LINE_NOISE
%   found the actual noise peak on the 20260708-1 dataset at ~13 cycles/line,
%   present at every z-plane and frame (a fixed-frequency electrical/timing
%   artifact, not fly signal).
%
%   movieOut = pb_remove_line_noise(movieIn, [f1 f2]) uses a custom stopband,
%   in cycles per line (i.e. cycles across the full width of one scan line).
%   f1, f2 are rounded to the nearest integer -- see METHOD below for why only
%   integer cycles/line are meaningful here.
%
%   movieIn may be [Y,X,T], [Y,X,Z,T], or any array with X as dimension 2; every
%   1D "line" (fixed Y, Z, T, ...) is treated independently and identically.
%
%   [movieOut, filt] = ... also returns the exact bins zeroed (filt.bins) and
%   the stopband used, so what was actually done is always inspectable.
%
%   METHOD (this replaces an earlier FIR/filtfilt bandstop -- see below for why):
%   Each line is exactly nX samples, so integer-cycle sine/cosine components
%   (1 cycle/line, 2 cycles/line, ...) are precisely the DFT basis for that
%   line. Zeroing the FFT bins for cycles f1:f2 (and their conjugate mirror
%   bins) removes exactly and only that frequency content, per line, then
%   ifft reconstructs the rest untouched. This is algebraically identical to a
%   least-squares fit-and-subtract of those sinusoids, just computed via FFT
%   instead of a matrix solve.
%
%   WHY NOT THE OLD FIR-FILTER APPROACH: an FIR/filtfilt bandstop filter (this
%   function's previous implementation, and remove_scannoise.m's cx_fft_filter_1d
%   before that) is a LOCAL convolution with a finite-length kernel. Narrow
%   enough to isolate ~13 cyc/line, that kernel has non-negligible spatial
%   extent (time-frequency uncertainty: a sharper notch needs a LONGER kernel,
%   not shorter), and since a bright point-like pixel (a hot pixel, or genuine
%   punctate signal) has energy at every frequency including the notch band,
%   convolving with that kernel leaves a visible ringing "streak" trailing off
%   in X from every bright point. This was confirmed on a synthetic
%   point-source test and on this dataset's real frames (see conversation/
%   session history around 2026-07-10) -- it is a real artifact of local
%   filtering, not motion correction and not a bug in filtfilt.
%   Per-line FFT-bin zeroing has no such issue in practice: because it acts on
%   the WHOLE line via its exact orthogonal Fourier basis rather than a local
%   sliding kernel, subtracting the same amount of frequency content leaves no
%   visible ringing around point sources on this dataset, is mathematically
%   exact (no window-design tradeoffs), and is ~2 orders of magnitude faster
%   (FFT/IFFT vs. a 40-tap filtfilt pass).
%
%   Name-value options:
%     verbose (1,1) logical = true
%
%   This only targets the periodic scan-line artifact. For the separate,
%   unstructured per-pixel shot noise (Poisson noise from low photon counts --
%   the sparse bright single-pixel puncta visible in raw frames), see
%   PB_REMOVE_SHOT_NOISE; the two are different noise sources and are handled
%   as two separate, independently-inspectable steps rather than one
%   do-everything filter.
%
%   See also PB_DIAGNOSE_LINE_NOISE, PB_REMOVE_SHOT_NOISE, PB_REGISTER.

arguments
    movieIn {mustBeNumeric}
    stopband (1,2) double = [7 23]
    options.verbose (1,1) logical = true
end

warning('pb_remove_line_noise:deprecated', ...
    ['This spatial notch filter was found to smear real glomerulus structure on real data ' ...
     '(the noise band overlaps real signal spatial frequencies -- see this function''s help ' ...
     'text for the evidence). Use PB_REMOVE_SCAN_NOISE instead unless you specifically need ' ...
     'to reproduce/inspect this abandoned approach.']);

sz = size(movieIn);
if numel(sz) < 2
    error('pb_remove_line_noise:badInput', 'movieIn must have at least 2 dimensions [Y,X,...].');
end
nX = sz(2);

f1 = round(stopband(1));
f2 = round(stopband(2));
if f1 <= 0 || f2 >= nX/2 || f1 > f2
    error('pb_remove_line_noise:badStopband', ...
        'stopband must satisfy 0 < f1 <= f2 < Nyquist (%.1f cyc/line for a %d-px line); got [%.2f %.2f].', ...
        nX/2, nX, stopband(1), stopband(2));
end
bins = f1:f2;

originalClass = class(movieIn);

if options.verbose
    fprintf('pb_remove_line_noise: zeroing %d cyc/line bins %d-%d (and mirrors) on %d lines of length %d\n', ...
        numel(bins), f1, f2, numel(movieIn)/nX, nX);
end

tic
Xf = fft(double(movieIn), [], 2);
Xf(:, bins+1, :, :, :, :) = 0;      % positive-frequency bins (1-indexed: bin k+1 <-> k cycles/line)
Xf(:, nX-bins+1, :, :, :, :) = 0;   % conjugate/mirror bins, keeps the result real-valued
movieOut = real(ifft(Xf, [], 2));
if options.verbose
    fprintf('pb_remove_line_noise: done in %.3f s\n', toc);
end

if any(strcmp(originalClass, {'int8','int16','int32','int64','uint8','uint16','uint32','uint64'}))
    warning('pb_remove_line_noise:castToSingle', ...
        'Input was %s; output is single because the filtered signal is not integer-valued (can go negative). Cast explicitly if you really want integers.', ...
        originalClass);
    movieOut = single(movieOut);
else
    movieOut = cast(movieOut, originalClass);
end

filt.stopband = [f1 f2];
filt.bins = bins;
filt.method = 'exact per-line FFT bin zeroing';

end
