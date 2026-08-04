function movieOut = pb_remove_shot_noise(movieIn, options)
%PB_REMOVE_SHOT_NOISE Spatial median filter for per-pixel shot noise.
%
%   movieOut = pb_remove_shot_noise(movieIn) applies a 2D median filter
%   (default 3x3) to every [Y,X] frame independently.
%
%   Scan noise (see PB_REMOVE_LINE_NOISE) is a structured, deterministic
%   sinusoidal artifact shared by every pixel on a line -- it can be surgically
%   removed in the frequency domain. Shot noise is the opposite: at the low
%   photon counts typical of single-frame GCaMP imaging, each pixel's count is
%   an approximately independent Poisson draw, which is exactly what the sparse
%   bright single-pixel "puncta" in a raw frame look like. There is no
%   frequency or template to notch out -- it is unstructured, so the
%   appropriate tool is a classical denoising filter, not a spectral one. This
%   repo's own denoise_batch_pablo.m / fft_scratchpad.m already reached for
%   medfilt2 for the same reason after their FFT-based step.
%
%   HONEST TRADEOFF: a median filter suppresses any pixel that disagrees with
%   its neighborhood, whether that pixel is noise or a genuinely sharp one- to
%   two-pixel-wide feature (the PB glomeruli in this dataset are themselves
%   only a few pixels across at this resolution/zoom). Expect peak brightness
%   of small real features to drop along with the noise. Run this on a test
%   chunk and compare against the un-median-filtered output before deciding to
%   use it -- it is not applied by default anywhere in this pipeline.
%
%   Name-value options:
%     kernelSize (1,2) double = [3 3]   median filter neighborhood (rows, cols)
%     verbose    (1,1) logical = true
%
%   See also PB_REMOVE_LINE_NOISE.

arguments
    movieIn {mustBeNumeric}
    options.kernelSize (1,2) double {mustBePositive, mustBeInteger} = [3 3]
    options.verbose (1,1) logical = true
end

sz = size(movieIn);
nFrames = prod(sz(3:end));
movieR = reshape(movieIn, sz(1), sz(2), nFrames);
originalClass = class(movieIn);

if options.verbose
    fprintf('pb_remove_shot_noise: median filter [%d %d] on %d frames\n', ...
        options.kernelSize(1), options.kernelSize(2), nFrames);
end

tic
out = zeros(size(movieR), 'like', double(movieR(1)));
for f = 1:nFrames
    out(:,:,f) = medfilt2(double(movieR(:,:,f)), options.kernelSize);
end
if options.verbose
    fprintf('pb_remove_shot_noise: done in %.2f s\n', toc);
end

movieOut = cast(reshape(out, sz), originalClass);

end
