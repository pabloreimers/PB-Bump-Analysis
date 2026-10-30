function movie = pb_lowrank_movie(dn, volumes)
%PB_LOWRANK_MOVIE Reconstruct the denoised z-summed photon-rate movie from pb_lowrank_denoise.
%
%   movie = pb_lowrank_movie(dn) returns Y x X x T single, photons per pixel
%   per volume (all valid planes summed). Inside dn.mask it is the low-rank
%   rate; outside it is the static registered mean image (no dynamics are
%   modelled there -- it is background for display only). Negative rates
%   (possible where mu is tiny) are clipped to 0.
%
%   movie = pb_lowrank_movie(dn, volumes) reconstructs only those volumes,
%   e.g. 1:500, so you don't need the whole 600+ MB movie in memory.
%
%   See also PB_LOWRANK_DENOISE.

arguments
    dn (1,1) struct
    volumes (1,:) double = 1:size(dn.V, 1)
end

[H, W] = size(dn.mask);
rate = dn.mu + sqrt(dn.mu) .* (dn.U * diag(dn.S) * dn.V(volumes, :)');
movie = repmat(single(dn.meanImage), 1, 1, numel(volumes));
movie = reshape(movie, H*W, []);
movie(dn.mask(:), :) = single(max(rate, 0));
movie = reshape(movie, H, W, []);
end
