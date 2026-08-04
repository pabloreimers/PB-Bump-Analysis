function pb_play_movie(movie, options)
%PB_PLAY_MOVIE Play a [Y,X,T] movie in a figure window.
%
%   pb_play_movie(movie) plays the 3-D array movie (Y x X x T) frame by frame
%   in a single imagesc axis, contrast-scaled to the whole-movie intensity
%   range so the display does not flicker.
%
%   Name-value options:
%     fps    (1,1) double = 30     playback rate (frames per second)
%     clim   (1,2) double = []     [lo hi] display range; default is the
%                                   2nd/99.8th percentiles of the whole movie
%     cmap   (1,:) char   = 'gray' colormap name
%     loop   (1,1) logical = false replay from the start when the end is reached
%
%   Press any key while the figure is focused to stop playback.
%
%   Example:
%     movie = pb_read_tiff_movie('...\denoised_srdtrans.tif');
%     pb_play_movie(movie, 'fps', 20, 'loop', true);
%
%   See also PB_READ_TIFF_MOVIE.

arguments
    movie (:,:,:)
    options.fps  (1,1) double {mustBePositive} = 30
    options.clim (1,2) double = [0 0]
    options.cmap (1,:) char = 'gray'
    options.loop (1,1) logical = false
end

nFrames = size(movie, 3);

if isequal(options.clim, [0 0])
    lo = prctile(movie(:), 0.2);
    hi = prctile(movie(:), 99.8);
    if hi <= lo; hi = lo + 1; end
    options.clim = [lo hi];
end

fig = figure('Name', 'pb_play_movie', 'NumberTitle', 'off', 'Color', 'k');
% Stop playback on any keypress by flagging the figure.
setappdata(fig, 'stop', false);
set(fig, 'KeyPressFcn', @(~,~) setappdata(fig, 'stop', true));

ax  = axes('Parent', fig);
im  = imagesc(movie(:, :, 1), 'Parent', ax);
axis(ax, 'image', 'off');
colormap(ax, options.cmap);
caxis(ax, options.clim);
ttl = title(ax, sprintf('frame 1/%d', nFrames), 'Color', 'w');

dt = 1 / options.fps;
while isvalid(fig)
    for i = 1:nFrames
        if ~isvalid(fig) || getappdata(fig, 'stop'); return; end
        set(im, 'CData', movie(:, :, i));
        set(ttl, 'String', sprintf('frame %d/%d', i, nFrames));
        drawnow;
        pause(dt);
    end
    if ~options.loop; break; end
end

end
