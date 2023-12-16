function plot_gif_fast(inp, ncol, flipdim, filename, titopt)

szo = size(inp);
numdims = ndims(inp);

if ~exist('flipdim', 'var')
    flipdim = 0;
end

if numdims==2
    sznew = [size(inp) 1];
elseif numdims==3
    sznew = size(inp);
elseif numdims>3
    if flipdim
        inp = permute(inp, [1 2 4 3]);
        szo = size(inp);
    end
    inp = reshape(inp, size(inp,1), size(inp,2), []);
    sznew = size(inp);
end

h = figure;
for i = 1:sznew(end)

    if i==1
        himg = imshow(inp(:,:,i), 'InitialMagnification', 'fit');
    else
        himg.CData = inp(:,:,i);
    end
    axis off; axis image;

    sgtitle(titopt, 'FontSize', 6)

    frame = getframe(h);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im,ncol);

    if i == 1
        imwrite(imind,cm,filename, 'DelayTime', 0, 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'DelayTime', 0,'WriteMode','append');
    end
end

end