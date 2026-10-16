function fig = pb_spotcheck_gui(imgData, fly, opts)
% PB_SPOTCHECK_GUI  Scrub through a registered movie and the matching bump data.
%
%   fig = pb_spotcheck_gui(imgData, fly)
%   fig = pb_spotcheck_gui(imgData, fly, 'VideoFile', 'path\to\fictrac-raw-*.avi')
%
% Optional FicTrac video: with 'VideoFile' the raw FicTrac .avi is loaded,
% cropped to the imaging period, and shown in a third image panel next to the
% imaging frame. Cropping follows example_videos_script.m: the camera/sync
% light makes the summed intensity of the top 'SyncRows' rows of each video
% frame jump above its mean while imaging runs, so the first/last frame above
% the mean mark the imaging window. The cropped video is then assumed to span
% xb(1)..xb(end) evenly (same assumption as example_videos_script.m, which
% interpolates imaging onto linspace(xb(1),xb(end),nVideoFrames)), and each
% imaging frame k shows the video frame nearest in time to xb(k). Video is
% stored as uint8 grayscale (channel 1), so memory is ~H*W*nFrames bytes.
%
% Interactive replacement for the frame-by-frame animation loop in
% PB_spotcheck.m. Instead of playing through the movie, a slider (plus
% arrow keys / play button / frame box) picks one frame, and every panel
% shows the data for that frame:
%
%   top-left      registered image frame (scaled to the movie's 99th pct)
%   top-right     bump profile z(:,i) vs glomerulus angle, with the fitted
%                 bump position mu(i) as a red dot
%   bottom        full-trial heatmap of z (y = unwrapped alpha) with the cue
%                 trace in cyan and a vertical line at the current time
%
% Inputs
%   imgData  [H x W x T] registered movie (e.g. imgData_smooth_reg)
%   fly      one element of all_data, i.e. all_data(fly_num); needs
%            fly.ft.xb, fly.ft.xf, fly.ft.cue, fly.im.z, fly.im.alpha, fly.im.mu
%
% Controls
%   slider / left,right arrows   step 1 frame
%   shift + left,right           step 10 frames
%   home / end                   first / last frame
%   space or Play button         play / pause (real-time/3, like the script)
%   frame box                    type a frame number
%
% Why a function and not a script: the original script overwrote CData in a
% loop, so it could not be paused or rewound. Here the movie is only
% indexed per frame (no scaled copy is made), so it stays cheap on big movies.

arguments
    imgData
    fly
    opts.VideoFile (1,:) char = ''
    opts.SyncRows (1,1) double = 150
    opts.Parent = []        % draw into this uipanel/figure instead of a new figure
end

xb  = fly.ft.xb(:);
xf  = fly.ft.xf(:);
cue = fly.ft.cue(:);
z   = fly.im.z;
mu  = fly.im.mu(:);
nT  = min([size(imgData,3), size(z,2), numel(xb), numel(mu)]);

% Original script: tmp = 200*img/p99, shown with image() -> colormap indices
% clip at 64. Equivalent to CLim [0 64*p99/200]. Subsample frames for speed.
p100  = max(reshape(imgData(:,:,1:max(1,floor(nT/200)):nT), 1, []));
clim = [0, 64*p100/200];

hasVid = ~isempty(opts.VideoFile);
if hasVid
    vid = loadFictracVideo(opts.VideoFile, opts.SyncRows, numel(xf));
    nV  = size(vid,3);
    % video frame shown for each imaging frame (nearest in time)
    tv     = linspace(xb(1), xb(end), nV);
    vidIdx = round(interp1(tv, 1:nV, xb, 'nearest', 'extrap'));
    vidIdx = max(1, min(nV, vidIdx));
    pos1 = [0.03 0.58 0.28 0.38];
    posV = [0.34 0.58 0.28 0.38];
    pos3 = [0.69 0.58 0.28 0.36];
end
if ~hasVid
    pos1 = [0.05 0.58 0.38 0.38];
    pos3 = [0.56 0.58 0.40 0.36];
end

ownFig = isempty(opts.Parent);
if ownFig
    fig = figure('Name','PB spotcheck','NumberTitle','off','Color','w', ...
        'Units','normalized','Position',[0.1 0.1 0.8 0.8]);
    parent = fig;
else
    parent = opts.Parent;
    fig = ancestor(parent,'figure');
end
state.i = 1; state.playing = false;

% --- panels (leave bottom 10% for controls) ---
ax1 = axes('Parent',parent,'Position',pos1);
hImg = imagesc(ax1, imgData(:,:,1), clim);
axis(ax1,'image'); axis(ax1,'off'); colormap(ax1,'parula');
hTitle = title(ax1,'');

if hasVid
    axV = axes('Parent',parent,'Position',posV);
    hVid = imagesc(axV, vid(:,:,vidIdx(1)), [0 255]);
    axis(axV,'image'); axis(axV,'off'); colormap(axV,'gray');
    hVidTitle = title(axV,'');
end

ax3 = axes('Parent',parent,'Position',pos3);
hProf = plot(ax3, fly.im.alpha, z(:,1)); hold(ax3,'on');
hProf.XData(round(end/2)) = nan;      % break the line between hemispheres
hMu = plot(ax3, mu(1), 0, 'ro');
ylim(ax3,[-2,3]); xlim(ax3,[-pi,pi]);
xlabel(ax3,'\alpha (rad)'); ylabel(ax3,'z');

ax2 = axes('Parent',parent,'Position',[0.05 0.16 0.91 0.34]);
imagesc(ax2, xb, unwrap(fly.im.alpha), z); hold(ax2,'on');
hCue = plot(ax2, xf, -cue, 'c', 'LineWidth', 2);
hCue.YData(abs(diff(hCue.YData)) > pi) = nan;   % no wrap lines
hMuTrace = plot(ax2, xb, fly.im.mu, 'w', 'LineWidth', 2);
hMuTrace.YData(abs(diff(hMuTrace.YData)) > pi) = nan;   % no wrap lines
hLine = xline(ax2, xb(1), 'w', 'LineWidth', 1.5);
xlim(ax2,[xb(1) xb(end)]);
xlabel(ax2,'time (s)'); ylabel(ax2,'unwrapped \alpha');

% --- controls ---
sld = uicontrol(parent,'Style','slider','Units','normalized', ...
    'Position',[0.12 0.04 0.70 0.04],'Min',1,'Max',nT,'Value',1, ...
    'SliderStep',[1/max(nT-1,1), 100/max(nT-1,1)]);
addlistener(sld,'ContinuousValueChange',@(~,~) setFrame(round(sld.Value)));
uicontrol(parent,'Style','pushbutton','String','Play','Units','normalized', ...
    'Position',[0.03 0.04 0.07 0.04],'Callback',@(~,~) togglePlay());
btn = parent.Children(1);
edt = uicontrol(parent,'Style','edit','Units','normalized', ...
    'Position',[0.84 0.04 0.06 0.04],'String','1', ...
    'Callback',@(s,~) setFrame(round(str2double(s.String))));
uicontrol(parent,'Style','text','Units','normalized','Position',[0.90 0.04 0.07 0.04], ...
    'String',sprintf('/ %d',nT),'BackgroundColor','w','HorizontalAlignment','left');

fig.KeyPressFcn = @onKey;
if ownFig
    fig.CloseRequestFcn = @(s,~) onClose(s);
end
setFrame(1);

    function setFrame(k)
        if isnan(k), k = state.i; end
        k = max(1, min(nT, k));
        state.i = k;
        hImg.CData  = imgData(:,:,k);
        hProf.YData = z(:,k);
        hMu.XData   = mu(k);
        hLine.Value = xb(k);
        sld.Value   = k;
        edt.String  = num2str(k);
        hTitle.String = sprintf('frame %d / %d   t = %.2f s', k, nT, xb(k));
        if hasVid
            hVid.CData = vid(:,:,vidIdx(k));
            hVidTitle.String = sprintf('FicTrac video frame %d / %d', vidIdx(k), nV);
        end
        drawnow limitrate
    end

    function togglePlay()
        state.playing = ~state.playing;
        btn.String = ternary(state.playing, 'Pause', 'Play');
        dt = mean(diff(xb)) / 3;
        while state.playing && isvalid(parent)
            if state.i >= nT, state.i = 0; end
            setFrame(state.i + 1);
            pause(dt);
        end
    end

    function onKey(~, e)
        step = 1 + 9*any(strcmp(e.Modifier,'shift'));
        switch e.Key
            case 'rightarrow', setFrame(state.i + step);
            case 'leftarrow',  setFrame(state.i - step);
            case 'home',       setFrame(1);
            case 'end',        setFrame(nT);
            case 'space',      togglePlay();
        end
    end

    function onClose(s)
        state.playing = false;
        delete(s);
    end

    function out = ternary(c, a, b)
        if c, out = a; else, out = b; end
    end
end

function vid = loadFictracVideo(file, syncRows, nExpected)
% Read a FicTrac .avi as uint8 grayscale [H x W x T], cropped to imaging.
% Frames are read sequentially (readFrame) and only channel 1 is kept, which
% avoids the 4-D RGB array that read(vidObj) would allocate.
vr = VideoReader(file);
vid = zeros(vr.Height, vr.Width, max(1, ceil(vr.Duration*vr.FrameRate)) + 10, 'uint8');
s = zeros(size(vid,3), 1);
r = 1:min(syncRows, vr.Height);
n = 0;
while hasFrame(vr)
    f = readFrame(vr);
    n = n + 1;
    if n > size(vid,3)              % duration estimate was short; grow
        vid(:,:,end+1000) = 0; s(end+1000) = 0; %#ok<AGROW>
    end
    vid(:,:,n) = f(:,:,1);
    s(n) = sum(double(f(r,:,1)), 'all');
end
vid = vid(:,:,1:n); s = s(1:n);

above = s > mean(s);
i0 = find(above, 1, 'first');
i1 = find(above, 1, 'last');
vid = vid(:,:,i0:i1);
fprintf(['FicTrac video: %d frames read, cropped to %d..%d (%d frames); ' ...
    '%d fictrac samples in ft.xf\n'], n, i0, i1, size(vid,3), nExpected);
end
