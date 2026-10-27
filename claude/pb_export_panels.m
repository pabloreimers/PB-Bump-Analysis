function [files, rects] = pb_export_panels(fig, out_dir, prefix, options)
%PB_EXPORT_PANELS Save each panel of a multi-panel figure as its own PNG.
%
%   files = pb_export_panels(fig, out_dir, prefix) finds every axes in fig,
%   groups axes that physically abut into one panel, and writes one PNG per
%   panel to out_dir, named <prefix>_r<row>c<col>.png (rows top to bottom,
%   columns left to right, both counted in panel units, not tile units).
%
%   Why render-once-and-crop:
%     copying axes into a fresh figure loses their layout-computed size and
%     legend/colorbar placement, and hiding everything outside a panel makes
%     tiledlayouts reflow (hidden axes free up their label space, so nested
%     layouts like Figure 3's get squashed). Instead the whole figure is
%     rendered ONCE (print -RGBImage at the figure's own size, i.e. exactly
%     what the full-figure export shows) and each panel is cropped out of
%     that image. Crop lines run through the widest blank band between a
%     panel and each neighbouring panel (one overlapping it in x or y), then
%     each crop is trimmed to its content plus a small margin. The top-level
%     title (sgtitle / root layout title) is recoloured to the background for
%     the render (not removed, so the layout doesn't move).
%
%   Grouping: two axes are one panel when their inner boxes are within
%   gapPx of each other along one dimension and overlap by >= 50% along the
%   other -- e.g. a trace axes stacked on a heatmap with TileSpacing 'none'.
%   Ordinary tile spacing ('compact' and up) keeps panels separate.
%
%   Name-value options:
%     Resolution (1,1) double = 300   render dpi
%     GapPx      (1,1) double = 4     max gap (pixels) for abutting axes
%     Groups     cell = {}            skip the automatic grouping: each cell is
%                                     an array of axes exported as one panel
%                                     (e.g. several figure rows together)
%     Names      cell = {}            file names (no extension) for Groups;
%                                     default <prefix>_p<k>
%
%   Returns the full paths of the written files, in panel order, and the
%   crop rectangles [x1 x2 y1 y2] in the full render (rows of rects, same order).

arguments
    fig (1,1) matlab.ui.Figure
    out_dir (1,:) char
    prefix (1,:) char
    options.Resolution (1,1) double = 300
    options.GapPx (1,1) double = 4
    options.Groups cell = {}
    options.Names cell = {}
end

if ~isfolder(out_dir); mkdir(out_dir); end
drawnow

% axes that draw anything (axis-off image axes still count)
axs = findall(fig,'Type','axes');
keep = false(size(axs));
for i = 1:numel(axs)
    keep(i) = strcmp(axs(i).Visible,'on') || ...
        any(strcmp(get(findall(axs(i).Children,'flat'),'Visible'),'on'));
end
axs = axs(keep);
n = numel(axs);
pos = zeros(n,4);
for i = 1:n; pos(i,:) = getpixelposition(axs(i),true); end

if isempty(options.Groups)
    % union-find over abutting axes
    lbl = 1:n;
    for i = 1:n
        for j = i+1:n
            if abuts(pos(i,:),pos(j,:),options.GapPx)
                lbl(lbl==lbl(j)) = lbl(i);
            end
        end
    end
    [~,~,lbl] = unique(lbl);
    n_p = max(lbl);
    box = zeros(n_p,4); % [left bottom right top]
    for k = 1:n_p
        p = pos(lbl==k,:);
        box(k,:) = [min(p(:,1)), min(p(:,2)), max(p(:,1)+p(:,3)), max(p(:,2)+p(:,4))];
    end

    % rows: walk panels top-down; a panel joins the current row when it overlaps
    % the row's first panel vertically by >= 30% of the shorter height
    [~,ord] = sort(box(:,4),'descend');
    row = zeros(n_p,1); r = 0; ref = [];
    for k = ord'
        if isempty(ref) || vfrac(box(k,:),ref) < 0.3
            r = r + 1; ref = box(k,:);
        end
        row(k) = r;
    end
    col = zeros(n_p,1);
    for r = 1:max(row)
        in_r = find(row==r);
        [~,o] = sort(box(in_r,1));
        col(in_r(o)) = 1:numel(in_r);
    end
    names = arrayfun(@(k) sprintf('%s_r%dc%d',prefix,row(k),col(k)), (1:n_p)', 'UniformOutput', false);
else
    n_p = numel(options.Groups);
    lbl = zeros(n,1);
    for k = 1:n_p; lbl(ismember(axs, options.Groups{k})) = k; end
    axs = axs(lbl>0); pos = pos(lbl>0,:); lbl = lbl(lbl>0);
    row = (1:n_p)'; col = ones(n_p,1);
    names = options.Names(:);
    if isempty(names); names = arrayfun(@(k) sprintf('%s_p%d',prefix,k), (1:n_p)', 'UniformOutput', false); end
end

% render once, with the figure title(s) recoloured to the background
bg = fig.Color;
ttl = findall(fig,'Type','subplottext');
for L = findall(fig,'Type','tiledlayout')'
    if isa(L.Parent,'matlab.ui.Figure'); ttl = [ttl; L.Title; L.Subtitle]; end %#ok<AGROW>
end
ttl_col = get(ttl,'Color'); if ~iscell(ttl_col); ttl_col = {ttl_col}; end
ppm = fig.PaperPositionMode; ihc = fig.InvertHardcopy;
fig.PaperPositionMode = 'auto'; fig.InvertHardcopy = 'off';
if ~isempty(ttl); set(ttl,'Color',bg); end
try
    % exportgraphics with Padding 'figure' (R2025a+) renders at the size in
    % fig.Position even when the on-screen canvas was clamped to the screen;
    % print (the fallback for older releases) renders the canvas as it is
    tmp_png = [tempname '.png'];
    try
        exportgraphics(fig, tmp_png, 'Resolution', options.Resolution, 'Padding', 'figure');
        img = imread(tmp_png); delete(tmp_png);
    catch
        img = print(fig,'-RGBImage',sprintf('-r%d',options.Resolution));
    end
catch err
    restore(); rethrow(err)
end
restore();

% axes boxes -> image coordinates [x1 x2 y1 y2] (y from the top)
fpos = fig.Position;
if ~strcmp(fig.Units,'pixels'); fpos = getpixelposition(fig); end
sx = size(img,2)/fpos(3); sy = size(img,1)/fpos(4);
fprintf('pb_export_panels: figure %dx%d px, render %dx%d\n', round(fpos(3)), round(fpos(4)), size(img,2), size(img,1));
if abs(sx/sy - 1) > 0.02
    % seen with a figure taller than the screen once another window had been
    % shown in the session: Position keeps the requested size but the canvas
    % (and so the render) is clamped to the screen
    error('pb_export_panels:clampedCanvas', ['render is %dx%d px but the figure is %dx%d -- ' ...
        'its canvas was probably clamped to the screen size; export before any figure window is shown'], ...
        size(img,2), size(img,1), round(fpos(3)), round(fpos(4)));
end
ibox = zeros(n_p,4);
for k = 1:n_p
    p = pos(lbl==k,:);
    l = min(p(:,1)); b = min(p(:,2)); r = max(p(:,1)+p(:,3)); t = max(p(:,2)+p(:,4));
    ibox(k,:) = round([l*sx, r*sx, (fpos(4)-t)*sy, (fpos(4)-b)*sy]);
end
ibox = max(ibox,1); ibox(:,[1 2]) = min(ibox(:,[1 2]),size(img,2)); ibox(:,[3 4]) = min(ibox(:,[3 4]),size(img,1));
blank = all(img == reshape(uint8(round(255*bg)),1,1,3), 3);
pad = round(0.1*options.Resolution); % 0.1 in margin around each trimmed panel

files = cell(n_p,1); rects = zeros(n_p,4);
[~,porder] = sortrows([row col]);
for k = porder'
    X1 = 1; X2 = size(img,2); Y1 = 1; Y2 = size(img,1);
    a = ibox(k,:);
    for j = setdiff(1:n_p,k)
        o = ibox(j,:);
        xo = [max(a(1),o(1)) min(a(2),o(2))]; yo = [max(a(3),o(3)) min(a(4),o(4))];
        if xo(2) > xo(1) && yo(2) <= yo(1)       % vertical neighbour
            if o(4) <= a(3); Y1 = max(Y1, cut_line(blank, o(4), a(3), xo, 2));
            else;            Y2 = min(Y2, cut_line(blank, a(4), o(3), xo, 2)); end
        elseif yo(2) > yo(1) && xo(2) <= xo(1)   % horizontal neighbour
            if o(2) <= a(1); X1 = max(X1, cut_line(blank, o(2), a(1), yo, 1));
            else;            X2 = min(X2, cut_line(blank, a(2), o(1), yo, 1)); end
        end
    end
    rects(k,:) = [X1 X2 Y1 Y2];
    crop = img(Y1:Y2, X1:X2, :);
    nb = ~blank(Y1:Y2, X1:X2);
    rr = find(any(nb,2)); cc = find(any(nb,1));
    if ~isempty(rr)
        crop = crop(max(rr(1)-pad,1):min(rr(end)+pad,end), max(cc(1)-pad,1):min(cc(end)+pad,end), :);
    end
    files{k} = fullfile(out_dir, [names{k} '.png']);
    imwrite(crop, files{k})
end
files = files(porder); rects = rects(porder,:);

    function restore()
        fig.PaperPositionMode = ppm; fig.InvertHardcopy = ihc;
        if ~isempty(ttl); set(ttl, {'Color'}, ttl_col(:)); end
    end
end

function c = cut_line(blank, lo, hi, span, dim)
% centre of the longest run of fully blank lines between lo and hi (rows if
% dim==2, columns if dim==1), blank across span; midpoint if there is none
lo = max(lo,1); hi = max(hi,lo);
if dim == 2
    is_b = all(blank(lo:hi, span(1):span(2)), 2)';
else
    is_b = all(blank(span(1):span(2), lo:hi), 1);
end
d = diff([0 is_b 0]); s0 = find(d==1); s1 = find(d==-1) - 1;
if isempty(s0); c = round((lo+hi)/2); return; end
[~,m] = max(s1 - s0);
c = lo - 1 + round((s0(m)+s1(m))/2);
end

function tf = abuts(a, b, gap)
dx = max(0, max(a(1),b(1)) - min(a(1)+a(3),b(1)+b(3)));
dy = max(0, max(a(2),b(2)) - min(a(2)+a(4),b(2)+b(4)));
ox = min(a(1)+a(3),b(1)+b(3)) - max(a(1),b(1));
oy = min(a(2)+a(4),b(2)+b(4)) - max(a(2),b(2));
tf = (dy <= gap && ox >= 0.5*min(a(3),b(3))) || ...
     (dx <= gap && oy >= 0.5*min(a(4),b(4)));
end

function f = vfrac(a, b) % a, b as [left bottom right top]
f = (min(a(4),b(4)) - max(a(2),b(2))) / min(a(4)-a(2), b(4)-b(2));
end
