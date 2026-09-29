function pb_stim_diff_grid(flies, varargin)
% PB_STIM_DIFF_GRID  Grid of in-stim minus out-of-stim mean images.
%
%   pb_stim_diff_grid(flies, 'Name', value, ...)
%
% Rows = (fly, condition), columns = LED intensity; each tile is the mean
% background-subtracted frame during the LED pulses minus the mean outside
% them (res.glom(i).diff_img from pb_glom_delta_fly), averaged over the
% trials in that cell if there are several. Red = brighter during the pulse,
% blue = dimmer, white = no change. Works for one fly (the per-fly figure
% pb_glom_delta_fly writes) or many (the dataset-wide figure the scripts
% write).
%
% Color scale: always symmetric about zero, set from the climPct-th percentile
% of |diff|, computed either
%   climMode 'fly'   (default) once per FLY over all its trials, printed in the
%                    row label -- conditions within a fly are comparable, flies
%                    are not (different brightness), by design; or
%   climMode 'panel' per panel, printed inside each tile -- every panel uses
%                    its full dynamic range, so weak conditions become visible
%                    but amplitudes no longer compare across panels.
%
% flies   struct array with .label and .res (pb_glom_delta_fly output)
% options
%   intensityOrder ({'low','medium','high'})
%   condOrder      ({})   row order within a fly; conds not listed follow first-seen
%   climPct        (99)
%   climMode       ('fly') 'fly' | 'panel', see above
%   figNum         (17)
%   outFile        ('')   png to write; '' = just draw

p = inputParser;
p.addParameter('intensityOrder', {'low', 'medium', 'high'});
p.addParameter('condOrder', {});
p.addParameter('climPct', 99);
p.addParameter('climMode', 'fly', @(s) any(strcmpi(s, {'fly', 'panel'})));
p.addParameter('figNum', 17);
p.addParameter('outFile', '');
p.parse(varargin{:});
o = p.Results;
nInt = numel(o.intensityOrder);

% rows: (fly, cond) in fly order, conds ordered by condOrder then first-seen
rows = struct('fi', {}, 'cond', {});
for fi = 1:numel(flies)
    conds = unique({flies(fi).res.glom.cond}, 'stable');
    [~, ord] = ismember(o.condOrder, conds); ord = ord(ord > 0);
    conds = conds([ord, setdiff(1:numel(conds), ord, 'stable')]);
    for c = 1:numel(conds)
        rows(end+1) = struct('fi', fi, 'cond', conds{c}); %#ok<AGROW>
    end
end
nRows = numel(rows);

% per-fly color scale
cFly = zeros(1, numel(flies));
for fi = 1:numel(flies)
    d = cat(3, flies(fi).res.glom.diff_img);
    cFly(fi) = prctile(abs(double(d(:))), o.climPct);
end

half = 128;
rwb = [[linspace(0,1,half)', linspace(0,1,half)', ones(half,1)]; [ones(half,1), linspace(1,0,half)', linspace(1,0,half)']];

[imH, imW] = size(flies(1).res.glom(1).diff_img);
tileW = 300; tileH = round(tileW * imH / imW) + 28;
figure(o.figNum); clf
figPx = [260 + tileW*nInt, 60 + tileH*nRows];
set(gcf, 'Position', [50 50 figPx], 'Color', 'w')
setappdata(gcf, 'pb_target_px', figPx) % MATLAB clamps tall figures to the screen; pb_export_png prints at this size instead
tl = tiledlayout(nRows, nInt, 'TileSpacing', 'compact', 'Padding', 'loose');
for r = 1:nRows
    R = flies(rows(r).fi).res;
    for c = 1:nInt
        ax = nexttile(tl, (r-1)*nInt + c);
        idx = find(strcmp({R.glom.cond}, rows(r).cond) & strcmp({R.glom.intensity}, o.intensityOrder{c}));
        if isempty(idx)
            text(ax, 0.5, 0.5, 'missing', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'Color', [0.5 0.5 0.5])
        else
            d = mean(cat(3, R.glom(idx).diff_img), 3);
            if strcmpi(o.climMode, 'panel')
                cTile = prctile(abs(double(d(:))), o.climPct);
            else
                cTile = cFly(rows(r).fi);
            end
            imagesc(ax, d); clim(ax, [-1 1] * cTile); axis(ax, 'image')
            colormap(ax, rwb)
            if strcmpi(o.climMode, 'panel')
                text(ax, 0.02, 0.04, sprintf('%c%.0f', char(177), cTile), 'Units', 'normalized', ...
                    'HorizontalAlignment', 'left', 'FontSize', 7, 'FontWeight', 'bold', 'Color', [0.2 0.2 0.2])
            end
            if numel(idx) > 1
                text(ax, 0.98, 0.04, sprintf('mean of %d', numel(idx)), 'Units', 'normalized', ...
                    'HorizontalAlignment', 'right', 'FontSize', 7, 'Color', [0.3 0.3 0.3])
            end
        end
        % ticks off but axes kept ON: 'axis off' would also hide the ylabel used as
        % the row label below, and a free text object outside the axes gets clipped
        % because tiledlayout reserves no room for it
        set(ax, 'XTick', [], 'YTick', [], 'Box', 'off', 'XColor', 'none', 'YColor', 'none')
        if r == 1; title(ax, [o.intensityOrder{c} ' intensity'], 'FontSize', 10); end
        if c == 1
            if strcmpi(o.climMode, 'panel')
                rowLbl = {flies(rows(r).fi).label, rows(r).cond};
            else
                rowLbl = {flies(rows(r).fi).label, rows(r).cond, sprintf('%c%.0f', char(177), cFly(rows(r).fi))};
            end
            yl = ylabel(ax, rowLbl, 'Rotation', 0, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
                'Interpreter', 'none', 'FontSize', 8, 'FontWeight', 'bold');
            yl.Color = 'k'; % YColor 'none' would otherwise hide the label too
        end
    end
end
% Figure title as an annotation, not title(tl): the layout title's font scales
% with the figure when pb_export_png prints a clamped figure at its intended
% (much taller) size, and comes out several times too big; an annotation with
% a point-sized font does not.
if strcmpi(o.climMode, 'panel')
    scaleNote = 'color scale per PANEL (\pm value in each panel; amplitudes not comparable across panels)';
else
    scaleNote = 'color scale per fly (\pm value in row label)';
end
annotation(gcf, 'textbox', [0 0.985 1 0.015], 'String', ...
    {'in-stim minus out-of-stim mean image (background-subtracted F)', ...
     ['red = brighter during the pulse, blue = dimmer; ' scaleNote]}, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'EdgeColor', 'none', 'FontSize', 10, 'FitBoxToText', 'off')

if ~isempty(o.outFile)
    pb_export_png(gcf, o.outFile);
end
end
