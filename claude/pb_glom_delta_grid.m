function pb_glom_delta_grid(flies, varargin)
% PB_GLOM_DELTA_GRID  Cross-fly grid of hemisphere-folded stim responses.
%
%   pb_glom_delta_grid(flies, 'Name', value, ...)
%
% Same content and styling as pb_glom_delta_fly's per-fly 1x3 overlay
% (columns = LED intensity, one line per trial, color = condition), stacked
% so every fly is a row. Each row keeps its own pairing (its own consensus
% shift P*) on the x-axis and the row label carries that P*. By default the
% y-axes are linked across ALL panels so amplitudes compare across flies;
% yLink 'row' links within each fly, and 'none' lets every panel autoscale to
% its own dynamic range (weak flies become readable, amplitudes no longer
% compare by eye). The legend is the union of conditions seen in any fly,
% ordered by condOrder then first-seen.
%
% flies   struct array, fields:
%   .label   row label (e.g. the fly name with the genotype suffix stripped)
%   .res     output of pb_glom_delta_fly
% options
%   intensityOrder ({'low','medium','high'})
%   condColors     (Map)  containers.Map cond -> [r g b]
%   unknownColor   ([.5 .5 .5])
%   condOrder      ({})   legend order
%   condLabel      ('condition')
%   titleStr       ('')   figure title; default built from condLabel
%   yLink          ('all') 'all' | 'row' | 'none' -- see above
%   fold           (true) true = the hemisphere-folded line (res.glom(i).folded, x = aligned
%                         pair); false = the raw per-glomerulus line (res.glom(i).delta,
%                         x = glomerulus 1..nClusters along the arch, dashed midline)
%   outFile        ('')   png to write; '' = just draw
% condColors / condOrder / condLabel default to what the flies were analysed
% with (flies(1).res.params), so a grid can be rebuilt from saved results.

p = inputParser;
p.addParameter('intensityOrder', {'low', 'medium', 'high'});
p.addParameter('condColors', containers.Map());
p.addParameter('unknownColor', [0.5 0.5 0.5]);
p.addParameter('condOrder', {});
p.addParameter('condLabel', 'condition');
p.addParameter('titleStr', '');
p.addParameter('yLink', 'all', @(s) any(strcmpi(s, {'all', 'row', 'none'})));
p.addParameter('fold', true);
p.addParameter('outFile', '');
p.addParameter('unitLabel', ''); % y-axis quantity; default taken from flies(1).res.unitLabel ('\Delta dF/F' or '\Delta z(F)')
p.parse(varargin{:});
o = p.Results;
if isempty(o.unitLabel)
    if isfield(flies(1).res, 'unitLabel'); o.unitLabel = flies(1).res.unitLabel; else; o.unitLabel = '\Delta dF/F'; end
end
if isfield(flies(1).res, 'params') % fall back to the settings the flies were analysed with
    prm = flies(1).res.params;
    if o.condColors.Count == 0 && isfield(prm, 'condColors'); o.condColors = prm.condColors; end
    if isempty(o.condOrder) && isfield(prm, 'condOrder');     o.condOrder  = prm.condOrder;  end
    if any(strcmp(p.UsingDefaults, 'condLabel')) && isfield(prm, 'condLabel'); o.condLabel = prm.condLabel; end
end
if isempty(o.titleStr)
    switch lower(o.yLink)
        case 'all';  yNote = 'shared y-axis';
        case 'row';  yNote = 'y-axis shared within each row';
        case 'none'; yNote = 'independent y-axis per panel';
    end
    if o.fold; what = 'hemisphere-folded'; else; what = 'per-glomerulus (unfolded)'; end
    o.titleStr = sprintf('%s %s per trial: rows = flies, cols = LED intensity, color = %s (%s)', what, o.unitLabel, o.condLabel, yNote);
end

nFlies = numel(flies);
nInt   = numel(o.intensityOrder);

figure(15); clf
figPx = [1500, 270*nFlies + 80];
set(gcf, 'Position', [50 50 figPx], 'Color', 'w')
setappdata(gcf, 'pb_target_px', figPx) % MATLAB clamps tall figures to the screen; pb_export_png prints at this size instead
tl = tiledlayout(nFlies, nInt, 'TileSpacing', 'compact', 'Padding', 'compact');
hLeg = gobjects(0); legNames = {};
for fi = 1:nFlies
    R = flies(fi).res;
    pairsF = R.align.pairs; nPairsF = size(pairsF, 1);
    nGlom = numel(R.glom(1).delta);
    rowAxes = gobjects(1, nInt);
    for c = 1:nInt
        ax = nexttile(tl, (fi-1)*nInt + c); hold(ax, 'on')
        rowAxes(c) = ax;
        yline(ax, 0, 'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off');
        if ~o.fold
            xline(ax, nGlom/2 + 0.5, '--', 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off'); % midline between hemispheres
        end
        for i = find(strcmp({R.glom.intensity}, o.intensityOrder{c}))
            cond = R.glom(i).cond;
            if isKey(o.condColors, cond); col = o.condColors(cond); else; col = o.unknownColor; end
            if o.fold
                h = plot(ax, 1:nPairsF, R.glom(i).folded, '-o', 'Color', col, 'LineWidth', 2, ...
                    'MarkerSize', 4, 'MarkerFaceColor', col, 'DisplayName', cond);
            else
                h = plot(ax, 1:nGlom, R.glom(i).delta, '-o', 'Color', col, 'LineWidth', 1.5, ...
                    'MarkerSize', 3, 'MarkerFaceColor', col, 'DisplayName', cond);
            end
            if ~ismember(cond, legNames)
                hLeg(end+1) = h; legNames{end+1} = cond; %#ok<AGROW>
            end
        end
        if o.fold
            xlim(ax, [0.5 nPairsF + 0.5]); xticks(ax, 1:nPairsF)
            xticklabels(ax, arrayfun(@(a,b) sprintf('%d/%d', a, b), pairsF(:,1), pairsF(:,2), 'UniformOutput', false))
            set(ax, 'FontSize', 7, 'XTickLabelRotation', 90)
        else
            xlim(ax, [0.5 nGlom + 0.5]); xticks(ax, 0:5:nGlom)
            set(ax, 'FontSize', 7)
        end
        if fi == 1; title(ax, [o.intensityOrder{c} ' intensity'], 'FontSize', 11); end
        if c == 1
            if o.fold
                if isfield(R.align, 'Psource') && strcmp(R.align.Psource, 'default')
                    pStr = sprintf('(P* = %d, default)', R.align.P); % shift not determinable from the data, see pb_glom_delta_fly
                else
                    pStr = sprintf('(P* = %d)', R.align.P);
                end
                ylabel(ax, {flies(fi).label, pStr, [o.unitLabel ' (folded)']}, 'FontSize', 9)
            else
                ylabel(ax, {flies(fi).label, o.unitLabel}, 'FontSize', 9)
            end
        end
        if fi == nFlies
            if o.fold
                xlabel(ax, 'aligned glomerulus pair (1st-half glomerulus / its 2nd-half partner)', 'FontSize', 8);
            else
                xlabel(ax, 'glomerulus (along the PB arch; dashed = midline)', 'FontSize', 8);
            end
        end
    end
    if strcmpi(o.yLink, 'row'); linkaxes(rowAxes, 'y'); end
end
if strcmpi(o.yLink, 'all'); linkaxes(findobj(gcf, 'Type', 'Axes'), 'y'); end

[~, ord] = ismember(o.condOrder, legNames); ord = ord(ord > 0);
ord = [ord, setdiff(1:numel(legNames), ord, 'stable')];
lgd = legend(hLeg(ord), legNames(ord), 'Interpreter', 'none', 'Box', 'off');
lgd.Layout.Tile = 'east';
% annotation rather than title(tl): the layout title's font scales with the
% figure when pb_export_png prints a clamped figure at its intended size
annotation(gcf, 'textbox', [0 0.985 1 0.015], 'String', o.titleStr, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', 'EdgeColor', 'none', 'FontSize', 11, 'FitBoxToText', 'off')

if ~isempty(o.outFile)
    pb_export_png(gcf, o.outFile);
end
end
