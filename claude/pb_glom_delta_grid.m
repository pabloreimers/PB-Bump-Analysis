function pb_glom_delta_grid(flies, varargin)
% PB_GLOM_DELTA_GRID  Cross-fly grid of hemisphere-folded stim responses.
%
%   pb_glom_delta_grid(flies, 'Name', value, ...)
%
% Same content and styling as pb_glom_delta_fly's per-fly 1x3 overlay
% (columns = LED intensity, one line per trial, color = condition), stacked
% so every fly is a row. Each row keeps its own pairing (its own consensus
% shift P*) on the x-axis and the row label carries that P*. y-axes are
% linked across ALL panels so amplitudes compare across flies (set sharedY
% false to link per row instead). The legend is the union of conditions seen
% in any fly, ordered by condOrder then first-seen.
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
%   sharedY        (true)
%   outFile        ('')   png to write; '' = just draw

p = inputParser;
p.addParameter('intensityOrder', {'low', 'medium', 'high'});
p.addParameter('condColors', containers.Map());
p.addParameter('unknownColor', [0.5 0.5 0.5]);
p.addParameter('condOrder', {});
p.addParameter('condLabel', 'condition');
p.addParameter('titleStr', '');
p.addParameter('sharedY', true);
p.addParameter('outFile', '');
p.parse(varargin{:});
o = p.Results;
if isempty(o.titleStr)
    o.titleStr = sprintf('hemisphere-folded \\DeltadF/F per trial: rows = flies, cols = LED intensity, color = %s', o.condLabel);
end

nFlies = numel(flies);
nInt   = numel(o.intensityOrder);

figure(15); clf
set(gcf, 'Position', [50 50 1500 270*nFlies + 80], 'Color', 'w')
tl = tiledlayout(nFlies, nInt, 'TileSpacing', 'compact', 'Padding', 'compact');
hLeg = gobjects(0); legNames = {};
for fi = 1:nFlies
    R = flies(fi).res;
    pairsF = R.align.pairs; nPairsF = size(pairsF, 1);
    rowAxes = gobjects(1, nInt);
    for c = 1:nInt
        ax = nexttile(tl, (fi-1)*nInt + c); hold(ax, 'on')
        rowAxes(c) = ax;
        yline(ax, 0, 'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off');
        for i = find(strcmp({R.glom.intensity}, o.intensityOrder{c}))
            cond = R.glom(i).cond;
            if isKey(o.condColors, cond); col = o.condColors(cond); else; col = o.unknownColor; end
            h = plot(ax, 1:nPairsF, R.glom(i).folded, '-o', 'Color', col, 'LineWidth', 2, ...
                'MarkerSize', 4, 'MarkerFaceColor', col, 'DisplayName', cond);
            if ~ismember(cond, legNames)
                hLeg(end+1) = h; legNames{end+1} = cond; %#ok<AGROW>
            end
        end
        xlim(ax, [0.5 nPairsF + 0.5]); xticks(ax, 1:nPairsF)
        xticklabels(ax, arrayfun(@(a,b) sprintf('%d/%d', a, b), pairsF(:,1), pairsF(:,2), 'UniformOutput', false))
        set(ax, 'FontSize', 7, 'XTickLabelRotation', 90)
        if fi == 1; title(ax, [o.intensityOrder{c} ' intensity'], 'FontSize', 11); end
        if c == 1
            if isfield(R.align, 'Psource') && strcmp(R.align.Psource, 'default')
                pStr = sprintf('(P* = %d, default)', R.align.P); % shift not determinable from the data, see pb_glom_delta_fly
            else
                pStr = sprintf('(P* = %d)', R.align.P);
            end
            ylabel(ax, {flies(fi).label, pStr, '\Delta dF/F (folded)'}, 'FontSize', 9)
        end
        if fi == nFlies; xlabel(ax, 'aligned glomerulus pair (1st-half glomerulus / its 2nd-half partner)', 'FontSize', 8); end
    end
    if ~o.sharedY; linkaxes(rowAxes, 'y'); end
end
if o.sharedY; linkaxes(findobj(gcf, 'Type', 'Axes'), 'y'); end

[~, ord] = ismember(o.condOrder, legNames); ord = ord(ord > 0);
ord = [ord, setdiff(1:numel(legNames), ord, 'stable')];
lgd = legend(hLeg(ord), legNames(ord), 'Interpreter', 'none', 'Box', 'off');
lgd.Layout.Tile = 'east';
title(tl, o.titleStr)

if ~isempty(o.outFile)
    outDir = fileparts(o.outFile);
    if ~isempty(outDir) && ~isfolder(outDir); mkdir(outDir); end
    exportgraphics(gcf, o.outFile, 'Resolution', 200);
    fprintf('saved %s\n', o.outFile);
end
end
