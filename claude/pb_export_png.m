function outFile = pb_export_png(fig, outFile, resolution)
% PB_EXPORT_PNG  Write a figure to png at its INTENDED size, robust to locked files.
%
%   outFile = pb_export_png(fig, outFile [, resolution])
%
% Two problems this works around, both hit by the long batch runs in this repo:
%
% 1. MATLAB silently clamps a figure's on-screen size to the screen (1080 px
%    tall here), even in -batch and even for invisible figures. A 27-row grid
%    asked to be 4000 px tall comes back ~1000 px, tiledlayout squeezes every
%    tile to a sliver, and exportgraphics faithfully renders the sliver. print
%    with an explicit PaperPosition is NOT clamped (MATLAB re-lays the figure
%    out at paper size while printing), so if the caller recorded the size it
%    wanted -- setappdata(fig, 'pb_target_px', [W H]) right after
%    set(fig, 'Position', ...) -- and the figure ended up smaller than that,
%    this prints at the intended size instead. Otherwise exportgraphics is
%    used as before (tighter cropping).
% 2. On Windows a png that is open in a viewer (or being synced/scanned)
%    can't be overwritten and the export errors with "PNG library failed:
%    Could not open file" -- at the end of an hour-long run that throws away
%    everything since the last save. So: retry once after a pause, then fall
%    back to <name>_<yyyymmdd-HHMMSS>.png with a warning. Returns the path
%    actually written.

if nargin < 3; resolution = 200; end
outDir = fileparts(outFile);
if ~isempty(outDir) && ~isfolder(outDir); mkdir(outDir); end

target = getappdata(fig, 'pb_target_px');
drawnow; % Position still reports the REQUESTED size until the window manager has actually clamped it
actual = get(fig, 'Position'); actual = actual(3:4);
usePrint = ~isempty(target) && any(target > actual + 2);
if usePrint
    set(fig, 'PaperUnits', 'inches', 'PaperPositionMode', 'manual', 'PaperPosition', [0 0 target / 96]);
end

try
    write_once(fig, outFile, resolution, usePrint);
catch firstErr
    pause(2);
    try
        write_once(fig, outFile, resolution, usePrint);
    catch
        [d, n, e] = fileparts(outFile);
        alt = fullfile(d, sprintf('%s_%s%s', n, datestr(now, 'yyyymmdd-HHMMSS'), e)); %#ok<TNOW1,DATST>
        warning('pb_export_png:fallback', 'Could not write %s (%s); writing %s instead.', outFile, firstErr.message, alt);
        write_once(fig, alt, resolution, usePrint);
        outFile = alt;
    end
end
fprintf('saved %s\n', outFile);
end

function write_once(fig, outFile, resolution, usePrint)
if usePrint
    print(fig, outFile, '-dpng', sprintf('-r%d', resolution));
else
    exportgraphics(fig, outFile, 'Resolution', resolution);
end
end
