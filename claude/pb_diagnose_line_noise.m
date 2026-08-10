function stats = pb_diagnose_line_noise(movie, options)
%PB_DIAGNOSE_LINE_NOISE Plot the fast-scan-axis spectrum to find/verify a noise stopband.
%
%   stats = pb_diagnose_line_noise(movie) averages the 1D power spectrum along
%   dimension 2 (X, the fast-scan axis) over every row/plane/frame in movie and
%   plots it, with the single largest non-DC peak annotated. Run this on raw
%   data to pick a stopband for PB_REMOVE_LINE_NOISE, and again on the filtered
%   output to confirm the peak is gone -- do not just trust the default [10 20]
%   band, since this is one fly/one rig/one day's electrical environment and it
%   can drift.
%
%   Name-value options:
%     maxFrames    (1,1) double = 200     subsample to at most this many frames
%                                          (taken evenly across dim 3+) for speed
%     figureVisible(1,1) logical = true
%     title        (1,:) char    = ''     prepended to the figure title, e.g.
%                                          'raw' or 'after filtering'
%
%   stats fields: freq_cyc_per_line, magnitude (both half-spectrum, DC first),
%   peakFreq, peakMagnitude, peakProminenceRatio (peak / local baseline -- values
%   well above ~2-3 indicate a real narrowband artifact worth notching).
%
%   See also PB_REMOVE_LINE_NOISE.

arguments
    movie {mustBeNumeric}
    options.maxFrames (1,1) double {mustBePositive} = 200
    options.figureVisible (1,1) logical = true
    options.title (1,:) char = ''
end

sz = size(movie);
nY = sz(1);
nX = sz(2);
nTrailing = prod(sz(3:end));
movieR = reshape(movie, nY, nX, nTrailing);

nUse = min(options.maxFrames, nTrailing);
frameIdx = round(linspace(1, nTrailing, nUse));

lines = reshape(permute(double(movieR(:,:,frameIdx)), [2 1 3]), nX, []); % X x (Y*nUse)
lines = lines - mean(lines, 1);
Lf = abs(fft(lines, [], 1));
Lf_avg = mean(Lf, 2);

half = 1:(floor(nX/2) + 1);
freq = (0:nX-1)'; % FFT bin k = k cycles across the full nX-px line
freq = freq(half);
mag = Lf_avg(half);

% find largest peak excluding DC and its immediate skirt
searchIdx = 3:numel(mag)-1; % skip DC (index 1) and index 2
[peakMag, relIdx] = max(mag(searchIdx));
peakIdx = searchIdx(relIdx);
peakFreq = freq(peakIdx);
baseline = median(mag(searchIdx));

stats.freq_cyc_per_line = freq;
stats.magnitude = mag;
stats.peakFreq = peakFreq;
stats.peakMagnitude = peakMag;
stats.peakProminenceRatio = peakMag / baseline;

if options.figureVisible
    figHandle = figure('Name', 'pb_diagnose_line_noise');
else
    figHandle = figure('Visible', 'off');
end
figure(figHandle);
plot(freq, mag, 'LineWidth', 1.2);
hold on
plot(peakFreq, peakMag, 'rv', 'MarkerFaceColor', 'r');
text(peakFreq, peakMag, sprintf('  %.1f cyc/line (%.1fx baseline)', peakFreq, stats.peakProminenceRatio), ...
    'VerticalAlignment', 'bottom');
grid on
xlabel('spatial frequency (cycles per line)')
ylabel('avg |FFT| magnitude (DC excluded from axis scaling)')
ylim([0, max(mag(searchIdx)) * 1.3])
titleStr = 'fast-scan-axis spectrum';
if ~isempty(options.title)
    titleStr = sprintf('%s -- %s', options.title, titleStr);
end
title(titleStr, 'Interpreter', 'none')

fprintf('pb_diagnose_line_noise: largest non-DC peak at %.1f cyc/line, %.1fx the median baseline magnitude.\n', ...
    peakFreq, stats.peakProminenceRatio);
if stats.peakProminenceRatio < 2
    fprintf('  (no strong narrowband peak found -- spectrum looks close to flat/broadband)\n');
end

end
