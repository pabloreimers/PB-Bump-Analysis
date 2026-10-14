function s = pb_h5_sync(h5Path)
% PB_H5_SYNC  Timing / stimulus channels from a berg4 WaveSurfer .h5.
%
%   s = pb_h5_sync(h5Path)
%
% Returns, on the fine DAQ clock (s.fs, s.t_fine) and the imaging volume
% clock (s.t_volume = si_volumeclock rising edges):
%   s.led        logical, stim LED on, per DAQ sample
%   s.led_on / s.led_off   row vectors of pulse onset / offset times (s)
%   s.ledSource  'DI led_stim_feedback' or 'AI vr_led_stim_copy' (see below)
%   s.bar_rad    bar angle on the LED display (rad; vr_ao0_copy is 1 V/rad)
%   s.vr_start / s.vr_end  first / last time the bar channel is > 0.2 V
%
% The stim LED has been logged two ways on berg4: as the digital line
% led_stim_feedback (trials up to Oct 5 2026) or, from Oct 12 2026, only as
% the analog command copy vr_led_stim_copy (0 / 5 V). This reads whichever
% exists (digital preferred; analog thresholded at 2.5 V). If neither exists
% s.led is all false and led_on/led_off are empty.
%
% Same conventions as obs_h5_sync in data scripts\overshoot_bar_stim_script.m
% (vr_ao0_copy calibration: slope 0.998, offset 0.011 V vs the VR's own log).

info = h5info(h5Path);
sweepNames = {info.Groups.Name};
sweepNames = sweepNames(~strcmp(sweepNames, '/header'));
assert(numel(sweepNames) == 1, 'pb_h5_sync: expected one sweep group in %s, found %d', h5Path, numel(sweepNames));
sw = sweepNames{1};

s.fs = double(h5read(h5Path, '/header/AcquisitionSampleRate'));
diNames = strtrim(string(h5read(h5Path, '/header/DIChannelNames')));
aiNames = strtrim(string(h5read(h5Path, '/header/AIChannelNames')));
coef = double(h5read(h5Path, '/header/AIScalingCoefficients'));   % c0..c3 per AI channel, counts -> volts

digi = int32(h5read(h5Path, [sw '/digitalScans']));
ana = h5read(h5Path, [sw '/analogScans']);
if size(ana, 1) < size(ana, 2), ana = ana'; end
toVolts = @(ch) polyval(flipud(coef(:, ch)), double(ana(:, ch)));

volCh = find(diNames == "si_volumeclock", 1);
barCh = find(aiNames == "vr_ao0_copy", 1);
assert(~isempty(volCh) && ~isempty(barCh), 'pb_h5_sync: missing si_volumeclock / vr_ao0_copy in %s (DI: %s; AI: %s)', ...
    h5Path, strjoin(diNames, ', '), strjoin(aiNames, ', '));

N = numel(digi);
s.t_fine = (0:N-1)' / s.fs;
s.t_volume = s.t_fine(find(diff(bitget(digi, volCh)) > 0) + 1);

ledDi = find(diNames == "led_stim_feedback", 1);
ledAi = find(aiNames == "vr_led_stim_copy", 1);
if ~isempty(ledDi)
    s.led = bitget(digi, ledDi) > 0; s.ledSource = 'DI led_stim_feedback';
elseif ~isempty(ledAi)
    s.led = toVolts(ledAi) > 2.5; s.ledSource = 'AI vr_led_stim_copy (> 2.5 V)';
else
    s.led = false(N, 1); s.ledSource = 'none';
end
onI  = find(diff([0; double(s.led)]) == 1);
offI = find(diff([double(s.led); 0]) == -1);
m = min(numel(onI), numel(offI));
s.led_on  = s.t_fine(onI(1:m))';
s.led_off = s.t_fine(offI(1:m))';

s.bar_rad = toVolts(barCh);                                  % volts == radians for this channel
s.vr_start = s.t_fine(find(s.bar_rad > 0.2, 1));
s.vr_end   = s.t_fine(find(s.bar_rad > 0.2, 1, 'last'));
end
