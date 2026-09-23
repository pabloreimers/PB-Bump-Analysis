function s = pb_h5_stim_timing(h5Path)
% PB_H5_STIM_TIMING  LED stim logical + per-volume timestamps from a WaveSurfer h5.
%
%   s = pb_h5_stim_timing(h5Path)
%
% Reads Janelia WaveSurfer's HDF5 layout (/header + one /sweep_000N group;
% the sweep name is found dynamically since a restarted acquisition can be
% /sweep_0002) and pulls two digital lines out of /sweep_000N/digitalScans:
%   si_volumeclock     one pulse per completed ScanImage volume -- its rising
%                      edges ARE the volume timestamps on the DAQ clock (see
%                      data scripts/led_stim_berg4_script.m's header for the
%                      verification against the tif's own frame timestamps)
%   led_stim_feedback  the actual LED TTL
%
% Returns
%   s.fs          DAQ sample rate (Hz)
%   s.durationSec real recording duration (N samples / fs)
%   s.t_volume    [nVol x 1] volume timestamps (s, DAQ clock)
%   s.stims       [nVol x 1] logical, LED on/off downsampled onto the volume
%                 clock (nearest sample) -- the only resolution the glomerulus
%                 analyses need
%
% Shared by data scripts/overshoot_glom_delta_script.m and
% led_stim_berg4_glom_delta_script.m (each of the older scripts in data scripts/
% carries its own private copy of this logic; this is the same code).

info = h5info(h5Path);
sweepNames = {info.Groups.Name};
sweepNames = sweepNames(~strcmp(sweepNames, '/header'));
if numel(sweepNames) ~= 1
    error('pb_h5_stim_timing:sweep', 'Expected exactly one sweep group in %s, found %d.', h5Path, numel(sweepNames));
end

fs      = double(h5read(h5Path, '/header/AcquisitionSampleRate'));
diNames = strtrim(string(h5read(h5Path, '/header/DIChannelNames')));
digi    = int32(h5read(h5Path, [sweepNames{1} '/digitalScans']));

volCh = find(diNames == "si_volumeclock", 1);
ledCh = find(diNames == "led_stim_feedback", 1);
if isempty(volCh) || isempty(ledCh)
    error('pb_h5_stim_timing:chan', ...
        'Could not find si_volumeclock/led_stim_feedback in %s DIChannelNames (found: %s).', h5Path, strjoin(diNames, ', '));
end

volBit = bitget(digi, volCh);
ledBit = bitget(digi, ledCh);
N      = numel(digi);
t_fine = (0:N-1)' / fs;

s.fs          = fs;
s.durationSec = N / fs;
s.t_volume    = t_fine(find(diff(volBit) > 0) + 1);
s.stims       = interp1(t_fine, double(ledBit), s.t_volume, 'nearest', 'extrap') > 0.5;
end
