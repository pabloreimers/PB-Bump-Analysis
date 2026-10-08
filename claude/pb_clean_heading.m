function [hClean, info] = pb_clean_heading(h, opts)
% PB_CLEAN_HEADING  Remove single-frame fictrac heading glitches (spikes and steps).
%
%   [hClean, info] = pb_clean_heading(h)          h: unwrapped heading (rad), uniformly sampled
%   [hClean, info] = pb_clean_heading(h, 'minThreshRad', 0.25, 'nMad', 8, 'growFrames', 1)
%
% Fictrac occasionally loses the ball for a frame or two and re-acquires at a
% different heading. In the unwrapped heading this looks like either a SPIKE
% (jump out and straight back, e.g. 160 deg at the start of 20230109 fly 1
% trial 1 while fictrac locks on) or a STEP (jump that never comes back, e.g.
% 65 deg in one frame in 20221221 fly 3 trial 1). Real turns never move the
% heading by more than a few degrees per 60-Hz frame (a fast saccade peaks
% near 1000 deg/s = 17 deg/frame and ramps over several frames), so a
% single-frame increment far above the trial's own spread is not behaviour.
%
% Method: work on the increments dh = diff(h). An increment is a glitch if
% |dh| > max(minThreshRad, nMad * 1.4826 * MAD(dh)); flagged frames are grown
% by growFrames on each side and grouped into runs. Each run is then either
%   "wrap": its net excursion is within wrapTolRad of k*2*pi (k ~= 0) -- a
%           2*pi wrap that unwrap() missed because it landed on noisy frames
%           (20230109 fly 1 trial 1, 5 s: ~350 deg in a few frames). Exactly
%           k turns are subtracted from the largest increment; the residual
%           motion is kept.
%   "loss": anything else -- fictrac lost the ball. The increments in the run
%           are linearly interpolated from their neighbours, so a spike's
%           out-and-back cancels and a step is removed entirely; everything
%           after a removed step shifts by the step size, which is intended:
%           the heading stays continuous.
% The heading is then re-integrated from the cleaned increments.
%
% Limitation: a genuine saccade peaking above ~minThreshRad per frame
% (~860 deg/s at 60 Hz) can be clipped. In the 2026-10-07 epg_7f batch this
% affected only 20230616 fly 2 (18-25 deg/frame events, < 33 deg total change).
%
% info: thresh (rad/frame), nBad (frames replaced), segs (n x 2 sample
% indices of each run), kind ("wrap"/"loss" per run), netRemoved (rad, total
% heading change taken out), madRad.

arguments
    h (1,:) double
    opts.minThreshRad (1,1) double = 0.25   % ~14 deg/frame = 860 deg/s at 60 Hz
    opts.nMad (1,1) double = 8
    opts.growFrames (1,1) double = 1
    opts.wrapTolRad (1,1) double = 0.6      % a run whose net excursion is within this of k*2*pi is a missed wrap
end

dh = diff(h);
madRad = 1.4826 * median(abs(dh - median(dh, 'omitnan')), 'omitnan');
thresh = max(opts.minThreshRad, opts.nMad * madRad);
bad = abs(dh) > thresh;
if opts.growFrames > 0
    bad = conv(double(bad), ones(1, 2 * opts.growFrames + 1), 'same') > 0;
end
dhClean = dh;
d = diff([false bad false]); s = find(d == 1); e = find(d == -1) - 1;
kind = strings(numel(s), 1);
good = ~bad & ~isnan(dh);
for r = 1:numel(s)
    run = s(r):e(r);
    net = sum(dh(run));
    k = round(net / (2*pi));
    if k ~= 0 && abs(net - 2*pi*k) < opts.wrapTolRad
        % a missed 2*pi wrap (the wrap landed on noisy frames, so unwrap did
        % not catch it): take exactly k turns out of the largest increment
        % and keep the small residual motion
        [~, im] = max(abs(dh(run))); dhClean(run(im)) = dh(run(im)) - 2*pi*k;
        kind(r) = "wrap";
    else
        % tracking loss (spike or step): the heading cannot have moved this
        % much in a frame; interpolate the increments across the run
        dhClean(run) = interp1(find(good), dh(good), run, 'linear', 'extrap');
        kind(r) = "loss";
    end
end
hClean = h(1) + [0 cumsum(dhClean)];

info = struct('thresh', thresh, 'nBad', nnz(bad), 'segs', [s(:) e(:)], 'kind', kind, ...
    'netRemoved', sum(dh(bad) - dhClean(bad)), 'madRad', madRad);
end
