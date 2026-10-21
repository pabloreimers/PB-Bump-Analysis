%% fig_mobility_example_candidates
% Scratch: pick the example fly + time windows for Figure 3's bump-mobility
% example row (closed-loop pre-ATP tracking, then the same fly's dark
% pre- and post-ATP trials). Same idea as Figure_3_scratch.m section 10
% (figure E), but scored on Figure_3_claude.m's own row-6 numbers
% (walk_fly_pre/post, same bouts) and with automatic window picks:
%
%  CLOSED LOOP (pre-ATP): sliding win_s windows; score = bump-vs-heading
%    VELOCITY correlation (0.5 s bins -- tracking of turns, not of a shared
%    slow drift), only windows with >= min_bouts walking bouts and mean rho
%    above the trial median. Heading = unwrap(-ft.cue).
%  DARK pre: same velocity-correlation score as closed loop -- the point
%    of the panel is "the bump follows self-motion in the dark before ATP".
%    (A first version matched the window's path-length ratio to the trial
%    slope instead; that picked windows where bump JITTER made up the path
%    length while the bump didn't follow the fly's turns at all -- path
%    length alone isn't a visual tracking measure.)
%  DARK post: windows where the fly turns at least min_turn rad in total
%    (so the ratio means something), picking the one whose own bump/heading
%    path-length ratio (Figure 3's metric, summed over the window's bouts)
%    is closest to that TRIAL's pooled sqrt(dur)-weighted slope --
%    representative, not the most dramatic. Heading = cumsum(r_speed)*dt.
%
% Walking bouts use the Figure 3 definition (walk_trial_bouts parameters).
% Traces are re-zeroed at the window start (mean of its first second), so
% a dark-trial drift between bump and heading stays visible.
%
% REQUIRES the Figure_3_claude.m workspace (walk_* variables).

%% check workspace
need = {'walk_data','walk_fly_ix','walk_fly_list','walk_trial_num','walk_has_pulse','walk_is_dark','walk_fly_group','walk_fly_pre','walk_fly_post','walk_fly_rho_pre','walk_fly_rho_post','walk_smooth_window','walk_turn_thresh','walk_max_gap_frames','walk_min_walking_frames'};
missing = {}; for v = need; if exist(v{1},'var')~=1; missing{end+1} = v{1}; end; end %#ok<AGROW>
if ~isempty(missing); error('run Figure_3_claude.m first -- missing: %s', strjoin(missing,', ')); end
repo_root  = fileparts(fileparts(fileparts(mfilename('fullpath'))));
export_dir = fullfile(repo_root,'ugly_figures','exports');

%% candidates + window parameters
cand_tags = {'20260708/fly 1','20260723/fly 2','20260727/fly 4','20260728/fly 9'}; % lpsp, ranked by CL tracking + representative dark change (see chat 2026-10-20)
win_s = 40; min_bouts = 4; step_s = 1; min_turn = 2*pi;

%% score windows
C = struct('tag',{},'f',{},'cl',{},'dpre',{},'dpost',{});
for k = 1:numel(cand_tags)
    f = find(contains(walk_fly_list, cand_tags{k})); assert(isscalar(f), 'fly %s not found', cand_tags{k});
    ti = find(walk_fly_ix==f); [~,o] = sort(walk_trial_num(ti)); ti = ti(o);
    p = find(walk_has_pulse(ti),1);
    pre = ti(1:p-1); post = ti(p+1:end); post = post(~walk_has_pulse(post));
    c = struct('tag',cand_tags{k},'f',f);
    % closed loop: best window over all pre-ATP CL trials
    c.cl = best_window(walk_data, pre(~walk_is_dark(pre)), 'cl', win_s, step_s, min_bouts, walk_smooth_window, walk_turn_thresh, walk_max_gap_frames, walk_min_walking_frames, min_turn);
    c.dpre  = best_window(walk_data, pre(walk_is_dark(pre)),   'darkpre', win_s, step_s, min_bouts, walk_smooth_window, walk_turn_thresh, walk_max_gap_frames, walk_min_walking_frames, min_turn);
    c.dpost = best_window(walk_data, post(walk_is_dark(post)), 'darkpost', win_s, step_s, min_bouts, walk_smooth_window, walk_turn_thresh, walk_max_gap_frames, walk_min_walking_frames, min_turn);
    for cn = {'dpre','dpost'} % fly never turned min_turn in any window -> relax the turning gate rather than drop the panel
        if isnan(c.(cn{1}).i)
            warning('%s %s: no window turns >= %.1f rad -- relaxing to 0', c.tag, cn{1}, min_turn);
            kinds = struct('dpre','darkpre','dpost','darkpost'); src = struct('dpre',pre(walk_is_dark(pre)),'dpost',post(walk_is_dark(post)));
            c.(cn{1}) = best_window(walk_data, src.(cn{1}), kinds.(cn{1}), win_s, step_s, min_bouts, walk_smooth_window, walk_turn_thresh, walk_max_gap_frames, walk_min_walking_frames, 0);
        end
    end
    C(k) = c;
    fprintf('%-16s CL trial %d t0=%4.0fs vel-r=%.2f (%d bouts) | dark pre trial %d t0=%4.0fs vel-r=%.2f ratio %.2f vs trial %.2f (%d bouts) | dark post trial %d t0=%4.0fs ratio %.2f vs trial %.2f (%d bouts)\n', ...
        c.tag, walk_trial_num(c.cl.i), c.cl.t0, c.cl.score, c.cl.nb, ...
        walk_trial_num(c.dpre.i), c.dpre.t0, c.dpre.score, c.dpre.ratio, c.dpre.trial_slope, c.dpre.nb, ...
        walk_trial_num(c.dpost.i), c.dpost.t0, c.dpost.ratio, c.dpost.trial_slope, c.dpost.nb);
end

%% gallery: one row per candidate fly; CL pre | dark pre | dark post
fig = figure('color','w','Visible','off','Position',[40 40 1800 330*numel(C)]);
tl = tiledlayout(fig,numel(C),3,'TileSpacing','compact','Padding','compact');
cond = {'cl','dpre','dpost'}; cond_lab = {'closed loop, pre-ATP','dark, pre-ATP','dark, post-ATP'};
for k = 1:numel(C)
    for j = 1:3
        w = C(k).(cond{j}); ax = nexttile(tl); hold(ax,'on')
        plot_window(ax, walk_data(w.i), w, win_s, walk_smooth_window, walk_turn_thresh, walk_max_gap_frames, walk_min_walking_frames)
        if j==1
            ttl = sprintf('%s -- %s, trial %d, t=%.0f-%.0f s\nbump-vs-heading velocity r = %.2f, %d bouts', C(k).tag, cond_lab{j}, walk_trial_num(w.i), w.t0, w.t0+win_s, w.score, w.nb);
        else
            ttl = sprintf('%s, trial %d, t=%.0f-%.0f s\nvel r %.2f, window bump/heading path %.2f (trial slope %.2f), %d bouts, rho %.2f', cond_lab{j}, walk_trial_num(w.i), w.t0, w.t0+win_s, w.vel_r, w.ratio, w.trial_slope, w.nb, w.rho);
        end
        title(ax, ttl, 'FontSize', 8, 'Interpreter','none')
        if k==numel(C); xlabel(ax,'time in window (s)'); end
    end
    title(tl, sprintf('Figure 3 mobility-row example candidates (lpsp>P2X2), %d s windows. black = bump (im.mu), blue = heading, grey = walking bouts. Fly mobility pre -> post per row in the command window.', win_s), 'FontSize', 10)
end
exportgraphics(fig, fullfile(export_dir,'fig_mobility_example_candidates.png'), 'Resolution', 130)
fig.Visible = 'on';

%% local functions
function [is_walking, r_s, dt, xf] = walk_mask(tr, sm, thr, max_gap, min_walk)
    % Figure 3 walking-bout mask (walk_trial_bouts), on the ft.xf clock
    xf = tr.ft.xf(:); dt = median(diff(xf));
    r_s = smoothdata(tr.ft.r_speed(:),'gaussian',sm);
    is_walking = abs(r_s) > thr;
    is_walking = imclose(is_walking,ones(max_gap+1,1));
    is_walking = bwareaopen(is_walking,min_walk);
end

function [head, mu] = window_traces(tr, kind, xf, dt)
    % heading and unwrapped bump on the xf clock (not yet re-zeroed)
    if strcmp(kind,'cl'); head = unwrap(-tr.ft.cue(:)); else; head = cumsum(tr.ft.r_speed(:))*dt; end
    mu = interp1(tr.ft.xb(:), unwrap(tr.im.mu(:)), xf, 'linear', 'extrap');
end

function W = best_window(data, trials, kind, win_s, step_s, min_bouts, sm, thr, max_gap, min_walk, min_turn)
    W = struct('i',nan,'t0',nan,'score',-inf,'nb',0,'ratio',nan,'trial_slope',nan,'rho',nan,'vel_r',nan);
    best_cost = inf;
    for i = trials(:)'
        tr = data(i);
        [iw, r_s, dt, xf] = walk_mask(tr, sm, thr, max_gap, min_walk);
        [head, mu] = window_traces(tr, kind(1:min(2,end)), xf, dt);
        mu_s = smoothdata(mu,'gaussian',sm);
        rho_f = interp1(tr.ft.xb(:), tr.im.rho(:), xf, 'linear', 'extrap');
        d = diff([false; iw; false]); b0 = find(d==1); b1 = find(d==-1)-1;
        bump_path = arrayfun(@(b) sum(abs(diff(mu_s(b0(b):b1(b))))), 1:numel(b0))';
        head_path = arrayfun(@(b) sum(abs(r_s(b0(b):b1(b))))*dt, 1:numel(b0))';
        dur = (b1-b0+1)*dt;
        trial_slope = (sqrt(dur).*head_path) \ (sqrt(dur).*bump_path);
        rho_med = median(rho_f,'omitnan');
        nwin = round(win_s/dt); bin = round(.5/dt);
        for s0 = 1:round(step_s/dt):numel(xf)-nwin
            s1 = s0+nwin-1;
            inb = b0>=s0 & b1<=s1; nb = sum(inb);
            if nb < min_bouts; continue; end
            rho_w = mean(rho_f(s0:s1),'omitnan');
            hb = head(s0:bin:s1); mb = mu_s(s0:bin:s1);
            r = corr(diff(hb), diff(mb)); % bump-vs-heading velocity correlation, 0.5 s bins
            if strcmp(kind,'darkpre') && sum(head_path(inb)) < min_turn; continue; end % tracking needs turns to show
            if any(strcmp(kind,{'cl','darkpre'}))
                if rho_w < rho_med; continue; end
                if r > W.score
                    W = struct('i',i,'t0',xf(s0),'score',r,'nb',nb,'ratio',sum(bump_path(inb))/sum(head_path(inb)),'trial_slope',trial_slope,'rho',rho_w,'vel_r',r);
                end
            else
                if sum(head_path(inb)) < min_turn; continue; end
                ratio = sum(bump_path(inb))/sum(head_path(inb));
                cost = abs(ratio - trial_slope) - .01*nb; % representative; ties -> more bouts
                if cost < best_cost
                    best_cost = cost;
                    W = struct('i',i,'t0',xf(s0),'score',-cost,'nb',nb,'ratio',ratio,'trial_slope',trial_slope,'rho',rho_w,'vel_r',r);
                end
            end
        end
    end
end

function plot_window(ax, tr, w, win_s, sm, thr, max_gap, min_walk)
    [iw, ~, dt, xf] = walk_mask(tr, sm, thr, max_gap, min_walk);
    kind = 'dark'; if contains(tr.ft.pattern,'bar') && ~contains(tr.ft.pattern,'background'); kind = 'cl'; end
    [head, mu] = window_traces(tr, kind, xf, dt);
    k = xf>=w.t0 & xf<=w.t0+win_s; t = xf(k)-w.t0;
    h = head(k); m = mu(k); z = t<=1;
    h = h - mean(h(z)); m = m - mean(m(z));
    d = diff([false; iw(k); false]); b0 = find(d==1); b1 = find(d==-1)-1;
    yl = max([abs(h); abs(m); pi])*1.1;
    for b = 1:numel(b0); patch(ax, t([b0(b) b1(b) b1(b) b0(b)]), [-yl -yl yl yl], [.9 .9 .9], 'EdgeColor','none'); end
    plot(ax, t, h, '-', 'Color', [0 .2 1], 'LineWidth', 1.3)
    plot(ax, t, m, 'k-', 'LineWidth', 1.3)
    xlim(ax,[0 win_s]); ylim(ax,[-yl yl])
    step = pi*max(1, ceil(yl/(3*pi))); kk = -floor(yl/step):floor(yl/step);
    yticks(ax, kk*step); yticklabels(ax, arrayfun(@(q) pi_lab(q*step/pi), kk, 'UniformOutput', false))
    set(ax,'FontSize',8)
end

function s = pi_lab(m)
    if m == 0; s = '0'; elseif m == 1; s = '\pi'; elseif m == -1; s = '-\pi'; else; s = sprintf('%d\\pi',m); end
end
