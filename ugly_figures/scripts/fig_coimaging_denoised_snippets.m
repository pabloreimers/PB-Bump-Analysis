%% fig_coimaging_denoised_snippets
% Quality check of epg_coimage_denoised_20261118.mat (EPG jRGECO1a + EPG
% GRAB(DA2m) co-imaging, photon-count denoising + HMM bump phase): one 30 s
% snippet per trial (the window with the most turning), 10 trials per page.
% Per trial, top to bottom:
%   EPG jRGECO1a denoised z (32 PB bins: half 1 = rows 1-16, half 2 = 17-32)
%   EPG GRAB(DA) denoised z
%   ball rotation (rad/s)
% Overlays on the heatmaps: white = that channel's HMM bump phase (.mu),
% drawn on half 1 and, shifted by the trial's left-right offset, on half 2;
% red dashed (closed loop only) = cue-predicted bump position, cue(t + lag)
% with sign, lag and offset fitted to the jRGECO HMM phase. Panel title:
% trial #, id, condition, cue tracking |r| of the jRGECO HMM phase.

%% load
if exist('.data','dir') && isfile(fullfile('.data','epg_coimage_denoised_20261118.mat')), data_dir = '.data'; else, data_dir = 'data'; end
S = load(fullfile(data_dir,'epg_coimage_denoised_20261118.mat'),'all_data'); A = S.all_data(:);
a16 = linspace(-pi, pi-2*pi/16, 16)'; Wsec = 30; perPage = 10; nPage = ceil(numel(A)/perPage);
if ~exist(fullfile('ugly_figures','exports'),'dir'), mkdir(fullfile('ugly_figures','exports')); end

%% pages
for pg = 1:nPage
    figure('Color','w','Position',[30 30 1900 1150]); tl = tiledlayout(6, 5, 'TileSpacing','tight','Padding','compact');
    ks = (pg-1)*perPage + 1 : min(pg*perPage, numel(A));
    for q = 1:numel(ks)
        k = ks(q); a = A(k); xb = a.ft.xb(:)'; xf = a.ft.xf(:); n = numel(xb); fr = (n-1)/(xb(end)-xb(1));
        dark = contains(string(a.ft.pattern), "background");
        vr = interp1(xf, smoothdata(a.ft.r_speed(:),'gaussian',30), xb, 'linear', 'extrap');
        cue = interp1(xf, unwrap(a.ft.cue(:)), xb, 'linear', 'extrap');
        % cue -> bump fit on the jRGECO HMM phase (closed loop)
        mu = a.geco.mu(:)'; best = 0; lag = 0; sg = 1; off = 0;
        if ~dark
            for tt = -round(2.5*fr):round(1*fr)
                idx = 1:n; idx = idx(idx+tt >= 1 & idx+tt <= n);
                for s_ = [-1 1], m_ = mean(exp(1i*(mu(idx) + s_*cue(idx+tt)))); if abs(m_) > best, best = abs(m_); lag = tt; sg = s_; off = angle(m_); end, end
            end
        end
        % window with the most turning
        L = round(Wsec*fr); cs = cumsum(abs(vr)); tot = cs(L+1:end) - cs(1:end-L); [~, i0] = max(tot(1:max(1,end-1))); w = i0 + (0:L-1); w = w(w <= n); tw = (w - w(1))/fr;
        toRow = @(p) mod(p + pi, 2*pi)/(2*pi)*16 + 0.5; col = mod(q-1, 5) + 1;
        for c = 1:2
            ch = {'geco','grab'}; s = a.(ch{c}); Z = s.z(:, w);
            Z0 = s.z; Z0(isnan(Z0)) = 0; vL = exp(1i*a16).'*Z0(1:16,:); vR = exp(1i*a16).'*Z0(17:32,:); dLR = angle(mean(vL .* conj(vR)));
            ax = nexttile(tl, ((q > 5)*3 + c - 1)*5 + col); hold(ax,'on');
            imagesc(ax, tw, 1:32, Z); set(ax,'YDir','reverse'); clim(ax, prctile(Z(:),[1 99.5])); colormap(ax, parula); axis(ax,'tight');
            ph = s.mu(w)'; y1 = toRow(ph); y2 = toRow(ph - dLR) + 16;
            y1(abs([0 diff(y1)]) > 8) = nan; y2(abs([0 diff(y2)]) > 8) = nan;
            plot(ax, tw, y1, 'w.', 'MarkerSize', 3); plot(ax, tw, y2, 'w.', 'MarkerSize', 3);
            if ~dark
                ix = w + lag; ok = ix >= 1 & ix <= n; ref = nan(size(w)); ref(ok) = off - sg*cue(ix(ok)); r1 = toRow(ref); r1(abs([0 diff(r1)]) > 8) = nan;
                if c == 2, r1 = toRow(ref + angle(mean(exp(1i*(a.grab.mu(:)' - a.geco.mu(:)'))))); r1(abs([0 diff(r1)]) > 8) = nan; end
                plot(ax, tw, r1, 'r--', 'LineWidth', 1);
            end
            yline(ax, 16.5, 'w-', 'LineWidth', 0.5); set(ax, 'XTick', [], 'YTick', [8 24], 'YTickLabel', {'L','R'}, 'FontSize', 7);
            if c == 1
                if dark, cond = 'dark'; trk = ''; else, cond = 'CL'; trk = sprintf(', track %.2f', best); end
                title(ax, sprintf('#%d %s (%s%s)', k, regexp(a.meta,'\d{8}-\d+','match','once'), cond, trk), 'FontSize', 8, 'FontWeight', 'normal');
                ylabel(ax, 'jRGECO', 'FontSize', 7);
            else, ylabel(ax, 'GRAB', 'FontSize', 7); end
        end
        ax = nexttile(tl, ((q > 5)*3 + 2)*5 + col); plot(ax, tw, vr(w), 'k-', 'LineWidth', 1); yline(ax, 0, ':');
        xlim(ax, [0 Wsec]); set(ax, 'Box', 'off', 'FontSize', 7); ylabel(ax, 'v_r (rad/s)', 'FontSize', 7); xlabel(ax, 'time (s)', 'FontSize', 7);
    end
    title(tl, sprintf('Denoised EPG jRGECO1a / GRAB(DA) co-imaging, 30 s snippets with the most turning (page %d/%d) | white = HMM bump phase, red = cue-predicted', pg, nPage), 'FontSize', 10);
    exportgraphics(gcf, fullfile('ugly_figures','exports', sprintf('fig_coimaging_denoised_snippets_p%d.png', pg)), 'Resolution', 110);
end
