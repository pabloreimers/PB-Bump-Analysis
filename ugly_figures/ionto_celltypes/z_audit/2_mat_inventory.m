% Step 2: one row per trial copy in every .data/lpsp_p2x2*.mat (meta path,
% stims, frame rate, pattern) -> p2x2_inventory.csv in the folder this is
% started from. Batch: matlab -batch "run('<path>\2_mat_inventory.m')"
work_dir = pwd;
repo_dir = 'C:\Users\ReimersPabloAlejandr\Documents\GitHub\LPsP_2p\MelData\PB-Bump-Analysis';
fs = dir(fullfile(repo_dir,'.data','lpsp_p2x2*.mat'));
inv = {};
for k = 1:numel(fs)
    t = load(fullfile(fs(k).folder,fs(k).name),'all_data'); ad = t.all_data;
    fprintf('\n=== %s: %d trials, fields: %s\n',fs(k).name,numel(ad),strjoin(fieldnames(ad),','));
    if isfield(ad,'ft') && ~isempty(ad(1).ft); fprintf('  ft fields: %s\n',strjoin(fieldnames(ad(1).ft),',')); end
    for i = 1:numel(ad)
        hasatp = isfield(ad(i),'atp') && isstruct(ad(i).atp) && isfield(ad(i).atp,'d') && ~isempty(ad(i).atp.d) && ~all(isnan(ad(i).atp.d(:)));
        nst = nan; slen = nan; fr = nan; pat = '';
        if isfield(ad(i),'ft') && isstruct(ad(i).ft)
            if isfield(ad(i).ft,'stims') && ~isempty(ad(i).ft.stims) && max(ad(i).ft.stims)>0
                st = ad(i).ft.stims(:); s = st>0.5*max(st); on = find(diff([false;s])==1); off = find(diff([s;false])==-1);
                nst = numel(on); slen = median(ad(i).ft.xf(off)-ad(i).ft.xf(on));
            end
            if isfield(ad(i).ft,'xb'); fr = 1/median(diff(ad(i).ft.xb)); end
            if isfield(ad(i).ft,'pattern'); pat = ad(i).ft.pattern; end
        end
        nrow = nan; if isfield(ad(i),'im') && isfield(ad(i).im,'d'); nrow = size(ad(i).im.d,1); end
        inv(end+1,:) = {fs(k).name, ad(i).meta, hasatp, nst, slen, fr, pat, nrow}; %#ok<SAGROW>
    end
    clear t ad
end
T = cell2table(inv,'VariableNames',{'file','meta','hasatp','nstims','stimlen','fr','pattern','nglom'});
writetable(T,fullfile(work_dir,'p2x2_inventory.csv'));
fprintf('INVENTORY DONE\n');
