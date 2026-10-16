function fig = pb_spotcheck_browser(all_data)
% PB_SPOTCHECK_BROWSER  Pick a fly/trial from all_data and spot-check it.
%
%   fig = pb_spotcheck_browser(all_data)
%   fig = pb_spotcheck_browser()          % asks for a dataset .mat file
%
% Adds a fly dropdown and a trial dropdown on top of pb_spotcheck_gui. On
% "Load" it finds the trial's files from all_data(i).meta, loads them, and
% draws pb_spotcheck_gui into the panel below.
%
% Layout assumed (what the p2x2 / lpsp pipelines write):
%   all_data(i).meta = <date>\fly N\<trial folder>\registration
%   <trial folder>\registration\imgData_smooth_reg.mat   registered movie
%                                                         (falls back to imgData_reg.mat)
%   <trial folder>\FicTracData\fictrac-raw*.avi          FicTrac video (optional)
%
% A "fly" is the <date>\fly N path prefix (same definition as trial_fly_id in
% lpsp_p2x2_claude.m); each all_data entry is one trial of a fly. The registered
% movie and the video are read from disk (typically the Z: drive), so a trial
% whose files are not reachable reports that in the status line and keeps the
% previous view.

if nargin < 1
    [f, p] = uigetfile('*.mat', 'Select dataset containing all_data');
    if isequal(f, 0), fig = []; return; end
    S = load(fullfile(p, f), 'all_data');
    all_data = S.all_data;
end

n = numel(all_data);
fly_id = cell(n,1);
for i = 1:n, fly_id{i} = fly_of(all_data(i).meta); end
[fly_ids, ~, fly_ix] = unique(fly_id, 'stable');
fly_labels = cellfun(@fly_label, fly_ids, 'UniformOutput', false);

fig = figure('Name','PB spotcheck browser','NumberTitle','off','Color','w', ...
    'Units','normalized','Position',[0.1 0.1 0.8 0.8]);

uicontrol(fig,'Style','text','String','Fly:','Units','normalized', ...
    'Position',[0.01 0.955 0.03 0.03],'BackgroundColor','w','HorizontalAlignment','right');
ddFly = uicontrol(fig,'Style','popupmenu','String',fly_labels,'Units','normalized', ...
    'Position',[0.045 0.955 0.22 0.035],'Callback',@(~,~) onFly());
uicontrol(fig,'Style','text','String','Trial:','Units','normalized', ...
    'Position',[0.28 0.955 0.04 0.03],'BackgroundColor','w','HorizontalAlignment','right');
ddTrial = uicontrol(fig,'Style','popupmenu','String',{' '},'Units','normalized', ...
    'Position',[0.325 0.955 0.28 0.035]);
uicontrol(fig,'Style','pushbutton','String','Load','Units','normalized', ...
    'Position',[0.62 0.955 0.06 0.035],'Callback',@(~,~) onLoad());
chkVid = uicontrol(fig,'Style','checkbox','String','FicTrac video','Value',1, ...
    'Units','normalized','Position',[0.69 0.955 0.10 0.035],'BackgroundColor','w');
status = uicontrol(fig,'Style','text','String','Pick a fly and trial, then Load.', ...
    'Units','normalized','Position',[0.79 0.95 0.20 0.04],'BackgroundColor','w', ...
    'HorizontalAlignment','left');

pnl = uipanel(fig,'Position',[0 0 1 0.94],'BorderType','none','BackgroundColor','w');
trial_idx = [];     % all_data indices for the fly currently in the trial dropdown
onFly();

    function onFly()
        trial_idx = find(fly_ix == ddFly.Value);
        labels = cell(numel(trial_idx),1);
        for j = 1:numel(trial_idx)
            labels{j} = sprintf('%s  (#%d)', trial_name(all_data(trial_idx(j)).meta), trial_idx(j));
        end
        ddTrial.String = labels;
        ddTrial.Value = 1;
    end

    function onLoad()
        ti = trial_idx(ddTrial.Value);
        trial_dir = fileparts(all_data(ti).meta);
        say('Loading movie...');
        try
            imgData = load_movie(all_data(ti).meta);
        catch err
            say(['Could not load movie: ' err.message]); return
        end

        vid_file = '';
        if chkVid.Value
            v = dir(fullfile(trial_dir, 'FicTracData', 'fictrac-raw*.avi'));
            if isempty(v)
                say('No FicTrac video found; showing imaging only.');
            else
                if numel(v) > 1, warning('%d FicTrac videos in %s; using the first.', numel(v), trial_dir); end
                vid_file = fullfile(v(1).folder, v(1).name);
                say('Loading FicTrac video...');
            end
        end

        drawnow
        delete(pnl.Children)
        try
            pb_spotcheck_gui(imgData, all_data(ti), 'VideoFile', vid_file, 'Parent', pnl);
        catch err
            say(['Could not build view: ' err.message]); return
        end
        fig.Name = ['PB spotcheck: ' trial_dir];
        say(sprintf('Loaded #%d.', ti));
    end

    function say(msg)
        status.String = msg;
        drawnow
    end
end

function imgData = load_movie(meta)
% Registered movie for one trial; imgData_smooth_reg.mat preferred.
for nm = ["imgData_smooth_reg", "imgData_reg"]
    f = fullfile(meta, nm + ".mat");
    if isfile(f)
        S = load(f);
        if isfield(S, nm), imgData = S.(nm); else
            fn = fieldnames(S); imgData = S.(fn{1});
        end
        return
    end
end
error('no imgData_smooth_reg.mat / imgData_reg.mat in %s', meta);
end

function id = fly_of(meta)
% <date folder>\fly N prefix; falls back to the trial path if there is no "fly N".
parts = strsplit(meta, filesep);
parts(cellfun(@isempty, parts)) = [];
k = find(~cellfun(@isempty, regexpi(parts, '^fly\s*\d+', 'once')), 1);
if isempty(k), id = meta; else, id = strjoin(parts(1:k), filesep); end
end

function s = fly_label(id)
parts = strsplit(id, filesep);
s = strjoin(parts(max(1, end-1):end), ' / ');   % "20260727 / fly 3"
end

function s = trial_name(meta)
parts = strsplit(fileparts(meta), filesep);
s = parts{end};                                  % trial folder name
end
