
%% Adding toolboxes

% EEGLAB and VHTP toolbox paths
eeglab_dir     = '/srv/TOOLKITS/eeglab-2022.0'; % https://tinyurl.com/59h6ksjs
vhtp_dir       = '/srv/Preprocessing/dependencies/vhtp'; % https://tinyurl.com/3fcbexp8

% Load toolkits - reset matlab paths
restoredefaultpath;
addpath(fullfile( eeglab_dir ));
addpath(genpath(fullfile( vhtp_dir )));

% Run eeglab
try eeglab nogui; catch, error('Check EEGLAB install'); end



%% Get list of EEG files

% raw data folder
rawdata_folder = '/srv/RAWDATA/Sri_Projects/entrain_bblong/old_triggers/';
rawdata_filelist = dir(strcat(rawdata_folder, '*.raw'));

% preprocessed (completed) data folder big epochs
preprocessed_data_folder = '/srv/Analysis/Sri_Projects/entrain_bblong_09-04-2024/Preprocessed_files/old_BBLong_odd_only/';
%%
count = 1;
for i = 1:length(rawdata_filelist)

    disp(i)
    raw_file_path = fullfile(rawdata_filelist(i).folder, rawdata_filelist(i).name);
    eeg_examp = pop_readegi(raw_file_path);
    events = {eeg_examp.event.type};
    unique_event = unique(events);
    display(unique_event);
    
end


%% check one file and test epoching fxn


% loop through file list
for i = 1 %: length(rawdata_filelist)

% load 1st file
raw_file_path = fullfile(rawdata_filelist(i).folder, rawdata_filelist(i).name);
eeg_examp = pop_readegi(raw_file_path);

% event_conversion_struct = struct(...
%     'DI71', 'rest', ...
%     'DI61', 'sham', ...
%     'DIN7', '7_hz_stimulus', ...
%     'DIN8', '8_hz_stimulus', ...
%     'DIN9', '9_hz_stimulus', ...
%     'DI10', '10_hz_stimulus', ...
%     'DI11', '11_hz_stimulus', ...
%     'DI12', '12_hz_stimulus', ...
%     'DI13', '13_hz_stimulus');

event_conversion_struct = struct(...
    'DIN4', 'sham', ...    
    'DIN6', 'rest', ...
    'DIN7', '7_hz_stimulus', ...
    'DIN9', '9_hz_stimulus', ...
    'DI11', '11_hz_stimulus', ...
    'DI13', '13_hz_stimulus');

in = eeg_examp;
out = select_and_rename_events(in, event_conversion_struct);

end

%% Preprocessing loop

% loop through file list
parfor i = 1 : length(rawdata_filelist)

    % get file name and pathway
    raw_file_path = fullfile(rawdata_filelist(i).folder, rawdata_filelist(i).name);

    % load eeg data
    eeg_raw = pop_readegi(raw_file_path);

    % set eeg set name and file name
    [~, eeg_raw.setname, ~] = fileparts(raw_file_path);
    [~, eeg_raw.filename, ~] = fileparts(raw_file_path);


    % resample
    eeg_resamp = eeg_htpEegResampleDataEeglab(eeg_raw, 'srate', 1000);


    % filter
    eeg_filt = eeg_htpEegFilterEeglab(eeg_resamp, 'method', 'highpass', 'highpassfilt', 1, 'dynamicfiltorder', true);
    eeg_filt = eeg_htpEegFilterEeglab(eeg_filt, 'method', 'lowpass', 'lowpassfilt', 80, 'dynamicfiltorder', true);
    eeg_filt = eeg_htpEegFilterEeglab(eeg_filt, 'method', 'notch', 'notchfilt', [55, 65], 'dynamicfiltorder', true);


    % ASR cleaning (20 sd threshold)
    eeg_asr = eeg_htpEegAsrCleanEeglab(eeg_filt, ...
        'asrmode', 5, ...
        'asrflatline', 'off', ...
        'asrhighpass', 'off', ...
        'asrchannel', 'off', ...
        'asrnoisy', 'off', ...
        'asrburst', 10, ...
        'asrwindow', 'off');


    % this is where I would put ica if I wanted to


    % compare asr cleaned data
    %vis_artifacts(EEG_asr, EEG_filt)




    % select events and rename
    % event_conversion_struct = struct(...
    %     'DI71', 'rest', ...
    %     'DI61', 'sham', ...
    %     'DIN7', '7_hz_stimulus', ...
    %     'DIN8', '8_hz_stimulus', ...
    %     'DIN9', '9_hz_stimulus', ...
    %     'DI10', '10_hz_stimulus', ...
    %     'DI11', '11_hz_stimulus', ...
    %     'DI12', '12_hz_stimulus', ...
    %     'DI13', '13_hz_stimulus');
    event_conversion_struct = struct(...
        'DIN4', 'sham', ...    
        'DIN6', 'rest', ...
        'DIN7', '7_hz_stimulus', ...
        'DIN9', '9_hz_stimulus', ...
        'DI11', '11_hz_stimulus', ...
        'DI13', '13_hz_stimulus');
    eeg_asr = select_and_rename_events(eeg_asr, event_conversion_struct);




    % epoch into big blocks
    trigger_names = {"rest", "sham", "7_hz_stimulus", "8_hz_stimulus", "9_hz_stimulus", "10_hz_stimulus", "11_hz_stimulus", "12_hz_stimulus", "13_hz_stimulus"};
    epoch_window = [0 60];
    eeg_epoch = pop_epoch(eeg_asr, trigger_names, epoch_window);
    eeg_epoch = eeg_checkset(eeg_epoch);




    % remove remaining bad epochs (voltage above a certain threshold)


    % make new file name
    preprocessed_name = strcat(regexp(eeg_epoch.filename, '^\d+', 'match'), '_BBLong_preprocessed.set');
    disp(preprocessed_name)
    preprocessed_fullPath = char(fullfile(preprocessed_data_folder, preprocessed_name));


    % save the eeg structure
    pop_saveset(eeg_epoch, 'filename', preprocessed_fullPath);

end


 

%% Functions

% rename epochs function
function EEG = select_and_rename_events(EEG, event_conversion_struct)

    % delete events not specified in struct
    events_to_keep = fieldnames(event_conversion_struct);
    EEG = pop_selectevent(EEG, ...
        'type', events_to_keep, ...
        'deleteepochs', 'off', ...
        'deleteevents', 'on');

    % rename events using struct
    for i = 1:length(EEG.event)
        old_type = EEG.event(i).type;  
        EEG.event(i).type = event_conversion_struct.(old_type);
    end
    
    EEG = eeg_checkset(EEG);  

end



