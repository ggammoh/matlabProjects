
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
rawdata_folder = '/srv/RAWDATA/Sri_Projects/entrain_bblong/';
rawdata_filelist = dir(strcat(rawdata_folder, '*.raw'));

% preprocessed (completed) data folder
preprocessed_data_folder = '/srv/Analysis/Sri_Projects/entrain_bblong_09-04-2024/Preprocessed_files/';

%% check out one file


% loop through file list
i = 1 

% get file name and pathway
raw_file_path = fullfile(rawdata_filelist(i).folder, rawdata_filelist(i).name);

% load eeg data
eeg_raw = pop_readegi(raw_file_path);



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
        'asrburst', 20, ...
        'asrwindow', 'off');


    % this is where I would put ica if I wanted to


    % compare asr cleaned data
    %vis_artifacts(EEG_asr, EEG_filt)


    % epoch into big blocks
    trigger_names = {"DI71", "DI61", "DIN7", "DIN8", "DIN9", "DI10", "DI11", "DI12", "DI13"};
    epoch_window = [0 60];
    eeg_epoch = pop_epoch(eeg_filt, trigger_names, epoch_window);
    eeg_epoch = eeg_checkset(eeg_epoch);

    % epoch big blocks further into 2 sec epochs
    event_conversion_struct = struct(...
        'DI71', 'rest', ...
        'DI61', 'sham', ...
        'DIN7', '7_hz_stimulus', ...
        'DIN8', '8_hz_stimulus', ...
        'DIN9', '9_hz_stimulus', ...
        'DI10', '10_hz_stimulus', ...
        'DI11', '11_hz_stimulus', ...
        'DI12', '12_hz_stimulus', ...
        'DI13', '13_hz_stimulus');
    eeg_epoch = make_sub_epochs_bblong(eeg_epoch, 2, event_conversion_struct);


    % remove remaining bad epochs (voltage above a certain threshold)


    % make new file name
    preprocessed_name = strcat(regexp(eeg_epoch.filename, '^\d+', 'match'), '_BBLong_preprocessed.set');
    disp(preprocessed_name)
    preprocessed_fullPath = char(fullfile(preprocessed_data_folder, preprocessed_name));


    % save the eeg structure
    pop_saveset(eeg_epoch, 'filename', preprocessed_fullPath);

end


%% TESTING 


event_conversion_struct = struct(...
    'DI64', 'stimtracker_sound_din', ...
    'DI51', 'inter_stimulus_interval', ...
    'DI71', 'rest', ...
    'DI61', 'sham', ...
    'DIN7', '7_hz_stimulus', ...
    'DIN8', '8_hz_stimulus', ...
    'DIN9', '9_hz_stimulus', ...
    'DI10', '10_hz_stimulus', ...
    'DI11', '11_hz_stimulus', ...
    'DI12', '12_hz_stimulus', ...
    'DI13', '13_hz_stimulus');



%% Check din count 

results_for_table = {};

for i = 3 %: length(rawdata_filelist)

    % get file name and pathway
    file_path = rawdata_filelist(i).folder;
    file_name = string(rawdata_filelist(i).name);
    full_file = fullfile(file_path, file_name);

    % load EEG data
    eeg_raw = pop_readegi(full_file);
    [~, eeg_raw.setname, ~] = fileparts(file_name);


    % save info
    results_for_table{i, 1} = eeg_raw.setname;

    if (isempty(eeg_raw.event))
        results_for_table{i, 2} = 0;
    else
        %results_for_table{i, 2} = size(unique({eeg_raw.event.type}), 2);
        results_for_table{i, 2} = unique({eeg_raw.event.type});

    end


end

unique_event_table = cell2table(results_for_table, 'VariableNames', {'file_name', 'dins'});




%% Test epoching function
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
% out = make_sub_epochs_bblong(eeg_epoch, 2, event_conversion_struct)



%% Functions


% 2 second epoch split function
function eeg_new_epoch = split_epochs(eeg_epoch,  new_epoch_length)
    % Inputs:
    %   eeg_epoch - EEG structure containing the original epoched data
    %   trigger_names - Cell array of trigger names to preserve
    %   new_epoch_length - Desired length of new epochs in seconds
    %
    % Output:
    %   eeg_new_epoch - EEG structure containing the new, shorter epochs
    
    % Convert epoched data back to continuous
    eeg_cont = eeg_epoch2continuous(eeg_epoch);
    eeg_cont = eeg_checkset(eeg_cont);
    
    % Create new epochs
    eeg_new_epoch = eeg_regepochs(eeg_cont, 'recurrence', new_epoch_length, 'limits', [0 new_epoch_length], 'rmbase', NaN);
    
    % Update the EEG structure
    eeg_new_epoch = eeg_checkset(eeg_new_epoch);


end


% Epoch big blocks into 2 second epochs
function EEG_epoch = make_sub_epochs_bblong(EEG, epoch_duration_sec, event_conversion_struct)
    
    
    og_epoch_start = EEG.xmin * EEG.srate; % miliseconds
    og_epoch_end = EEG.xmax * EEG.srate;
    epoch_duration_ms = epoch_duration_sec * EEG.srate;
    
    % loop through origina lepochs
    for og_epoch_idx = 1 : length(EEG.epoch)
    
        cur_epoch_type = EEG.event(og_epoch_idx).type;
        cur_epoch_latency = EEG.event(og_epoch_idx).latency;
    
        disp(cur_epoch_type)
    
        for new_epoch_time = og_epoch_start : epoch_duration_ms : og_epoch_end - epoch_duration_ms
    
            new_event_type = event_conversion_struct.(cur_epoch_type);
            new_event = struct('type', new_event_type, ...
                'latency', cur_epoch_latency + new_epoch_time, ...
                'urevent', length(EEG.event)+1, ...
                'epoch', og_epoch_idx);
    
            EEG.event(end+1) = new_event;
    
        end
    
    end
    
    
    % have to select 1 sample less to ensure epochs don't overlap
    EEG_epoch = pop_epoch(EEG, ...
        struct2cell(event_conversion_struct), ...
        [0 epoch_duration_sec - (1 / EEG.srate)], ...
        'newname', EEG.filename, ...
        'epochinfo', 'no');
    

    % remove original events
    EEG_epoch = pop_selectevent(EEG_epoch, 'type', struct2cell(event_conversion_struct), 'deleteevents', 'on', 'deleteepochs', 'off');
    
    EEG_epoch = eeg_checkset(EEG_epoch);
    
    
    disp('success yay')

end