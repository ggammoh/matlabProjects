%Load in EEG Structure

eeglab_dir     = '/srv/TOOLKITS/eeglab-2022.0'; % https://tinyurl.com/59h6ksjs
vhtp_dir       = '/srv/Preprocessing/dependencies/vhtp'; % https://tinyurl.com/3fcbexp8
restoredefaultpath;
addpath(fullfile( eeglab_dir ));
addpath(genpath(fullfile( vhtp_dir )));
try eeglab nogui; catch, error('Check EEGLAB install'); end
raw_file_path = fullfile('/srv/RAWDATA/Sri_Projects/entrain_bblong/0228_BBLong.raw');
eeg = pop_readegi(raw_file_path);
eeg.subject = 3416;
eeg.filename = "0228_BBLong.raw";
%% run power function without EEG structure
x = eeg_htpCalcRestPower(eeg);
%%run power function with EEG structure