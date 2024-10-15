%Load in EEG Structure

eeglab_dir     = '/srv/TOOLKITS/eeglab-2022.0'; % https://tinyurl.com/59h6ksjs
vhtp_dir       = '/srv/Preprocessing/dependencies/vhtp'; % https://tinyurl.com/3fcbexp8
restoredefaultpath;
addpath(fullfile( eeglab_dir ));
addpath(genpath(fullfile( vhtp_dir )));
try eeglab nogui; catch, error('Check EEGLAB install'); end
raw_file_path1 = fullfile('/srv/RAWDATA/Sri_Projects/entrain_bblong/0228_BBLong.raw');
raw_file_path2 = fullfile('/srv/RAWDATA/Sri_Projects/entrain_bblong/0706_BBLong.raw');
eeg1 = pop_readegi(raw_file_path1);
eeg1.subject = 0228;
eeg1.filename = "0228_BBLong.raw";

eeg2 = pop_readegi(raw_file_path2);
eeg2.subject = 0706;
eeg2.filename = "0706_BBLong.raw";
%% run power function with EEG structure
x = eeg_htpCalcRestPower(eeg1);
%% run pipeline function with EEG structure
y = eeg_htpEegAssessPipelineHAPPE(eeg1, eeg1);