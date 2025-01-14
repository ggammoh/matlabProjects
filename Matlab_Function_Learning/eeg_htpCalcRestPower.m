function [EEG, results] = eeg_htpCalcRestPower(EEG, varargin)
    % Description: calculate spectral power on continuous data.
    % Category: Analysis
    % ShortTitle: Resting Spectral Power Analysis
    % Tags: Power
    %      Power is calculated using MATLAB pWelch function. Key parameter is
    %      window length with longer window providing increased frequency
    %      resolution. Overlap is set at default at 50%. A hanning window is
    %      also implemented. Speed is greatly increased by GPU.
    %
    % Usage:
    %    >> [ EEG, results ] = eeg_htpCalcRestPower( EEG, varargin )
    %
    % Require Inputs:
    %     EEG       - EEGLAB Structure
    % Function Specific Inputs:
    %     gpuon     - [logical] use gpuArray. default: false
    %     duration  - [integer] duration to calculate on. default: 80 seconds
    %                 if duration is greater sample, will default to max size.
    %     offset    - [integer] start time in seconds. default: 0
    %     outputdir - default, same as EEG path
    %.    useParquet - [logical] use Parquet format on save
    % Common Visual HTP Inputs:
    %     'bandDefs'   - cell-array describing frequency band definitions
    %     {'delta', 2 ,3.5;'theta', 3.5, 7.5; 'alpha1', 8, 10; 'alpha2', 10.5, 12.5;
    %     'beta', 13, 30;'gamma1', 30, 55; 'gamma2', 65, 80; 'epsilon', 81, 120;}
    %     'outputdir' - path for saved output files (default: tempdir)
    %
    % Outputs:
    % Outputs:
    %     EEG       - EEGLAB Structure with modified .vhtp field
    %                 [table] summary_table: subject chan power_type_bandname
    %                 [table] spectro: channel average power for spectrogram
    %     results   - .vhtp structure
    %
    %  This file is part of the Cincinnati Visual High Throughput Pipeline,
    %  please see http://github.com/cincibrainlab
    %
    %  Contact: kyle.cullion@cchmc.org

    timestamp = datestr(now, 'yymmddHHMMSS'); % timestamp
    functionstamp = mfilename; % function name for logging/output

    [note] = htp_utilities();

    % Inputs: Function Specific
    defaultGpu = 0;%variable to set GPU acceleration on or off
    defaultDuration = 60; %possible used to determine length of EEG file
    defaultOffset = 0;%probably used tg set offsets for data indexing
    defaultWindow = 2;% window size for data
    defaultUseParquet = false;% determines if data will be set in parquet format

    % Inputs: Common across Visual HTP functions
    defaultOutputDir = '/srv/Analysis/Sri_Projects/Matlab_Function_Learning/eeg_htpCalcRestPower/';
    defaultBandDefs = {'delta', 2, 3.5; 'theta', 3.5, 7.5; 'alpha1', 8, 10;
                    'alpha2', 10.5, 12.5; 'beta', 13, 30; 'gamma1', 30, 55;
                    'gamma2', 65, 80; 'epsilon', 81, 120; };% cell array where each row defines
    %EEG frequency bands, with its name first and then its range in hz.

    % MATLAB built-in input validation
    ip = inputParser(); %creates input parser object used in addparameter function
    addRequired(ip, 'EEG', @isstruct);%checks to see if the EEG file is a structure
    addParameter(ip, 'gpuOn', defaultGpu, @mustBeNumericOrLogical);%all of these check if the corresponding value has the correct data type
    addParameter(ip, 'duration', defaultDuration, @isnumeric);
    addParameter(ip, 'offset', defaultOffset, @isnumeric);
    addParameter(ip, 'window', defaultWindow, @isnumeric);
    addParameter(ip, 'outputdir', defaultOutputDir, @isfolder);
    addParameter(ip, 'bandDefs', defaultBandDefs, @iscell);
    addParameter(ip, 'useParquet', defaultUseParquet, @islogical);

    parse(ip, EEG, varargin{:});%parses the 2 structures

    outputdir = '/srv/Analysis/Sri_Projects/Matlab_Function_Learning/eeg_htpCalcRestPower/';%pulling out bandefs and outputdir from the inputParser results
    bandDefs = ip.Results.bandDefs;
    fprintf("Finished Parsing")

    % File Management (create subfolder with function name)
    [~, basename, ~] = fileparts(EEG.filename);%acquires file name
    analysis_outputdir =  fullfile(outputdir, mfilename);%creates a file path for anaylysis of output directory
    if ~exist("analysis_outputdir", "dir")%checks to see if analysis_outputdir exists if not then it creates the directory
        mkdir(analysis_outputdir);
    end
    pow_file   = '_eeg_htpCalcRestPower_band.csv'; %fullfile(analysis_outputdir, [basename '_eeg_htpCalcRestPower_band.csv']); stores power band result
    spectro_file   = '_eeg_htpCalcRestPower_spectro.csv';%fullfile(analysis_outputdir, [basename '_eeg_htpCalcRestPower_spectro.csv']); store spectogram results
    qi_file   = '_eeg_htpCalcRestPower_qi.csv';%fullfile(analysis_outputdir, [basename '_eeg_htpCalcRestPower_qi.csv']);

    % START: Signal Processing

    % Key Parameters
    t = ip.Results.duration; % time in seconds
    fs = EEG.srate; % sampling rate
    win = ceil(ip.Results.window * fs); % window
    nfft = win; % FFT points--
    noverlap = .5 * win; % points overlap
    channo = EEG.nbchan;% number of EEG channels
    EEG.subject = EEG.setname;%stores the setnames in the subject identifier for EEG

    labels = bandDefs(:, 1);
    freq = cell2mat(bandDefs(:, 2:3));%extract requency band labels and their corresponding frequency ranges
    % dataset validation
    % is size sufficient for duration and offset?
    samples = t * fs; % if using time, number of samples
    start_sample = ip.Results.offset * fs; if start_sample == 0, start_sample = 1; end
    total_samples = EEG.pnts * EEG.trials;%calculates the number of samples to analyze and the starting sample
    
    if samples >= total_samples - start_sample %checks if the analysis duration go over the available data
        samples = total_samples;%if true then will adjust the analysis to use all available data and issues a warning
        start_samples = 1; % in samples
        warning("Insufficient Data, using max samples.")
    end

    % calculate power from first and last frequency from banddefs
    if ndims(EEG.data) > 2 %#ok<ISMAT> : Checks whether the EEG data has more than 2 dimensions
        %dat = permute(detrend3(permute(EEG.data, [2 3 1])), [3 1 2]);
        dat = EEG.data; %if it does the it reshapes the data into 2d data
        cdat = reshape(dat, size(dat, 1), size(dat, 2) * size(dat, 3));
    else
        cdat = EEG.data;
    end
    fprintf("Finished processing")
    % define final input data
    cdat = cdat(:, start_sample:end);
%all of this is to prepare for doing welch statistics
    % switch on gpu
    if ip.Results.gpuOn, cdat = gpuArray(cdat); end

    % power computation; pxx contain the power spectral density and f
    % contains the frequency
    [pxx, f] = pwelch(cdat', hanning(win), noverlap, freq(1, 1):.5:freq(end, 2), fs); %#ok<*ASGLU> ;computes power spectral density using Welch's method
    if size(pxx,1)==EEG.nbchan, pxx=pxx'; f=f'; end % added for single channel data (like mouse electrodes); ensures the PSD output is the correct orientation
    if ip.Results.gpuOn, pxx = gather(pxx); f = gather(f); end

    % power derivations
    pow_abs = pxx(1:end, :); % absolute power (V^2/Hz); copies all rows and collumns of pxx in pow_abs
    pow_db = 10 * log10(pow_abs); % absolute power dB/Hz; converts pwoer to decibels

    pow_rel = NaN * ones(size(pow_abs)); % relative power (unitless); creates a matrix for relative power filled with NaNs

    for chani = 1:size(pow_rel, 2)%loop uesed to calculate relative power for each channel
        pow_rel(:, chani) = pow_abs(:, chani) ./ sum(pow_abs(:, chani));%dividing power at each frequency by the total pwer across all frequencies for a channel
    end

    % band averaged power
    pow_prealloc = zeros(length(freq), channo);%just preallocating matrices
    pow_abs_band = pow_prealloc; pow_db_band = pow_prealloc; pow_rel_band = pow_prealloc;

    for bandi = 1:length(freq)%loop used to calculate the average pwoer within each frequency band
        current_band = freq(bandi, :);
        freqidx = [find(f == current_band(1)):find(f == current_band(2))];
        pow_abs_band(bandi, :) = squeeze(mean(pow_abs(freqidx, :), 1));
        pow_db_band(bandi, :) = squeeze(mean(pow_db(freqidx, :), 1));
        pow_rel_band(bandi, :) = squeeze(mean(pow_rel(freqidx, :), 1));
    end

    % create output table
    pow_abs_band = pow_abs_band';%unsure about the single quote
    pow_db_band = pow_db_band';
    pow_rel_band = pow_rel_band';
    % creating labels for each type of power measurement 
    abs_labels = cellfun(@(x) sprintf('abs_%s', x), labels', 'uni', 0);
    db_labels = cellfun(@(x) sprintf('db_%s', x), labels', 'uni', 0);
    rel_labels = cellfun(@(x) sprintf('rel_%s', x), labels', 'uni', 0);
    %combines all power measurements into a single table
    allbandlabels = [abs_labels db_labels rel_labels];
    powertable = array2table([pow_abs_band pow_db_band pow_rel_band], 'VariableNames', allbandlabels);
    % creates a table filled with the EEg metadata
    infocolumns = table(string(repmat(EEG.subject, channo, 1)), string(repmat(EEG.filename, channo, 1)), ...
        {EEG.chanlocs.labels}', 'VariableNames', {'eegid', 'filename', 'chan'});
    % combines the EEG and power data into one table 
    csvtable = [infocolumns, powertable];

    chan_names = {EEG.chanlocs.labels};
    slabel = @( measure_prefix ) cellfun(@(chan_name_cell) sprintf('%s_%s', ... % function slabel is used to generate channel specific names
        measure_prefix, chan_name_cell), chan_names, 'uni',0);

    abs_labels = slabel('abspow');% uses slabel function to generate labels for each type power measurment
    db_labels = slabel('dbpow');
    rel_labels = slabel('relpow');
    spectro_labels = [{'freq'} {abs_labels{:}} ...
         {db_labels{:}} {rel_labels{:}}];

    spectro_values = array2table([f pow_abs pow_db pow_rel ], ... % creates a table called spectro_values  that contain each of the power measurements for each channel
         'VariableNames', spectro_labels);
    spectro_info = table(repmat(EEG.subject, length(f), 1), ...
         repmat('spectrogram', length(f), 1), 'VariableNames', {'eegid', 'measure'});%creates a table of metadata for each row of the spectral data
% both spectro_values and spectro_info are designed to  work together 
    % deprecated: mean spectrogram across channels
%     spectro_values = array2table([f ...
%                                 mean(pow_abs, 2) mean(pow_db, 2) mean(pow_rel, 2)], ...
%         'VariableNames', {'freq', 'abspow', 'dbpow', 'relpow'});
% 
%     spectro_info = table(repmat(EEG.subject, length(f), 1), ...
%         repmat('mean', length(f), 1), 'VariableNames', {'eegid', 'chan'});

    % END: Signal Processing

    % QI Table
    qi_table = cell2table({EEG.setname, EEG.filename, functionstamp, timestamp, ... %creates another table containing, EEG dataset name, EEG filename, name of current function, timestamp, total EEG trials, number of time points, sampling rate, max and min time values
        EEG.trials, EEG.pnts, EEG.srate, EEG.xmin, EEG.xmax}, ...
    'VariableNames', {'eegid', 'filename', 'scriptname', 'timestamp', ...
    'trials', 'points', 'srate', 'xmin', 'xmax'});

    % Outputs:
    EEG.vhtp.eeg_htpCalcRestPower.summary_table = csvtable; %stores the table created earlier into the EEG structure
    EEG.vhtp.eeg_htpCalcRestPower.pow.spectro = [spectro_info, spectro_values]; %combines spectro_info and spectro_values into the EEG structures
    EEG.vhtp.eeg_htpCalcRestPower.qi_table = qi_table; %Stores qi_table into the EEG structure
    results = EEG.vhtp.eeg_htpCalcRestPower;% assigns all the sotred results to a variable called 'results'

    % file management
    if ~isempty(ip.Results.outputdir) %checks if output directory has been specified
        if ~ip.Results.useParquet% if parquet format is not used will do the following
            writetable(results.summary_table, pow_file);%writes summary to pow_file
            writetable(results.pow.spectro, spectro_file);% writes spectogram data to spectro file
            writetable(results.qi_table, qi_file);% writes the qi_table into qi_file

        else% if in parquet format it will do the following
            parquetwrite(strrep(pow_file,'.csv','.parquet'),results.summary_table);%writes the summary table and spectrogram data in  parquet format
            parquetwrite(strrep(spectro_file,'.csv','.parquet'),results.pow.spectro);
            writetable(results.qi_table, qi_file);%still the same either way

        end
        note(sprintf('%s saved in %s.\n', EEG.setname, ip.Results.outputdir))%used to provide information about where the results are saved
    else
        note('Specify output directory to save results.\n')
    end

    function note = htp_utilities()%defines htp_utilities function, which prints messages with the current function name as a prefix
            note        = @(msg) fprintf('%s: %s\n', mfilename, msg );
    end
fprintf("Finished")
end  