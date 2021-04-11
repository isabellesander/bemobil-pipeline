% bemobil_preprocess() - Preprocessing of EEG data: Fill EEG structure with
% ur-data. Remove unused electrodes of electrode arrays. Import channel
% locations from vicra file. Change channel names and declare EOG.
% Resample.
%
% Usage:
%   >>  [ ALLEEG EEG CURRENTSET ] = bemobil_preprocess(ALLEEG, EEG, CURRENTSET, channel_locations_filepath, channels_to_remove, eog_channels, resample_freq);
%   >>  [ ALLEEG EEG CURRENTSET ] = bemobil_preprocess(ALLEEG, EEG, CURRENTSET, channel_locations_filepath, channels_to_remove, eog_channels, resample_freq, out_filename, out_filepath);
%
% Inputs:
%   ALLEEG                  - complete EEGLAB data set structure
%   EEG                     - current EEGLAB EEG structure
%   CURRENTSET              - index of current EEGLAB EEG structure within ALLEEG
%   channel_locations_file  - channel_locations file (with path); OR []
%   channels_to_remove      - cell of all channels that should be thrown out
%       per se (e.g. {'N29' 'N30' 'N31'}); OR []
%   eog_channels            - cell of channels that should be declared as EOG
%       for later use (e.g. {'G16' 'G32'}); OR []
%   resample_freq           - Resample frequency (Hz), if [], no resampling will be applied
%   out_filename            - output filename (OPTIONAL ARGUMENT)
%   out_filepath            - output filepath (OPTIONAL ARGUMENT - File will only be saved on disk
%       if both a name and a path are provided)
%
% Outputs:
%   ALLEEG                  - complete EEGLAB data set structure
%   EEG                     - current EEGLAB EEG structure
%   Currentset              - index of current EEGLAB EEG structure within ALLEEG
%
%   .set data file of current EEGLAB EEG structure stored on disk (OPTIONALLY)
%
% See also:
%   EEGLAB, pop_eegfiltnew, pop_resample, pop_chanedit, pop_select
%
% Authors: Lukas Gehrke, 2017, Marius Klug, 2021

function [ ALLEEG EEG CURRENTSET ] = bemobil_preprocess(ALLEEG, EEG, CURRENTSET,...
    channel_locations_filepath, channels_to_remove, eog_channels, resample_freq,...
    out_filename, out_filepath, rename_channels, ref_channel, zaplineConfig)

% only save a file on disk if both a name and a path are provided
save_file_on_disk = (exist('out_filename', 'var') && ~isempty(out_filename) &&...
    exist('out_filepath', 'var') && ~isempty(out_filepath));

% check if file already exist and show warning if it does
if save_file_on_disk
    mkdir(out_filepath); % make sure that folder exists, nothing happens if so
    dir_files = dir(out_filepath);
    if ismember(out_filename, {dir_files.name})
        warning([out_filename ' file already exists in: ' out_filepath '. File will be overwritten...']);
    end
end

% fill/copy all ur_structures with raw data (e.g. copy event to urevent)
EEG = eeg_checkset(EEG, 'makeur');

% remove unused neck electrodes from file (if BeMoBIL layout is used as is)
if ~isempty(channels_to_remove)
    if all(ismember(channels_to_remove, {EEG.chanlocs.labels}))
        % b) remove not needed channels import "corrected" chanlocs file
        EEG = pop_select( EEG,'nochannel', channels_to_remove);
        EEG = eeg_checkset( EEG );
        disp(['Removed electrodes: ' channels_to_remove ' from the dataset.']);
    else
        error('Not all of the specified channels to remove were present als data channels!')
    end
else
    disp('No channels to remove specified, skipping this step.')
end

% Resample if frequency is provided
if ~isempty(resample_freq)
    if resample_freq ~= EEG.srate
        warndlg(['Resampling in BIDS processing did not work properly. Resampling to ' num2str(resample_freq), 'Hz.'])
        EEG = pop_resample(EEG, resample_freq);
        EEG = eeg_checkset( EEG );
    end
end

%% Clean line noise with ZapLine: de Cheveigne, A. (2020) ZapLine: a simple and effective method to remove power line
% artifacts. Neuroimage, 1, 1-13.
if exist('zaplineConfig','var')
    
    [EEG, plothandles] = bemobil_clean_data_with_zapline(EEG, zaplineConfig);
        
    if save_file_on_disk 
        disp('Saving ZapLine figures...')

        filenamesplit = strsplit(out_filename,'.set');

        for i_fig = 1:length(plothandles)
        
            savefig(plothandles(i_fig),fullfile(out_filepath,[filenamesplit{1}...
                '_' matlab.lang.makeValidName(['zapline_' num2str(zaplineConfig.linefreqs(i_fig))]) '.fig']))
            saveas(plothandles(i_fig),fullfile(out_filepath,[filenamesplit{1}...
                '_' matlab.lang.makeValidName(['zapline_' num2str(zaplineConfig.linefreqs(i_fig))]) '.png']))
            close(plothandles(i_fig))
            
        end

        disp('...done')
    end
    
end

%%

% rename channels if specified
if exist('rename_channels','var') && ~isempty(rename_channels)
    disp('Renaming channels...')
    
    if length(rename_channels) == 1
        
        for i_chan = 1:EEG.nbchan
            EEG.chanlocs(i_chan).labels = erase(EEG.chanlocs(i_chan).labels,rename_channels{1});
        end
        
    else
        
        for i_pair = 1:size(rename_channels,1)

            old_chanidx = find(strcmp({EEG.chanlocs.labels},rename_channels{i_pair,1}));

            if ~isempty(old_chanidx)
                EEG=pop_chanedit(EEG, 'changefield',{old_chanidx 'labels' rename_channels{i_pair,2}});
            else
                warning(['Did not find channel ' rename_channels{i_pair,1} '. Skipping...'])
            end

        end
    end
end

% 1b3) add ref channel as zero if specified
if exist('ref_channel','var') && ~isempty(ref_channel)
    disp('Adding reference channel with zeros...')
    
    EEG.nbchan = EEG.nbchan + 1;
    EEG.data(end+1,:) = zeros(1, EEG.pnts);
    EEG.chanlocs(end+1).labels = ref_channel;
    
    EEG = eeg_checkset(EEG);
    
    disp('...done.')
    
end

% 1c) import chanlocs and copy to urchanlocs
if ~isempty(channel_locations_filepath)
    EEG = pop_chanedit(EEG, 'load',...
        {channel_locations_filepath 'filetype' 'autodetect'});
    disp('Imported channel locations.');
    EEG.urchanlocs = EEG.chanlocs;
else
    standard_channel_locations_path =...
        fullfile(fileparts(which('dipfitdefs')),'standard_BESA','standard-10-5-cap385.elp');
    
    EEG = pop_chanedit(EEG,'lookup',standard_channel_locations_path);
end

% this has to happen after loading chanlocs because chanlocs are being completely overwritten in the process
if exist('ref_channel','var')
    
    disp('Declaring ref for all channels...')
    
    [EEG.chanlocs(:).ref] = deal(ref_channel);
    
end


% 1d) change channel types in standard MoBI montage declaring the EOG and EEG channels
for n = 1:length(EEG.chanlocs)
    if ismember(lower(EEG.chanlocs(n).labels), lower(eog_channels))
        EEG.chanlocs(n).type = strcat('EOG');
        disp(['Added channel type: ', EEG.chanlocs(n).labels, ' is EOG electrode!!']);
    elseif exist('ref_channel','var') && strcmpi(EEG.chanlocs(n).labels, ref_channel)
        EEG.chanlocs(n).type = strcat('REF');
        disp(['Added channel type: ', EEG.chanlocs(n).labels, ' is REF electrode!!']);
    else
        EEG.chanlocs(n).type = strcat('EEG');
        disp(['Added channel type: ', EEG.chanlocs(n).labels, ' is EEG electrode.']);
    end
end

EEG = eeg_checkset( EEG );


%% new data set in EEGLAB
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'gui', 'off');
EEG = eeg_checkset( EEG );

% save on disk
if save_file_on_disk
    EEG = pop_saveset( EEG, 'filename',out_filename,'filepath', out_filepath);
    disp('...done')
end

[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);