%addpath(genpath("/Users/izzys/Documents/matlab_analyses/bids-example-specification-main/"));
addpath(genpath("/Users/izzys/Documents/dataanalysis/bids-example-specification-main/"));

%% general metadata, recording environment and hardware specs
general = readstruct('general.json');
eeg = readstruct('eeg.json');
phys = readstruct('physio.json');
env = readstruct('bemobil.json');

subj = [20]
%% info about the streams of this particular subject as they appear in the XDF file to parse
subject = readstruct('subject.json');

%% parse the XDF file and write the BIDS files
subject.bids_target_folder    = strcat(general.TaskName{1}, env.bids_target_folder{1});
subject.task                  = general.TaskName{1};

subject.eeg.stream_name       = eeg.eeg_stream_name{1};
subject.eeg.channel_labels    = cellstr(eeg.eeg_chanloc_names);
subject.eeg.ref_channel       = eeg.EEGReference;
subject.eeg.chanloc           = [];

subject.phys.streams{1}.stream_name = 'EyeLogic';
phys.physio = phys;

subject.bids_parsemarkers_custom = '';
for s = subj
    tmp_subject_naming = strcat('vp', sprintf('%02d', s)); % fixed the loop variable and sprintf
    
    for session = subject.session_names
        
        subject.subject               = s;                     % setting subject.subject to current loop variable s
        subject.session               = session{1};              % optional
        subject.filename              = fullfile(env.bids_source_folder, tmp_subject_naming, strcat(tmp_subject_naming, '_', session, '.xdf'));
    
        % try
            if exist(subject.filename, 'file') == 2
                bemobil_xdf2bids(subject, ...
                    'general_metadata', general,...
                    'eeg_metadata', eeg, ...
                    'physio_metadata', phys ...
                );
            else
                fprintf('File %s does not exist, skipping...\n', subject.filename);
            end
        % catch ME
        %     fprintf('Error processing file %s: %s\n', subject.filename, ME.message);
        % end
    end
end
    
subject.bids_target_folder = strcat(general.TaskName{1}, env.bids_target_folder{1}, filesep);
subject.set_folder         = fullfile(subject.bids_target_folder, 'derivatives', env.raw_eeglab_folder{1});

subject.session_names = {'a', 'b'}; % fixed this here again for string compatibility

subject.other_data_types   = {env.other_data_types{1}};
subject.match_electrodes_channels = {eeg.eeg_chanloc_names};

for s = subj
    subject.subject            = s;
    bemobil_bids2set(subject);
end

