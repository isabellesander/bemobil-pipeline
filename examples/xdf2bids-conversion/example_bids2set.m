addpath(genpath("/Users/izzys/Documents/matlab_analyses/bids-example-specification-main/"));
%addpath(genpath("/Users/lukasgehrke/code/Izzy_conversion_BIDS/bids-example-specification-main"));

env = readstruct('bemobil.json');
subject = readstruct('subject.json');
eeg = readstruct('eeg.json');
% phys = readstruct('physio.json');
    
subject.bids_target_folder = strcat(general.TaskName{1}, env.bids_target_folder{1}, filesep);
subject.set_folder         = fullfile(subject.bids_target_folder, 'derivatives', env.raw_eeglab_folder{1});

subject.session_names = {'a', 'b'}; % fixed this here again for string compatibility

subject.other_data_types   = {env.other_data_types{1}};
subject.match_electrodes_channels = {eeg.eeg_chanloc_names};

for s = [20]
    subject.subject            = s;
    bemobil_bids2set(subject);
end
