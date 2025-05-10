clear;
close all;
%% define the subject and start diary
dtime = string(datetime);
diary_name = string(strcat('merge_blocks_',dtime,'.txt'));
diary_name = strrep(diary_name, ' ', '-');
diary_name = strrep(diary_name, ':', '-');
diary(diary_name)
disp(dtime)
addpath("eeglab2024.2\")
eeglab;
pre_processing_params_path = "pre_processing_parameters_final.xlsx"
pre_processing_params =  readtable(pre_processing_params_path)
params_all = pre_processing_params(strcmp(pre_processing_params.site,"All"),:)
epoch_max = params_all.epoch_max
baseline_min = params_all.baseline_min
epoch_window = [baseline_min epoch_max] %window for epoching in s
datapath_base = 'D:\REFTEP_ALL\EEG_preprocessing_data\'
datapath_raw = "D:\REFTEP_ALL\Data_raw\"
sites = {'Aalto', 'Tuebingen'};
for site=sites
     %go through all the sites and define specific paths
    site_char = char(site);
    directory_name_site = fullfile(datapath_base,strcat('Preprocessing_',site_char,"\"));
    rawsite = char(fullfile(datapath_raw,strcat(site_char,"_raw","\")))
    files_and_folders = dir(rawsite)
    is_subfolder = [files_and_folders.isdir]
    folders = files_and_folders(is_subfolder)
    names = {folders.name}
    subject_names = names(contains(names,"sub"))
    %trigger names for the site
    params_site = pre_processing_params(strcmp(pre_processing_params.site,site_char),:);
    trigger_names = eval(params_site.trigger_names{1})
    for index = 1:length(subject_names)
        subject_name = subject_names{index}
        subject_datapath_raw = char(fullfile(datapath_raw,strcat(site_char,"_raw","\"),strcat(subject_name,"\")))
        eeg_preprocessing_path_subject = char(fullfile(directory_name_site,subject_name))
        file_list = dir(subject_datapath_raw);
        fnames = {file_list.name}
        tep_files_set = file_list(contains(fnames,'task-tep') & contains(fnames,'.set'))
        setind = 1;
        for k = 1:length(tep_files_set)
            %filename = tep_files_set(k).name
            filename = char(strcat(subject_name,'_task-tep_run-0',string(k),'_eeg.set'))
            EEG = pop_loadset(filename, subject_datapath_raw);
            if ~isempty(EEG.event)
                EEG.event = EEG.event(ismember({EEG.event.type},trigger_names)); %select only trigger events
                epochs_name = strcat('epochs_block_',string(setind));
                EEG = pop_epoch(EEG, trigger_names, epoch_window, 'newname', epochs_name, 'epochinfo', 'yes')
                if EEG.trials > 0 & ~isempty(EEG.event)
                    if setind == 1
                        EEG_merged = EEG;
                    else
                        EEG_merged = pop_mergeset(EEG_merged,EEG);
                    end
                    setind = setind + 1;
                end
            end
        end
        mkdir(eeg_preprocessing_path_subject);
        EEG_merged = pop_editset(EEG_merged, 'setname', strcat(subject_name,'-teps_epochs_merged'));
        filetosave_merged = strcat(subject_name,'_task-tep_epochs_merged_eeg.set')
        pop_saveset(EEG_merged, 'filename',filetosave_merged,'filepath',eeg_preprocessing_path_subject);
    end
end