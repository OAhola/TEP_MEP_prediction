% run ICA for REFTEP subjects (manual selection of ocular artifacts)
clear;
close all;
dtime = string(datetime);
diary_name = string(strcat('preprocessing_diary_ICA_',dtime,'.txt'));
diary_name = strrep(diary_name, ' ', '-');
diary_name = strrep(diary_name, ':', '-');
diary(diary_name)
disp(dtime)
addpath("eeglab2024.2\")
eeglab;
pre_processing_params_path = "pre_processing_parameters_final.xlsx"
pre_processing_params =  readtable(pre_processing_params_path)
params_all = pre_processing_params(strcmp(pre_processing_params.site,"All"),:)
timtopo_timerange = eval(params_all.timtopo_timerange{1})
timtopo_topotimes = eval(params_all.timtopo_topotimes{1})
baseline_min = params_all.baseline_min
baseline_max = params_all.baseline_max
baseline_window = [baseline_min*1000 baseline_max*1000] %baseline window in ms
datapath_base = 'D:\REFTEP_ALL\EEG_preprocessing_data\'
%% define the subject and start diary
sites = {'Aalto'};
for site=sites
    %go through all the sites and define specific paths
    site_char = char(site);
    directory_name_site = fullfile(datapath_base,strcat('Preprocessing_',site_char,"\"));
    files_and_folders = dir(directory_name_site)
    is_subfolder = [files_and_folders.isdir]
    folders = files_and_folders(is_subfolder)
    names = {folders.name}
    subject_names = names(contains(names,"sub"))
    for index = 1:length(subject_names)
        reftep_subject = char(subject_names(index))
        directory_path = char(fullfile(directory_name_site,reftep_subject,"\"));
        filename_eeg = append(reftep_subject,'_EEG_with_ICA.set');
        EEG = pop_loadset(filename_eeg,directory_path); %load the data
        %% remove only ocular artifacts -- blinks and eye movement
        %EEG = pop_tesa_pcacompress(EEG, 'compVal', n_ica_components, 'plot','off' );
        %EEG = pop_tesa_fastica(EEG, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'off' );
        EEG = pop_tesa_compplot(EEG,'figSize','large','plotTimeX',[-200 300],'plotFreqX',[1 100], 'freqScale','log', 'saveWeights','on' );
        %% baseline correct and save the data after ICA, you could also look at the data'
        EEG = pop_rmbase(EEG, baseline_window ,[]);
        figname_ICA_pruned = strcat(reftep_subject,'_ICA_pruned_TEP_data');
        figure; pop_timtopo(EEG, timtopo_timerange, timtopo_topotimes, figname_ICA_pruned);
        saveas(gcf,append(directory_path,strcat(figname_ICA_pruned, '.png')));
        close(gcf);
        afterica_baseline_name = append(reftep_subject,'_EEG_after_ICA_baseline_corrected.set');
        pop_saveset(EEG, 'filename',afterica_baseline_name,'filepath',directory_path);
    end
end