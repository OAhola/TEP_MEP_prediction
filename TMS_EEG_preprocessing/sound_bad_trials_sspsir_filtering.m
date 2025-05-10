% remove bad trials, run SOUND and SSP-SIR and filter the data
clear;
close all;
dtime = string(datetime);
diary_name = string(strcat('preprocessing_diary_sound_trials_sspsir_filtering_',dtime,'.txt'));
diary_name = strrep(diary_name, ' ', '-');
diary_name = strrep(diary_name, ':', '-');
diary(diary_name)
disp(dtime)
addpath("eeglab2024.2\");
eeglab;
pre_processing_params_path = "pre_processing_parameters_final.xlsx"
pre_processing_params =  readtable(pre_processing_params_path)
params_all = pre_processing_params(strcmp(pre_processing_params.site,"All"),:)
timtopo_timerange = eval(params_all.timtopo_timerange{1})
timtopo_topotimes = eval(params_all.timtopo_topotimes{1})
baseline_min = params_all.baseline_min
baseline_max = params_all.baseline_max
baseline_window = [baseline_min*1000 baseline_max*1000] %baseline window in ms
baseline_window_in_seconds = [baseline_min baseline_max]
detection_window_r1 = eval(params_all.detection_window_r1{1})
detection_window_r2 = eval(params_all.detection_window_r2{1})
detection_window_r3 = [detection_window_r1(2) detection_window_r2(2)]
epoch_max = params_all.epoch_max
epoch_window = [baseline_min epoch_max] %window for epoching in s
sound_lambda = params_all.sound_lambda
sound_n_iters = params_all.sound_n_iters
sspsir_percentile = params_all.ssp_sir_percentile
pulse_artifact_window2 = eval(params_all.pulse_artifact_window2{1})
interpolation_window2 = eval(params_all.interpolation_window2{1})
resample_to = params_all.resample_to
powerline_range1 = eval(params_all.powerline_range1{1})
filter_range = eval(params_all.filter_range{1})
filter_order = params_all.filter_order
final_time_window_eeg = eval(params_all.final_time_window_eeg{1})
rejtrial_threshold_global_pre = params_all.rejtrial_threshold_global_pre
rejtrial_thresholds_local_pre = eval(params_all.rejtrial_thresholds_local_pre{1})
nums_over_rejtrial_thresholds_local_pre = eval(params_all.nums_over_rejtrial_thresholds_local_pre{1})
rejtrial_threshold_global_post = params_all.rejtrial_threshold_global_post
rejtrial_thresholds_local_post = eval(params_all.rejtrial_thresholds_local_post{1})
nums_over_rejtrial_thresholds_local_post = eval(params_all.nums_over_rejtrial_thresholds_local_post{1})
time_window_to_detect_bad_trials_pre = baseline_window_in_seconds
window_lims = linspace(time_window_to_detect_bad_trials_pre(1),time_window_to_detect_bad_trials_pre(2),4);
window_lims2 = [];
windows = {};
for i=1:length(window_lims)-1
    window = [window_lims(i) window_lims(i+1)];
    windows{end+1} = window;
    window_lims2(end+1) = mean(window);
end
windows{end+1} = [window_lims(1) window_lims2(1)];
for i=1:length(window_lims2)-1
    windows{end+1} = [window_lims2(i) window_lims2(i+1)];
end
windows{end+1} = [window_lims2(end) window_lims(end)];
windows
datapath_base = 'D:\REFTEP_ALL\EEG_preprocessing_data\'
sites = {'Aalto','Tuebingen'};
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
        reftep_subject = char(subject_names(index));
        directory_name = char(fullfile(directory_name_site,reftep_subject,"\")); %the directory of the EEG data to be loaded
        filename_eeg = char(strcat(reftep_subject,'_EEG_after_ICA_baseline_corrected.set')); %.set file name of the eeg data with good channels
        filename_eeg_allchans = char(strcat(reftep_subject,'_EEG_with_all_chans_epoched.set'));
        EEG = pop_loadset(filename_eeg, directory_name); %load the data
        %% apply SOUND
        %pop_eegplot(EEG,1,1,1)
        EEG = pop_reref(EEG,[]);
        EEG = pop_tesa_sound(EEG, 'lambdaValue', sound_lambda, 'iter', sound_n_iters, 'leadfieldInFile' ,[], 'leadfieldChansFile', [], 'replaceChans', char(fullfile(directory_name,filename_eeg_allchans)), 'multipleConds', []);
        EEG = pop_rmbase(EEG, baseline_window ,[]);
        %plot average TEP after SOUND
        figname_after_SOUND = strcat(reftep_subject,'_after_SOUND_TEP_data');
        figure; pop_timtopo(EEG, timtopo_timerange, timtopo_topotimes, figname_after_SOUND);
        saveas(gcf,append(directory_name,strcat(figname_after_SOUND, '.png')));
        close(gcf);
        %save data after SOUND
        filename_eeg_after_sound = char(strcat(reftep_subject,'_EEG_after_SOUND.set'));
        EEG.setname = 'epochs_after_sound';
        %pop_saveset(EEG,'filename',filename_eeg_after_sound,'filepath',directory_name)
        %% remove bad trials by detecting outlier epoched data from different windows
        %select windows of EEG data
        %pop_eegplot(EEG,1,1,1)
        epoch_inds = 1:size(EEG.data,3);
        EEG_baseline = pop_select(EEG, 'time', time_window_to_detect_bad_trials_pre);
        bad_trials = []; %init array
        bad_trials_baseline = find_bad_trials(EEG_baseline, rejtrial_threshold_global_pre,rejtrial_thresholds_local_pre,nums_over_rejtrial_thresholds_local_pre);
        bad_trials = union(bad_trials, bad_trials_baseline);
        formatted_array = strjoin(string(time_window_to_detect_bad_trials_pre),', ');
        fprintf('n bad trials from baseline window [%s]: %d\n',formatted_array, length(bad_trials_baseline));

        for i=1:length(windows)%go through pre-stim detection windows
            window = windows{i};
            EEG_window = pop_select(EEG, 'time', window);
            bad_trials_in_window = find_bad_trials(EEG_window, rejtrial_threshold_global_pre,rejtrial_thresholds_local_pre,nums_over_rejtrial_thresholds_local_pre);
            formatted_array = strjoin(string(window),', ');
            fprintf('n bad trials from window [%s]: %d\n',formatted_array, length(bad_trials_in_window));
            bad_trials = union(bad_trials, bad_trials_in_window);
        end
        detection_windows_post = {detection_window_r1,detection_window_r2,detection_window_r3};
        for i=1:length(detection_windows_post) %go through post-stim detection windows
            detection_window = detection_windows_post{i};
            EEG_window_post = pop_select(EEG, 'time', detection_window);
            bad_trials_in_window_post = find_bad_trials(EEG_window_post, rejtrial_threshold_global_post, rejtrial_thresholds_local_post, nums_over_rejtrial_thresholds_local_post);
            formatted_array = strjoin(string(detection_window),', ');
            fprintf('n bad trials from window [%s]: %d\n',formatted_array, length(bad_trials_in_window_post));
            bad_trials = union(bad_trials, bad_trials_in_window_post);
        end
        disp(reftep_subject)
        fprintf("trials removed: %.f \n", length(bad_trials))
        bad_trials_binary = ismember(epoch_inds,bad_trials);
        EEG = pop_rejepoch(EEG, bad_trials_binary,0); %reject the conbined bad trials
        EEG.setname = 'epochs_after_removing_bad_trials';
        %plot average TEP after removing bad trials
        figname_after_removing_bad_trials = strcat(reftep_subject,'_bad_trials_removed_TEP_data');
        figure; pop_timtopo(EEG, timtopo_timerange, timtopo_topotimes, figname_after_removing_bad_trials);
        saveas(gcf,append(directory_name,strcat(figname_after_removing_bad_trials, '.png')));
        close(gcf);
        filename_eeg_good_trials = char(strcat(reftep_subject,'_EEG_good_trials.set'));
        %pop_saveset(EEG,'filename',filename_eeg_good_trials,'filepath',directory_name)
    
        %% suppress muscle-artifacts with SSP-SIR
        EEG = pop_tesa_sspsir(EEG, 'artScale','automatic','PC', {'data',[sspsir_percentile]});
        figname_after_sspsir = strcat(reftep_subject,'_after_sspsir_',num2str(sspsir_percentile),'_TEP_data');
        figure; pop_timtopo(EEG, timtopo_timerange, timtopo_topotimes, figname_after_sspsir);
        saveas(gcf,append(directory_name,strcat(figname_after_sspsir, '.png')));
        close(gcf);
        
        %save the data after SSP-SIR
        filename_eeg_after_sspsir = char(strcat(reftep_subject,'_EEG_after_SSPSIR.set'));
        EEG.setname = 'epochs_after_sspsir';
        %pop_saveset(EEG,'filename',filename_eeg_after_sspsir,'filepath',directory_name)
        %%
        %Interpolate the artifact-recovered pulse artifact window
        EEG = pop_tesa_removedata(EEG, pulse_artifact_window2);
        EEG = pop_tesa_interpdata(EEG, 'cubic', interpolation_window2);
        %%
        %filter to frequencies of interest and filter out powerline noise
        EEG = pop_tesa_filtbutter(EEG, 1, filter_range(2), filter_order, 'bandpass');
        EEG = pop_tesa_filtbutter(EEG, powerline_range1(1), powerline_range1(2), filter_order, 'bandstop');
        
        %get rid of edge artifacts by cropping the epoch window
        EEG = pop_select(EEG, 'time', final_time_window_eeg);
    
        %downsample the data
        EEG = pop_resample(EEG, resample_to);
        baseline_window_new = [final_time_window_eeg(1)*1000 baseline_max*1000];
        EEG = pop_rmbase(EEG, baseline_window_new ,[]); %baseline-correct
        EEG = pop_reref(EEG, []); %finally set to average reference
        EEG.setname = 'epochs_cleaned_eeg_baseline_corrected_average_ref';
        figname_cleaned_data = strcat(reftep_subject,'_cleaned_TEP_data');
        figure; pop_timtopo(EEG, timtopo_timerange, timtopo_topotimes, figname_cleaned_data);
        saveas(gcf,append(directory_name,strcat(figname_cleaned_data, '.png')));
        close(gcf);
        %% save the final pre-processed data (bad emg trials have not yet been
        %removed here)
        eeg_cleaned_name = append(reftep_subject,'_EEG_cleaned.set');
        pop_saveset(EEG, 'filename',eeg_cleaned_name,'filepath',directory_name);
    end
end