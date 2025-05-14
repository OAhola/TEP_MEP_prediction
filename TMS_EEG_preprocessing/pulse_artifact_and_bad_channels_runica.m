clear;
close all;
%% define the subject and start diary
dtime = string(datetime);
diary_name = string(strcat('preprocessing_diary_channels_',dtime,'.txt'));
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
interpolation_window = eval(params_all.interpolation_window1{1})
baseline_min = params_all.baseline_min
baseline_max = params_all.baseline_max
baseline_window = [baseline_min*1000 baseline_max*1000] %baseline window in ms
baseline_window_in_seconds = [baseline_min baseline_max]
epoch_max = params_all.epoch_max
epoch_window = [baseline_min epoch_max] %window for epoching in s
rejchan_threshold_pre = params_all.rejchan_threshold_pre
rejchan_threshold_post = params_all.rejchan_threshold_post
frontal_threshold_shift = eval(params_all.frontal_threshold_shift{1})
bad_trial_thresh = params_all.bad_trial_thresh
detection_window_r1 = eval(params_all.detection_window_r1{1})
detection_window_r2 = eval(params_all.detection_window_r2{1})
pulse_artifact_window = eval(params_all.pulse_artifact_window1{1})
pulse_artifact_window_plot = eval(params_all.pulse_artifact_window_plot{1})
good_trial_perc_scale = params_all.good_trial_perc_scale
bad_trial_perc_scale = params_all.bad_trial_perc_scale
pre_peak_threshes_large = eval(params_all.pre_peak_threshes_large{1})
pre_peak_threshes_small = eval(params_all.pre_peak_threshes_small{1})
n_ica_components = params_all.n_ica_components
sites = {'Tuebingen','Aalto'};
datapath_base = 'D:\REFTEP_ALL\EEG_preprocessing_data\'
for site=sites
    %go through all the sites and define specific paths
    site_char = char(site);
    directory_name_site = fullfile(datapath_base,strcat('Preprocessing_',site_char,"\"));
    files_and_folders = dir(directory_name_site)
    is_subfolder = [files_and_folders.isdir]
    folders = files_and_folders(is_subfolder)
    names = {folders.name}
    subject_names = names(contains(names,"sub"))
    %trigger names for the site
    params_site = pre_processing_params(strcmp(pre_processing_params.site,site_char),:);
    trigger_names = eval(params_site.trigger_names{1})
    for index = 1:length(subject_names)
        reftep_subject = char(subject_names(index))
        %% define file paths and directories for the subject
        directory_path = char(fullfile(directory_name_site,reftep_subject,"\")); %the directory of the EEG data to be loaded
        filename_tmseeg = char(strcat(reftep_subject,'_task-tep_epochs_merged_eeg.set')); %.set file name of the eeg data
        %% Load data, and divide EEG and EMG:
        % Divide the data into epochs
        EEG_and_EMG_epoched = pop_loadset(filename_tmseeg,directory_path);
        %EEG_and_EMG.event = EEG_and_EMG.event(ismember({EEG_and_EMG.event.type},trigger_names)); %select only trigger events
        %EEG_and_EMG_epoched = pop_epoch(EEG_and_EMG, trigger_names, epoch_window, 'newname', 'epochs', 'epochinfo', 'yes');
        n_epochs = size(EEG_and_EMG_epoched.data,3);
        epoch_inds = 1:n_epochs;
        length(epoch_inds)
        if strcmp(reftep_subject,'sub-105') %no noise masking in some trials
            bad_epochs_noisemasking = 280:285;
            epoch_inds = get_epoch_inds_good(EEG_and_EMG_epoched,bad_epochs_noisemasking);
        end
        if strcmp(reftep_subject,'sub-120') %no noise masking in some trials
            bad_epochs_noisemasking = 278:290;
            epoch_inds = get_epoch_inds_good(EEG_and_EMG_epoched,bad_epochs_noisemasking);
        end
        if strcmp(reftep_subject,'sub-110') %no noise masking in some trials
            bad_epochs_noisemasking = 580:590;
            epoch_inds = get_epoch_inds_good(EEG_and_EMG_epoched,bad_epochs_noisemasking);
        end
        if strcmp(reftep_subject,'sub-118') %no noise masking in this trial
            bad_epochs_noisemasking = 286:287;
            epoch_inds = get_epoch_inds_good(EEG_and_EMG_epoched,bad_epochs_noisemasking);
        end
        %{
        %optimally someting like this should exist here if get_epoch_inds_good was used for no noise masking trials so that the next step does not overwrite epoch_inds
        %however, we did not observe that the condition if length(uniq_events) ~= length([EEG_and_EMG_epoched.event.urevent]) | length(unique([EEG_and_EMG_epoched.event.epoch])) ~= length([EEG_and_EMG_epoched.event.epoch])
        %was fulfilled with the subjects that had no noise masking trials
        if bad_noises
            EEG_merged  = pop_select(EEG_and_EMG_epoched, 'trial',epoch_inds);
            n_epochs = size(EEG_and_EMG_epoched.data,3);
            epoch_inds = 1:n_epochs;
        end
        %}
        length(epoch_inds)
        uniq_events = unique([EEG_and_EMG_epoched.event.urevent]);
        if length(uniq_events) ~= length([EEG_and_EMG_epoched.event.urevent]) | length(unique([EEG_and_EMG_epoched.event.epoch])) ~= length([EEG_and_EMG_epoched.event.epoch])
            [unique_vals, ~, idx] = unique([EEG_and_EMG_epoched.event.urevent]);
            counts = histcounts(idx,numel(unique_vals));
            duplicates_urevents = unique_vals(counts>1);
            [unique_vals, ~, idx] = unique([EEG_and_EMG_epoched.event.epoch]);
            counts = histcounts(idx,numel(unique_vals));
            duplicates_epochs = unique_vals(counts>1);
            formatted_array1 = strjoin(string(duplicates_urevents),', ');
            formatted_array2 = strjoin(string(duplicates_epochs),', ');
            fprintf('dropping duplicate urevents/epochs! [%s] / [%s]\n',formatted_array1, formatted_array2)
            eventurevents = [];
            for t=1:size(EEG_and_EMG_epoched.epoch,2)
                eventurevents(end+1) = EEG_and_EMG_epoched.epoch(t).eventurevent{1};
            end
            contaminated_trials_binary = ismember(eventurevents,duplicates_urevents);
            fprintf('dropping %d epochs\n', sum(contaminated_trials_binary))
            EEG_and_EMG_epoched = pop_rejepoch(EEG_and_EMG_epoched, contaminated_trials_binary, 0);
            n_epochs = size(EEG_and_EMG_epoched.data,3);
            epoch_inds = 1:n_epochs;
        end
        n_epochs = length(epoch_inds);
        %select at maximum 1200 epochs
        if n_epochs > 1200 %restrict to maximum 1200 trials
            n_epochs = 1200;
            disp("more than 1200 epochs so restricting to 1200 epochs")
            epoch_inds = epoch_inds(epoch_inds <= n_epochs)
        end
        EEG_and_EMG_epoched  = pop_select(EEG_and_EMG_epoched, 'trial',epoch_inds);
        EEG_and_EMG_epoched = eeg_checkset(EEG_and_EMG_epoched);
        channel_names_original = {EEG_and_EMG_epoched.chanlocs.labels};
        channel_names_original_lower = lower(channel_names_original);
        emg_channels = channel_names_original(contains(channel_names_original_lower,'emg') | ...
            contains(channel_names_original_lower,'fdi') | contains(channel_names_original_lower,'apb')); %select emg channels
        emg_channels
        if length(emg_channels) ~=2
            disp("more or less than 2 emg channels detected")
        end
        EMG_epoched = pop_select(EEG_and_EMG_epoched, 'channel', emg_channels);
        for i=1:EMG_epoched.nbchan
            EMG_epoched.chanlocs(i).type = 'EMG'; %ecg also marked as emg for now
        end
        raw_emg_name = append(reftep_subject,'_EMG_raw_epoched.set');
        pop_saveset(EMG_epoched, 'filename',raw_emg_name,'filepath',directory_path);%select emg channels
        eeg_channels = channel_names_original(~contains(channel_names_original_lower,'emg') & ...
            ~contains(channel_names_original_lower,'fdi') & ~contains(channel_names_original_lower,'apb'));
        eeg_channels
        EEG_epoched = pop_select(EEG_and_EMG_epoched, 'channel', eeg_channels); %select eeg channels
        for i=1:EEG_epoched.nbchan
            EEG_epoched.chanlocs(i).type = 'EEG';
        end
        all_eeg_channels_name_after_dropping = append(reftep_subject,'_EEG_with_all_chans_epoched.set'); %save data with all channels
        pop_saveset(EEG_epoched, 'filename',all_eeg_channels_name_after_dropping,'filepath',directory_path);
        EEG = pop_rmbase(EEG_epoched, baseline_window ,[]); %baseline correct the data
        %% Plot data
        % Visualize the raw TEP data
        figname_non_processed_TEP_data = strcat(reftep_subject,'_non-processed_TEP_data');
        figure; pop_timtopo(EEG, timtopo_timerange, timtopo_topotimes, figname_non_processed_TEP_data);
        saveas(gcf,append(directory_path,strcat(figname_non_processed_TEP_data, '.png')));
        close(gcf);
        % Visualize the pulse artifact
        figname_pulse_artifact = strcat(reftep_subject,'_pulse_artifact_visualization');
        figure; pop_timtopo(EEG, pulse_artifact_window_plot, [NaN], figname_pulse_artifact);
        saveas(gcf,append(directory_path,strcat(figname_pulse_artifact, '.png')));
        close(gcf);
        %% Take care of pulse artifact
        % Remove Pulse artefact
        EEG = pop_tesa_removedata(EEG, pulse_artifact_window);
        % Interpolate the missing pulse data:
        EEG = pop_tesa_interpdata(EEG, 'cubic', interpolation_window);
        
        %visualize the raw TEP data and and pulse artifact window after
        %interpolation
        figname_pulse_artifact_interpolated_TEP_data = strcat(reftep_subject,'_pulse_artifact_interpolated_TEP_data');
        figure; pop_timtopo(EEG, timtopo_timerange, timtopo_topotimes, figname_pulse_artifact_interpolated_TEP_data);
        saveas(gcf,append(directory_path,strcat(figname_pulse_artifact_interpolated_TEP_data, '.png')));
        close(gcf);
        figname_pulse_artifact_interpolated = strcat(reftep_subject,'_pulse_artifact_interpolated_visualization');
        figure; pop_timtopo(EEG, pulse_artifact_window_plot, [NaN], figname_pulse_artifact_interpolated);
        saveas(gcf,append(directory_path,strcat(figname_pulse_artifact_interpolated, '.png')));
        close(gcf);
        %% Detect bad channels from the data and drop them
        EEG_baseline = pop_select(EEG, 'time', baseline_window_in_seconds);
        EEG_response1 = pop_select(EEG, 'time', detection_window_r1);
        EEG_response2 = pop_select(EEG, 'time', detection_window_r2);
        n_epochs = size(EEG.data,3);
        good_trial_lim = good_trial_perc_scale*n_epochs;
        bad_trial_lim = bad_trial_perc_scale*n_epochs;
        bad_channels_pre = find_bad_channels(EEG_baseline,rejchan_threshold_pre, frontal_threshold_shift, bad_trial_thresh, bad_trial_lim, good_trial_lim,pre_peak_threshes_large, pre_peak_threshes_small, 'pre')
        bad_channels_post1 = find_bad_channels(EEG_response1,rejchan_threshold_post, frontal_threshold_shift, NaN, NaN, NaN, [], [], 'post')
        bad_channels_post2 = find_bad_channels(EEG_response2,rejchan_threshold_post, frontal_threshold_shift, NaN, NaN, NaN, [], [], 'post')
        union_bad_elecs = union(bad_channels_pre, union(bad_channels_post1, bad_channels_post2));
        chan_names = {EEG.chanlocs.labels};
        bad_chans = chan_names(union_bad_elecs)
        %pop_eegplot(EEG,1,1,1)
        EEG = pop_select(EEG,'nochannel',bad_chans);
        %plot the average response after removing bad channels
        figname_goodchannels = strcat(reftep_subject,'_good_channels_TEP_data');
        figure; pop_timtopo(EEG, timtopo_timerange, timtopo_topotimes, figname_goodchannels);
        saveas(gcf,append(directory_path,strcat(figname_goodchannels, '.png')));
        close(gcf);
        %save the data with good channels
        goodchans = append(reftep_subject,'_EEG_good_chans.set');
        pop_saveset(EEG, 'filename',goodchans,'filepath',directory_path);
        %% Do ICA (and remove comps in the next sript)
        EEG = pop_tesa_pcacompress(EEG, 'compVal', n_ica_components, 'plot','off' );
        EEG = pop_tesa_fastica(EEG, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'off' );
        with_ica_name = append(reftep_subject,'_EEG_with_ICA.set');
        pop_saveset(EEG, 'filename',with_ica_name,'filepath',directory_path);
    end
end