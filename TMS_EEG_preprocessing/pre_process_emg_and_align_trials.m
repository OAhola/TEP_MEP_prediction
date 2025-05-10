%pre process emg data
clear;
close all;
dtime = string(datetime);
diary_name = string(strcat('preprocessing_diary_eeg_and_emg_alignment_',dtime,'.txt'));
diary_name = strrep(diary_name, ' ', '-');
diary_name = strrep(diary_name, ':', '-');
diary(diary_name)
disp(dtime)
addpath("eeglab2024.2\")
eeglab;
pre_processing_params_path = "pre_processing_parameters_final.xlsx"
pre_processing_params =  readtable(pre_processing_params_path)
params_all = pre_processing_params(strcmp(pre_processing_params.site,"All"),:)
emg_high_pass = params_all.emg_high_pass
emg_preact_window = eval(params_all.emg_preact_window{1})
emg_preact_threshold_muv = params_all.emg_preact_threshold_muv
threshes_extra = emg_preact_threshold_muv + [10 20 30 40 50]
filter_order = params_all.filter_order
baseline_min = params_all.baseline_min
final_time_window_eeg = eval(params_all.final_time_window_eeg{1})
baseline_window = [baseline_min*1000 emg_preact_window(2)*1000] %baseline window in ms
final_time_window_emg = [emg_preact_window(1) final_time_window_eeg(2)]
datapath_base = 'D:\REFTEP_ALL\EEG_preprocessing_data\'
%% go through sites and subjects
sites = {'Tuebingen','Aalto'};
for site=sites
    %go through all the sites and define specific paths
    site_char = char(site);
    params_site = pre_processing_params(strcmp(pre_processing_params.site,site_char),:);
    trigger_names = eval(params_site.trigger_names{1})
    directory_name_site = fullfile(datapath_base,strcat('Preprocessing_',site_char,"\"));
    files_and_folders = dir(directory_name_site)
    is_subfolder = [files_and_folders.isdir]
    folders = files_and_folders(is_subfolder)
    names = {folders.name}
    subject_names = names(contains(names,"sub"))
    %trigger names for the site
    for index = 1:length(subject_names)
        reftep_subject = subject_names(index)
        directory_name = char(fullfile(directory_name_site,reftep_subject,"\")); %the directory of the EEG data to be loaded
        filename_emg = append(reftep_subject,'_EMG_raw_epoched.set');
        EMG = pop_loadset(filename_emg, directory_name); %load the data
        if EMG.nbchan ~=2
            EMG = pop_select(EMG, 'channel', {'EMG1','EMG2'});
            disp("selecting 2 emg channels")
        end
        if EMG.nbchan ~=2
            disp("not 2 emg channels!")
            return
        end
        EMG = pop_rmbase(EMG, baseline_window ,[]);
        EMG = pop_tesa_filtbutter(EMG, emg_high_pass, [], filter_order, 'highpass');
        EMG = pop_select(EMG, 'time', final_time_window_emg);
        EMG = pop_rmbase(EMG, emg_preact_window*1000 ,[]);
        window1 = emg_preact_window*1000 + [0 1]; %add buffer for detrend for pop_tesa_detrend
        window1_s = window1/1000; %in seconds
        EMG_pre = pop_select(EMG, 'time', window1_s);
        EMG_pre = pop_tesa_detrend(EMG_pre, 'linear',emg_preact_window*1000); %detrend the preact window
        EMG_pre = pop_select(EMG_pre, 'time', emg_preact_window); %select the preact window
        if size(EMG.data,3)~=size(EMG_pre.data,3)
            fprintf("check time range cropping")
            return
        end
        n_trials = size(EMG_pre.data,3);
        emg_trials = 1:n_trials;
        n_channels = size(EMG_pre.data,1);
        bad_emg_trials = [];
        %go through all trials
        for i=1:n_trials
            for chan_number=1:n_channels
                data_now = EMG_pre.data(chan_number,:,i);
                diff = max(data_now) - min(data_now);
                if diff > emg_preact_threshold_muv
                %check if the threshold is exceeded within this trial
                    bad_emg_trials(end+1) = i;
                end
            end
        end

        bad_emg_trials = unique(bad_emg_trials);
        fprintf("Found %d bad emg trials\n",length(bad_emg_trials))
        bad_emg_trials_binary = ismember(emg_trials, bad_emg_trials);
        thresh_ind = 1;
        %check how many trials are rejected
        while thresh_ind <= length(threshes_extra)
            if sum(bad_emg_trials_binary) > n_trials/3
                fprintf("Adjusting threshold to %d as %d bad trials were found...\n",threshes_extra(thresh_ind),sum(bad_emg_trials_binary))
                bad_emg_trials = [];
                %go through all trials
                for i=1:n_trials
                    for chan_number=1:n_channels
                        data_now = EMG_pre.data(chan_number,:,i);
                        diff = max(data_now) - min(data_now); %min max diff
                        if diff > threshes_extra(thresh_ind)
                        %check if the threshold is exceeded within this trial
                            bad_emg_trials(end+1) = i;
                        end
                    end
                end
                bad_emg_trials = unique(bad_emg_trials);
                bad_emg_trials_binary = ismember(emg_trials, bad_emg_trials);
            else
                break
            end
            thresh_ind = thresh_ind + 1; %adjust treshold index
        end
        EMG = pop_rejepoch(EMG, bad_emg_trials_binary,0);
        fprintf("Found %d bad emg trials\n",length(bad_emg_trials))
    
        %align EEG and EMG trials
        filename_eeg_cleaned = append(reftep_subject,'_EEG_cleaned.set');
        EEG = pop_loadset(filename_eeg_cleaned,directory_name);
        events_eeg = [EEG.event.urevent];
        events_emg = [EMG.event.urevent];
        %find existing bad trials and reject them
        bad_emg_trials_in_eeg = ~ismember(events_eeg, events_emg);
        bad_eeg_trials_in_emg = ~ismember(events_emg, events_eeg);
        EEG = pop_rejepoch(EEG, bad_emg_trials_in_eeg, 0);
        EMG = pop_rejepoch(EMG, bad_eeg_trials_in_emg, 0);
        events_eeg_cleaned = [EEG.event.urevent];
        events_emg_cleaned = [EMG.event.urevent];
        if size(EEG.data,3) ~= size(EMG.data,3) && ~isequal(events_eeg_cleaned,events_emg_cleaned)
            fprintf("bad trial alignment\n")
            return
        end
    
        %save the final results
        eeg_final_name = char(append(reftep_subject,'_EEG_aligned_final.set'));
        emg_final_name = char(append(reftep_subject,'_EMG_aligned_final.set'));
        time_inds = find(EMG.times > 10 & EMG.times < 50);
        %create a figure of the emg post stim
        figure;
        for ch_ind=1:n_channels
            subplot(n_channels,1,ch_ind)
            hold on
            mean_data = mean(squeeze(EMG.data(ch_ind,time_inds,:)),2);
            std_data = std(squeeze(EMG.data(ch_ind,time_inds,:)),0,2);
            fill([EMG.times(time_inds),fliplr(EMG.times(time_inds))],[mean_data + std_data; flipud(mean_data - std_data)],[0.8, 0.8, 0.8],'EdgeColor','none')
            plot([EMG.times(time_inds)], mean_data)
            xlabel("Time (ms)")
            ylabel("Voltage (uV)")
            legend('std','mean')
            title(string(EMG.chanlocs(ch_ind).labels))
        end
        figname_emg = strcat(reftep_subject,' - emg data');
        sgtitle(figname_emg)
        saveas(gcf,char(append(directory_name,strcat(figname_emg, '.png'))));
        close(gcf);
        %save the data
        pop_saveset(EEG, 'filename',eeg_final_name,'filepath',directory_name);
        pop_saveset(EMG, 'filename',emg_final_name,'filepath',directory_name);
    end
end