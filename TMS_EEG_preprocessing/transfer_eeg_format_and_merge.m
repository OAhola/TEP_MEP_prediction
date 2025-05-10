%%
clear all;
clc;
dtime = string(datetime);
diary_name = string(strcat('diary_change_data_format_',dtime,'.txt'));
diary_name = strrep(diary_name, ' ', '-');
diary_name = strrep(diary_name, ':', '-');
diary(diary_name)
disp(strcat('the datetime is ', datestr(now, 'dd/mm/yy-HH:MM')))
addpath("eeglab2024.2") %eeglab in the current folder
eeglab %start eeglab
electrode_location_path = 'D:\REFTEP_ALL\Digitization\standard_1005.elc'
%reftep_eeg_excel = "D:\REFTEP_ALL\Data_raw\REFTEP_ppEEG_TueAalto.xlsx";
reftep_eeg_excel = ""
eeg_table = readtable(reftep_eeg_excel); %the table with eeg information
subids = eeg_table.Subject_id;
subids = unique(subids);
for sub_ind = 1:length(subids)
    sub_id = char(subids{sub_ind}); %the current subject id
    if strcmp(sub_id(1),"0")
        site = 'Tuebingen'
    elseif strcmp(sub_id(1),"1")
        site = 'Aalto'
    else
        return
    end
    site_raw = strcat(site,'_raw');
    rootDir = char("D:\REFTEP_ALL\Data_raw");
    eeg_path_on_pc = char(fullfile(rootDir,site_raw));
    site_preprocessing = strcat('Preprocessing_',site);
    eeg_preprocessing_path = char(fullfile('D:\REFTEP_ALL\EEG_preprocessing_data',site_preprocessing));
    datapath_eeg_subject = fullfile(eeg_path_on_pc,['sub-' sub_id ]) %eeg datapath to put data into
    eeg_preprocessing_path_subject = fullfile(eeg_preprocessing_path,['sub-' sub_id ]) %eeg datapath to put merged data into
     %% load subject info
    %go through all the blocks
    indices_subjects_now_eeg = strcmp(eeg_table.Subject_id,sub_id);
    eeg_table_sub = eeg_table(indices_subjects_now_eeg,:);
    blocks_all = eeg_table_sub.Block;
    %blocks = blocks(contains(blocks,'tep')); %tep blocks
    %go through all blocks
    mkdir(datapath_eeg_subject)
    tep_block_nro = 0;
    for block_ind = 1:length(blocks_all)
        block_str = char(blocks_all{block_ind})  % Extract the string from the cell array
        indices_blocks_now = strcmp(eeg_table_sub.Block,block_str)
        eeg_table_block = eeg_table_sub(indices_blocks_now,:)
        session_num = eeg_table_block.EEG_session_index
        subject_in_reftep = append(['REFTEP_' sub_id]);
        if strcmp(site,'Aalto')
            subject_in_reftep = append(['REFTEP' sub_id]); %change
        end
        session_path = char(fullfile(rootDir,strcat("Raw_NeurOne_",site),subject_in_reftep,eeg_table_block.EEG_session_folder))
        EEG = custom_pop_readneurone_tue_noeloc(session_path, session_num);
        channel_names_original = {EEG.chanlocs.labels};
        channel_names_original_lower = lower(channel_names_original);
        channels = channel_names_original(~contains(channel_names_original_lower,'input')) %select proper channels, this one is bad in tuebingen
        EEG = pop_select(EEG, 'channel', channels);
        EEG = pop_chanedit(EEG,'lookup',electrode_location_path,'nosedir','+Y');
        EEG = pop_editset(EEG, 'setname', ['sub-' sub_id '-' block_str]);
        %remove "Input" channel if it is there
        EEG = eeg_checkset(EEG);
        %% save eeg data to eeglab format
        filetosave = append(['sub-' sub_id '_task-' block_str '_eeg.set']);
        pop_saveset( EEG, 'filename',filetosave,'filepath',datapath_eeg_subject);
        if contains(block_str,'tep')
            if tep_block_nro == 0
                tep_block_nro = tep_block_nro + 1;
                EEG_merged = EEG;
            else
                EEG_merged = pop_mergeset(EEG_merged, EEG);
            end
        end
    end
    mkdir(eeg_preprocessing_path_subject);
    EEG_merged = pop_editset(EEG_merged, 'setname', ['sub-' sub_id '-teps_merged']);
    filetosave_merged = append(['sub-' sub_id '_task-tep_all_eeg.set']);
    pop_saveset(EEG_merged, 'filename',filetosave_merged,'filepath',eeg_preprocessing_path_subject);
end
diary off