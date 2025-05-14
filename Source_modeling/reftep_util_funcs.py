import h5netcdf
import mne
import mne_connectivity
import os
import numpy as np
import sklearn
from sklearn.decomposition import PCA
import sklearn.metrics
import collections
from scipy.signal import find_peaks
from sklearn.cluster import KMeans
import xml.etree.ElementTree as ET
import scipy.spatial.transform
import matplotlib.pyplot as plt


#Functions for dipole fitting
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def lsq_dipole(forward, evoked, tmin=None, tmax=None):
    """
    forward: the forward solution
    evoked: evoked structure
    tmin, tmax = minimum and maximum time, can be tmin=tmax
    """
    L = forward['sol']['data']
    n_dips_x_oris = int(L.shape[1]) #number of dipoles x orientations
    L = L - np.mean(L,axis=0) #average reference the leadfield
    exp_var = -np.inf #initialize that the explained variance is -np.inf
    evoked_full_data = evoked.data - np.mean(evoked.data, axis=0) #average reference EEG data
    n_channels = evoked.data.shape[0] #number of channels
    #if tmin and tmax have not been specified, then use the whole time range of the evoked response (or single trial)
    if tmin is None and tmax is None:
        evoked_cropped = evoked
    else:
        evoked_cropped = evoked.copy().crop(tmin, tmax)
    for time_ind, time in enumerate(evoked_cropped.times): #go through the times
        y_measured = evoked_cropped.data[:,time_ind] - np.mean(evoked_cropped.data[:,time_ind], axis=0) #the measured topography at the current time in average reference
        for lf_ind, dipole_start_index in enumerate(np.arange(0,int(n_dips_x_oris),3)): #go through lead field triplets, i.e. positions with orientations
            pickcol = L[:,dipole_start_index:dipole_start_index+3]#picked triplet here
            free_orientation_stc_now = np.matmul(np.linalg.pinv(pickcol),y_measured) #q_hat = pinv(L)y
            y_predicted = np.matmul(pickcol,free_orientation_stc_now) #y_hat = Lq_hat
            exp_var_new = sklearn.metrics.explained_variance_score(y_measured, y_predicted) #expvar between y and y_hat
            if exp_var_new > exp_var: #if the best explained variance was exceeded, then save the dipole info as the best one
                exp_var = exp_var_new #best new expvar
                #update best position
                best_pos_ind = lf_ind
                best_match = dipole_start_index
                best_topo = y_predicted #save the best (so far) topography
                best_free_ori_stc_now = free_orientation_stc_now #dipole orientation
                best_time = time #the time of the fit


    evoked_fit = mne.EvokedArray(data=best_topo.reshape(n_channels,1), info=evoked.info, tmin=best_time, nave=1, kind='average', baseline=None) #the dipolar acitity
    #project out the dipolar activity from the data and return the resulting residual
    lf_cols = L[:,best_match:best_match+3] #best matching source
    tc_of_source = np.matmul(np.linalg.pinv(lf_cols),evoked_full_data) #full time course of that source
    source_projected = np.matmul(lf_cols,tc_of_source) #time course of the source projected to sensor space
    full_residual_data = evoked_full_data - source_projected #calculate the residual (dipolar activity removed)
    full_residual_data = full_residual_data - np.mean(full_residual_data,axis=0) #average reference
    evoked_residual = mne.EvokedArray(data=full_residual_data, info=evoked.info, tmin=evoked.times[0], nave=1, kind='average', baseline=None) #residual data

    #extract dipole information and create a dipole object
    dipole_amplitude = np.linalg.norm(best_free_ori_stc_now) #dipole amplitude
    dipole_pos = forward['source_rr'][best_pos_ind] #dipole position in mri
    dipole_ori = best_free_ori_stc_now / dipole_amplitude #dipole orientation normalized (q/amplitude)
    dipole = mne.Dipole([best_time], [dipole_pos], [dipole_amplitude], [dipole_ori], [exp_var*100], name=f"dipole_{best_pos_ind}") #create dipole object (mne.Dipole)
    return dipole, best_free_ori_stc_now, evoked_fit, best_match, best_pos_ind, evoked_residual

def lsq_dipole_to_pos(forward, evoked, tmin, tmax, pos_ind, n_times, ori_fixed=None):
    """
    This function fits a dipole to a pre-defined location and gives it a free orientation
    """
    start_of_triplet_ind = pos_ind*3
    n_channels = evoked.data.shape[0] #number of channels
    L = forward['sol']['data']
    L = L - np.mean(L,axis=0) #average reference the leadfield
    exp_var = -np.inf #initialize that the explained variance is -np.inf
    evoked_cropped = evoked.copy().crop(tmin, tmax)
    if len(evoked_cropped.times) != n_times:
        print(len(evoked_cropped.times), n_times)
        return False
    for time_ind, time in enumerate(evoked_cropped.times):
        y_measured = evoked_cropped.data[:,time_ind] - np.mean(evoked_cropped.data[:,time_ind], axis=0) #the measured topography at the current time in average reference
        pickcol = L[:,start_of_triplet_ind:start_of_triplet_ind+3] #picked triplet here
        if ori_fixed is None:
            orientation_stc_now = np.matmul(np.linalg.pinv(pickcol),y_measured) #q_hat = pinv(L)y
        else:
            pickcol_projected = np.matmul(pickcol,ori_fixed)[:, np.newaxis] #project the leadfield to the fixed orientation
            amplitude = np.matmul(np.linalg.pinv(pickcol_projected),y_measured)[0] #amplitude of the fixed source
            orientation_stc_now = ori_fixed*amplitude #the source strengths in all dimensions
        y_predicted = np.matmul(pickcol,orientation_stc_now) #y_hat = Lq_hat
        exp_var_new = sklearn.metrics.explained_variance_score(y_measured, y_predicted) #expvar between y and y_hat
        if exp_var_new > exp_var: #if the best explained variance was exceeded, then save the dipole info as the best one
            exp_var = exp_var_new
            best_topo = y_predicted
            best_ori_stc_now = orientation_stc_now
            best_time = time

    #extract dipole information and create a dipole object
    dipole_amplitude = np.linalg.norm(best_ori_stc_now) #dipole amplitude
    #print(dipole_amplitude)
    evoked_fit = mne.EvokedArray(data=best_topo.reshape(n_channels,1), info=evoked.info, tmin=best_time, nave=1, kind='single_epoch', baseline=None)
    dipole_pos = forward['source_rr'][pos_ind] #dipole position in mri
    dipole_ori = best_ori_stc_now / dipole_amplitude #dipole orientation normalized (q/amplitude)
    dipole = mne.Dipole([best_time], [dipole_pos], [dipole_amplitude], [dipole_ori], [exp_var*100], name=f"dipole_{pos_ind}") #create dipole object (mne.Dipole)
    return dipole, best_ori_stc_now, evoked_fit



def define_fitting_ranges(center_times, percent_increase, min_time):
    percent_decrease = 1-percent_increase
    center_times = np.array(list(center_times.values()))
    time_edges_pos = []
    for i in range(len(center_times)-1):
        t1  = center_times[i]
        t2 = center_times[i+1]
        t_edge = t1 + percent_increase*(t2-t1)
        time_edges_pos.append(t_edge)
    last_edge = center_times[-1] + percent_increase*(center_times[-1]-center_times[-2])
    time_edges_pos.append(last_edge)


    time_edges_neg = []
    first_edge = center_times[0] - percent_increase*(center_times[1]-center_times[0])
    time_edges_neg.append(first_edge)
    for i in range(len(center_times)-1):
        t1  = center_times[i]
        t2 = center_times[i+1]
        t_edge = t1 + percent_decrease*(t2-t1)
        time_edges_neg.append(t_edge)
    time_edges_neg_adjusted = [time if time >= min_time else min_time for time in time_edges_neg]
    return time_edges_neg_adjusted, time_edges_pos

    
def get_mep_ptps(emg_epochs,tmin,tmax):
    data = emg_epochs.copy().crop(tmin=tmin, tmax=tmax).get_data(copy=True)
    mins = np.min(data, axis=2)
    maxs = np.max(data, axis=2)
    min_max_diffs = maxs - mins
    mean_min_max_diffs = np.mean(min_max_diffs,axis=0)
    maxind = np.argmax(mean_min_max_diffs) #which one has on average a higher amplitude
    return min_max_diffs[:,maxind] #return those amplitudes (min-max differences)
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


def get_psd_over_freqs(epochs,fmin,fmax,bw_scale):
    sfreq = epochs.info["sfreq"] #sampling frequency
    data = epochs.get_data(copy=True) #epoch data
    n_times = data.shape[-1] #number of time points
    bw = bw_scale * (sfreq / n_times) #used bandwidth
    psds, freqs = mne.time_frequency.psd_array_multitaper(data, sfreq, fmin=fmin, fmax=fmax, bandwidth=bw, output='power')
    return psds, freqs

def get_individual_alpha_peak(epochs,bw_scale,fmin,fmax):
    psd, freqs = get_psd_over_freqs(epochs,fmin=fmin,fmax=fmax,bw_scale=bw_scale) #psd from fmin to fmax
    mean_psd = np.mean(psd, axis=0)
    mean_psd = np.mean(mean_psd, axis=0) #mean psd over epochs and channels
    peaks, _ = find_peaks(mean_psd) #find peaks in the psd
    if len(peaks) == 0:
        return False, freqs #no individual peak found
    if len(peaks) == 1:
        max_peak_ind = 0
        alpha_peak = freqs[peaks[max_peak_ind]] #peak found!
    else:
        peaks_vals = mean_psd[peaks]
        max_peak_ind = np.argmax(peaks_vals)
        alpha_peak = freqs[peaks[max_peak_ind]]
        print("more than one peak found, picked one with highest value")
    fig, axs  = plt.subplots()
    axs.plot(freqs,mean_psd)
    axs.axvline(alpha_peak, mean_psd[max_peak_ind])
    plt.show()
    return alpha_peak, freqs

def get_freq_ranges_based_on_alpha_peak(alpha_peak, freq_range_names):
    if len(freq_range_names) == 7:
        alpha = [alpha_peak-2.5, alpha_peak+2.5]
        min_alpha = np.min(alpha)
        theta = [min_alpha-4, min_alpha]
        max_alpha = np.max(alpha)
        low_beta = [max_alpha, max_alpha+10]
        max_low_beta = np.max(low_beta)
        high_beta = [max_low_beta, max_low_beta+10]
        max_high_beta = np.max(high_beta)
        low_gamma = [max_high_beta, max_high_beta+20]
        max_low_gamma = np.max(low_gamma)
        mid_gamma  = [max_low_gamma, max_low_gamma+30]
        max_mid_gamma = np.max(mid_gamma)
        high_gamma  = [max_mid_gamma, max_mid_gamma+50]
        frequency_ranges = [theta, alpha, low_beta, high_beta, low_gamma, mid_gamma, high_gamma]
    elif len(freq_range_names) == 6:
        alpha = [alpha_peak-2.5, alpha_peak+2.5]
        min_alpha = np.min(alpha)
        theta = [min_alpha-4, min_alpha]
        max_alpha = np.max(alpha)
        low_beta = [max_alpha, max_alpha+10]
        max_low_beta = np.max(low_beta)
        high_beta = [max_low_beta, max_low_beta+10]
        max_high_beta = np.max(high_beta)
        low_gamma = [max_high_beta, max_high_beta+35]
        max_low_gamma = np.max(low_gamma)
        high_gamma  = [max_low_gamma, max_low_gamma+50]
        frequency_ranges = [theta, alpha, low_beta, high_beta, low_gamma, high_gamma]
    elif len(freq_range_names) == 5:
        alpha = [alpha_peak-2.5, alpha_peak+2.5]
        min_alpha = np.min(alpha)
        theta = [min_alpha-4, min_alpha]
        max_alpha = np.max(alpha)
        low_beta = [max_alpha, max_alpha+10]
        max_low_beta = np.max(low_beta)
        high_beta = [max_low_beta, max_low_beta+10]
        max_high_beta = np.max(high_beta)
        gamma = [max_high_beta, max_high_beta+40]
        frequency_ranges = [theta, alpha, low_beta, high_beta, gamma]
    elif len(freq_range_names) == 4:
        alpha = [alpha_peak-2.5, alpha_peak+2.5]
        min_alpha = np.min(alpha)
        theta = [min_alpha-4, min_alpha]
        max_alpha = np.max(alpha)
        beta = [max_alpha, max_alpha+20]
        max_beta = np.max(beta)
        gamma = [max_beta, max_beta+50]
        frequency_ranges = [theta, alpha, beta, gamma]
    if len(frequency_ranges) != len(freq_range_names):
        return False
    freq_range_dict = {name:freq_range for name, freq_range in zip(freq_range_names, frequency_ranges)}
    return frequency_ranges, freq_range_dict


def get_intersection_of_good_trial_inds(basepath, response_types, feature_path, subject_path, reftep_subject_name, subject, n_trials, gof_thresh_average, gof_thresh_single, mep_thresh_average, mep_thresh_single):
    #initialize lists
    inds_under_thresh = []
    gofs_all = []
    gof_threshes = []
    checkers = []
    meps_filename = feature_path + subject + "/" + reftep_subject_name + "emg_ptps.npy"
    for response in response_types: #go through responses (p30, n45, p60 and emg)
        if response == 'emg':
            #load emg peak to peaks
            meps_all = np.load(meps_filename)
            meps_mean = np.mean(meps_all)
            if mep_thresh_single is not None:
                mep_inds_under_thresh = set(i for i, mep_amplitude in enumerate(meps_all) if mep_amplitude < mep_thresh)
            else:
                mep_inds_under_thresh = set()
            if mep_thresh_average is not None: #restrict on the average gof
                if meps_mean <= mep_thresh_average:
                    checkers.append(True)
                else:
                    checkers.append(False)
            else:
                checkers.append(True)
        else:
            #load dipole amplitudes
            dip_path = os.path.join(subject_path, f"{reftep_subject_name}_dipoles", f"dipole_{response}/")
            mean_dipole_gof = mne.read_dipole(basepath + subject + "/" + reftep_subject_name + "peakfit_dipole_" + response.lower(), verbose=False).gof[0]
            if gof_thresh_average is not None: #restrict on the average gof
                if mean_dipole_gof >= gof_thresh_average:
                    checkers.append(True)
                else:
                    checkers.append(False)
            else:
                checkers.append(True)
            if gof_thresh_single is not None:
                gof_threshes.append(mean_dipole_gof-gof_thresh_single)
                gofs_all.append([mne.read_dipole(dip_path + subject + response + f"fixed_peakfit_trial_free_ori{False}{trial}", verbose=False).gof[0] for trial in range(n_trials)])
    inds_under_thresh = [set(i for i, gof in enumerate(gofs) if gof < gof_thresh) for gofs, gof_thresh in zip(gofs_all, gof_threshes)]
    inds_under_thresh.append(mep_inds_under_thresh) #add the bad mep trial indices
    union_of_indices = set.union(*inds_under_thresh)
    bad_inds = list(union_of_indices)
    good_inds = [ind for ind in range(n_trials) if ind not in bad_inds]
    return good_inds, checkers
#------------------------------------------------------------------------------------------------------------------------------------------------------------------     









#Functions for getting estimates for good peak fitting times
#------------------------------------------------------------------------------------------------------------------------------------------------------------------

def get_aucs_for_chs_around_times(evoked, channels, times, t_shift):
    aucs = []
    for channel, time in zip(channels,times):
        tmin, tmax = time-t_shift, time+t_shift
        evoked_now = evoked.copy().crop(tmin, tmax).pick(channel) #channel data at a time range
        data_absed = np.abs(evoked_now.data[0]) #absolute values of the data
        auc_full_range = np.trapz(data_absed) #auc of the time range
        y_start, y_end = data_absed[0], data_absed[-1] #starting and ending value in the time range
        y = np.array([np.min(data_absed) for _ in range(len(data_absed))])
        #y = np.linspace(start=y_start,stop=y_end,num=len(data_absed), endpoint=True) #line from y_start to y_end (including y_end)
        auc_quadrilateral = np.trapz(y) #auc of the quadrilateral
        difference = auc_full_range-auc_quadrilateral #calculate the difference
        if difference < 0: #sanity check
            return False
        else:
            aucs.append(difference)
    return np.array(aucs)


def get_peak_amplitudes_for_times(evoked,channels,times):
    amplitudes = []
    for channel, time in zip(channels,times):
        evoked_now = evoked.copy().crop(time, time).pick(channel) #channel data at a time range
        data_absed = np.abs(evoked_now.data[0]) #absolute values of the data
        amplitudes.append(data_absed[0]) #value at time
    return np.array(amplitudes)


def get_weights_on_time_counts(times):
    n_times = len(times) #number of time points in the cluster
    n_bins = int(np.ceil(np.sqrt(n_times))) #define the number of bins
    _, bin_edges = np.histogram(times, bins = n_bins) #edges of the bins
    bin_indices = np.digitize(times, bin_edges) #bin the data
    counts = collections.Counter(bin_indices)
    #number of samples in the bin where the time point is in
    counts_for_times = np.array([counts[ind] for ind in bin_indices])
    return counts_for_times

def get_weights(inputs):
    weights = np.ones(shape=inputs[0].shape)
    for arr in inputs:
        arr_scaled = scale_with_max(arr)
        weights = weights*arr_scaled
    return weights

def scale_with_max(arr):
    #scales all the values
    return arr/np.max(arr)

def apply_weighting(times, inputs):
    weights = get_weights(inputs)
    weighted_times = times*weights #weighted_times[i] = times[i]*weights[i]
    scale = np.sum(weights) #normalizer
    weighted_mean =  np.sum(weighted_times)/scale
    return weighted_times, weights, weighted_mean

def clustered_peak_times(evoked, tmin, tmax, n_clusters, dropped_peaks=[], prev_peak_indices_list=[], picks=None):
    if picks is not None:
        evokedp = evoked.copy().pick(picks) #define channel space to detect peaks from
    else:
        evokedp = evoked.copy() #no picks
    evoked2 = evokedp.copy().crop(tmin,tmax)
    times = evoked2.times
    #evoked_data_absed = np.abs(evoked2.data)
    ch_names = evoked2.info['ch_names']

    peak_indices_list = [] #store peak indices here
    peak_times_list = [] #store peak times here
    origin_info_list = [] #track which array each peak belongs to

    for ind, arr in enumerate(evoked2.data):
        if len(prev_peak_indices_list) == 0:
            peaks_max, properties_max = find_peaks(arr)
            peaks_min, properties_max = find_peaks(-arr)
            abs_arr = np.abs(arr)
            peaks = np.concatenate((peaks_max, peaks_min)) #combine peaks of maxima and minima
            #peak_amplitudes = abs_arr[peaks_conc]
            #n_components = n_clusters + 1
            #top_n_inds = np.argsort(peak_amplitudes)[-n_components:]
            #peaks = np.sort(peaks_conc[top_n_inds]) #sort the values from smallest to largest
        else:
            peaks = prev_peak_indices_list[ind] #don't need to recompute the peaks if they already have been detected
        peaks_cleaned = np.array([peak for peak in peaks if (ind, peak) not in dropped_peaks]) #take those peaks that have not been dropped
        peak_indices_list.append(peaks_cleaned)
        origin_info_list.extend([ind]*len(peaks_cleaned))
        if len(peaks_cleaned) > 0:
            peak_times_list.append(times[peaks_cleaned])
        else:
            peak_times_list.append([])

    #flatten for clustering
    all_peaks_inds = np.concatenate(peak_indices_list).reshape(-1,1)
    all_peaks = np.concatenate(peak_times_list).reshape(-1,1)
    all_origins = np.array(origin_info_list)


    # Output peak indices for each array
    #for idx, peaks in enumerate(peak_indices_list):
        #print(f"Peaks in array {ch_names[idx]}: {peaks}")

    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10) #k means with 3 clusters and random state 42
    kmeans.fit(all_peaks)
    clusters = kmeans.predict(all_peaks) #get cluster labels

    #group by clusters
    cluster_dict = {}
    for i, label in enumerate(clusters):
        if label not in cluster_dict: #does not exists yet so create a key
            cluster_dict[label] = []
        cluster_dict[label].append(all_peaks[i][0])

    #sort the clusters by mean and apply weighting
    sorted_weighted_cluster_means = []
    for cluster_id in range(n_clusters):
        cluster_peaks_times = all_peaks[clusters==cluster_id][:,0] #peaks for this cluster
        cluster_peaks_inds = all_peaks_inds[clusters==cluster_id][:,0] #indices of peaks for this cluster
        N_c = len(cluster_peaks_times) #number of elements in the cluster
        cluster_origins = all_origins[clusters==cluster_id] #origins of peaks for this cluster
        chs = [ch_names[origin] for origin in cluster_origins] #channels of the peaks in the cluster
        # get weighting values
        #aucs_for_times = get_aucs_for_chs_around_times(evoked_longer_range, chs, cluster_peaks_times, t_shift) #aucs around peak times with t_shift
        amplitudes = get_peak_amplitudes_for_times(evoked2, chs, cluster_peaks_times)
        counts_for_times = get_weights_on_time_counts(cluster_peaks_times) #weights for time distribution
        _ , weights, w_mean = apply_weighting(cluster_peaks_times, [amplitudes, counts_for_times]) #apply weighting
        normalized_weights = weights/(np.linalg.norm(weights)*N_c) #normalize in-class weights for comparisons across cluster and give less weight for larger clusters.
        sorted_weighted_cluster_means.append((w_mean,cluster_peaks_inds.flatten(),cluster_peaks_times.flatten(),cluster_origins,normalized_weights))

    #sort by the mean
    sorted_weighted_cluster_means.sort(key=lambda x: x[0]) #sort by the first value, i.e. the cluster mean
    weighted_means = np.array([sorted_weighted_cluster_means[cluster_id][0] for cluster_id in range(n_clusters)]) #weighted means for each cluster
    cluster_peak_times = [sorted_weighted_cluster_means[cluster_id][2] for cluster_id in range(n_clusters)] #values of peak times in each cluster


    recompute=False #init parameter for telling the script to recompute or not after this run
    min_weight = np.inf #initialize value for comparison
    for _, (_, cluster_peaks, _, cluster_origins, normalized_cluster_weights) in enumerate(sorted_weighted_cluster_means):
        #check if peaks from the same channel are in the same cluster
        duplicates = [] 
        origin_counts = collections.Counter(cluster_origins)
        for channel_index in origin_counts:
            if origin_counts[channel_index] > 1: #more than one peak from one channel in a cluster
                recompute=True #need to recompute
                duplicates.extend([i for i, value in enumerate(cluster_origins) if value==channel_index])
        #remove the "worst" duplicate out of all
        for duplicate_index in duplicates:
            weight_in_cluster = normalized_cluster_weights[duplicate_index]
            if weight_in_cluster < min_weight:
                worst_one = (cluster_origins[duplicate_index],cluster_peaks[duplicate_index])
    if recompute: #add the worst one to be removed
        dropped_peaks.append(worst_one)


    return weighted_means, cluster_peak_times, recompute, peak_indices_list, dropped_peaks, kmeans.inertia_
#------------------------------------------------------------------------------------------------------------------------------------------------------------------