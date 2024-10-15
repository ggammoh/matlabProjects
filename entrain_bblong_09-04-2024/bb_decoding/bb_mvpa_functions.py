"""
This module contains custom functions for the binaural beats decoding/mvpa analysis
"""

# load libraries
import numpy as np
import pandas as pd
import mne 
import matplotlib.pyplot as plt
import copy

# libraries for decoding/mvpa/cross validation
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedKFold, cross_val_score



# divide eeg file into equal epochs and label using event annotations
def epoch_bb_raw(
    raw, 
    epoch_len_sec, 
    event_mappings
):

    # get events from annotations
    (events_from_annot, event_dict) = mne.events_from_annotations(raw, event_mappings)

    # create dictionary of time ranges for stimuli
    evlim = copy.copy(events_from_annot)

    # insert ranges for before 0 and 12 
    evlim = np.insert(evlim, 0, [0, 0, 0], axis = 0)
    last_time = int(raw.times[-1] * 500)
    evlim = np.insert(evlim, len(evlim), [last_time, 0 , 12], axis = 0)

    # get time ranges from original events
    stimulus_time_ranges = {evlim[i, 2]: (evlim[i, 0], evlim[i+1, 0]) for i in range(len(evlim) - 1)}

    # create equally spaced epochs
    epoch = mne.make_fixed_length_epochs(
        raw, 
        duration = epoch_len_sec,
        preload = True)

    # assign epoch events based on time ranges
    for ev_i, ev in enumerate(epoch.events):
        for key, value in stimulus_time_ranges.items():
            if value[0] <= ev[0] <= value[1]:
                ev[2] = key
                break
    
    # return epoched eeg
    return epoch
            

def select_roi(
    epoch,
    tfr, 
    stim1, 
    stim2, 
    roi
):

    # select only epochs belonging to stimuli of interest
    epoch_mask = (epoch.events[:, 2] == stim1) | (epoch.events[:, 2] == stim2)

    # select roi
    chan_list = epoch.info.ch_names
    chan_mask = np.asarray([i in roi for i in chan_list])

    
    # epoch_data is raw data, y is labels
    tfr_sel = tfr[epoch_mask, :]
    tfr_sel = tfr_sel[:, chan_mask, :]
    
    labels = epoch.events[:, 2]  
    labels = labels[epoch_mask]


    return tfr_sel, labels




# calculate tfr, pick epochs matching stimuli, return avg tfr 
def calc_bb_tfr(
    epoch, 
    freqs, 
    cycles, 
    dimensions_2_average = None
):
        
    # get data from epochs
    epoch_data = epoch.get_data()
    
    # get sampling rate from eeg 
    sfreq = epoch.info['sfreq']

    # calculate tfr morlet 
    tfr = mne.time_frequency.tfr_array_morlet(
        epoch_data = epoch_data, 
        sfreq = sfreq, 
        freqs = freqs,
        n_cycles = cycles, 
        zero_mean = False, 
        use_fft = True, 
        decim = 1, 
        output = 'power', 
        n_jobs = 10, 
        verbose = None
    )
    
    # average across specified dimensions (default average time)
    if dimensions_2_average is None:
        avg_tfr = tfr
    else:
        avg_tfr = np.mean(tfr, axis = dimensions_2_average)

    # avg_tfr is the features, y is the labels
    return avg_tfr



# classify AT EACH FREQUENCY, CHANNELS ARE USED AS FEATURES
def calc_crossval_score_keepfreq(
    avg_tfr, 
    y, 
    freqs, 
    n_splits, 
    svm_kernel = 'rbf'
):

    # cross validation scheme
    cv = StratifiedKFold(n_splits = n_splits, shuffle=True, random_state=1)

    # svm classifier
    clf = SVC(kernel = svm_kernel)

    # output array of scores at each frequency
    scores = np.zeros(freqs.shape)
    scores_shuffled = np.zeros(freqs.shape)

    # loop through and classify at each frequency
    for freq_i, freq in enumerate(freqs):
            
        # pick data from current frequency 
        X = avg_tfr[..., freq_i].reshape(len(y), -1)

        # shuffle labels to get a null/random for comparison        
        y_shuf = copy.copy(y)
        np.random.shuffle(y_shuf)

        # save mean scores over folds for each frequency 
        scores[freq_i] = np.mean(
            cross_val_score(estimator=clf, X=X, y=y, scoring="roc_auc", cv=cv), axis=0
        )
        
        # save shuffled scores
        scores_shuffled[freq_i] = np.mean(
            cross_val_score(estimator=clf, X=X, y=y_shuf, scoring="roc_auc", cv=cv), axis=0
        )
    
    # return scores and 'null' shuffled scores
    return (scores, scores_shuffled)



# classify AT EACH CHANNEL, FREQUENCIES ARE USED AS FEATURES
def calc_crossval_score_keepchan(
    avg_tfr, 
    y, 
    chan_list, 
    n_splits, 
    svm_kernel = 'rbf'
):
    
    # cross validation scheme
    cv = StratifiedKFold(n_splits = n_splits, shuffle=True, random_state=1)

    # svm classifier
    clf = SVC(kernel = svm_kernel)

    # output array of scores at each frequency
    scores = np.zeros((len(chan_list), 1))
    scores_shuffled = np.zeros(scores.shape)

    # loop through and classify at each frequency
    for chan_i, chan in enumerate(chan_list):
            
        # pick data from current frequency 
        X = avg_tfr[:, chan_i, :].reshape(len(y), -1)

        # shuffle labels to get a null/random for comparison        
        y_shuf = copy.copy(y)
        np.random.shuffle(y_shuf)

        # save mean scores over folds for each frequency 
        scores[chan_i] = np.mean(
            cross_val_score(estimator=clf, X=X, y=y, scoring="roc_auc", cv=cv), axis=0
        )
        
        # save shuffled scores
        scores_shuffled[chan_i] = np.mean(
            cross_val_score(estimator=clf, X=X, y=y_shuf, scoring="roc_auc", cv=cv), axis=0
        )

    # return scores and 'null' shuffled scores
    return (scores, scores_shuffled)
