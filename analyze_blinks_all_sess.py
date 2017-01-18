import os
import glob
import json as js
import matplotlib.pyplot as plt
from ptsa.data.readers import BaseEventReader
import numpy as np
from scipy.signal import butter, filtfilt


def detect_blinks(subj, exp, sys='EGI', sample_rate=500):
    if sys == 'EGI':
        eog_chans = [('025', '127'), ('008', '126')]
        blink_thresh = 100
        gain = .0762963
    elif sys == 'Biosemi':
        eog_chans = [('EXG3', 'EXG1'), ('EXG4', 'EXG2')]
        blink_thresh = 1500
        gain = 1
    else:
        print 'ERROR: Cannot run on unknown system.'
        return

    if exp == 'ltpFR':
        dur = 4500 * sample_rate / 1000
        stim_off_time = 2000 * sample_rate / 1000  # actually 3000 ms, but only counts towards stat if blink in first 2k
    elif exp == 'ltpFR2':
        dur = 4000 * sample_rate / 1000
        stim_off_time = 1600 * sample_rate / 1000
    else:
        print 'ERROR: Unknown experiment.'
        return

    sessions = '/Users/jessepazdera/rhino_mount/data/eeg/scalp/ltp/%s/%s/session_*' % (exp, subj)
    sessions = glob.glob(sessions)

    pres_count = 0
    blink_count = np.zeros(dur)
    # beginning_blink_counts = np.zeros(dur)
    # no_beginning_blink_counts = np.zeros(dur)

    beginning_period_dur1 = 100 * sample_rate / 1000
    beginning_period_dur2 = 200 * sample_rate / 1000
    beginning_period_dur3 = 300 * sample_rate / 1000
    blink_rate = 0
    blink_rate_1= 0
    blink_rate_2 = 0
    blink_rate_3 = 0

    for sess in sessions:
        print sess

        event_path = sess + '/events.mat'
        if not os.path.isfile(event_path):
            continue
        try:
            events = BaseEventReader(filename=event_path).read()

            reref = glob.glob(sess + '/eeg/eeg.reref/*')
            if len(reref) == 0:
                continue
            eeg_path = os.path.splitext(reref[0])[0]

            artifact_mask = None
            for i in range(len(eog_chans)):
                ch = eog_chans[i]
                if isinstance(ch, tuple):
                    eeg1 = np.fromfile(eeg_path + '.' + ch[0], 'int16') * gain
                    eeg2 = np.fromfile(eeg_path + '.' + ch[1], 'int16') * gain
                    eeg = eeg1 - eeg2
                else:
                    eeg = np.fromfile(eeg_path + '.' + ch, 'int16') * gain

                if artifact_mask is None:
                    artifact_mask = np.empty((len(eog_chans), len(eeg)))

                artifact_mask[i] = find_blinks(eeg, blink_thresh)

            blinks = np.any(artifact_mask, 0)

            for i in range(len(events)):
                if events[i].type != 'WORD' or events[i].eegoffset < 0:
                    continue

                pres_count += 1
                ev_blink = blinks[events[i].eegoffset:events[i].eegoffset + dur]
                blink_count += ev_blink

                if np.any(ev_blink[:stim_off_time]):
                    blink_rate += 1
                    if np.any(ev_blink[beginning_period_dur1:stim_off_time]):
                        blink_rate_1 += 1
                    if np.any(ev_blink[beginning_period_dur2:stim_off_time]):
                        blink_rate_2 += 1
                    if np.any(ev_blink[beginning_period_dur3:stim_off_time]):
                        blink_rate_3 += 1
                '''
                if np.any(ev_blink[:beginning_period_dur]):
                    beginning_blink_counts += ev_blink
                else:
                    no_beginning_blink_counts += ev_blink
                '''
        except:
            continue
    if pres_count > 0:
        blink_rate /= float(pres_count)
        blink_rate_1 /= float(pres_count)
        blink_rate_2 /= float(pres_count)
        blink_rate_3 /= float(pres_count)
    else:
        blink_rate = -999
        blink_rate1 = -999
        blink_rate2 = -999
        blink_rate3 = -999
    return blink_count, pres_count, blink_rate, blink_rate_1, blink_rate_2, blink_rate_3

def find_blinks(data, thresh):
    """
    Locates the blinks in an EEG signal. It does so by maintaining two running averages - one that changes quickly
    and one that changes slowly.

    The "slow" running average gives each new sample a weight of .025 and the current running average a weight of
    .975, then adds them to produce the new average. This average is used as a type of baseline measure against
    which the "fast" running average is calculated.

    The fast average tracks the divergence of the signal voltage from the slow/baseline average, in order to detect
    large spikes in voltage. In calculating each new average, it gives the new sample a weight of .5 and the
    current fast running average a weight of .5, allowing it to be heavily influenced by rapid voltage changes.

    Blinks are identified as occurring on any sample where the fast average exceeds the given threshold, i.e.
    whenever voltage displays a large and rapid divergence from baseline.

    :param data: A numpy array containing the EEG signal that will be searched for blinks.
    :param thresh: The mV threshold used for determining when a blink has occurred.
    :return: A numpy array containing one boolean per EEG sample, indicating whether that sample contains a blink.
    """
    # Create an array for each of the two running averages
    num_samples = len(data)
    fast = np.empty(num_samples)
    slow = np.empty(num_samples)

    # Calculate the starting "fast" and "slow" averages
    start_mean = np.mean(data[0:10])
    fast[0] = .5 * (data[0] - start_mean)
    slow[0] = .975 * start_mean + .025 * data[0]

    # Track the running averages across all samples
    for i in range(1, num_samples):
        fast[i] = .5 * fast[i - 1] + .5 * (data[i] - slow[i - 1])
        slow[i] = .975 * slow[i - 1] + .025 * data[i]

    # Mark whether each sample's "fast average" exceeded the threshold
    return abs(fast) >= thresh


def butter_filt(data, freq_range, sample_rate=500, filt_type='bandstop', order=4):
    """
    Designs and runs an Nth order digital butterworth filter on an array of data.

    NOTE: In order to match the original MATLAB implementation of the filtfilt function, the padlen argument must be
    set to 3*(max(len(Ab), len(Bb))-1)). Default padlen in SciPy is 3*max(len(a), len(b)), which will cause it to return
    different values from our MATLAB scripts if not manually set (i.e. default padlen for SciPy is always 3 greater than
    the default padlen for MATLAB).

    :param data: An array containing the data to be filtered
    :param freq_range: A single frequency to filter on or a 2D matrix where each row is a frequency range to filter on
    :param sample_rate: The sampling rate of the EEG recording
    :param filt_type: The type of filter to run - can be 'bandstop', 'highpass', 'lowpass', or 'bandpass'
    :param order: The order of the filter
    :return: The filtered data
    """
    # Calculate Nyquist frequency
    nyq = sample_rate/2.
    # Convert the frequency range to a 2D matrix if it was input as an integer
    freq_range = np.array([[freq_range]]) if isinstance(freq_range, (int, float)) else np.array(freq_range)
    freq_range = freq_range / nyq
    # Get the Butterworth values and run the filter for zero phase distortion
    for i in range(np.shape(freq_range)[0]):
        Bb, Ab = butter(order, freq_range[i, :], btype=filt_type)
        pad = 3 * (max(len(Ab), len(Bb)) - 1)  # calculate padlen like MATLAB instead of SciPy
        data = filtfilt(Bb, Ab, data, padlen=pad)
    return data

if __name__ == "__main__":
    subj = 'LTP249'  # raw_input('Enter subject ID: ')
    exp = 'ltpFR'  # raw_input('Enter experiment name: ')
    sys = 'EGI'
    blink_count, pres_count, br, br1, br2, br3 = detect_blinks(subj, exp, sys=sys)

    json = {}
    json['pres_count'] = pres_count
    json['blink_count'] = blink_count.tolist()
    json['blink_rate'] = br
    json['blink_rate_1'] = br1
    json['blink_rate_2'] = br2
    json['blink_rate_3'] = br3

    outfile = '/Users/jessepazdera/Desktop/Blink Stats/blink_count_%s.json' % subj
    with open(outfile, 'w') as f:
        js.dump(json, f)
    plt.plot(blink_count/pres_count*100)
    plt.show()
