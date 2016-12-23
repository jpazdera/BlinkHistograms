import numpy as np
from ptsa.data.readers import BaseEventReader
import glob
import matplotlib.pyplot as plt
import os
import matplotlib.backends.backend_pdf
import json
import traceback

subjs = glob.glob('/Users/jessepazdera/rhino_mount/data/eeg/scalp/ltp/ltpFR/behavioral/events/events_fr_LTP*.mat')
low = []

try:
    with open('/Users/jessepazdera/Desktop/low_blinks.json') as f:
        low = json.load(f)['low']
except:
    traceback.print_exc()

for s in subjs:
    try:
        if s.endswith('backup.mat') or s in low:
            continue

        snum = os.path.basename(s)[10:-4]
        print snum
        r = BaseEventReader(filename=s, use_reref_eeg=True)
        e = r.read()
        pres = e[np.where(e.type == 'WORD')[0]]

        if len(pres) == 0:
            print 'No event data. Skipping...'
            low.append(s)
            continue

        blink_trials = pres[np.where(pres.artifactMS >= 0)[0]]
        blink_rate = float(len(np.where(blink_trials.artifactMS <= 2000)[0])) / float(len(pres))
        print 'Blink Rate: ', blink_rate
        if blink_rate < .2:
            print 'Skipping...'
            low.append(s)
            continue

        print 'High blink rate. Creating histogram...'

        plt.figure(facecolor='white')
        plt.subplot(311)
        plt.hist(blink_trials.artifactMS, range(0, 4401, 40), color='b', normed=True)
        plt.axvline(3000, color='k', linewidth=2)
        plt.title(snum + ' (' + str(int(blink_rate*100)) + '%) - Blink Onset Time')
        plt.xlabel('MS after word onset')
        plt.ylabel('Frequency')

        plt.subplot(312)
        plt.hist(blink_trials.artifactMeanMS, range(0, 4401, 40), color='r', normed=True)
        plt.axvline(3000, color='k', linewidth=2)
        plt.title('Avg Time of Blink Samples')
        plt.xlabel('MS after word onset')
        plt.ylabel('Frequency')

        plt.subplot(313)
        hist3 = plt.hist(blink_trials.artifactFrac * 100, range(0, 101, 1), color='g', normed=True)
        plt.title('Percent of Samps w/Blinks')
        plt.xlabel('Percent')
        plt.ylabel('Frequency')

        plt.tight_layout()

    except:
        traceback.print_exc()

try:
    j = dict()
    j['low'] = low
    with open('/Users/jessepazdera/Desktop/low_blinks.json', 'w') as f:
        json.dump(j, f)
except:
    traceback.print_exc()

pdf = matplotlib.backends.backend_pdf.PdfPages('/Users/jessepazdera/Desktop/Blink Stats.pdf')
for fig in xrange(1, plt.figure().number):
    pdf.savefig(fig)
pdf.close()

