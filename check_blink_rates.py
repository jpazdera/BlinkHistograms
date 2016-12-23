import numpy as np
from ptsa.data.readers import BaseEventReader
import glob
import json
import os

THRESH = .5
subjs = glob.glob('/Users/jessepazdera/rhino_mount/data/eeg/scalp/ltp/ltpFR/behavioral/events/events_fr_LTP*.mat')
low = []

for s in subjs:
    print os.path.basename(s)[10:-4]
    if s.endswith('backup.mat'):
        continue

    r = BaseEventReader(filename=s, use_reref_eeg=True)
    e = r.read()
    pres = e[np.where(e.type == 'WORD')[0]]

    if len(pres) == 0:
        print 'No event data. Skipping...'
        low.append(s)
        continue

    blink_trials = pres[np.where(pres.artifactMS >= 0)[0]]
    blink_rate = float(len(np.where(blink_trials.artifactMS <= 2000)[0])) / float(len(pres))

    if blink_rate < THRESH:
        low.append(s)

j = dict()
j['low'] = low
with open('/Users/jessepazdera/Desktop/low_blinks_' + str(THRESH) + '.json', 'w') as f:
    json.dump(j, f)