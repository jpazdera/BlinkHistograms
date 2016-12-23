import numpy as np
from ptsa.data.readers import BaseEventReader
import glob
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import json
import traceback

subjs = glob.glob('/Users/jessepazdera/rhino_mount/data/eeg/scalp/ltp/ltpFR/behavioral/events/events_fr_LTP*.mat')
low = []
ams = None
amms = None
af = None

try:
    with open('/Users/jessepazdera/Desktop/low_blinks.json') as f:
        low = json.load(f)['low']
except:
    traceback.print_exc()

for s in subjs:
    try:
        if s not in low:
            continue

        r = BaseEventReader(filename=s, use_reref_eeg=True)
        e = r.read()

        if ams is None:
            ams = e[np.where(e.type == 'WORD')[0]].artifactMS
            amms = e[np.where(e.type == 'WORD')[0]].artifactMeanMS
            af = e[np.where(e.type == 'WORD')[0]].artifactFrac
        else:
            ams = np.concatenate((ams, e[np.where(e.type == 'WORD')[0]].artifactMS))
            amms = np.concatenate((amms, e[np.where(e.type == 'WORD')[0]].artifactMeanMS))
            af = np.concatenate((af, e[np.where(e.type == 'WORD')[0]].artifactFrac))

        print len(ams), len(amms), len(af)
    except:
        traceback.print_exc()

plt.figure(facecolor='white')
plt.subplot(311)
plt.hist(ams, range(0, 4401, 40), color='b', normed=True)
plt.axvline(3000, color='k', linewidth=2)
plt.title('Blink Onset Time')
plt.xlabel('MS after word onset')
plt.ylabel('Frequency')

plt.subplot(312)
plt.hist(amms, range(0, 4401, 40), color='r', normed=True)
plt.axvline(3000, color='k', linewidth=2)
plt.title('Avg Time of Blink Samples')
plt.xlabel('MS after word onset')
plt.ylabel('Frequency')

plt.subplot(313)
hist3 = plt.hist(af * 100, range(0, 101, 1), color='g', normed=True)
plt.title('Percent of Samps w/Blinks')
plt.xlabel('Percent')
plt.ylabel('Frequency')

plt.tight_layout()

pdf = matplotlib.backends.backend_pdf.PdfPages('/Users/jessepazdera/Desktop/Average Blink Stats (<20%).pdf')
for fig in xrange(1, plt.figure().number):
    pdf.savefig(fig)
pdf.close()