"""
Script to make figure comparing noisy and quiet period coherence and 
transfer function
"""

import numpy as np
from obspy import read
from obspy.core import UTCDateTime
from obspy.imaging.cm import obspy_sequential
from obspy.signal.tf_misfit import cwt
import math
from scipy import signal
import sys

# import matplotlib as mpl
import matplotlib.pyplot as plt

deckfile = 'Deck_taurus_0274_20171005_000000.seed'
reffile = 'Ref_taurus_2241_20171005_000000.seed'
outfile = 'tf_coherence_comp.png'

tstarts = [UTCDateTime(2017, 10, 7, 15), UTCDateTime(2017, 10, 7, 15),
           UTCDateTime(2017, 10, 8, 3)]
lengths = [int(86400*1.0), int(86400*0.5), int(86400*0.5)]
colors = ['black', 'green', 'red']
component = 'Z'

# First pull out a long data segment with an extra hour on either end
tstart = tstarts[0] - 3600
tend = tstart + lengths[0] + 7200
st_deck = read(deckfile, starttime=tstart, endtime=tend)
st_ref = read(reffile, starttime=tstart, endtime=tend)

st_deck_comp = st_deck.select(component=component)
st_deck_comp.merge(method=1)

st_ref_comp = st_ref.select(component=component)
st_ref_comp.merge(method=1)


# Combine into single stream, decimate and plot
st_combined = st_ref_comp + st_deck_comp
st_combined.decimate(5)
st_combined.plot(outfile=component + '_long_raw.png')
# st_combined.spectrogram(outfile='rawspec.png')

# Bandpass and plot
f1 = 0.001
f2 = 5.0
st_bp = st_combined.copy()
st_bp.detrend()
st_bp.detrend('demean')
st_bp.filter('bandpass', freqmin=f1, freqmax=f2, corners=2, zerophase=True)

st_bp.plot(outfile=component + '_long_filtered.png')

fig = plt.figure(figsize=(10, 10))
ax1 = plt.subplot2grid((17,1), (0,0), rowspan=5)
ax2 = plt.subplot2grid((17,1), (6,0), rowspan=5)
ax3 = plt.subplot2grid((17,1), (12,0), rowspan=5)

dt = st_bp[0].stats.delta
npts = st_bp[0].stats.npts
t = np.linspace(0, dt*npts, npts)

ax1.plot(t,st_bp[0].data, color='k')
ax1.plot(t,st_bp[1].data, color='b')

# Make boxes for smaller segments
t1 = tstarts[1] - tstart
t2 = t1 + lengths[1]
t3 = tstarts[2] - tstart
t4 = t3 + lengths[2]

ax1.axvspan(t1, t2, alpha=0.5, color=colors[1])
ax1.axvspan(t3, t4, alpha=0.5, color=colors[2])

plt.sca(ax1)
plt.xlim([0, dt*npts])
plt.xlabel("Time (s)")
plt.ylabel("Counts")

# Loop over windows and compute coherence and transfer function
for i, tstart in enumerate(tstarts):
    tend = tstart + lengths[i]
    st_deck = read(deckfile, starttime=tstart, endtime=tend)
    st_ref = read(reffile, starttime=tstart, endtime=tend)

    st_deck_comp = st_deck.select(component=component)
    st_deck_comp.merge(method=1)

    st_ref_comp = st_ref.select(component=component)
    st_ref_comp.merge(method=1)


    # Combine into single stream and decimate
    st_combined = st_ref_comp + st_deck_comp
    st_combined.decimate(5)

    # Bandpass
    st_bp = st_combined.copy()
    st_bp.detrend()
    st_bp.detrend('demean')
    st_bp.filter('bandpass', freqmin=f1, freqmax=f2, corners=2, zerophase=True)

    fs = st_bp[0].stats['sampling_rate']
    refdata = st_bp[0].data
    deckdata = st_bp[1].data

    freq, Cxy = signal.coherence(refdata, deckdata, fs, nperseg=8192)
    plt.sca(ax2)
    plt.loglog(freq, Cxy, colors[i])
    plt.ylim([0.001, 1.0])
    plt.xlim([0.005, 10.0])
    # plt.xlabel("Frequency (Hz)")
    plt.ylabel("Coherence")
    # plt.plot([0.001, 10.0],[0.7, 0.7],'r--')

    # Now do transfer function
    lamda = 1.0e6
    dt = st_bp[0].stats['delta']
    seglength = int(1000*st_bp[0].stats['sampling_rate'])

    i1 = 0
    i2 = i1 + seglength
    tf_list = []
    while (i2 < len(st_bp[0].data)):
        refdata = st_bp[0].data[i1:i2]
        reffft = np.fft.fft(refdata)
        deckdata = st_bp[1].data[i1:i2]
        deckfft = np.fft.fft(deckdata)
        freq = np.fft.fftfreq(len(reffft), dt)

        # Regularize the spectral division
        transfer = (np.conjugate(reffft)*deckfft/
                    (np.conjugate(reffft)*reffft + lamda))
        print(len(transfer))
        tf_list.append(transfer)
        i1 += seglength
        i2 += seglength

    num_tf = len(tf_list)
    print('Finished calculating ' + str(num_tf) + ' transfer functions')
    tf_array = np.array(tf_list)
    tf_average = np.average(tf_array, axis=0)

    plt.sca(ax3)
    plt.plot(freq[freq>0], np.abs(tf_average[freq>0]), color=colors[i])
    plt.xlim([0.005, 10.0])
    plt.xscale('log')
    plt.ylim([0.01, 100.0])
    plt.yscale('log')
    plt.ylabel('Transfer magnitude')
    plt.xlabel('Frequency (Hz)')
    
    

fig.savefig(outfile)

