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
import csv

import matplotlib as mpl
import matplotlib.pyplot as plt

deckfile = 'Deck_taurus_0274_20171005_000000.seed'
reffile = 'Ref_taurus_2241_20171005_000000.seed'
deck_temp_file = 'Deck_taurus_0274_SOH_20171006_000000.csv'
ref_temp_file = 'Ref_taurus_2241_SOH_20171006_000000.csv'
outfile = 'temp_motion_comp.png'

tstart = UTCDateTime(2017, 10, 6, 19)
length = 68.*3600
component = 'Z'

# First pull out the seismic data
tend = tstart + length
st_deck = read(deckfile, starttime=tstart, endtime=tend)
st_ref = read(reffile, starttime=tstart, endtime=tend)

st_deck_comp = st_deck.select(component=component)
st_deck_comp.merge(method=1, fill_value=0.)

st_ref_comp = st_ref.select(component=component)
st_ref_comp.merge(method=1, fill_value=0)


# Combine into single stream and plot
st_combined = st_ref_comp + st_deck_comp
st_combined.plot(outfile=component + '_long_raw.png')

# Bandpass and plot
f1 = 1.0
f2 = 50.0
st_bp = st_combined.copy()
st_bp.detrend()
st_bp.detrend('demean')
st_bp.filter('bandpass', freqmin=f1, freqmax=f2, corners=2, zerophase=True)

st_bp.plot(outfile=component + '_long_filtered.png')

# Read in temperature data
deck_temp_t = []
deck_temp = []
with open(deck_temp_file) as csvfile:
    reader = csv.reader(csvfile)
    next(reader, None) # skip the header row
    for row in reader:
        deck_temp_t.append(UTCDateTime(row[1])-tstart)
        deck_temp.append(float(row[2]))
deck_temp_t = np.array(deck_temp_t)
deck_temp = np.array(deck_temp)
        
ref_temp_t = []
ref_temp = []
with open(ref_temp_file) as csvfile:
    reader = csv.reader(csvfile)
    next(reader, None) # skip the header row
    for row in reader:
        ref_temp_t.append(UTCDateTime(row[1])-tstart)
        ref_temp.append(float(row[2]))
ref_temp_t = np.array(deck_temp_t)
ref_temp = np.array(deck_temp)

# Make figure
fig = plt.figure(figsize=(10, 10))
ax1 = plt.subplot2grid((17,1), (0,0), rowspan=5)
ax2 = plt.subplot2grid((17,1), (6,0), rowspan=5)
ax3 = plt.subplot2grid((17,1), (12,0), rowspan=5)

# Plot up temperature data in ax1
ax1.plot(ref_temp_t/3600., ref_temp, color='k', label='Reference')
ax1.plot(deck_temp_t/3600., deck_temp+1.0, color='b', label='Deck')
plt.sca(ax1)
plt.ylabel(r'Temperature ($^{\circ}$C)')
plt.xlim([0, length/3600.])
plt.ylim([22, 28])
plt.legend(loc=4)

# Create time array for seismic data and plot raw in ax2 and bp in ax3
dt = st_combined[1].stats.delta
npts = st_combined[1].stats.npts
t = np.linspace(0, dt*npts, npts)/3600.
ax2.plot(t,st_combined[1].data, color='b')

dt = st_combined[0].stats.delta
npts = st_combined[0].stats.npts
t = np.linspace(0, dt*npts, npts)/3600.
ax2.plot(t,st_combined[0].data, color='k')

ax2.ticklabel_format(style='sci', axis='y')
ax2.yaxis.major.formatter.set_powerlimits((0,0))
plt.sca(ax2)
plt.ylabel('Counts')
plt.xlim([0, length/3600.])
plt.ylim([-1e5, 1e5])

dt = st_bp[1].stats.delta
npts = st_bp[1].stats.npts
t = np.linspace(0, dt*npts, npts)/3600.
ax3.plot(t,st_bp[1].data, color='b')

dt = st_bp[0].stats.delta
npts = st_bp[0].stats.npts
t = np.linspace(0, dt*npts, npts)/3600.
ax3.plot(t,st_bp[0].data, color='k')

ax3.ticklabel_format(style='sci', axis='y')
ax3.yaxis.major.formatter.set_powerlimits((0,0))
plt.sca(ax3)
plt.ylabel('Counts')
plt.xlabel('Time (hours)')
plt.xlim([0, length/3600.])
plt.ylim([-1e4, 1e4])

# Output outfile
mpl.rcParams['agg.path.chunksize'] = 10000
fig.savefig(outfile)

