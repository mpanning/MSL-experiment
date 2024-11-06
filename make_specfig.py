"""
Script to do basic processing of data collected on the deck of testbed MSL and
reference seismometer on ground
"""

import numpy as np
from obspy import read
from obspy.core import UTCDateTime
from obspy.imaging.cm import obspy_sequential
from obspy.signal.tf_misfit import cwt
import math
from scipy import signal

import matplotlib.pyplot as plt

deckfile = 'Deck_taurus_0274_20171005_000000.seed'
reffile = 'Ref_taurus_2241_20171005_000000.seed'
calc_spec = True

tstart = UTCDateTime(2017, 10, 8, 8, 30)
length = 86400
tend = tstart + length

# Pull out the data segment for all components
st_deck = read(deckfile, starttime=tstart, endtime=tend)
st_ref = read(reffile, starttime=tstart, endtime=tend)

# components = ['Z', 'N', 'E']
components = ['Z']
for i, component in enumerate(components):
    # Extract component and merge segments into single trace
    st_deck_comp = st_deck.select(component=component)
    st_deck_comp.merge(method=1)

    st_ref_comp = st_ref.select(component=component)
    st_ref_comp.merge(method=1)

    # Combine into single stream, decimate and plot
    st_combined = st_ref_comp + st_deck_comp
    st_combined.decimate(5)
    st_combined.plot(outfile=component + '_raw.png')

    # Bandpass and plot
    f1 = 0.001
    f2 = 5.0
    st_bp = st_combined.copy()
    st_bp.detrend()
    st_bp.detrend('demean')
    st_bp.filter('bandpass', freqmin=f1, freqmax=f2, corners=2, zerophase=True)

    st_bp.plot(outfile=component + '_filtered.png')

    # Do tighter bandpass for display seismograms
    f1 = 0.05
    f2 = 1.0
    st_bp2 = st_combined.copy()
    st_bp2.detrend()
    st_bp2.detrend('demean')
    st_bp2.filter('bandpass', freqmin=f1, freqmax=f2, corners=2, zerophase=True)
    
    # Do spectrogram with scalogram
    if calc_spec:
        axnum = 210
        fig = plt.figure(figsize=(9,9))
        # Make 3 axes, with the top one for seismograms
        axs = [plt.subplot2grid((25, 1), (0, 0), rowspan=2)]
        axs.append(plt.subplot2grid((25, 1), (4, 0), rowspan=8, sharex=axs[0]))
        axs.append(plt.subplot2grid((25, 1), (14, 0), rowspan=8, sharex=axs[0]))
        axs.append(plt.subplot2grid((25, 1), (24, 0)))
        axnum = 0

        # Plot display seismograms
        dt = st_bp2[0].stats.delta
        npts = st_bp2[0].stats.npts
        t = np.linspace(0, dt*npts, npts)
        axs[0].plot(t, st_bp2[0].data, color='k')
        axs[0].plot(t, st_bp2[1].data, color='b')
        plt.sca(axs[0])
        plt.ylim([-2500., 2500.])
        
        for stn in ['REF', 'DECK']:
            print('Scalogram working on station: ' + stn + ' and component: '
                  + component)
            tr_comb = st_bp.select(station=stn)[0]

            npts = tr_comb.stats.npts
            dt = tr_comb.stats.delta
            t = np.linspace(0, dt * npts, npts)
            f_min = 1./150.
            f_max = 5.

            scalogram = cwt(tr_comb.data, dt, 6, f_min, f_max, nf=100)
   
            axnum += 1
            ax = axs[axnum]

            x, y = np.meshgrid(t,
                               np.logspace(np.log10(f_min), 
                                           np.log10(f_max), 
                                           scalogram.shape[0]))

            m = ax.pcolormesh(x, y, np.log10((scalogram)), 
                              cmap='viridis',
                              vmin=0, vmax=3)
            ax.set_title(stn)
            ax.set_ylabel("Frequency [Hz]")
            ax.set_yscale('log')
    
            ax.set_ylim(f_min, f_max)
            ax.set_xlim(0,length)
            if axnum == 2:
                ax.set_xlabel('time / sec')


        cbar = plt.colorbar(mappable=m, cax=axs[3], orientation='horizontal',
                            ticks=[0, 1, 2, 3])
        cbar.ax.set_xticklabels([r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$'])
        fig.savefig(component + '_specfig.png')

