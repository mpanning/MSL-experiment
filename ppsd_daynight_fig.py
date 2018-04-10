"""
Script to plot up ppsd figures for reference and deck data
"""

import numpy as np
from obspy import read, UTCDateTime
from obspy.signal import PPSD
# from obspy.io.xseed import Parser
from obspy.signal.spectral_estimation import get_nhnm, get_nlnm
import matplotlib as mpl
import matplotlib.pyplot as plt

deckfile = 'Deck_taurus_0274_20171005_000000.seed'
reffile = 'Ref_taurus_2241_20171005_000000.seed'


# Trying with a paz dictionary (not sure if this is right, but the numbers
# look right)
TrillC = {'gain': 749.1,
          'poles': [complex(-3.691000e-02,3.712000e-02),
                    complex(-3.691000e-02,-3.712000e-02),
                    complex(-3.712000e+02,0.0),
                    complex(-3.739000e+02,4.755000e+02),
                    complex(-3.739000e+02,-4.755000e+02),
                    complex(-5.884000e+02,1.508000e+03),
                    complex(-5.884000e+02,-1.508000e+03)],
          # 'sensitivity': 4.91313873354e15,
          'sensitivity': 4.344928E+17,
          'zeros': [0, 0, -4.341E+02],
          'normalization_frequency': 1.0
}

dayStartTimes = [UTCDateTime(2017, 10, 6, 20), UTCDateTime(2017, 10, 7, 15),
                 UTCDateTime(2017, 10, 8, 15), UTCDateTime(2017, 10, 9, 15)] 
dayLengths = [7.*3600, 12.*3600, 12.*3600, 3.*3600]
nightStartTimes = [UTCDateTime(2017, 10, 7, 3), UTCDateTime(2017, 10, 8, 3),
                   UTCDateTime(2017, 10, 9, 3)]
nightLengths = [12.*3600, 12.*3600, 12.*3600]

# Pull out the data segments for all components for day and night time
day_deck_sts = []
day_ref_sts = []
for i, starttime in enumerate(dayStartTimes):
    day_deck_sts.append(read(deckfile, starttime=starttime,
                             endtime=starttime + dayLengths[i]))
    day_ref_sts.append(read(reffile, starttime=starttime,
                            endtime=starttime + dayLengths[i]))
night_deck_sts = []
night_ref_sts = []
for i, starttime in enumerate(nightStartTimes):
    night_deck_sts.append(read(deckfile, starttime=starttime,
                               endtime=starttime + nightLengths[i]))
    night_ref_sts.append(read(reffile, starttime=starttime,
                              endtime=starttime + nightLengths[i]))

components = ['Z', 'N', 'E']

day_deck_mean_psd = []
night_deck_mean_psd = []
day_ref_mean_psd = []
night_ref_mean_psd = []
day_deck_mean_pd = []
night_deck_mean_pd = []
day_ref_mean_pd = []
night_ref_mean_pd = []
for i, component in enumerate(components):
    day_deck_comp_sts = []
    day_ref_comp_sts = []
    night_deck_comp_sts = []
    night_ref_comp_sts = []
    for j, st in enumerate(day_deck_sts):
        day_deck_comp_sts.append(st.select(component=component))
        day_deck_comp_sts[j].merge(method=1)
    for j, st in enumerate(day_ref_sts):
        day_ref_comp_sts.append(st.select(component=component))
        day_ref_comp_sts[j].merge(method=1)
    for j, st in enumerate(night_deck_sts):
        night_deck_comp_sts.append(st.select(component=component))
        night_deck_comp_sts[j].merge(method=1)
    for j, st in enumerate(night_ref_sts):
        night_ref_comp_sts.append(st.select(component=component))
        night_ref_comp_sts[j].merge(method=1)

    day_ppsd_deck = PPSD(day_deck_comp_sts[0][0].stats, TrillC,
                         ppsd_length=600.0)
    for st in day_deck_comp_sts:
        day_ppsd_deck.add(st)
    # plotfile = 'day_deck_ppsd_' + component + '.png'
    # day_ppsd_deck.plot(plotfile, show_coverage=False)
    (meanpd, meanpsd) = day_ppsd_deck.get_mean()
    day_deck_mean_pd.append(meanpd)
    day_deck_mean_psd.append(meanpsd)

    night_ppsd_deck = PPSD(night_deck_comp_sts[0][0].stats, TrillC,
                           ppsd_length=600.0)
    for st in night_deck_comp_sts:
        night_ppsd_deck.add(st)
    # plotfile = 'night_deck_ppsd_' + component + '.png'
    # night_ppsd_deck.plot(plotfile, show_coverage=False)
    (meanpd, meanpsd) = night_ppsd_deck.get_mean()
    night_deck_mean_pd.append(meanpd)
    night_deck_mean_psd.append(meanpsd)

    day_ppsd_ref = PPSD(day_ref_comp_sts[0][0].stats, TrillC,
                        ppsd_length=600.0)
    for st in day_ref_comp_sts:
        day_ppsd_ref.add(st)
    # plotfile = 'day_ref_ppsd_' + component + '.png'
    # day_ppsd_ref.plot(plotfile, show_coverage=False)
    (meanpd, meanpsd) = day_ppsd_ref.get_mean()
    day_ref_mean_pd.append(meanpd)
    day_ref_mean_psd.append(meanpsd)

    night_ppsd_ref = PPSD(night_ref_comp_sts[0][0].stats, TrillC,
                          ppsd_length=600.0)
    for st in night_ref_comp_sts:
        night_ppsd_ref.add(st)
    # plotfile = 'night_ref_ppsd_' + component + '.png'
    # night_ppsd_ref.plot(plotfile, show_coverage=False)
    (meanpd, meanpsd) = night_ppsd_ref.get_mean()
    night_ref_mean_pd.append(meanpd)
    night_ref_mean_psd.append(meanpsd)

# Plot up Earth noise models and mean psd values
nlnm_pd, nlnm = get_nlnm()
nhnm_pd, nhnm = get_nhnm()

fig = plt.figure(figsize=(6, 9))

# First the vertical component
ax = plt.subplot2grid((2, 8), (0, 1), colspan=7)

ax.semilogx(nhnm_pd, nhnm, linewidth=2, color='darkgrey',
            label='Earth noise')
ax.semilogx(nlnm_pd, nlnm, linewidth=2, color='darkgrey')

ax.semilogx(day_ref_mean_pd[0], day_ref_mean_psd[0], color='red')
ax.semilogx(day_deck_mean_pd[0], day_deck_mean_psd[0], color='red', ls='--')
ax.semilogx(night_ref_mean_pd[0], night_ref_mean_psd[0], color='green')
ax.semilogx(night_deck_mean_pd[0], night_deck_mean_psd[0], color='green',
            ls='--')
plt.xlim([0.01, 100])
plt.ylim([-200, -50])
plt.ylabel(r'Power ($m^2/s^4/Hz$) (dB)')

# Average the horizontals
ax = plt.subplot2grid((2, 8), (1, 1), colspan=7)


ax.semilogx(nhnm_pd, nhnm, linewidth=2, color='darkgrey',
            label='Earth noise')
ax.semilogx(nlnm_pd, nlnm, linewidth=2, color='darkgrey')

ax.semilogx(day_ref_mean_pd[0], 0.5*(day_ref_mean_psd[1] + day_ref_mean_psd[2]),
            color='red')
ax.semilogx(day_deck_mean_pd[0],
            0.5*(day_deck_mean_psd[1] + day_deck_mean_psd[2]), color='red',
            ls='--')
ax.semilogx(night_ref_mean_pd[0],
            0.5*(night_ref_mean_psd[1] + night_ref_mean_psd[2]), color='green')
ax.semilogx(night_deck_mean_pd[0],
            0.5*(night_deck_mean_psd[1] + night_deck_mean_psd[2]),
            color='green', ls='--')
plt.xlim([0.01, 100])
plt.ylim([-200, -50])
plt.xlabel('Period (s)')
plt.ylabel(r'Power ($m^2/s^4/Hz$) (dB)')


plt.savefig('Summary_psd.png')
