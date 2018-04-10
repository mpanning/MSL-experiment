# MSL-experiment
Archive for data and tools associated with study of the MSL seismic response

```
Deck_taurus_0274_20171005_000000.seed
Ref_taurus_2241_20171005_000000.seed
```

These are the miniseed format files containing the seismic data for the seismometer on the deck and on the ground below MSL, respectively.

```
Deck_taurus_0274_SOH_20171006_000000.csv
Ref_taurus_2241_SOH_20171006_000000.csv
```

These contain the temperature and state of health data for the seismometers in csv format.

Python scripts for processing and plotting the data to reproduce the figures in the paper are also included.  These have been verified with Python 2.7, and require the package ObsPy 1.0.3 (some errors show up in the spectrogram plotting when using the more recent ObsPy version 1.1.0)

```
temp_fig.py
make_specfig.py
find_events.py
ppsd_daynight_fig.py
```

These are example python codes to use and process the data.  They are the source code for making figures 3, 4, 5, and 7, respectively.
