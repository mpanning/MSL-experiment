import numpy as np
import matplotlib.pyplot as plt
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

client = Client('IRIS')

starttime = UTCDateTime(2017, 10, 8, 8, 30)
length = 86400
endtime = starttime + length
minmag = 6.0
maxmag = 9.0
localmin = 4.5
localrad = 30.0
lat = 34.2001
lon = -118.176


big_events = client.get_events(starttime=starttime, endtime=endtime,
                               minmagnitude=minmag, maxmagnitude=maxmag)
print(big_events)

local_events = client.get_events(starttime=starttime, endtime=endtime,
                                 minmagnitude=localmin, maxmagnitude=maxmag,
                                 latitude=lat, longitude=lon,
                                 maxradius=localrad)
print(local_events)

all_events = big_events + local_events

# Get location of PASC which is closest station to Mars Yard
inventory = client.get_stations(network='CI', station='PASC',
                                starttime=starttime)

# Plot
fig = inventory.plot(label=False, show=False, projection='global')
all_events.plot(color='date', outfile='events.png', fig=fig)
