# load metadata about event
import obspy
import yaml
import numpy as np
import copy
with open('metadata.yml') as f:
    metadata = yaml.load(f, Loader=yaml.FullLoader)

# location_visible_start
# location_visible_end
# location_impact
type_of_loc = 'location_impact'

ev_lat = metadata[type_of_loc][0]
ev_lon = metadata[type_of_loc][1]

time_in_video = obspy.UTCDateTime(metadata['event_time_visual_video'])

# read downloaded metadata
import os
inv = obspy.Inventory()
for fn in os.listdir('stations'):
    if '.xml' in fn:
        inv += obspy.read_inventory(f'stations/{fn}', format='STATIONXML')

# read downloaded data and write distances into .stats
st = obspy.read('waveforms/*HHZ*.mseed')

# REMOVE RESPONSE
pre_filt = (0.005, 0.006, 250.0, 400.0)
st.remove_response(inventory=inv, output='DISP', pre_filt=pre_filt)

st.trim(
    starttime=obspy.UTCDateTime('2020-04-06T13:30:00.0Z'),
    endtime=obspy.UTCDateTime('2020-04-06T13:50:00.0Z')
)

for tr in st:
    selected_station_meta = inv.select(
        network=tr.stats.network,
        station=tr.stats.station
        )
    
    if len(selected_station_meta) == 0:
        print(f'WARNING: NO metadata found for {tr.stats.station}')
        continue
    
    sta_lat = selected_station_meta[0][0].latitude
    sta_lon = selected_station_meta[0][0].longitude

    dist, az, baz = obspy.geodetics.gps2dist_azimuth(
        lat1=sta_lat, 
        lon1=sta_lon, 
        lat2=ev_lat, 
        lon2=ev_lon
        )
    tr.stats.distance = dist
    tr.coords = {}
    tr.coords['latitude'] = sta_lat
    tr.coords['longitude'] = sta_lon

max_dist = np.max([tr.stats.distance for tr in st])

# st.plot(type='section', outfile='test.png')

import pylab as plt
import matplotlib.dates as mdates
from sven_utils import suArray

from obspy.signal.filter import highpass
# %matplotlib inline

filter_pairs = (
    (1, 4),
    (4, 8),
    (8, 16),
    (16, 32),
    )

fig, axs = plt.subplots(1, len(filter_pairs), sharey=True, sharex=True, figsize=(8, 4))
fig.subplots_adjust(wspace=0, left=.075, right=.975, bottom=.05, top=.925)
for filter_pair, ax in zip(filter_pairs, axs.flatten()):
    for tr in st:
        if tr.stats.sampling_rate < 2*filter_pair[1]:
            print(f'INFO: Skipping {tr.stats.station}. Sampling rate too low for Filter Pair.')
            continue

        if tr.coords['latitude'] > ev_lat:
            print(f'INFO: Skipping {tr.stats.station}. South of event.')
            continue

        # create fresh copy to work with
        data = copy.deepcopy(tr.data)
        # filter
        data = suArray.bandpass(
            array=data,
            freqmin=filter_pair[0],
            freqmax=filter_pair[1],
            df=tr.stats.sampling_rate
            )

        # cut start and end to avoid filter artifacts
        data = data[100:-100]
        times = tr.times('matplotlib')[100:-100]

        # normalize and scale up for visibility
        data /= np.max(np.abs(data))
        scale_factor = 8
        data *= scale_factor

        ax.plot(
            times, 
            data + tr.stats.distance/1000,
            alpha=.25,
            lw=.1,
            color='k')
        
        ax.set_title(f'{filter_pair}Hz')

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(8) 
    
    if ax == axs.flatten()[0]:
        ax.set_ylabel(f'Distance from {type_of_loc}')
    
    # ax.axvline(time_in_video.matplotlib_date, ls=':', lw='.5', color='k')
    
    ax.plot([time_in_video.matplotlib_date, times[-1]], [0, 0.3*(tr.stats.endtime-time_in_video)], ls=':', c='#1f77b4', lw='.25', label='0.3km/s')
    ax.plot([time_in_video.matplotlib_date, times[-1]], [0, 3*(tr.stats.endtime-time_in_video)], ls='--', c='#1f77b4', lw='.25', label='3km/s')

    ax.set_ylim(0, max_dist/1000)


ax.legend()

ax.xaxis_date()
ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
fig.autofmt_xdate()

# ax.set_xlabel(f"UTC on {metadata['event_time_approx'].split('T')[0]}")
# ax.set_ylabel("Distance from visible endpoint [km]")
fig.savefig(f'out_{type_of_loc}.png', dpi=300)
plt.close(fig)