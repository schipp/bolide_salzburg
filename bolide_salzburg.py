# load metadata about event
import obspy
import yaml
import numpy as np
import copy
import pylab as plt
import matplotlib.dates as mdates
from bandpass import bandpass

# valid locs: location_impact, location_visible_end, location_visible_start
type_of_loc = 'location_visible_start'

# filter pairs to look at
filter_pairs = (
    (1, 10),
    )

# read from event metadata
with open('metadata.yml') as f:
    metadata = yaml.load(f, Loader=yaml.FullLoader)
ev_lat = metadata[type_of_loc][0]
ev_lon = metadata[type_of_loc][1]
time_in_video = obspy.UTCDateTime(metadata['event_time_visual_video'])

# read downloaded metadata
import os
inv = obspy.Inventory()
for fn in os.listdir('stations'):
    if '.xml' in fn:
        inv += obspy.read_inventory(f'stations/{fn}', format='STATIONXML')

# read downloaded data
st = obspy.read('waveforms/deconv/HHZ.mseed')

# # REMOVE RESPONSE
# pre_filt = (0.005, 0.006, 250.0, 400.0)
# st.remove_response(inventory=inv, output='DISP', pre_filt=pre_filt)

# trim to reasonable window
st.trim(
    starttime=obspy.UTCDateTime('2020-04-06T13:32:00.0Z'),
    endtime=obspy.UTCDateTime('2020-04-06T13:55:00.0Z')
)

# compute distances
dists = []
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

    dists.append(dist)

max_dist = np.max([tr.stats.distance for tr in st])

st_sorted = [y for x, y in sorted(zip(dists, st))]

fig, axs = plt.subplots(len(st_sorted), 1, sharey=True, figsize=(6, 6))
fig.subplots_adjust(wspace=0, hspace=0, left=.075, right=.975, bottom=0, top=.95)
# axs = [axs]

filter_pair = filter_pairs[0]

stations_used_by_ZAMG = (
    'LESA',
    'SOSA',
    'KBA',
    'A038A',
    'ABTA'
)

for tr, ax in zip(st_sorted[::-1], axs):
    if tr.stats.sampling_rate < 2*filter_pair[1]:
        print(f'INFO: Skipping {tr.stats.station}. Sampling rate too low for Filter Pair.')
        continue

    # create fresh copy to work with
    data = copy.deepcopy(tr.data)

    # filter
    data = bandpass(
        array=data,
        freqmin=filter_pair[0],
        freqmax=filter_pair[1],
        df=tr.stats.sampling_rate
        )

    # cut start and end to avoid filter artifacts
    data = data[500:-500]
    times = tr.times('matplotlib')[500:-500]

    # normalize and scale up for visibility
    data /= np.max(np.abs(data))
    # scale_factor = 8
    # data *= scale_factor
    col = 'k'
    if tr.stats.station in stations_used_by_ZAMG:
        col = '#1f77b4'

    if tr.stats.network == 'Z3':
        col = 'grey'
    
    ax.plot(
        times, 
        data,
        alpha=1,
        lw=.1,
        color=col)
    
    # ax.set_title(f'{filter_pair}Hz')
    ax.set_frame_on(False)
    if ax != axs[-1]:
        ax.set_xticks([])
    ax.set_yticks([])

    # ax.axvline((time_in_video + tr.stats.distance/3000).matplotlib_date, ls='-', c='#2ca02c', lw='1')
    # ax.axvline((time_in_video-1*60 + tr.stats.distance/300).matplotlib_date, ls='-', c='#1f77b4', lw='1')
    ax.scatter((time_in_video + tr.stats.distance/330).matplotlib_date, 0, s=2, c='#ff7f0e')
    # ax.axvline((time_in_video+1*60 + tr.stats.distance/500).matplotlib_date, ls='-', c='#ff7f0e', lw='1')

    ax.set_xlim(times[0], times[-1])

    ax.text(0.01, 0.95, f'{tr.stats.distance/1000:0.1f}km', fontsize=4, va='top', ha='left', transform=ax.transAxes)


for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(4) 

if ax == axs[0]:
    ax.set_ylabel(f'Distance from {type_of_loc}')

# ax.plot([time_in_video.matplotlib_date, times[-1]], [0, 0.3*(tr.stats.endtime-time_in_video)], ls=':', c='#1f77b4', lw='.25', label='0.3km/s')
# ax.plot([time_in_video.matplotlib_date, times[-1]], [0, 3*(tr.stats.endtime-time_in_video)], ls='--', c='#1f77b4', lw='.25', label='3km/s')

# ax.set_ylim(0, max_dist/1000)


ax.xaxis_date()
ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
fig.autofmt_xdate()

fig.savefig(f'out_{type_of_loc}.png', dpi=300)
plt.close(fig)