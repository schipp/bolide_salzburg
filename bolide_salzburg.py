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
    (8, 15),
    )

v_sound = 330
zoom_window = 45

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
    endtime=obspy.UTCDateTime('2020-04-06T13:50:00.0Z')
)

# compute distances
dists = []
sta_lats = []
pythagoras_dists = []
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

    # dist from visible event locations
    dist, az, baz = obspy.geodetics.gps2dist_azimuth(
        lat1=sta_lat, 
        lon1=sta_lon, 
        lat2=ev_lat, 
        lon2=ev_lon
        )
    
    # dist from trajectory line
    p1 = np.array([metadata['location_visible_start'][1], metadata['location_visible_start'][0], obspy.geodetics.kilometers2degrees(80)])
    p2 = np.array([metadata['location_impact'][1], metadata['location_impact'][0], 0])
    # p2 = np.array([metadata['location_visible_end'][1], metadata['location_visible_end'][0], obspy.geodetics.kilometers2degrees(40)])
    p3 = np.array([sta_lon, sta_lat, 0])
    d = np.linalg.norm(np.cross(p2-p1, p1-p3))/np.linalg.norm(p2-p1)
    dist = obspy.geodetics.degrees2kilometers(d) * 1000

    tr.stats.distance = dist
    tr.stats.coords = {}
    tr.stats.coords['latitude'] = sta_lat
    tr.stats.coords['longitude'] = sta_lon

    sta_lats.append(sta_lat)
    dists.append(dist)
    # pythagoras_dists.append(pythagoras_dist)

max_dist = np.max([tr.stats.distance for tr in st])

st_sorted = [y for x, y in sorted(zip(dists, st))]
# st_sorted = [y for x, y in sorted(zip(sta_lats, st))]

fig, axs = plt.subplots(len(st_sorted), 1, sharey=True, figsize=(6, 6))
fig.subplots_adjust(wspace=0, hspace=0, left=.075, right=.975, bottom=0, top=.95)
# axs = [axs]

fig_zoom, axs_zoom = plt.subplots(len(st_sorted), 1, figsize=(6, 6))
fig_zoom.subplots_adjust(wspace=0, hspace=0, left=.075, right=.975, bottom=0, top=.95)

filter_pair = filter_pairs[0]

stations_used_by_ZAMG = (
    'LESA',
    'SOSA',
    'KBA',
    'A038A',
    'ABTA'
)

for tr, ax, ax_zoom in zip(st_sorted[::-1], axs, axs_zoom):
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

    # # AGC
    # (data, D, E) = tf_agc(data, tr.stats.sampling_rate, t_scale=30, f_scale=1.0, causal_tracking=True)

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

    zoom_time_start = time_in_video + tr.stats.distance/v_sound - zoom_window
    zoom_time_end = time_in_video + tr.stats.distance/v_sound + zoom_window

    # ax.axvline((time_in_video + tr.stats.distance/3000).matplotlib_date, ls='-', c='#2ca02c', lw='1')
    # ax.axvline((time_in_video-1*60 + tr.stats.distance/300).matplotlib_date, ls='-', c='#1f77b4', lw='1')
    # ax.scatter((time_in_video + tr.stats.distance/v_sound).matplotlib_date, 0, s=2, c='#ff7f0e')
    
    ax.axvline((time_in_video + tr.stats.distance/v_sound).matplotlib_date, lw=.75, c='#ff7f0e', zorder=0)
    # ax.axvline(zoom_time_start.matplotlib_date, lw=.5, c='#ff7f0e', zorder=0)
    # ax.axvline(zoom_time_end.matplotlib_date, lw=.5, c='#ff7f0e', zorder=0)
    
    # la_dist, _, _ = obspy.geodetics.gps2dist_azimuth(lat1=47.42, lon1=13, lat2=tr.stats.coords.latitude, lon2=13)
    # ax.scatter((time_in_video+2*60 + la_dist/v_sound).matplotlib_date, 0, s=2, c='#ff7f0e')
    # ax.axvline((time_in_video+1*60 + tr.stats.distance/500).matplotlib_date, ls='-', c='#ff7f0e', lw='1')

    ax.set_xlim(times[0], times[-1])

    ax.text(0.01, 0.95, f'{tr.stats.distance/1000:0.1f}km', fontsize=4, va='top', ha='right', transform=ax.transAxes)
    ax.text(0.01, 0.95, f', {tr.stats.coords.longitude:0.1f}°E', fontsize=4, va='top', ha='left', transform=ax.transAxes)
    # ax.text(0.01, 0.95, f'{tr.stats.coords.latitude}°', fontsize=4, va='top', ha='left', transform=ax.transAxes)

    # ZOOM IN VIEW
    ax_zoom.plot(
        times, 
        data,
        alpha=1,
        lw=.25,
        color=col)
    
    # ax.set_title(f'{filter_pair}Hz')
    ax_zoom.set_frame_on(False)
    if ax_zoom != axs[-1]:
        ax_zoom.set_xticks([])
    ax_zoom.set_yticks([])

    # ax.axvline((time_in_video + tr.stats.distance/3000).matplotlib_date, ls='-', c='#2ca02c', lw='1')
    # ax.axvline((time_in_video-1*60 + tr.stats.distance/300).matplotlib_date, ls='-', c='#1f77b4', lw='1')
    # ax.scatter((time_in_video + tr.stats.distance/v_sound).matplotlib_date, 0, s=2, c='#ff7f0e')
    ax_zoom.axvline((time_in_video + tr.stats.distance/v_sound).matplotlib_date, lw=.5, c='#ff7f0e', zorder=0)
    # la_dist, _, _ = obspy.geodetics.gps2dist_azimuth(lat1=47.42, lon1=13, lat2=tr.stats.coords.latitude, lon2=13)
    # ax.scatter((time_in_video+2*60 + la_dist/v_sound).matplotlib_date, 0, s=2, c='#ff7f0e')
    # ax.axvline((time_in_video+1*60 + tr.stats.distance/500).matplotlib_date, ls='-', c='#ff7f0e', lw='1')

    ax_zoom.set_xlim(zoom_time_start.matplotlib_date, zoom_time_end.matplotlib_date)
    # print()
    data_in_xlims = data[np.argmin(np.abs(zoom_time_start.matplotlib_date - times)):np.argmin(np.abs(zoom_time_end.matplotlib_date - times))]
    ax_zoom.set_ylim(-1.1*np.max(np.abs(data_in_xlims)), 1.1*np.max(np.abs(data_in_xlims)))

    ax_zoom.text(0.01, 0.95, f'{tr.stats.distance/1000:0.1f}km', fontsize=4, va='top', ha='right', transform=ax_zoom.transAxes)
    ax_zoom.text(0.01, 0.95, f', {tr.stats.coords.longitude:0.1f}°E', fontsize=4, va='top', ha='left', transform=ax_zoom.transAxes)
    
    if ax_zoom == axs_zoom[-1]:
        ax_zoom.set_xticks([])


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


ax_zoom.xaxis_date()
ax_zoom.xaxis.set_major_locator(mdates.SecondLocator(interval=10))
ax_zoom.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))
fig_zoom.autofmt_xdate()

for tick in ax_zoom.xaxis.get_major_ticks():
    tick.label.set_fontsize(4) 

fig.savefig(f'out_{type_of_loc}.png', dpi=300)
fig_zoom.savefig(f'out_zoom_{type_of_loc}.png', dpi=300)
plt.close(fig)
plt.close(fig_zoom)