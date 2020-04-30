import yaml
with open('metadata.yml') as f:
    metadata = yaml.load(f, Loader=yaml.FullLoader)

# download data near trajectory (use visible end location for now)
import obspy
from obspy.clients.fdsn import Client
client = Client('ORFEUS')
starttime = obspy.UTCDateTime(metadata['event_time_approx']) - 30*60
endtime = obspy.UTCDateTime(metadata['event_time_approx']) + 30*60
from obspy.clients.fdsn.mass_downloader import CircularDomain, \
    Restrictions, MassDownloader

client.set_eida_token('eidatoken.txt')

domain = CircularDomain(
    latitude=metadata['location_visible_end'][0],
    longitude=metadata['location_visible_end'][1],
    minradius=None,
    maxradius=2
    )

restrictions = Restrictions(
    starttime=starttime,
    endtime=endtime,
    reject_channels_with_gaps=True,
    minimum_length=0.95,
    minimum_interstation_distance_in_m=5E3,
    channel_priorities=["HH[Z]", "BH[Z]"],
    # exclude_networks=['Z3']
    )

# No specified providers will result in all known ones being queried.
mdl = MassDownloader([client])
# The data will be downloaded to the ``./waveforms/`` and ``./stations/``
mdl.download(domain, restrictions, mseed_storage="waveforms",
             stationxml_storage="stations")
