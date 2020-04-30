import obspy
import os
inv = obspy.Inventory()
for fn in os.listdir('stations'):
    if '.xml' in fn:
        inv += obspy.read_inventory(f'stations/{fn}', format='STATIONXML')

# read downloaded data
st = obspy.read('waveforms/*HHZ*.mseed')

# REMOVE RESPONSE
pre_filt = (0.005, 0.006, 250.0, 400.0)
st.remove_response(inventory=inv, output='DISP', pre_filt=pre_filt)

st.write('waveforms/deconv/HHZ.mseed')