import os
import cdflib
from functions import cdf_times_to_np_date, compute_unique_times_list

# Open file
cdf_dir = 'donnees/cdf_files_basic_sync_functions_201'
swarm_file = os.path.join(cdf_dir, 'SW_OPER_VOBS_4M_2__20140301T000000_20210701T000000_0201_basic_sync_functions.cdf')
ground_file = os.path.join(cdf_dir, 'GObs_4M_19970301T000000_20211101T000000_0103_basic_sync_functions.cdf')
cdf_read = cdflib.CDF(swarm_file)
info = cdf_read.cdf_info()
print('info', info)
zvars = info['zVariables']
alldata = {name: cdf_read.varget(name) for name in zvars}
print('observations keys: ', alldata.keys())

# Plot the localisations
times = alldata['Timestamp']
unique_times = cdf_times_to_np_date(compute_unique_times_list(times))
all_latitudes = alldata['Latitude']
all_longitudes = alldata['Longitude']
N_VO = int(len(all_latitudes) / len(unique_times))



