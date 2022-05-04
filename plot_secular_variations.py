import numpy as np
import os
from functions import cdf_times_to_np_date, compute_unique_times_list, compute_MF_list
import cdflib
import matplotlib.pyplot as plt


cdf_dir = 'donnees/cdf_files_basic_sync_functions_201'
data_file = os.path.join(cdf_dir, 'SW_OPER_VOBS_4M_2__20140301T000000_20210701T000000_0201_basic_sync_functions.cdf')
cdf_read = cdflib.CDF(data_file)
info = cdf_read.cdf_info()
zvars = info['zVariables']
alldata = {name: cdf_read.varget(name) for name in zvars}
print('>' ,alldata.keys())

# Choice of the date and component
times_list = alldata['Timestamp_SV']
print(times_list, np.shape(times_list))
for k in range(len(times_list)-1 , -1, -1):
    if np.isnan(times_list[k]):
        times_list = np.delete(times_list, k)
print(times_list, np.shape(times_list))
unique_times_list = compute_unique_times_list(cdf_times_to_np_date(times_list))
print(unique_times_list, np.shape(unique_times_list))

date = '2014-09' # Arbitrary choice of the date to plot
index_temps = np.argwhere((cdf_times_to_np_date(alldata['Timestamp_SV'])) == np.datetime64(date))[0][0]
print(index_temps)
N_VO = 300
choice = 0

# Creation of the different lists
latitudes = alldata['Latitude'][index_temps:index_temps+N_VO]
longitudes = alldata['Longitude'][index_temps:index_temps+N_VO]
MF = compute_MF_list(alldata['B_SV'][index_temps:index_temps+N_VO], choice)
print(np.shape(latitudes), np.shape(longitudes), np.shape(MF))

fig = plt.figure()
'''ax1 = fig.add_subplot(221, projection=ccrs.Robinson(central_longitude=0))'''
ax1=fig.add_subplot(221)
map_cf = ax1.scatter(longitudes, latitudes, c=MF, s=500, alpha=0.3,cmap='bwr')
'''ax1.coastlines()
ax1.set_global()'''
ax1.set_title("Geomagnetic field's component chosen")
plt.colorbar(map_cf, ax=ax1)
plt.show()

