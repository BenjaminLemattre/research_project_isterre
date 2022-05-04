import os
from functions import *


#  dates ['2014-05', '2014-09', '2015-01', '2015-05', '2015-09, '2016-01', '2016-05', '2016-09', '2017-01', '2017-05' \\
#  '2017-09', 2018-01', '2018-05', '2018-09', '2019-01', '2019-05', '2019-09', '2020-01', '2020-05', '2020-09', '2021-01', '2021-05']

#  Directories
N_VO = 300
L_max = 13
date = '2017-05'
cm_prior = np.loadtxt('prior_cov_matrix.txt')
cdf_dir = 'donnees/cdf_files_basic_sync_functions_201'
swarm_file = os.path.join(cdf_dir, 'SW_OPER_VOBS_4M_2__20140301T000000_20210701T000000_0201_basic_sync_functions.cdf')
cdf_read = cdflib.CDF(swarm_file)
info = cdf_read.cdf_info()
zvars = info['zVariables']
alldata = {name: cdf_read.varget(name) for name in zvars}
times = cdf_times_to_np_date(alldata['Timestamp'])
index_temps = np.argwhere(times == np.datetime64(date))
unique_times = compute_unique_times_list(times)
X_CF = compute_gauss_coefs_vector(swarm_file, '2017-05', 13, 300, 'B_CF', 'sigma_CF', cm_prior)
X_OB = compute_gauss_coefs_vector(swarm_file, '2017-05', 13, 300, 'B_OB', 'sigma_OB', cm_prior)
