from functools import lru_cache
import numpy as np
from scipy.special import lpmv
import math
import cdflib
import h5py
import matplotlib.pyplot as plt

r_earth = 6371.2
r_iono = 60000


def compute_direct_obs_operator(positions, max_degree, internal_source=True, eps = 10e-7):
    """
    Computes the direct observation operator H of a measure of given max degree for a list of given positions.

    :param positions: List of the positions stored under the form [r1, th1, ph1, r2, th2, ph2 ...]
    :type positions: list or 1D numpy.ndarray
    :param max_degree: Maximal degree of the observed measure
    :type max_degree: int
    :param internal_source: whether the Gauss coefficient represent an internal or an external source (a/r or r/a dependency)
    :type internal_source: bool
    :return: the direct observation operator H
    :rtype: 2D numpy.ndarray (nb_observations x nb_coefs)
    """
    # Define radius
    r_earth = 6371.2
    r_iono = 60000
    # Number of observation parameters
    Nobs = len(positions)
    # Number of Gauss coefficients
    Ncoefs = max_degree * (max_degree + 2)

    H = np.zeros([Nobs, Ncoefs])

    for i in range(0, Nobs, 3):
        r = positions[i]
        th = positions[i + 1]
        ph = positions[i + 2]
        H[i:i + 3, :] = compute_lines_direct_obs_operator(r, th, ph, internal_source, eps, max_degree)
    return H


@lru_cache(maxsize=2048)
def compute_lines_direct_obs_operator(r, th, ph, internal_source, eps, max_degree):
    """
    Convenience function, as most of the the observations positions are similar through time (GVOs are unchanged for each satellite era),
    a cache is used to avoid computing the same coefficients for each time step.
    """
    global a_on_b, sign_cos_th
    Ncoefs = max_degree * (max_degree + 2)
    r_th_ph_coefs = np.zeros(shape=(3, Ncoefs))
    # Compute prefactors
    if internal_source:
        radius_to_earth = r_earth / r
    else:
        radius_to_earth = r / r_earth
        a_on_b = r_earth / r_iono
    sin_th = math.sin(th)
    if abs(sin_th) == 0:  # theta is either O or pi
        sin_th_0 = True
        sign_cos_th = np.sign(math.cos(th))
    else:
        sin_th_0 = False

    for k in range(0, Ncoefs):
        n, m, coef = get_degree_order(k)
        plm = schmidt_semi_normalized_plm(n, m, th)
        dplm = derivativeThetaAssocLegendrePoly(n, m, th)
        cos_mph = math.cos(m * ph)
        sin_mph = math.sin(m * ph)
        if internal_source:
            if coef == "g":
                r_th_ph_coefs[0, k] = (n + 1) * radius_to_earth ** (n + 2) * plm * cos_mph
                r_th_ph_coefs[1, k] = - radius_to_earth ** (n + 2) * dplm * cos_mph
                if sin_th_0:
                    #  assert sin_mph == 0 or plm == 0, 'One of these values should be 0, {}, {}'.format(sin_mph, plm)
                    r_th_ph_coefs[2, k] = radius_to_earth ** (n + 2) * m * sign_cos_th * dplm * sin_mph
                else:
                    r_th_ph_coefs[2, k] = radius_to_earth ** (n + 2) * m * plm * sin_mph / sin_th
            else:
                r_th_ph_coefs[0, k] = (n + 1) * radius_to_earth ** (n + 2) * plm * sin_mph
                r_th_ph_coefs[1, k] = -radius_to_earth ** (n + 2) * dplm * sin_mph
                if sin_th_0:
                    #  assert cos_mph == 0 or plm == 0, 'One of these values should be 0, {}, {}'.format(cos_mph, plm)
                    r_th_ph_coefs[2, k] = -radius_to_earth ** (n + 2) * m * sign_cos_th * dplm * cos_mph
                else:
                    r_th_ph_coefs[2, k] = -radius_to_earth ** (n + 2) * m * plm * cos_mph / sin_th
        else:  # if the source is external to the position, then the dependence in the radius is as (r / a)^l-1 instead of (a / r)^l+2 shown in Sabaka 2010
            if coef == "g":
                r_th_ph_coefs[0, k] = (n + 1) * a_on_b ** (2 * n + 1) * radius_to_earth ** (n - 1) * plm * cos_mph
                r_th_ph_coefs[1, k] = (n + 1) / n * a_on_b ** (2 * n + 1) * radius_to_earth ** (n - 1) * dplm * cos_mph
                assert not sin_th_0
                r_th_ph_coefs[2, k] = - (n + 1) / n * a_on_b**(2*n+1) * radius_to_earth ** (n - 1) * m * plm * sin_mph / sin_th
            else:
                r_th_ph_coefs[0, k] = (n + 1) * a_on_b**(2*n+1) * radius_to_earth ** (n - 1) * plm * sin_mph
                r_th_ph_coefs[1, k] = (n + 1) / n * a_on_b**(2*n+1) * radius_to_earth ** (n - 1) * dplm * sin_mph
                assert not sin_th_0
                r_th_ph_coefs[2, k] = (n + 1) / n * a_on_b**(2*n+1) * radius_to_earth ** (n - 1) * m * plm * cos_mph / sin_th
    return r_th_ph_coefs


def get_degree_order(k):
    """
    Gets the degree, order and coefficient name ("g" or "h") of a coefficient.
    :param k: index of the coefficient
    :type k: int
    :return: degree, order and coef name of the coefficient
    :rtype: int, int, str
    :rq: !! still need to compute the coef q and s !!
    """

    assert (int(k) == k and k >= 0)

    # For m = 0, k = l**2 -1.
    floating_sqrt_result = math.sqrt(k + 1)
    degree = int(floating_sqrt_result)
    if degree == floating_sqrt_result:
        return degree, 0, "g"

    # We need now to find m verifying:
    #    for g : k = l**2 + 2m - 2 => 2m = k - l**2 + 2
    # OR for h : k = l**2 + 2m - 1 => 2m = k - l**2 + 1

    twice_order = k - degree * degree + 2
    if twice_order % 2 == 0:
        coef = "g"
    else:
        coef = "h"
        # If twice_order is not divisible by 2, it means that the 2m is in fact (twice_order - 1)
        twice_order = twice_order - 1

    order = twice_order // 2
    return degree, order, coef


def norm_coef(n, m):
    if m == 0:
        return 1
    else:
        return np.sqrt(2 * ((math.factorial(n - m)) / (math.factorial(n + m))))


def full_norm(n, m):
    return np.sqrt(2 * n + 1) * np.sqrt(np.math.factorial(n - m) / np.math.factorial(n + m))


def assocLegendrePoly(n, m, theta):
    Plm = lpmv(m, n, np.cos(theta)) * ((-1) ** m)
    return Plm


def derivativeThetaAssocLegendrePoly(n, m, theta):
    Plm = (n + 1) * (-np.cos(theta) / np.sin(theta)) * lpmv(m, n, np.cos(theta)) + \
          (n - m + 1) * (1 / np.sin(theta)) * lpmv(m, n + 1, np.cos(theta))
    return Plm * ((-1) ** m) * norm_coef(n, m)


def secondDerivativeThetaAssocLegendrePoly(n, m, theta):
    Plm = (n + 1) * ((n + 1) * lpmv(m, n, np.cos(theta)) * np.cos(theta) -
                     (n - m + 1) * lpmv(m, n + 1, np.cos(theta))) * np.cos(theta) + \
          (n + 1) * lpmv(m, n, np.cos(theta)) * (np.sin(theta) ** 2) + \
          (n + 1) * lpmv(m, n, np.cos(theta)) * (np.cos(theta) ** 2) + \
          ((n + 2) * lpmv(m, n + 1, np.cos(theta)) * np.cos(theta) -
           (n - m + 2) * lpmv(m, n + 2, np.cos(theta))) * \
          (-n + m - 1) - (n - m + 1) * lpmv(m, n + 1, np.cos(theta)) * np.cos(theta)
    return Plm * ((-1) ** m) / (np.sin(theta) ** 2)



def compute_grid_positions(r, theta, phi, N):
    """
    params:
    r: radius in km
    theta: list of colatitude in deg
    phi: list of longitudes in deg
    N: number of colatitudes/longitudes in the grid
    return:
    list of lengh 3*N^2 positions [r, th1, ph1, r, th1, ph2 ...]
    """
    pos = []
    for i in range(N):
        for j in range(N):
            pos.append(r)
            pos.append(np.deg2rad(theta[i]))
            pos.append(np.deg2rad(phi[j]))
    return pos


def schmidt_semi_normalized_plm(n, m, th):
    return assocLegendrePoly(n, m, th) * norm_coef(n, m)


def compute_list_coords_obs_at_t(d, data, times_list, N_VO):
    """
    params:
    d: date chosen
    data: dict of data from file
    times_list: list of times from the dict data
    N_VO: number of observatories
    return: 
    list of coords of observatories of lengh 3*N_VO [r1, th1, ph1, r2 ...]
    """
    index_t = np.argwhere(times_list == np.datetime64(d))[0][0]
    slice_at_time = slice(N_VO*index_t, N_VO*(index_t + 1))
    radius_list_at_t = data['Radius'][slice_at_time] / 1000
    th_list_at_t = data['Latitude'][slice_at_time]
    ph_list_at_t = data['Longitude'][slice_at_time]
    list_coords_at_t = []
    for k in range(len(radius_list_at_t)):
        list_coords_at_t.append(radius_list_at_t[k])
        list_coords_at_t.append(np.deg2rad(90-th_list_at_t[k]))
        list_coords_at_t.append(np.deg2rad(ph_list_at_t[k]))
    return radius_list_at_t, th_list_at_t, ph_list_at_t, list_coords_at_t


def compute_cov_matrix_obs_at_t(list_sigma):
    """
    Takes the list of standard deviation errors at locations and form a
    covariance matrix N_VO x N_VO
    """
    cm_obs = np.zeros((len(list_sigma), len(list_sigma)))
    for k in range(0, len(list_sigma), 3):
        cm_obs[k, k] = list_sigma[k]**2
        cm_obs[k+1, k+1] = list_sigma[k+1]**2
        cm_obs[k+2, k+2] = list_sigma[k+2]**2
    return cm_obs


def compute_obs_vector_at_t(data, d, N_VO, key, times_list):
    """
    params:
    data: dict of data from file
    d: date chosen 
    N_VO: number of observatories
    key: key chosen from data
    times_list: times list from data 
    returns: 
    observations vector with shape(N_VO,3) [[BR1, Bth1, Bph1],[BR2, Bth2, Bph2]...]
    """
    index_t = np.argwhere(times_list == np.datetime64(d))[0][0]
    obs = data[key][N_VO*index_t:N_VO*(index_t + 1), :]
    obs_vec = []
    for k in range(len(obs)):
        obs_vec.append(obs[k, 0])
        obs_vec.append(obs[k, 1])
        obs_vec.append(obs[k, 2])
    return obs_vec


def delete_nan_values(observations, positions, list_sigma):
    """
    params: 
    observations: list of MF observations of lengh 3*N_VO
    positions: list of observatories positions of lengh 3*N_VO
    list_sigma: list of sigma of lengh 3*N_VO
    returns:
    observations, positions, list_sigma with nan values deleted
    """
    for k in range(len(observations)-1, -1, -3):
        L = [np.isnan(observations[k]), np.isnan(observations[k-1]), np.isnan(observations[k-2])]
        if np.any(np.array(L)):
            observations = np.delete(observations, [k, k-1, k-2], axis=None)
            positions = np.delete(positions, [k, k-1, k-2], axis=None)
            list_sigma = np.delete(list_sigma, [k, k-1, k-2], axis=None)
    return observations, positions, list_sigma


def compute_gauss_coefs_vector(file, date_chosen, max_degree, N_VO, measure_MF, measure_sigma, cm_prior):
    """
    params: 
    file: file chosen (swarm for instance, with os.path.join and cdf_dir)
    date_chosen: date chosen
    max_degree: max degree chosen
    N_VO: number of observatories
    measure_MF: list of MF observations with lengh 3*N_VO
    measure_sigma: list of sigma with lengh 3*N_VO
    cm_prior: prior covariance matrix with shape (3*N_VO^2, 3*N_VO^2)
    returns:
    gauss coefficient vector computed with the parameters chosen
    """
    cdf_read = cdflib.CDF(file)
    info = cdf_read.cdf_info()
    zvars = info['zVariables']
    alldata = {name: cdf_read.varget(name) for name in zvars}
    times = cdf_times_to_np_date(alldata['Timestamp'])
    unique_times = compute_unique_times_list(times)
    list_sigma_obs = compute_obs_vector_at_t(alldata, date_chosen, N_VO, measure_sigma, unique_times)
    obs_vec_at_t = compute_obs_vector_at_t(alldata, date_chosen, N_VO, measure_MF, unique_times)
    r, th, ph, list_coords_obs = compute_list_coords_obs_at_t(date_chosen, alldata, unique_times, N_VO)
    obs_vec_at_t, list_coords_obs, list_sigma_obs = delete_nan_values(obs_vec_at_t, list_coords_obs,
                                                                            list_sigma_obs)
    cm_obs = compute_cov_matrix_obs_at_t(list_sigma_obs)
    H = compute_direct_obs_operator(list_coords_obs, max_degree)
    mean_prior = np.zeros(max_degree * (max_degree + 2))
    sigma_obs_inv = np.linalg.inv(cm_obs)
    sigma_prior_inv = np.linalg.inv(cm_prior)
    Kalman_gain = np.linalg.inv(H.T @ sigma_obs_inv @ H + sigma_prior_inv)
    return Kalman_gain @ (sigma_prior_inv @ mean_prior + H.T @ sigma_obs_inv @ obs_vec_at_t)


def cdf_times_to_np_date(times):
    """
    Transform the times in cdf epoch into numpy dates, rounding month to nearest.
    """
    dates = []
    for t in times:
        if t != t:
            # if time is a nan, date should be a nan as well
            dates.append(np.float('nan'))
            continue
        year, month, day = cdflib.cdfepoch.breakdown_epoch(t)[:3]
        if day > 15:
            if month == 12:
                year += 1
                month = 1
            else:
                month += 1
        dates.append(np.datetime64('{}-{:02}'.format(year, month)))
    return dates


def read_observatories_cdf_data(cdf_filename, measure_type, obs, remove_bias_crust=True, huge_sigma=99999):
    """
    Reads the data of observatories file in cdf format.

    :param huge_sigma: Replace Nans in sigma by this big value
    :type cdf_filename: str
    :param cdf_filename: name of the .cdf file
    :param measure_type: measure type MF = Main field, SV = Secular variation
    :param obs: observation GO = Ground Observatories, CH = Champ, SW = Swarm, OR = Oersted, CR = Cryosat , C0= Composite
    """
    global times, Bs, sigma2
    possible_obs = np.array(['GO', 'CH', 'SW', 'OR', 'CR', 'CO'])

    if not np.any(obs == possible_obs):  # error test
        raise ValueError('Got {} as obs value, possible values are {}'.format(obs, possible_obs))
    if not np.any(measure_type == np.array(['MF', 'SV'])):
        raise ValueError('Got {} as measure_type, possible values are {}'.format(measure_type, ('MF', 'SV')))

    cdf_file = cdflib.CDF(cdf_filename)

    thetas = 90 - cdf_file.varget("Latitude")  # convert in colatitude

    if np.any(thetas > 180) or np.any(thetas < 0):  # error test
        raise ValueError('angle th {} represent the colatitude, did you use latitudes instead?'.format(thetas))

    thetas = np.deg2rad(thetas)  # convert in rad
    phis = np.deg2rad(cdf_file.varget("Longitude"))  # longitudes in rad
    radii = cdf_file.varget('Radius') / 1000  # to convert in kms
    assert np.amax(radii) < 10000 and np.amin(radii) > 6000, 'Are you sure radii are stored in meters?'  # error test

    if measure_type == 'MF':
        times = cdf_file.varget('Timestamp')
        Bs, sigma2 = cdf_file.varget('B_CF'), cdf_file.varget('sigma_CF')**2
        if remove_bias_crust and obs == 'GO':
            Bs -= cdf_file.varget('bias_crust')
    elif measure_type == 'SV':
        times = cdf_file.varget('Timestamp_SV')
        Bs, sigma2 = cdf_file.varget('B_SV'), cdf_file.varget('sigma_SV')**2

    sigma2[np.where(np.isnan(sigma2))] = huge_sigma # replace Nans

    err = 'At least some invalid values of sigma correspond to valid values in B, file, measure and obs {}, {}, {}'
    err = err.format(cdf_filename.split('/')[-1], measure_type, obs)

    if obs == 'GO':
        locs = [''.join([l[0] for l in loc]) for loc in cdf_file.varget('Obs')]
        N_VO = len(list(dict.fromkeys(zip(thetas, phis, locs))))  # small trick that counts the number of unique pairs
        # of theta/phi/locs
    else:
        N_VO = len(dict.fromkeys(zip(thetas, phis)).keys())  # small trick that counts the number of unique pairs of
        # theta/phi
    dates = cdf_times_to_np_date(times)
    return dates, Bs


def slice_for_obs_ids(cdf_filename, radius=False):
    """
    Params:
    cdf_filename: .cdf file
    Returns:
    alldata: dict that contains all data from cdf file 'file', the keys are
              alldata[zvar_name] = zvar_data
    unique_obs_ids: list of the ordered locations of the observatories with their latitude and longitude
    mask_obs_id: dict. Given a location (th, ph) of an observatory, returns the array indices
                 such that data[mask_obs_id[obs_id]] = data from observatories at different times,
                 where data is e.g. B_CF.
    """
    global unique_observatories_ids, full_obs_ids
    cdf_read = cdflib.CDF(cdf_filename)
    info = cdf_read.cdf_info()

    zvars = info['zVariables']

    #  V0 = Virtual observatories GVO = Geomagnetic virtual observatories
    obs_type = 'GO' if 'GObs' in cdf_filename else 'VO'
    err = 'Allowed values for obs_type are GO or VO, got {} instead'.format(obs_type)
    assert (obs_type == 'GO' or obs_type == 'VO'), err
    #  -------------------------------------

    alldata = {name: cdf_read.varget(name) for name in zvars}  # alldata is a dictionary : alldata['zVariables'] to get those

    thetas = 90 - alldata['Latitude']  # list of thetas
    phis = alldata['Longitude']  # list of phis

    if obs_type == 'VO':
        if not radius:
            full_obs_ids = np.array([(th, ph) for th, ph in zip(thetas, phis)])
            unique_observatories_ids = list(dict.fromkeys(zip(thetas, phis)))
            print('unique_observatories_ids')
            print(unique_observatories_ids)
            print('\n')
        else:  # This part is never taken in account because radius = False
            radiuss = alldata['Radius']
            full_obs_ids = np.array([(r, th, ph) for r, th, ph in zip(radiuss, thetas, phis)])
            unique_observatories_ids = list(dict.fromkeys(zip(radiuss, thetas, phis)))

    elif obs_type == 'GO':
        # shape locs from numpy array to list of strings containing locations
        locs = [''.join([l[0] for l in loca]) for loca in alldata['Obs']]
        if not radius:
            full_obs_ids = np.array([(th, ph, loca) for th, ph, loca in zip(thetas, phis, locs)])
            unique_observatories_ids = list(dict.fromkeys(zip(thetas, phis, locs)))
        else:  # This part is never taken in account because radius = False
            radiuss = alldata['Radius']
            full_obs_ids = np.array([(r, th, ph, loca) for r, th, ph, loca in zip(radiuss, thetas, phis, locs)])
            unique_observatories_ids = list(dict.fromkeys(zip(radiuss, thetas, phis, locs)))

    N_VO = len(unique_observatories_ids)  # Number of observatories
    print('Number of observatories')
    print(N_VO)
    print('\n')
    print('Number of points')
    print(len(full_obs_ids))
    print('\n')

    # test is too strong since order of obs_ids does not matter
    assert(np.all(unique_observatories_ids == full_obs_ids[:N_VO]))
    # not true for GO since it is composed of unequal series
    assert(len(full_obs_ids) % N_VO == 0)

    # dictionary whose keys are lats and lons.
    # values of the dict returns all index of data corresponding to the location given by the key
    # e.g. if (th_0, ph_0) correspond to the first VO, mask_obs_id[(th_0, ph_0)] = array(0, 300, 600, ...)
    # (th_, ph_1) to the second VO, mask_obs_id[(th_0, ph_0)] = array(1, 301, 601, ...)
    # since there is 300 VOs.
    # therefore, B_CF[masks_obs_id[(0.0, 5.0)] ] returns the MF at latitude = 0 and longitude = 5
    # at all available times
    # For GO obs_id also contains the name of the observatory
    mask_obs_loc = {'MF': {}, 'SV': {}}
    for obs_id in unique_observatories_ids:
        mask_obs_loc['MF'][obs_id] = (full_obs_ids == obs_id).all(axis=1).nonzero()[0]

    # for new brut series, the timestamp_SV has not the same length as the timestamp.
    # This calls for another mask in that case
    if alldata['Timestamp'].shape != alldata['Timestamp_SV'].shape:
        # loop will stop at the end of timestamp_SV
        full_obs_ids_SV = np.array([(th, ph) for th, ph, t_SV in zip(thetas, phis, alldata['Timestamp_SV'])])
        for obs_id in unique_observatories_ids:
            mask_obs_loc['SV'][obs_id] = (full_obs_ids_SV == obs_id).all(axis=1).nonzero()[0]
    else:
        mask_obs_loc['SV'] = mask_obs_loc['MF']
    return alldata, unique_observatories_ids, mask_obs_loc


def compute_unique_times_list(times_list):
    """
    params: times list from data 
    returns: unique time list of data observations
    """
    list_t = []
    for t in times_list:
        if t not in list_t:
            list_t.append(t)
    return list_t


def compute_MF_at_loc(r, th, ph, gauss_coefs, index_temps: int, max_degree = 20):
    H = compute_direct_obs_operator([r, th, ph], max_degree)
    B = np.dot(H, gauss_coefs)
    return B[0], B[1], B[2]


def compute_grid_positions(r, theta, phi, N):
    """
    r: radius in km
    theta: list of colatitude in deg
    phi: list of longitudes in deg
    """
    pos = []
    for i in range(N):
        for j in range(N):
            pos.append(r)
            pos.append(np.deg2rad(theta[i]))
            pos.append(np.deg2rad(phi[j]))
    return pos


def create_graph(fig, zone, world, list_th, list_ph, B_choice, title):
    """
    params: 
    fig: figure
    zone: zone for subplot
    world: map of the world with geopandas
    list_th, list_ph: array of colatitudes/ longitudes previously meshed with lengh (N^2,N^2)
    """
    levels_B = lambda B_choice: np.linspace(np.min(B_choice), np.max(B_choice), 300)
    ax1 = fig.add_subplot(zone)
    world.plot(ax=ax1)
    map_cf = ax1.contourf(list_ph, list_th, B_choice, levels=levels_B(B_choice), alpha=0.3, cmap='bwr')
    plt.colorbar(map_cf, ax=ax1)
    ax1.set_title(title)   


def compute_B_meshgrid(Magnetic_Field_at_time, choice, Nloc):
    """
    Magnetic_Field_at_time: list of Magnetic field at a given time index
    choix: 0 for BR, 1 for Bth, 2 for Bphi
    Nloc: length of latitudes list and longitudes list
    """
    B_meshgrid = np.zeros((Nloc, Nloc))
    B_choice = []
    for i in range(choice, Magnetic_Field_at_time.shape[0], 3):
        B_choice.append(Magnetic_Field_at_time[i])
    index_loc = 0
    for i in range(Nloc):
        for j in range(Nloc):
            B_meshgrid[i][j] = B_choice[index_loc]
            index_loc += 1
    return B_meshgrid


def compute_B(pos, N_loc, max_degree, X):
    """
    params: 
    pos: list of positions with lengh 3*N_loc^2 (radius, colatitude, longitudes)
    N_loc: number of colatitude/longitude for the grid 
    max_degree: max_degree chosen
    X: gauss coefficient vector computed with previous parameters
    returns:
    estimated list of MF at lengh 3*N_loc^2
    """
    MF = np.zeros(3 * N_loc * N_loc)
    for k in range(0, len(pos), 3):
        H = compute_direct_obs_operator([pos[k], np.deg2rad(90) - pos[k + 1], pos[k + 2]], max_degree)
        BR, Bth, Bphi = np.dot(H, X)
        MF[k], MF[k + 1], MF[k + 2] = BR, Bth, Bphi
    return MF


def extract_MF_list(Magnetic_field, choice):
    """
    params: 
    Magnetic_field: list of magnetic field from data
    choice: choice of component to extract
    returns: 
    list of magnetic field according to component chosen with lengh N_VO
    """
    B_choice = []
    for k in range(len(Magnetic_field)):
        B_choice.append(Magnetic_field[k][choice])
    return B_choice


def delete_nan_values_t(times_list):
    """
    params: 
    times_list: list of times from data
    returns: 
    times_list with nan values deleted
    """
    for k in range(len(times_list)-1, -1):
        if math.isnan(times_list[k]):
            print(np.isnan(times_list[k]))
            times_list = np.delete(times_list, k, axis=None)
    return times_list


def compute_gauss_coefs_matrix(file, max_degree, N_VO, measure_MF, measure_sigma, cm_prior, unique_time_list):
    """
    params: same of compute_gauss_coefs_vector apart from date
    returns: gauss coefficients matrix with shape(max_degree*(max_degree+2), len(unique_time_list))
    """
    gauss_matrix = np.zeros((max_degree*(max_degree+2),len(unique_time_list)))
    j = 0
    for date in unique_time_list:
        X_column = compute_gauss_coefs_vector(file, date, max_degree, N_VO, measure_MF, measure_sigma, cm_prior)
        for k in range(len(X_column)):
            gauss_matrix[k, j] = X_column[k]
        j += 1
    return gauss_matrix


def create_graph_compare_gauss_coefs(i, zone, unique_times, X_matrix_OB, X_matrix_CF):
    """
    params: 
    i: choice of the index of gauss coefficient to plot
    zone: zone to subplot
    unique_times: unique times list from our file
    X_matrix_OB: gauss coefficients matrix computed with data 'B_OB'
    X_matrix_CF: gauss coefficients matrix computed with data 'B_CF'
    returns:
    graph comparing the gauss coefficients calculated with 'B_OB' and 'B_CF' chosen
    """
    g_ob_chosen = np.zeros(len(unique_times))
    g_cf_chosen = np.zeros(len(unique_times))
    for k in range(len(unique_times)): 
        g_ob_chosen[k] = X_matrix_OB[i, k]
        g_cf_chosen[k] = X_matrix_CF[i, k]
    ax1 = plt.subplot(zone)
    ax2 = plt.subplot(zone)
    ax1.plot(unique_times, g_ob_chosen,'b', label='g_ob')
    ax2.plot(unique_times, g_cf_chosen,'r', label='g_cf')
    ax1.legend(fontsize='xx-large')
    ax2.legend(fontsize='xx-large')
    plt.grid()
    plt.title(label='Comparison between OB and CF data in gauss coefficient %i'%i, fontsize='xx-large')

def create_hist(fig, i, zone, alldata, L_max, unique_times_list, date, file, cm_prior, title):
    N_VO = 300
    index_temps = np.argwhere((cdf_times_to_np_date(alldata['Timestamp'])) == np.datetime64(date))[0][0]
    radius = alldata['Radius'][index_temps:index_temps+N_VO]*10**(-3)
    latitudes = alldata['Latitude'][index_temps:index_temps+N_VO]
    longitudes = alldata['Longitude'][index_temps:index_temps+N_VO]
    r, th, ph, list_coords_obs = compute_list_coords_obs_at_t(date, alldata, unique_times_list, N_VO)
    B_CF_real_choice = extract_MF_list(alldata['B_CF'][index_temps:index_temps+N_VO], i)
    list_sigma_choice = extract_MF_list(alldata['sigma_CF'][index_temps:index_temps+N_VO], i)
    
    X_CF = compute_gauss_coefs_vector(file, date, L_max, N_VO, 'B_CF', 'sigma_CF', cm_prior)
    H = compute_direct_obs_operator(list_coords_obs, L_max)
    B_CF_estimated = H @ X_CF

    B_CF_estimated_choice = np.zeros(N_VO)
    for k in range(N_VO): 
        B_CF_estimated_choice[k] = B_CF_estimated[i+3*k]
    B_CF_diff_choice = np.zeros(N_VO)
    for k in range(N_VO):
        B_CF_diff_choice[k] = (B_CF_real_choice[k] - B_CF_estimated_choice[k]) / list_sigma_choice[k]
    B_CF_real_choice, B_CF_estimated_choice, list_sigma_choice = delete_nan_values(B_CF_real_choice, B_CF_estimated_choice, list_sigma_choice)
    m = sum(B_CF_real_choice - B_CF_estimated_choice)/len(B_CF_real_choice - B_CF_estimated_choice)
    
    list_index = []
    for k in range(len(B_CF_real_choice)):
        if abs(B_CF_real_choice[k] - B_CF_estimated_choice[k]) >3:
            list_index.append(k)
    print('> number of value too far in component %i : '%i, len(list_index))
    i = 0
    for k in range(len(list_index)):
        if abs(latitudes[list_index[k]]) > 50:
            i+=1
    print('> number of those values with |lats| > 60 : ', i)
            

    ## Lines plotting the geomagnetic field 
    ax1=fig.add_subplot(zone)
    ax1.set_title(title)
    ax1.hist(B_CF_diff_choice, bins=50, label='mean value = %1.4f'%m)
    ax1.legend()