import matplotlib.pyplot as plt
import numpy as np
from functions import compute_grid_positions, compute_MF_at_loc, compute_B_meshgrid, create_graph, times_chaos, gauss_coeffs


if __name__ == '__main__':
    radius = 6861
    eps = 10e-7
    N = 20
    L_max = 20
    latitudes = np.linspace(0, 180, N)
    longitudes = np.linspace(-180, 180, N)
    th = 90 - latitudes + eps
    ph = longitudes
    date = 2018.
    positions = compute_grid_positions(radius, th, ph, N)
    index_temps = np.argwhere(times_chaos == date)
    B = np.zeros(3*N*N)
    for k in range(0, len(positions), 3):
        B[k], B[k+1], B[k+2] = compute_MF_at_loc(radius, np.deg2rad(90)-positions[k+1], positions[k+2], gauss_coeffs, index_temps)
    #  Lines plotting it
    fig = plt.figure()
    ph, th = np.meshgrid(ph, th)
    levels_B = lambda B: np.linspace(np.min(B), np.max(B), N)
    BR_meshgrid = compute_B_meshgrid(B, 0, N)
    Bth_meshgrid = compute_B_meshgrid(B, 1, N)
    Bphi_meshgrid = compute_B_meshgrid(B, 2, N)
    create_graph(fig, 221, th, ph, BR_meshgrid, "Geomagnetic field's radial component", levels_B)
    create_graph(fig, 222, th, ph, Bth_meshgrid, "Geomagnetic field's polar component", levels_B)
    create_graph(fig, 223, th, ph, Bphi_meshgrid, "Geomagnetic field's azimuthal component", levels_B)
    plt.show()
