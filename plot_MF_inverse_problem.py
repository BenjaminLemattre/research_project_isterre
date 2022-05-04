import numpy as np
import matplotlib.pyplot as plt
from functions import compute_grid_positions, compute_B, compute_B_meshgrid, create_graph
from solving_inverse_problem import X_CF


if __name__ == '__main__':
    '''
    make sure the date on "solving_inverse_problem" and on "plot_MF_inverse_problem" correspond
    '''
    date = '2017-05'

    # Plotting parameters
    radius = 6861
    eps = 10e-7
    N = 20
    th = np.linspace(-90+eps, 90+eps, N)
    ph = np.linspace(-180, 180, N)
    positions = compute_grid_positions(radius, th, ph, N)

    # Compute B
    L_max = 13
    B = compute_B(positions, N, L_max, X_CF)

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