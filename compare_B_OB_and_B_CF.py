from solving_inverse_problem import X_CF, X_OB
import matplotlib.pyplot as plt
import numpy as np
from plot_MF_inverse_problem import compute_grid_positions, compute_B_meshgrid, create_graph, compute_B


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
    B_CF = compute_B(positions, N, L_max, X_CF)
    B_OB = compute_B(positions, N, L_max, X_OB)
    B_difference = B_OB - B_CF

    #  Lines plotting it
    fig = plt.figure()
    ph, th = np.meshgrid(ph, th)
    levels_B = lambda B: np.linspace(np.min(B_difference), np.max(B_difference), N)
    BR_meshgrid = compute_B_meshgrid(B_difference, 0, N)
    Bth_meshgrid = compute_B_meshgrid(B_difference, 1, N)
    Bphi_meshgrid = compute_B_meshgrid(B_difference, 2, N)

    create_graph(fig, 221, th, ph, BR_meshgrid, "Difference OB-CF radial component", levels_B)
    create_graph(fig, 222, th, ph, Bth_meshgrid, "Difference OB-CF polar component", levels_B)
    create_graph(fig, 223, th, ph, Bphi_meshgrid, "Difference OB_CF azimuthal component", levels_B)
    plt.show()
