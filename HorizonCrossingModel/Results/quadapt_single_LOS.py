# This file analyzes the step sizes used by adaptive quadrature along a single LOS

import matplotlib.pyplot as plt
plt.rc('text', usetex=True)

from AnalyzeCrossing import AnalyzeCrossing

# Define constants
ES = AnalyzeCrossing(cb="Earth", H=420)
E_kev = 4.0   # keV
t = 50.0  # sec
tol = 1e-12
tol_range = [1e-7, 1e-8]

def main():

    for tol_i in tol_range:
        tau, dx_list, x_midpoints = ES.tau_adaptive_simpson(t, tol_i)
        
        plt.plot(x_midpoints, dx_list, "-o", label=fr"tol={tol_i}")
        plt.legend()
    plt.title("Step Sizes Identified in Adaptive Quadrature")
    plt.xlabel('$x$ (km) along the line of sight')
    plt.ylabel("$dx$ step size (km)")
    plt.show()
    return 0


if __name__ == '__main__':
    main()
