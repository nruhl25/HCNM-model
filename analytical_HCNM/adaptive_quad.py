# Author: Nathaniel Ruhl
# This script performs adaptive integration of the gamma vs x curve with Simpson's rule using an algorithm inspired by pseudo code on p. 642 of "Numerical Methods for Engineers, 7th ed" by Steven C. Chapra and Raymond P. Canale

import numpy as np
import matplotlib.pyplot as plt

# import local libraries
from AnalyzeCrossing import AnalyzeCrossing

ISS = AnalyzeCrossing(cb="Earth", H=420)

# Recursive function that calculates the area under gamma from a to b, with midpoint c.
# gamma curve is defined for E_kev and t
def qstep(a, b, tol, E_kev, t):
    h1 = b - a
    h2 = h1/2
    c = (a+b)/2
    d = (a+c)/2
    e = (c+b)/2
    # evaluate gamma at the desired points
    ga = ISS.gamma_vs_x(E_kev, a, t)
    gb = ISS.gamma_vs_x(E_kev, b, t)
    gc = ISS.gamma_vs_x(E_kev, c, t)
    gd = ISS.gamma_vs_x(E_kev, d, t)
    ge = ISS.gamma_vs_x(E_kev, e, t)

    # Evaluate Integrals
    I1 = (h1/6)*(ga + 4*gc + gb)
    I2 = (h2/6)*(ga + 4*gd + 2*gc + 4*ge + gb)

    # Euler-Mclaurin error
    epsilon = (I2-I1)/15

    if abs(epsilon)<=tol:
        I = I2 + epsilon   # better estimate of the integral
        return I, a, c, b
    else:
        # Cut b,the upper bound of the integral, in half
        b = c
        return qstep(a, b, tol, E_kev, t)

# This is the main function that does the adaptive qudrature and returns the total optical depth
def adaptive_simpson(E_kev, t, tol):

    # Lists for step size and distance along the LOS
    dx_list = []
    x_midpoints = []
    tau_list = []   # keep track of tau along the LOS, will sum at the end

    # First integral is over the entire domain
    a_moving = 0.0  # will update this in the while loop
    b_fixed = ISS.d_tot(t)/2

    tau = 0
    # Each iteration of the loop goes through one round of adaptive quadrature
    while a_moving < b_fixed:
        tau_i, x_lower, x_mid, x_upper = qstep(a_moving, b_fixed, tol, E_kev, t)
        x_midpoints.append(x_mid)
        dx_list.append(x_upper - x_mid)
        tau_list.append(tau_i)
        a_moving = x_upper   # Upper bound of last integral is lower bound of next integral
    tau = 2*sum(tau_list)
    return tau, dx_list, x_midpoints

def main():
    plt.rc('text', usetex=True)
    E_kev = 4.0   # keV
    t = 40.0  # sec
    tol = 1e-12
    tol_range = [1e-8, 1e-10, 1e-12, 1e-14]
    # perform adaptive quadrature with multiple tolerance goals
    for tol_i in tol_range:
        tau, dx_list, x_midpoints = adaptive_simpson(E_kev, t, tol_i)

        plt.plot(x_midpoints, dx_list, label=fr"tol={tol_i}")
        plt.legend()
    plt.xlabel('$x$ (km) along the LOS')
    plt.ylabel("$dx$ step size (km)")
    plt.title("Adaptive Simpson's rule when integration a telescopic LOS")
    plt.show()
    return 0

if __name__ == '__main__':
    main()


