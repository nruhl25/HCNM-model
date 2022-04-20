# Author: Nathaniel Ruhl
# This script analyzes the convergence of the adaptive quadrature method when integrating the LOS

# Import local modules
from AnalyzeCrossing import AnalyzeCrossing

# import standard libraries
import numpy as np
import time
import matplotlib.pyplot as plt
plt.rc("text", usetex=True)

# Define constants
E_kev = 4.0
ES = AnalyzeCrossing(cb="Earth", H=420)   # used for all analysis
time_array = np.arange(0, ES.time_final+1, 1)

# Calculates a transmittance curve with adaptive quadrature for a given tolerance, and returns the transmittance array, run-time

def calculate_transmit_quadapt(tol):

    transmit_adapt = np.zeros_like(time_array)

    # Time the full horizon crossing (over 300 evaluations of the integral)
    start_time = time.time()

    for i, t in enumerate(time_array):
        tau_adaptive, dx_list, x_midpoints = ES.tau_adaptive_simpson(t, tol=tol)
        transmit_adapt[i] = np.exp(-tau_adaptive)
    run_time = time.time() - start_time
    return transmit_adapt, run_time

# Calculates a transmittance curve with N steps via gaussian quadrature, and returns the transmittance array and run-time
def calculate_transmit_gauss(N):

    transmit_gauss = np.zeros_like(time_array)

    # Time the full horizon crossing (over 300 evaluations of the integral)
    start_time = time.time()

    for i, t in enumerate(time_array):
        transmit_gauss[i] = np.exp(-ES.tau_gauss(t, N))
    run_time = time.time() - start_time

    return transmit_gauss, run_time


def main():
    # define a dictionary containing N (int) as the keys and the transmittance array as the values
    transmit_dict = {}

    tol_list = [1e-5, 1e-6, 1e-7, 1e-8, 1e-9]
    # list of run-times corresponding to transmit_dict
    run_time_list = []
    error_list = []   # difference from the "truth value"

    # Take this as the "truth value"
    transmit_best, run_time_best = calculate_transmit_gauss(N=100)

    comp_range = np.where((transmit_best > 0.01) & (transmit_best < 0.99))[0]

    for tol in tol_list:
        transmit_array, run_time = calculate_transmit_quadapt(tol)
        transmit_dict[tol] = transmit_array
        run_time_list.append(run_time)
        error_list.append(
            np.sum(np.abs(transmit_array[comp_range] - transmit_best[comp_range])/transmit_best[comp_range])/len(comp_range))

    # Plot the results
    for tol in tol_list:
        plt.figure(1)
        plt.plot(time_array, transmit_dict[tol], label=f"tol={tol}")

    plt.figure(1)
    plt.title("Transmittance curves calculated with Adaptive Quadrature")
    plt.xlabel(r"Time since $t_0$ (seconds)")
    plt.xlim([45, 80])
    plt.ylabel("Transmittance")
    plt.legend()

    plt.figure(2)
    plt.title("Run time of Adaptive Quadrature")
    plt.plot(tol_list, run_time_list, "-o")
    plt.xscale('log')
    plt.xlabel("Optical Depth Tolerance")
    plt.ylabel("Run time (seconds)")

    plt.figure(3)
    plt.title("Convergence of Adaptive Quadrature")
    plt.plot(tol_list, error_list, '-o')
    plt.xscale('log')
    plt.xlabel("Optical Depth Tolerance")
    plt.ylabel("Fractional Difference from Guassian Quadrature with N=100")

    plt.show()

    return 0


if __name__ == '__main__':
    main()
