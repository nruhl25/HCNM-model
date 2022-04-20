# Author: Nathaniel Ruhl
# This file analyzes the converges and runtime of gaussian quadrature when integrating the LOS

import numpy as np
import time
import matplotlib.pyplot as plt
plt.rc("text", usetex=True)

from AnalyzeCrossing import AnalyzeCrossing

# Define constants
E_kev = 4.0
ES = AnalyzeCrossing(cb="Earth", H=420)   # used for all analysis
time_array = np.arange(0, ES.time_final+1, 1)

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

    N_list = np.arange(6, 13, 1)
    # list of run-times corresponding to transmit_dict
    run_time_list = []
    error_list = []   # difference from the "truth value"
    
    # Take this as the "truth value"
    transmit100, run_time_100 = calculate_transmit_gauss(100)

    comp_range = np.where((transmit100 > 0.01) & (transmit100 < 0.99))[0]
    
    for i, N in enumerate(N_list):
        transmit_array, run_time = calculate_transmit_gauss(N)
        transmit_dict[N] = transmit_array
        run_time_list.append(run_time)
        error_list.append(np.sum(np.abs(transmit_array[comp_range] - transmit100[comp_range])/transmit100[comp_range])/len(comp_range))

    
    # Plot the results
    for i, N in enumerate(N_list):
        plt.figure(1)
        plt.plot(time_array, transmit_dict[N], label = f"N={N}")

    plt.figure(1)
    plt.title("Transmittance curves calculated with Gaussian Quadrature")
    plt.xlabel(r"Time since $t_0$")
    plt.ylabel("Transmittance")
    plt.legend()

    plt.figure(2)
    plt.title("Run time of Gaussian Quadrature")
    plt.plot(N_list, run_time_list, "-o")
    plt.xlabel("Number of Slices")
    plt.ylabel("Run time (seconds)")

    plt.figure(3)
    plt.title("Convergence of Gaussian Quadrature")
    plt.plot(N_list, error_list, "-o")
    plt.xlabel("Number of Slices")
    plt.ylabel("Fractional Difference from Gaussian Quadrature with N = 100")

    plt.show()

    return 0


if __name__ == '__main__':
    main()
