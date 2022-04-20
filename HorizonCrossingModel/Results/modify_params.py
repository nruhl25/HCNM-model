# Author: Nathaniel Ruhl

# This script checks to see how the changes in cross section, surface density, and scale height change a transmission curve

import numpy as np
import matplotlib.pyplot as plt

from AnalyzeCrossing import AnalyzeCrossing

# Function to calculate a transmittance array for a given SAT/orbital parameters
def calc_transmit(SAT, time_array):
    transmit_array = np.zeros_like(time_array)
    for i, t in enumerate(time_array):
        transmit_array[i] = np.exp(-SAT.tau_gauss(t, N=10))
    return transmit_array

# The functions below change parameters and make plots

def change_sigma():
    plt.figure()
    plt.title("Transmittance vs Time for Standard LEO")
    SAT = AnalyzeCrossing(cb="Earth", H=420, E_kev=4)
    time_array = np.arange(0, SAT.time_final + 1, 1, dtype=float)
    alpha_list = [0.8, 0.9, 1.0, 1.1, 1.2]    # list of scaling factors for cross section
    for alpha in alpha_list:
        SAT.sigma = alpha*SAT.sigma
        transmit_array = calc_transmit(SAT, time_array)
        plt.plot(time_array, transmit_array,
                 label=f"alpha={alpha}, sigma = {SAT.sigma}")
        SAT.sigma = SAT.reset_sigma()
    plt.legend()
    return 0

def change_rho0():
    plt.figure()
    plt.title("Transmittance vs Time for Standard LEO")
    SAT = AnalyzeCrossing(cb="Earth", H=420, E_kev=4)
    time_array = np.arange(0, SAT.time_final + 1, 1, dtype=float)
    # list of scaling factors for rho0
    alpha_list = [0.8, 0.9, 1.0, 1.1, 1.2]
    for alpha in alpha_list:
        SAT.rho0 = alpha*SAT.rho0
        transmit_array = calc_transmit(SAT, time_array)
        plt.plot(time_array, transmit_array,
                 label=f"alpha={alpha}, rho0 = {SAT.rho0}")
        SAT.rho0 = SAT.reset_rho0()
    plt.legend()
    return 0

def change_scale_height():
    plt.figure()
    plt.title("Transmittance vs time with different scale heights")
    plt.ylabel("Transmittance")
    plt.xlabel("Time")
    plt.xlim([0, 200])
    SAT = AnalyzeCrossing(cb="Earth", H=420, E_kev=4)
    time_array = np.arange(0, SAT.time_final + 1, 1, dtype=float)
    # list of scaling factors for scale_height
    alpha_list = [0.8, 0.9, 1.0, 1.1, 1.2]
    for alpha in alpha_list:
        SAT.scale_height = alpha*SAT.scale_height
        transmit_array = calc_transmit(SAT, time_array)
        plt.plot(time_array, transmit_array,
                 label=fr"$\beta$={alpha}, scale_height = {SAT.scale_height:.2f} km")
        SAT.scale_height = SAT.reset_scale_height()
    plt.legend()
    return 0

if __name__ == '__main__':
    # plt.rc("text", usetex=True)
    # change_sigma()
    # change_rho0()
    change_scale_height()
    plt.show()




