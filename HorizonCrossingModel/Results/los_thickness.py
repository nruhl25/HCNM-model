# Author: Nathaniel Ruhl

# This script makes plots that provide motivation for adaptive quadrature

import numpy as np
import matplotlib.pyplot as plt
plt.rc("text", usetex=True)

from AnalyzeCrossing import AnalyzeCrossing

def main():
    ISS = AnalyzeCrossing(cb="Earth", H=420, E_kev=4.0)

    # define the array of times in a horizon crossing
    time_array = np.arange(0, ISS.time_final + 1, 1)
    # km, tangent altitudes during the horizon crossing
    h_array = ISS.tan_alt(time_array)
    # Transmission at each time during the horizon crossing
    transmit_array = np.zeros_like(time_array)

    for i, t in enumerate(time_array):
        transmit_array[i] = np.exp(-ISS.tau_gauss(t, N=10))

    plt.figure(1)
    # This is for an ISS-like orbit
    plt.title(f"Transmittance of {ISS.E_kev} keV photons vs tangent altitude (km)")
    plt.plot(h_array, transmit_array)
    plt.ylabel("Transmittance")
    plt.xlabel("Tangent altitude (km)")

    # Based on this plot above, identify the time range in which transmission curve rises
    t_range = np.arange(50, 65, 2)
    plt.figure(2)
    plt.title(fr"Optical depth per km vs $x$ coordinate on the line of sight (km)")
    for t_i in t_range:
        # Define quantities for plotting
        dtot = ISS.d_tot(t_i)
        h_i = ISS.tan_alt(t_i)
        x_array_i = np.arange(0, dtot+1, 1)
        gamma_array_i = ISS.gamma_vs_x(x_array_i, t_i)
        transmit_i = np.exp(-ISS.tau_gauss(t_i, N=10))  # transmission coefficient at the time t_i
        plt.plot(x_array_i, gamma_array_i, label=fr"t={t_i:.1f} sec, h = {h_i:.1f} km, T = {transmit_i:.1f}")

    plt.ylabel(r"Optical depth per km, $\gamma$")
    plt.xlabel("$x$ coordinate on the telescopic line of sight")
    plt.legend(prop={'size': 8})

    plt.show()
    return 0

if __name__=='__main__':
    main()
