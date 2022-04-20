# This file plots percent contribution to absorption along the LOS. It uses two modes of comparison, contribution as defined by contribution to total optical depth and contribution as defined by contribution to total absorption

import numpy as np
import matplotlib.pyplot as plt

# import local libraries
from AnalyzeCrossing import AnalyzeCrossing
from gaussxw import gaussxwab

# Define the satellite to be used
SAT = AnalyzeCrossing(cb="Earth", H=420, E_kev=4.0)

# Calculates total optical depth along the LOS at time t. Called from the function that calculates percent contribution
def calc_total_optical_depth(t):
    dtot = SAT.d_tot(t)
    N = 500
    a = 0.0
    b = dtot/2
    h = (b-a)/N
    xlist, wlist = gaussxwab(N, a, b)
    gamma_array = SAT.gamma_vs_x(xlist, t)   # Optical depth per km
    # Integrate gamma vs x with gaussian quadrature
    tau_total = 2*np.sum(wlist*gamma_array)
    return tau_total

# Note that tau_total is calculated in the above function and used as an input to this to save computation time. "comp_string" determines if we're comparing optical depth or absorptio
def calc_percent_contribution(t, x1, tau_total, comp_string):
    dtot = SAT.d_tot(t)
    N = 500
    a = 0.0
    b = dtot/2
    h = (b-a)/N
    # Calculate the area under the specified region of the curve
    xlist, wlist = gaussxwab(N, x1, b)
    gamma_region = SAT.gamma_vs_x(xlist, t)
    # optical depth within the specified portion of the curve
    tau1 = 2*np.sum(wlist*gamma_region)
    if comp_string == "tau":
        contribution = tau1/tau_total
    elif comp_string == "absorption":
        A1 = 1 - np.exp(-tau1)   # Absorption within specified range
        Atot = 1 - np.exp(-tau_total)   # Absorption on total LOS
        contribution = A1/Atot
    else:
        raise RuntimeError("Invalid Argument: 'comp_string' must be either 'tau' or 'absorption'")
    return contribution

# This function looks at percent contribution on different parts of the LOS vs time.
# comp_string="tau" or "absorptio"
def contribution_vs_time(comp_string):
    time_list = np.arange(50, 65, 2)
    for ti in time_list:
        dtot_i = SAT.d_tot(ti)
        x1_list_i = np.linspace(dtot_i/4, dtot_i/2, 100)
        tau_total = calc_total_optical_depth(ti)  # used as an input to calc_percent_contribution()
        contribution_list = []   # contribution list for a single time
        for x1 in x1_list_i:
            contribution_list.append(calc_percent_contribution(ti, x1, tau_total, comp_string))
        # X axes for plots
        xlist = (dtot_i/2)-x1_list_i
        zlist = SAT.x_to_z(xlist, ti)
        plt.figure(1)
        plt.plot(xlist, contribution_list, label=fr"t={ti} sec, h={SAT.tan_alt(ti):.2f}")
        plt.figure(2)
        plt.plot(zlist, contribution_list, label=fr"t={ti} sec, h={SAT.tan_alt(ti):.2f}")

    plt.figure(1)
    plt.ylabel("Contribution")
    plt.xlabel(r"$\pm$ x km from the tangent point")
    plt.title(
        fr"Relative contribution to total {comp_string} (E={SAT.E_kev} keV), scale height = {SAT.scale_height} km)")
    plt.legend()

    plt.figure(2)
    plt.ylabel("Contribution")
    plt.xlabel(r"Altitude above Earth, $z$ (km)")
    plt.title(fr"Relative contribution to total {comp_string} (E={SAT.E_kev} keV, scale height = {SAT.scale_height} km)")
    plt.legend()

    plt.show()
    return 0

# In order for this curve to be smooth, we must select a time that is in the "interesting time range" for all scale heights
def contribution_vs_scale_height(comp_string):
    t = 50
    dtot = SAT.d_tot(t)
    x1_list = np.linspace(dtot/4, dtot/2, 100)
    # Lists for plot x axes
    xlist = (dtot/2)-x1_list
    zlist = SAT.x_to_z(x1_list, t)
    for scale_height in [6, 7, 8, 9]:
        SAT.scale_height = scale_height
        # used as an input to calc_percent_contribution()
        tau_total = calc_total_optical_depth(t)
        contribution_list = []   # contribution list for a single time
        for x1 in x1_list:
            contribution_list.append(
                calc_percent_contribution(t, x1, tau_total, comp_string))
        plt.figure(1)
        plt.plot(xlist, contribution_list,
                    label=fr"scale height = {SAT.scale_height} km")
        plt.figure(2)
        plt.plot(zlist, contribution_list,
                 label=fr"scale height = {SAT.scale_height} km")
    plt.figure(1)
    plt.ylabel("Contribution")
    plt.xlabel(r"$\pm$ x km from the tangent point")
    plt.title(
        fr"Relative Contribution to total {comp_string} (E={SAT.E_kev} keV, t={t} sec, h={SAT.tan_alt(t):.2f})")
    plt.legend()

    plt.figure(2)
    plt.ylabel("Contribution")
    plt.xlabel(r"Altitude, $z$ (km)")
    plt.title(
        fr"Relative Contribution to total {comp_string} (E={SAT.E_kev} keV, t={t} sec, h={SAT.tan_alt(t):.2f})")
    plt.legend()
    plt.show()
    return 0

if __name__ == '__main__':
    plt.rc("text", usetex=True)
    # contribution_vs_scale_height("absorption")
    contribution_vs_time("absorption")
