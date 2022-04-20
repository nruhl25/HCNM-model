# Author: Nathaniel Ruhl

# This script uses Newton's method to solve for sigma, rho0, or scale_height from a horizon crossing

import numpy as np
import random
import matplotlib.pyplot as plt
plt.rc("text", usetex=True) # uncover only for plots for paper

from AnalyzeCrossing import AnalyzeCrossing

STD = 0.05  # Standard deviation of normal distribution from which noise is generated for generating the transmittance 'data'

COMP_RANGE = [0.01, 0.9]  # range in which noise is added to the model and in which rho0 is solved

# This function generates noisy data for a horizon crossing
def generate_crossing(SAT):
    time_array = np.arange(0, SAT.time_final + 1, 1, dtype=float)
    transmit_model = np.zeros_like(time_array)
    transmit_data = np.zeros_like(time_array)
    for i, t in enumerate(time_array):
        transmit_model[i] = np.exp(-SAT.tau_gauss(t, N=10))

    for i, t in enumerate(time_array):
        if (transmit_model[i] > COMP_RANGE[0]) & (transmit_model[i] < COMP_RANGE[1]):
            transmit_data[i] = transmit_model[i]*np.random.normal(1, STD)
        else:
            transmit_data[i] = transmit_model[i]
    '''
    plt.figure()
    plt.title("Simulated horizon crossing data")
    plt.ylabel("Transmittance")
    plt.xlabel("Time (s)")
    plt.plot(time_array, transmit_data, ".")
    '''
    return transmit_data

# This is the function used for Newton's method
def f(SAT, transmit_data_i, time_i, rho0):
    return np.log(transmit_data_i) + 2*SAT.sigma*rho0*SAT.exp_integral(time_i)

# This function employs Newton's method to solve for the rho0 for which the transmission equation is valid, rho0_guess is the first guess for rho0. transmit_model_i is needed for chisq analysis
# Returns the root, and rho0 and chisq lists that track the iterations of Newton's method

def find_root(SAT, transmit_data_i, time_i, rho0_guess, transmit_model_i):
    rho0 = rho0_guess
    rho0_last = rho0_guess - 0.1
    delta = 1.0
    accuracy = 1e-9

    rho0_error_list = []
    chisq_list = []
    while abs(delta) > accuracy:
        rho0_error_list.append(abs(rho0 - SAT.rho0))
        b = f(SAT, transmit_data_i, time_i, rho0)
        m = (f(SAT, transmit_data_i, time_i, rho0) -
             f(SAT, transmit_data_i, time_i, rho0_last))/(rho0-rho0_last)
        delta = b/m
        chisq = ((transmit_data_i-transmit_model_i)/STD)**2
        chisq_list.append(chisq)
        rho0_last = rho0
        rho0 -= delta
    return rho0, np.array(rho0_error_list), np.array(chisq_list)

# Function to solve rho0 if cross section and scale height are known
def solve_rho0(SAT, transmit_data, plot_bool):
    # Calculate tha model in order to specify the solution range
    time_array = np.arange(0, SAT.time_final + 1, 1, dtype=float)
    transmit_model = np.zeros_like(time_array)
    for i, t in enumerate(time_array):
        transmit_model[i] = np.exp(-SAT.tau_gauss(t, N=10))
    # Index range in which to solve for rho0
    sol_range = np.where((transmit_model > COMP_RANGE[0]) & (transmit_model < COMP_RANGE[1]))[0]

    # Use Newton's method to solve for rho0
    rho0_list = []  # measured rho0 for each data point in the desired range
    for indx in sol_range:
        rho0, rho0_error_list, chisq_list = find_root(
            SAT, transmit_data[indx], time_array[indx], 100*SAT.rho0, transmit_model[indx])
        rho0_list.append(rho0)
        # Plots to look more into the convergence of Newton's method
        # plt.figure()
        # plt.title("Iterations of Newton's Method")
        # plt.ylabel(r"$\chi^2$")
        # plt.xlabel(r"|$\rho_{{0,guess}}-\rho_0$|")
        # plt.plot(chisq_list, rho0_error_list, '.', label=fr"Transmittance = {transmit_model[indx]:.2f}")
        # plt.legend()
        # print(f"num iterations = {len(chisq_list)}")
    
    rho0_list = np.array(rho0_list)
    if plot_bool is True:
        plt.figure()
        plt.title(f"Surface-level density of {SAT.cb} measured from a horizon crossing")
        plt.ylabel(r"Density (g/cm$^3$)")
        plt.xlabel("Transmittance of Data Point")
        plt.plot(transmit_data[sol_range], rho0_list,
                label=fr'Mean $\rho_0$={np.mean(rho0_list):.6f} g/cm$^3$')
        plt.axhline(y=SAT.rho0, xmin=0, xmax=1, linestyle='--', color='r', label=fr"Expected $\rho_0=$ {SAT.rho0} g/cm$^3$")
        plt.legend()

        # Plot percent fractional difference from expected value
        diff = abs(SAT.rho0-rho0_list)/SAT.rho0
        plt.figure()
        plt.title(r"Fractional error of $\rho_0$ at different moments in the horizon crossing")
        plt.xlabel("Transmittance of Data Point")
        plt.ylabel(r"Fractional Error of $\rho_0$")
        plt.plot(transmit_data[sol_range], diff)
    
    return np.mean(rho0_list)

# This input determines if we're looking at a single horizon crossing or 50 in order to get determine the mean rho0 value calcuted
def main(mean_bool):
    SAT = AnalyzeCrossing(cb="Earth", H=420, E_kev=4.0)
    transmit_data = generate_crossing(SAT)
    if mean_bool is True:
        rho0_mean_list = []
        for i in range(50):
            rho0_mean = solve_rho0(SAT, transmit_data, plot_bool=False)
            rho0_mean_list.append(rho0_mean)
        
        rho0_mean_list = np.array(rho0_mean_list)
        print(f"rho0 mean over 50 HC's: {np.mean(rho0_mean_list)} g/cm^3")
    else:
        rho0_mean = solve_rho0(SAT, transmit_data, plot_bool=True)
        print(f"rho0 mean over 1 HC: {rho0_mean} g/cm^3")
        plt.show()
    return 0

if __name__ == '__main__':
    main(mean_bool=False)
