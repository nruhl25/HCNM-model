# Author: Nathaniel Ruhl

# This script uses Newton's method to solve for atmospheric scale height, L

from AnalyzeCrossing import AnalyzeCrossing
import numpy as np
import random
import matplotlib.pyplot as plt
# plt.rc("text", usetex=True)  # uncover only for plots for paper


STD = 0.05  # Standard deviation of normal distribution from which noise is generated for generating the transmittance 'data'

# range in which noise is added to the model and in which rho0 is solved
COMP_RANGE = [0.01, 0.9]

# This function generates noisy data for a horizon crossing


def generate_crossing(SAT, plot_bool):
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
    
    if plot_bool is True:
        plt.figure()
        plt.title("Simulated horizon crossing data")
        plt.ylabel("Transmittance")
        plt.xlabel("Time (s)")
        plt.plot(time_array, transmit_data, ".")

    return transmit_data

# This is the function used for Newton's method

def f(SAT, transmit_data_i, time_i, L_guess):
    return np.log(transmit_data_i) + 2*SAT.sigma*SAT.rho0*SAT.exp_integral(time_i, scale_height=L_guess)

# This function employs Newton's method to solve for the scale height (L) for which the transmission equation is valid, L_guess is the first guess for L. transmit_model_i is needed for chisq analysis
# Returns the root, L, and chisq list track the iterations of Newton's method

def find_root(SAT, transmit_data_i, time_i, L_guess, transmit_model_i):
    L = L_guess
    L_last = L_guess - 0.1
    delta = 1.0
    accuracy = 1e-4

    L_error_list = []
    chisq_list = []
    while abs(delta) > accuracy:
        L_error_list.append(abs(L - SAT.scale_height))
        b = f(SAT, transmit_data_i, time_i, L)
        m = (f(SAT, transmit_data_i, time_i, L) -
             f(SAT, transmit_data_i, time_i, L_last))/(L-L_last)
        delta = b/m
        chisq = ((transmit_data_i-transmit_model_i)/STD)**2
        chisq_list.append(chisq)
        L_last = L
        L -= delta
    return L, np.array(L_error_list), np.array(chisq_list)

# Function to solve rho0 if cross section and scale height are known


def solve_L(SAT, transmit_data, L0_guess, crossing_plot_bool, chisq_plot_bool):
    # Calculate tha model in order to specify the solution range
    time_array = np.arange(0, SAT.time_final + 1, 1, dtype=float)
    transmit_model = np.zeros_like(time_array)
    for i, t in enumerate(time_array):
        transmit_model[i] = np.exp(-SAT.tau_gauss(t, N=10))
    # Index range in which to solve for rho0
    sol_range = np.where((transmit_model > COMP_RANGE[0]) & (
        transmit_model < COMP_RANGE[1]))[0]

    # Use Newton's method to solve for scale height, L
    L_list = []  # measured rho0 for each data point in the desired range
    for indx in sol_range:
        L, L_error_list, chisq_list = find_root(
            SAT, transmit_data[indx], time_array[indx], L0_guess, transmit_model[indx])
        L_list.append(L)
        print(f"num iterations = {len(chisq_list)}")
        if chisq_plot_bool is True:
            # Plots to look more into the convergence of Newton's method
            plt.figure()
            plt.title("Iterations of Newton's Method")
            plt.ylabel(r"$\chi^2$")
            plt.xlabel(r"$|L_{{guess}}-L|$")
            plt.plot(chisq_list, L_error_list, '.', label=fr"Transmittance = {transmit_model[indx]:.2f}")
            plt.legend()

    L_list = np.array(L_list)
    if crossing_plot_bool is True:
        plt.figure()
        plt.title(
            f"Scale height of {SAT.cb} measured from a horizon crossing")
        plt.ylabel(r"Scale height (km)")
        plt.xlabel("Transmittance of Data Point")
        plt.plot(transmit_data[sol_range], L_list,
                 label=fr'Mean $L$={np.mean(L_list):.4f} km')
        plt.axhline(y=SAT.scale_height, xmin=0, xmax=1, linestyle='--',
                    color='r', label=fr"Expected $L=$ {SAT.scale_height} km")
        plt.legend()

        # Plot percent fractional difference from expected value
        diff = abs(SAT.scale_height-L_list)/SAT.scale_height
        plt.figure()
        plt.title(
            r"Fractional error of scale height at different moments in the horizon crossing")
        plt.xlabel("Transmittance of Data Point")
        plt.ylabel(r"Fractional Error of scale height")
        plt.plot(transmit_data[sol_range], diff,
                 label=fr"Mean Calculated $L=${np.mean(L_list):.4f} km")
        plt.plot([], [], label=fr"Expected $L=${SAT.scale_height} km")
        plt.legend()

    return np.mean(L_list)

def main():
    SAT = AnalyzeCrossing(cb="Earth", H=420, E_kev=4.0)
    transmit_data = generate_crossing(SAT, plot_bool=True)
    L_mean = solve_L(SAT, transmit_data, L0_guess=SAT.scale_height+1, crossing_plot_bool=True, chisq_plot_bool=False)
    plt.show()
    return 0


if __name__ == '__main__':
    main()
