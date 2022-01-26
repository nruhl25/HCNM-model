# Author: Nathaniel Ruhl
# Date: 1/26/2022
# This script calculates transmittance curves for different orbital velocities in a circular orbit, based on the ISS.
# The results enable us to quantify the "allowable error" in velocity for HCNM to a certain degree of accuracy.

# import local modules
from Python_Scripts.HCNM import TransmitOrbitModel as model
from Observations.V4641_Feb3 import feb3_dict as f3

# import standard libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


def main():
    # Using a circular model orbit for the ISS based on the V4641 Sgr horizon crossing observation on 02/03/2020
    num_decimals = 6        # 1 m and 1 m/s precision. These floats must also work as dictionary keys
    R0 = round(np.linalg.norm(f3.feb3_dict["r0_mkf_ellipsoid"]), num_decimals)    # km, radius of orbit at r0
    v0 = round(np.linalg.norm(f3.feb3_dict["v0_mkf_ellipsoid"]), num_decimals)    # km/sec, orbital velocity at r0
    omega_best = v0 / R0

    # Our goal is to vary v0 and calculate multiple transmittance curves
    v0_step = 0.1  # km/sec
    v0_array = np.arange(v0 - 0.5, v0 + 0.5 + v0_step, v0_step)
    v0_array = np.round(v0_array, num_decimals)
    omega_array = v0_array / R0   # rad/sec
    transmit_dict = {}   # key will correspond to linear velocity, value will be transmittance arrays
    tan_alt_dict = {}    # key will correspond to linear velocity, value will be tangent altitude arrays
    time_dict = {}    # same format. Time during horizon crossing

    e_band = np.array([2.0, 3.0])

    for i, omega in enumerate(omega_array):
        # Note that this simulates a "fake" Earth orbit since radius and period don't correspond via Kepler's law
        # Can make model step sizes smaller before getting final results
        model.CalculateModelTransmitVsTime.set_angular_velocity(omega)
        model.CalculateModelTransmitVsTime.set_ds_km_step(1)  # km steps along the line of sight
        model.CalculateModelTransmitVsTime.set_dx_steps_effective_transmit(1)  # keV steps within energy band
        if i == 1000:  # set to 0 if we want to see the data
            model_obj = model.CalculateModelTransmitVsTime(observation_dict=f3.feb3_dict,
                                                           binning_args=(e_band, 1),
                                                           parameter_string_args=('circle', 'ellipsoid', 'HEASARC'),
                                                           use_pymsis=True, pymsis_version=00, bin_data=True)
            transmit_dict[v0_array[i]] = model_obj.effective_transmit_3d
            transmit_data = model_obj.band_transmittance
            transmit_data_times = model_obj.band_times_binned
        else:
            model_obj = model.CalculateModelTransmitVsTime(observation_dict=f3.feb3_dict,
                                                           binning_args=(e_band, 1),
                                                           parameter_string_args=('circle', 'ellipsoid', 'HEASARC'),
                                                           use_pymsis=True, pymsis_version=00, bin_data=False)
            transmit_dict[v0_array[i]] = model_obj.effective_transmit_3d
            tan_alt_dict[v0_array[i]] = model_obj.tangent_geodetic_altitudes_pymap
            time_dict[v0_array[i]] = model_obj.t_prime
            print(f'v0={R0*model_obj.get_angular_velocity()} km/sec')   # Confirm the velocity is changing correctly

    del model_obj
    print("Finished calculating models")

    # Identify time and tangent altitude at the 50% transmission point, then compare them to the "best" model
    # First we need to reduce the domain for transmittance so it is a 1-1 function that can be interpolated
    interp_range_indices = np.where((transmit_dict[v0] > 0.01) & (transmit_dict[v0] < 0.99))[0]
    transmit_short = transmit_dict[v0][interp_range_indices]
    time_short = time_dict[v0][interp_range_indices]
    tan_alt_short = tan_alt_dict[v0][interp_range_indices]

    time_vs_transmit_best = interp1d(x=transmit_short, y=time_short, kind='cubic')
    alt_vs_transmit_best = interp1d(x=transmit_short, y=tan_alt_short, kind='cubic')
    t50_best = time_vs_transmit_best(0.5)
    alt50_best = alt_vs_transmit_best(0.5)

    t50_array = np.empty_like(v0_array)   # indices will correspond to v0_array
    alt50_array = np.empty_like(v0_array)
    # key corresponds to the linear velocity, which is the dict key
    for i, key in enumerate(transmit_dict):
        interp_range_indices = np.where((transmit_dict[key] > 0.01) & (transmit_dict[key] < 0.99))[0]
        time_vs_transmit = interp1d(x=transmit_dict[key][interp_range_indices],
                                    y=time_dict[key][interp_range_indices], kind='cubic')
        altitude_vs_transmit = interp1d(x=transmit_dict[key][interp_range_indices],
                                        y=tan_alt_dict[key][interp_range_indices], kind='cubic')
        t50_array[i] = time_vs_transmit(0.5)
        alt50_array[i] = altitude_vs_transmit(0.5)

    plot_simulations(transmit_dict, tan_alt_dict, time_dict)
    plot_t50_comparison(v0, t50_best, alt50_best, v0_array, t50_array, alt50_array)
    plt.show()

    return 0


def plot_simulations(transmit_dict, tan_alt_dict, time_dict):
    plt.figure(1)
    plt.title("Transmittance vs. Time for different velocities")
    plt.ylabel('Transmittance')
    plt.xlabel(r"Time from t_{0,e} (sec)")

    for velocity, transmit_curve in transmit_dict.items():
        plt.plot(time_dict[velocity], transmit_curve, label=f"v={velocity:.2f} km/sec")

    plt.legend()

    plt.figure(2)
    plt.title("Transmittance vs. Tangent Altitude for different velocities")
    plt.ylabel('Transmittance')
    plt.xlabel(r"Tangent Altitude (km)")

    for velocity, tan_alt in tan_alt_dict.items():
        plt.plot(tan_alt, transmit_dict[velocity], label=f"v={velocity:.2f} km/sec")

    plt.legend()

    plt.figure(3)
    plt.title("Tangent Altitude vs. Time for different velocities")
    plt.ylabel('Tangent Altitude')
    plt.xlabel(r"Time from t_{0,e}")

    for velocity, tan_alt_curve in tan_alt_dict.items():
        plt.plot(tan_alt_curve, label=f"v={velocity:.2f} km/sec")

    plt.legend()

    return


def plot_t50_comparison(v0, t50_best, alt50_best, v0_array, t50_array, alt50_array):
    plt.figure(4)
    plt.plot(v0_array - v0, t50_array - t50_best)
    plt.title(r"Error in t_{50} vs Error in linear velocity")
    plt.ylabel(r"Error in t_{50}")
    plt.xlabel("Error in velocity (km/sec)")

    plt.figure(5)
    plt.plot(v0_array - v0, alt50_array - alt50_best)
    plt.title(r"Error in alt_{50} vs Error in linear velocity")
    plt.ylabel(r"Error in alt_{50} (km)")
    plt.xlabel("Error in velocity (km/sec)")

    plt.figure(6)
    plt.plot(v0_array - v0, v0 * (t50_array - t50_best))
    plt.title(r"Simulated in-track error vs error in linear velocity")
    plt.ylabel(r"In-track error (km)")
    plt.xlabel("Error in velocity (km/sec)")

    return


if __name__ == "__main__":
    import time
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))

