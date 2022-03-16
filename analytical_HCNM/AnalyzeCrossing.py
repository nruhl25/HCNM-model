# Author: Nathaniel Ruhl
# This class assembles all the "tools" methods to analyze a horizon crossing

import numpy as np

# import local libraries
from Orbit import Orbit
from xsects import BCM

class AnalyzeCrossing(Orbit):

    def __init__(self, cb="Earth", H=420):
        Orbit.__init__(self, cb, H) # instansiates both Orbit and Planet classes

    ## Define functions relevant to the 2d geometry below:

    # Tangent altitude (km) as a function of elevation angle (rad)
    def tan_alt(self, t):
        h = self.R_orbit*np.sin(self.theta+self.elevation(t))-self.R
        return h

    # Relationship between the total length of the line of sight (km) and elevation angle (rad)


    def d_tot(self, t):
        dtot = 2*np.sqrt(self.R_orbit**2 - (self.R+self.tan_alt(t))**2)
        return dtot

    # Relationship between elevation angle (rad) and angular velocity (rad/sec)


    def elevation(self, t):
        epsilon = self.omega*t
        return epsilon

    # Define a functions to convert between a point at distance x on the line of sight and an altitude above Earth, z (km).


    def x_to_z(self, x_km, t):
        z = np.sqrt((self.R+self.tan_alt(t))**2+((self.d_tot(t)/2)-x_km)**2)-self.R
        return z


    def z_to_x(self, z_km, t):
        x = (self.d_tot(t)/2) - np.sqrt((self.R+z_km)**2-(self.R+self.tan_alt(t))**2)
        return x

    # Evaluate the density at a single x distance on the LOS

    def rho_vs_x(self, x_km, t):
        z_km = self.x_to_z(x_km, t)  # radial altitude above Earth
        rho = self.rho_vs_z(z_km, t)   # g/cm^3, mass density
        return rho

    # Exponential density as a function of altitude (km)
    def rho_vs_z(self, z_km, t):
        z0 = 0.0   # km, reference point for rho0
        rho0 = 0.001225  # g/cm^3, density at z0
        rho = rho0*np.exp(-(z_km-z0)/self.scale_height)
        return rho

    # Define functions for integrating a single line of sight at the time t
    # E is a single value in keV

    # This function returns an array of gamma = optical depth per km along the LOS, corresponding to the input a_array
    def gamma_vs_x(self, E_kev, x_array_km, t):
        sigma = BCM.get_total_xsect(E_kev, mix_N=self._mix_N, mix_O=self._mix_O, mix_Ar=self._mix_Ar)  # cm^2/g, cross section
        gamma_array = sigma*self.rho_vs_x(x_array_km, t)  # cm^-1
        gamma_array = gamma_array*10**5  # km^-1
        return gamma_array

    # This function calculates optical depth for a line of sight at the time t with a numerical integral
    def tau(self, E_kev, t):
        dtot_km = self.d_tot(t)   # km, total length of the los

        # Evaluate the integral over the LOS with Simpson's rule
        # tau = sigma*dx_cm*np.sum(self.rho_vs_x(x_array_km[:-1], t)) is left-hand reimann sum
        N = 10000
        a = 0
        b = dtot_km/2
        dx_km = (b-a)/N
        dx_cm = dx_km * 10 ** 5
        x_array_km = np.arange(0, dtot_km+dx_km, dx_km)

        density_array = self.rho_vs_x(x_array_km, t)    # array of densities along the LOS, g/cm^3
        sigma = BCM.get_total_xsect(E_kev, mix_N=self._mix_N, mix_O=self._mix_O, mix_Ar=self._mix_Ar)  # cm^2/g, cross section

        s_odd = 0
        s_even = 0
        for i in range(len(x_array_km)):
            if i % 2 == 0:
                s_even += density_array[i]
            else:
                s_odd += density_array[i]
        tau = (dx_cm/3)*sigma*(density_array[0] + density_array[-1] +
                        4*s_odd + 2*s_even)   # value of integral
        return 2*tau

# Code to test the class
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    ISS = AnalyzeCrossing(cb="Earth", H=420)
    time_array = np.arange(0, ISS.time_final+1, 1)
    transmit_array = np.zeros_like(time_array)
    for i, t in enumerate(time_array):
        transmit_array[i] = np.exp(-ISS.tau(4, t))

    plt.plot(time_array, transmit_array)
    plt.show()
