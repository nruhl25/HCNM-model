# Author: Nathaniel Ruhl

Mars = {
    "Mass": 0.64169 * 10 ** 24,  # kg
    "Radius": 3396.2,   # km, equatorial radius
    "surface_gravity": 3.71,  # m/s^2, acceleration due to gravity at surface level

    # Density Profile
    "surface_density": 0.00002,  # g/cm^3, 
    "scale_height": 11.1,  # km, scale height for density exponential

    # Volumetric mix of the atmosphere, include 0.951 CO2, (0.06 CO not included)
    "mix_N": 0.03,
    "mix_O": 0.63,
    "mix_Ar": 0.02,
    "mix_C": 0.32
}

if __name__ == '__main__':
    tot = Mars["mix_N"]+Mars["mix_O"]+Mars["mix_Ar"]+Mars["mix_C"]
    print(tot)
