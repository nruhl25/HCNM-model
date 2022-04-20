# Author: Nathaniel Ruhl

Venus = {
    "Mass": 4.8673 * 10 ** 24,  # kg
    "Radius": 6051.8,   # km, equatorial radius
    "surface_gravity": 8.87,  # m/s^2, acceleration due to gravity at surface level

    # Density Profile
    "surface_density": 0.065,  # g/cm^3
    "scale_height": 15.9,  # km, scale height for density exponential

    # Volumetric mix of the atmosphere, (included an extra 0.05% O)
    "mix_N": 0.035,
    "mix_O": 0.645,
    "mix_Ar": 0.0,
    "mix_C": 0.32,
}

if __name__ == '__main__':
    print(Venus["mix_N"]+Venus["mix_O"]+Venus["mix_Ar"]+Venus["mix_C"])
