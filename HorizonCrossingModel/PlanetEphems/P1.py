# Author: Nathaniel Ruhl
# This is a made-up planet with the size of Mars, but the atmosphere of Earth

P1 = {
    # General Info
    "Mass": 0.64169 * 10 ** 24,  # kg
    "Radius": 3396.2,   # km, equatorial radius
    "surface_gravity": 3.71,  # m/s^2, acceleration due to gravity at surface level

    # Density Profile
    "surface_density": 0.001225,  # g/cm^3
    "scale_height": 8.5,  # km, scale height for density exponential

    # Volumetric mix of the atmosphere
    "mix_N": 0.78,
    "mix_O": 0.21,
    "mix_Ar": 0.01,
    "mix_C": 0.0
}

if __name__ == '__main__':
    print(P1)
