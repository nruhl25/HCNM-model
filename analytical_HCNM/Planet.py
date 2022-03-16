# Author: Nathaniel Ruhl
# The class in this script defines parameters of the central body
# Can also contain information about the atmosphere

G = 6.6743*10**(-11)     # Nm^2/kg^2, Gravitational constant

# For now the default central body "cb" is Earth, and others can be added in the future
class Planet():
    def __init__(self, cb="Earth"):
        if cb != "Earth":
            raise RuntimeError("Planet class is only defined for cb='Earth'")
        self.cb = cb
        self.M = 5.972 * 10 ** 24   # kg, mass of Earth
        self.R = 6378.137  # km, semi-major axis of Earth (equatorial radius)

        # Default properties of the Planet's atmosphere defined below... can be changed internally with the setter methods
        self._mix_N = 0.78  # default mixing ratios of atmospheric consituents
        self._mix_O = 0.21 
        self._mix_Ar = 0.01
        self._scale_height = 8  # km, default atmospheric scale height for exponential density model

    @property
    def mix_N(self):
        return self._mix_N

    @mix_N.setter
    def mix_N(self, mix_N):
        self._mix_N = mix_N

    @property
    def mix_O(self):
        return self._mix_O

    @mix_O.setter
    def mix_O(self, mix_O):
        self._mix_O = mix_O

    @property
    def mix_Ar(self):
        return self._mix_Ar

    @mix_Ar.setter
    def mix_Ar(self, mix_Ar):
        self._mix_Ar = mix_Ar

    @property
    def scale_height(self):
        return self._scale_height

    @scale_height.setter
    def scale_height(self, scale_height):
        self._scale_height = scale_height


if __name__ == '__main__':
    Earth = Planet(cb="Earth")
    Earth.mix_N = 0.5
    print(Earth.mix_N)
    print(Earth.scale_height)
    Earth.scale_height = 6
    print(Earth.scale_height)


