from .spectrum import *


class Blank(Spectrum):
    """
    A very boring empty Spectrum full of nans.
    """

    def __init__(self):
        Spectrum.__init__(
            self,
            wavelength=default_wavelength_grid,
            flux=np.nan * np.ones(len(default_wavelength_grid)) * u.W / u.m**2 / u.nm,
        )
