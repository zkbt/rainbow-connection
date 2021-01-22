from .spectrum import *


class Blank(Spectrum):
    def __init__(self):
        Spectrum.__init__(
            self,
            wavelength=np.linspace(400, 700, 100) * u.nm,
            flux=np.nan * np.ones(100),
        )
