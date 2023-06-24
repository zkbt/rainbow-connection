from .spectrum import *


class Sun(Spectrum):
    def __init__(self):
        self.radius = 1 * u.Rsun
        self.teff = (
            (u.Lsun / (con.sigma_sb * 4 * np.pi * u.Rsun**2)) ** (1.0 / 4.0)
        ).to(u.K)

        filename = os.path.join(data_directory, "solarspectrum.txt")
        d = ascii.read(filename, comment="#")
        self._wavelength = d["wavelength"].data * u.nm
        factor = ((1 * u.au) ** 2 / (1 * u.Rsun) ** 2).decompose()
        self._flux = d["irradiance"].data * u.W / u.nm / u.m**2 * factor
        self.wavelength = self._wavelength

    def integrate(self, lower=None, upper=None):
        """
        Integrate the spectrum over wavelength.

        It gives a number with units identical to the results of
        `.spectrum()` but without the wavelength (W or W/m**2).

        Parameters
        ----------
        lower : astropy.units.quantity.Quantity
            The lower wavelength limit.

        upper : astropy.units.quantity.Quantity
            The lower wavelength limit.

        Returns
        -------
        integral : astropy.units.quantity.Quantity
            The integral over wavelength.
        """

        # if wavelength limits are used, revert back to the numerical integral
        if (lower is not None) or (upper is not None):
            return super().integrate(lower=lower, upper=upper)

        # if there are infinite wavelength limits, use bolometric
        return 1 * u.Lsun / self.normalization()
