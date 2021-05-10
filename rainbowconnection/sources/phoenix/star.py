from ..spectrum import *
from .library import read_phoenix


class Star(Spectrum):
    def __init__(
        self,
        teff=5800 * u.K,
        radius=1 * u.Rsun,
        mass=1 * u.Msun,
        metallicity=0.0,
        R=None,
        extend_wavelengths=False
    ):

        self.teff = teff
        self.radius = radius
        self.mass = mass
        self.logg = np.log10((con.G * self.mass / self.radius ** 2).to("cm/s**2").value)
        self.metallicity = metallicity

        w, f = read_phoenix(
            self.teff.to("K").value,
            logg=self.logg,
            metallicity=self.metallicity,
            R=R,
            photons=False,
            extend_wavelengths=extend_wavelengths
        )
        self._wavelength = w * u.nm
        self._flux = (f * u.erg / u.s / u.cm ** 2 / u.nm).to(u.W / u.nm / u.m ** 2)
        self.default_wavelengths = self._wavelength

        # FIXME -- check if it's surface flux or something else!??!?!

    """
    def integrate(self, lower=None, upper=None):
        '''
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
        '''

        # if wavelength limits are used, revert back to the numerical integral
        if (lower is not None) or (upper is not None):
            return super().integrate(lower=lower, upper=upper)

        # if there are infinite wavelength limits, use bolometric
        return 1*u.Lsun/self.normalization()
    """
