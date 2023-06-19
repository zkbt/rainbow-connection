from .spectrum import *
from chromatic import get_phoenix_photons


class Star(Spectrum):
    def __init__(
        self,
        teff=5800 * u.K,
        radius=1 * u.Rsun,
        mass=1 * u.Msun,
        metallicity=0.0,
        R=100,
        wavelength=None,
    ):

        self.teff = teff
        self.radius = radius
        self.mass = mass
        self.logg = np.log10((con.G * self.mass / self.radius**2).to("cm/s**2").value)
        self.metallicity = metallicity
        if wavelength is not None:
            R = None

        model_wavelength, model_photon_flux = get_phoenix_photons(
            temperature=self.teff.to_value("K"),
            logg=self.logg,
            metallicity=self.metallicity,
            R=R,
            wavelength=wavelength,
        )
        self._wavelength = model_wavelength.to("nm")
        energy_per_photon = (con.h * con.c / model_wavelength).to("J") / u.photon
        surface_flux = model_photon_flux * energy_per_photon
        self._flux = surface_flux
        self.wavelength = self._wavelength

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
