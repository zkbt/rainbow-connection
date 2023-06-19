from .spectrum import *
from .thermal import Thermal

__all__ = ["LightBulb", "Incandescent", "Sodium", "EmissionLines", "WhiteLED"]


def gauss(x, x0, sigma):
    """
    Simple Gaussian.
    """
    return 1.0 / np.sqrt(2 * np.pi) / sigma * np.exp(-0.5 * ((x - x0) / sigma) ** 2)


class LightBulb(Spectrum):
    def __repr__(self):
        """
        How should this object appear as a string?

        Returns
        -------
        s : str
            A simple string representation.
        """

        basic = f"{self.__class__.__name__}LightBulb ({self.power})"
        try:
            assert self.distance is not None
            return basic + f" at {self.distance}"
        except (AssertionError, AttributeError):
            return basic

    def __init__(self, wavelength=None):
        if wavelength is not None:
            self.wavelength = wavelength
        else:
            self.wavelength = default_wavelength_grid


class Incandescent(LightBulb, Thermal):
    def __init__(self, teff=3500 * u.K, power=100 * u.W, wavelength=None):

        Thermal.__init__(self, teff=teff, radius=1 * u.mm)
        LightBulb.__init__(self, wavelength=wavelength)

        # renormalize that radius to match the light output
        self.set_power(power)


class Sodium(LightBulb):
    def __init__(self, power=100 * u.W, wavelength=None):

        LightBulb.__init__(self, wavelength=wavelength)

        # start with a ridiculous radius
        self.radius = 1 * u.m

        # renormalize that radius to match the light output
        self.set_power(power)

    def surface_flux(self, wavelength):
        """
        This function calculates a cartoon spectrum of a sodium lamp.

            Inputs:
                wavelength = numpy array of wavelengths (with astropy units)
                temperature = a single number, the temperature (with astropy units)

            Outputs:
                Returns an array of thermal emission fluxes,
                in astropy units of W/(m^2*micron). This is a flux, which has
                already been integrated over solid angle.
        """

        # make sure a wavelength grid is defined
        w = self.get_wavelength(wavelength)
        flux = np.zeros(np.shape(w)) / u.nm

        # set the two lines
        lines = np.array([588.9950, 589.5924]) * u.nm
        width = 2.0 * u.nm

        # give a normalization
        brightness = 1 * u.W / u.m**2

        # add the two lines to the spectrum
        for l in lines:
            flux += gauss(w, x0=l, sigma=width)

        # return the flux, in convenient units
        return (brightness * flux).to("W/(nm * m**2)")


class EmissionLines(LightBulb):
    def __init__(
        self,
        line_centers,
        line_amplitudes,
        line_width=2 * u.nm,
        power=10 * u.W,
        wavelength=None,
    ):
        """
        Parameters
        ----------
        line_centers : array (with astropy units)
            The wavelengths of the centers of the lines.

        line_amplitudes : array (with or without units)
            The amplitudes of the lines.

        line_width : float or array (with astropy units)
            The line width, which can either be just
            one number, or an array indicating a separate
            linewidth for each line.

        power : float (with astropy units)
            The total power of the light bulb,
            integrated over all wavelengths.
        """

        LightBulb.__init__(self, wavelength=wavelength)

        # store the information about the lines
        self.line_centers = np.atleast_1d(line_centers)
        n_lines = len(self.line_centers)
        self.line_amplitudes = (
            np.atleast_1d(line_amplitudes) * u.Unit("") * np.ones(n_lines)
        )
        self.line_widths = np.atleast_1d(line_width) * np.ones(n_lines)

        # start with a ridiculous radius
        self.radius = 1 * u.m

        # renormalize that radius to match the light output
        self.set_power(power)

    def surface_flux(self, wavelength):
        """
        This function calculates a cartoon spectrum of an emission line lamp.

        Parameters
        ----------
        wavelength : array of wavelengths (with astropy units)

        Returns
        -------
        flux : array
            Emission line flux at each wavelength.
        """

        # make sure a wavelength grid is defined
        w = self.get_wavelength(wavelength)
        flux_unit = self.line_amplitudes.unit / u.nm
        flux = np.zeros(np.shape(w)) * flux_unit

        # give a normalization
        brightness = 1 * u.W / u.m**2

        # add the two lines to the spectrum
        for center, amplitude, width in zip(
            self.line_centers, self.line_amplitudes, self.line_widths
        ):
            flux += amplitude * gauss(w, x0=center, sigma=width)

        # return the flux, in convenient units
        return brightness * flux  # .to("W/(nm * m**2)")


class WhiteLED(EmissionLines):
    def __init__(self, power=100 * u.W, wavelength=None):
        """
        Parameters
        ----------
        power : float (with astropy units)
            The total power of the light bulb,
            integrated over all wavelengths.
        """
        EmissionLines.__init__(
            self,
            line_centers=np.array([450, 560]) * u.nm,
            line_amplitudes=[0.2, 1],
            line_width=np.array([15, 80]) * u.nm,
            power=power,
            wavelength=wavelength,
        )
