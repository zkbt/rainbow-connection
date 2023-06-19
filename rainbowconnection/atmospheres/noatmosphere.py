from .atmosphere import *


class NoAtmosphere(Atmosphere):
    """
    A fake atmosphere that represents a completely transparent atmosphere.
    This might be useful as a placeholder, if you want to make
    plots/animations that directly compare the spectrum of a
    source with or without extinction by an atmosphere.
    """

    def __repr__(self):
        """
        How should this object appear as a string?

        Returns
        -------
        s : str
            A simple string representation.
        """
        # H/R={scale},
        return f"{self.__class__.__name__}Atmosphere"

    def __init__(self, zenith_angle=0.0 * u.deg, **kwargs):
        """
        Initialize the fake atmosphereless atmosphere.
        """

        # read the transmission spectrum data
        self.wavelength = default_wavelength_grid

        # set the zenith angle (or fall back to the current setting)
        self.set_zenith_angle(zenith_angle)

    def transmission(self, wavelength=None, zenith_angle=None, **kw):
        """
        Calculate the transmission through the atmosphere.

        Parameters
        ----------
        wavelength : astropy.units.quantity.Quantity
            The wavelengths on which we want the transmission.

        zenith_angle : astropy.units.quantity.Quantity
            The angle from zenith along which transmission
            should be calculated.

        Returns
        -------
        transmission : numpy.ndarray
            The fractional transmission through the atmosphere.
        """

        # update the zenith angle, if necessary
        if zenith_angle is not None:
            self.set_zenith_angle(zenith_angle)

        # make sure at least some grid of wavelengths is defined
        w = self.get_wavelength(wavelength)

        # return 100% transmission at all wavelengths
        return np.ones(np.shape(w))
