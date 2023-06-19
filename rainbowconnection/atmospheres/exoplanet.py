from .atmosphere import *


class Exoplanet(DiscreteAtmosphere):
    def __init__(
        self,
        wavelength,
        planet_radius,
        H=20 * u.km,
        altitude=0.0,
        zenith_angle=0.0 * u.deg,
        **kwargs
    ):
        """
        Initialize a generic exoplanet atmosphere from a model exoplanet
        transmission spectrum, specified by the inputs `wavelength`
        and `planet_radius`

        Parameters
        ----------
        wavelength : astropy.units.quantity.Quantity
            The wavelengths of the transmission spectrum,
            in units equivalent to [microns].
        planet_radius : astropy.units.quantity.Quantity
            The wavelength-dependent radius of the planet,
            in units equivalent to [m].
        H : astropy.units.quantity.Quantity
            The scale height of the atmosphere. This is used
            only for estimating how to transform between
            slant and zenith viewing geometries.
        altitude : float
            The altitude at which we should be floating
            in the atmosphere, in units of scale heights.
        zenith_angle : astropy.units.quantity.Quantity
            The default angle from zenith along which transmission
            should be calculated (this can be changed).

        """

        # define the wavelength grid
        self._wavelength = wavelength.to("nm")

        # store the effective transit radius
        self._transit_radius = planet_radius

        # pick an (arbitrary-ish) reference radius
        reference_radius = np.min(self._transit_radius)

        # define geometry of the atmosphere (how spherical?)
        self.H = H
        self.radius = reference_radius

        # calculate the optical depth at zenith
        z_over_H = ((self._transit_radius - reference_radius) / self.H).decompose()
        tau_slant = np.exp(z_over_H)
        self._tau_zenith_reference = (tau_slant / self.fortney_factor()).decompose()

        # define the default wavelength grid to use
        self.wavelength = self._wavelength

        # set the zenith angle (or fall back to the current setting)
        self.set_zenith_angle(zenith_angle)

        # set the altitude at which we're hovering
        self.set_altitude(altitude)
