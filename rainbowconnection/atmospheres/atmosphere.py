from ..imports import *
from ..units import determine_quantity
from ..plottingtools import setup_axes_with_rainbow
from .sunset import Sunset


class Atmosphere:
    def __repr__(self):
        """
        How should this object appear as a string?

        Returns
        -------
        s : str
            A simple string representation.
        """
        scale = (self.H / self.radius).decompose()
        # H/R={scale},
        return f"{self.__class__.__name__}Atmosphere ({self.zenith_angle} from zenith, {self.altitude} scale heights from reference)"

    # FIXME (maybe both spectrum and atmosphere should inherit from the same thing?)
    def get_wavelength(self, wavelength=None):
        """
        A wrapper to ensure at least some grid of wavelengths
        gets defined. A default grid will be assumed, unless
        any wavelength array at all is passed.
        """

        # make sure at least some wavelengths are defined
        if wavelength is None:
            wavelength = self.wavelength
        return wavelength.to("nm")

    def plot(
        self,
        ax=None,
        wavelength=None,
        rainbow=True,
        color=None,
        style="dark_background",
        figsize=None,
        **kwargs,
    ):
        """
        A quick tool to plot the transmission through the atmosphere.

        Parameters
        ----------
        ax : matplotlib.axes._subplots.AxesSubplot
            Specify an axes object into which
            this plot should be drawn. This allows
            overplotting multiple spectra, for example
            with a structure like
                ```
                ax = Thermal(teff=5800*u.K, radius=1*u.Rsun).plot()
                ax = Sun.plot(ax)
                ```
            To overplot on current axes use `ax=plt.gca()`.

        wavelength : astropy.units.quantity.Quantity
            A grid of wavelengths on which the spectrum should
            be plotted. If None, the function defaults to
            covering visible wavelengths at 1nm resolution.

        rainbow : bool
            Should we add an extra rainbow above the plot,
            to indicate how wavelengths match to visible light?

        color : str
            The color for drawing the spectrum.
            'auto' represents the actual visible color.

        """

        # set up to use a dark background for the plot; make sure units match
        with plt.style.context(style), quantity_support():

            # setup the basic axes
            ax = setup_axes_with_rainbow(ax=ax, rainbow=rainbow, figsize=figsize)

            # make sure at least some wavelengths are defined
            w = self.get_wavelength(wavelength)

            # pull out the spectrum
            t = self.transmission(w)

            plt.plot(w, t * 100, color=color, label=self, **kwargs)

            # add the axis labels
            wunit = w.unit.to_string("latex_inline")
            plt.xlabel(f"Wavelength ({wunit})")
            plt.ylabel(f"Transmission (%)")

            plt.xlim(np.min(w).value, np.max(w).value)

        return ax

    def set_elevation_angle(self, elevation=90 * u.deg):
        """
        Set the elevation angle above the horizon.

        Parameters
        ----------
        elevation : astropy.units.quantity.Quantity
            The angle above the horizon along which
            the transmission of the atmosphere should
            be calculated.
        """
        self.set_zenith_angle(90 * u.deg - elevation)

    def set_zenith_angle(self, zenith_angle=0 * u.deg):
        """
        Set the angle from zenith.

        Parameters
        ----------
        zenith_angle : astropy.units.quantity.Quantity
            The angle away from zenith along which
            the transmission of the atmosphere should
            be calculated.
        """

        # set the zenith angle and airmass
        self.zenith_angle = np.minimum(zenith_angle, 90 * u.deg)
        self.airmass = 1 / np.cos(self.zenith_angle)

    def set_altitude(self, altitude=0):
        """
        Set how many atmospheric scale heights we are
        above or below the reference radius of this
        atmosphere.

        Parameters
        ----------
        altitude : float
            The altitude at which we should be floating
            in the atmosphere, in units of scale heights.
            For example, +1.0 is one scale height above
            the reference radius and therefore a more
            transparent atmosphere, and -1.0 is one scale
            height deeper than the reference radius and
            therefore a more opaque atmosphere.
            Negative altitudes are probably a weird concept
            if the reference radius is the surface of a
            rocky planet, but they're quite reasonable
            on gas-dominated planets.
        """
        self.altitude = altitude
        self.tau_zenith = self._tau_zenith_reference * np.exp(-altitude)

        # guess what is scattering vs. extinction
        self.guess_scattering()

    def transmit(self, spectrum):
        """
        Calculate the spectrum of light that results from
        transmitting a known light source through this
        atmosphere.

        Parameters
        ----------
        spectrum : Spectrum
            Any light source that you want to transmit.
        """

        # create a composite object, linking the source and the atmosphere
        return Sunset(source=spectrum, atmosphere=self)


class DiscreteAtmosphere(Atmosphere):
    """
    The DiscreteAtmosphere represents an atmosphere
    whose transmission as a function of wavelength
    can be represented as a grid of zenithward optical
    depths. Most pre-calculated models will fall
    into this category.
    """

    def __init__(self, zenith_angle=0.0 * u.deg, altitude=0.0, **kwargs):
        """
        In inherited classes, this initialization relies on the
        existence of a method called ".read_transmission"
        that will create the attributes:
            ._wavelengths (with units of wavelength)
            ._tau_zenith_reference (unitless)
            .H (with units of length)
            .radius (with units of length)
        """

        # read the transmission spectrum data
        self.read_transmission(**kwargs)
        self.wavelength = self._wavelength

        # set the zenith angle (or fall back to the current setting)
        self.set_zenith_angle(zenith_angle)

        # set the altitude at which we're hovering
        self.set_altitude(altitude)

    def guess_scattering(self, visualize=False):
        """
        Make a guess at the Rayleigh scattering component of
        the atmospheric extinction, for use in estimating
        the sky color.

        This should be called once per every new altitude.
        """

        # figure out the appropriate normalization
        self._rayleigh_normalization = np.min(self.tau_zenith * self.wavelength**4)

        w = self.wavelength
        self.tau_zenith_scatter = self._rayleigh_normalization / w**4
        self.tau_zenith_absorb = self.tau_zenith - self.tau_zenith_scatter

        if visualize:

            plt.plot(w, self.tau_zenith, zorder=100, label="total")
            plt.plot(
                w,
                self.tau_zenith_scatter,
                alpha=0.5,
                label="scattering",
            )
            plt.plot(
                w,
                self.tau_zenith_absorb,
                alpha=0.5,
                label="absorbtion",
            )
            plt.legend()
            plt.xscale("log")
            plt.yscale("log")
            plt.xlabel("Wavelength ({})".format(w.unit.to_string("latex_inline")))
            plt.ylabel(r"$\tau_{zenith}$")

    """
    def tau_zenith_scatter(self):
        '''
        Estimate the optical depth to *absorption*.
        (at the default wavelengths of the grid)
        '''
        w = self.wavelength
        return self._rayleigh_normalization/w**4

    def tau_zenith_absorb(self):
        '''
        Estimate the optical depth to *absorption*.
        (at the default wavelengths of the grid)
        '''
        # calculate the extinction from absorption
        return self.tau_zenith - self.tau_zenith_scatter()
    """

    def fortney_factor(self):
        """
        The ratio of optical depth between a horizontal slant path
        through an atmosphere to a vertical path, accounting
        for the sphericity of the planet. This is a simple approach
        to connect exoplanet transmission spectra to vertical
        optical depths through an atmosphere (see Fortney 2005).

        Returns
        -------
        fortney_factor : float
            The ratio of slant to vertical optical depths.
        """

        # calculate the fortney factor = ratio of slant/vertical optical depths
        return np.sqrt(2 * np.pi * self.radius / self.H).decompose()

    def transit_radius(self, wavelength=None):
        """
        Calculate the effective transit radius of the planet's atmosphere,
        as seen from infinitely far away.

        Parameters
        ----------
        wavelength : astropy.units.quantity.Quantity
            The wavelengths on which we want the transit radius.

        Returns
        -------
        transit_radius : astropy.units.quantity.Quantity
            The projected radius of the planet, with units.
        """
        # make sure at least some grid of wavelengths is defined
        w = self.get_wavelength(wavelength)

        # figure out the slant optical depth at the reference radius
        effective_airmass = self.fortney_factor()
        tau = self._tau_zenith_reference * effective_airmass

        # bin this spectrum to the particular wavelength grid
        # KLUDGE -- this should really be in transmission space!
        neww, tau_slant = bintogrid(
            self._wavelength.to("nm").value,
            tau,
            newx=w.to("nm").value,
            drop_nans=False,
        )

        # make sure the wavelengths match up
        assert np.all(neww == w.to("nm").value)

        # convert back to tau
        z_over_H = np.log(tau_slant)

        return self.radius + self.H * z_over_H

    def transmission(self, wavelength=None, zenith_angle=None, altitude=None):
        """
        Calculate the transmission through the atmosphere.

        Parameters
        ----------
        wavelength : astropy.units.quantity.Quantity
            The wavelengths on which we want the transmission.

        zenith_angle : astropy.units.quantity.Quantity
            The angle from zenith along which transmission
            should be calculated.

        altitude : float
            The altitude at which we should be floating
            in the atmosphere, in units of scale heights.

        Returns
        -------
        transmission : numpy.ndarray
            The fractional transmission through the atmosphere.
        """

        # update the zenith angle, if necessary
        if zenith_angle is not None:
            self.set_zenith_angle(zenith_angle)

        # update the altitude, if necessary
        if altitude is not None:
            self.set_altitude(altitude)

        # make sure at least some grid of wavelengths is defined
        w = self.get_wavelength(wavelength)

        # figure out the transmission at this altitude
        # FIXME -- this is a major kludge! do the integral!
        effective_airmass = np.minimum(self.fortney_factor() / 2.0, self.airmass)
        tau = self.tau_zenith * effective_airmass

        # bin this spectrum to the particular wavelength grid
        # FIXME: binning choice gets real scary with transmission
        neww, newt = bintogrid(
            self._wavelength.to("nm").value,
            np.exp(-tau),
            newx=w.to("nm").value,
            drop_nans=False,
        )

        # make sure the wavelengths match up
        assert np.all(neww == w.to("nm").value)

        # return the binned transmission
        return newt
