from ..imports import *
from ..sources import Spectrum, Thermal
from colour.phenomena.rayleigh import rayleigh_optical_depth


class Sky(Spectrum):
    def __init__(self, sunset):
        """
        Initialize this composite object, connecting a light source
        and an atmosphere together into a diffuse Sky calculator.

        Parameters
        ----------
        sunset : rainbowconnection.Sunset
            Any object with both a source and an atmosphere.
            Calculations related to the diffuse sky brightness
            require both, so it's simplest just to pass
            a sunset object in which they are already linked.
        """

        self.sunset = sunset
        self.source = sunset.source
        self.atmosphere = sunset.atmosphere

        try:
            self.B = Thermal(self.atmosphere.T).intensity
        except AttributeError:

            def f(wavelength):
                return np.zeros(np.shape(wavelength)) * u.W / u.m**2 / u.nm * u.sr

            self.B = f

        self.set_zenith_angle(0 * u.deg)
        # self.is_earth = 'Earth' in self.atmosphere.__class__.__name__

    def set_zenith_angle(self, zenith_angle=0 * u.deg):
        """
        Set the angle from zenith (for a slice of the sky.)

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

    def spectrum(self, wavelength=None):
        """
        Calculate intensity of the diffuse sky at a given zenith angle.
        Currently this accounts only for Rayleigh scattering (assuming
        perfect 1/wavelength**4) and for thermal emission (which is
        usually negligible below about 1700K). The relative contributions
        of scattering and aborption/emission are estimated in a very
        cartoonish fashion from the overall extinction, and assumes
        there is at least one wavelength in the spectrum that is
        dominated by scattering.

        Parameters
        ----------
        wavelength : astropy.units.quantity.Quantity
            The wavelengths on which we want the spectrum.
        """

        # figure out the airmass along this line of sight
        # airmass = 1/np.cos(self.zenith_angle)
        effective_airmass = np.minimum(
            self.atmosphere.fortney_factor() / 2.0, self.airmass
        )

        # calculate the optical depths at this angle
        tau_scatter = self.atmosphere.tau_zenith_scatter * effective_airmass
        tau_absorb = self.atmosphere.tau_zenith_absorb * effective_airmass

        # estimate the single
        albedo = tau_scatter / (tau_absorb + tau_scatter)

        # calculate the mean intensity field
        mean_intensity = self.sunset.mean_intensity()

        # calculate the intensity of the illuminated sky (away from the disk)
        thermal_intensity = (
            (1 - albedo)
            * self.B(self.atmosphere.wavelength)
            * (1 - np.exp(-tau_absorb))
        )
        scattering_intensity = albedo * mean_intensity * (1 - np.exp(-tau_scatter))
        sky_intensity = thermal_intensity + scattering_intensity

        # bin this intensity onto the desired wavelengths
        w = self.get_wavelength(wavelength)
        unit = sky_intensity.unit
        neww, newi = bintogrid(
            self.atmosphere.wavelength,
            sky_intensity.value,
            newx=w.to("nm").value,
            drop_nans=False,
        )

        # make sure the wavelengths match up
        assert np.all(neww == w.to("nm").value)

        # return this sky intensity
        return newi * unit
