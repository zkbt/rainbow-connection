from ..imports import *
from ..sources import Spectrum
from colour.phenomena.rayleigh import rayleigh_optical_depth

class Sky(Spectrum):

    def __init__(self, sunset):
        '''
        Initialize this composite object, connecting a light source
        and an atmosphere together into a diffuse Sky calculator.
        '''

        self.sunset = sunset
        self.source = sunset.source
        self.atmosphere = sunset.atmosphere

    def set_zenith_angle(self, z=0*u.deg):
        '''
        Set the angle from zenith.

        Parameters
        ----------
        zenith : astropy.units.quantity.Quantity
            The angle away from zenith along which
            the transmission of the atmosphere should
            be calculated.
        '''

        # set the zenith angle and airmass
        self.zenith_angle = np.minimum(z, 90*u.deg)
        self.airmass = 1/np.cos(self.zenith_angle)

    def tau_rayleigh_scatter(self, wavelength=None):
        '''
        (Crudely) estimate the optical depth for Rayleigh scattering.

        This scales the optical depth with the airmass and with the altitude
        of the viewing location within the atmosphere, but it is ultimately
        assuming an Earth-like composition.

        The goal is to provide an approximate means to estimate the
        diffuse brightness of the sky, away from the disk of
        any particular light source.
        '''

        # convert the wavelength to cm (without units)
        w = self.atmosphere.wavelength(wavelength)
        w_cm = w.to('cm').value

        # calculate the optical depth to Rayleigh scattering vertically, from sea level
        tau_zenith_sealevel = rayleigh_optical_depth(w_cm)

        # convert to the current altitude
        tau_zenith = tau_zenith_sealevel*np.exp(-self.atmosphere.altitude)

        # figure out the airmass along this line of sight
        # airmass = 1/np.cos(self.zenith_angle)
        effective_airmass = np.minimum(self.atmosphere.fortney_factor()/2.0,
                                       self.airmass)

        # actual optical depth along this line of sight
        tau = tau_zenith*effective_airmass


        '''
        An implicit assumption in this modeling is that the mean intensity
        field is constant everywhere. That's blatantly not true at sunset,
        where different zenith angle (and altitudes integrated across) are
        illuminated with very different spectra. However, hopefully this is
        a reasonable-ish start?

        It's also assuming exactly Earth-like scattering. That's
        definitely not coolsies.
        '''

        return tau


    def spectrum(self, wavelength=None):
        '''
        Calculate intensity of the diffuse sky at a given zenith angle.
        '''

        # calculate the optical depth at this wavelength/angle
        tau_scattering = self.tau_rayleigh_scatter(wavelength=wavelength)

        # (very crudely) assume the atmosphere is all scattering
        albedo = 1.0
        # (this could instead be estimated from the actual transmission)

        # calculate the mean intensity field
        mean_intensity = self.sunset.mean_intensity(wavelength)

        # calculate the intensity of the illuminated sky (away from the disk)
        sky_intensity = albedo*mean_intensity*(1 - np.exp(-tau_scattering))

        # return this sky intensity
        return sky_intensity
