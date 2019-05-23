from ..imports import *
from ..sources import Spectrum
from ..animatetools import *

class Sunset(Spectrum):

    def __init__(self, source, atmosphere):
        '''
        Initialize this composite object, connecting a light source
        and an atmosphere together into a Sunset calculator.
        '''

        self.source = source
        self.atmosphere = atmosphere

    def __repr__(self):
        '''
        How should this be represented as a string?
        '''

        s = repr(self.source)
        t = repr(self.atmosphere)
        return f'<{s} through {t}>'

    def set_zenith_angle(self, *args, **kwargs):
        '''
        Pass zenith commands up to the atmosphere.
        '''
        self.atmosphere.set_zenith_angle(*args, **kwargs)

    def set_altiude(self, *args, **kwargs):
        '''
        Pass altitude commands up to the atmosphere.
        '''
        self.atmosphere.set_altitude(*args, **kwargs)

    def spectrum(self, wavelength=None, zenith_angle=None, altitude=None):
        '''
        The spectrum of the light source viewed through the atmosphere.

        Parameters
        ----------
        wavelength : astropy.units.quantity.Quantity
            The wavelengths on which we want the spectrum.

        Returns
        -------
        spectrum : astropy.units.quantity.Quantity
            The luminosity (W/nm) or flux (W/nm/m**2).
        '''

        # the transmission wavelength sets the grid
        w = self.atmosphere.wavelength(wavelength)

        # get the transmission a
        t = self.atmosphere.transmission(wavelength=w,
                                         zenith_angle=zenith_angle,
                                         altitude=altitude)

        # get the original spectrum
        s = self.source.spectrum(wavelength=w)

        # return the transmitted spectrum
        return s*t


    def plot_sunset(self, overhead=45, set=95):

        zenith_angles = np.arange(overhead, set, 1)*u.deg
        elevation_angles = 90*u.deg - zenith_angles
        #azimuth_angles = np.zeros_like(zenith_angles)

        colors = []

        with plt.style.context('dark_background'), quantity_support():
            ax = plt.subplot()
            plt.axis('scaled')
            plt.xlim(-5*u.deg, 5*u.deg)
            plt.ylim(0*u.deg, overhead)
            for z in tqdm(zenith_angles):
                self.set_zenith_angle(z)
                #colors.append(self.to_color())
                self.disk(z, 0.0)


    def animate_sunset(self, maxelevation=30, minelevation=-5, filename='sunset.mp4'):
        wri = get_writer(filename)
        zenith_angles = np.arange(90-maxelevation, 90-minelevation, 0.5)*u.deg
        #azimuth_angles = np.zeros_like(zenith_angles)

        colors = []

        with plt.style.context('dark_background'), quantity_support():
            fi, ax = plt.subplots(1,1)

            with wri.saving(fi, filename, 400):

                for z in tqdm(zenith_angles):
                    plt.cla()
                    self.set_zenith_angle(z)
                    #colors.append(self.to_color())
                    self.disk(z, 0.0)
                    plt.axis('scaled')
                    plt.xlim(-5*u.deg, 5*u.deg)
                    plt.ylim(0*u.deg, maxelevation*u.deg + self.source.angular_size()*2)

                    wri.grab_frame()

    def disk(self, zenith_angle, azimuth_angle=0*u.deg):
        '''
        Create a simple cartoon of the star,
        with appropriate color, angular size, and position.
        '''

        # set the zenith angle
        self.set_zenith_angle(zenith_angle)

        # figure out the color of the star
        rgb = self.to_color()

        c = plt.Circle(xy=[azimuth_angle, 90*u.deg - zenith_angle],
                       radius=self.source.angular_size(),
                       facecolor=rgb)
        plt.gca().add_patch(c)
        return c
