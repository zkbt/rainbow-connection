from ..imports import *
from ..units import determine_quantity
from ..plottingtools import setup_axes_with_rainbow
from .transmitted import TransmittedSpectrum

class Atmosphere:
    def __repr__(self):
        '''
        How should this object appear as a string?

        Returns
        -------
        s : str
            A simple string representation.
        '''
        return f'{self.__class__.__name__}'

    # FIXME (maybe both spectrum and atmosphere should inherit from the same thing?)
    def wavelength(self, wavelength=None):
        '''
        A wrapper to ensure at least some grid of wavelengths
        gets defined. A default grid will be assumed, unless
        any wavelength array at all is passed.
        '''

        # make sure at least some wavelengths are defined
        if wavelength is None:
            wavelength = self.default_wavelengths
        return wavelength.to('nm')

    def plot(self,  ax=None,
                    wavelength=None,
                    rainbow=True,
                    color=None,
                    style='dark_background',
                    **kwargs):
        '''
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

        '''

        # set up to use a dark background for the plot; make sure units match
        with plt.style.context(style), quantity_support():

            # setup the basic axes
            ax = setup_axes_with_rainbow(ax=ax, rainbow=rainbow)

            # make sure at least some wavelengths are defined
            w = self.wavelength(wavelength)

            # pull out the spectrum
            t = self.transmission(w)

            plt.plot(w, t*100, color=color, label=self, **kwargs)

            # add the axis labels
            wunit = w.unit.to_string('latex_inline')
            plt.xlabel(f'Wavelength ({wunit})')
            plt.ylabel(f'Transmission (%)')

            plt.xlim(np.min(w).value, np.max(w).value)

        return ax

    def set_elevation_angle(self, elevation=90*u.deg):
        '''
        Set the elevation angle above the horizon.

        Parameters
        ----------
        elevation : astropy.units.quantity.Quantity
            The angle above the horizon along which
            the transmission of the atmosphere should
            be calculated.
        '''
        self.set_zenith_angle(90*u.deg - elevation)

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
        self.zenith_angle = z
        self.airmass = 1/np.cos(z)

    def transmit(self, spectrum, wavelength=None, zenith_angle=None):
        '''
        '''

        # create
        return TransmittedSpectrum(source=spectrum, atmosphere=self)
