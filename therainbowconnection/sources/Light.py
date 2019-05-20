import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import copy

from scipy.integrate import quad

import astropy.units as u
from astropy import constants

from colour.plotting import plot_visible_spectrum, plot_single_sd
from colour import SpectralDistribution

from ..colortools import plot_rainbow, CMF

class Light:
    # the default grid of wavelengths
    default_wavelengths = np.arange(200, 1000)*u.nm

    def __repr__(self):
        '''
        How should this object appear as a string?
        '''
        return f'<{self.__class__.__name__}>'

    def at(self, distance=1*u.au):
        '''
        Return a new Spectrum in which we're viewing
        the current source from some distance.

        Parameters
        ----------
        distance : astropy.units.quantity.Quantity
            The distance at which we're viewing this source.
        '''

        # create a copy of the current spectrum
        new = copy.deepcopy(self)

        # update this copy's distance and return
        new.distance = distance
        return new

    def surface_flux(self, wavelength):
        '''
        The surface flux of the light source,
        in units of W/nm/m**2.
        '''

        # make a boring spectrum
        shape = np.shape(wavelength)
        return np.ones(shape)*u.W/u.nm/u.m**2

    def surface_area(self):
        '''
        The surface area of the light source,
        in units of m**2.
        '''
        return 1*u.m**2

    def normalization(self):
        '''

        '''
        try:
            assert(self.distance is not None)
            return 4*np.pi*self.distance**2
        except (AttributeError, AssertionError):
            return 1.0

    def spectrum(self, wavelength):
        '''
        The spectrum of the light source,
        as specific luminosity (W/nm)
        or specific flux (W/nm/m**2).

        Parameters
        ----------
        wavelength : numpy.ndarray
            The wavelengths at which we want the spectrum.
        '''

        # simplify the factor as best you can
        factor = (self.surface_area()/self.normalization()).decompose()
        return factor*self.surface_flux(wavelength)


    def to_sd(self):
        w = CMF.wavelengths
        f = self.spectrum(w).value
        sd = SpectralDistribution(dict(zip(w.value, f.value)))
        return sd

    def plot(self, ax=None, wavelength=None, rainbow=True,
                   color='white', **kwargs):
        '''
        A quick tool to plot a spectrum,
        in units of W/nm (specific luminosity).

        Parameters
        ----------
        ax : matplotlib.axes._subplots.AxesSubplot
            Specify an axes object into which
            this plot should be drawn. For example,
            ax=plt.gca() to overplot directly on
            the current plot.

        wavelength : numpy.ndarray
            A grid of wavelengths on which
            the spectrum should be plotted.

        rainbow : bool
            Should we add an extra rainbow above
            the plot, to indicate visible light?

        color : str
            The color for drawing the spectrum.
            'auto' should represent

        '''

        # set up to use a dark background for the plot
        with plt.style.context('dark_background'):

            # create new ax(s), unless we're supposed to over plot on one
            if ax is None:
                # decide whether to add a horizontal rainbow cartoon
                if rainbow:
                    # create a two-part grid
                    gs = GridSpec(  2, 1,
                                    height_ratios=[0.1, 1],
                                    hspace=0.0)


                    # plot a cartoon rainbow, in a box above
                    ax_rainbow = plt.subplot(gs[0])
                    plt.sca(ax_rainbow)
                    plot_rainbow(axes=ax_rainbow)
                    plt.axis('off')

                    # create the main ax
                    ax = plt.subplot(gs[1], sharex=ax_rainbow)


                else:
                    # simply create one simple ax
                    ax = plt.subplots(1, 1)


            # make sure we point back at the first plot
            plt.sca(ax)

            # make sure at least some wavelengths are defined
            w = (wavelength or self.default_wavelengths).to('nm')

            # pull out the spectrum
            f = self.spectrum(w)

            # plot the spectrum
            plt.plot(w, f, color=color, label=self, **kwargs)

            # add the axis labels
            wunit = w.unit.to_string('latex')
            funit = f.unit.to_string('latex')
            plt.xlabel(f'Wavelength ({wunit})')
            plt.ylabel(f'Flux ({funit})')

            plt.xlim(np.min(w).value, np.max(w).value)

            # set the ylimits
            #plt.ylim(1e-20*max(self.flux), 10*max(self.flux))

            #plt.xscale('log')
            #plt.yscale('log')





        return ax

    def integrated_surface_flux(self, lower=None, upper=None):
        '''
        Integrate the surface flux spectrum,
        returning units of W/m**2.

        Parameters
        ----------

        lower : astropy.units.quantity.Quantity
            The lower wavelength limit.

        upper : astropy.units.quantity.Quantity
            The lower wavelength limit.
        '''

        if lower is not None:
            raise NotImplementedError('Wavelength limits not yet OK.')

        w = self.default_wavelengths
        f = self.surface_flux(w)

        return np.trapz(f, w)

        #np.trapz(f.value, w.value)*f.unit*w.unit
        #wlower = lower or self.default_wavelengths[0]
        #wupper = upper or self.default_wavelengths[-1]
        #return quad(self.spectrum, wlower, wupper)

    '''
    def normalize(self, power=100*u.W):

        #Renormalize the
        total = scipy.integrate.trapz(self.flux, self.wavelength)
        self.power = power
        self.flux = self.flux/total*power
    '''
