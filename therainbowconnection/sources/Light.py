import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import astropy.units as u
from astropy import constants

from colour.plotting import plot_visible_spectrum, plot_single_sd
from colour import SpectralDistribution

# the default grid of wavelengths
default_wavelengths = np.arange(200, 1000)*u.nm

class Light:

    def __repr__(self):
        '''
        How should this object appear as a string?
        '''
        return f'<{self.__class__.__name__}>'


    def flux(self, wavelength):
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

    def spectrum(self, wavelength):
        '''
        The specific luminosity of the light source,
        in units of W/nm.
        '''
        return self.surface_area()*self.flux(wavelength)

    def plot(self, ax=None, wavelength=None, rainbow=True):
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
        '''

        # set up to use a dark background for the plot
        with plt.style.context('dark_background'):

            # create an ax, unless we're supposed to over plot on one
            if ax is None:
                if rainbow:
                    _, both_ax = plt.subplots(2, 1,
                                sharex=True,
                                gridspec_kw=dict(height_ratios=[0.2, 1]))

                    ax_rainbow, ax = both_ax

                    # plot a cartoon rainbow
                    plt.sca(ax_rainbow)
                    #plot_visible_spectrum(axes=ax_rainbow)
                    #plt.axis('off')
                    # make sure we point back at
                    plt.sca(ax)

                else:
                    ax = plt.subplots(1, 1)




            # make sure at least some wavelengths are defined
            w = (wavelength or default_wavelengths).to('nm')

            # pull out the spectrum
            f = self.spectrum(w).to('W/nm')

            # plot the spectrum
            plt.plot(w, f, color='white', label=self)

            if rainbow:
                sd = SpectralDistribution(dict(zip(w.value, f.value)))
                plot_single_sd(sd)
            else:
                plt.plot(w, f)

            # add the axis labels
            wunit = w.unit.to_string('latex')
            funit = w.unit.to_string('latex')
            plt.xlabel(f'Wavelength ({wunit})')
            plt.ylabel(f'Flux ({funit})')

            plt.xlim(np.min(w).value, np.max(w).value)

            # set the ylimits
            #plt.ylim(1e-20*max(self.flux), 10*max(self.flux))

            #plt.xscale('log')
            #plt.yscale('log')

        return ax


    def integrate(self, limits=[]):
        '''
        Integrate the spectrum,
        returning units of W.

        Parameters
        ----------

        lower : astropy.u.quantity.Quantity
            The lower wavelength limit.

        upper : astropy.u.quantity.Quantity
            The lower wavelength limit.
        '''


        wlower = lower or default_wavelengths[0]
        wupper = upper or default_wavelengths[-1]

        return scipy.integrate.quad(self.spectrum,
                                    wlower,
                                    wupper)
    '''
    def normalize(self, power=100*u.W):

        #Renormalize the
        total = scipy.integrate.trapz(self.flux, self.wavelength)
        self.power = power
        self.flux = self.flux/total*power
    '''
