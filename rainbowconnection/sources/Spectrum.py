from ..imports import *
from ..colortools import plot_rainbow, CMF, SpectralDistribution
from ..units import determine_quantity

class Spectrum:
    '''
    The Spectrum class is a generic representation of the light from
    some emitting object, particularly for spherically symmetric emission.

    By default, the `.spectrum(wavelength)` method will return the
    spectral luminosity of the object. This is a quantity with units
    like W/nm, and it can be integrated over wavelength to provide
    the total luminosity of the light-emitting object, in W.

    The `.at(distance)` method creates an object representing the spectral
    flux from the object if seen from some particular distance. This is
    a quantity with units like W/nm/m**2, and it can be integrated over
    wavelength to provide the total flux of the light, in W/m**2.

    Classes that inherit from this will likely modify (at least) the
    surface_flux and surface_area methods.
    '''

    # the default grid of wavelengths
    default_wavelengths = np.arange(200, 1000)*u.nm

    def surface_flux(self, wavelength=None):
        '''
        The surface flux of the light source,
        in units like W/nm/m**2.

        (This should likely by overwritten by inheriting classes.)

        Parameters
        ----------
        wavelength : astropy.units.quantity.Quantity
            The wavelengths on which we want the flux.

        Returns
        -------
        surface_flux : astropy.units.quantity.Quantity
            The flux from the surface, usually in W/nm/m**2.
        '''

        # make sure a wavelength grid is defined
        w = self.wavelength(wavelength)

        # make a boring spectrum
        shape = np.shape(w)
        return np.ones(shape)*u.W/u.nm/u.m**2

    def surface_area(self):
        '''
        The surface area of the light source,
        in units like m**2.

        (This should likely by overwritten by inheriting classes.)

        Returns
        -------
        surface_area : astropy.units.quantity.Quantity
            The emitting area of the surface, usually in m**2.
        '''

        # make a very simple surface area
        return 1*u.m**2

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

    def spectrum(self, wavelength=None):
        '''
        The spectrum of the light source, as spectral luminosity (W/nm)
        or if a distance is defined as spectral flux (W/nm/m**2).

        Parameters
        ----------
        wavelength : astropy.units.quantity.Quantity
            The wavelengths on which we want the spectrum.

        Returns
        -------
        spectrum : astropy.units.quantity.Quantity
            The luminosity (W/nm) or flux (W/nm/m**2).
        '''

        # simplify the factor as best we can
        factor = (self.surface_area()/self.normalization()).decompose()

        # return the surface flux with appropriate normalization
        return factor*self.surface_flux(wavelength)

    def normalization(self):
        '''
        The normalization by which this Spectrum
        will be divided. It's flexible, to allow
        either spectral luminosity or spectral flux.

        Returns
        -------
        norm : astropy.units.quantity.Quantity
            A normalization
        '''
        try:
            # if there's a distance, return a flux
            assert(self.distance is not None)
            return 4*np.pi*self.distance**2
        except (AttributeError, AssertionError):
            # by default, simply return a luminosity
            return 1.0

    def at(self, distance=1*u.au):
        '''
        Create a new Spectrum representing the current
        light source viewed from some distance
        (assuming spherical symmetry).

        Parameters
        ----------
        distance : astropy.units.quantity.Quantity
            The distance at which we're viewing this source.

        Returns
        -------
        flux : Spectrum
            A new Spectrum, with the distance attached.
        '''

        # create a copy of the current spectrum
        new = copy.deepcopy(self)

        # update this copy's distance and return
        new.distance = distance
        return new

    # FIXME -- for analytic functions, it'd help to define some kind of
    # a bounding box in wavelength space, so this integral could be done
    # analytically or with scipy.integrate.quad
    def integrated_surface_flux(self, lower=None, upper=None):
        '''
        Integrate the spectrum over wavelength.

        It gives a number with units identical to the results of
        `.spectrum()` but without the wavelength (W or W/m**2).

        Parameters
        ----------
        lower : astropy.units.quantity.Quantity
            The lower wavelength limit.

        upper : astropy.units.quantity.Quantity
            The lower wavelength limit.

        Returns
        -------
        integral : astropy.units.quantity.Quantity
            The integral over wavelength
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


    def plot(self,  ax=None,
                    wavelength=None,
                    rainbow=True,
                    color='white',
                    style='dark_background',
                    **kwargs):
        '''
        A quick tool to plot a spectrum.

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

        # set up to use a dark background for the plot
        with plt.style.context(style):

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

            plt.ylabel(f'{determine_quantity(f.unit)} ({funit})')

            plt.xlim(np.min(w).value, np.max(w).value)

            # set the ylimits
            #plt.ylim(1e-20*max(self.flux), 10*max(self.flux))

            #plt.xscale('log')
            #plt.yscale('log')





        return ax

    def to_sd(self):
        w = CMF.wavelengths
        f = self.spectrum(w).value
        sd = SpectralDistribution(dict(zip(w.value, f.value)))
        return sd



    def __repr__(self):
        '''
        How should this object appear as a string?

        Returns
        -------

        '''
        try:
            assert(self.distance is not None)
            return f'<{self.__class__.__name__} at {self.distance}>'
        except (AssertionError, AttributeError):
            return f'<{self.__class__.__name__}>'

    '''
    def normalize(self, power=100*u.W):

        #Renormalize the
        total = scipy.integrate.trapz(self.flux, self.wavelength)
        self.power = power
        self.flux = self.flux/total*power
    '''
