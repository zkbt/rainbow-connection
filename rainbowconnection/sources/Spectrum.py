from ..imports import *
from ..colortools import plot_rainbow, rainbow_spectrum, CMFs, SpectralDistribution
import colour
from ..units import determine_quantity
from ..plottingtools import setup_axes_with_rainbow


bgkw = dict(color='gray', alpha=0.5, zorder=-100)

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

        Returns
        -------
        surface_area : astropy.units.quantity.Quantity
            The emitting area of the surface, usually in m**2.
        '''
        return 4*np.pi*self.radius**2


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

    def angular_size(self):
        '''
        The angular size of the light source (if viewed from a distance).
        '''

        try:
            # if there's a distance, return an angular size
            assert(self.distance is not None)
            return np.arctan(self.radius/self.distance).to('deg')
        except (AttributeError, AssertionError):
            # complain if no distance is defined
            raise ValueError('''
            This Spectrum has no .distance attribute.
            Please consider using `.at(distance)` to
            create a new light source as viewed from
            a distance.
            ''')

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
    def integrate(self, lower=None, upper=None):
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

        w = self.default_wavelengths
        f = self.spectrum(w)


        ok = np.ones(np.shape(w)).astype(np.bool)
        if lower is not None:
            ok *= w >= lower

        if upper is not None:
            ok *= w <= upper
            #raise NotImplementedError('Wavelength limits not yet OK.')


        return np.trapz(f[ok], w[ok])

        #np.trapz(f.value, w.value)*f.unit*w.unit
        #wlower = lower or self.default_wavelengths[0]
        #wupper = upper or self.default_wavelengths[-1]
        #return quad(self.spectrum, wlower, wupper)



    def to_sd(self):
        '''
        Create a `colour` SpectralDistribution from this object,
        solely covering the visible range. This can be used
        to estimate the true color of this spectrum.
        '''

        # pick wavelengths directly from the color-matching functions
        w = CMFs.wavelengths*u.nm

        # normalize only within the visible range
        ok = (w < 700*u.nm) & (w > 390*u.nm)
        w = w[ok]
        # (note: this is a kludge to avoid making spectra
        #  with most of their luminosity outside the visible
        #  appear as really dark and dull. there is probably
        #  a much cleverer way of making sure the normalization
        #  does something reasonable.)

        # calculate the spectrum at those wavelengths
        f = self.spectrum(w)

        # create the spectral distribution
        sd = SpectralDistribution(dict(zip(w.value, f.value/np.max(f.value))))

        # return it
        return sd

    def to_color(self):
        '''
        Determine the RGB color of this spectrum.

        Returns
        -------
        rgb : numpy.ndarray
            3-element array containing RGB values
            that can be fed into matplotlib.
        '''

        # create a colour SpectralDistribution
        sd = self.to_sd()

        # convert to XYZ
        # (maybe use `k=` option for relative scaling between sources?)
        XYZ = colour.sd_to_XYZ(sd)
        if (np.min(XYZ) < 0) or np.max(XYZ) > 100:
            print(f'XYZ={XYZ} is outside of [0, 100]!')

        # convert to RGB
        RGB = colour.XYZ_to_sRGB(XYZ / 100)
        if (np.min(RGB) < 0) or np.max(RGB) > 1:
            pass
            #print(f'RGB={RGB} is outside of [0, 1]!')

        # trim out underenderable colors (this is sneaky!)
        #clipped_RGB = np.maximum(0, np.minimum(1, RGB))
        clipped_RGB = np.maximum(0, RGB)
        # instead of clip, should we add to everything until we get to zero?

        # kludge to maximize brightness, for every color!
        clipped_RGB /= np.max(clipped_RGB)
        return clipped_RGB

    def __repr__(self):
        '''
        How should this object appear as a string?

        Returns
        -------
        s : str
            A simple string representation.
        '''
        try:
            assert(self.distance is not None)
            return f'{self.__class__.__name__} at {self.distance}'
        except (AssertionError, AttributeError):
            return f'{self.__class__.__name__}'

    def set_power(self, power=100*u.W):
        '''
        Change the radius of the object so that it will
        emit a specified power, with the same spectral shape.

        Parameters
        ----------
        power : astropy.units.quantity.Quantity
            The total power we want the object to emit.
        '''

        # calculate the total luminosity of this object
        total = self.integrate()

        # make sure we're dealing with an actual luminosity
        assert(total.unit.is_equivalent('W'))

        # calculation a new normalization
        normalization = (power/total).decompose()

        # change the radius of this object
        self.radius *= np.sqrt(normalization)
        self.power = power

    def mean_intensity(self, wavelength=None):
        '''
        Calculate the mean intensity field created by the source.

        The mean intensity field represents the intensity
        of the source, smeared over a full 4pi steradians.

        (This makes sense only for spectra viewed from a distance.)
        '''

        # what is the intensity of the disk
        F_disk = self.spectrum(wavelength)


        # make sure we're dealing with a flux
        assert(F_disk.unit.is_equivalent(u.W/u.m**2/u.nm))

        # calculate the mean intensity
        solid_angle = np.pi*self.angular_size()**2

        # calculate the mean intensity field
        J = F_disk/4/np.pi/u.sr

        return J

    def disk_intensity(self, wavelength=None):
        '''
        Calculate the intensity of the disk of the source.

        The disk intensity represents the intensity of staring
        directly at the disk. The flux of the source is its
        intensity integrated over solid angle, so two sources
        with the same intensities can have very different
        fluxes, based on their apparent angular sizes.
        For example, the Sun seen `.at` different distances
        will have different fluxes but the same disk intensity.

        (This makes sense only for spectra viewed from a distance.)
        '''

        # what is the intensity of the disk
        F_disk = self.spectrum(wavelength)


        # make sure we're dealing with a flux
        assert(F_disk.unit.is_equivalent(u.W/u.m**2/u.nm))

        # calculate the mean intensity
        solid_angle = np.pi*self.angular_size()**2
        I_disk = F_disk/solid_angle

        return I_disk

    # ===================================
    #
    #        PLOTTING METHODS
    #
    # ===================================

    def plot(self,  ax=None,
                    wavelength=None,
                    rainbow=True,
                    color='auto',
                    style='dark_background',
                    figsize=(5, 2.5),
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

        # set up to use a dark background for the plot; make sure units match
        with plt.style.context(style), quantity_support():

            # setup the basic axes
            ax = setup_axes_with_rainbow(ax=ax, rainbow=rainbow, figsize=figsize)

            # make sure at least some wavelengths are defined
            w = self.wavelength(wavelength)

            # pull out the spectrum
            f = self.spectrum(w)

            # plot the spectrum
            if color == 'auto':
                color = self.to_color()
                background = ax.get_facecolor()[0:3]
                if np.max(color - background) < 0.05:
                    print(f'''
                    The inferred color {color} might be a little
                    too close to {background} to be visible. Consider
                    plotting without the `color='auto'` option.
                    ''')
            plt.plot(w, f, color=color, label=self, **kwargs)

            # add the axis labels
            wunit = w.unit.to_string('latex_inline')
            funit = f.unit.to_string('latex_inline')
            plt.xlabel(f'Wavelength ({wunit})')
            plt.ylabel(f'{determine_quantity(f.unit)} ({funit})')

        return ax


    def plot_rgb(self,  ax=None,
                        wavelength=None,
                        rainbow=True,
                        foreground=True,
                        style='dark_background',
                        figsize=(5, 2.5),
                        **kwargs):
        '''
        A quick tool to plot a spectrum, just as RGB bars.

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
            ax = setup_axes_with_rainbow(ax=ax, rainbow=rainbow, figsize=figsize)

            # make sure at least some wavelengths are defined
            w = self.wavelength(wavelength)

            rgb = self.to_color()

            blue = [400, 495]
            green = [495, 590]
            red = [590, 685]

            centers = [np.mean(c) for c in [red, green, blue]]
            widths = [c[1] - c[0] for c in [red, green, blue]]

            if foreground:
                kw = dict(color = [np.array([1, 0, 0]),
                                   np.array([0, 1, 0]),
                                   np.array([0, 0, 1])],
                          edgecolor='white',
                          zorder = 0)
            else:

                kw = dict(**bgkw)
                kw['edgecolor'] = 'none'
            plt.bar(centers, rgb*100, widths, **kw) # linewidth=2,


            # add the axis labels
            wunit = w.unit.to_string('latex_inline')
            plt.xlabel(f'Wavelength ({wunit})')
            plt.ylabel(f'Brightness (%)')

        return ax

    def plot_as_rainbow(self,  ax=None,
                                rainbow=True,
                                color='auto',
                                style='dark_background',
                                foreground=True,
                                ylim=[0, 120],
                                figsize=(5,2.5),
                                **kwargs):
        '''
        A quick tool to plot a spectrum, just as RGB bars.

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
            ax = setup_axes_with_rainbow(ax=ax, rainbow=rainbow, figsize=figsize)

            # make sure at least some wavelengths are defined
            w = self.wavelength(np.arange(330, 760, 1)*u.nm)
            f = self.spectrum(w)
            # KLUDGE?
            norm = np.max(f.value[(w > 400*u.nm) & (w < 685*u.nm)]) / 100
            if foreground:
                rainbow_spectrum(axes=ax, wavelength=w.to('nm').value, flux=f.value/norm,
                                 rainbowtop=np.max(ylim))
                plt.plot(w, f/norm, color='white') #, linewidth=2
            else:
                plt.fill_between(w, f/norm, linewidth=0, **bgkw)
            plt.ylim(*ylim)

            # add the axis labels
            wunit = w.unit.to_string('latex_inline')
            plt.xlabel(f'Wavelength ({wunit})')
            plt.ylabel(f'Brightness (%)')

        return ax

    def cartoon_sun(self):
        '''
        Create a simple cartoon of the star,
        with appropriate color and angular size.
        '''

        # figure out the color of the star
        rgb = self.to_color()

        with plt.style.context('dark_background'):
            fi = plt.figure(figsize=(5, 2.5))
            ax = fi.add_axes([0, 0, 1, 1])
            plt.scatter(0, 0, c=rgb, s=8000)
            plt.xticks([])
            plt.yticks([])
            plt.axis('off')
