from ..imports import *
from ..sources import Spectrum
from ..animatetools import *
from .sky import *

class Sunset(Spectrum):

    def __init__(self, source, atmosphere):
        '''
        Initialize this composite object, connecting a light source
        and an atmosphere together into a Sunset calculator.
        '''

        self.source = source
        self.atmosphere = atmosphere
        self.sky = Sky(self)

    # KLUDGE - there's got to be a better way to organize this?!
    @property
    def radius(self):
        return self.source.radius

    @property
    def distance(self):
        try:
            return self.source.distance
        except AttributeError:
            return None

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




    def plot_disk(self, ax=None, zenith_angle=88*u.deg, azimuth_angle=0*u.deg):
        '''
        Plot the disk of the star.

        Parameters
        ----------
        zenith_angle : astropy.units.quantity.Quantity
            The angle away from zenith.
        azimuth_angle : astropy.units.quantity.Quantity
            The azimuth angle along the horizon.
        '''

        if ax is not None:
            plt.sca(ax)

        # set the zenith angle
        self.set_zenith_angle(zenith_angle)

        # figure out the color of the star
        rgb = self.to_color()

        c = plt.Circle(xy=[azimuth_angle, 90*u.deg - zenith_angle],
                       radius=self.source.angular_size(),
                       facecolor=rgb, zorder=1000)
        plt.gca().add_patch(c)
        return c

    def plot_sky(self, maxelevation=20*u.deg,
                       minelevation=0*u.deg,
                       skyresolution=0.5*u.deg,
                       skynormalization=0.7):
        '''
        Plot a diffuse sky behind the disk of the star.

        Parameters
        ----------
        maxelevation : astropy.units.quantity.Quantity
            The maximum elevation above the horizon.
        minelevation : astropy.units.quantity.Quantity
            The minimum elevation above the horizon.
        skyresolution : astropy.units.quantity.Quantity
            The vertical extent of constant-color stripes
            when rendering the atmosphere. Stripes larger
            than ~1/20 the vertical extent will be very
            noticeable.
        skynormalization : float
            By what factor do we down-weight the intensity
            of the sky. To be quantitatively correct, this
            should be something like the ratio of the
            stellar disk intensity to the diffuse sky
            brightness, but that factor is so tiny that
            the sky would always seem to be black.
        '''

        # set the elevations where we will calculate sky colors
        elevations = np.arange(minelevation.to('deg').value,
                               (maxelevation + skyresolution).to('deg').value,
                               skyresolution.to('deg').value)*u.deg

        # calculate a list of sky colors, one for each elevation
        colors = []
        for e in elevations:

            # point at a particular stripe of the sky
            self.sky.set_zenith_angle(90*u.deg - e)

            # estimate the color of that stripe
            colors.append(self.sky.to_color()*skynormalization)

        # draw horizontal bars with these colors
        plt.barh(y=elevations,
                 width=20,
                 height=skyresolution,
                 left=-10,
                 zorder=-100,
                 color=colors)

    def plot_sunset(self, ax=None,
                          zenith_angle=88*u.deg,
                          azimuth_angle=0*u.deg,
                          maxelevation=20*u.deg,
                          minelevation=-5*u.deg,
                          skyresolution=0.5*u.deg,
                          skynormalization=0.7):

        '''
        Plot the disk of the star and the diffuse sky.

        Parameters
        ----------
        ax : matplotlib.axes._subplots.AxesSubplot
            The axes into which this plot should be drawn.
        zenith_angle : astropy.units.quantity.Quantity
            The angle away from zenith.
        azimuth_angle : astropy.units.quantity.Quantity
            The azimuth angle along the horizon.
        maxelevation : astropy.units.quantity.Quantity
            The maximum elevation above the horizon.
        minelevation : astropy.units.quantity.Quantity
            The minimum elevation above the horizon.
        skyresolution : astropy.units.quantity.Quantity
            The vertical extent of constant-color stripes
            when rendering the atmosphere. Stripes larger
            than ~1/20 the vertical extent will be very
            noticeable.
        skynormalization : float
            By what factor do we down-weight the intensity
            of the sky. To be quantitatively correct, this
            should be something like the ratio of the
            stellar disk intensity to the diffuse sky
            brightness, but that factor is so tiny that
            the sky would always seem to be black.
        '''

        # make sure we're using a dark background and astropy units
        with plt.style.context('dark_background'), quantity_support():

            if ax is not None:
                plt.sca(ax)

            # put the sun at a particular elevation
            self.set_zenith_angle(zenith_angle)

            # plot the diffuse sky
            self.plot_sky(maxelevation=maxelevation,
                          minelevation=minelevation)

            # plot the disk of the sun
            self.plot_disk(zenith_angle=zenith_angle,
                           azimuth_angle=azimuth_angle)

            # fuss with the axis labels
            plt.axis('scaled')
            plt.xlim(-5*u.deg, 5*u.deg)
            plt.ylim(0*u.deg, maxelevation)
            plt.ylabel('Angle above Horizon')
            ax = plt.gca()
            ax.get_xaxis().set_visible(False)
        return ax

    def animate_sunset(self, filename='sunset.mp4',
                             maxelevation=30*u.deg,
                             skyresolution=0.5*u.deg,
                             skynormalization=0.7,
                             motionresolution=0.5*u.deg):
        '''
        Animate the sun setting.

        Parameters
        ----------
        filename : string
            The filename where this animation should be saved
        maxelevation : astropy.units.quantity.Quantity
            The maximum elevation above the horizon.
        skyresolution : astropy.units.quantity.Quantity
            The vertical extent of constant-color stripes
            when rendering the atmosphere. Stripes larger
            than ~1/20 the vertical extent will be very
            noticeable.
        skynormalization : float
            By what factor do we down-weight the intensity
            of the sky. To be quantitatively correct, this
            should be something like the ratio of the
            stellar disk intensity to the diffuse sky
            brightness, but that factor is so tiny that
            the sky would always seem to be black.
        motionresolution : astropy.units.quantity.Quantity
            The angular step the star moves between
            frames of the animation.
        '''


        # get the appropriate animation writer
        wri = get_writer(filename)

        # calculate the grid of zenith angles for the star
        min_zenith = 90*u.deg - maxelevation - self.source.angular_size()
        minelevation = -(self.source.angular_size() + motionresolution)
        max_zenith = 90*u.deg - minelevation
        zenith_angles = np.arange(min_zenith.to('deg').value,
                                  max_zenith.to('deg').value,
                                  motionresolution.to('deg').value)*u.deg


        # make a figure with a dark background
        with plt.style.context('dark_background'), quantity_support():

            # create a new figure
            fi, ax = plt.subplots(1,1)

            # save to frames of an animation
            with wri.saving(fi, filename, 400):

                # loop over zenith angles
                for z in tqdm(zenith_angles):

                    # clear whatever was in the previous axes
                    ax.cla()

                    # plot the sunset at this current stellar zenith angle
                    self.plot_sunset(ax=ax,
                                     zenith_angle=z,
                                     maxelevation=maxelevation,
                                     minelevation=minelevation)

                    # grab this frame in the animation
                    wri.grab_frame()

    def animate_everything(self, filename='everything-sunset.mp4',
                                 maxelevation=90*u.deg,
                                 skyresolution=0.5*u.deg,
                                 skynormalization=0.7,
                                 motionresolution=0.5*u.deg):

        # get the appropriate animation writer
        wri = get_writer(filename)

        # calculate the grid of zenith angles for the star
        min_zenith = 90*u.deg - maxelevation - self.source.angular_size()
        minelevation = -(self.source.angular_size() + motionresolution)
        max_zenith = 90*u.deg - minelevation
        zenith_angles = np.arange(min_zenith.to('deg').value,
                                  max_zenith.to('deg').value,
                                  motionresolution.to('deg').value)*u.deg

        with plt.style.context('dark_background'), quantity_support():

            fi = self.plot_everything(zenith_angles[0])

            # save to frames of an animation
            with wri.saving(fi, filename, 160):

                # loop over zenith angles
                for z in tqdm(zenith_angles):

                    # clear whatever was in the previous axes
                    fi.clf()

                    # plot the sunset at this current stellar zenith angle
                    self.plot_everything(fi=fi, zenith_angle=z)

                    # grab this frame in the animation
                    wri.grab_frame()


    def plot_everything(self, zenith_angle=85*u.deg, fi=None):

        with plt.style.context('dark_background'), quantity_support():

            if fi is None:
                xpixels, ypixels = 1920, 1080
                fi = plt.figure(figsize=(12,6.75))
                dpi = 1920/12 # 160

            # set up the basic columns
            gs = GridSpec(  1, 3,
                            width_ratios=[0.4, 0.2, 1],
                            wspace=0.35)

            ax = {}

            # set up the contextualizing plots on the left
            gs_geometry = GridSpecFromSubplotSpec(2, 1,
                                                  height_ratios=[1, 2],
                                                  subplot_spec=gs[0])

            #ax['cartoon'] = plt.subplot(gs_geometry[0])
            ax['sky-zoom'] = plt.subplot(gs_geometry[0])

            # set up the main panel showing the sunset over the full sky
            ax['sky'] = plt.subplot(gs[1])

            # set up the spectrum panels
            gs_spectra = GridSpecFromSubplotSpec(2, 1,
                                                 height_ratios=[1, 1],
                                                 subplot_spec=gs[2],
                                                 hspace=0.5)
            ax['rgb'] = plt.subplot(gs_spectra[0])
            ax['spectrum'] = plt.subplot(gs_spectra[1],
                                         sharex=ax['rgb'], sharey=ax['rgb'])

            # actually make the plots
            self.plot_sunset(ax=ax['sky'], zenith_angle=zenith_angle, maxelevation=90*u.deg)

            self.plot_disk(ax=ax['sky-zoom'], zenith_angle=zenith_angle)
            plt.axis('scaled')
            plt.axis('off')

            #self.plot_sunset(ax=ax['sky-zoom'], zenith_angle=zenith_angle, maxelevation=90*u.deg)
            width = 1.1*self.source.angular_size()
            plt.xlim(-width, width)
            plt.ylim(90*u.deg - zenith_angle - width, 90*u.deg - zenith_angle + width)
            ax['sky-zoom'].get_yaxis().set_visible(False)

            self.plot_rgb(ax=ax['rgb'])
            self.plot_as_rainbow(ax['spectrum'])
            plt.ylim(0, 120)
            plt.xlim(350*u.nm, 750*u.nm)

        return fi
