from ..imports import *
from ..sources import Spectrum
from ..animatetools import *
from .sky import *


class Sunset(Spectrum):
    def __init__(self, source, atmosphere):
        """
        Initialize this composite object, connecting a light source
        and an atmosphere together into a Sunset calculator.

        Parameters
        ----------
        source : rainbowconnection.Spectrum
            The light source from above the atmosphere.
        atmosphere : rainbowconnection.Atmosphere
            The atmosphere through which the light source propagates.
        """

        self.source = source
        self.atmosphere = atmosphere
        self.sky = Sky(self)

    # KLUDGE? (to make plots work)
    def get_wavelength(self, *args, **kwargs):
        return self.atmosphere.get_wavelength(*args, **kwargs)

    # KLUDGE - there's got to be a better way to organize this?!
    @property
    def radius(self):
        return self.source.radius

    # KLUDGE - there's got to be a better way to organize this?!
    @property
    def distance(self):
        try:
            return self.source.distance
        except AttributeError:
            return None

    def __repr__(self):
        """
        How should this be represented as a string?
        """

        s = repr(self.source)
        t = repr(self.atmosphere)
        return f"{s}\n through\n{t}"

    def set_zenith_angle(self, *args, **kwargs):
        """
        Pass zenith commands up to the atmosphere.
        """
        self.atmosphere.set_zenith_angle(*args, **kwargs)

    def set_altiude(self, *args, **kwargs):
        """
        Pass altitude commands up to the atmosphere.
        """
        self.atmosphere.set_altitude(*args, **kwargs)

    def spectrum(self, wavelength=None, zenith_angle=None, altitude=None):
        """
        The spectrum of the light source viewed through the atmosphere.

        Parameters
        ----------
        wavelength : astropy.units.quantity.Quantity
            The wavelengths on which we want the spectrum.

        zenith_angle : astropy.units.quantity.Quantity
            The angle away from zenith along which
            the transmission of the atmosphere should
            be calculated.

        altitude : float
            The altitude at which we should be floating
            in the atmosphere, in units of scale heights.

        Returns
        -------
        spectrum : astropy.units.quantity.Quantity
            The luminosity (W/nm) or flux (W/nm/m**2).
        """

        # the transmission wavelength sets the grid
        w = self.atmosphere.get_wavelength(wavelength)

        # get the transmission a
        t = self.atmosphere.transmission(
            wavelength=w,
            zenith_angle=zenith_angle,
            altitude=altitude,
        )

        # get the original spectrum
        s = self.source.spectrum(wavelength=w)

        # return the transmitted spectrum
        return s * t

    def transit_depth(self, wavelength=None):
        """
        The transit depth of the planet, passing in front of the light source.

        Parameters
        ----------
        wavelength : astropy.units.quantity.Quantity
            The wavelengths on which we want the spectrum.

        Returns
        -------
        transitdepth : np.array
            The fraction of starlight blocked, unitless.
        """

        # the transmission wavelength sets the grid
        w = self.atmosphere.get_wavelength(wavelength)

        # get the transmission at those wavelengths
        rp = self.atmosphere.transit_radius(wavelength=w)

        # get the original spectrum
        rs = self.source.radius

        return (rp / rs).decompose() ** 2

    def plot_disk(
        self,
        ax=None,
        zenith_angle=89 * u.deg,
        azimuth_angle=0 * u.deg,
    ):
        """
        Plot the disk of the star.

        Parameters
        ----------
        ax : matplotlib.axes._subplots.AxesSubplot
            The axes into which this plot should be drawn.
        zenith_angle : astropy.units.quantity.Quantity
            The angle away from zenith.
        azimuth_angle : astropy.units.quantity.Quantity
            The azimuth angle along the horizon.
        """

        with plt.style.context("dark_background"), quantity_support():

            if ax is not None:
                plt.sca(ax)

            # set the zenith angle
            self.set_zenith_angle(zenith_angle)

            # figure out the color of the star
            rgb = self.to_color()

            c = plt.Circle(
                xy=[azimuth_angle, 90 * u.deg - zenith_angle],
                radius=self.source.angular_size(),
                facecolor=rgb,
                zorder=1000,
            )
            plt.gca().add_patch(c)

            size = self.source.angular_size() * 1.1
            plt.xlim(azimuth_angle - size, azimuth_angle + size)
            plt.ylim(
                90 * u.deg - zenith_angle - size,
                90 * u.deg - zenith_angle + size,
            )
            plt.axis("scaled")
            return c

    def plot_sky(
        self,
        maxelevation=20 * u.deg,
        minelevation=0 * u.deg,
        skyresolution=0.5 * u.deg,
        skynormalization=0.7,
    ):
        """
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
        """

        if self.atmosphere.__class__.__name__ == "NoAtmosphere":
            return

        # set the elevations where we will calculate sky colors
        elevations = (
            np.arange(
                minelevation.to("deg").value,
                (maxelevation + skyresolution).to("deg").value,
                skyresolution.to("deg").value,
            )
            * u.deg
        )

        # calculate a list of sky colors, one for each elevation
        colors = []
        for e in elevations:

            # point at a particular stripe of the sky
            self.sky.set_zenith_angle(90 * u.deg - e)

            # estimate the color of that stripe
            colors.append(self.sky.to_color() * skynormalization)

        # draw horizontal bars with these colors
        plt.barh(
            y=elevations,
            width=20,
            height=skyresolution,
            left=-10,
            zorder=-100,
            color=colors,
        )

    def plot_sunset(
        self,
        ax=None,
        zenith_angle=89 * u.deg,
        azimuth_angle=0 * u.deg,
        maxelevation=20 * u.deg,
        minelevation=-5 * u.deg,
        skyresolution=0.5 * u.deg,
        skynormalization=0.7,
    ):

        """
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
        """

        # make sure we're using a dark background and astropy units
        with plt.style.context("dark_background"), quantity_support():

            if ax is not None:
                plt.sca(ax)

            # put the sun at a particular elevation
            self.set_zenith_angle(zenith_angle)

            # plot the diffuse sky
            self.plot_sky(
                maxelevation=maxelevation,
                minelevation=minelevation,
            )

            # plot the disk of the sun
            self.plot_disk(
                zenith_angle=zenith_angle,
                azimuth_angle=azimuth_angle,
            )

            # fuss with the axis labels
            plt.axis("scaled")
            plt.xlim(-5 * u.deg, 5 * u.deg)
            plt.ylim(0 * u.deg, maxelevation)
            plt.ylabel("Angle above Horizon")
            ax = plt.gca()
            ax.get_xaxis().set_visible(False)
        return ax

    def animate_sunset(
        self,
        filename="sunset.mp4",
        maxelevation=30 * u.deg,
        skyresolution=0.5 * u.deg,
        skynormalization=0.7,
        motionresolution=0.5 * u.deg,
        overwrite=False,
    ):
        """
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
        """

        if os.path.exists(filename):
            if overwrite is False:
                message = f"""
                Filename {filename} already exists.
                Not recreating the movie, unless
                you rerun with `overwrite=True`.
                """
                warnings.warn(message)
                return

        # get the appropriate animation writer
        wri = get_writer(filename)

        # calculate the grid of zenith angles for the star
        min_zenith = 90 * u.deg - maxelevation - self.source.angular_size()
        minelevation = -(self.source.angular_size() + motionresolution)
        max_zenith = 90 * u.deg - minelevation
        zenith_angles = (
            np.arange(
                min_zenith.to("deg").value,
                max_zenith.to("deg").value,
                motionresolution.to("deg").value,
            )
            * u.deg
        )

        # make a figure with a dark background
        with plt.style.context("dark_background"), quantity_support():

            # create a new figure
            fi, ax = plt.subplots(1, 1)

            # save to frames of an animation
            with wri.saving(fi, filename, 400):

                # loop over zenith angles
                for z in tqdm(zenith_angles):

                    # clear whatever was in the previous axes
                    ax.cla()

                    # plot the sunset at this current stellar zenith angle
                    self.plot_sunset(
                        ax=ax,
                        zenith_angle=z,
                        maxelevation=maxelevation,
                        minelevation=minelevation,
                    )

                    # grab this frame in the animation
                    wri.grab_frame()

    def animate_everything(
        self,
        filename="everything-sunset.mp4",
        ingredients=["sky", "sky-zoom", "rgb", "spectrum"],
        maxelevation=50 * u.deg,
        skyresolution=0.5 * u.deg,
        skynormalization=0.7,
        motionresolution=0.5 * u.deg,
        overwrite=False,
        **kwargs,
    ):

        """
        Make an animation that summarizes
        lots of information about a sunset.

        Parameters
        ----------
        filename : string
            The filename where this animation should be saved
        ingredients : list
            The components to draw. These are:
                'sky' = a slice of many elevations
                'sky-zoom' = just the disk of the star
                'rgb' = a plot of relative RGB colors
                'spectrum' = a plot of the spectrum
            Any subset of these is allowed.
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
        **kwargs : dict
            All extra keyword arguments will be passed to
            plot_everything.
        """

        # get the appropriate animation writer
        wri = get_writer(filename)

        # calculate the grid of zenith angles for the star
        min_zenith = 90 * u.deg - maxelevation - self.source.angular_size()
        minelevation = -(self.source.angular_size() + motionresolution)
        max_zenith = 90 * u.deg - minelevation
        zenith_angles = (
            np.arange(
                min_zenith.to("deg").value,
                max_zenith.to("deg").value,
                motionresolution.to("deg").value,
            )
            * u.deg
        )

        if os.path.exists(filename):
            if overwrite is False:
                message = f"""
                Filename {filename} already exists.
                Not recreating the movie, unless
                you rerun with `overwrite=True`.
                """
                warnings.warn(message)
                return

        with plt.style.context("dark_background"), quantity_support():

            fi = self.plot_everything(zenith_angles[0], **kwargs)

            # save to frames of an animation
            with wri.saving(fi, filename, fi.get_dpi()):

                # loop over zenith angles
                for z in tqdm(zenith_angles):

                    # clear whatever was in the previous axes
                    fi.clf()

                    # plot the sunset at this current stellar zenith angle
                    self.plot_everything(
                        fi=fi,
                        zenith_angle=z,
                        ingredients=ingredients,
                        **kwargs,
                    )

                    # grab this frame in the animation
                    wri.grab_frame()

    def plot_everything(
        self,
        ingredients=["sky", "sky-zoom", "rgb", "spectrum"],
        zenith_angle=85 * u.deg,
        maxelevation=50 * u.deg,
        fi=None,
        pixels=[1920, 1080],
        width=8,
        filename=None,
    ):
        """
        Make an animation that summarizes
        lots of information about a sunset.

        Parameters
        ----------
        ingredients : list
            The components to draw. These are:
                'sky' = a slice of many elevations
                'sky-zoom' = just the disk of the star
                'rgb' = a plot of relative RGB colors
                'spectrum' = a plot of the spectrum
            Any subset of these is allowed.
        zenith_angle : astropy.units.quantity.Quantity
            The angle away from zenith.
        maxelevation : astropy.units.quantity.Quantity
            The maximum elevation above the horizon
            (for setting the plot limits.)
        aspect : list
            The aspect ratio [width, height] of the
            size of the plot.
        """

        with plt.style.context("dark_background"), quantity_support():
            xpixels, ypixels = pixels
            dpi = xpixels / width
            scale = 12 / width
            if fi is None:

                fi = plt.figure(
                    figsize=(width, ypixels * width / xpixels),
                    dpi=dpi,
                )

            # set up the basic columns
            gs = GridSpec(
                1,
                4,
                width_ratios=[0.18, 0.4, 0.09 * scale, 0.9],
                wspace=0.0, figure=fi, 
            )

            ax = {}

            # set up the contextualizing plots on the left
            gs_geometry = GridSpecFromSubplotSpec(
                2,
                1,
                height_ratios=[1, 1],
                hspace=0.5 * scale,
                subplot_spec=gs[1],
            )

            # ax['cartoon'] = plt.subplot(gs_geometry[0])
            if "sky-zoom" in ingredients:
                ax["sky-zoom"] = plt.subplot(gs_geometry[0])
                self.plot_disk(ax=ax["sky-zoom"], zenith_angle=zenith_angle)
                width = 1.1 * self.source.angular_size()
                plt.xlim(-width, width)
                plt.ylim(
                    90 * u.deg - zenith_angle - width,
                    90 * u.deg - zenith_angle + width,
                )
                ax["sky-zoom"].get_xaxis().set_visible(False)
                ax["sky-zoom"].get_yaxis().set_visible(False)
                plt.axis("scaled")
                plt.axis("off")

            # set up the main panel showing the sunset over the full sky
            if "sky" in ingredients:
                ax["sky"] = plt.subplot(gs[0])
                # actually make the plots
                self.plot_sunset(
                    ax=ax["sky"],
                    zenith_angle=zenith_angle,
                    maxelevation=maxelevation,
                )

            #
            # set up the spectrum panels
            gs_spectra = GridSpecFromSubplotSpec(
                2,
                1,
                height_ratios=[1, 1],
                subplot_spec=gs[3],
                hspace=0.5 * scale,
            )
            if "rgb" in ingredients:
                ax["rgb"] = plt.subplot(gs_spectra[0])
                self.plot_rgb(ax=ax["rgb"])
                plt.ylim(0, 120)
                plt.xlim(360 * u.nm, 740 * u.nm)

            if "spectrum" in ingredients:
                ax["spectrum"] = plt.subplot(
                    gs_spectra[1],
                    sharex=ax["rgb"],
                    sharey=ax["rgb"],
                )
                self.plot_as_rainbow(ax["spectrum"])
                plt.ylim(0, 120)
                plt.xlim(360 * u.nm, 740 * u.nm)

        if filename is not None:
            plt.savefig(filename, facecolor="black")
        return fi

    def plot_rgb(self, ax=None, **kwargs):
        """
        Wrapper to add the unextincted spectrum
        in dark gray in the background for context.
        """
        ax = Spectrum.plot_rgb(self, ax=ax, **kwargs)
        self.source.plot_rgb(foreground=False, ax=ax)

    def plot_as_rainbow(self, ax=None, **kwargs):
        """
        Wrapper to add the unextincted spectrum
        in dark gray in the background for context.
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            ax = Spectrum.plot_as_rainbow(self, ax=ax, **kwargs)
            warnings.simplefilter("ignore")
            self.source.plot_as_rainbow(foreground=False, ax=ax)
