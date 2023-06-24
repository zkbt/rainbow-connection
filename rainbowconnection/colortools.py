"""
These tools are pulled *very* strongly from the great
ingredients in the `colour` package. They've been
modified and simplified here, in order to work better
with the design of rainbowconnection.
"""


import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

import numpy as np
from matplotlib.patches import Polygon

# colour-science throws off a lot of warnings; let's mostly ignore those
import warnings

with warnings.catch_warnings():
    from colour.utilities import (
        first_item,
        tstack,
        ColourWarning, ColourRuntimeWarning, ColourUsageWarning
    )

    warnings.simplefilter("ignore", category=ColourWarning)
    warnings.simplefilter("ignore", category=ColourRuntimeWarning)
    warnings.simplefilter("ignore", category=ColourUsageWarning)


    from colour.algebra import LinearInterpolator, normalise_maximum
    from colour.colorimetry import (
        CCS_ILLUMINANTS,
        SDS_ILLUMINANTS,
        LIGHTNESS_METHODS,
        LUMINANCE_METHODS,
        SpectralShape,
        sd_blackbody,
        sd_ones,
        sd_to_XYZ,
        wavelength_to_XYZ,
    )
    from colour.plotting import (
        ColourSwatch,
        CONSTANTS_COLOUR_STYLE,
        XYZ_to_plotting_colourspace,
        artist,
        filter_passthrough,
        filter_cmfs,
        filter_illuminants,
        override_style,
        render,
        plot_single_colour_swatch,
        plot_multi_functions,
        plot_single_sd,
    )


    from colour import SpectralDistribution

    # define a standard set of color-matching functions
    CMFs = first_item(filter_cmfs("CIE 2015 2 Degree Standard Observer").values())

    def plot_simple_rainbow(
        ax=None,
        flux=None,
        **kwargs
    ):
        """
        Plot a simple horizontal rainbow,
        based off colour's plot_single_sd.

        Parameters
        ----------
        ax : matplotlib.axes._subplots.AxesSubplot
            The ax into which the rainbow should be drawn.

        Returns
        -------
        ax : matplotlib.axes._subplots.AxesSubplot
            The ax into which the rainbow was drawn.
        """

        if ax is None:
            ax = plt.gca()


        # create a grid of wavelengths (at which CMFss are useful)
        wavelength = CMFs.wavelengths

        # create colors at theose wavelengths
        colours = XYZ_to_plotting_colourspace(
            wavelength_to_XYZ(wavelength, CMFs),
            CCS_ILLUMINANTS["CIE 1931 2 Degree Standard Observer"]["E"],
        )

        # normalize the colors to their maximum?
        #colours = CONSTANTS_COLOUR_STYLE.colour.colourspace.cctf_encoding(
        #    normalise_maximum(colours)
        #)
        colours = np.maximum(0, colours / np.max(colours))


        # create y values that will be plotted (these could be spectrum)
        if flux is None:
            flux = np.ones_like(wavelength)
        x_min, x_max = min(wavelength), max(wavelength)
        y_min, y_max = 0, max(flux) + max(flux) * 0.05

        # create a polygon to define the top of the bars?
        polygon = Polygon(
            np.vstack(
                [
                    (x_min, 0),
                    tstack([wavelength, flux]),
                    (x_max, 0),
                ]
            ),
            facecolor="none",
            edgecolor="none",
            clip_on=True
        )
        ax.add_patch(polygon)

        # draw bars, with the colors at each vertical stripe
        padding = 0.1
        ax.bar(
            x=wavelength - padding / 2,
            height=max(flux),
            width=1 + padding,
            color=colours,
            align="edge",
            clip_on=True
        )

        return ax

    def plot_with_rainbow_fill(
        ax=None,
        wavelength=None,
        flux=None,
        rainbowtop=np.inf,
        **kwargs
    ):
        """
        (This is still *real* blarg-y*.)
        Plot a spectrum, with a rainbow underneath.

        Parameters
        ----------
        ax : matplotlib.axes._subplots.AxesSubplot
            The ax into which the rainbow should be drawn.
        cmfs : string
            The color matching function(s?) to use.

        Returns
        -------
        ax : matplotlib.axes._subplots.AxesSubplot
            The ax into which the rainbow was drawn.
        """

        if ax is None:
            ax = plt.gca()

        if wavelength is None:
            # create a grid of wavelengths (at which CMFss are useful)
            wavelength = CMFs.wavelengths

        # create y values that will be plotted (these could be spectrum)
        if flux is None:
            flux = np.ones_like(wavelength)


        # get the XYZ for the valid wavelengths
        XYZ_on_CMF_grid = wavelength_to_XYZ(CMFs.wavelengths, CMFs)

        # create colors at those wavelengths
        colours_on_CMF_grid = XYZ_to_plotting_colourspace(XYZ_on_CMF_grid)

        # create imaginary wavelengths as needed
        if wavelength is None:
            wavelength = np.arange(300, 1000)

        # interpolate colors to actual wavelengths
        colour_interpolator = interp1d(CMFs.wavelengths, colours_on_CMF_grid, axis=0, bounds_error=False, fill_value=0)
        colours = colour_interpolator(wavelength)
        colours = np.maximum(0, colours / np.max(colours))

        # create intensities that will be plotted (these could be spectrum)
        if flux is None:
            flux = np.ones_like(wavelength)



  

        x_min, x_max = min(wavelength), max(wavelength)
        y_min, y_max = 0, max(flux) + max(flux) * 0.05

        # create a polygon to define the top of the bars?
        polygon = Polygon(
            np.vstack(
                [
                    (x_min, 0),
                    tstack([wavelength, flux]),
                    (x_max, 0),
                ]
            ),
            facecolor="none",
            edgecolor="none",
            clip_on=True
        )
        ax.add_patch(polygon)

        # draw bars, with the colors at each vertical stripe
        padding = 0.1
        dw = np.mean(np.diff(wavelength))
        ax.bar(
            x=wavelength,
            height=flux,
            width=dw,
            color=colours,
            clip_path=polygon,
            align="edge",
            clip_on=True,
        )

        return ax

    def plot_as_slit_rainbow(
        ax=None,
        wavelength=None,
        flux=None,
        **kwargs
    ):
        """
        (This is still *real* blarg-y*.)
        Plot a spectrum as a light source would be seen through
        a slit spectrometer, with vertical bands of light that
        are brighter or fainter depending on the intensity
        of the spectrum at that particular wavelength.

        Parameters
        ----------
        ax : matplotlib.axes._subplots.AxesSubplot
            The ax into which the rainbow should be drawn.

        wavelength : array
            The wavelengths to include in the spectrum.
            In units of nm, but not as astropy units.

        flux : array
            The fluxes to include in the spectrum.
            In units of whatever, but not as astropy units.


        Returns
        -------
        ax : matplotlib.axes._subplots.AxesSubplot
            The ax into which the rainbow was drawn.
        """

        # make sure our plotting ax is defined
        if ax is None:
            ax = plt.gca()


        # get the XYZ for the valid wavelengths
        XYZ_on_CMF_grid = wavelength_to_XYZ(CMFs.wavelengths, CMFs)

        # create colors at those wavelengths
        colours_on_CMF_grid = XYZ_to_plotting_colourspace(XYZ_on_CMF_grid)

        # create imaginary wavelengths as needed
        if wavelength is None:
            wavelength = np.arange(300, 1000)

        # interpolate colors to actual wavelengths
        colour_interpolator = interp1d(CMFs.wavelengths, colours_on_CMF_grid, axis=0, bounds_error=False, fill_value=0)
        colours = colour_interpolator(wavelength)

        # create intensities that will be plotted (these could be spectrum)
        if flux is None:
            flux = np.ones_like(wavelength)


        # normalize the colors to their maximum?
        # colours = CONSTANTS_COLOUR_STYLE.colour.colourspace.cctf_encoding(
        #    normalise_maximum(colours))

        # normalize the brightness by the flux
        colours *= flux[:, np.newaxis]/np.max(flux)

        # normalize to the brightest line
        colours = np.maximum(0, colours / np.max(colours))

        # draw as an RGB color image with imshow
        dw = np.gradient(wavelength)
        wavelength_edges = np.hstack([wavelength - dw/2, wavelength[-1] + dw[-1]/2])
        slit_edges = np.array([0,1])
        ax.pcolormesh(
            wavelength_edges,
            slit_edges,
            colours[np.newaxis, :, :],
        )
        plt.yticks([])
        return ax
