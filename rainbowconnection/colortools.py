"""
These tools are pulled *very* strongly from the great
ingredients in the `colour` package. They've been
modified and simplified here, in order to work better
with the design of rainbowconnection.
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon
from six.moves import reduce

from colour.algebra import LinearInterpolator
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
from colour.utilities import (
    domain_range_scale,
    first_item,
    normalise_maximum,
    tstack,
    ColourRuntimeWarning,
)

from colour import SpectralDistribution

# define a standard set of color-matching functions
CMFs = first_item(filter_cmfs("CIE 1931 2 Degree Standard Observer").values())


def plot_simple_rainbow(
    ax=None,
    wavelength=None,
    flux=None,
    cmfs="CIE 1931 2 Degree Standard Observer",
    **kwargs
):
    """
    Plot a simple horizontal rainbow,
    based off colour's plot_single_sd.

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

    # pull out the CMFss
    cmfs = first_item(filter_cmfs(cmfs).values())

    if wavelength is None:
        # create a grid of wavelengths (at which CMFss are useful)
        wavelength = cmfs.wavelengths

    ok = (wavelength >= np.min(cmfs.wavelengths)) & (
        wavelength <= np.max(cmfs.wavelengths)
    )

    # create colors at theose wavelengths
    colours = XYZ_to_plotting_colourspace(
        wavelength_to_XYZ(wavelength[ok], cmfs),
        CCS_ILLUMINANTS["CIE 1931 2 Degree Standard Observer"]["E"],
    )

    # normalize the colors to their maximum?
    colours = CONSTANTS_COLOUR_STYLE.colour.colourspace.cctf_encoding(
        normalise_maximum(colours)
    )

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
                tstack([wavelength[ok], flux[ok]]),
                (x_max, 0),
            ]
        ),
        facecolor="none",
        edgecolor="none",
    )
    ax.add_patch(polygon)

    # draw bars, with the colors at each vertical stripe
    padding = 0.2
    ax.bar(
        x=wavelength[ok] - padding / 2,
        height=max(flux[ok]),
        width=1 + padding,
        color=colours,
        align="edge",
    )

    return ax


def plot_with_rainbow_fill(
    ax=None,
    wavelength=None,
    flux=None,
    cmfs="CIE 1931 2 Degree Standard Observer",
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

    # pull out only the values that *can* be converted to colors
    ok = (wavelength >= np.min(CMFs.wavelengths)) & (
        wavelength <= np.max(CMFs.wavelengths)
    )

    w, f = wavelength[ok], flux[ok]

    # clip the top of the box?
    f = np.minimum(f, rainbowtop)

    XYZ = wavelength_to_XYZ(w)

    # create colors at theose wavelengths
    colours = XYZ_to_plotting_colourspace(XYZ)

    # normalize the colors to their maximum?
    # colours = CONSTANTS_COLOUR_STYLE.colour.colourspace.cctf_encoding(
    #    normalise_maximum(colours))
    colours = np.maximum(0, colours / np.max(colours))

    x_min, x_max = min(w), max(w)
    y_min, y_max = 0, max(f) + max(f) * 0.05

    # create a polygon to define the top of the bars?
    polygon = Polygon(
        np.vstack(
            [
                (x_min, 0),
                tstack([w, f]),
                (x_max, 0),
            ]
        ),
        facecolor="none",
        edgecolor="none",
    )
    ax.add_patch(polygon)

    # draw bars, with the colors at each vertical stripe
    padding = 0.0
    dw = np.mean(np.diff(w))
    ax.bar(
        x=w,
        height=f,
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
    cmfs="CIE 1931 2 Degree Standard Observer",
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

    cmfs : string
        The color matching function(s?) to use.

    vector : bool


    Returns
    -------
    ax : matplotlib.axes._subplots.AxesSubplot
        The ax into which the rainbow was drawn.
    """

    # make sure our plotting ax is defined
    if ax is None:
        ax = plt.gca()

    # make sure we have a grid of wavelengths defined
    if wavelength is None:
        # create a grid of wavelengths (at which CMFss are useful)
        wavelength = CMFs.wavelengths

    # create y values that will be plotted (these could be spectrum)
    if flux is None:
        flux = np.ones_like(wavelength)

    # pull out only the values that *can* be converted to colors
    ok = (wavelength >= np.min(CMFs.wavelengths)) & (
        wavelength <= np.max(CMFs.wavelengths)
    )
    w, f = wavelength[ok], flux[ok]

    # get the XYZ for the wavelengths
    XYZ = wavelength_to_XYZ(w)

    # create colors at those wavelengths
    colours = XYZ_to_plotting_colourspace(XYZ)

    # normalize the colors to their maximum?
    # colours = CONSTANTS_COLOUR_STYLE.colour.colourspace.cctf_encoding(
    #    normalise_maximum(colours))

    # normalize the brightness by the flux
    colours *= f[:, np.newaxis]

    # normalize to the brightest line
    colours = np.maximum(0, colours / np.max(colours))

    # draw as an RGB color image with imshow
    ax.imshow(
        colours[np.newaxis, :, :],
        aspect="auto",
        extent=[np.min(w), np.max(w), 0, 1],
        interpolation="nearest",
    )

    return ax
