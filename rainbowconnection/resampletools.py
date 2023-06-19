"""Tools for resampling array from grid of independent variables to another."""

import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt


def binsizes(x):
    """If x is an array of bin centers, calculate what their sizes are.
    (assumes outermost bins are same size as their neighbors)"""

    binsize = np.zeros_like(x)
    binsize[0:-1] = x[1:] - x[0:-1]
    binsize[-1] = binsize[-2]
    return binsize


def plotboxy(x, y, **kwargs):
    """
    Plot with boxes, to show the left and right edges of a box. This is useful,
    for example, to plot flux associated with pixels, in case you are trying to
    do a sub-pixel resample or interpolation or shift.

    (kwargs are passed on to plt.plot)
    """

    # what are the edges of the bins (making a guess for those on the ends)
    xbinsize = binsizes(x)
    xleft = x - xbinsize / 2.0
    xright = x + xbinsize / 2.0

    # create a array that doubles up the y values, and interleaves the edges
    plot_x = np.vstack((xleft, xright)).reshape((-1,), order="F")
    plot_y = np.vstack((y, y)).reshape((-1,), order="F")

    # plot those constructed arrays
    plt.plot(plot_x, plot_y, **kwargs)


def fluxconservingresample(
    xin_unsorted,
    yin_unsorted,
    xout,
    visualize=False,
    demo=False,
    treatnanas=0.0,
):
    """
    Starting from some initial x and y, resample onto a different grid
    (either higher or lower resolution), while conserving total flux.

    When including the entire range of xin, sum(yout) == sum(yin) should be true.
    """

    # sort to make sure x is strictly increasing
    s = np.argsort(xin_unsorted)
    xin = xin_unsorted[s]
    yin = yin_unsorted[s]

    # set up the bins, to calculate cumulative distribution of y?
    xinbinsize = binsizes(xin)
    xinleft = xin - xinbinsize / 2.0
    xinright = xin + xinbinsize / 2.0

    # the first element should be the left edge of the first pixel
    # last element will be right edge of last pixel
    xinforcdf = np.hstack([xinleft, xinright[-1]])

    # to the left of the first pixel, assume flux is zero
    yinforcdf = np.hstack([0, yin])

    # correct for non-finite
    bad = np.isnan(yinforcdf)
    yinforcdf[bad] = treatnanas

    # calculate the cumulative distribution function of the flux (at pixel edge locations)
    cdfin = np.cumsum(yinforcdf)

    # create an interpolator for that
    cdfinterpolator = scipy.interpolate.interp1d(
        xinforcdf,
        cdfin,
        kind="linear",
        bounds_error=False,
        fill_value=(0.0, np.sum(yin)),
    )

    # calculate bin edges (of size len(xout)+1)
    xoutbinsize = binsizes(xout)
    xoutleft = xout - xoutbinsize / 2.0
    xoutright = xout + xoutbinsize / 2.0
    xoutcdf = np.hstack([xoutleft, xoutright[-1]])

    # interpolate the CDF onto those bin edges
    cdfout = cdfinterpolator(xoutcdf)

    # take the derivative of the CDF to get the flux per resampled bin
    # (xout is the center of the bin, and yout is the flux in that bin)
    yout = np.diff(cdfout)

    if visualize:
        fi, (ax_cdf, ax_pdf) = plt.subplots(2, 1, sharex=True, figsize=(9, 6), constrained_layout=True)
        inkw = dict(
            color="black",
            alpha=1,
            linewidth=3,
            marker=".",
            markeredgecolor="none",
        )
        outkw = dict(
            color="darkorange",
            alpha=1,
            linewidth=1,
            marker=".",
            markersize=8,
            markeredgecolor="none",
        )

        legkw = dict(
            fontsize=10,
            frameon=False,
            loc="upper left",
            bbox_to_anchor=(1, 1),
        )

        # plot the PDFs
        plt.sca(ax_pdf)
        plt.ylabel("Flux per (Original) Pixel")
        plt.xlabel("Pixel")
        # plot the original pixels (in df/dpixel to compare with resampled)
        plotboxy(xin, yin / xinbinsize, label="Original Pixels", **inkw)

        # what would a bad interpolation look like?
        badinterpolation = scipy.interpolate.interp1d(
            xin,
            yin / xinbinsize,
            kind="linear",
            bounds_error=False,
            fill_value=0.0,
        )
        plt.plot(
            xout,
            badinterpolation(xout),
            color="cornflowerblue",
            alpha=1,
            linewidth=1,
            marker=".",
            markersize=8,
            markeredgecolor="none",
            label="Silly Simple Interpolation",
        )

        # plot the flux-conserving resampled data (again, in df/d"pixel")
        plt.plot(
            xout, yout / xoutbinsize, label="Flux-Conserving Interpolation", **outkw
        )
        plt.legend(**legkw)

        # plot the CDFs
        plt.sca(ax_cdf)
        plt.ylabel("Cumulative Flux (from left)")

        # plot the original CDF
        plt.plot(xinforcdf, cdfin, label="Original Pixels", **inkw)

        # plot the interpolated CDF
        plt.plot(xoutcdf, cdfout, label="Flux-Conserved Resample", **outkw)
        # plt.legend(**legkw)
        if demo:
            a = raw_input(
                "Pausing a moment to check on interpolation; press return to continue."
            )

        print("{:>6} = {:.5f}".format("Actual", np.sum(yin)))
        print(
            "{:>6} = {:.5f}".format(
                "Silly",
                np.sum(badinterpolation(xout) * xoutbinsize),
            )
        )
        print("{:>6} = {:.5f}".format("CDF", np.sum(yout)))

    # return the resampled y-values
    return yout


def bintogrid(
    x,
    y,
    unc=None,
    newx=None,
    dx=None,
    weighting="inversevariance",
    drop_nans=True,
):
    """
    x = the independent variable (wavelength)
    y = the measurement (transit depth)
    unc = the uncertainty on the measurement (sigmas on the transit depths)
    newx = a (linearly) uniformly spaced grid onto which we should bin
    dx = if newx isn't set, then create a grid with this fixed spacing

    Right now, this works only on 1-dimenions arrays,
    but in the long term we should update it to work in
    for arbitrary N-dimensional arrays.

    Arrays must be ordered from decreasing to increasing x.
    """

    # make up a grid, if one wasn't specified
    if newx is None:
        newx = np.arange(np.min(x), np.max(x) + dx, dx)

    # don't complain about zero-divisions in here (just make uncertainties -> infinity)
    with np.errstate(divide="ignore", invalid="ignore"):

        # resample the sums onto that new grid (in logarithmic space)
        if unc is None:
            weights = np.ones_like(x)
        else:
            if weighting == "inversevariance":
                weights = 1 / unc**2
        numerator = fluxconservingresample(x, y * weights, newx, visualize=False)
        denominator = fluxconservingresample(x, weights, newx, visualize=False)

        # the binned weighted means on the new grid
        newy = numerator / denominator

        # the standard error on the means, for those bins
        newunc = np.sqrt(1 / denominator)

    if drop_nans:
        ok = np.isfinite(newy)
    else:
        ok = np.ones_like(newx).astype(bool)
    # return binned x, y, unc
    if unc is None:
        return newx[ok], newy[ok]
    else:
        return newx[ok], newy[ok], newunc[ok]


def bintoR(x, y, unc=None, R=50, xlim=None, weighting="inversevariance"):
    """
    x = the independent variable (wavelength)
    y = the measurement (transit depth)
    unc = the uncertainty on the measurement (sigmas on the transit depths)
    R = the resolution (lambda/dlambda) to bin everything to
    xlim = the limits of the new binned grid
                    [if None, will be created from data]

    Right now, this works only on 1-dimenions arrays,
    but in the long term we should update it to work in
    for arbitrary N-dimensional arrays.

    Arrays must be ordered from decreasing to increasing x.
    """

    # create a new grid of x at the given resolution
    lnx = np.log(x)
    dnewlnx = 1.0 / R

    # set the limits of the new xgrid (in log space)
    if xlim is None:
        # use the input grid to set the limits
        lnxbottom, lnxtop = np.min(lnx), np.max(lnx)
    else:
        # use the custom xlim to set the limits
        lnxbottom, lnxtop = xlim

    # create a new, log-uniform grid of x values
    newlnx = np.arange(lnxbottom, lnxtop + dnewlnx, dnewlnx)

    # now do the binning on a uniform grid of lnx

    if unc is None:
        blnx, by = bintogrid(lnx, y, unc, newx=newlnx, weighting=weighting)
        return np.exp(blnx), by
    else:
        blnx, by, bunc = bintogrid(lnx, y, unc, newx=newlnx, weighting=weighting)
        return np.exp(blnx), by, bunc


# swiped from stack overflow
def find_nearest(array, value, verbose=False):
    """
    Find the value of the element in the array that is closest
    to the provided value.
    """
    idx = (np.abs(np.array(array) - value)).argmin()
    if verbose:
        print("{0} --> {1}".format(value, array[idx]))
    return array[idx]


# modified from above
def find_two_nearest(array, value, verbose=False):
    """
    Find the two values of the two elements in the array
    that most closely span the provided value.
    """
    # assumes ordered arrays and that value falls between the min and max of the array
    offset = value - np.array(array)
    signs = np.sign(offset)

    # handle 1-element arrays
    if len(array) == 1:
        return [array[0], array[0]]

    if (signs == -1).all() | (signs == 1).all():
        # value is below the minimum of the array
        m = np.argmin(np.abs(offset))
        left, right = m, m
    else:
        # the value is somewhere in the bounds between the array's min and max
        left = (signs[1:] - signs[0:-1]).nonzero()[0][0]
        right = left + 1

    nearest = [array[left], array[right]]

    if verbose:
        print("{0} is closest to {1}".format(value, nearest))
    return nearest


def interpolation_weights(bounds, value, verbose=True):
    """
    Figure out interpolation weights to apply to
    some value that falls somewhere between two bounds.
    """
    if bounds[0] == bounds[1]:
        return 1.0, 0.0
    assert (value >= np.min(bounds)) * (value <= np.max(bounds))
    span = float(bounds[1] - bounds[0])
    weights = [
        (bounds[1] - value) / span,
        (value - bounds[0]) / span,
    ]
    return weights
