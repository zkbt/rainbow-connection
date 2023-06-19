from .imports import *
from .colortools import plot_simple_rainbow

colors_to_wavelengths = dict(
    R=625,
    O=595,
    Y=570,
    G=520,
    B=465,
    I=448,
    V=430,
)


def setup_axes_with_rainbow(
    ax=None, rainbow=True, figsize=None, height_ratio=0.1, roygbiv=False
):
    # create new ax(s), unless we're supposed to over plot on one
    if ax is None:
        fi = plt.figure(figsize=figsize)
        # decide whether to add a horizontal rainbow cartoon
        if rainbow:
            # create a two-part grid
            gs = GridSpec(2, 1, height_ratios=[height_ratio, 1], hspace=0.0, figure=fi)

            # plot a cartoon rainbow, in a box above
            ax_rainbow = plt.subplot(gs[0])
            plt.sca(ax_rainbow)
            plot_simple_rainbow(ax=ax_rainbow)
            if roygbiv:
                for k, v in colors_to_wavelengths.items():
                    plt.text(v, 1, k, ha="center", va="bottom", color="gray")

            ax_rainbow.get_yaxis().set_visible(False)
            ax_rainbow.get_xaxis().set_visible(False)
            ax_rainbow.set_facecolor("black")
            # (kludge to ensure black background for rainbow)
            if plt.gcf().get_facecolor() == (0.0, 0.0, 0.0, 1.0):
                plt.axis("off")

            # create the main ax
            ax = plt.subplot(gs[1], sharex=ax_rainbow)

        else:
            # simply create one simple ax
            ax = plt.subplot()

    # make sure we point back at the first plot
    plt.sca(ax)
    return ax
