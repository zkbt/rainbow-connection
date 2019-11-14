from .imports import *
from .colortools import plot_rainbow
def setup_axes_with_rainbow(ax=None, rainbow=True, figsize=(5, 2.5)):
    # create new ax(s), unless we're supposed to over plot on one
    if ax is None:
        plt.figure(figsize=figsize)
        # decide whether to add a horizontal rainbow cartoon
        if rainbow:
            # create a two-part grid
            gs = GridSpec(  2, 1,
                            height_ratios=[0.05, 1],
                            hspace=0.0)


            # plot a cartoon rainbow, in a box above
            ax_rainbow = plt.subplot(gs[0])
            plt.sca(ax_rainbow)
            plot_rainbow(axes=ax_rainbow)

            ax_rainbow.get_yaxis().set_visible(False)
            ax_rainbow.get_xaxis().set_visible(False)
            ax_rainbow.set_facecolor('black')
            # (kludge to ensure black background for rainbow)
            if plt.gcf().get_facecolor() == (0.0, 0.0, 0.0, 1.0):
                plt.axis('off')

            # create the main ax
            ax = plt.subplot(gs[1], sharex=ax_rainbow)

        else:
            # simply create one simple ax
            ax = plt.subplot()


    # make sure we point back at the first plot
    plt.sca(ax)
    return ax
