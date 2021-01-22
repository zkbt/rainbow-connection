from ...imports import *


def stringify_metallicity(Z):
    """
    Convert a metallicity into a PHOENIX-style string.

    Parameters
    ----------
    Z : float
        [Fe/H]-style metallicity (= 0.0 for solar)
    """
    if Z <= 0:
        return "-{:03.1f}".format(np.abs(Z))
    else:
        return "+{:03.1f}".format(Z)
