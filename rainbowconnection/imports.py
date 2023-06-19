from IPython import get_ipython
try:
    if get_ipython().__class__.__name__ == 'ZMQInteractiveShell':
        get_ipython().magic(u'matplotlib inline')
except AttributeError:
    pass

import copy, inspect, os, glob, warnings
package_directory = os.path.dirname(inspect.getfile(inspect.currentframe()))
data_directory = os.path.join(package_directory, 'data')
package_name = "rainbow-connection"

import numpy as np
import matplotlib.pyplot as plt

plt.matplotlib.rcParams["figure.figsize"] = (8, 3)
plt.matplotlib.rcParams["figure.dpi"] = 100

from astropy.io import ascii
import astropy.units as u
import astropy.constants as con
from astropy.visualization import quantity_support
from astropy.utils.data import (
    download_file,
    cache_contents,
    check_download_cache,
    zipfile,
)


from .resampletools import bintogrid

from tqdm import tqdm

from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from scipy.integrate import quad


default_wavelength_grid = np.linspace(100, 1000, 901) * u.nm

