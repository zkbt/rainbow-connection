import copy, pkg_resources, os, glob, warnings

data_directory = pkg_resources.resource_filename("rainbowconnection", "data")
package_name = "rainbow-connection"

import numpy as np
import matplotlib.pyplot as plt

plt.matplotlib.rcParams["figure.figsize"] = (8, 3)
plt.matplotlib.rcParams["figure.dpi"] = 300

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
