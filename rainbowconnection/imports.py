import numpy as np
import matplotlib.pyplot as plt
import copy, pkg_resources, os, glob, warnings

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


data_directory = pkg_resources.resource_filename("rainbowconnection", "data")

package_name = "rainbow-connection"
