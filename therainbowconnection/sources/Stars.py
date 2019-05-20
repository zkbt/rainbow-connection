from .Thermal import *
from craftroom.resample import bintogrid
from astropy.io.ascii import read
import pkg_resources, os

data_directory = pkg_resources.resource_filename('therainbowconnection', 'data')

class Sun(Thermal):

    def __init__(self):
        self.radius = 1*u.Rsun
        self.teff = 5777*u.K

        filename = os.path.join(data_directory, 'solarspectrum.txt')
        d = read(filename, comment='#')
        self._wavelength = d['wavelength'].data
        self._flux = d['irradiance'].data
        self.default_wavelengths = self._wavelength*u.nm



    def surface_flux(self, wavelength):

        w, f = bintogrid(self._wavelength,
                         self._flux,
                         newx=wavelength.to('nm').value)

        # convert back to the surface flux from the Sun
        factor = ((1*u.au)**2/(1*u.Rsun)**2).decompose()
        return factor*f*u.W/u.nm/u.m**2
