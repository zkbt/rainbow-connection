from .thermal import *

class Sun(Thermal):

    def __init__(self):
        self.radius = 1*u.Rsun
        self.teff = ((u.Lsun/(con.sigma_sb*4*np.pi*u.Rsun**2))**(1.0/4.0)).to(u.K)

        filename = os.path.join(data_directory, 'solarspectrum.txt')
        d = ascii.read(filename, comment='#')
        self._wavelength = d['wavelength'].data
        self._flux = d['irradiance'].data
        self.default_wavelengths = self._wavelength*u.nm



    def surface_flux(self, wavelength=None):

        # make sure at least some grid of wavelengths is defined
        w = self.wavelength(wavelength)

        # bin this spectrum to the particular wavelength grid
        neww, newf = bintogrid(self._wavelength,
                         self._flux,
                         newx=w.to('nm').value,
                         drop_nans=False)

        # make sure the wavelengths match up
        assert(np.all(neww == w.to('nm').value))

        # convert back to the surface flux from the Sun
        factor = ((1*u.au)**2/(1*u.Rsun)**2).decompose()
        return factor*newf*u.W/u.nm/u.m**2
