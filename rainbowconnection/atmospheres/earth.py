from .atmosphere import *

class Earth(Atmosphere):
    def __init__(self):

        # load the vertical transmission through Earth's atmosphere
        filename = os.path.join(data_directory, 'earthtransmission.txt')

        # store th
        d = ascii.read(filename, comment='#')
        self._wavelength = d['wavelength'].data*1e3
        self._transmission = d['transmission'].data
        self._tau_zenith = -np.log(self._transmission)
        self.default_wavelengths = self._wavelength*u.nm
        self.H = 8*u.km
        self.radius = 1*u.Rearth
        self.fortney_factor = np.sqrt(2*np.pi*self.radius/self.H).decompose()

        self.set_zenith_angle()

    def transmission(self, wavelength, zenith_angle=None):

        # update the zenith angle, if necessary
        if zenith_angle is not None:
            self.set_zenith_angle(zenith_angle)

        # make sure at least some grid of wavelengths is defined
        w = self.wavelength(wavelength)

        # figure out the transmission at this altitude
        # FIXME -- this is a major kludge! do the integral!
        effective_airmass = np.minimum(self.fortney_factor, self.airmass)
        tau = self._tau_zenith*effective_airmass


        # bin this spectrum to the particular wavelength grid
        # FIXME: binning choice gets real scary with transmission
        neww, newt = bintogrid(self._wavelength, tau,
                         newx=w.to('nm').value,
                         drop_nans=False)

        # make sure the wavelengths match up
        assert(np.all(neww == w.to('nm').value))

        # return the binned transmission
        return np.exp(-newt)
