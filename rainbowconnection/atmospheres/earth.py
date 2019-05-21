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
        self.set_zenith_angle()

    def transmission(self, wavelength, zenith_angle=None):

        # update the zenith angle, if necessary
        if zenith_angle is not None:
            self.set_zenith_angle(zenith_angle)

        # make sure at least some grid of wavelengths is defined
        w = self.wavelength(wavelength)

        # figure out the transmission at this altitude
        tau = self._tau_zenith*self.airmass


        # bin this spectrum to the particular wavelength grid
        # FIXME: binning choice gets real scary with transmission
        neww, newt = bintogrid(self._wavelength, tau,
                         newx=w.to('nm').value,
                         drop_nans=False)

        # make sure the wavelengths match up
        assert(np.all(neww == w.to('nm').value))

        # return the binned transmission
        return np.exp(-newt)
