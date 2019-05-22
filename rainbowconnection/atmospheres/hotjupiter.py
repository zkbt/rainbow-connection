from .earth import *

class HotJupiter(Earth):
    def __init__(self, altitude=0):

        # load the vertical transmission through Earth's atmosphere
        filename = os.path.join(data_directory, 'fortney/lambda_1500K_g10_noTiOVO.dat')

        # store th
        d = ascii.read(filename, comment='#')
        self._wavelength = d['wavelength'].data*1e3
        radiusinkm = d['radiusinkm'].data*u.km
        reference_radius = np.min(radiusinkm)
        T = 1500*u.K
        mu = 2.3
        g = 10*u.m/u.s**2
        self.H = (con.k_B*T/mu/g/con.m_p).to('km')
        self.radius = reference_radius
        self.fortney_factor = np.sqrt(2*np.pi*self.radius/self.H).decompose()

        z_over_H = ((radiusinkm - reference_radius)/self.H).decompose() - altitude
        self._tau_zenith = (np.exp(z_over_H)*np.sqrt(self.H/2/np.pi/self.radius)).decompose()

        self._transmission = np.exp(-self._tau_zenith)

        #self._tau_zenith = -np.log(self._transmission)
        self.default_wavelengths = self._wavelength*u.nm

        self.set_zenith_angle()
