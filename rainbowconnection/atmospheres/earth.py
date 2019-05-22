from .atmosphere import *

class Earth(DiscreteAtmosphere):
    '''
    Earth's transmission spectrum.
    '''

    def read_transmission(self, **kwargs):
        '''
        Read Earth's transmission spectrum in from a file
        and define the required attributes:

            ._wavelengths (with units of wavelength)
            ._tau_zenith_reference (unitless)
            .H (with units of length)
            .radius (with units of length)
        '''

        # load the ESO vertical transmission through Earth's atmosphere
        filename = os.path.join(data_directory, 'earthtransmission.txt')
        d = ascii.read(filename, comment='#')

        # calculate the optical depth at zenith
        self._wavelength = d['wavelength'].data*1e3*u.nm
        self._transmission_zenith = d['transmission'].data
        self._tau_zenith_reference = -np.log(self._transmission_zenith)

        # define geometry of the atmosphere (how spherical?)
        mu = 29
        T = 273*u.K
        g = 9.8*u.m/u.s**2
        self.H = (con.k_B*T/mu/g/con.m_p).to('km')
        self.radius = 1*u.Rearth
