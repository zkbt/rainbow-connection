from .spectrum import *
from .thermal import Thermal

def gauss(x, x0, sigma):
    '''
    Simple Gaussian.
    '''
    return 1.0/np.sqrt(2*np.pi)/sigma*np.exp(-0.5*((x-x0)/sigma)**2)

class LightBulb(Spectrum):
    def __repr__(self):
        '''
        How should this object appear as a string?

        Returns
        -------
        s : str
            A simple string representation.
        '''

        basic = f'{self.__class__.__name__}LightBulb ({self.power})'
        try:
            assert(self.distance is not None)
            return basic + f' at {self.distance}'
        except (AssertionError, AttributeError):
            return basic

class Incandescent(LightBulb, Thermal):
    def __init__(self, power=100*u.W):

        Thermal.__init__(self, teff=3500*u.K, radius=1*u.mm)

        # renormalize that radius to match the light output
        self.set_power(power)


class Sodium(LightBulb):
    def __init__(self, power=100*u.W):

        # start with a ridiculous radius
        self.radius = 1*u.m

        # renormalize that radius to match the light output
        self.set_power(power)

    def surface_flux(self, wavelength):
        '''
        This function calculates a cartoon spectrum of a sodium lamp.

            Inputs:
                wavelength = numpy array of wavelengths (with astropy units)
                temperature = a single number, the temperature (with astropy units)

            Outputs:
                Returns an array of thermal emission fluxes,
                in astropy units of W/(m^2*micron). This is a flux, which has
                already been integrated over solid angle.
        '''

        # make sure a wavelength grid is defined
        w = self.wavelength(wavelength)
        flux = np.zeros(np.shape(w))/u.nm

        # set the two lines
        lines = np.array([588.9950, 589.5924])*u.nm
        width = 2.0*u.nm

        # give a normalization
        brightness = 1*u.W/u.m**2

        # add the two lines to the spectrum
        for l in lines:
            flux += gauss(w, x0=l, sigma=width)


        # return the flux, in convenient units
        return (brightness*flux).to('W/(nm * m**2)')
