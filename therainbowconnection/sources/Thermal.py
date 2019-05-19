from .Light import *

class Thermal(Light):
    def __init__(self, teff=5800*u.K, radius=1*u.Rsun):

        self.teff = teff
        self.radius = radius

    def surface_area(self):
        return 4*np.pi*self.radius**2

    def intensity(self, wavelength):
        '''
        This function calculates the thermal emission intensity spectrum of a surface.

            Inputs:
                wavelength = numpy array of wavelengths (with astropy units)

            Outputs:
                Returns an array of thermal emission intensities,
                in astropy units of W/(m^2*micron*sr). This is a flux, which has
                already been integrated over solid angle.
        '''

        temperature = self.teff

        # define variables as shortcut to the constants we need
        h = constants.h
        k = constants.k_B
        c = constants.c

        # this is the thing that goes into the exponent (it's units better cancel!)
        up = h*c/(wavelength*k*temperature)

        # calculate the intensity from the Planck function
        intensity = (2*h*c**2/wavelength**5/(np.exp(up) - 1))/u.steradian

        # return the intensity
        return intensity.to('W/(m**2*nm*sr)')

    def flux(self, wavelength):
        '''
        This function calculates the thermal emission flux spectrum of a surface.

            Inputs:
                wavelength = numpy array of wavelengths (with astropy units)
                temperature = a single number, the temperature (with astropy units)

            Outputs:
                Returns an array of thermal emission fluxes,
                in astropy units of W/(m^2*micron). This is a flux, which has
                already been integrated over solid angle.
        '''

        # calculate the flux, knowing the angle integral will be pi steradians (for isotropic emission)
        flux = self.intensity(wavelength)*np.pi*u.steradian

        # return the flux, in convenient units
        return flux.to('W/(m**2*nm)')
