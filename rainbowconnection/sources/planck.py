from .spectra import *



        
class Thermal(Spectrum):
    def __init__(self, teff=5800 * u.K, radius=1 * u.Rsun):

        self.teff = teff
        self.radius = radius
        self.wavelength = np.logspace(2, 3, 1000) * u.nm

    def intensity(self, wavelength):
        """
        This function calculates the thermal emission intensity spectrum of a surface.

            Inputs:
                wavelength = numpy array of wavelengths (with astropy units)

            Outputs:
                Returns an array of thermal emission intensities,
                in astropy units of W/(m^2*micron*sr). This is a flux, which has
                already been integrated over solid angle.
        """

        temperature = self.teff

        # define variables as shortcut to the constants we need
        h = con.h
        k = con.k_B
        c = con.c

        # this is the thing that goes into the exponent (it's units better cancel!)
        up = h * c / (wavelength * k * temperature)

        # calculate the intensity from the Planck function
        intensity = (2 * h * c**2 / wavelength**5 / (np.exp(up) - 1)) / u.steradian

        # return the intensity
        return intensity.to("W/(m**2*micron*sr)")


    def surface_flux(self, wavelength):
        """
        This function calculates the thermal emission flux spectrum of a surface.

            Inputs:
                wavelength = numpy array of wavelengths (with astropy units)
                temperature = a single number, the temperature (with astropy units)

            Outputs:
                Returns an array of thermal emission fluxes,
                in astropy units of W/(m^2*micron). This is a flux, which has
                already been integrated over solid angle.
        """

        # calculate the flux, knowing the angle integral will be pi steradians (for isotropic emission)
        flux = self.intensity(wavelength) * np.pi * u.steradian

        # return the flux, in convenient units
        return flux.to("W/(micron * m**2)")

    def __repr__(self):
        """
        How should this object appear as a string?

        Returns
        -------
        s : str
            A simple string representation.
        """

        basic = f"{self.__class__.__name__} ({self.teff:.0f}, {self.radius})"
        try:
            assert self.distance is not None
            return basic + f" at {self.distance}"
        except (AssertionError, AttributeError):
            return basic

    def integrate(self, lower=None, upper=None):
        """
        Integrate the spectrum over wavelength.

        It gives a number with units identical to the results of
        `.spectrum()` but without the wavelength (W or W/m**2).

        Parameters
        ----------
        lower : astropy.units.quantity.Quantity
            The lower wavelength limit.

        upper : astropy.units.quantity.Quantity
            The lower wavelength limit.

        Returns
        -------
        integral : astropy.units.quantity.Quantity
            The integral over wavelength.
        """

        # if wavelength limits are used, revert back to the numerical integral
        if (lower is not None) or (upper is not None):
            return super().integrate(lower=lower, upper=upper)

        # if there are infinite wavelength limits, do the integral analytically
        surface_flux = con.sigma_sb * self.teff**4
        factor = (self.surface_area() / self.normalization()).decompose()
        return factor * surface_flux


class SumOfThermal(Thermal):
    def intensity(self,wavelength,n_grid=90,plots=False):
        # THANKS PAT!!
        

        temperature = self.teff
        T_dmax = temperature

        # define variables as shortcut to the constants we need
        h = con.h
        k = con.k_B
        c = con.c
        
        # define a 2D grid of latitude and longitude    
        dlat = np.pi / n_grid
        dlon = np.pi / n_grid
        
        lats = np.arange(-np.pi/2., np.pi/2., dlat)
        lons = np.arange(-np.pi, np.pi, dlon)
        lats_mesh, lons_mesh = np.meshgrid(lats, lons, indexing='ij')
        
        # define inclination
        inc = np.pi/2.

        # calculate substellar temperature (f=1) 
        T_sub = T_dmax*(3./2.)**(1./4.)

        # calculate cosine from substellar and temperature at each point in 2D map
        cos_angle_hotspot = np.cos(lats_mesh) * np.cos(lons_mesh)
        T_points = T_sub * cos_angle_hotspot**0.25
        T_points[cos_angle_hotspot <= 0] = 0 # set nightside to 0
        if plots:
            plt.imshow(T_points.value)
            plt.colorbar()
            plt.show()

        # calculate spectrum at each point
        planet_spectrum = (2*h*c**2 / wavelength**5 / (np.exp(h*c/(wavelength * k * T_points[:,:,np.newaxis])) - 1)) / u.steradian
        cos_angle_obs = np.sin(lats_mesh) * np.sin(np.pi/2 - inc) + np.cos(lats_mesh) * np.cos(np.pi/2 - inc) * np.cos(lons_mesh)
        planet_spectrum[cos_angle_obs <= 0] = 0
        unit_solid_angle = np.cos(lats_mesh) * dlat * dlon
        
        # calculate intensity and sum them
        intensity = planet_spectrum * (unit_solid_angle * cos_angle_obs)[:,:,np.newaxis]
        sum_intensity = intensity.sum(axis=(0,1))/np.pi
        if plots:
            plt.imshow(intensity[:,:,0].value)
            plt.colorbar()
            plt.show()

        return sum_intensity.to("W/(m**2*nm*sr)")