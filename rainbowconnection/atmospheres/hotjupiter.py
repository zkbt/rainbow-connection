from .atmosphere import *


class HotJupiter(DiscreteAtmosphere):
    def read_transmission(self, TiO=False, **kwargs):
        """
        Read a Hot Jupiter's transmission spectrum in from a file
        and define the required attributes:

            ._wavelengths (with units of wavelength)
            ._tau_zenith_reference (unitless)
            .H (with units of length)
            .radius (with units of length)
        """

        # load a transmission spectrum from Fortney et al. (2010)
        if TiO:
            choice = "lambda_1500K_g10_wTiOVO.dat"
        else:
            choice = "lambda_1500K_g10_noTiOVO.dat"
        filename = os.path.join(data_directory, f"fortney/{choice}")
        d = ascii.read(filename, comment="#")

        # define the wavelength grid
        self._wavelength = d["wavelength"].data * 1e3 * u.nm

        # store the effective transit radius
        radiusinkm = d["radiusinkm"].data * u.km
        self._transit_radius = radiusinkm.to("Rjup")

        # pick an (arbitrary-ish) reference radius
        reference_radius = np.min(radiusinkm)

        # define geometry of the atmosphere (how spherical?)
        mu = 2.3
        self.T = 1500 * u.K
        g = 10 * u.m / u.s**2
        self.H = (con.k_B * self.T / mu / g / con.m_p).to("km")
        self.radius = reference_radius

        # calculate the optical depth at zenith
        z_over_H = ((radiusinkm - reference_radius) / self.H).decompose()
        tau_slant = np.exp(z_over_H)
        self._tau_zenith_reference = (tau_slant / self.fortney_factor()).decompose()

        # set the default wavelengths for simple plots
        self.wavelength = self._wavelength
