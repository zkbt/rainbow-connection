"""
A few shortcuts for dealing with units.
"""
import astropy.units as u

# define some useful units
spectral_luminosity_unit = u.W / u.nm
spectral_flux_unit = u.W / u.nm / u.m**2
spectral_intensity_unit = u.W / u.nm / u.m**2 / u.sr

luminosity_unit = u.W
flux_unit = u.W / u.m**2
intensity_unit = u.W / u.m**2 / u.sr

unit2name = {
    spectral_flux_unit: "Flux",
    spectral_luminosity_unit: "Luminosity",
    flux_unit: "Flux",
    luminosity_unit: "Luminosity",
    intensity_unit: "Intensity",
    spectral_intensity_unit: "Intensity",
}
# name2unit = {v.si:k for k, v in unit2name.items()}


def determine_quantity(unit):
    """
    Determine the name of a particular unit.

    Parameters
    ----------
    unit : astropy.units.core.Unit, astropy.units.core.CompositeUnit

    Returns
    -------
    """

    # try all the listed units
    for k in unit2name:
        if unit.is_equivalent(k):
            return unit2name[k]

    # if none match, return "?"
    return "?"
