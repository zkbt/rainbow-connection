from ..sources import Spectrum

class TransmittedSpectrum(Spectrum):

    def __init__(self, source, atmosphere):
        self.source = source
        self.atmosphere = atmosphere

    def __repr__(self):
        s = repr(self.source)
        t = repr(self.atmosphere)
        return f'<{s} through {t}>'

    def set_elevation_angle(self, *args, **kwargs):
        self.atmosphere.set_elevation_angle(*args, **kwargs)

    def set_zenith_angle(self, *args, **kwargs):
        self.atmosphere.set_zenith_angle(*args, **kwargs)

    def spectrum(self, wavelength=None, zenith_angle=None):

        # the transmission wavelength dominates over the source
        w = self.atmosphere.wavelength(wavelength)

        # get the transmission
        t = self.atmosphere.transmission(wavelength=w,
                                         zenith_angle=zenith_angle)

        # get the original spectrum
        s = self.source.spectrum(wavelength=w)

        # return the transmitted spectrum
        return s*t
