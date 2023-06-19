import rainbowconnection as rc
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from .directory import *


def test_spectrum():
    N = 100
    w = np.linspace(1, 1000, N) * u.nm
    f = np.random.normal(100, 1, N) * u.W / u.m**2 / u.nm
    s = rc.Spectrum(w, f)
    s.plot()
    save("random.png")


def test_thermal():
    s = rc.Thermal(teff=4000 * u.K, radius=1 * u.m)
    s.plot()
    save("thermal.png")


def test_sun():
    s = rc.Sun()
    s.plot()
    save("sun.png")


def test_lamps():
    ax = rc.Sodium().plot()
    rc.Incandescent().plot(ax)
    rc.WhiteLED().plot(ax)
    plt.yscale("log")
    plt.ylim(1e-5, 1e2)
    save("lamps.png")


def test_phoenix_spectra():
    rc.Star().plot()
    rc.Star(teff=10000 * u.K, R=500).at(1 * u.au).plot()


if __name__ == "__main__":
    outputs = {k.split("_")[-1]: v() for k, v in locals().items() if "test_" in k}
