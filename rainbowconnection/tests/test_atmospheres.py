import rainbowconnection as rc
import matplotlib.pyplot as plt
import astropy.units as u
from .directory import *


def test_earth():
    a = rc.Earth()
    a.plot()
    save("earth-atmosphere.png")


def test_hotjupiter():
    a = rc.HotJupiter()
    a.plot()
    plt.xlim(300 * u.nm, 2000 * u.nm)
    save("hotjupiter-atmosphere.png")


def test_transmission():
    a = rc.Earth()
    s = rc.Sun()
    a.transmit(s).plot()
    save("sun-through-earth-atmosphere.png")


if __name__ == "__main__":
    outputs = {k.split("_")[-1]: v() for k, v in locals().items() if "test_" in k}
