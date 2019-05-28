import rainbowconnection as rc
import matplotlib.pyplot as plt
import astropy.units as u
from directory import *

plt.ion()

def test_earth():
    a = rc.Earth()
    a.plot()

    save('earth-atmosphere.pdf')


def test_hotjupiter():
    a = rc.HotJupiter()
    a.plot()
    plt.xlim(300*u.nm, 2000*u.nm)
    save('hotjupiter-atmosphere.pdf')

def test_transmission():
    a = rc.Earth()
    s = rc.Sun()
    a.transmit(s).plot()
    save('sun-through-earth-atmosphere.pdf')

if __name__ == '__main__':
    test_earth()
    test_hotjupiter()
    test_transmission()