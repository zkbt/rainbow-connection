import rainbowconnection as rc
import matplotlib.pyplot as plt
import astropy.units as u
from directory import *
plt.ion()

def test_sunset():
    s = rc.Sun().at(1*u.au)
    e = rc.Earth()
    t = e.transmit(s)
    t.animate_sunset(d('earth-sun-quick.mp4'), motionresolution=10*u.deg)
    return t

if __name__ == '__main__':
    test_sunset()
