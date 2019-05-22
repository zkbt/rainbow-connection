import rainbowconnection as rc
import matplotlib.pyplot as plt
import astropy.units as u

plt.ion()

def test_earth():
    a = rc.Earth()
    a.plot()
    plt.show()


def test_hotjupiter():
    a = rc.HotJupiter()
    a.plot()
    plt.show()

def test_transmission():
    a = rc.Earth()
    s = rc.Sun()
    a.transmit(s).plot()
    plt.show()

if __name__ == '__main__':
    test_earth()
    test_hotjupiter()
    test_transmission()
