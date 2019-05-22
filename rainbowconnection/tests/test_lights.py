import rainbowconnection as rc
import matplotlib.pyplot as plt
import astropy.units as u

plt.ion()

def test_thermal():
    s = rc.Thermal(teff=4000*u.K, radius=1*u.m)
    s.plot()
    plt.show()


def test_sun():
    s = rc.Sun()
    s.plot()
    plt.show()

def test_lamps():
    ax = rc.Sodium().plot()
    rc.Incandescent().plot(ax)
    plt.yscale('log')
    plt.ylim(1e-5, 1e2)
    plt.show()

if __name__ == '__main__':
    test_sun()
    test_thermal()
    test_lamps()
