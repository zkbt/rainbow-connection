import rainbowconnection as trc
import matplotlib.pyplot as plt
import astropy.units as u

plt.ion()

def test_thermal():
    s = trc.Thermal(teff=4000*u.K, radius=1*u.m)
    s.plot()
    plt.show()


def test_sun():
    s = trc.Sun()
    s.plot()
    plt.show()

def test_lamps():
    ax = trc.Sodium().plot()
    trc.Incandescent().plot(ax)
    plt.yscale('log')
    plt.ylim(1e-5, 1e2)
    plt.show()

if __name__ == '__main__':
    test_sun()
    test_thermal()
    test_lamps()
