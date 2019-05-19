import therainbowconnection as trc
import matplotlib.pyplot as plt

def test_thermal():
    l = trc.Thermal()
    l.plot()
    plt.show()


def test_sun():
    l = trc.Sun()
    l.plot()
    plt.show()

if __name__ == '__main__':
    test_sun()
    #test_thermal()
