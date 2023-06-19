import rainbowconnection as rc
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from .directory import *

def test_plots():
    s = rc.Sun()
    s.plot()
    save("sun-plot.png")
    s.plot_as_rainbow()
    save("sun-plot-as-rainbow.png")
    s.plot_as_slit_rainbow()
    save("sun-plot-as-slit-rainbow.png")

    xlim = (300,1000)
    s = rc.Sun()
    s.plot()
    plt.xlim(*xlim)
    save("sun-plot-zoom.png")
    s.plot_as_rainbow()
    plt.xlim(*xlim)
    save("sun-plot-as-rainbow-zoom.png")
    s.plot_as_slit_rainbow()
    plt.xlim(*xlim)
    save("sun-plot-as-slit-rainbow-zoom.png")
