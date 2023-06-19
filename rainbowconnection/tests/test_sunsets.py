import rainbowconnection as rc
import matplotlib.pyplot as plt
import astropy.units as u
from .directory import *


def test_sunset():
    s = rc.Sun().at(1 * u.au)
    e = rc.Earth()
    t = e.transmit(s)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t.animate_sunset(
            d("earth-sun-quick.mp4"),
            motionresolution=15 * u.deg,
        )


def test_everything():
    s = rc.Sun().at(1 * u.au)
    e = rc.Earth()
    t = e.transmit(s)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t.animate_everything(
            d("earth-sun-everything-quick.mp4"),
            motionresolution=15 * u.deg,
        )


if __name__ == "__main__":
    outputs = {k.split("_")[-1]: v() for k, v in locals().items() if "test_" in k}
