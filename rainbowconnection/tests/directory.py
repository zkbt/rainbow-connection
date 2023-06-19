import os, warnings
import matplotlib.pyplot as plt

directory = "examples"
try:
    os.mkdir(directory)
except:
    pass


def d(s):
    return os.path.join(directory, s)


def save(s):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        plt.tight_layout()
        plt.savefig(d(s), facecolor="black")
