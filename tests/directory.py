import os
import matplotlib.pyplot as plt

directory = "examples"
try:
    os.mkdir(directory)
except:
    pass


def d(s):
    return os.path.join(directory, s)


def save(s):
    plt.tight_layout()
    plt.savefig(d(s), facecolor="black")
