from rainbowconnection.resampletools import *


def test_grid(dx=1.0, N=500, snr=1000, seed=42, irregular=False):
    """
    a simple test to show how binning to a grid works
    """

    print(
        "Testing binning {} {} data binned to dx={}\n".format(
            N,
            {True: "irregularly gridded", False: "gridded"}[irregular],
            dx,
        )
    )

    # set up fake arrays
    np.random.seed(seed)
    if irregular:
        x = np.sort(np.random.uniform(0.5, 25, N))
    else:
        x = np.linspace(0.5, 25, N)
    unc = np.ones_like(x) / snr
    model = np.zeros_like(x) + np.exp(-0.5 * (x - 5) ** 2 / 2**2) * 0.01  # x*0.0002
    y = np.random.normal(model, unc, N)

    bx, by, bunc = bintogrid(x, y, unc, dx=dx)

    # plot the demonstrations
    plt.errorbar(
        x,
        y,
        unc,
        color="gray",
        linewidth=0,
        elinewidth=1,
        alpha=0.1,
    )
    plt.errorbar(bx, by, bunc, color="black", linewidth=0, elinewidth=3)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(
        "{} {} data binned to dx={}".format(
            N,
            {True: "irregularly gridded", False: "gridded"}[irregular],
            dx,
        )
    )


def test_R(R=10, N=500, snr=1000, seed=42, irregular=False):
    """
    a simple test to show how binning to constant R works
    """

    print(
        "Testing binning {} {} data binned to R={}\n".format(
            N,
            {True: "irregularly gridded", False: "gridded"}[irregular],
            R,
        )
    )

    # set up fake arrays
    np.random.seed(seed)
    if irregular:
        x = np.sort(np.random.uniform(0.5, 25, N))
    else:
        x = np.linspace(0.5, 25, N)
    unc = np.ones_like(x) / snr
    model = np.zeros_like(x) + np.exp(-0.5 * (x - 5) ** 2 / 2**2) * 0.01  # x*0.0002
    y = np.random.normal(model, unc, N)

    bx, by, bunc = bintoR(x, y, unc, R=R)

    # plot the demonstrations
    plt.errorbar(
        x,
        y,
        unc,
        color="gray",
        linewidth=0,
        elinewidth=1,
        alpha=0.1,
    )
    plt.errorbar(bx, by, bunc, color="black", linewidth=0, elinewidth=3)
    plt.xscale("log")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(
        "{} {} data binned to R={}".format(
            N,
            {True: "irregularly gridded", False: "gridded"}[irregular],
            R,
        )
    )


def test_FCR(supersample=True):
    """this function tests out the resampling code

    supersample=True
            means that there will be multiple new pixels per original pixel

    supersample=False
            means that there will be fewer new pixels than original pixels
    """

    print(
        "Testing flux-conserving resampling, when {}supersampling the original grid\n".format(
            {True: "", False: "not "}[supersample]
        )
    )
    xinitial = np.arange(39, 47)
    yinitial = np.random.uniform(0.0, 0.1, len(xinitial))
    if supersample:
        xresample = np.linspace(np.min(xinitial) - 1.5, np.max(xinitial) + 1.5, 50)
    else:
        xresample = np.linspace(np.min(xinitial) - 1.5, np.max(xinitial) + 1.5, 5)
    yresample = fluxconservingresample(xinitial, yinitial, xresample, visualize=True)


if __name__ == "__main__":
    testFCR()
    testFCR(supersample=False)
    testR(R=10)
    testR(R=10, irregular=True)
    testgrid(dx=0.5)
    testgrid(dx=0.5, irregular=True)
