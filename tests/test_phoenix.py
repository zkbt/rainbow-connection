from rainbowconnection.sources.phoenix import *


def test_download():
    print(
        """
    This might be a slightly annoylingly slow test
    for the first time that it is run, requiring a download
    of a 600MB zip file. However, the second time you
    run this test it should go a lot faster, since the
    file will have already been downloaded.
    """
    )
    return unzip_metallicity(0.0)


def test_phoenix_spectra():
    Star().plot()
    Star(teff=10000 * u.K, R=500).at(1 * u.au).plot()
    Star(teff=10000 * u.K, R=500, extend_wavelengths=True).at(1 * u.au).plot()
