"""
Set up what's needed to access a directory full
of PHOENIX model spectra.
"""

__all__ = [
    "get_metallicity_directory",
    "unzip_metallicity",
    "setup_metallicities",
    "read_exact_phoenix",
    "read_phoenix",
]

from ...imports import *
from .utils import *
from ...resampletools import *
from astropy.io import fits
from ..thermal import Thermal

preloaded = {}

# specify where to get the data from online
base_url = "http://phoenix.astro.physik.uni-goettingen.de/data/"
online_library_directory = "MedResFITS/R10000FITS/"

# FIXME -- make the local storage directory more flexible?
base_directory = os.path.join(os.getenv("HOME"), ".rainbow-connection")
directory_template = "PHOENIX-ACES-AGSS-COND-2011_R10000FITS_Z{metallicity}"

# define a more flexible template string to catch all files
file_template = (
    "lte{Teff:05.0f}-{logg:04.2f}{metallicity}.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits"
)
flexible_file_template = "lte*-*{metallicity}.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits"


def get_metallicity_directory(Z=0.0):
    """
    Figure out model directory for a particular metallicity.

    Parameters
    ----------
    Z : float
        [Fe/H]-style metallicity (= 0.0 for solar)
    """
    return directory_template.format(metallicity=stringify_metallicity(Z))


def get_metallicity_url(Z=0.0):
    """
    Figure out the URL to download a particular
    metallicity subset of models.

    Parameters
    ----------
    Z : float
        [Fe/H]-style metallicity (= 0.0 for solar)
    """
    model_directory = get_metallicity_directory(Z)
    url = f"{base_url}{online_library_directory}{model_directory}.zip"
    return url


def download_metallicity(Z=0.0):
    """
    Make sure the zip file of a particular metallicity has
    been downloaded.

    Parameters
    ----------
    Z : float
        [Fe/H]-style metallicity (= 0.0 for solar)
    """

    # what's the URL of the file to download
    url = get_metallicity_url(Z)

    print(
        f"""
    Checking for a cached version of
    {url}
    or downloading it for the first time.
    """
    )
    filename = download_file(url, cache=True, show_progress=True, pkgname=package_name)

    return filename


def unzip_metallicity(Z=0.0):
    """
    Unzip the downloaded PHOENIX model zip file.
    (Figure out the URL from the metallicity,
    make sure astropy has downloaded it into the
    cache directory, and then unzip into a more
    friendly directory structure.)

    Parameters
    ----------
    Z : float
        [Fe/H]-style metallicity (= 0.0 for solar)
    """

    # figure out the local
    url = get_metallicity_url(Z)
    local_path_to_zip_file = download_metallicity(Z)

    # find the base directory (something like /Users/zkbt/.rainbow-connection/)
    i_base = local_path_to_zip_file.find("cache/download")
    assert i_base > 0
    local_base_directory = local_path_to_zip_file[:i_base]
    local_directory = os.path.join(local_base_directory, get_metallicity_directory(Z))

    # check if unzipped files already exist
    pattern = os.path.join(local_directory, "*.fits")
    unzipped_files = glob.glob(pattern)
    N = len(unzipped_files)
    print(
        f"""
    Checking the file pattern
    {pattern}
    for existing models that have already been extracted.
    """
    )
    if N > 0:
        print(
            f"""
    The file pattern
    {pattern}
    already matches {N} files.
    If you want to re-extract these files,
    please delete the existing ones.
        """
        )
    else:
        # update what's happening
        print(
            f"""
        The file downloaded from
        {url}
        was stored locally as
        {local_path_to_zip_file}
        and is now being unzipped to
        {local_directory}
        It may take a while...
        """
        )
        with zipfile.ZipFile(local_path_to_zip_file, "r") as zip_ref:
            zip_ref.extractall(local_directory)

    return local_directory


def setup_metallicities(Zs=[-0.5, 0.0, 0.5]):
    """
    Download and unzip a set of metallicities.

    Parameters
    ----------
    Z : list, array
        [Fe/H]-style metallicities (= 0.0 for solar)
    """
    for Z in Zs:
        unzip_metallicity(Z)


def get_downloaded_models():
    """
    Identify what models have already been downloaded,
    in terms of metallicity, Teff, logg.

    Returns
    -------
    available : dict
        A dictionary with (float) metallicities as keys.
        Each value is a pair of arrays (Teff, logg) to
        indicate the values of effective temperature and
        log(surface gravity) that are available for that
        surface gravity value.
    """

    # figure out what metallicity directories are avaialble
    available_metallicity_directories = glob.glob(
        os.path.join(
            base_directory,
            directory_template.format(metallicity="*"),
        )
    )

    available_metallicities = {
        float(d.split("/")[-1].split("_Z")[-1]): d
        for d in available_metallicity_directories
    }

    availability_by_metallicity = {}

    if len(available_metallicities) == 0:
        print(
            """
        It is totally fine to proceed using `rainbow-connection` as is, without
        downloading any addition files. However, if you want to have access
        access to stellar spectra of arbitrary temperatures, you will need
        to download some PHOENIX model spectra to your computer. It appears
        that no spectra have yet been downloaded.

        To download a set of model spectra for a particular metallicity, run
        ```
        from rainbowconnection.sources.phoenix import unzip_metallicity
        unzip_metallicity(Z=0.0)
        ```
        replacing the `Z=` argument with whatever (log-solar) metallicity
        you want to download. Available options are in increments of 0.5 dex.

        This will download a very large file for each metallicty, so it may
        take a fairly long time. Please be patient. It will download and
        unzip these files into a (hidden) directory called `.rainbow-connection`
        in your $HOME directory.
        """
        )
    else:
        print(
            "The following PHOENIX stellar spectra have been downloaded\n"
            "and are available to use via `Star()` source objects:\n"
        )
        for metallicity in available_metallicities:
            Z = stringify_metallicity(metallicity)
            this_metallicity = os.path.join(
                base_directory,
                directory_template,
                flexible_file_template,
            ).format(metallicity=stringify_metallicity(metallicity))

            # for now, just interpolate in temperature
            Teffs, loggs = [], []
            for t in np.sort(glob.glob(this_metallicity)):
                x = t.split("/lte")[-1].split("-")[0:2]
                Teffs.append(float(x[0]))
                loggs.append(float(x[1].split("+")[0]))
            Teffs = np.array(Teffs)
            loggs = np.array(loggs)

            availability_by_metallicity[metallicity] = (
                Teffs,
                loggs,
            )
            print(
                f"Metallicity [Fe/H] = {metallicity} has {len(Teffs)} spectra, stored in..."
            )
            print(f" {this_metallicity}")

    return availability_by_metallicity


availability_by_metallicity = get_downloaded_models()


def read_exact_phoenix(Teff=5800, logg=4.5, metallicity=0.0, photons=True, R=None, extend_wavelengths=False):

    key = (Teff, logg, metallicity, photons, R)

    try:
        wave, flux = preloaded[key]
    except KeyError:
        filename = os.path.join(
            base_directory, directory_template, file_template
        ).format(
            Teff=Teff,
            logg=logg,
            metallicity=stringify_metallicity(metallicity),
        )

        hdus = fits.open(filename)
        # flux units are erg/s/cm^2/cm
        flux_cgs = hdus[0].data
        # flux units are erg/s/cm^2/nm
        flux_nm = flux_cgs / 1e7
        header = hdus[0].header

        # wavelength units are nm
        wave = np.exp(header["CRVAL1"] + header["CDELT1"] * np.arange(header["NAXIS1"])) / 10

        if extend_wavelengths:
            # include kludge to extend to wavelengths beyond 2500nm
            t = Thermal(Teff*u.K, radius=1*u.Rsun).at(1*u.Rsun)
            extra_wave = np.linspace(np.max(wave), 6000)*u.nm
            extra_flux = t.surface_flux(extra_wave).to('erg/(s * cm**2 * nm)')

            wave = np.hstack([wave, extra_wave.value])
            flux_nm = np.hstack([flux_nm, extra_flux.value])

        # energy = (con.h*con.c/w/u.ph).to('J/ph')
        # f = (t.surface_flux(w)/energy).to('ph/(s * cm**2 * nm)')
        """
        k_B = 1.38e-16       # erg/K
        if extend_wavelengths:

            # define extra wavelength grid
            original_R = 1/header["CDELT1"]
            w_extra = np.max(wave)*np.exp(np.arange(int(original_R))/original_R)
            w_extra_in_cm = w_extra/1e7

            # calculate the thing in the
            up = h * c / (w_extra_in_cm * k_B * Teff)

            # calculate the intensity from the Planck function
            flux_extra = np.pi * (2 * h * c ** 2 / w_extra_in_cm ** 5 / (np.exp(up) - 1))

            wave = np.hstack([wave, w_extra])
            flux_nm = np.hstack([flux_nm, flux_extra])
            print('extending wavelengths')"""

        if photons:
            h = 6.6260755e-27	 # erg s
            c = 2.99792458e10    # cm/s
            w_cm = wave / 1e7
            photon_energy = h * c / w_cm
            flux = flux_nm / photon_energy
        else:
            flux = flux_nm

        if R is not None:
            wave, flux = bintoR(wave, flux, R=R)
        preloaded[key] = wave, flux

    return wave, flux


def read_phoenix(Teff=5800, logg=4.5, metallicity=0.0, photons=True, R=None, extend_wavelengths=False):
    try:
        w, f = read_exact_phoenix(
            Teff=Teff,
            logg=logg,
            metallicity=metallicity,
            photons=photons,
            R=R,
            extend_wavelengths=extend_wavelengths
        )
    except FileNotFoundError:
        # print('The exact parameters {} were not found in {}'.format(locals(), phoenix_directory))

        Teffs, loggs = availability_by_metallicity[metallicity]
        # figure out the closest gravity
        best_gravity = find_nearest(loggs, logg)
        """if best_gravity != logg:
            print(
                "nudged logg from {} to {}".format(
                    logg, best_gravity
                )
            )"""

        Teffs = Teffs[loggs == best_gravity]
        close_temperatures = find_two_nearest(Teffs, Teff, verbose=False)
        weights = interpolation_weights(close_temperatures, Teff)

        w_one, f_one = read_exact_phoenix(
            Teff=close_temperatures[0],
            logg=best_gravity,
            metallicity=metallicity,
            photons=photons,
            R=R,
            extend_wavelengths=extend_wavelengths
        )
        w_other, f_other = read_exact_phoenix(
            Teff=close_temperatures[1],
            logg=best_gravity,
            metallicity=metallicity,
            photons=photons,
            R=R,
            extend_wavelengths=extend_wavelengths
        )
        assert (w_one == w_other).all()
        w = w_one
        # do logarithmic interpolation
        f = np.exp(weights[0] * np.log(f_one) + weights[1] * np.log(f_other))

    return w, f
