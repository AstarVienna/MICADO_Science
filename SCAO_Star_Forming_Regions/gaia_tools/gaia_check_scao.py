from astropy.io import ascii, fits
from astropy.table import Table, Column
import numpy as np

GAIA_FILE  = "./J_A+A_673_A114_gaia_stars.dat.fits.zip"
MAG_LIMIT  = 16.0   # faintest guide star to consider
DIST_LIMIT = 10.0     # arcsec radius from cluster to include

# check how many clisters do NOT have a suitable guide star
def get_no_scao(hdulist, mag_lim = 16.0, dist_lim = 50.0):
    no_scao = 0
    dist_lim /= 3600  # arcsec to deg
    for hdu in hdulist[1:]:
        # print(f"{hdu.header["EXTNAME"]}")
        if hdu.data['source_id'][0] == "DR3--" \
            or all(np.float64(list(hdu.data['phot_rp_mean_mag'])) > mag_lim) \
            or all(np.float64(list(hdu.data['dist'])) > dist_lim):
            no_scao += 1
    return no_scao


if __name__ == "__main__":
    with fits.open(GAIA_FILE) as hdulist:
        num_clusters = len(hdulist) - 1 # exclude primary header
        print(f"Searching {num_clusters} clusters in catalogue...")
        print("Search parameters:")
        print(f"\tGuide star brightness: Rp mag <= {MAG_LIMIT}")
        print(f"\tGuide star distance: <= {DIST_LIMIT} arcsec")

        bad_fields = get_no_scao(hdulist, MAG_LIMIT, DIST_LIMIT)
        num_good_fields = num_clusters - bad_fields

        print(f"{num_good_fields = }")
    
