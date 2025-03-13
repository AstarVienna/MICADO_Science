from astropy.io import ascii, fits
from astropy.table import Table, Column
import numpy as np

GAIA_TABLE = "./gaia_data.txt"

# Read the ascii table and convert to a multi-extension fits file
# each extension contains:
#   - header keywords with the cluster information
#   - data = fits table with list of GAIA objects
def read_gaia_table(fname):
    table = ascii.read(fname, format="tab")
    cols = table.colnames
    c_cols = cols[0:3]  # cluster info
    g_cols = cols[3:]   # guide star info

    # collapse table in to list of unique clusters
    hdict = {}
    for row in table:
        hdict[row["c_name"]] = 1

    # build hdulist - one extension per cluster
    hdulist = fits.HDUList()
    for name in list(hdict.keys()):
        # add list of GAIA stars
        mask = table["c_name"] == name
        hdulist.append(fits.table_to_hdu(table[g_cols][mask]))
        # add cluster info to header
        hdr = hdulist[-1].header
        hdr["EXTNAME"] = name
        for kw in c_cols[1:]:
            hdr[kw] = table[kw][mask][0]
    
    hdulist.writeto("J_A+A_673_A114_gaia_stars.dat.fits", overwrite=True)

if __name__ == "__main__":
    read_gaia_table(GAIA_TABLE)
