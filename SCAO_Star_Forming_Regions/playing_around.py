from matplotlib import pyplot as plt
from astropy.io import fits, ascii
import numpy as np

mamajek_table = ascii.read("EEM_dwarf_UBVIJHK_colors_Teff.txt")
Jmag = mamajek_table["M_J"]
Msun = mamajek_table["Msun"]

with fits.open("J_A+A_673_A114_clusters.dat.gz.fits") as hdulist:
    tbl = hdulist[1].data
    ra = tbl["RAdeg"]
    dec = tbl["DEdeg"]
    dist = tbl["dist50"]        # Unsure what the 16 and 84% estimates mean
    radius = tbl["rc"]          #

# Assumption is that clusters need to be at least 5 deg above
# ELT Alt limit of 60 deg from zenith
LAT_ARMAZONES = -25    # degrees South
dec_mask = dec < LAT_ARMAZONES + 55
dist_mask = dist < 18000

mask = dec_mask * dist_mask
print(f"Number of open clusters visible to MICADO: {sum(mask)}")

# plt.hist(dist[mask], bins=100)
# plt.semilogy()
# plt.show()


