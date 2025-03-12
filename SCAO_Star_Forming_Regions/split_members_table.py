from matplotlib import pyplot as plt
from astropy.io import fits, ascii
from astropy.table import Table, Column
import numpy as np


# Get the names of the clusters and their Hunt+23 IDs
hdulist_cluster = fits.open("J_A+A_673_A114_clusters.dat.gz.fits")
tbl = hdulist_cluster[1]
names = tbl.data["Name"]
hunt_ids = tbl.data["ID"]

# # split the members file into individual clusters
# with fits.open("J_A+A_673_A114_members.dat.fits") as hdulist:
#     col_names = ["ID", "GaiaDR3", "Prob",
#                  "RAdeg", "DEdeg", "pmRA", "pmDE",
#                  "Gmag", "Rpmag", "BP-G", "G-RP"]
#     cols = [Column(hdulist_cluster[1].data[col_name], name=col_name)
#             for col_name in col_names]
#     reduced_tbl = Table(cols, names=col_names)
#
#     new_hdulist = fits.HDUList([hdulist_cluster[0]])
#     for hunt_id, name in zip(hunt_ids, names):
#         mask = reduced_tbl["ID"] == hunt_id
#         new_hdulist.append(fits.table_to_hdu(reduced_tbl[mask]))
#         new_hdulist[-1].header["NAME"] = name
#         print(f"Made table for {hunt_id}: {name}. Length: {sum(mask)}")
#
#     new_hdulist.writeto("reduced_J_A+A_673_A114_members.dat.fits",
#                         overwrite=True)

clusters_tbl = Table([Column(tbl.data[colname]) for colname in tbl.data.names], names=tbl.data.names)
hdulist = fits.open("reduced_J_A+A_673_A114_members.dat.fits", mode="update")
for ii in hunt_ids:
    hdu = hdulist[ii+1]
    if not clusters_tbl[ii]["ID"] == ii:
        raise ValueError(f'{clusters_tbl[ii]["ID"]}  {ii}')
    hdu.header.update(dict(clusters_tbl[ii]))

hdulist.flush()


