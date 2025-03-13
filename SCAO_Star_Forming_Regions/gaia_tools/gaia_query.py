import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia

# Query GAIA DR3 for objects within a suare region
#   ra:             RA [deg]
#   dec:            DEC [deg]
#   region_size:    square field size [arcsec]
def gaia_query(ra, dec, region_size):
    coord = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame='icrs')
    width = u.Quantity(region_size, u.arcsecond)
    height = u.Quantity(region_size, u.arcsecond)
    r = Gaia.query_object_async(coordinate=coord, width=width, height=height)
    return r

if __name__ == "__main__":
    result = gaia_query(49.85574296, 28.39919379, 50)
    # result.pprint(max_lines=12, max_width=130)
    for r in result:
        print(f"{r["source_id"]}\t{r["ra"]}\t{r["dec"]}\t{r["phot_rp_mean_mag"]}")
