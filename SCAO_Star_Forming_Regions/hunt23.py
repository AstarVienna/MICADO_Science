from path import Path
from matplotlib import pyplot as plt
from astropy.io import fits, ascii
from astropy.table import Table, Column
import numpy as np


CATALOGUE_FILE = Path(__file__).parent / "J_A+A_673_A114_clusters.dat.gz.fits"
MEMBERS_FILE = Path(__file__).parent / "reduced_J_A+A_673_A114_members.dat.fits"
README_FILE = Path(__file__).parent / "Hunt23_ReadMe.txt"

readme = ascii.read(README_FILE)


class Cluster():
    def __init__(self, cluster_id):
        self.cluster_id = cluster_id
        with fits.open(MEMBERS_FILE) as hdulist:
            self.hdu = hdulist[cluster_id]

            self.meta = dict(self.hdu.header)
            self.table = Table([Column(self.hdu.data[colname])
                                for colname in self.hdu.data.names],
                               names=self.hdu.data.names)

        self.ra = self.table["RAdeg"]
        self.dec = self.table["DEdeg"]
        self.G_mag = self.table["Gmag"]
        self.Rp_mag = self.table["Rpmag"]

        self.meta["GS_coverage"] = {"cluster_radius": "R50",
                                    "gs_radius": 25,
                                    "Rp_limit": 16,
                                    }

    @property
    def n_guide_stars(self):
        return sum(self.Rp_mag <= self.meta["GS_coverage"]["Rp_limit"])


    def calculate_fraction_covered_by_guide_stars(self):
        """
        Calculate the percentage of the cluster to a certain cluster radius
        covered by potential guide stars

        Parameters
        ----------
        cluster_radius : str
            Radius to calculate : ["RC", R50", "RT"]
        gs_radius : float, int
            [arcsec] Acceptable distance to guide star. Default=25"

        Returns
        -------
        frac_overlaps
            min(1, sum_overlaps)
        """

        ra_med, dec_med = self.meta["RADEG"], self.meta["DEDEG"]

        mask = self.Rp_mag < self.meta["GS_coverage"]["Rp_limit"]
        dras = (self.ra[mask] - ra_med) * 3600
        ddecs = (self.dec[mask] - dec_med) * 3600

        gs_radius = self.meta["GS_coverage"]["gs_radius"]
        radius_str = self.meta["GS_coverage"]["cluster_radius"]
        radius = self.meta[radius_str] * 3600

        if radius > 1024:
            f = 512 / radius
            radius *= f
            gs_radius *= f

        grid = add_circles_to_grid(int(2 * radius), dras + radius, ddecs + radius, gs_radius)
        grid_mask = add_circles_to_grid(int(2 * radius), [radius], [radius], radius)
        frac_overlaps = np.sum(grid * grid_mask) / (4 * radius**2 * 0.7854)

        return frac_overlaps

    def plot(self):
        ra_med, dec_med = self.meta["RADEG"], self.meta["DEDEG"]
        mask = self.Rp_mag < 16

        dras, ddecs = (self.ra - ra_med) * 3600, (self.dec - dec_med) * 3600

        plt.scatter(dras[np.invert(mask)], ddecs[np.invert(mask)], c="k")
        plt.scatter(dras[mask], ddecs[mask], c="r")

        gs_radius = self.meta["GS_coverage"]["gs_radius"]
        for dra, ddec in zip(dras[mask], ddecs[mask]):
            patch = plt.Circle((dra, ddec), gs_radius, fill=True, alpha=0.3, fc="r")
            plt.gca().add_artist(patch)

        rc = plt.Circle((0, 0), self.meta["RC"] * 3600, fill=False, ec="g", ls="--")
        r50 = plt.Circle((0, 0), self.meta["R50"] * 3600, fill=False, ec="b")
        rt = plt.Circle((0, 0), self.meta["RT"] * 3600, fill=False, ec="orange", ls=":")

        for patch in [rc, r50, rt]:
            plt.gca().add_artist(patch)

        gs_coverage = self.calculate_fraction_covered_by_guide_stars()
        plt.title(f'Coverage by Rp<{self.meta["GS_coverage"]["Rp_limit"]} '
                  f'Guide stars out to {self.meta["GS_coverage"]["cluster_radius"]}: '
                  f'{round(gs_coverage, 2)}')

        plt.colorbar()
        plt.show()

    def __getattr__(self, item):
        return getattr(self.table, item)

    def __getitem__(self, item):
        return self.table[item]


class Catalogue(object):
    def __init__(self):
        with fits.open(CATALOGUE_FILE) as hdulist:
            self.hdu = hdulist[1]

            self.meta = dict(self.hdu.header)
            self.table = Table([Column(self.hdu.data[colname])
                                for colname in self.hdu.data.names],
                               names=self.hdu.data.names)

    def __getattr__(self, item):
        return getattr(self.table, item)

    def __getitem__(self, item):
        return self.table[item]


def circles_fraction_overlap(circle1, circle2):
    """
    The fraction of circle1 that is covered by circle2

    Parameters
    ----------
    circle1, circle2 : tuple
        (x,y,r)

    Returns
    -------
    fraction_overlap : float

    """

    x1, y1, r1 = circle1
    x2, y2, r2 = circle2
    d = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

    if d == 0:
        return min(1, (r2 / r1)**2)
    elif d >= r1 + r2:
        return 0

    # x, y are the intercepts
    x = (r1 ** 2 - r2 ** 2 + d ** 2) / (2 * d)
    y = np.sqrt(r1 ** 2 - x ** 2)
    ang_1 = np.arctan2(y, x)
    ang_2 = np.arctan2(y, d-x)

    # The overlapping area is the sum of the D-shapes from each circle intercept
    # made by subtracting the triangles (x*y/2) from the chords (ang/[2*pi*r])
    overlap_area_1 = ang_1 * r1**2 - x * y
    overlap_area_2 = ang_2 * r2 ** 2 - (d - x) * y

    fraction_overlap = (overlap_area_1 + overlap_area_2) / (np.pi * r1**2)
    return fraction_overlap


def circle_mask(grid_size, cx, cy, radius):
    """
    Creates a 2D array where elements inside a given circle are set to 1, others remain 0.

    Parameters:
        grid_size (int, tuple): Dimensions of the 2D array (height, width).
        cx, cy (float): Center of the circle.
        radius (float): Radius of the circle.

    Returns:
        np.ndarray: A 2D array with ones inside the circle and zeros elsewhere.
    """
    if isinstance(grid_size, int):
        grid_size = (grid_size, grid_size)
    height, width = grid_size
    y, x = np.ogrid[:height, :width]

    mask = (x - cx) ** 2 + (y - cy) ** 2 <= radius ** 2

    return mask.astype(int)


def add_circles_to_grid(grid_size, x_positions, y_positions, radius):
    """
    Adds multiple circles to a grid based on given x and y positions.

    Parameters:
        grid_size (tuple): Dimensions of the 2D array (height, width).
        x_positions (list): List of x-coordinates for circle centers.
        y_positions (list): List of y-coordinates for circle centers.
        radius (float): Radius of all circles.

    Returns:
        np.ndarray: A 2D array with ones where circles are added and zeros elsewhere.
    """
    if isinstance(grid_size, int):
        grid_size = (grid_size, grid_size)
    grid = np.zeros(grid_size, dtype=int)
    for cx, cy in zip(x_positions, y_positions):
        grid += circle_mask(grid_size, cx, cy, radius)

    return np.clip(grid, 0, 1)


for i in np.linspace(1, 1766, 20):
    cl = Cluster(int(i))
    print(int(i), cl.n_guide_stars, cl.calculate_fraction_covered_by_guide_stars())
