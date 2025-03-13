from astropy.io import ascii
import numpy as np

cat_name = "Hunt23_clusters"

# Read the fixed-width table
# -> Returns an astropy.table.Table
def read_hunt(fname = cat_name):
    hdr = ascii.read(fname+".txt", format="fixed_width_two_line")
    names = []
    col_begin = []
    col_end = []
    for row in hdr:
        if type(row[0]) == np.str_:
            names.append(str(row[3]))
            col_range = row[0].split("-")
            cb = int(col_range[0])-1
            col_begin.append(cb)
            if len(col_range) > 1:
                col_end.append(int(col_range[1])-1)
            else:
                col_end.append(cb+1)

    # print(names)
    # print(col_begin)
    # print(col_end)

    data = ascii.read(fname+".dat", format="fixed_width_no_header",
                    col_starts=col_begin,
                    col_ends=col_end,
                    names=names)

    return data

if __name__ == "__main__":
    clusters = read_hunt()