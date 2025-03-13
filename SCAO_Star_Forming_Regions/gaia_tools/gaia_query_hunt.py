import read_hunt as catalogue
import gaia_query as gq
import time

SAVE_TO_FILE = True  # Save results to file?
MAGLIM = 16.0        # Limiting Rp mag to include
REGION_SIZE = 50     # Size of square region, N arcsec on a side
DATA_CHUNK_N = 0     # data chunk number within cluster catalogue
DATA_CHUNK_SZ = 1000 # number of clusters in chunk

# standalone script to query one 'chunk' of clusters from the Hunt23 catalogue
# from the GAIA archive and sace the list of guide stars to a text file
if __name__ == "__main__":
    clusters = catalogue.read_hunt()

    total = len(clusters)
    print(f"Querying GAIA DR3 for {DATA_CHUNK_SZ} cluster locations in catalogue...")

    header = "c_name\tc_ra[deg]\tc_dec[deg]\tsource_id\tdist\tra[deg]\tdec[deg]\tref_epoch\tpmra[mas/yr]\tpmdec[mas/yr]\tphot_rp_mean_mag\n"
    print(header)

    if SAVE_TO_FILE:
        ofile = open(f"results{DATA_CHUNK_N}.txt", "w")
        ofile.write(header)

    chunk_start = DATA_CHUNK_N*DATA_CHUNK_SZ
    chunk_end = chunk_start + DATA_CHUNK_SZ
    num_scao_clusters = 0
    for i, c in enumerate(clusters[chunk_start:chunk_end]):
        result = None
        retry = 0
        while result == None and retry < 3:
            try:
                result = gq.gaia_query(float(c["RAdeg"]), float(c["DEdeg"]), REGION_SIZE)
            except:
                print("***Query exception, retrying...***")
                time.sleep(1)
                retry = retry + 1

        c_name = c["Name"]

        lines = []
        for obj in result:
            if float(obj["phot_rp_mean_mag"]) <= MAGLIM:
                lines.append(f"{c_name}\t{c["RAdeg"]}\t{c["DEdeg"]}\tDR3_{obj["source_id"]}\t{obj["dist"]}\t{obj["ra"]}\t{obj["dec"]}\t{obj["ref_epoch"]}\t{obj["pmra"]}\t{obj["pmdec"]}\t{obj["phot_rp_mean_mag"]}")
        if len(lines) == 0:
            lines = [f"{c_name}\t{c["RAdeg"]}\t{c["DEdeg"]}\tDR3--\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0"]
        else:
            num_scao_clusters += 1

        print(f"{i}/{DATA_CHUNK_SZ}:")
        for l in lines:
            print(l)
            if SAVE_TO_FILE:
                ofile.write(l)
                ofile.write("\n")
        if SAVE_TO_FILE:
            ofile.flush()

    l = f"{num_scao_clusters = }\n"
    print(l)

