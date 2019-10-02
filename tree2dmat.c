#include "tree.h"
#include "utils.h"
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <stdio.h>
#include <unistd.h>
#include <stdint.h>
#include <ctype.h>
#include <string.h>
#include "hdf5.h"

static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage: tree2dmat  [options] <in.nwk>\n\n");
	fprintf(stderr, "Options: -o STR        the HDF5 file to write the distance matrix to [distmat.h5]\n");
	fprintf(stderr, "         -d STR        the dataset to write the distance matrix to [distances]\n");
	fprintf(stderr, "         -n STR        the dataset to write the leaf names to [leaf_names]\n");
	fprintf(stderr, "         -v            verbose output\n");
	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char *argv[])
{
    FILE *file;
    int n_leaves, i, c, verbose = 0;
    char *code, *nwk_path,
         *h5_outpath = "distmat.h5",
         *h5_dset_path = "distances",
         *h5_names = "leaf_names";
    size_t n = 0;

	while ((c = getopt(argc, argv, "o:d:n:v")) >= 0) {
		switch (c) {
		case 'o': h5_outpath = optarg; break;
		case 'd': h5_dset_path = optarg; break;
		case 'n': h5_names = optarg; break;
		case 'v': verbose = 1; break;
        }
    }
    if (argc - optind < 1) return usage();

    nwk_path = argv[optind];
    if (verbose) logmsg("reading Newick string from %s", nwk_path);
    file = fopen(nwk_path, "r");

    if (file == NULL) {
        return 1; //could not open file
    }
    // compute file size
    fseek(file, 0, SEEK_END);
    long f_size = ftell(file);
    fseek(file, 0, SEEK_SET);
    code = (char *) malloc(f_size);

    // read in newick string
    int internodes_n = 0;
    while ((c = fgetc(file)) != EOF) {
        if (c == 10 || c == 32 || c == 9) {
            continue;
        }
        if (c == 44)
            internodes_n++;
        code[n++] = (char)c;
    }
    code[n] = '\0';

    n_leaves = internodes_n+1;
    if (verbose) logmsg("found %d leaves", n_leaves);

    double * curr_dist = (double *) malloc(sizeof(double)*(n_leaves));
    double * weights = (double *) malloc(sizeof(double)*((n_leaves*(n_leaves-1))/2));
    if (verbose) logmsg("computing %d weights", ((n_leaves*(n_leaves-1))/2));
    char ** names = (char **) malloc(sizeof(char*)*(n_leaves));
    for (i = 0; i < n_leaves; i++){
        curr_dist[i] = 0.0;
    }
    for (i = i+1; i < (n_leaves*(n_leaves-1)/2); i++){
        weights[i] = 0.0;
    }

    i = 0;
    int * id = &i;
    char * nwk_ptr = code;
    if (verbose) logmsg("reading branch lengths and computing distances");
    read_tree(&nwk_ptr, curr_dist, weights, names, id, n_leaves, 1);

    hid_t dataset, dataspace, h5file, dtype;
    hsize_t dimsf[] = { ((n_leaves*(n_leaves-1))/2) };

    if (verbose) logmsg("writing file to %s", h5_outpath);
    h5file = H5Fcreate (h5_outpath, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    dataspace = H5Screate_simple(1, dimsf, NULL);

    if (verbose) logmsg("writing dataset to %s", h5_dset_path);
    dataset = H5Dcreate (h5file, h5_dset_path, H5T_IEEE_F64LE, dataspace, H5P_DEFAULT,
                H5P_DEFAULT, H5P_DEFAULT);

    H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, weights);

    H5Dclose(dataset);
    H5Sclose(dataspace);

    dimsf[0] = n_leaves;

    dtype = H5Tcopy (H5T_C_S1);
    H5Tset_size (dtype, H5T_VARIABLE);

    /*
     * Create dataspace.  Setting maximum size to NULL sets the maximum
     * size to be the current size.
     */
    dataspace = H5Screate_simple (1, dimsf, NULL);

    /*
     * Create the dataset and write the variable-length string data to
     * it.
     */
    dataset = H5Dcreate (h5file, h5_names, dtype, dataspace, H5P_DEFAULT, H5P_DEFAULT,
                H5P_DEFAULT);
    H5Dwrite (dataset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, names);


    H5Fclose(h5file);

    free(curr_dist);
    free(weights);
    free(code);
    free(names);
    if (verbose) logmsg("done");
}
