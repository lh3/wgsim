#include "tree.h"
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
    /*
	fprintf(stderr, "Options: -e FLOAT      base error rate [%.3f]\n", ERR_RATE);
	fprintf(stderr, "         -d INT        outer distance between the two ends [270]\n");
	fprintf(stderr, "         -s INT        standard deviation [50]\n");
	fprintf(stderr, "         -N INT        number of read pairs [1000000]\n");
	fprintf(stderr, "         -1 INT        length of the first read [150]\n");
	fprintf(stderr, "         -2 INT        length of the second read [150]\n");
	fprintf(stderr, "         -r FLOAT      rate of mutations [%.4f]\n", MUT_RATE);
	fprintf(stderr, "         -R FLOAT      fraction of indels [%.2f]\n", INDEL_FRAC);
	fprintf(stderr, "         -X FLOAT      probability an indel is extended [%.2f]\n", INDEL_EXTEND);
	fprintf(stderr, "         -S INT        seed for random generator [-1]\n");
	fprintf(stderr, "         -A FLOAT      discard if the fraction of ambiguous bases higher than FLOAT [%.2f]\n", MAX_N_RATIO);
	fprintf(stderr, "         -h            haplotype mode\n");
	fprintf(stderr, "         -f STR        the directory containing input FastA files [.]\n");
	fprintf(stderr, "         -o STR        the directory to write files output files to [.]\n");
	fprintf(stderr, "         -I            interlaved output\n");
	fprintf(stderr, "\n");
    */
	return 1;
}

int main(int argc, char *argv[])
{
    char *nwk_path = argv[1];
    if (argc < 2){
        return usage();
    }

    FILE *file = fopen(nwk_path, "r");
    int n_leaves, i, c;
    char *code;
    size_t n = 0;


    if (file == NULL) return 1; //could not open file
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

    double * curr_dist = (double *) malloc(sizeof(double)*(n_leaves));
    double * weights = (double *) malloc(sizeof(double)*((n_leaves*(n_leaves-1))/2));
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
    read_tree(&nwk_ptr, curr_dist, weights, names, id, n_leaves, 1);

    int j;
    for (i = 0; i < n_leaves; i++) {
        for (j = i+1; j < n_leaves; j++) {
        printf("%0.6f\n", weights[i*n_leaves + j - (i+1)*(i+2)/2]);
        }
    }

    hid_t dataset, datatype, dataspace, h5file;
    herr_t status;

    h5file = H5Fcreate("test.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    const hsize_t dimsf[] = {((n_leaves*(n_leaves-1))/2)};
    dataspace = H5Screate_simple(1, dimsf, NULL);

    datatype = H5Tcopy(H5T_IEEE_F64LE);
    status = H5Tset_order(datatype, H5T_ORDER_LE);

    /* Same as above, but use the property list */
    dataset = H5Dcreate1(h5file, "/distances", datatype, dataspace, H5P_DEFAULT);

    status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, weights);

    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Tclose(datatype);
    H5Fclose(h5file);

    free(curr_dist);
    free(weights);
    free(code);
    free(names);

}
