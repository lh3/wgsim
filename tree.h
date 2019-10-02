#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "xrand.h"


typedef struct __tree_sampler {
    char ** names;
    double * cum_weights;
    int * node_ids;
    int rem_leaves;

} tree_sampler;


void read_tree(char ** nwk_ptr, double * curr_dist, double * weights, char ** names, int * id, int nleaves, int full_dmat) {
    char * nwk = *nwk_ptr;
    char *name = malloc(256), *blen = malloc(256);
    size_t name_n = 0, blen_n = 0;
    int s, m;   // s is starting index of nodes in the left tree, m is starting index of nodes in the right tree

    int read_blen = 0;
    char c = *nwk;
    s = *id;
    m = *id;
    while (c != '\0') {
        if (c == '(') { // read left node
            nwk++;
            read_tree(&nwk, curr_dist, weights, names, id, nleaves, full_dmat);
            m = *id;
        } else if (c == ')') { // done reading branch length
            blen[blen_n] = '\0';
            read_blen = 0;
            nwk++;
            break;
       } else if (c == ','){
            if (read_blen){ // add branch length to ret
                blen[blen_n] = '\0';
                break;
            } else {       // read right node
                nwk++;
                read_tree(&nwk, curr_dist, weights, names, id, nleaves, full_dmat);
            }
        } else if (c == ':') {  // read branch length
            if (name_n > 0){   // we were reading a leaf node
                name[name_n] = '\0';
                names[(*id)++] = name;
            }
            read_blen = 1;
            nwk++;
        } else {           // reading branch lenght or name
            if (read_blen) {
                blen[blen_n++] = c;
            } else {
                name[name_n++] = c;
            }
            nwk++;
        }
        c = *nwk;
    }
    int e = *id;    // the end index of nodes in the right tree
    assert( s <= m );
    assert( m < e );
    assert( m < nleaves );
    (*nwk_ptr) = nwk;
    int i, j;
    double tmp;
    if (full_dmat) {
        for (i = s; i < m; i++){
            for (j = m; j < e; j++){
                tmp = curr_dist[i] + curr_dist[j];
                weights[i * nleaves + j - ((i+1)*(i+2))/2] = tmp;
            }
        }
    } else {
        for (i = s; i < m; i++){
            for (j = m; j < e; j++){
                tmp = curr_dist[i] + curr_dist[j];
                weights[i] += tmp;
                weights[j] += tmp;
            }
        }
    }
    tmp = atof(blen);
    for (i = s; i < m; i++) {
        curr_dist[i] += tmp;
    }
    for (i = m; i < e; i++) {
        curr_dist[i] += tmp;
    }
    free(blen);
}

int binsearch(double * weights, int l, int r, double val){
    if (l == r)
        return l;
    int mid = l + (r-l)/2;
    if (val < weights[mid]){
        return binsearch(weights, l, mid, val);
    } else {
        return binsearch(weights, mid+1, r, val);
    }
}

void shift_weights(double * cum, int idx, int len) {
    int i;
    double w = -cum[idx];
    if (idx > 0)
        w += cum[idx-1]; // weight of idx
    if (len > 1)
        w += cum[len-1] - cum[len-2]; // weight of last element
    for (i = idx; i < len-1; i++)
        cum[i] = cum[i] + w;
    cum[len-1] = -1.0;
}

char ** sample_tree(char * nwk_path, int nsamples) {

    FILE *file = fopen(nwk_path, "r");
    char *code;
    size_t n = 0;
    int c;

    if (file == NULL) return NULL; //could not open file
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

    // allocate distance and weights
    int nleaves = internodes_n+1;
    if (nsamples > nleaves){
        printf("nsamples must be fewer than the number of leaves (%d)\n", nleaves);
        return NULL;
    }
    double * curr_dist = (double *) malloc(sizeof(double)*(nleaves));
    double * weights = (double *) malloc(sizeof(double)*(nleaves));
    char ** names = (char **) malloc(sizeof(char*)*(nleaves));
    int i;
    for (i = 0; i <= internodes_n; i++){
        curr_dist[i] = 0.0;
        weights[i] = 0.0;
    }

    // parse newick string and compute weights
    i = 0;
    int * id = &i;
    char * nwk_ptr = code;
    read_tree(&nwk_ptr, curr_dist, weights, names, id, nleaves, 0);

    // compute cumulative weights
    double * cum_weights = (double*) malloc(sizeof(double)*nleaves);
    double sum = 0.0;
    for (i = 0; i < nleaves; i++){
        sum += weights[i];
        cum_weights[i] = sum;
    }

    // allocate node ids to keep track of samples
    int * node_ids = (int*) malloc(nleaves*sizeof(int));
    for (i = 0; i < nleaves; i++){
        node_ids[i] = i;
    }

    // sample leaves
    int idx;
    int rem_leaves = nleaves;
    double U;
    int tmp;
    double r;
    for (i = 0; i < nsamples; i++){
        r = xdrand();
        U = cum_weights[rem_leaves-1]*r;  // draw random number
        idx = binsearch(cum_weights, 0, rem_leaves, U); // get index
        shift_weights(cum_weights, idx, rem_leaves);    // update weights
        tmp = node_ids[idx];
        node_ids[idx] = node_ids[rem_leaves-1];
        node_ids[rem_leaves-1] = tmp;
        rem_leaves--;
    }

    tmp = rem_leaves;
    for (i = 0; i < tmp; i++){
        free(names[node_ids[i]]);
    }
    char ** ret = (char **) malloc(nsamples * sizeof(char *));
    for (i = tmp; i < nleaves; i++){
        ret[i-tmp] = names[node_ids[i]];
    }
    free(curr_dist);
    free(weights);
    free(cum_weights);
    free(node_ids);
    free(code);
    free(names);
    return ret;
}

tree_sampler * sampler_init(char * nwk_path){
    FILE *file = fopen(nwk_path, "r");
    char *code;
    size_t n = 0;
    int c;

    if (file == NULL) return NULL; //could not open file
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

    // allocate distance and weights
    int nleaves = internodes_n+1;
    double * curr_dist = (double *) malloc(sizeof(double)*(nleaves));
    double * weights = (double *) malloc(sizeof(double)*(nleaves));
    char ** names = (char **) malloc(sizeof(char*)*(nleaves));
    int i;
    for (i = 0; i <= internodes_n; i++){
        curr_dist[i] = 0.0;
        weights[i] = 0.0;
    }

    // parse newick string and compute weights
    i = 0;
    int * id = &i;
    char * nwk_ptr = code;
    read_tree(&nwk_ptr, curr_dist, weights, names, id, nleaves, 0);
    free(curr_dist);
    free(code);

    // compute cumulative weights
    double * cum_weights = (double*) malloc(sizeof(double)*nleaves);
    double sum = 0.0;
    for (i = 0; i < nleaves; i++){
        sum += weights[i];
        cum_weights[i] = sum;
    }
    free(weights);

    // allocate node ids to keep track of samples
    int * node_ids = (int*) malloc(nleaves*sizeof(int));
    for (i = 0; i < nleaves; i++){
        node_ids[i] = i;
    }

    tree_sampler * ret = (tree_sampler *) malloc(sizeof(tree_sampler));
    ret->cum_weights = cum_weights;
    ret->node_ids = node_ids;
    ret->names = names;
    ret->rem_leaves = nleaves;;

    return ret;
}

int next_leaf(tree_sampler * ts){
    double * cum_weights = ts->cum_weights;
    int rem_leaves = ts->rem_leaves;
    int * node_ids = ts->node_ids;
    double U = cum_weights[rem_leaves-1]*xdrand();       // draw random number
    int idx = binsearch(cum_weights, 0, rem_leaves, U);  // get index
    shift_weights(cum_weights, idx, rem_leaves);         // update weights
    int tmp = node_ids[idx];
    node_ids[idx] = node_ids[rem_leaves-1];
    node_ids[rem_leaves-1] = tmp;
    ts->rem_leaves--;
    return tmp;
}



double * calc_abund(int nsamples) {
    double * abundance = (double *) malloc(sizeof(double)*nsamples);
    double sum = 0.0;
    int i;
    for (i = 0; i < nsamples; i++){
        abundance[i] = -log(-xdrand() + 1);
        sum += abundance[i];
    }
    for (i = 0; i < nsamples; i++)
        abundance[i] = abundance[i]/sum;
    return abundance;
}

