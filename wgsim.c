/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).
                 2011 Heng Li <lh3@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* This program is separated from maq's read simulator with Colin
 * Hercus' modification to allow longer indels. */


#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <stdio.h>
#include <unistd.h>
#include <stdint.h>
#include <ctype.h>
#include <string.h>
#include <zlib.h>
#include "kseq.h"
#include "tree.h"
#include "xrand.h"
#include "utils.h"
KSEQ_INIT(gzFile, gzread)

#define PACKAGE_VERSION "0.1"

const uint8_t nst_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

/* Simple normal random number generator, copied from genran.c */

double ran_normal()
{
	static int iset = 0;
	static double gset;
	double fac, rsq, v1, v2;
	if (iset == 0) {
		do {
			v1 = 2.0 * xdrand() - 1.0; // drand48
			v2 = 2.0 * xdrand() - 1.0; // drand48
			rsq = v1 * v1 + v2 * v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac = sqrt(-2.0 * log(rsq) / rsq);
		gset = v1 * fac;
		iset = 1;
		return v2 * fac;
	} else {
		iset = 0;
		return gset;
	}
}

/* wgsim */

enum muttype_t {NOCHANGE = 0, INSERT = 0x1000, SUBSTITUTE = 0xe000, DELETE = 0xf000};
typedef unsigned short mut_t;
static mut_t mutmsk = (mut_t)0xf000;

typedef struct {
	int l, m; /* length and maximum buffer size */
	mut_t *s; /* sequence */
} mutseq_t;

typedef struct {
    char * id;
    const char * fn; /* the fasta file path */
    uint64_t tot_len;
    uint64_t wlen;
    int n_ref;
} ref_t;

static double ERR_RATE = 0.02;
static double MUT_RATE = 0.001;
static double INDEL_FRAC = 0.15;
static double INDEL_EXTEND = 0.3;
static double MAX_N_RATIO = 0.05;

void wgsim_mut_diref(const kseq_t *ks, int is_hap, mutseq_t *hap1, mutseq_t *hap2)
{
	int i, deleting = 0;
	mutseq_t *ret[2];

	ret[0] = hap1; ret[1] = hap2;
	ret[0]->l = ks->seq.l; ret[1]->l = ks->seq.l;
	ret[0]->m = ks->seq.m; ret[1]->m = ks->seq.m;
	ret[0]->s = (mut_t *)calloc(ks->seq.m, sizeof(mut_t));
	ret[1]->s = (mut_t *)calloc(ks->seq.m, sizeof(mut_t));
	for (i = 0; i != ks->seq.l; ++i) {
		int c;
		c = ret[0]->s[i] = ret[1]->s[i] = (mut_t)nst_nt4_table[(int)ks->seq.s[i]];
        if (deleting) {
            if (xdrand() < INDEL_EXTEND) { // drand48
                if (deleting & 1) ret[0]->s[i] |= DELETE;
                if (deleting & 2) ret[1]->s[i] |= DELETE;
                continue;
            } else deleting = 0;
        }
		if (c < 4 && xdrand() < MUT_RATE) { // mutation
			if (xdrand() >= INDEL_FRAC) { // substitution // drand48
				double r = xdrand(); // drand48
				c = (c + (int)(r * 3.0 + 1)) & 3;
				if (is_hap || xdrand() < 0.333333) { // hom // drand48
					ret[0]->s[i] = ret[1]->s[i] = SUBSTITUTE|c;
				} else { // het
					ret[xdrand()<0.5?0:1]->s[i] = SUBSTITUTE|c; // drand48
				}
			} else { // indel
				if (xdrand() < 0.5) { // deletion // drand48
					if (is_hap || xdrand() < 0.333333) { // hom-del // drand48
						ret[0]->s[i] = ret[1]->s[i] = DELETE;
                        deleting = 3;
					} else { // het-del
                        deleting = xdrand()<0.5?1:2; // drand48
						ret[deleting-1]->s[i] = DELETE;
					}
				} else { // insertion
                    int num_ins = 0, ins = 0;
                    do {
                        num_ins++;
                        ins = (ins << 2) | (int)(xdrand() * 4.0); // drand48
                    } while (num_ins < 4 && xdrand() < INDEL_EXTEND); // drand48

					if (is_hap || xdrand() < 0.333333) { // hom-ins // drand48
						ret[0]->s[i] = ret[1]->s[i] = (num_ins << 12) | (ins << 4) | c;
					} else { // het-ins
						ret[xdrand()<0.5?0:1]->s[i] = (num_ins << 12) | (ins << 4) | c; // drand48
					}
				}
			}
		}
	}
}
void wgsim_print_mutref(FILE * out, char* id, const char *name, const kseq_t *ks, mutseq_t *hap1, mutseq_t *hap2)
{
	int i, j = 0; // j keeps the end of the last deletion
	for (i = 0; i != ks->seq.l; ++i) {
		int c[3];
		c[0] = nst_nt4_table[(int)ks->seq.s[i]];
		c[1] = hap1->s[i]; c[2] = hap2->s[i];
		if (c[0] >= 4) continue;
		if ((c[1] & mutmsk) != NOCHANGE || (c[2] & mutmsk) != NOCHANGE) {
			if (c[1] == c[2]) { // hom
				if ((c[1]&mutmsk) == SUBSTITUTE) { // substitution
					fprintf(out, "%s\t%s\t%d\t%c\t%c\t-\n", id, name, i+1, "ACGTN"[c[0]], "ACGTN"[c[1]&0xf]);
				} else if ((c[1]&mutmsk) == DELETE) { // del
					if (i >= j) {
						fprintf(out, "%s\t%s\t%d\t", id, name, i+1);
						for (j = i; j < ks->seq.l && hap1->s[j] == hap2->s[j] && (hap1->s[j]&mutmsk) == DELETE; ++j)
							fputc("ACGTN"[nst_nt4_table[(int)ks->seq.s[j]]], out);
						fprintf(out, "\t-\t-\n");
					}
				} else if (((c[1] & mutmsk) >> 12) <= 4) { // ins
					fprintf(out, "%s\t%s\t%d\t-\t", id, name, i+1);
                    int n = (c[1]&mutmsk) >> 12, ins = c[1] >> 4;
                    while (n > 0) {
                        fputc("ACGTN"[ins & 0x3], out);
						ins >>= 2;
                        n--;
                    }
                    fprintf(out, "\t-\n");
				} // else: deleted base in a long deletion
			} else { // het
				if ((c[1]&mutmsk) == SUBSTITUTE || (c[2]&mutmsk) == SUBSTITUTE) { // substitution
					fprintf(out, "%s\t%s\t%d\t%c\t%c\t+\n", id, name, i+1, "ACGTN"[c[0]], "XACMGRSVTWYHKDBN"[1<<(c[1]&0x3)|1<<(c[2]&0x3)]);
				} else if ((c[1]&mutmsk) == DELETE) {
					if (i >= j) {
						fprintf(out, "%s\t%s\t%d\t", id, name, i+1);
						for (j = i; j < ks->seq.l && hap1->s[j] != hap2->s[j] && (hap1->s[j]&mutmsk) == DELETE; ++j)
							fputc("ACGTN"[nst_nt4_table[(int)ks->seq.s[j]]], out);
						fprintf(out, "\t-\t-\n");
					}
				} else if ((c[2]&mutmsk) == DELETE) {
					if (i >= j) {
						fprintf(out, "%s\t%s\t%d\t", id, name, i+1);
						for (j = i; j < ks->seq.l && hap1->s[j] != hap2->s[j] && (hap2->s[j]&mutmsk) == DELETE; ++j)
							fputc("ACGTN"[nst_nt4_table[(int)ks->seq.s[j]]], out);
						fprintf(out, "\t-\t-\n");
					}
				} else if (((c[1] & mutmsk) >> 12) <= 4 && ((c[1] & mutmsk) >> 12) > 0) { // ins1
					fprintf(out, "%s\t%s\t%d\t-\t", id, name, i+1);
                    int n = (c[1]&mutmsk) >> 12, ins = c[1] >> 4;
                    while (n > 0) {
                        fputc("ACGTN"[ins & 0x3], out);
						ins >>= 2;
                        n--;
                    }
                    fprintf(out, "\t+\n");
				} else if (((c[2] & mutmsk) >> 12) <= 4 || ((c[2] & mutmsk) >> 12) > 0) { // ins2
					fprintf(out, "%s\t%s\t%d\t-\t", id, name, i+1);
                    int n = (c[2]&mutmsk) >> 12, ins = c[2] >> 4;
                    while (n > 0) {
                        fputc("ACGTN"[ins & 0x3], out);
                        ins >>= 2;
                        n--;
                    }
                    fprintf(out, "\t+\n");
				} // else: deleted base in a long deletion
			}
		}
	}
}

double perc_nuc(kseq_t * ks){
    uint64_t len = strlen(ks->seq.s);
    char * s = ks->seq.s;
    uint64_t counts[128] = {
	    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
    };
    uint64_t i;
    for (i = 0; i < len; i++) {
        counts[(int) *(s++)]++;
    }
    int idx[] = {65, 67, 71, 84};
    double n_nucs = 0;
    for (i = 0; i < 4; i++){
        n_nucs += counts[idx[i]] + counts[idx[i]+32];
    }
    return n_nucs/len;
}

int check_seq(kseq_t *ks, double min_amb, int min_len){
    int ret = 0;
    int l = ks->seq.l;
    double pn;
    pn = perc_nuc(ks);
    if (pn < min_amb){
        ret += 1;
    }
	if (l < min_len) {
        ret += 2;
	}
    return ret;
}

void read_seq(const char *fn, ref_t * ref, int min_len){
	gzFile fp_fa;
    int l, n_ref;
    uint64_t tot_len;
    // ks->seq.s should be the sequence
	kseq_t *ks;

	fp_fa = gzopen(fn, "r");
	ks = kseq_init(fp_fa);
	tot_len = n_ref = 0;
    logmsg("[%s] calculating the total length of the reference sequence %s...", __func__, fn);
	while ((l = kseq_read(ks)) >= 0) {
        if (check_seq(ks, 0.999, min_len)) continue;
		tot_len += l;
		++n_ref;
	}

    ref->fn = fn;
    ref->tot_len = tot_len;
    ref->n_ref = n_ref;

	logmsg("[%s] %d useable sequences in %s, total length: %llu", __func__, n_ref, fn, (long long)tot_len);
	kseq_destroy(ks);
	gzclose(fp_fa);
}

uint64_t wgsim_core(FILE *fpout1, FILE *fpout2, FILE *mutout, ref_t * ref, int is_hap, uint64_t N, int dist, int std_dev, int size_l, int size_r)
{
	kseq_t *ks;
    mutseq_t rseq[2];
	gzFile fp_fa;
	uint64_t ii;
	int i, l;
    uint64_t actual_n_pairs = 0;
	char *qstr;
	int /*size[2],*/ Q, max_size;
	uint8_t *tmp_seq[2];
    mut_t *target;

	l = size_l > size_r? size_l : size_r;
	qstr = (char*)calloc(l+1, 1);
	tmp_seq[0] = (uint8_t*)calloc(l+2, 1);
	tmp_seq[1] = (uint8_t*)calloc(l+2, 1);
	//size[0] = size_l; size[1] = size_r;
	max_size = size_l > size_r? size_l : size_r;

    // make quality scores high
	// Q = (ERR_RATE == 0.0)? 'I' : (int)(-10.0 * log(ERR_RATE) / log(10.0) + 0.499) + 33;
    Q = 'I';


	fp_fa = gzopen(ref->fn, "r");
	ks = kseq_init(fp_fa);
    int ok_bit;
    int min_len = dist + 3 * std_dev;
    double min_nuc = 0.999;
	while ((l = kseq_read(ks)) >= 0) {
        ok_bit = check_seq(ks, min_nuc, min_len);
        if (ok_bit) {
            if (ok_bit % 2) logmsg("[%s] skip sequence '%s' (length = %d) from %s -- fraction of noniambiguous nucleotides less than %0.5f", __func__, ks->name.s, l, ref->id, min_nuc);
            if (ok_bit / 2) logmsg("[%s] skip sequence '%s' from %s -- sequence is shorter than %d bases", __func__, ks->name.s, ref->id, min_len);
            continue;
        }
		uint64_t n_pairs = (uint64_t)((long double)l / ref->tot_len * N + 0.5);
		// generate mutations and print them out
		wgsim_mut_diref(ks, is_hap, rseq, rseq+1);
		wgsim_print_mutref(mutout, ref->id, ks->name.s, ks, rseq, rseq+1);

	    FILE *fpo[] = { fpout1, fpout2 };
        int s[] = { size_l, size_r };
		for (ii = 0; ii != n_pairs; ++ii) { // the core loop
			double ran;
			int d, pos, /*s[2],*/ order[2], is_flip = 0;
			int n_sub[2], n_indel[2], n_err[2], ext_coor[2], j, k;
			//FILE *fpo[2];

			do { // avoid boundary failure
				ran = ran_normal();
				ran = ran * std_dev + dist;
				d = (int)(ran + 0.5);
				d = d > max_size? d : max_size;
				pos = (int)((l - d + 1) * xdrand()); // drand48
			} while (pos < 0 || pos >= ks->seq.l || pos + d - 1 >= ks->seq.l);

			// flip or not
			if (xdrand() < 0.5) { // drand48
                order[0] = 0;
                order[1] = 1;
			} else {
                order[0] = 1;
                order[1] = 0;
				is_flip = 1;
			}

			// generate the read sequences
			target = rseq[xdrand()<0.5?0:1].s; // haplotype from which the reads are generated // drand48
			n_sub[0] = n_sub[1] = n_indel[0] = n_indel[1] = n_err[0] = n_err[1] = 0;

#define __gen_read(x, start, iter) do {									\
				for (i = (start), k = 0, ext_coor[x] = -10; i >= 0 && i < ks->seq.l && k < s[x]; iter) {	\
					int c = target[i], mut_type = c & mutmsk;			\
					if (ext_coor[x] < 0) {								\
						if (mut_type != NOCHANGE && mut_type != SUBSTITUTE) continue; \
						ext_coor[x] = i;								\
					}													\
					if (mut_type == DELETE) ++n_indel[x];				\
					else if (mut_type == NOCHANGE || mut_type == SUBSTITUTE) { \
						tmp_seq[x][k++] = c & 0xf;						\
						if (mut_type == SUBSTITUTE) ++n_sub[x];			\
					} else {											\
						int n, ins;										\
						++n_indel[x];									\
						tmp_seq[x][k++] = c & 0xf;						\
						for (n = mut_type>>12, ins = c>>4; n > 0 && k < s[x]; --n, ins >>= 2) \
							tmp_seq[x][k++] = ins & 0x3;				\
					}													\
				}														\
				if (k != s[x]) ext_coor[x] = -10;						\
			} while (0)

			__gen_read(0, pos, ++i);
			__gen_read(1, pos + d - 1, --i);
			for (k = 0; k < s[1]; ++k) tmp_seq[1][k] = tmp_seq[1][k] < 4? 3 - tmp_seq[1][k] : 4; // complement
			if (ext_coor[0] < 0 || ext_coor[1] < 0) { // fail to generate the read(s)
				--ii;
				continue;
			}

			// generate sequencing errors
			for (j = 0; j < 2; ++j) {
				int n_n = 0;
				for (i = 0; i < s[j]; ++i) {
					int c = tmp_seq[j][i];
					if (c >= 4) { // actually c should be never larger than 4 if everything is correct
						c = 4;
						++n_n;
					} else if (xdrand() < ERR_RATE) { // drand48
						// c = (c + (int)(xdrand() * 3.0 + 1)) & 3; // random sequencing errors // drand48
						c = (c + 1) & 3; // recurrent sequencing errors
						++n_err[j];
					}
					tmp_seq[j][i] = c;
				}
				if ((double)n_n / s[j] > MAX_N_RATIO) break;
			}
			if (j < 2) { // too many ambiguous bases on one of the reads
				--ii;
				continue;
			}

			// print
			for (k = 0; k < 2; ++k) {
                j = order[k];
				for (i = 0; i < s[j]; ++i) qstr[i] = Q;
				qstr[i] = 0;
				fprintf(fpo[j], "@%s:%s_%u_%u_%d:%d:%d_%d:%d:%d_%llx/%d\n", ref->id, ks->name.s, ext_coor[0]+1, ext_coor[1]+1,
						n_err[0], n_sub[0], n_indel[0], n_err[1], n_sub[1], n_indel[1],
						(long long)ii, j==0? is_flip+1 : 2-is_flip);
				for (i = 0; i < s[j]; ++i)
					fputc("ACGTN"[(int)tmp_seq[j][i]], fpo[j]);
				fprintf(fpo[j], "\n+\n%s\n", qstr);
			}
            actual_n_pairs++;
		}
		free(rseq[0].s); free(rseq[1].s);
	}
	kseq_destroy(ks);
	gzclose(fp_fa);
	free(qstr);
	free(tmp_seq[0]); free(tmp_seq[1]);
    return actual_n_pairs;
}

static int simu_usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: wgsim (short read simulator)\n");
	fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
	fprintf(stderr, "Contact: Andrew Tritt <ajtritt@lbl.gov>\n\n");
	fprintf(stderr, "Usage:   wgsim [options] <in.nwk> <num_taxa>\n\n");
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
	return 1;
}

int main(int argc, char *argv[])
{
	int64_t N;
	int dist, std_dev, c, size_l, size_r, is_hap = 0, ileave = 0;
	FILE *fpout1, *fpout2, *about, *mutout;
	int seed = -1;
    int n_taxa, i;
    char * indir = ".";
    char * outdir = ".";

	N = 1000000; dist = 270; std_dev = 50;
	size_l = size_r = 150;
	while ((c = getopt(argc, argv, "e:d:s:N:1:2:r:R:hX:S:A:f:o:I")) >= 0) {
		switch (c) {
		case 'd': dist = atoi(optarg); break;
		case 's': std_dev = atoi(optarg); break;
		case 'N': N = atoi(optarg); break;
		case '1': size_l = atoi(optarg); break;
		case '2': size_r = atoi(optarg); break;
		case 'e': ERR_RATE = atof(optarg); break;
		case 'r': MUT_RATE = atof(optarg); break;
		case 'R': INDEL_FRAC = atof(optarg); break;
		case 'X': INDEL_EXTEND = atof(optarg); break;
		case 'A': MAX_N_RATIO = atof(optarg); break;
		case 'S': seed = atoi(optarg); break;
		case 'h': is_hap = 1; break;
		case 'f': indir = optarg; break;
		case 'o': outdir = optarg; break;
		case 'I': ileave = 1; break;
		}
	}

    struct stat sb;
    if (!(stat(outdir, &sb) == 0 && S_ISDIR(sb.st_mode))){
        mkdir(outdir, 0755);
    }


	if (argc - optind < 2) return simu_usage();
    n_taxa = atoi(argv[optind+1]);
    char outpath[128];
    if (ileave == 1) {
        logmsg("[wgsim] Interleaving output");
        snprintf(outpath, 128, "%s/reads.fastq", outdir);
	    fpout1 = fopen(outpath, "w");
        fpout2 = fpout1;
    } else {
        snprintf(outpath, 128, "%s/reads_1.fastq", outdir);
	    fpout1 = fopen(outpath, "w");
        snprintf(outpath, 128, "%s/reads_2.fastq", outdir);
	    fpout2 = fopen(outpath, "w");
    }
    snprintf(outpath, 128, "%s/mutations.txt", outdir);
    mutout = fopen(outpath, "w");

	if (!fpout1 || !fpout2) {
		logmsg("[wgsim] file open error");
		return 1;
	}
	//if (seed <= 0) seed = time(0)&0x7fffffff;
	if (seed <= 0) seed = time(NULL);
    srand(seed);
	xsrand(seed);
	logmsg("[wgsim] seed = %d", seed);
    logmsg("[wgsim] sampling tree %s, getting sequence from %s", argv[optind], indir);
    logmsg("[wgsim] fragment size mean = %d, fragment size stdev = %d", dist, std_dev);
    char ** leaves = sample_tree(argv[optind], n_taxa);
    logmsg("[wgsim] calculating abundances");
    double * abund = calc_abund(n_taxa);
    double Ltot = 0.0;
    ref_t * refs = (ref_t*) malloc(n_taxa*sizeof(ref_t));
    ref_t * tmp_ref = refs;
    char * fn;
    for (i = 0; i < n_taxa; i++) {
        fn = (char *) malloc(128*sizeof(char));
        snprintf(fn, 128, "%s/%s.fasta", indir, leaves[i]);
        logmsg("[wgsim] Looking for %s in %s.", leaves[i], fn);
        tmp_ref->id = leaves[i];
        read_seq(fn, tmp_ref, dist + 3 * std_dev);
        tmp_ref->wlen = abund[i] * tmp_ref->tot_len;
        Ltot += tmp_ref->wlen;
        tmp_ref++;
    }
    tmp_ref = refs;
    snprintf(outpath, 128, "%s/abundance.txt", outdir);
    about = fopen(outpath, "w");
    uint64_t n_pairs;
    uint64_t actual;
    uint64_t actual_total = 0;
    for (i = 0; i < n_taxa; i++){
        n_pairs = (uint64_t) tmp_ref->wlen*N/Ltot;
        logmsg("[wgsim] Sampling %lld (%0.12f) pairs from %s", (long long) n_pairs, abund[i], tmp_ref->fn);
	    actual = wgsim_core(fpout1, fpout2, mutout, tmp_ref, is_hap, n_pairs, dist, std_dev, size_l, size_r);
        fprintf(about, "%s\t%s\t%0.12f\t%lld\n", tmp_ref->id, tmp_ref->fn, abund[i], (long long) actual);
        actual_total += actual;
        tmp_ref++;
    }
    logmsg("[wgsim] Done sampling reads, simulated %lld actual reads", (long long) actual_total);

    fclose(about);
    fclose(mutout);
	fclose(fpout1);
    if (ileave != 1) {
        fclose(fpout2);
    }
    logmsg("[wgsim] Closed outputs");
	return 0;
}
