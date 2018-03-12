#include <stdio.h>  //fprintf
#include <stdlib.h> //free qsort
#include <zlib.h>
#include <inttypes.h>
#include <math.h> // pow()
#include <stdint.h>

// klib (H.li)
#include "kseq.h"
#include "kvec.h"
#include "kstring.h"

KSEQ_INIT(gzFile, gzread)

#define VERSION "0.0.1"

typedef kvec_t(double) double_array;

double log2(double n) {
  // log(n)/log(2) is log2.
  return log(n) / log(2);
}

double sgn(double x) {
  if (x > 0) return 1;
  if (x < 0) return -1;
  return 0;
}

double mean(double *counts, int n) {
  double sum = 0;
  int i = 0;
  for (; i < n; i++) {
    sum += counts[i];
  }
  return sum / n;
}

double logmean(int *counts, int n) {
  double sum = 0;
  int i = 0;
  for (; i < n; i++) {
    sum += log(counts[i]);
  }
  return sum / n;
}

double sd(double *counts, int n, double mean) {
  double sum = 0;
  int i = 0;
  for (; i < n; i++) {
    sum += pow(counts[i] - mean, 2);
  }
  return sqrt(sum / (n-1));
}

int cmp_reverse_double_pointers(const void * a, const void * b) {
  const double aa = **(const double **)a;
  const double bb = **(const double **)b;
  if (aa > bb) {
    return -1;
  }
  if (bb > aa) {
    return  1;
  }
  return 0;
}

int compare_double(const void *a,const void *b) {
  double *x = (double *) a;
  double *y = (double *) b;
  if (*x < *y) return -1;
  else if (*x > *y) return 1; return 0;
}

double dmin(double a, double b) {
  if (a > b) {
    return b;
  } else {
    return a;
  }
}

int main(int argc, char *argv[])
{
  char *counts_file;
  double sampling_rate = 1;

  int c;
  while ((c = getopt(argc, argv, "s:")) >= 0) {
    switch (c) {
      case 's': sampling_rate = atof(optarg); break;
    }
  }

  if(sampling_rate <= 0 || sampling_rate > 1) {
    fprintf(stderr, "Invalid value for sampling rate [%.2f], must be comprised between 0 and 1.\n", sampling_rate);
    return 1;
  }

  if ((optind) >= argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   computeNF [options] <counts.tsv>\n\n");
		fprintf(stderr, "Options: -s FLOAT  sampling rate [%.2f]\n", sampling_rate);
		fprintf(stderr, "\n");
		return 1;
	}

  counts_file     = argv[optind++];

  gzFile fp;
	kstream_t *ks;
	kstring_t *str;
  kvec_t(char*) samples;
  kvec_t(double_array) norm_counts;
  int dret = 0, sampling_modulo = (int) 1 / sampling_rate;
  size_t j, line = 0, nb_samples = 0;

  str = (kstring_t*)calloc(1, sizeof(kstring_t));
  kv_init(samples);
  kv_init(norm_counts);

  // 1. Get samples names from counts_file

  fp = gzopen(counts_file, "r");
  if(!fp) { fprintf(stderr, "Failed to open %s\n", counts_file); exit(EXIT_FAILURE); }
  ks = ks_init(fp);

  ks_getuntil(ks, KS_SEP_SPACE, str, &dret);
  while(ks_getuntil(ks, KS_SEP_SPACE, str, &dret) >= 0) {
    kv_push(char*, samples, ks_release(str));
    // Create a new array for the counts of this samples
    double_array sample_counts;
    kv_init(sample_counts);
    kv_push(double_array, norm_counts, sample_counts);
    if (dret == '\n') break;
  }
  ks_destroy(ks);
  gzclose(fp);

  nb_samples = kv_size(samples);

  // 1. Get samples conditions indicies and Normalization factors

  // Open counts file
  fp = gzopen(counts_file, "r");
  if(!fp) { fprintf(stderr, "Failed to open %s\n", counts_file); exit(EXIT_FAILURE); }
  ks = ks_init(fp);

  int *counts = (int*) malloc(sizeof(int) * nb_samples);

  // skip header line
  ks_getuntil(ks, KS_SEP_LINE, str, &dret);
  while(ks_getuntil(ks, KS_SEP_SPACE, str, &dret) >= 0) {
    line++;

    // Go to next line depending of the sampling rate
    if(line % sampling_modulo != 0) {
      // Skip the rest of the line
      ks_getuntil(ks, KS_SEP_LINE, str, &dret);
      // Go to next loop iteration
      continue;
    }

    char *kmer = ks_release(str);

    // load counts
    j = 0;
    while(j < nb_samples && ks_getuntil(ks, KS_SEP_SPACE, str, &dret) >= 0) {
      //fprintf(stderr, "TOTO\n");
      counts[j] = atoi(str->s);
      j++;
    }

    // Initial code from compute_norm_factors.R
    //loggeomeans <- rowMeans(log(selected_kmers_counts[,2:ncol(selected_kmers_counts)]))
    //function(cnts) { exp(median((log(cnts) - loggeomeans)[is.finite(loggeomeans) & cnts > 0]))})

    // Compute raw mean
    double log_row_mean = logmean(counts, nb_samples);
    double value = 0;

    for(j = 0; j < nb_samples; j++) {
      if(counts[j] > 0 && isfinite(log_row_mean)) {
        value = log(counts[j]) - log_row_mean;
        kv_push(double,kv_A(norm_counts,j),value);
      }
    }

    if (dret != '\n') {
      fprintf(stderr, "inconsistent number of column (line %zu)\n", line);
      exit(EXIT_FAILURE);
    }

    // Free k-mer
    free(kmer);
  }
  ks_destroy(ks);
  gzclose(fp);

  fprintf(stdout, "sample\tnormalization_factor\n");

  for(j = 0; j < nb_samples; j++) {
    double_array a = kv_A(norm_counts,j);
    qsort(a.a, kv_size(a), sizeof(double), compare_double);
    double median = kv_A(a,int(kv_size(a)/2));
    double nf = exp(median);
    fprintf(stdout, "%s\t%f\n", kv_A(samples,j), nf);
    kv_destroy(a);
  }

  return 0;
}
