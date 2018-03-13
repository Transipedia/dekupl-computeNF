#include <stdio.h>  //fprintf
#include <stdlib.h> //free qsort
#include <math.h> // log()
#include <zlib.h>

// klib (H.li)
#include "kseq.h"
#include "kvec.h"
#include "kstring.h"

KSEQ_INIT(gzFile, gzread)

#define VERSION "0.0.1"

typedef kvec_t(double) double_array;
typedef kvec_t(double_array) double_array_collection;

double logmean(int *counts, int n) {
  double sum = 0;
  int i = 0;
  for (; i < n; i++) {
    sum += log(counts[i]);
  }
  return sum / n;
}

int compare_double(const void *a,const void *b) {
  double *x = (double *) a;
  double *y = (double *) b;
  if (*x < *y) {
    return -1;
  } else if (*x > *y) {
    return 1;
  }
  return 0;
}

double* compute_nf(double_array_collection norm_counts, size_t nb_samples) {
  double * normalization_factors = (double*) malloc(sizeof(double) * nb_samples);
  for(size_t j = 0; j < nb_samples; j++) {
    double_array a = kv_A(norm_counts,j);
    qsort(a.a, kv_size(a), sizeof(double), compare_double);
    double median = kv_A(a,int(kv_size(a)/2));
    normalization_factors[j] = exp(median);
  }
  return normalization_factors;
}

int main(int argc, char *argv[])
{
  char *counts_file;
  int nb_kmers_per_step = 1000000;
  double min_error_per_sample = 0.001;

  if ((optind) >= argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   computeNF [options] <counts.tsv>\n\n");
		return 1;
	}

  counts_file     = argv[optind++];

  gzFile fp;
	kstream_t *ks;
	kstring_t *str;
  kvec_t(char*) samples;
  double_array_collection norm_counts;
  int dret = 0;
  size_t j, line = 0, nb_samples = 0, nb_step = 0;
  double* normalization_factors;
  double min_error;

  str = (kstring_t*)calloc(1, sizeof(kstring_t));
  kv_init(samples);
  kv_init(norm_counts);

  /* 1. Get samples names from counts_file */

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
  normalization_factors = (double*) malloc(sizeof(double) * nb_samples);
  min_error = min_error_per_sample * nb_samples;
  /* 2. For each sample store counts normalized with log row mean */

  // Open counts file
  fp = gzopen(counts_file, "r");
  if(!fp) { fprintf(stderr, "Failed to open %s\n", counts_file); exit(EXIT_FAILURE); }
  ks = ks_init(fp);

  int *counts = (int*) malloc(sizeof(int) * nb_samples);

  // skip header line
  ks_getuntil(ks, KS_SEP_LINE, str, &dret);
  while(ks_getuntil(ks, KS_SEP_SPACE, str, &dret) >= 0) {
    line++;

    char *kmer = ks_release(str);

    // load counts
    j = 0;
    while(j < nb_samples && ks_getuntil(ks, KS_SEP_SPACE, str, &dret) >= 0) {
      counts[j] = atoi(str->s);
      j++;
    }

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

    // We have reached a step, we compute the NF
    if(line % nb_kmers_per_step == 0) {
      fprintf(stderr, "compute normalization_factors\n");
      // Create a new array to store the NF
      double * new_nf = compute_nf(norm_counts, nb_samples);
      // We already have done one step, we compare the values
      if(nb_step > 0) {
        double sum = 0;
        for(j = 0; j < nb_samples; j++) {
          sum += abs(normalization_factors[j] - new_nf[j]);
        }
        fprintf(stderr, "step %zu - error with previous step : %f\n", nb_step, sum);
        // New factors replace the old ones
        free(normalization_factors);
        normalization_factors = new_nf;
        // We have found a good approximation we stop the computation
        if(sum < min_error) {
          break;
        }
      }
      nb_step++;
    }

    free(kmer);
  }
  ks_destroy(ks);
  gzclose(fp);

  // We have never calculated the NF
  if(nb_step == 0) {
    normalization_factors = compute_nf(norm_counts, nb_samples);
  }

  fprintf(stdout, "sample\tnormalization_factor\n");

  /* 3. Compute median for each sample and print normalization factor */
  for(j = 0; j < nb_samples; j++) {
    fprintf(stdout, "%s\t%f\n", kv_A(samples,j), normalization_factors[j]);
    kv_destroy(kv_A(norm_counts,j));
  }

  return 0;
}
