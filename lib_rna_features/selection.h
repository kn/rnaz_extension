#include "lrfutils.h"

extern void classify_base_pairs(int numBps, bp_info* bps, int len, int* classes);
extern void get_avg_sub_rate_for_classes(int numSeq, int len, char** seqs, int *classes,
										double *avg1, double *avg2, double *avg3, double *avg13);