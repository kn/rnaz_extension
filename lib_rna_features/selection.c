/*
		selection.c
	Author: Katsuya Noguchi

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"
#include "dnapars.h"
#include "selection.h"

#define PUBLIC
#define PRIVATE static

#define NUM_CLASS 3
#define MIN_STEM_LEN 4
#define MIN_LOOP_SIZE 3
#define IGNORE_CLASS_UP_TO 2
#define TRUE 1
#define FALSE 0

PUBLIC void get_avg_sub_rate_for_classes(int numSeq, int len, char** seqs, int *classes,
										double *avg1, double *avg2, double *avg3, double *avg13);
PUBLIC void classify_base_pairs(int numBps, bp_info* bps, int len, int *classes);
PRIVATE void separate_sequences_by_class(int numSeq, int len, int* classes, char** seqs,
										int *len1, char** seqC1, int *len2,
										char** seqC2, int *len3, char** seqC3, int *len13, char** seqC13);
PRIVATE int annotate_bp(tuple bp, tuple loop, int stemLen, int class);

/**
 * Determines the class of bases in give base pairs.
 */
PUBLIC void classify_base_pairs(int numBps, bp_info* bps, int len, int* classes) {
	int i, j, k;
	tuple outLoop, inLoop;
	
	for (i = 0; i < len; i++)
		classes[i] = 0; //NUM_CLASS + 1;
	
	i = 0;
	while (i < numBps) {
		j = 1;
		// Determine stem length
		while (i + j < numBps && bps[i + j].bp.i == bps[i].bp.i + j && bps[i + j].bp.j == bps[i].bp.j - j)
			j++;
		
		// Determine loop type on outside of stem
		outLoop = bps[i].loop;
		
		// Determine loop type on inside of stem
		inLoop = bps[i + j - 1].loop;
		
		// Annotate base pairs in current stem
		for (k = 0; k < j / 2; k++) {
			classes[bps[i + k].bp.i] = classes[bps[i + k].bp.j] = annotate_bp(bps[i + k].bp, outLoop, j, k);
			classes[bps[i + j - 1 - k].bp.i] = classes[bps[i + j - 1 - k].bp.j] = annotate_bp(bps[i + j - 1 - k].bp, inLoop, j, k);
		}
		// Annotate middle base pair in stem of odd length
		if (j % 2 == 1) {
			classes[bps[i + j / 2].bp.i] = annotate_bp(bps[i + j / 2].bp, outLoop, j, j / 2);
			if (classes[bps[i + j / 2].bp.i]  != -1)
				classes[bps[i + j / 2].bp.i] = annotate_bp(bps[i + j / 2].bp, inLoop, j, j / 2);
			classes[bps[i + j / 2].bp.j] = classes[bps[i + j / 2].bp.i];
			
		}
		i += j;
	}
}

/**
 * Calculates the average substitution rate for each class (I, II and III).
 */
PUBLIC void get_avg_sub_rate_for_classes(int numSeq, int len, char** seqs, int *classes,
										double *avg1, double *avg2, double *avg3, double *avg13) {
	int i = 0;
	int *len1 = (int*)malloc(sizeof(int));
	int *len2 = (int*)malloc(sizeof(int));
	int *len3 = (int*)malloc(sizeof(int));
	int *len13 = (int*)malloc(sizeof(int));
	char** seqC1 = (char**)malloc(sizeof(char*) * numSeq);
	char** seqC2 = (char**)malloc(sizeof(char*) * numSeq);
	char** seqC3 = (char**)malloc(sizeof(char*) * numSeq);
	char** seqC13 = (char**)malloc(sizeof(char*) * numSeq);
	
	for (i = 0; i < numSeq; i++) {
		seqC1[i] = (char*)malloc(sizeof(char) * len);
		seqC2[i] = (char*)malloc(sizeof(char) * len);
		seqC3[i] = (char*)malloc(sizeof(char) * len);
		seqC13[i] = (char*)malloc(sizeof(char) * len);
	}
	separate_sequences_by_class(numSeq, len, classes, seqs,
								len1, seqC1, len2,
								seqC2, len3, seqC3, len13, seqC13);
	(*avg1) = get_min_parsimony_score(numSeq, *len1, seqC1) / (double)(*len1);
	(*avg2) = get_min_parsimony_score(numSeq, *len2, seqC2) / (double)(*len2);
	(*avg3) = get_min_parsimony_score(numSeq, *len3, seqC3) / (double)(*len3);
	(*avg13) = get_min_parsimony_score(numSeq, *len13, seqC13) / (double)(*len13);
}

/**
 * Separates given set of sequences into sets of sequences by class (1, 2 and 3). 
 */
PRIVATE void separate_sequences_by_class(int numSeq, int len, int* classes, char** seqs,
										int *len1, char** seqC1, int *len2,
										char** seqC2, int *len3, char** seqC3,
										int *len13, char** seqC13) {
	int c1, c2, c3, c13, i, j;
	c1 = c2 = c3 = c13 = 0;
	for (i = 0; i < len; i++) {
		if (classes[i] == 1) {
			for (j = 0; j < numSeq; j++) {
				seqC1[j][c1] = seqs[j][i];
				seqC13[j][c13] = seqs[j][i];
			}
			c13++;
			c1++;
		} else if (classes[i] == 2) {
			for (j = 0; j < numSeq; j++)
				seqC2[j][c2] = seqs[j][i];
			c2++;
		} else if (classes[i] == 3) {
			for (j = 0; j < numSeq; j++) {
				seqC3[j][c3] = seqs[j][i];
				seqC13[j][c13] = seqs[j][i];
			}
			c13++;
			c3++;
		}
	}
	for (j = 0; j < numSeq; j++) {
		seqC1[j][c1] = '\0';
		seqC2[j][c2] = '\0';
		seqC3[j][c3] = '\0';
		seqC13[j][c13] = '\0';
	}
	(*len1) = c1;
	(*len2) = c2;
	(*len3) = c3;
	(*len13) = c13;
}

/**
 * Determine the class of given base pair.
 */
PRIVATE int annotate_bp(tuple bp, tuple loop, int stemLen, int class) {
	class = class < NUM_CLASS ? class + 1 : NUM_CLASS;
	if (stemLen < MIN_STEM_LEN) {
		// stem length is less than minimum
		return -1;
	} else if (loop.i != 0 || loop.j != 0) {
		// loop is internal loop
		if (loop.i + loop.j >= MIN_LOOP_SIZE || class > IGNORE_CLASS_UP_TO) {
			return class;
		} else {
			return -1;
		}
	} else {
		// loop is exterior or other loop
		return class;
	}
}