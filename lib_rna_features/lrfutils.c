/*
			utils.c
	Author: Katsuya Noguchi
	
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "lrfutils.h"

#define PRIVATE static
#define PUBLIC

PUBLIC void *allocate(unsigned int size, char* name);
PUBLIC void find_base_pairs(const char* structure, int* numBps, bp_info* bps);
PRIVATE void sort_base_pairs(int len, bp_info* bps);
PRIVATE void quick_sort(bp_info* bps, int low, int high);
PRIVATE int partition(bp_info* bps, int low, int high, int pivot);
PRIVATE void swap(bp_info* bps, int i, int j);

PUBLIC void *allocate(unsigned int size, char* name) {
	void *pointer;
	
	if ((pointer = (void*) malloc(size)) == NULL) {
		fprintf(stderr, "Failed to allocate memory for %s.", name);
		exit(1);
	}
	return pointer;
}

/**
 * Determine base pairs and types of loops they close,
 * given its parenthesis structure.
 */
PUBLIC void find_base_pairs(const char* structure, int* numBps, bp_info* bps) {
	int i, j, c, l;
	int len = strlen(structure);
	int size = len/2 + 1;
	int *stack = (int *)allocate(sizeof(int) * len, "stack");
	int *numItems = (int *)allocate(sizeof(int) * len, "numItem");
	tuple t;
	tuple loops[size][size];
	
	for (i = 0; i < size; i++)
		numItems[i] = 0;
	
	j = c = 0;
	l = 1;
	for (i = 0; i < len; i++) {
		if (structure[i] == '(') {
			stack[j++] = i;
			l++;
		} else if (structure[i] == ')') {
			if (j == 0) {
				fprintf(stderr, "Unbalanced structure is given.\n");
				exit(1);
			}
			bps[c].bp.i = stack[--j];
			bps[c].bp.j = i;
			t = loops[--l][0];
			// Check if current base pair closes an interior loop.
			if (numItems[l] == 1 && (t.i != bps[c].bp.i + 1 || t.j != bps[c].bp.j - 1)) {
				t.i = t.i - bps[c].bp.i - 1;
				t.j = bps[c].bp.j - t.j - 1;
				bps[c].loop.i = t.i > t.j ? t.j : t.i;
				bps[c].loop.j = t.i > t.j ? t.i : t.j;
			}
			numItems[l] = 0;
			loops[l - 1][numItems[l - 1]++] = bps[c].bp;
			c++;			
		}
	}
	if (j != 0) {
		fprintf(stderr, "Unbalanced structure is given.\n");
		exit(1);
	}
	sort_base_pairs(c, bps);
	for (i = c - 1; i >= 0; i--) {
		if (bps[i].loop.i != 0 || bps[i].loop.j != 0) {
			bps[i + 1].loop.i = bps[i].loop.i;
			bps[i + 1].loop.j = bps[i].loop.j;
		}
	}
	(*numBps) = c;
}

/**
 * Sorts given array of base pairs by index number of the first base.
 **/
PRIVATE void sort_base_pairs(int len, bp_info* bps) {
	srand ( time(NULL) );
	quick_sort(bps, 0, len - 1);
}

PRIVATE void quick_sort(bp_info* bps, int low, int high) {
	int len = high - low + 1;
	if (len <= 1)
		return;
	int pivot = rand() % len + low;
	pivot = partition(bps, low, high, pivot);	
	quick_sort(bps, low, pivot - 1);
	quick_sort(bps, pivot + 1, high);
}

PRIVATE int partition(bp_info* bps, int low, int high, int pivot) {
	int i, storeIndex;
	storeIndex = low;
	swap(bps, pivot, high);
	pivot = high;
	for (i = low; i < high; i++) {
		if (bps[i].bp.i < bps[pivot].bp.i) {
			swap(bps, i, storeIndex++);
		}
	}
	swap(bps, storeIndex, pivot);
	return storeIndex;
}

PRIVATE void swap(bp_info* bps, int i, int j) {
	bp_info tmp = bps[i];
	bps[i] = bps[j];
	bps[j] = tmp;
}