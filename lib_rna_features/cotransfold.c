/*
								cotransfold.c
							Author: Katsuya Noguchi
	
	Calculates Cis and Trans of a given RNA sequence and its structure.
	
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"
#include "cotransfold.h"

#define PUBLIC
#define PRIVATE static

#define MIN_STEM_LEN 9 /* The minimum stem length for alternative helix. */
#define TRUE 1
#define FALSE 0
#define NO_PAIR 0
#define AU 1
#define UA 2
#define CG 3
#define GC 4
#define GU 5
#define UG 6

PUBLIC void cotransfold(char *seq, int numBps, bp_info *bps, double *cis, double *trans, int is_plain);

PRIVATE double calculate_threeCis_or_threeTrans(int i, int c, int len, double energy, int is_plain);
PRIVATE double calculate_fiveCis_or_fiveTrans(int i, int c, double energy, int is_plain);
PRIVATE int get_stem_len(double *energy, const char *seq, int len, int b1, int b2, int is_forward);
PRIVATE int has_at_least_two_adjacent_bp(int numBps, bp_info *bps, int p);
PRIVATE int is_alt_helix(double *energy, const char *seq, int len, int b1, int b2);
PRIVATE int which_base_pair(char b1, char b2);
PRIVATE double get_stacking_energy(char pre_b1, char pre_b2, char b1, char b2);
PRIVATE double get_au_stacking_energy(char b1, char b2);
PRIVATE double get_cg_stacking_energy(char b1, char b2);
PRIVATE double get_gc_stacking_energy(char b1, char b2);
PRIVATE double get_gu_stacking_energy(char b1, char b2);
PRIVATE double get_ua_stacking_energy(char b1, char b2);
PRIVATE double get_ug_stacking_energy(char b1, char b2);


void cotransfold(char *seq, int numBps, bp_info *bps, double *cis, double *trans, int is_plain) {
	int len = (int) strlen(seq);
	int p, c;
	double threeCis = 0, threeTrans = 0, fiveCis = 0, fiveTrans = 0;
	int *pairs = (int *)malloc(sizeof(int) * len);
	double *energy = (double *)malloc(sizeof(double));

	for (p = 0; p < numBps; p++) {
		if (has_at_least_two_adjacent_bp(numBps, bps, p) == TRUE) {
			for (c = 0; c < bps[p].bp.i - 3; c++) {
				if (is_alt_helix(energy, seq, len, c, bps[p].bp.i) == TRUE) {
					fiveCis += calculate_fiveCis_or_fiveTrans(bps[p].bp.i, c, *energy, is_plain);
				}
			}
   			for (c = bps[p].bp.j + 1; c < len; c++) {
				if (is_alt_helix(energy, seq, len, bps[p].bp.i, c) == TRUE) {
					threeTrans += calculate_threeCis_or_threeTrans(bps[p].bp.j, c, len, *energy, is_plain);
				}
			}
			for (c = 0; c < bps[p].bp.i; c++) {
				if (is_alt_helix(energy, seq, len, c, bps[p].bp.j) == TRUE) {
					fiveTrans += calculate_fiveCis_or_fiveTrans(bps[p].bp.i, c, *energy, is_plain);					
				}
			}
			for (c = bps[p].bp.j + 4; c < len; c++) {
				if (is_alt_helix(energy, seq, len, bps[p].bp.j, c) == TRUE) {
					threeCis += calculate_threeCis_or_threeTrans(bps[p].bp.j, c, len, *energy, is_plain);
				}
			}
		}
	}

	*cis = fiveCis - threeCis;
	*trans = threeTrans - fiveTrans;
}

PRIVATE double calculate_threeCis_or_threeTrans(int i, int c, int len, double energy, int is_plain) {
	double d = (double) c - i;
	double l = (double) len - i;
	return (is_plain ? 1 : abs(energy)) / (d * log(l));
}

PRIVATE double calculate_fiveCis_or_fiveTrans(int i, int c, double energy, int is_plain) {
	double d = (double) i - c;
	double l = (double) i + 1;
	return (is_plain ? 1 : abs(energy)) / (d * log(l));
}

/**
 * Returns true if given base pairs have at least two adjacent base pairs; o/w false.
 */
PRIVATE int has_at_least_two_adjacent_bp(int numBps, bp_info *bps, int p) {
	int count = 0;
	int k = 1;
	while (p - k >= 0 && count < 2 && bps[p - k].bp.i + k == bps[p].bp.i && bps[p - k].bp.j - k == bps[p].bp.j) {
		count++;
		k++;
	}
	k = 1;
	while (p + k < numBps && count < 2 && bps[p + k].bp.i - k == bps[p].bp.i && bps[p + k].bp.j + k == bps[p].bp.j) {
		count++;
		k++;
	}
	if (count < 2)
		return FALSE;
	return TRUE;
}

/**
 * Returns true if given bases are qualified to be alternative helix; o/w false.
 */
PRIVATE int is_alt_helix(double *energy, const char *seq, int len, int b1, int b2) {
	int stem_len = 0;
	*energy = 0;
	stem_len += get_stem_len(energy, seq, len, b1, b2, TRUE);
	stem_len += get_stem_len(energy, seq, len, b1, b2, FALSE);
	stem_len--;
	if (stem_len >= MIN_STEM_LEN)
		return TRUE;
	return FALSE;
}

/**
 * Recursively calculate stem length in specified direction and its stacking energy.
 * Returns the stem length.
 */
PRIVATE int get_stem_len(double *energy, const char *seq, int len, int b1, int b2, int is_forward) {
	if (b1 < 0 || b2 >= len || b1 >= b2 || which_base_pair(seq[b1], seq[b2]) == NO_PAIR)
		return 0;
	int stem_len;
	if (is_forward == TRUE) {
		stem_len = get_stem_len(energy, seq, len, b1 + 1, b2 - 1, is_forward) + 1;
		if (stem_len > 1)
			*energy += get_stacking_energy(seq[b1], seq[b2], seq[b1 + 1], seq[b2 - 1]);
	} else {
		stem_len = get_stem_len(energy, seq, len, b1 - 1, b2 + 1, is_forward) + 1;
		if (stem_len > 1)
			*energy += get_stacking_energy(seq[b1 - 1], seq[b2 + 1], seq[b1], seq[b2]);
	}
	return stem_len;
}

/**
 * Returns a type of base pair of given bases.
 */
PRIVATE int which_base_pair(char b1, char b2) {
	switch (b1) {
	case 'A': case 'a':
		if (b2 == 'U' || b2 == 'u')
			return AU;
		break;
	case 'C': case 'c':
		if (b2 == 'G' || b2 == 'g')
			return CG;
		break;	
	case 'G': case 'g':
		if (b2 == 'C' || b2 == 'c')
			return GC;
		if (b2 == 'U' || b2 == 'u')
			return GU;
		break;
	case 'U': case 'u':
		if (b2 == 'A' || b2 == 'a')
			return UA;
		if (b2 == 'G' || b2 == 'g')
			return UG;
		break;
	}
	return NO_PAIR;
}

/**
 * Returns stacking energy of given base stacking.
 */
PRIVATE double get_stacking_energy(char pre_b1, char pre_b2, char b1, char b2) {
	switch (which_base_pair(pre_b1, pre_b2)) {
	case AU:
		return get_au_stacking_energy(b1, b2);
	case CG:
		return get_cg_stacking_energy(b1, b2);
	case GC:
		return get_gc_stacking_energy(b1, b2);
	case GU:
		return get_gu_stacking_energy(b1, b2);
	case UA:
		return get_ua_stacking_energy(b1, b2);
	case UG:
		return get_ug_stacking_energy(b1, b2);
	}
	return 0;
}

PRIVATE double get_au_stacking_energy(char b1, char b2) {
	switch (which_base_pair(b1, b2)) {
	case AU:
		return -0.9;
	case CG:
		return -2.2;
	case GC:
		return -2.1;
	case GU:
		return -0.6;
	case UA:
		return -1.1;
	case UG:
		return -1.4;
	}
	return 0;
}

PRIVATE double get_cg_stacking_energy(char b1, char b2) {
	switch (which_base_pair(b1, b2)) {
	case AU:
		return -2.1;
	case CG:
		return -3.3;
	case GC:
		return -2.4;
	case GU:
		return -1.4;
	case UA:
		return -2.1;
	case UG:
		return -2.1;
	}
	return 0;
}

PRIVATE double get_gc_stacking_energy(char b1, char b2) {
	switch (which_base_pair(b1, b2)) {
	case AU:
		return -2.4;
	case CG:
		return -3.4;
	case GC:
		return -3.3;
	case GU:
		return -1.5;
	case UA:
		return -2.2;
	case UG:
		return -2.5;
	}
	return 0;
}

PRIVATE double get_gu_stacking_energy(char b1, char b2) {
	switch (which_base_pair(b1, b2)) {
	case AU:
		return -1.3;
	case CG:
		return -2.5;
	case GC:
		return -2.1;
	case GU:
		return -0.5;
	case UA:
		return -1.4;
	case UG:
		return 1.3;
	}
	return 0;
}

PRIVATE double get_ua_stacking_energy(char b1, char b2) {
	switch (which_base_pair(b1, b2)) {
	case AU:
		return -1.3;
	case CG:
		return -2.4;
	case GC:
		return -2.1;
	case GU:
		return -1.0;
	case UA:
		return -0.9;
	case UG:
		return -1.3;
	}
	return 0;
}

PRIVATE double get_ug_stacking_energy(char b1, char b2) {
	switch (which_base_pair(b1, b2)) {
	case AU:
		return -1.0;
	case CG:
		return -1.5;
	case GC:
		return -1.4;
	case GU:
		return 0.3;
	case UA:
		return -0.6;
	case UG:
		return -0.5;
	}
	return 0;
}