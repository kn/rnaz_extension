#ifndef __VIENNA_RNA_PACKAGE_PART_FUNC_H__
#define __VIENNA_RNA_PACKAGE_PART_FUNC_H__

#include "data_structures.h"

#define FLT_OR_DBL double

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#else
#define DEPRECATED(func) func
#endif


/**
 *  \file boltzmann.h
 * 
 *  \brief Partition function of single RNA sequences
 * 
 *  This file includes (almost) all function declarations within the <b>RNAlib</b> that are related to
 *  Partion function folding...
 */

/**
 *  a flag indicating that auxilary arrays are needed throughout the computations which are necessary for stochastic backtracking
 */
extern  int st_back;

/*
#################################################
# PARTITION FUNCTION COMPUTATION                #
#################################################
*/

/**
 *  \brief Compute the partition function \f$Q\f$ of an RNA sequence
 * 
 *  If \a structure is not a NULL pointer on input, it contains on
 *  return a string consisting of the letters " . , | { } ( ) " denoting
 *  bases that are essentially unpaired, weakly paired, strongly paired without
 *  preference, weakly upstream (downstream) paired, or strongly up-
 *  (down-)stream paired bases, respectively.
 *  If #fold_constrained is not 0, the \a structure string is
 *  interpreted on input as a list of constraints for the folding. The
 *  character "x" marks bases that must be unpaired, matching brackets " ( ) "
 *  denote base pairs, all other characters are ignored. Any pairs
 *  conflicting with the constraint will be forbidden. This is usually sufficient
 *  to ensure the constraints are honored.
 *  If #do_backtrack has been set to 0 base pairing probabilities will not
 *  be computed (saving CPU time), otherwise #pr will contain the probability
 *  that bases \a i and \a j pair.
 *  \note The global array #pr is deprecated and the user who wants the computed
 *  base pair probabilities for further computations is advised to use the function export_bppm()
 * 
 *  \see pf_circ_fold(), bppm_to_structure(), export_bppm()
 * 
 *  \param sequence   The RNA sequence to be computed
 *  \param structure  A pointer to a char array where a base pair probability information might be stored in a pseudo-dot-bracket notation (might be NULL, too)
 *  \returns          The Gibbs free energy of the ensemble (\f$G = -RT \cdot \log(Q) \f$) in kcal/mol
 */
float   boltzmann(const char *sequence, char *structure, double mfe, FLT_OR_DBL *Q, FLT_OR_DBL *X, FLT_OR_DBL *Y, double ali_mfe, double *mfe_prob);

//float   pf_fold(const char *sequence, char *structure);

/**
 *  \brief Free arrays from pf_fold()
 */
void    free_pf_arrays(void);

#endif
