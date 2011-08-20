/*********************************************************************
 *                                                                   *
 *                              zscore.h                             *
 *                                                                   *
 *	Compute a z-score to assess significance of a predicted MFE  *
 *                                                                   *
 *	          c Stefan Washietl, Ivo L Hofacker                  *
 *                                                                   *
 *	   $Id: zscore.h,v 1.3 2004/09/19 13:31:42 wash Exp $        *
 *                                                                   *
 *********************************************************************/

void regression_svm_init();

void regression_svm_free();

double mfe_zscore(const char *seq, double mfe, int *type, int avoid_shuffle, char* warning_string);
