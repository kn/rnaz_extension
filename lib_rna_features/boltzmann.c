/*
							boltzmann.c
 Caslculates expected energy and variance of boltzmann distribution

                  			Katsuya Noguchi
*/

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>    /* #defines FLT_MAX ... */
#include "utils.h"
#include "energy_par.h"
#include "fold_vars.h"
#include "pair_mat.h"
#include "params.h"
#include "loop_energies.h"
#include "boltzmann.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define ISOLATED  256.0

/*
#################################
# GLOBAL VARIABLES              #
#################################
*/
PUBLIC  int     st_back=0;

/*
#################################
# PRIVATE VARIABLES             #
#################################
*/
PRIVATE FLT_OR_DBL  *q, *qb=NULL, *qm, *qm1, *qqm, *qqm1, *qq, *qq1;
PRIVATE FLT_OR_DBL  *x, *xb=NULL, *xm, *xm1, *xxm, *xxm1, *xx, *xx1;
PRIVATE FLT_OR_DBL  *y, *yb=NULL, *ym, *ym1, *yym, *yym1, *yy, *yy1;
PRIVATE FLT_OR_DBL  *probs, *prml, *prm_l, *prm_l1, *q1k, *qln;
PRIVATE FLT_OR_DBL  *scale;
PRIVATE FLT_OR_DBL  *expMLbase;
PRIVATE FLT_OR_DBL  qo, qho, qio, qmo, *qm2;
PRIVATE int         *jindx;
PRIVATE int         init_length;  /* length in last call to init_pf_fold() */
PRIVATE int         circular=0;
PRIVATE char        *pstruc;
PRIVATE char        *sequence;
PRIVATE char        *ptype;       /* precomputed array of pair types */
PRIVATE pf_paramT   *pf_params;   /* the precomputed Boltzmann weights */
PRIVATE short       *S, *S1;

#ifdef _OPENMP

/* NOTE: all variables are assumed to be uninitialized if they are declared as threadprivate
         thus we have to initialize them before usage by a seperate function!
         OR: use copyin in the PARALLEL directive!
         e.g.:
         #pragma omp parallel for copyin(pf_params)
*/
#pragma omp threadprivate(q, qb, qm, qm1, qqm, qqm1, qq, qq1, prml, prm_l, prm_l1, q1k, qln,\
                          probs, scale, expMLbase, qo, qho, qio, qmo, qm2, jindx, init_length,\
                          circular, pstruc, sequence, ptype, pf_params, S, S1)

#endif

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/
PRIVATE void  init_partfunc(int length);
PRIVATE void  scale_pf_params(unsigned int length);
PRIVATE void  get_arrays(unsigned int length);
PRIVATE void  make_ptypes(const short *S, const char *structure);
PRIVATE void  pf_linear(const char *sequence, char *structure);
PRIVATE void  pf_create_bppm(const char *sequence, char *structure);
PRIVATE void  backtrack(int i, int j);
PRIVATE void  backtrack_qm(int i, int j);
PRIVATE void  backtrack_qm1(int i,int j);
PRIVATE void  backtrack_qm2(int u, int n);

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PRIVATE void init_partfunc(int length){
  if (length<1) nrerror("init_pf_fold: length must be greater 0");

#ifdef _OPENMP
/* Explicitly turn off dynamic threads */
  omp_set_dynamic(0);
  free_pf_arrays(); /* free previous allocation */
#else
  if (init_length>0) free_pf_arrays(); /* free previous allocation */
#endif

#ifdef SUN4
  nonstandard_arithmetic();
#else
#ifdef HP9
  fpsetfastmode(1);
#endif
#endif
  make_pair_matrix();
  get_arrays((unsigned) length);
  scale_pf_params((unsigned) length);

  init_length = length;
}

PRIVATE void get_arrays(unsigned int length){
  unsigned int size;

  size  = sizeof(FLT_OR_DBL) * ((length+1)*(length+2)/2);

  q     = (FLT_OR_DBL *) space(size);
  qb    = (FLT_OR_DBL *) space(size);
  qm    = (FLT_OR_DBL *) space(size);
  qm1   = (st_back || circular) ? (FLT_OR_DBL *) space(size) : NULL;
  x     = (FLT_OR_DBL *) space(size);
  xb    = (FLT_OR_DBL *) space(size);
  xm    = (FLT_OR_DBL *) space(size);
  xm1   = (st_back || circular) ? (FLT_OR_DBL *) space(size) : NULL;
  y     = (FLT_OR_DBL *) space(size);
  yb    = (FLT_OR_DBL *) space(size);
  ym    = (FLT_OR_DBL *) space(size);
  ym1   = (st_back || circular) ? (FLT_OR_DBL *) space(size) : NULL;
  probs = (FLT_OR_DBL *) space(size);

  ptype     = (char *) space(sizeof(char)*((length+1)*(length+2)/2));
  q1k       = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+1));
  qln       = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+2));
  qq        = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+2));
  qq1       = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+2));
  qqm       = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+2));
  qqm1      = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+2));
  xx        = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+2));
  xx1       = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+2));
  xxm       = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+2));
  xxm1      = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+2));
  yy        = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+2));
  yy1       = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+2));
  yym       = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+2));
  yym1      = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+2));
  prm_l     = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+2));
  prm_l1    = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+2));
  prml      = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+2));
  expMLbase = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+1));
  scale     = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+1));

  iindx     = get_iindx(length);
  jindx     = get_indx(length);
}

/**
*** Allocate memory for all matrices and other stuff
**/
PUBLIC void free_pf_arrays(void){
  if(q)         free(q);
  if(qb)        free(qb);
  if(qm)        free(qm);
  if(qm1)       free(qm1);
  if(x)         free(x);
  if(xb)        free(xb);
  if(xm)        free(xm);
  if(xm1)       free(xm1);
  if(y)         free(y);
  if(yb)        free(yb);
  if(ym)        free(ym);
  if(ym1)       free(ym1);
  if(qm2)       free(qm2);
  if(ptype)     free(ptype);
  if(qq)        free(qq);
  if(qq1)       free(qq1);
  if(qqm)       free(qqm);
  if(qqm1)      free(qqm1);
  if(xx)        free(xx);
  if(xx1)       free(xx1);
  if(xxm)       free(xxm);
  if(xxm1)      free(xxm1);
  if(yy)        free(yy);
  if(yy1)       free(yy1);
  if(yym)       free(yym);
  if(yym1)      free(yym1);
  if(q1k)       free(q1k);
  if(qln)       free(qln);
  if(probs)     free(probs);
  if(prm_l)     free(prm_l);
  if(prm_l1)    free(prm_l1);
  if(prml)      free(prml);
  if(expMLbase) free(expMLbase);
  if(scale)     free(scale);
  if(iindx)     free(iindx);
  if(jindx)     free(jindx);
  if(S)         free(S);
  if(S1)        free(S1);
  if(pr)        free(pr);

  S = S1 = NULL;
  q = x = y = pr = probs = qb = qm = qm1 = xb = xm = xm1 = yb = ym = ym1 = qm2 = qq = qq1 = qqm = qqm1 = xx = xx1 = xxm = xxm1 = yy = yy1 = yym = yym1 = q1k = qln = prm_l = prm_l1 = prml = expMLbase = scale = NULL;
  iindx = jindx = NULL;
  ptype = NULL;

#ifdef SUN4
  standard_arithmetic();
#else
#ifdef HP9
  fpsetfastmode(0);
#endif
#endif

  init_length = 0;
}

/*-----------------------------------------------------------------*/
PUBLIC float boltzmann(const char *sequence, char *structure, double mfe, FLT_OR_DBL *Q, FLT_OR_DBL *X, FLT_OR_DBL *Y, double ali_mfe, double *mfe_prob){

  int         n = (int) strlen(sequence);

#ifdef _OPENMP
  /* always init everything since all global static variables are uninitialized when entering a thread */
  init_partfunc(n);
#else
  if (n > init_length) init_partfunc(n);
#endif

  S   = encode_sequence(sequence, 0);
  S1  = encode_sequence(sequence, 1);

  make_ptypes(S, structure);

  /* do the linear pf fold and fill all matrices  */
  pf_linear(sequence, structure);

  *Q = q[iindx[1]-n];
  *X = x[iindx[1]-n];
  *Y = y[iindx[1]-n];

  *mfe_prob += log10(exp(-ali_mfe * 1000 / pf_params->kT) / *Q);

  return log10(exp(-mfe * 1000 / pf_params->kT) / *Q);
}

PRIVATE void pf_linear(const char *sequence, char *structure){

  int n, i,j,k,l, ij, u,u1,d,ii, type, type_2, tt;
  FLT_OR_DBL temp, tempx, tempy, Qmax=0;
  FLT_OR_DBL qbt1, xbt1, ybt1, *tmp, *tmpx, *tmpy;
  FLT_OR_DBL expEnergy, energy;
  FLT_OR_DBL xi, qi, yi;

  FLT_OR_DBL  expMLclosing = pf_params->expMLclosing;
  double      max_real;
  double kT = pf_params->kT / -1000;

  max_real = (sizeof(FLT_OR_DBL) == sizeof(float)) ? FLT_MAX : DBL_MAX;

  n = (int) strlen(sequence);

  /*array initialization ; qb,qm,q
    qb,qm,q (i,j) are stored as ((n+1-i)*(n-i) div 2 + n+1-j */

  for (d=0; d<=TURN; d++)
    for (i=1; i<=n-d; i++) {
      j=i+d;
      ij = iindx[i]-j;
      q[ij]=1.0; //*scale[d+1];
	  x[ij]=0.0;
	  y[ij]=0.0;
      qb[ij]=qm[ij]=0.0;
	  xb[ij]=xm[ij]=0.0;
	  yb[ij]=ym[ij]=0.0;
    }

  for (i=1; i<=n; i++)
    qq[i]=qq1[i]=qqm[i]=qqm1[i]=xx[i]=xx1[i]=xxm[i]=xxm1[i]=yy[i]=yy1[i]=yym[i]=yym1[i]=0;

  for (j=TURN+2;j<=n; j++) {
    for (i=j-TURN-1; i>=1; i--) {
      /* construction of partition function of segment i,j*/
      /*firstly that given i binds j : qb(i,j) */
      u = j-i-1; ij = iindx[i]-j;
      type = ptype[ij];
      if (type!=0) {
        /*hairpin contribution*/
        if (((type==3)||(type==4))&&no_closingGU) qbt1 = 0;
        else  {
			expEnergy = exp_E_Hairpin(u, type, S1[i+1], S1[j-1], sequence+i-1, pf_params);
			energy = log(expEnergy) * kT;
			qbt1 = expEnergy; //* scale[u+2];
			xbt1 = expEnergy * energy;
			ybt1 = xbt1 * energy;
        }
		/* interior loops with interior pair k,l */
        for (k=i+1; k<=MIN2(i+MAXLOOP+1,j-TURN-2); k++) {
          u1 = k-i-1;
          for (l=MAX2(k+TURN+1,j-1-MAXLOOP+u1); l<j; l++) {
            type_2 = ptype[iindx[k]-l];
            if (type_2) {
              type_2 = rtype[type_2];
			  expEnergy = exp_E_IntLoop(u1, j-l-1, type, type_2,
										S1[i+1], S1[j-1], S1[k-1], S1[l+1], pf_params);
			  energy = log(expEnergy) * kT;
              qbt1 += qb[iindx[k]-l] * (expEnergy);
			  xbt1 += expEnergy * (xb[iindx[k]-l] + energy * qb[iindx[k]-l]);
			  ybt1 += (yb[iindx[k]-l] * expEnergy + expEnergy * energy * energy * qb[iindx[k]-l]
					+ 2.0 * xb[iindx[k]-l] * expEnergy * energy);
            }
          }
        }
        /*multiple stem loop contribution*/
        ii = iindx[i+1]; /* ii-k=[i+1,k-1] */
        temp = tempx = tempy = 0.0;
        for (k=i+2; k<=j-1; k++) {
			temp += qm[ii-(k-1)] * qqm1[k];
			tempx += (xm[ii-(k-1)]*qqm1[k] + xxm1[k]*qm[ii-(k-1)]);
			tempy += (ym[ii-(k-1)]*qqm1[k] + yym1[k]*qm[ii-(k-1)] + 2.0*xm[ii-(k-1)]*xxm1[k]);
		}
        tt = rtype[type];
		
		expEnergy = expMLclosing;
		energy = log(expEnergy) * kT;
		qi = temp * expEnergy;
		xi = expEnergy * (tempx + energy * temp);
		yi = tempy * expEnergy + expEnergy*energy*energy * temp + 2.0 * tempx * expEnergy*energy;
		
		expEnergy = exp_E_MLstem(tt, S1[j-1], S1[i+1], pf_params);
		energy = log(expEnergy) * kT;
		
		qbt1 += qi * expEnergy; //* scale[2];
		xbt1 += expEnergy * (xi + energy * qi);
		ybt1 += yi * expEnergy + expEnergy*energy*energy * qi + 2.0 * xi * expEnergy*energy;
        qb[ij] = qbt1;
		xb[ij] = xbt1;
		yb[ij] = ybt1;
      } /* end if (type!=0) */
      else {
		qb[ij] = 0.0;
		xb[ij] = 0.0;
		yb[ij] = 0.0;
	  }

      /* construction of qqm matrix containing final stem
         contributions to multiple loop partition function
         from segment i,j */
	  expEnergy = expMLbase[1];
	  energy = log(expEnergy) * kT;
      qqm[i] = qqm1[i]*expEnergy;
	  xxm[i] = expEnergy * (xxm1[i] + energy * qqm1[i]);
	  yym[i] = yym1[i] * expEnergy + expEnergy*energy*energy * qqm1[i] + 2.0 * xxm1[i] * expEnergy*energy;
      if (type) {
		expEnergy = exp_E_MLstem(type, ((i>1) || circular) ? S1[i-1] : -1, ((j<n) || circular) ? S1[j+1] : -1, pf_params);
		energy = log(expEnergy) * kT;
		qbt1 = qb[ij] * expEnergy;
		xbt1 = expEnergy * (xb[ij] + energy * qb[ij]);
		ybt1 = yb[ij] * expEnergy + expEnergy*energy*energy * qb[ij] + 2.0 * xb[ij] * expEnergy*energy;
        qqm[i] += qbt1;
		xxm[i] += xbt1;
		yym[i] += ybt1;
      }
      if (qm1) qm1[jindx[j]+i] = qqm[i]; /* for stochastic backtracking and circfold */

      /*construction of qm matrix containing multiple loop
        partition function contributions from segment i,j */
      temp = tempx = tempy = 0.0;
      ii = iindx[i];  /* ii-k=[i,k-1] */
      for (k=j; k>i; k--) {
		expEnergy = expMLbase[k-i];
		energy = log(expEnergy) * kT;
		qi = qm[ii-(k-1)] + expEnergy;
		temp += qi*qqm[k];
		xi = xm[ii-(k-1)] + expEnergy * energy;
		tempx += (xi * qqm[k] + xxm[k] * qi);
		yi = ym[ii-(k-1)] + expEnergy * energy * energy;
		tempy += (yi * qqm[k] + yym[k] * qi + 2.0 * xi * xxm[k]);
	  }
      qm[ij] = (temp + qqm[i]);
   	  xm[ij] = (tempx + xxm[i]);
	  ym[ij] = (tempy + yym[i]);

      /*auxiliary matrix qq for cubic order q calculation below */
      qbt1 = qb[ij];
	  xbt1 = xb[ij];
	  ybt1 = yb[ij];
      if(type) {
		expEnergy = exp_E_ExtLoop(type, ((i>1) || circular) ? S1[i-1] : -1, ((j<n) || circular) ? S1[j+1] : -1, pf_params);
		energy = log(expEnergy) * kT;
		ybt1 = ybt1 * expEnergy + expEnergy*energy*energy * qbt1 + 2.0 * xbt1 * expEnergy*energy;
		xbt1 = expEnergy * (xbt1 + energy * qbt1);
		qbt1 *= expEnergy; 
	  }

      qq[i] = qq1[i] + qbt1; //*scale[1] + qbt1;
	  xx[i] = xx1[i] + xbt1;
	  yy[i] = yy1[i] + ybt1;

      /*construction of partition function for segment i,j */
	  temp = 1.0 + qq[i]; //*scale[1+j-i] + qq[i];
	  tempx = xx[i];
	  tempy = yy[i];
      for (k=i; k<=j-1; k++) {
		temp += q[ii-k]*qq[k+1];
		tempx += (x[ii-k] * qq[k+1] + xx[k+1] * q[ii-k]);
		tempy += (y[ii-k] * qq[k+1] + yy[k+1] * q[ii-k] + 2.0 * x[ii-k] * xx[k+1]);
	  }
      q[ij] = temp;
	  x[ij] = tempx;
	  y[ij] = tempy;
      if (temp>Qmax) {
        Qmax = temp;
        if (Qmax>max_real/10.)
          fprintf(stderr, "Q close to overflow: %d %d %g\n", i,j,temp);
      }
      if (temp>=max_real) {
        PRIVATE char msg[128];
        sprintf(msg, "overflow in pf_fold while calculating q[%d,%d]\n"
                     "use larger pf_scale", i,j);
        nrerror(msg);
      }
    }
    tmp = qq1;  qq1 =qq;  qq =tmp;
    tmp = qqm1; qqm1=qqm; qqm=tmp;
	tmpx = xx1;  xx1 =xx;  xx =tmpx;
    tmpx = xxm1; xxm1=xxm; xxm=tmpx;
	tmpy = yy1;  yy1 =yy;  yy =tmpy;
    tmpy = yym1; yym1=yym; yym=tmpy;
  }
}

PRIVATE void scale_pf_params(unsigned int length){
  unsigned int i;
  double  kT;

  if(pf_params) free(pf_params);
  pf_params = get_scaled_pf_parameters();

  kT = pf_params->kT;   /* kT in cal/mol  */

   /* not using scale factors anymore
	scaling factors (to avoid overflows)
  if (pf_scale == -1) { /* mean energy for random sequences: 184.3*length cal
    pf_scale = exp(-(-185+(pf_params->temperature-37.)*7.27)/kT);
    if (pf_scale<1) pf_scale=1;
  }
  scale[0] = 1.;
  scale[1] = 1./pf_scale;
  */
  expMLbase[0] = 1;
  expMLbase[1] = pf_params->expMLbase;
  for (i=2; i<=length; i++) {
    //scale[i] = scale[i/2]*scale[i-(i/2)];
	expMLbase[i] = pow(pf_params->expMLbase, (double)i);// * scale[i];
  }
}

/*---------------------------------------------------------------------------*/
PRIVATE void make_ptypes(const short *S, const char *structure){
  int n,i,j,k,l;

  n=S[0];
  for (k=1; k<n-TURN; k++)
    for (l=1; l<=2; l++) {
      int type,ntype=0,otype=0;
      i=k; j = i+TURN+l; if (j>n) continue;
      type = pair[S[i]][S[j]];
      while ((i>=1)&&(j<=n)) {
        if ((i>1)&&(j<n)) ntype = pair[S[i-1]][S[j+1]];
        if (noLonelyPairs && (!otype) && (!ntype))
          type = 0; /* i.j can only form isolated pairs */
        qb[iindx[i]-j] = 0.;
        ptype[iindx[i]-j] = (char) type;
        otype =  type;
        type  = ntype;
        i--; j++;
      }
    }

  if (fold_constrained && (structure != NULL))
    constrain_ptypes(structure, (unsigned int)n, ptype, NULL, TURN, 1);
}

/*
  stochastic backtracking in pf_fold arrays
  returns random structure S with Boltzman probabilty
  p(S) = exp(-E(S)/kT)/Z
*/
char *pbacktrack(char *seq){
  double r, qt;
  int i,j,n, start;
  sequence = seq;
  n = strlen(sequence);

  if (init_length<1)
    nrerror("can't backtrack without pf arrays.\n"
            "Call pf_fold() before pbacktrack()");
  pstruc = space((n+1)*sizeof(char));

  for (i=0; i<n; i++) pstruc[i] = '.';

  start = 1;
  while (start<n) {
  /* find i position of first pair */
    for (i=start; i<n; i++) {
      r = urn() * qln[i];
      if (r > qln[i+1]*scale[1])  break; /* i is paired */
    }
    if (i>=n) break; /* no more pairs */
    /* now find the pairing partner j */
    r = urn() * (qln[i] - qln[i+1]*scale[1]);
    for (qt=0, j=i+1; j<=n; j++) {
      int type;
      type = ptype[iindx[i]-j];
      if (type) {
        double qkl;
        qkl = qb[iindx[i]-j];
        if (j<n) qkl *= qln[j+1];
        qkl *=  exp_E_ExtLoop(type, (i>1) ? S1[i-1] : -1, (j<n) ? S1[j+1] : -1, pf_params);
        qt += qkl;
        if (qt > r) break; /* j is paired */
      }
    }
    if (j==n+1) nrerror("backtracking failed in ext loop");
    start = j+1;
    backtrack(i,j);
  }

  return pstruc;
}
char *pbacktrack_circ(char *seq){
  double r, qt;
  int i, j, k, l, n;
  FLT_OR_DBL  expMLclosing      = pf_params->expMLclosing;

  sequence = seq;
  n = strlen(sequence);

  if (init_length<1)
    nrerror("can't backtrack without pf arrays.\n"
      "Call pf_circ_fold() before pbacktrack_circ()");
  pstruc = space((n+1)*sizeof(char));

  /* initialize pstruct with single bases  */
  for (i=0; i<n; i++) pstruc[i] = '.';

  qt = 1.0*scale[n];
  r = urn() * qo;

  /* open chain? */
  if(qt > r) return pstruc;

  for(i=1; (i < n); i++){
    for(j=i+TURN+1;(j<=n); j++){

      int type, u;
      /* 1. first check, wether we can do a hairpin loop  */
      u = n-j + i-1;
      if (u<TURN) continue;

      type = ptype[iindx[i]-j];
      if (!type) continue;

      type=rtype[type];

      char loopseq[10];
      if (u<7){
        strcpy(loopseq , sequence+j-1);
        strncat(loopseq, sequence, i);
      }

      qt += qb[iindx[i]-j] * exp_E_Hairpin(u, type, S1[j+1], S1[i-1],  loopseq, pf_params) * scale[u];
      /* found a hairpin? so backtrack in the enclosed part and we're done  */
      if(qt>r){ backtrack(i,j); return pstruc;}

      /* 2. search for (k,l) with which we can close an interior loop  */
      for(k=j+1; (k < n); k++){
        int ln1, lstart;
        ln1 = k - j - 1;
        if(ln1+i-1>MAXLOOP) break;

        lstart = ln1+i-1+n-MAXLOOP;
        if(lstart<k+TURN+1) lstart = k + TURN + 1;
        for(l=lstart; (l <= n); l++){
            int ln2, type2;
            ln2 = (i - 1) + (n - l);
            if((ln1+ln2) > MAXLOOP) continue;

            type2 = ptype[iindx[k]-l];
            if(!type) continue;
            type2 = rtype[type2];
            qt += qb[iindx[i]-j] * qb[iindx[k]-l] * exp_E_IntLoop(ln2, ln1, type2, type, S1[l+1], S1[k-1], S1[i-1], S1[j+1], pf_params) * scale[ln1 + ln2];
            /* found an exterior interior loop? also this time, we can go straight  */
            /* forward and backtracking the both enclosed parts and we're done      */
            if(qt>r){ backtrack(i,j); backtrack(k,l); return pstruc;}
        }
      } /* end of kl double loop */
    }
  } /* end of ij double loop  */
  {
    /* as we reach this part, we have to search for our barrier between qm and qm2  */
    qt = 0.;
    r = urn()*qmo;
    for(k=TURN+2; k<n-2*TURN-3; k++){
      qt += qm[iindx[1]-k] * qm2[k+1] * expMLclosing;
      /* backtrack in qm and qm2 if we've found a valid barrier k  */
      if(qt>r){ backtrack_qm(1,k); backtrack_qm2(k+1,n); return pstruc;}
    }
  }
  /* if we reach the actual end of this function, an error has occured  */
  /* cause we HAVE TO find an exterior loop or an open chain!!!         */
  nrerror("backtracking failed in exterior loop");
  return pstruc;
}

PRIVATE void backtrack_qm(int i, int j){
  /* divide multiloop into qm and qm1  */
  double qmt, r;
  int k;
  while(j>i){
    /* now backtrack  [i ... j] in qm[] */
    r = urn() * qm[iindx[i] - j];
    qmt = qm1[jindx[j]+i]; k=i;
    if(qmt<r)
      for(k=i+1; k<=j; k++){
        qmt += (qm[iindx[i]-(k-1)]+expMLbase[k-i])*qm1[jindx[j]+k];
        if(qmt >= r) break;
      }
    if(k>j) nrerror("backtrack failed in qm");

    backtrack_qm1(k,j);

    if(k<i+TURN) break; /* no more pairs */
    r = urn() * (qm[iindx[i]-(k-1)] + expMLbase[k-i]);
    if(expMLbase[k-i] >= r) break; /* no more pairs */
    j = k-1;
  }
}

PRIVATE void backtrack_qm1(int i,int j){
  /* i is paired to l, i<l<j; backtrack in qm1 to find l */
  int ii, l, type;
  double qt, r;
  r = urn() * qm1[jindx[j]+i];
  ii = iindx[i];
  for (qt=0., l=i+TURN+1; l<=j; l++) {
    type = ptype[ii-l];
    if (type)
      qt +=  qb[ii-l] * exp_E_MLstem(type, S1[i-1], S1[l+1], pf_params) * expMLbase[j-l];
    if (qt>=r) break;
  }
  if (l>j) nrerror("backtrack failed in qm1");
  backtrack(i,l);
}

PRIVATE void backtrack_qm2(int k, int n){
  double qom2t, r;
  int u;
  r= urn()*qm2[k];
  /* we have to search for our barrier u between qm1 and qm1  */
  for (qom2t = 0.,u=k+TURN+1; u<n-TURN-1; u++){
    qom2t += qm1[jindx[u]+k]*qm1[jindx[n]+(u+1)];
    if(qom2t > r) break;
  }
  if(u==n-TURN) nrerror("backtrack failed in qm2");
  backtrack_qm1(k,u);
  backtrack_qm1(u+1,n);
}

PRIVATE void backtrack(int i, int j){
  do {
    double r, qbt1;
    int k, l, type, u, u1;

    pstruc[i-1] = '('; pstruc[j-1] = ')';

    r = urn() * qb[iindx[i]-j];
    type = ptype[iindx[i]-j];
    u = j-i-1;
    /*hairpin contribution*/
    if (((type==3)||(type==4))&&no_closingGU) qbt1 = 0;
    else
      qbt1 = exp_E_Hairpin(u, type, S1[i+1], S1[j-1], sequence+i-1, pf_params)*
        scale[u+2]; /* add scale[u+2] */

    if (qbt1>=r) return; /* found the hairpin we're done */

    for (k=i+1; k<=MIN2(i+MAXLOOP+1,j-TURN-2); k++) {
      u1 = k-i-1;
      for (l=MAX2(k+TURN+1,j-1-MAXLOOP+u1); l<j; l++) {
        int type_2;
        type_2 = ptype[iindx[k]-l];
        if (type_2) {
          type_2 = rtype[type_2];
          /* add *scale[u1+u2+2] */
          qbt1 += qb[iindx[k]-l] * (scale[u1+j-l+1] *
            exp_E_IntLoop(u1, j-l-1, type, type_2,
                          S1[i+1], S1[j-1], S1[k-1], S1[l+1], pf_params));
        }
        if (qbt1 > r) break;
      }
      if (qbt1 > r) break;
    }
    if (l<j) {
      i=k; j=l;
    }
    else break;
  } while (1);

  /* backtrack in multi-loop */
  {
    double r, qt;
    int k, ii, jj;

    i++; j--;
    /* find the first split index */
    ii = iindx[i]; /* ii-j=[i,j] */
    jj = jindx[j]; /* jj+i=[j,i] */
    for (qt=0., k=i+1; k<j; k++) qt += qm[ii-(k-1)]*qm1[jj+k];
    r = urn() * qt;
    for (qt=0., k=i+1; k<j; k++) {
      qt += qm[ii-(k-1)]*qm1[jj+k];
      if (qt>=r) break;
    }
    if (k>=j) nrerror("backtrack failed, can't find split index ");

    backtrack_qm1(k, j);

    j = k-1;
    backtrack_qm(i,j);
  }
}