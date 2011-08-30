/*********************************************************************                
 *                                                                   *
 *                          svm_helper.c                             *
 *                                                                   *
 *   Functions relating to SVM regression/classification and         *
 *    for interacting with the SVMLIB libraries                      *
 *                                                                   *
 *	     c Stefan Washietl, Ivo L Hofacker, Katsuya Noguchi          *
 *                                                                   *
 *                                                                   *
 *********************************************************************/

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "svm.h"
#include "svm_helper.h"
#include "rnaz_utils.h"
#include "utils.h"

#include "mfe_stdv.inc"
#include "mfe_avg.inc"
#include "decision.inc"
#include "decision_dinucleotide.inc"
#include "decision_dinucleotide_structural.inc"
#include "decision_rna_features.inc"


/* Here we read in the models as strings */
#include "L50-200_GC20-30_avg.inc"
#include "L50-200_GC20-30_stdv.inc"
#include "L50-200_GC30-36_avg.inc"
#include "L50-200_GC30-36_stdv.inc"
#include "L50-200_GC36-40_avg.inc"
#include "L50-200_GC36-40_stdv.inc"
#include "L50-200_GC40-46_avg.inc"
#include "L50-200_GC40-46_stdv.inc"
#include "L50-200_GC46-50_avg.inc"
#include "L50-200_GC46-50_stdv.inc"
#include "L50-200_GC50-56_avg.inc"
#include "L50-200_GC50-56_stdv.inc"
#include "L50-200_GC56-60_avg.inc"
#include "L50-200_GC56-60_stdv.inc"
#include "L50-200_GC60-66_avg.inc"
#include "L50-200_GC60-66_stdv.inc"
#include "L50-200_GC66-70_avg.inc"
#include "L50-200_GC66-70_stdv.inc"
#include "L50-200_GC70-80_avg.inc"
#include "L50-200_GC70-80_stdv.inc"



/* All values in the SVMs are normalized to mean 0 and standard
   deviation 1. Therefore the mean and standard deviation of the
   training set must be known to scale the values appropriately.*/

/*                              Mean  Std.Dev. */
double scale_regression[][2]={{0.5,0.1581213081},    /* GC-content  */
			      {0.5,0.1581213081},    /* "A-content" */
			      {0.5,0.1581213081},    /* "C-content" */
			      {225.0,114.569772373}};/* length      */

double scale_decision[][2]= {{-7.87,2.76 }, /* z   */
			     {0,1.23 },     /* SCI */
			     {52.35,99.55}, /* ID  */
			     {2,6}};        /* k   */

double scale_decision_dinucleotide[][2]= {{-8.15,2.01 }, /* z   */
					  {0,1.29 },     /* SCI */
					  {0,1.28718}}; /* norm. Shannon entropy  */

double scale_decision_dinucleotide_structural[][2]= {{-8.13,2.01 }, /* z   */
						     {0,1.31 },     /* SCI */
						     {0,1.38353}}; /* norm. Shannon entropy  */
						
double scale_decision_rna_features[][2] = {
							{-6.18, 2.73},	/* z */
							{0, 1.28},	/* SCI */
							{0.00619, 0.82296},	/* norm. Shannon entropy */
							{-2.36296, 2.71455},	/* CIS */
							{-97.76376999999999, 26.64685},	/* TRANS */
							{-5.233691, -0.213344},	/* MFE log probability */
							{0.010874, 0.050835},	/* The difference between MFE and the Boltzmann expected energy (per base) */
							{0.000163, 0.001514}};	/* The Boltzmann variance (per base)

/* Also y-values are scaled in the regression SVMs, so also here the
   mean and standard deviation are needed */

double scale_avg[]={-58.60276,45.24618};
double scale_stdv[]={4.098457,1.107606};


/* Scales a SVM-lib "node" according to the scale_arrays defined
   above. Index 1 of the node is scaled with values of index 0 of the
   array.*/

void scale_node(struct svm_node* node,double scale[][2]){

  int i=0;
  double mu,sigma;
  while (1){
	mu=scale[node[i].index-1][0];
	sigma=scale[node[i].index-1][1];
	node[i].value=(node[i].value-mu)/sigma;
	++i;
	if (node[i].index==-1) break;
  }
}

void scale_node_linear(struct svm_node* node,double scale[][2]){

  int i=0;
  double min,max;
  double from,to;
  from=-1.0;
  to=1.0;
  while (1){
	min=scale[node[i].index-1][0];
	max=scale[node[i].index-1][1];
	node[i].value=from+(to-from)*(node[i].value-min)/(max-min);
	++i;
	if (node[i].index==-1) break;
  }
}



/* Scales x-data for the z-score regression */

void scale_regression_node(struct svm_node* node){
  scale_node(node,scale_regression);
}


void scale_decision_node(struct svm_node* node, int decision_model_type){
  if (decision_model_type == 1) {
    scale_node_linear(node,scale_decision);
  }
  if (decision_model_type == 2) {
    scale_node_linear(node,scale_decision_dinucleotide);
  }
  if (decision_model_type == 3) {
    scale_node_linear(node,scale_decision_dinucleotide_structural);
  }
  if (decision_model_type == 4) {
	scale_node_linear(node,scale_decision_rna_features);
  }
}



/* Back-scales y-data to original dimensions */

void backscale_regression(double* avg, double* stdv){
  *avg=(*avg)*scale_avg[1]+scale_avg[0];
  *stdv=(*stdv)*scale_stdv[1]+scale_stdv[0];
}


/* Loads both models for average and standard deviation from files
   given by a common basename. If no name is given, default models
   hard-coded in this file are used */

void get_regression_models(struct svm_model** avg_model,
			   struct svm_model** stdv_model,
			   int type){

  if (type == -1) {
    *avg_model=svm_load_model_string(default_avg_model_string);
    *stdv_model=svm_load_model_string(default_stdv_model_string);
  }
  if (type == 0) {
    *avg_model=svm_load_model_string(avg_20_30_string);
    *stdv_model=svm_load_model_string(stdv_20_30_string);
  }
  if (type == 1) {
    *avg_model=svm_load_model_string(avg_30_36_string);
    *stdv_model=svm_load_model_string(stdv_30_36_string);
    }
  if (type == 2) { 
    *avg_model=svm_load_model_string(avg_36_40_string);
    *stdv_model=svm_load_model_string(stdv_36_40_string);
  }
  if (type == 3) {
    *avg_model=svm_load_model_string(avg_40_46_string);
    *stdv_model=svm_load_model_string(stdv_40_46_string);
  }
  if (type == 4) {
    *avg_model=svm_load_model_string(avg_46_50_string);
    *stdv_model=svm_load_model_string(stdv_46_50_string);
  }
  if (type == 5) {
    *avg_model=svm_load_model_string(avg_50_56_string);
    *stdv_model=svm_load_model_string(stdv_50_56_string);
  }
  if (type == 6) {
    *avg_model=svm_load_model_string(avg_56_60_string);
    *stdv_model=svm_load_model_string(stdv_56_60_string);
  }
  if (type == 7) {
      *avg_model=svm_load_model_string(avg_60_66_string);
      *stdv_model=svm_load_model_string(stdv_60_66_string);
  }
  if (type == 8) {
    *avg_model=svm_load_model_string(avg_66_70_string);
    *stdv_model=svm_load_model_string(stdv_66_70_string);
  }
  if (type == 9) {
    *avg_model=svm_load_model_string(avg_70_80_string);
    *stdv_model=svm_load_model_string(stdv_70_80_string);
  }
}


struct svm_model* get_decision_model(char *basefilename, int decision_model_type){

  char fn[256];

  struct svm_model* model=NULL;
  
  if (basefilename!=NULL){
    strcpy(fn, basefilename); strcat(fn, "/decision.model");
    model=svm_load_model(fn);

  } else {
    if (decision_model_type == 1) {
	model=svm_load_model_string(default_decision_model_string);
    }
    if (decision_model_type == 2) {
	model=svm_load_model_string(decision_dinucleotide_string);
    }
    if (decision_model_type == 3) {
	model=svm_load_model_string(decision_dinucleotide_structural_string);
    }
	if (decision_model_type == 4) {
	model=svm_load_model_string(decision_rna_features_model_string);
	}
  }
  
  return model;
  
}


/*********************************************************************
 *********************************************************************

 The folllowing code enables us to read a libsvm model from a string
 rather than a file. This should be part of the library, but currently
 it is not. So here this function is implemented.

 Model files can easily be modified in C-strings using this perl
 one-liner (plus adding a "char* filename=" and a ";" at the end):

 cat svm.model |perl -ne 'chomp;print "\"$_\\n\"\n"' >svm_model.inc

 The new function needs the structure definition for "svm_model" and
 the two arrays "svm_type_table" and "kernel_type_table" which are
 coded in svm.cpp but are not exported throgh svm.h. Therefore they
 are here redefined.

*/


struct svm_model{
  struct svm_parameter param;
  int nr_class;
  int l;
  struct svm_node **SV;
  double **sv_coef;
  double *rho;
  double *probA;
  double *probB;
  int *label;
  int *nSV;
  int free_sv;
};

struct svm_model* svm_load_model_string(char *modelString){

  /* redefinition from svm.cpp */
  char *svm_type_table[]={"c_svc","nu_svc","one_class","epsilon_svr","nu_svr",NULL};
  char *kernel_type_table[]={"linear","polynomial","rbf","sigmoid",NULL};

  struct svm_model *model;
  char **lines, **fields;
  int i,j,k,l,m;
  char *key, *value, *field;
  char c;
  int dataStart, elements;
  int isColon;
  struct svm_node *x_space=NULL;
  
  model = (struct svm_model*)space(sizeof(struct svm_model));

  model->rho = NULL;
  model->probA = NULL;
  model->probB = NULL;
  model->label = NULL;
  model->nSV = NULL;


  /* Read header until support vectors start */
  lines=splitLines(modelString);
  i=0;
  while (strcmp(lines[i],"SV")!=0){
	fields=splitFields(lines[i]);

	key=fields[0];
	  
	if(strcmp(key,"svm_type")==0){
	  value=fields[1];
	  for(j=0;svm_type_table[j];j++){
		if(strcmp(svm_type_table[j],value)==0){
		  model->param.svm_type=j; 
		  break;
		}
	  }
	  if(svm_type_table[i] == NULL){
		fprintf(stderr,"unknown svm type.\n");
		free(model->rho);
		free(model->label);
		free(model->nSV);
		free(model);
		return NULL;
	  }
	} else
	  
	if(strcmp(key,"kernel_type")==0){
	  value=fields[1];
	  for(j=0;kernel_type_table[j];j++){
		if(strcmp(kernel_type_table[j],value)==0){
		  model->param.kernel_type=j; 
		  break;
		}
	  }
	  if(kernel_type_table[i] == NULL){
		fprintf(stderr,"unknown kernel type.\n");
		free(model->rho);
		free(model->label);
		free(model->nSV);
		free(model);
		return NULL;
	  }
	} else
 
	if (strcmp(key,"gamma")==0){
	  value=fields[1];
	  sscanf(value,"%lf",&model->param.gamma);
	}
	  
	if (strcmp(key,"degree")==0){
	  value=fields[1];
	  sscanf(value,"%lf",&model->param.degree);
	} else 

	if (strcmp(key,"coef0")==0){
	  value=fields[1];
	  sscanf(value,"%lf",&model->param.coef0);
	} else 
	if (strcmp(key,"nr_class")==0){
	  value=fields[1];
	  sscanf(value,"%d",&model->nr_class);
	} else 
	if (strcmp(key,"total_sv")==0){
	  value=fields[1];
	  sscanf(value,"%d",&model->l);
	} else

	if (strcmp(key,"rho")==0){
	  int n = model->nr_class * (model->nr_class-1)/2; 
	  model->rho = (double*)space(sizeof(double)*n);
	  for(j=0;j<n;j++){
		sscanf(fields[j+1],"%lf",&model->rho[j]);
	  }
	} else 

	if (strcmp(key,"nr_sv")==0){
	  int n = model->nr_class;
	  model->nSV = (int*)space(sizeof(int)*n);
	  for(j=0;j<n;j++){
		sscanf(fields[j+1],"%d",&model->nSV[j]);
	  }
	} else

	if (strcmp(key,"label")==0){
	  int n = model->nr_class;
	  model->label = (int*)space(sizeof(int)*n);
	  for(j=0;j<n;j++){
		sscanf(fields[j+1],"%d",&model->label[j]);
	  }
	} else 

	if (strcmp(key,"probA")==0){
	  int n = model->nr_class * (model->nr_class-1)/2;
	  model->probA = (double*)space(sizeof(double)*n);
	  for(j=0;j<n;j++){
		sscanf(fields[j+1],"%lf",&model->probA[j]);
	  }
	} else 

	if (strcmp(key,"probB")==0){
	  int n = model->nr_class * (model->nr_class-1)/2;
	  model->probB = (double*)space(sizeof(double)*n);
	  for(j=0;j<n;j++){
		sscanf(fields[j+1],"%lf",&model->probB[j]);
	  }
	}
	i++;
	freeFields(fields);
  }

  dataStart=i+1;
  elements=0;

  /* Count number of nodes (by counting colons) in advance to allocate
	 memory in one block */
  while (lines[i]!=NULL){
	j=0;
	while ((c=lines[i][j])!='\0'){
	  if (c==':'){
		elements++;
	  }
	  j++;
	}
	elements++;
	i++;
  }

  /* allocate memory for SVs and coefficients */
  m = model->nr_class - 1;
  l = model->l;
  model->sv_coef = (double**)space(sizeof(double*)*m);
  for(i=0;i<m;i++){
	model->sv_coef[i] = (double*)space(sizeof(double)*l);
  }
  model->SV = (struct svm_node**)space(sizeof(struct svm_node*)*l);
	
		
  if(l>0){
	x_space = (struct svm_node*)space(sizeof(struct svm_node)*(elements));
  }
	

  /* parse support vector data */
  
  j=0; 
  for(i=0;i<l;i++){
	fields=splitFields(lines[dataStart+i]);
	model->SV[i] = &x_space[j];
	  
	k=0;
	while ((field=fields[k])!=NULL){
	  if (k<m){
		sscanf(fields[k],"%lf",&model->sv_coef[k][i]);
	  } else {
		sscanf(fields[k],"%d:%lf",&(x_space[j].index),&(x_space[j].value));
		j++;
	  }
	  k++;
	}
	x_space[j++].index = -1;
	freeFields(fields);
  }
  
  freeFields(lines);

  model->free_sv = 1;
  
  return(model);
}


int print_model(const char *model_file_name, struct svm_model *model)
{

  struct svm_node *p;
  int nr_class, i,j, l;
  FILE *fp = stdout;

  char *svm_type_table[] =
	{
	  "c_svc","nu_svc","one_class","epsilon_svr","nu_svr",NULL
	};

  char *kernel_type_table[]=
	{
	  "linear","polynomial","rbf","sigmoid",NULL
	};

  
	if(fp==NULL) return -1;


	fprintf(fp,"svm_type %s\n", svm_type_table[model->param.svm_type]);
	fprintf(fp,"kernel_type %s\n", kernel_type_table[model->param.kernel_type]);

	if(model->param.kernel_type == POLY)
		fprintf(fp,"degree %g\n", model->param.degree);

	if(model->param.kernel_type == POLY || model->param.kernel_type == RBF || model-> param.kernel_type == SIGMOID)
		fprintf(fp,"gamma %g\n", model->param.gamma);

	if(model->param.kernel_type == POLY || model->param.kernel_type == SIGMOID)
		fprintf(fp,"coef0 %g\n", model->param.coef0);

	nr_class = model->nr_class;
	l = model->l;
	fprintf(fp, "nr_class %d\n", nr_class);
	fprintf(fp, "total_sv %d\n",l);
	
	{
		fprintf(fp, "rho");
		for(i=0;i<nr_class*(nr_class-1)/2;i++)
			fprintf(fp," %g",model->rho[i]);
		fprintf(fp, "\n");
	}
	
	if(model->label)
	{
		fprintf(fp, "label");
		for(i=0;i<nr_class;i++)
			fprintf(fp," %d",model->label[i]);
		fprintf(fp, "\n");
	}

	if(model->probA)
	{
		fprintf(fp, "probA");
		for(i=0;i<nr_class*(nr_class-1)/2;i++)
			fprintf(fp," %g",model->probA[i]);
		fprintf(fp, "\n");
	}
	if(model->probB)
	{
		fprintf(fp, "probB");
		for(i=0;i<nr_class*(nr_class-1)/2;i++)
			fprintf(fp," %g",model->probB[i]);
		fprintf(fp, "\n");
	}

	if(model->nSV)
	{
		fprintf(fp, "nr_sv");
		for(i=0;i<nr_class;i++)
			fprintf(fp," %d",model->nSV[i]);
		fprintf(fp, "\n");
	}


	fprintf(fp, "SV\n");

	for(i=0;i<l;i++)
	{
	  for(j=0;j<nr_class-1;j++)
		fprintf(fp, "%.16g ",model->sv_coef[j][i]);
	  
	  p = model->SV[i];
	  while(p->index != -1)
		{
		  fprintf(fp,"%d:%.8g ",p->index,p->value);
		  p++;
		}
	  fprintf(fp, "\n");
	}

	return 0;
}
