bin_PROGRAMS = RNAz

RNAz_SOURCES=RNAz.c rnaz_utils.c svm_helper.c rnaz_utils.h svm_helper.h zscore.c zscore.h cmdline.c cmdline.h strand.c strand.h  ../models/mfe_avg.inc  ../models/mfe_stdv.inc  ../models/strand.inc ../models/decision.inc 

nodist_RNAz_SOURCES=../libsvm-2.89/svm.cpp 

RNAz_LDADD = ../lib_phylip/libPHYLIP.a ../lib_rna_features/libRNAFeatures.a ../librna2/libRNA2.a
 
RNAz_LINK = $(CXX) -o $@
 
INCLUDES = -I../lib_phylip -I../lib_rna_features -I../librna2 -I../libsvm-2.89 -I../models

