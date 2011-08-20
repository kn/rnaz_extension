#ifndef LRFUTILS_H
#define LRFUTILS_H
typedef struct {
	int i;
	int j;
} tuple;

typedef struct {
	tuple bp;
	tuple loop;
} bp_info;

extern void *allocate(unsigned int size, char* name);
extern void find_base_pairs(const char* structure, int* numBps, bp_info* bps);
#endif