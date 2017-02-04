#ifndef __KERNEL_TABLE__
#define __KERNEL_TABLE__

typedef struct {
	double minsegma;
	double maxsegma;
	int bins;
	int tablehalfwidth;
	double * kerneltable;
} KernelTableStruct;

#endif