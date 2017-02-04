#ifndef __ColorModel_H
#define __ColorModel_H

#define MAX_COLORMODEL_SAMPLES 5000
#define UPDATE_PROFILE_RATE 0.9
#define MAX_TARGETS 25000

#include "KernelTable.h"
#include "FastGauss.h"


// computational methods
#define DIRECT_EVALUATION 1
#define FAST_GAUSS_3D 2
#define KERNEL_LUT 3
#define FAST_GAUSS_2D 4

class ColorModel {
 
private :

	unsigned int VecSize;
	unsigned char * Sample;
	unsigned int rows,cols;
	unsigned int init_profile_flag;

	unsigned int comp_method;

	double * DSample;
	double * Sources;
	double * SrcStrength;
	double * segma;
	double * FGresults;
	double * results;
	double * Targets;

	KernelTableStruct * KernelTable;


	FastGauss * pFastGauss;



	double * temp;

	double Profile_h;

public :
	double * YProfile;
	double * Prob_image;

	unsigned int SampleCount;
	unsigned int SubSampleCount;

	ColorModel(unsigned int OneSampleSize,unsigned int Rows,unsigned int Cols);
	~ColorModel();
	void ColorModel::Clear();
	void AddSample(unsigned char * Sample,double * DSampleVector);

	double ColorModel::EstimateProbability(unsigned char * Vec,KernelTableStruct * KernelTable);
	void ColorModel::WriteSamples(char * fname);

	void ColorModel::InitializeFastGaussComputation();

	void ColorModel::InitializeComputation(unsigned int compmethod);


	void ColorModel::SetKernelTable(KernelTableStruct * kerneltableptr){
		KernelTable = kerneltableptr;
	};

	void ColorModel::Kernel_LUTEvaluation(unsigned char * srcs,
											double * str,
											double * segma,
											unsigned char * trgs,
											double * result,
											unsigned int Nsources,
											unsigned int Ntargets,
											unsigned int Dim) ;


	float ColorModel::EstimateProbabilityFG(unsigned int FGPixelCount,
																			 unsigned int * FGPixelLocX,
																			 unsigned int * FGPixelLocY,
																			 unsigned int * FGPixelMap,
																			 double * FGtargets,
																			 unsigned char * FGQtargets,
																			 unsigned int Ntragets,
																			 double * Pimage);


	void ColorModel::ResetProbImage();


	void ColorModel::UpdateProfile(int * cnt,int ymin,int ymax,int top_bottom_sep);
	void SmoothProfile(double sigma,KernelTableStruct * KernelTable);


};


#endif __ColorModel_H
