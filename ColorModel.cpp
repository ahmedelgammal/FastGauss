#include "ColorModel.h"
#include <math.h>
#include <assert.h>
#include <stdio.h>

#define PI 3.14159

ColorModel::ColorModel(unsigned int OneSampleSize,unsigned int R,unsigned int C)
{
	rows=R;
	cols=C;


	Sample = new unsigned char [MAX_COLORMODEL_SAMPLES*OneSampleSize];
	DSample = new double [MAX_COLORMODEL_SAMPLES*OneSampleSize];
	//Sources = new double [MAX_COLORMODEL_SAMPLES*OneSampleSize];
	
	Sources = new double [MAX_COLORMODEL_SAMPLES*2];
	
	SrcStrength= new double [MAX_COLORMODEL_SAMPLES];
	segma= new double [OneSampleSize];

	FGresults= new double [MAX_TARGETS];
	results= new double [MAX_TARGETS];

	Targets= new double[MAX_TARGETS*2];

	SampleCount=0;
	SubSampleCount=0;
	VecSize=OneSampleSize;

	YProfile = new double [rows*2];
	temp= new double [rows*2];

	for( int y=0;y<rows*2;y++)
		YProfile[y]=0;

	init_profile_flag=1;

}


ColorModel::~ColorModel()
{
	
	delete Sample;
}

void ColorModel::Clear()
{

	SampleCount=0;

}



void ColorModel::AddSample(unsigned char * SampleVector, double * DSampleVector)
{
	for (unsigned int i=0;i<VecSize;i++) {
		Sample[3*SampleCount+i]=SampleVector[i];
		DSample[3*SampleCount+i]=DSampleVector[i];
	}

	SampleCount= (SampleCount  +1) % MAX_COLORMODEL_SAMPLES;
}


double ColorModel::EstimateProbability(unsigned char * Vec,KernelTableStruct * KernelTable)
{
	unsigned int i;
	int d;
	double sum,Pr;
	double k1,k2,k3;
	int idx;
	
	//
	double segmamin=KernelTable->minsegma;
	double segmamax=KernelTable->maxsegma;
	int segmabins=KernelTable->bins;
	int kerneltablewidth = 2*KernelTable->tablehalfwidth+1;
	
	double * kerneltable = KernelTable->kerneltable;

	double segma1=4.0;
	double segma2=1.0;
	double segma3=1.0;

	int segma_bin1;
	int segma_bin2;
	int segma_bin3;

	int kernelbase1;
	int kernelbase2;
	int kernelbase3;



	segma_bin1=(int) floor((segma1-segmamin)/(segmamax-segmamin)*segmabins +1 );
	segma_bin2=(int) floor((segma2-segmamin)/(segmamax-segmamin)*segmabins +1 );
	segma_bin3=(int) floor((segma3-segmamin)/(segmamax-segmamin)*segmabins +1 );

	kernelbase1=segma_bin1*kerneltablewidth;
	kernelbase2=segma_bin2*kerneltablewidth;
	kernelbase3=segma_bin3*kerneltablewidth;


	sum=0;
	for (i=0;i< SampleCount;i++) {


		idx=i*VecSize;

		d=(int) Vec[0]- (int) Sample[idx]+255;
		//k1=KernelTable[d];
		k1=kerneltable[kernelbase1+d];


		d=(int) Vec[1]-(int) Sample[idx+1]+255;
		k2=kerneltable[kernelbase2+d];

		d=(int) Vec[2]-(int) Sample[idx+2]+255;
		k3=kerneltable[kernelbase3+d];

		sum+=k1*k2*k3;
		//sum+=k2*k3;

	}

	Pr=(SampleCount > 0 ? sum/SampleCount : 0);

	return Pr;
}


void ColorModel::InitializeFastGaussComputation()
{
	unsigned int i;
	unsigned int r;


	SubSampleCount=200;


	pFastGauss = new FastGauss(3);



	segma[0]=0.05;
	segma[1]=0.01;
	segma[2]=0.01;

	// subsample
	if (SampleCount > SubSampleCount) 
		for(i=0;i<SubSampleCount;i++) {
			r = (unsigned int) (((double) rand() / (double) RAND_MAX) * SampleCount);

			Sources[3*i]  =DSample[3*r];
			Sources[3*i+1]=DSample[3*r+1];
			Sources[3*i+2]=DSample[3*r+2];
			
		}
	

		else {
			SubSampleCount=SampleCount;
			for(i=0;i<SubSampleCount;i++) {
				Sources[3*i]  =DSample[3*i];
				Sources[3*i+1]=DSample[3*i+1];
				Sources[3*i+2]=DSample[3*i+2];				
			}

		}


	double srcstrength = 1.0 / (double) SubSampleCount;


	for (i=0;i< SubSampleCount;i++) {
		SrcStrength[i] = srcstrength;	
	}



	pFastGauss->BuildSourceStructure(SubSampleCount,Sources,SrcStrength,segma,1e-2);

	pFastGauss->PreAllocateTargets(MAX_TARGETS);



}

void ColorModel::InitializeComputation(unsigned int compmethod)
{
	unsigned int i;
	unsigned int r;


	comp_method=compmethod;

	SubSampleCount=200;

	segma[0]=0.05;
	segma[1]=0.01;
	segma[2]=0.01;

	if (comp_method==FAST_GAUSS_2D) {
		pFastGauss = new FastGauss(2);
		segma[0]=0.01;
		segma[1]=0.01;	
	}
	else if (comp_method==FAST_GAUSS_3D)
		pFastGauss = new FastGauss(3);




	// subsample
	if (SampleCount > SubSampleCount) 
		for(i=0;i<SubSampleCount;i++) {
			r = (unsigned int) (((double) rand() / (double) RAND_MAX) * SampleCount);

			if (comp_method==FAST_GAUSS_2D) {
				Sources[2*i]  =DSample[3*r+1];
				Sources[2*i+1]=DSample[3*r+2];
			}				
			else {
				Sources[3*i]  =DSample[3*r];
				Sources[3*i+1]=DSample[3*r+1];
				Sources[3*i+2]=DSample[3*r+2];
			}

			
		}
	

		else {
			SubSampleCount=SampleCount;
			for(i=0;i<SubSampleCount;i++) 
				if (comp_method==FAST_GAUSS_2D) {
					Sources[2*i]  =DSample[3*i+1];
					Sources[2*i+1]=DSample[3*i+2];
				}				
				else {
					Sources[3*i]  =DSample[3*i];
					Sources[3*i+1]=DSample[3*i+1];
					Sources[3*i+2]=DSample[3*i+2];
				}

		}


	double srcstrength = 1.0 / (double) SubSampleCount;


	for (i=0;i< SubSampleCount;i++) {
		SrcStrength[i] = srcstrength;	
	}



	if (comp_method==FAST_GAUSS_2D ||
		comp_method==FAST_GAUSS_3D) {

		pFastGauss->BuildSourceStructure(SubSampleCount,Sources,SrcStrength,segma,1e-2);

		pFastGauss->PreAllocateTargets(MAX_TARGETS);
	}
}


void ColorModel::WriteSamples(char * fname)
{
	FILE * FP;

	FP=fopen(fname,"wb");

	fwrite(Sources,sizeof(double),SampleCount*VecSize,FP);

	fclose(FP);

}

void ColorModel::Kernel_LUTEvaluation(unsigned char * srcs,
											double * str,
											double * segma,
											unsigned char * trgs,
											double * result,
											unsigned int Nsources,
											unsigned int Ntargets,
											unsigned int Dim) 
{

	unsigned int i,j;
	int d;
	double sum;
	double k1,k2,k3;
	int srcidx,trgidx;

	double segmamin=KernelTable->minsegma;
	double segmamax=KernelTable->maxsegma;
	int segmabins=KernelTable->bins;
	int kerneltablewidth = 2*KernelTable->tablehalfwidth+1;
	
	double * kerneltable = KernelTable->kerneltable;

	double segma1=segma[0]*KernelTable->tablehalfwidth/sqrt(2);
	double segma2=segma[1]*KernelTable->tablehalfwidth/sqrt(2);
	double segma3=segma[2]*KernelTable->tablehalfwidth/sqrt(2);

	int segma_bin1;
	int segma_bin2;
	int segma_bin3;

	int kernelbase1;
	int kernelbase2;
	int kernelbase3;

	segma_bin1=(int) floor((segma1-segmamin)/(segmamax-segmamin)*segmabins +1 );
	segma_bin2=(int) floor((segma2-segmamin)/(segmamax-segmamin)*segmabins +1 );
	segma_bin3=(int) floor((segma3-segmamin)/(segmamax-segmamin)*segmabins +1 );

	kernelbase1=segma_bin1*kerneltablewidth;
	kernelbase2=segma_bin2*kerneltablewidth;
	kernelbase3=segma_bin3*kerneltablewidth;


	for (i=0;i<Ntargets;i++) {
		result[i]=0.0;
	}

	double srcstr=1.0 / (double) Nsources;


	for (i=0;i<Ntargets;i++) {
		trgidx=i*Dim;
		sum=0;
		for (j=0;j<Nsources;j++) {
			srcidx=j*Dim;
			d=(int) srcs[srcidx]- (int) trgs[trgidx]+255;
			k1=kerneltable[kernelbase1+d];

			d=(int) srcs[srcidx+1]- (int) trgs[trgidx+1]+255;
			k2=kerneltable[kernelbase2+d];

			d=(int) srcs[srcidx+2]- (int) trgs[trgidx+2]+255;
			k3=kerneltable[kernelbase3+d];


			sum+=srcstr*k1*k2*k3;

		}

		result[i]=sum;
	}

}

void DirectEvaluation(double * srcs,
											double * str,
											double * segma,
											double * trgs,
											double * result,
											unsigned int Nsources,
											unsigned int Ntargets,
											unsigned int Dim) 
{

	unsigned int i,j,k;
	double x,sum,s;

	for (i=0;i<Ntargets;i++) {
		result[i]=0.0;
	}

	for (i=0;i<Ntargets;i++)
		for (j=0;j<Nsources;j++) {
			sum=0;
			for (k=0;k<Dim;k++) {
				x =(srcs[j*Dim+k] - trgs[i*Dim+k]) / segma[k];
				sum+= x*x;
			}
			s=str[j] * exp(-sum);
			result[i]+=s;
			//result[i]+=str[j] * exp(-sum);

		}
}


void DirectEvaluationI(double * srcs,
											double * str,
											double * segma,
											double * trgs,
											double * result,
											unsigned int Nsources,
											unsigned int Ntargets,
											unsigned int Dim) 
{

	unsigned int i,j;
	double x,x1,x2,sum1,sum2;

	for (i=0;i<Ntargets;i++) {
		result[i]=0.0;
	}

	for (i=0;i<Ntargets;i++) {
		sum1=0;
		sum2=0;
		for (j=0;j<Nsources;j++) {
				
			// s
			x =(srcs[j*Dim] - trgs[i*Dim]) / segma[0];
			sum1+=str[j] * exp(-x*x);

			// rg
			x1=(srcs[j*Dim+1] - trgs[i*Dim+1]) / segma[1];
		
			x2=(srcs[j*Dim+2] - trgs[i*Dim+2]) / segma[2];

			sum2+= str[j] * exp (-x1*x1-x2*x2);


		}
		result[i]+=sum1*sum2;
		
	}
}


double eval_error(double * result1,
									double * result2,
									unsigned int N,
									double & maxerr,
									double & minerr)
{
	unsigned int i;
	double err=0,e;
		
	maxerr=0;
	minerr=100;

	for (i=0;i<N;i++) {
		e=fabs(result1[i]-result2[i]);
		err+=e;
		maxerr=(e>maxerr) ? e : maxerr ;
		minerr=(e<minerr) ? e : minerr ;
	}

	return err;
}



float ColorModel::EstimateProbabilityFG(unsigned int FGPixelCount,
																			 unsigned int * FGPixelLocX,
																			 unsigned int * FGPixelLocY,
																			 unsigned int * FGPixelMap,
																			 double * FGtargets,
																			 unsigned char * FGQtargets,
																			 unsigned int Ntargets,
																			 double * Pimage)
{

	unsigned int i;
	unsigned int x;
	unsigned int y;

	
	LARGE_INTEGER pht1,pht2,phf;
	float t1;



	QueryPerformanceCounter(&pht1);


	if (comp_method==DIRECT_EVALUATION)
		DirectEvaluationI(Sources,SrcStrength,segma,FGtargets,FGresults,SubSampleCount,Ntargets,3);
	else if (comp_method==FAST_GAUSS_3D)
		pFastGauss->EvaluateTargets(Ntargets,FGtargets,FGresults);
	else if (comp_method==KERNEL_LUT)
		Kernel_LUTEvaluation(Sample,SrcStrength,segma,FGQtargets,FGresults,SampleCount,Ntargets,3);
	else if (comp_method==FAST_GAUSS_2D) {

		for (i=0;i<Ntargets;i++) {
			Targets[i*2]=FGtargets[i*3+1];
			Targets[i*2+1]=FGtargets[i*3+2];
		}

		pFastGauss->EvaluateTargets(Ntargets,Targets,FGresults);
	
	}


	QueryPerformanceCounter(&pht2);
	QueryPerformanceFrequency(&phf);
	t1=(float) (1000.0*(float)(pht2.QuadPart - pht1.QuadPart) / phf.QuadPart); 

	
	// compare results
	//double sumerr,maxerr,minerr;
	//sumerr=eval_error(results,FGresults,Ntargets,maxerr,minerr);

	
	// fill in Pr. image



	for (i=0; i< FGPixelCount;i++) {
		x=FGPixelLocX[i];
		y=FGPixelLocY[i];

		Pimage[y*cols+x]=FGresults[FGPixelMap[i]];

	}

	return t1;

}


void ColorModel::ResetProbImage()
{

	memset(Prob_image,0,rows*cols*sizeof(double));

}

void ColorModel::UpdateProfile(int * cnt,int ymin,int ymax,int top_bottom_sep)
{
	int sum;
	int y,y1;
	double a=UPDATE_PROFILE_RATE;

	sum=0;
	for(y=ymin;y<=ymax;y++)
		sum+=cnt[y];


	//for (y=0;y<= (ymax-ymin);y++)
	//	YProfile[y]=a*YProfile[y]+(1-a)*((double) cnt[ymax-y]/(double) sum);

	if (init_profile_flag){
		for(y=ymin;y<=ymax;y++) {
			y1=rows+top_bottom_sep-y;
			YProfile[y1]=((double) cnt[y]/(double) sum);
		}
		Profile_h=ymax-ymin+1;
		init_profile_flag=0;
	}
	else
		for(y=ymin;y<=ymax;y++) {
			y1=rows+(int) ((top_bottom_sep-y) * (Profile_h / (double) (ymax-ymin+1)) + 0.5);
			YProfile[y1]=a*YProfile[y1]+(1-a)*((double) cnt[y]/(double) sum);
		}


}

void ColorModel::SmoothProfile(double sigma,KernelTableStruct * KernelTable)
{

	double segmamin=KernelTable->minsegma;
	double segmamax=KernelTable->maxsegma;
	int segmabins=KernelTable->bins;
	int kerneltablewidth = 2*KernelTable->tablehalfwidth+1;
	
	double * kerneltable = KernelTable->kerneltable;
	int y,y1,d;
	double k,sum;

	assert(sigma >=segmamin && sigma <= segmamax);

	int segma_bin=(int) floor((sigma-segmamin)/(segmamax-segmamin)*segmabins +1 );
	int kernelbase=segma_bin*kerneltablewidth;


	for (y=0;y<2*rows;y++) 
		temp[y]=0.0;


	sum=0;
	for (y=0;y<2*rows;y++) {
		if (YProfile[y] > 0.0) {
			for (y1=0;y1<2*rows;y1++) {

				d=(y1-y) / 2 +255;
				
				k= kerneltable[kernelbase+d];
				temp[y1] += YProfile[y]*k;
				sum += YProfile[y]*k;
			}
		}

	
	}

	for (y=0;y<2*rows;y++) 
		YProfile[y]=temp[y]/sum;

}
