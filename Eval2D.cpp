
#include "FastGauss.h"
#include <stdlib.h>
#include <math.h>
//#include <windows.h>
#include <stdio.h>
#include <iostream.h>


# define PI 3.141592654

/*--------------------macros-----------------------------*/


// TIMED_FUNCTION macro call the function "fn" and measure
// the time required to execute this function and add that
// to "x"

#define TIMED_FUNCTION(fn,x) \
{			\
	LARGE_INTEGER __t1,__t2,__f;  \
	QueryPerformanceFrequency(&__f);\
	QueryPerformanceCounter(&__t1);\
	(fn);\
	QueryPerformanceCounter(&__t2);	\
	(x)+=(float) (1000.0*(float)(__t2.QuadPart - __t1.QuadPart) / __f.QuadPart); \
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
				x =(srcs[j*Dim+k] - trgs[i*Dim+k]) / segma[j*Dim+k];
				sum+= x*x;
			}
			s=str[j] * exp(-sum);
			result[i]+=s;
			//result[i]+=str[j] * exp(-sum);

		}
}

/* 
  generate random number  0<= r < 1
	r is strictly less than 1
*/
double rand_1 ()
{
	double r;

	do {
		r = (double) rand() / (double) RAND_MAX;
	}
	while (r >= 1.0);

	return r;
}

void GenerateRandomST(unsigned int Nsources,
											unsigned int Ntargets,
											unsigned int Dim,
											double * sources,
											double * segma,
											double * str,
											double * targets, 
											double * fixedsegma)
{
	unsigned int i,j;
	double strsum;

	// Generate random sources
	strsum=0;
	for ( i=0;i<Nsources;i++) {
		for (j=0;j<Dim;j++) {
			sources[i*Dim+j]=rand_1();
			segma[i*Dim+j]=fixedsegma[j];
		}
		str[i]=rand_1();
		strsum+=str[i];
	}


	// normalize sum of str[] to 1
	for ( i=0;i<Nsources;i++) {
		str[i] = str[i] / strsum;
	}
	
	/*
	for ( i=0;i<Nsources;i++) {
		printf("%g %g \n",sources[2*i],sources[2*i+1]);
	}
	*/

	// Generate random targets
	for ( i=0;i<Ntargets;i++) {
		for (j=0;j<Dim;j++) {
			targets[i*Dim+j]=rand_1();
		}
	}

}


void GenerateSourcesinCircle(unsigned int Nsources,
													double * sources,
													double * segma,
													double * str,
													double * fixedsegma,
													double R,
													double & cx, double & cy)

{

	//double cx,cy;
	//double R;
	double r,t,x,y;
	double strsum;
	unsigned int i;

	// random center
	cx=rand_1()* 0.5 + 0.25;
	cy=rand_1()* 0.5 + 0.25;

	R=0.1;

	strsum=0.0;
	// Generate random sources
	for ( i=0;i<Nsources;i++) {
		r=rand_1()*R;
		t=rand_1()*2*PI - PI;;

		x=r * cos(t) +cx;
		y=r * sin(t) +cy;
		
		sources[i*2]=x;
		sources[i*2+1]=y;

		segma[i*2]=fixedsegma[0];
		segma[i*2+1]=fixedsegma[1];

		str[i]=rand_1();
		strsum+=str[i];

	}

		// normalize sum of str[] to 1
	for ( i=0;i<Nsources;i++) {
		str[i] = str[i] / strsum;
	}


}


void GenerateSourcesinCircle3D(unsigned int Nsources,
													double * sources,
													double * segma,
													double * str,
													double * fixedsegma,
													double R, 
													double & cx, double & cy, double & cz)
													
{

	//double cx,cy,cz;
	double r,t1,t2,x,y,z;
	double strsum;
	unsigned int i;

	// random center
	cx=rand_1()* 0.5 + 0.25;
	cy=rand_1()* 0.5 + 0.25;
	cz=rand_1()* 0.5 + 0.25;

	strsum=0;

	// Generate random sources
	for ( i=0;i<Nsources;i++) {
		r=rand_1()*R;
		t1=rand_1()*2*PI - PI;;
		t2=rand_1()*2*PI - PI;;


		x=r * cos(t2) * cos(t1) +cx;
		y=r * cos(t2) * sin(t1) +cy;
		z=r * sin(t2) + cz;

		sources[i*3]=x;
		sources[i*3+1]=y;
		sources[i*3+2]=z;

		segma[i*3]=fixedsegma[0];
		segma[i*3+1]=fixedsegma[1];
		segma[i*3+2]=fixedsegma[2];

		str[i]=rand_1();
		strsum+=str[i];

	}

	// normalize sum of str[] to 1
	for ( i=0;i<Nsources;i++) {
		str[i] = str[i] / strsum;
	}



}

void GenerateTargetsinCircle2D(unsigned int Ntargets,
													double * targets,
													double R,
													double cx, double cy)

{
	double r,t,x,y;
	unsigned int i;

	// random center
	//cx=rand_1()* 0.5 + 0.25;
	//cy=rand_1()* 0.5 + 0.25;

	//R=0.1;

	// Generate random targets
	for ( i=0;i<Ntargets;i++) {
		r=rand_1()*R;
		t=rand_1()*2*PI - PI;;

		x=r * cos(t) +cx;
		y=r * sin(t) +cy;
		
		targets[i*2]=x;
		targets[i*2+1]=y;

	}



}

void GenerateTargetsinCircle3D(unsigned int Ntargets,
															 double * targets,
															 double R,
															 double cx, double cy, double cz)

{
	//double cx,cy,cz;
	double r,t1,t2,x,y,z;
	unsigned int i;

	// random center
	//cx=rand_1()* 0.5 + 0.25;
	//cy=rand_1()* 0.5 + 0.25;
	//cz=rand_1()* 0.5 + 0.25;


	// Generate random targets
	for ( i=0;i<Ntargets;i++) {
		r=rand_1()*R;
		t1=rand_1()*2*PI - PI;;
		t2=rand_1()*2*PI - PI;;

		x=r * cos(t2) * cos(t1) +cx;
		y=r * cos(t2) * sin(t1) +cy;
		z=r * sin(t2) + cz;
		
		targets[i*3]=x;
		targets[i*3+1]=y;
		targets[i*3+2]=z;

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


//int main()

int main(int argc, char * argv[])
{

	unsigned int Nsources;
	unsigned int Ntargets;
	unsigned int Dim=3;
	unsigned int i,j;


	//Nsources=atoi(argv[1]);

	//Ntargets=atoi(argv[1]);
	//Nsources=1200;
	
	//cout << "Enter number of source & targets :";
	//cin >> Nsources;

	/*
	cout << "Enter number of source :";
	cin >> Nsources;
	
	cout << "Enter number of targets :";
	cin >> Ntargets;
	*/


	// read sources and targets from files
	char fname[150];
	FILE * FP;
	unsigned int nsourcesread,ntargetsread;

	double * sources_readbuffer=new double[10000*Dim];
	double * targets_readbuffer=new double[10000*Dim];


	FP=fopen("C:\\users\\Ahmed\\Research\\code\\FastGauss\\MultiDScale2Phase\\Debug\\Obj1_Head.txt","rb");

	nsourcesread=fread(sources_readbuffer,sizeof(double)*Dim,10000,FP);

	fclose(FP);


	FP=fopen("C:\\users\\Ahmed\\Research\\code\\FastGauss\\MultiDScale2Phase\\Debug\\targets89.txt","rb");

	ntargetsread=fread(targets_readbuffer,sizeof(double)*Dim,10000,FP);

	fclose(FP);


	Nsources=nsourcesread;
	Ntargets=ntargetsread;

	


	cout << "Number of source  = " << Nsources << endl;
	
	cout << "Number of targets = " << Ntargets << endl;


	//Ntargets=Nsources;
	
	


	double * fixedsegma = new double[Dim+1];
	
	float t1=0.0;
	float t2=0.0;
	float t3=0.0;

	fixedsegma[0]=0.05;
	fixedsegma[1]=0.01;
	fixedsegma[2]=0.01;

	double tolerance=1e-6;

	FastGauss * pFastGauss = new FastGauss(Dim);

	// allocate memory for sources and targets
	double * sources = new double [Nsources*Dim];
	double * str = new double [Nsources];
	double * segma = new double [Nsources*Dim];
	double * targets = new double [Ntargets*Dim];
	double * results = new double [Ntargets];
	double * FFresults;
	// double * FFresults =new double [Ntargets];
 


	// fill in the data
	
	//GenerateRandomST(Nsources,Ntargets,Dim,sources,segma,str,targets,fixedsegma);

	
	{
		memcpy(sources,sources_readbuffer,Nsources*Dim*sizeof(double));

		memcpy(targets,targets_readbuffer,Ntargets*Dim*sizeof(double));

		unsigned int i,j;

		for ( i=0;i<Nsources;i++) {
			str[i] = 1.0 / (double) Nsources;
			for (j=0;j<Dim;j++) {
				segma[i*Dim+j]=fixedsegma[j];
			}

		}
	}
	

	double cx,cy,cz;
	
	
	
	/*
	if (Dim==2)
		GenerateSourcesinCircle(Nsources,sources,segma,str,fixedsegma,0.1,cx,cy);
	else if (Dim==3)
		GenerateSourcesinCircle3D(Nsources,sources,segma,str,fixedsegma,0.1,cx,cy,cz);
	
 
	
	if (Dim==2)
		GenerateTargetsinCircle2D(Ntargets,targets,0.1,cx,cy);
	else if (Dim==3)
	  GenerateTargetsinCircle3D(Ntargets,targets,0.1,cx,cy,cz);
	
	
	*/

	
	/*
	
	// Generate random targets
	for ( i=0;i<Ntargets;i++) {
		for (j=0;j<Dim;j++) {
			targets[i*Dim+j]=rand_1();
		}
	}
	
	*/
	

	/*
	sources[0]=0.5;
	sources[1]=0.5;

	targets[0]=0.4;
	targets[1]=0.4;

	str[0]=1.0;
	*/

	/*
	printf("\n");

	for ( i=0;i<Ntargets;i++) {
		printf("%g %g \n",targets[2*i],targets[2*i+1]);
	}
	*/


	TIMED_FUNCTION(DirectEvaluation(sources,str,segma,targets,results,Nsources,Ntargets,Dim),t1);


	cout << "N sources = " << Nsources << endl;
	cout << "N targets = " << Ntargets << endl;
	cout << "Dim	= " << Dim << endl;
	cout << "segma= " << fixedsegma[0]	<< " , " 
										<< fixedsegma[1]	<< " , " 
										<< fixedsegma[2]	<< " , " 
										<< endl;

	FFresults=NULL;

	//TIMED_FUNCTION(pFastGauss->FixedScaleEvaluation(sources,str,segma,targets,FFresults,Nsources,Ntargets,tolerance),t2);
	//printf("FG = %g msec \n",t);

	
	TIMED_FUNCTION(pFastGauss->BuildSourceStructure(Nsources,sources,str,segma,tolerance),t2);

	TIMED_FUNCTION(pFastGauss->EvaluateTargets(Ntargets,targets,FFresults),t3);


	// compare results
	double sumerr,maxerr,minerr;
	sumerr=eval_error(results,FFresults,Ntargets,maxerr,minerr);


	cout << "avg error = " << (sumerr / Ntargets) << endl;
	cout << "max error = " << maxerr << endl;
	cout << "min error = " << minerr << endl;

	cout << " Direct evaluation     = " << t1 << " msec " << endl;
	cout << " Fast Gauss evaluation (sources) = " << t2 << " msec " << endl;
	cout << " Fast Gauss evaluation (targets) = " << t3 << " msec " << endl;
	cout << " Fast Gauss evaluation  " << t2+t3 << " msec " << endl;


	return 1;
}

