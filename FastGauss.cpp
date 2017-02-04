#include "FastGauss.h"
#include <stdlib.h>
#include <math.h>
#include <iostream.h>
#include <assert.h>

FastGauss::FastGauss(unsigned int D)
{

	Dim=D;
	PreAllocateTargetFlag=false;

}



int FastGauss::Allocate()
{


	/*
	sources = new double [Nsources*Dim];
	str = new double [Nsources];
	segma = new double [Nsources*Dim];
	targets = new double [Ntargets*Dim];
	*/

	// source boxes
	SrcsInBoxCount = new int [NBoxes];
	SrcsInBox = new int [NBoxes];
	SrcNext = new int [Nsources];

	BoxCenter = new double [NBoxes*Dim];

	Box2Subscript = new int [NBoxes*Dim];


	SrcBoxFFexp =new double * [NBoxes];

	// target boxes
	TrgsInBoxCount = new int [NBoxes];
	TrgsInBox = new int [NBoxes];
	UsedTrgsBox = new int[NBoxes];

	TrgBoxLocalExp = new double * [NBoxes];
	TrgBoxLocalExpFlag = new unsigned char [NBoxes];

	NbrList = new int [NBoxes];


	//TrgNext = new int [Ntargets];


	// Nbr boxes
	NbrDisp = new int [Ndisp*Dim];
	NbrRelativeIndex=  new int [Ndisp];

	// expansion work space;
	//LocalExp = new double [Nterms2pDim];
	
	FFexp = new double [Nterms2pDim*MAXFFEXP];
	

	// temp vectors
	tmp=new int [Dim];
	Dtmp= new double[Dim];

	/*
	return (sources != (double * ) NULL)  && 
		(str != (double *) NULL) &&
		(segma != (double *) NULL) &&
		(targets != (double *) NULL) ;
	*/

	return 1;

}

int FastGauss::AllocateTargets()
{

	// result buffer 
	if (TargetField== NULL)
		TargetField = new double [Ntargets];
	TrgNext = new int [Ntargets];

	return 1;
}

int FastGauss::PreAllocateTargets(int MaxTargets)
{

	TrgNext = new int [MaxTargets];
	PreAllocateTargetFlag = true;

	return 1;
}



#define SRC(i,j) sources[(i)*Dim+(j)]
#define TRG(i,j) targets[(i)*Dim+(j)]

#define CENT(i,j)  BoxCenter[(i)*Dim+(j)]
//#define SEGMA(i,j)	segma[(i)*Dim+(j)]

void FastGauss::GenerateBoxes()
{
	// find boxes center

	unsigned int i,j;
	double * ptr;
	int * iptr;
	int f;
	int d;

	double * h = new double [Dim];
	double * hh = new double [Dim];

	for (i=0;i<Dim;i++){
		h[i] = 1.0 / Nboxes_eachside[i];
		hh[i] = 0.5 * h[i];
	}

	ptr=BoxCenter;
	iptr=Box2Subscript;

	for (i=0;i<NBoxes;i++) {
		f=i;
		for (j=0;j<Dim;j++) {
			d = f % Nboxes_eachside[j];
			f = f / Nboxes_eachside[j];

			// SrcBoxCenter(i,j)=d*h+hh;
			*ptr++=d*h[j]+hh[j];
			*iptr++=d;
		}

	}

	delete [] h;
	delete [] hh;

}

void FastGauss::AssignSources()
{

	unsigned int i;
	int j;
	unsigned int box_index;
	int * dim_index = new int [Dim];
	double * h = new double [Dim];

	for (i=0;i<Dim;i++){
		h[i] = 1.0 / Nboxes_eachside[i];
	}


	// initialize source boxes' count 
	// initialize Sources in box link
	for (i=0;i<NBoxes;i++) {
		SrcsInBoxCount[i]=0;
		SrcsInBox[i]=-1;
	}


	// assign sources to boxes
	for (i=0;i<Nsources;i++) {
		// find the box that contains source i 
		box_index=0;
		for (j=Dim-1;j>=0;j--) {
			dim_index[j]=(int) (SRC(i,j) / h[j]);
			//box_index+=dim_index[j]*Nboxes_powers[j];
			box_index = box_index * Nboxes_eachside[j] + dim_index[j];
		}

		SrcsInBoxCount[box_index]++;

		SrcNext[i]=SrcsInBox[box_index];
		SrcsInBox[box_index]=i;

	}


	delete [] dim_index;
	delete [] h;
}


void FastGauss::AssignTargets()
{

	unsigned int i;
	int j;
	unsigned int box_index;
	//int * dim_index = new int [Dim];
	int dim_index[MAXDIM];
	//double * h = new double [Dim];
	double h[MAXDIM];

	for (i=0;i<Dim;i++){
		h[i] = 1.0 / Nboxes_eachside[i];
	}

	// initialize target boxes' count 
	// initialize targets in box link
	for (i=0;i<NBoxes;i++) {
		TrgsInBoxCount[i]=0;
		TrgsInBox[i]=-1;
	}

	UsedTrgsBoxCount=0;

	for (i=0;i<Ntargets;i++) {
		// find the box that contains source i 
		box_index=0;
		for (j=Dim-1;j>=0;j--) {
			dim_index[j]=(int) (TRG(i,j) / h[j]);
			//box_index+=dim_index[j]*Nboxes_powers[j];
			box_index = box_index * Nboxes_eachside[j] + dim_index[j];
		}

		if (TrgsInBoxCount[box_index] == 0) {			
			UsedTrgsBox[UsedTrgsBoxCount]=box_index;
			UsedTrgsBoxCount++;
		}

		TrgsInBoxCount[box_index]++;

		TrgNext[i]=TrgsInBox[box_index];
		TrgsInBox[box_index]=i;


	}

	//delete [] dim_index;
	//delete [] h;

}

void FastGauss::GenerateNbrs(int boxidx,int * nbrlist,unsigned int & NbrCount)
{

	int i;
	unsigned int j;
	int outofboundryflag;
	int d;

	NbrCount=0;

	/*
	// factorize box index
	int f=boxidx;
	
	for (j=0;j<Dim;j++) {
		tmp[j] = f % Nboxes_eachside[j];
		f = f / Nboxes_eachside[j];
	}
	*/

	int * boxsubscript=Box2Subscript+boxidx*Dim;

	// tmp[j] is box in dim j

	for (i=0;i<Ndisp;i++) {
		
		// check inclusion of relative nbr [i]
		outofboundryflag=0;
		for (j=0;j<Dim;j++) {
			//d=tmp[j]+NbrDisp[i*Dim+j];
			d=boxsubscript[j]+NbrDisp[i*Dim+j];
			if ( d < 0 || d >= Nboxes_eachside[j]) {
				outofboundryflag=1;
				break;
			}
		}

		if (! outofboundryflag) {
			nbrlist[NbrCount]=NbrRelativeIndex[i]+boxidx;
			NbrCount++;			
		}
	}



}

void FastGauss::GenerateDisplacements()
{
	int i;
	int j;
	int k,w,f,d;
	int idx;

	k=NbrSize;
	w=2*k+1;
	//Ndisp=pow(w,Dim);


	int * NbrDispptr=NbrDisp;

	for (i=0;i< Ndisp;i++) {
		f=i;
		idx=0;
		for (j=Dim-1;j>=0;j--) {
			d = f % w - k;
			f = f / w;
	
			//idx+=d*Nboxes_powers[j];
			idx= idx * Nboxes_eachside[j] + d;
			// NbrsDim(i,j)=d 
			//*NbrDispptr++=d;
			NbrDisp[i*Dim+j]=d;
		}

		NbrRelativeIndex[i]=idx;

	}

}

/*
	Evaluate all targets at box (trgbox) at all sources at srcbox
*/
void FastGauss::DirectEvaluation(int srcbox,int trgbox)
{
	int src,trg;
	unsigned int j;
	double d,dd;
	double sum=0;
	
	trg=TrgsInBox[trgbox];
	while (trg >= 0) {
		sum=0;

		src=SrcsInBox[srcbox];
		while (src >= 0 ) {

			dd=0;
			for (j=0;j<Dim;j++) {
				d=(SRC(src,j)-TRG(trg,j)) / fixedsegma[j];
				dd+=d*d;
			}

			sum+=exp(- dd) * str[src];
			src=SrcNext[src];
		}

		TargetField[trg]+=sum;
		trg=TrgNext[trg];
	}


}


void	TabulateHermitPoly1D(unsigned int nterms,
													 double x,
													 double * hermitX)
{

		double x2,facx;
		unsigned int i;
		int i2;

		x2=2.0 * x;
		
		facx = exp(-x*x);

		hermitX[0]=facx;

		hermitX[1]=x2*facx;

		for (i=1;i<nterms;i++) {
			i2=2*i;
			hermitX[i+1] = x2 * hermitX[i] - i2 * hermitX[i-1];
		}

		

}

void	TabulateHermitPoly2D(unsigned int nterms,
													 double x,
													 double y,
													 double * hermitX,
													 double * hermitY)
{

		double x2;
		double y2;
		double facx,facy;
		unsigned int i;
		int i2;

		x2=2.0 * x;
		y2=2.0 * y;

		facx = exp(-x*x);
		facy = exp(-y*y);

		hermitX[0]=facx;
		hermitY[0]=facy;

		hermitX[1]=x2*facx;
		hermitY[1]=y2*facy;

		for (i=1;i<nterms;i++) {
			i2=2*i;

			hermitX[i+1] = x2 * hermitX[i] - i2 * hermitX[i-1];
			hermitY[i+1] = y2 * hermitY[i] - i2 * hermitY[i-1];

		}

		

}

void	TabulateHermitPoly3D(unsigned int nterms,
													 double x,
													 double y,
													 double z,
													 double * hermitX,
													 double * hermitY,
													 double * hermitZ)
{

		double x2;
		double y2;
		double z2;
		double facx,facy,facz;
		unsigned int i;
		int i2;

		x2=2.0 * x;
		y2=2.0 * y;
		z2=2.0 * z;


		facx = exp(-x*x);
		facy = exp(-y*y);
		facz = exp(-z*z);


		hermitX[0]=facx;
		hermitY[0]=facy;
		hermitZ[0]=facz;

		hermitX[1]=x2*facx;
		hermitY[1]=y2*facy;
		hermitZ[1]=z2*facz;

		for (i=1;i<nterms;i++) {
			i2=2*i;

			hermitX[i+1] = x2 * hermitX[i] - i2 * hermitX[i-1];
			hermitY[i+1] = y2 * hermitY[i] - i2 * hermitY[i-1];
			hermitZ[i+1] = z2 * hermitZ[i] - i2 * hermitZ[i-1];

		}

		

}

/*
	TaylorExpansion2D :
	2D taylor expansion of source at loc src[] with strenght str
	about location cent[].
	results is added to local expansion: locexp
*/

void TaylorExpansion2D(double * src,
											 double srcstr,
											 double * cent,
											 double * locexp,
											 unsigned int Nterms,
											 double * segma)	
{


	double hermitX[MAXTERMS+1], hermitY[MAXTERMS+1];
	double x,y,fac;
	unsigned int i,j;
	double * expptr;
	double hx;

	x=(cent[0]-src[0]) / segma[0];
	y=(cent[1]-src[1]) / segma[1];

	TabulateHermitPoly2D(Nterms,x,y,hermitX,hermitY);

	fac=1.0;
	for (i=0;i<Nterms;i++) {
		hermitX[i] *= fac;
		hermitY[i] *= fac;
		fac = - fac / double (i+1);
		
	}

	expptr = locexp;
	for (i=0;i<Nterms; i++) {
		hx= hermitX[i]*srcstr;
		for (j=0;j<Nterms;j++) {
			*expptr++ += hx * hermitY[j];
		}
	}
}

/*
	TaylorExpansion3D :
	3D taylor expansion of source at loc src[] with strenght str
	about location cent[].
	results is added to local expansion: locexp
*/

void TaylorExpansion3D(double * src,
											 double srcstr,
											 double * cent,
											 double * locexp,
											 unsigned int Nterms,
											 double * segma)	
{


	double hermitX[MAXTERMS+1], 
				hermitY[MAXTERMS+1],
				hermitZ[MAXTERMS+1];

	double x,y,z,fac;
	unsigned int i,j,k;
	double * expptr;
	double hx,hxy;

	x=(cent[0]-src[0]) / segma[0];
	y=(cent[1]-src[1]) / segma[1];
	z=(cent[2]-src[2]) / segma[2];

	//TabulateHermitPoly2D(Nterms,x,y,hermitX,hermitY);
	//TabulateHermitPoly1D(Nterms,x,hermitX);
	//TabulateHermitPoly1D(Nterms,y,hermitY);
	//TabulateHermitPoly1D(Nterms,z,hermitZ);
	TabulateHermitPoly3D(Nterms,x,y,z,hermitX,hermitY,hermitZ);


	fac=1.0;
	for (i=0;i<Nterms;i++) {
		hermitX[i] *= fac;
		hermitY[i] *= fac;
		hermitZ[i] *= fac;
		fac = - fac / double (i+1);
		
	}

	expptr = locexp;
	for (i=0;i<Nterms; i++) {
		hx= hermitX[i]*srcstr;
		for (j=0;j<Nterms;j++) {
			hxy= hx * hermitY[j];
			for (k=0;k<Nterms;k++) {
				*expptr++ += hxy * hermitZ[k];
			}
		}
	}
}


/*
	BoxTaylorExpansion2D:
	convert all sources in box (srcbox) into taylor exp expansion about
	the center of target box (trgbox)
*/

void FastGauss::BoxTaylorExpansion2D(int srcbox,int trgbox)
{
	double * cent;
	double src2D[2];
	double * locexp;
	int src;
	unsigned int i;

	// cent[] : pointer to target box center
	cent=BoxCenter+trgbox*Dim;
	//locexp=LocalExp+trgbox*	Nterms2pDim;		
	//locexp=LocalExp;
	locexp=TrgBoxLocalExp[trgbox];


	// reset Local exp
	// reset all coeff to 0;
	//for (i=0;i<Nterms2pDim;i++)
	//	locexp[i]=0.0;


	src=SrcsInBox[srcbox];
	while (src >=0) {
	
		src2D[0]=SRC(src,0);
		src2D[1]=SRC(src,1);

		TaylorExpansion2D(src2D,str[src],cent,locexp,Nterms,fixedsegma);	

		src=SrcNext[src];
	}


}

/*
	BoxTaylorExpansion3D:
	convert all sources in box (srcbox) into taylor exp expansion about
	the center of target box (trgbox)
*/

void FastGauss::BoxTaylorExpansion3D(int srcbox,int trgbox)
{
	double * cent;
	double src3D[3];
	double * locexp;
	int src;
	unsigned int i;

	// cent[] : pointer to target box center
	cent=BoxCenter+trgbox*Dim;
	//locexp=LocalExp+trgbox*	Nterms2pDim;		
	
	locexp=TrgBoxLocalExp[trgbox];

	// reset Local exp
	// reset all coeff to 0;
	//for (i=0;i<Nterms2pDim;i++)
	//	locexp[i]=0.0;


	src=SrcsInBox[srcbox];
	while (src >=0) {
	
		src3D[0]=SRC(src,0);
		src3D[1]=SRC(src,1);
		src3D[2]=SRC(src,2);

		TaylorExpansion3D(src3D,str[src],cent,locexp,Nterms,fixedsegma);	

		src=SrcNext[src];
	}


}


/* 
	FarFieldExpansion2D
	Generate far feild expansion at the center of source box i 
	for all sources in that bos

*/
void FastGauss::FarFieldExpansion2D(int srcbox, double * ffexp)
{
	unsigned int i,j,k;
	double x,y,xp;
	double yp[MAXTERMS];
	double * FFexpptr;
	int src;

	// reset all coeff to 0;
	for (i=0;i<Nterms2pDim;i++)
		ffexp[i]=0.0;

	src=SrcsInBox[srcbox];
	while (src >= 0 ) {

		//for (i=0;i<Dim;i++) {
		//	Dtmp[i]=(SRC(src,i)-CENT(srcbox,i)) / SEGMA(src,i)
		//}
			
		x=(SRC(src,0)-CENT(srcbox,0)) / fixedsegma[0];
		y=(SRC(src,1)-CENT(srcbox,1)) / fixedsegma[1];

		yp[0]=1.0;
		for (k=1;k<Nterms;k++) {
			yp[k]=yp[k-1]*y / (double) k;
		}

		FFexpptr=ffexp;
		xp=str[src];
		for (i=0;i<Nterms;i++) {
			for (j=0;j<Nterms;j++) {
				// FFexp(i,j) +=xp(i) * yp(j)
				(*FFexpptr)+= xp * yp[j];
				FFexpptr++;
			}
			xp*= x / (double) (i+1);
		}
		src=SrcNext[src];
	}


}

/* 
	FarFieldExpansion3D
	Generate far feild expansion at the center of source box i 
	for all sources in that bos

*/
void FastGauss::FarFieldExpansion3D(int srcbox,double * ffexp)
{
	unsigned int i,j,k;
	double x,y,z,xp,xyp;
	double yp[MAXTERMS],zp[MAXTERMS];
	double * FFexpptr;
	int src;

	// reset all coeff to 0;
	for (i=0;i<Nterms2pDim;i++)
		ffexp[i]=0.0;

	src=SrcsInBox[srcbox];
	while (src >= 0 ) {

		//for (i=0;i<Dim;i++) {
		//	Dtmp[i]=(SRC(src,i)-CENT(srcbox,i)) / SEGMA(src,i)
		//}
			
		x=(SRC(src,0)-CENT(srcbox,0)) / fixedsegma[0];
		y=(SRC(src,1)-CENT(srcbox,1)) / fixedsegma[1];
		z=(SRC(src,2)-CENT(srcbox,2)) / fixedsegma[2];

		zp[0]=1.0;
		yp[0]=1.0;
		for (k=1;k<Nterms;k++) {
			zp[k]=zp[k-1]*z / (double) k;
			yp[k]=yp[k-1]*y / (double) k;
		}

		FFexpptr=ffexp;
		xp=str[src];
		//yp=1.0;
		for (i=0;i<Nterms;i++) {
			for (j=0;j<Nterms;j++) {
				xyp=xp*yp[j];
				// FFexp(i,j) +=xp(i) * yp(j)
				//xyp=xp * yp;
				for (k=0;k<Nterms;k++) {
					(*FFexpptr)+= xyp* zp[k];
					FFexpptr++;
				}
			}
			xp*= x / (double) (i+1);
		}
		src=SrcNext[src];
	}


}

/*
	EvalFarFieldExp_atTarget2D
	Evaluate far field exp about cent[] FFexp at target position trg[]
*/
double FastGauss::EvalFarFieldExp_atTarget2D(double * trg,double * cent,double * ffexp)
{

	unsigned int i,j;
	double x,y;
	double sum;
	double hermitX[MAXTERMS+1], hermitY[MAXTERMS+1];
	double * ffexpptr;

	x= (trg[0]-cent[0]) / fixedsegma[0];
	y= (trg[1]-cent[1]) / fixedsegma[1];

	TabulateHermitPoly2D(Nterms,x,y,hermitX,hermitY);

	sum=0;
	ffexpptr=ffexp;
	for (i=0;i<Nterms;i++) {
		for (j=0;j<Nterms;j++)
			// 			sum+=hermitX[i]*hermitY[j]* FFexp(i,j);
			sum+=hermitX[i]*hermitY[j]* (*ffexpptr++);
	}

	return sum;
}


/*
	EvalFarFieldExp_atTarget3D
	Evaluate far field exp about cent[] FFexp at target position trg[]
*/
double FastGauss::EvalFarFieldExp_atTarget3D(double * trg,double * cent,double * ffexp)
{

	unsigned int i,j,k;
	double x,y,z;
	double sum;
	double hermitX[MAXTERMS+1], 
				 hermitY[MAXTERMS+1],
				 hermitZ[MAXTERMS+1];
	double * ffexpptr;
	double hx,hxy;

	x= (trg[0]-cent[0]) / fixedsegma[0];
	y= (trg[1]-cent[1]) / fixedsegma[1];
	z= (trg[2]-cent[2]) / fixedsegma[2];

	TabulateHermitPoly3D(Nterms,x,y,z,hermitX,hermitY,hermitZ);

	sum=0;
	ffexpptr=ffexp;
	for (i=0;i<Nterms;i++) {
		hx=hermitX[i];
		for (j=0;j<Nterms;j++) {
			hxy=hx*hermitY[j];
			for (k=0;k<Nterms;k++) {
			// 			sum+=hermitX[i]*hermitY[j]* FFexp(i,j);
				sum+= hxy* hermitZ[k] * (*ffexpptr++);
			}
		}
	}

	return sum;
}


/*
	EvalFarFieldExp
	Evaluate far field exp FFexp at each target of box trgbox
*/
void FastGauss::EvalFarFieldExp(int trgbox,double * cent,double * ffexp)	
{
	double * trgptr;
	int trg;

	trg=TrgsInBox[trgbox];

	if (Dim==2) {
		while (trg >= 0) {
			trgptr=& (TRG(trg,0));
			TargetField[trg]+=EvalFarFieldExp_atTarget2D(trgptr,cent,ffexp);
			trg=TrgNext[trg];
		}
	}
	else if (Dim==3) {
		while (trg >= 0) {
			trgptr=& (TRG(trg,0));
			TargetField[trg]+=EvalFarFieldExp_atTarget3D(trgptr,cent,ffexp);			
			trg=TrgNext[trg];
		}
	}

}

void ShiftExpansion2D(double * srcboxcent, 
											double * trgboxcent,
											double * ffexp,
											double * localexp,
											unsigned int Nterms,
											double * segma)
{
	double x,y;
	double hermitX[2*MAXTERMS+1];
	double hermitY[2*MAXTERMS+1];
	double tempexp[MAXTERMS*MAXTERMS];
	double fac[MAXTERMS];


	unsigned int i,k,l,m;
	int j;

	double * ptr, *kptr, *jptr, *tempexpptr , * locexpptr;
	double sum;

	// reset all coeff to 0;
	//for (i=0;i<Nterms*Nterms;i++)
	//	localexp[i]=0.0;

	x=(trgboxcent[0]-srcboxcent[0]) / segma[0];
	y=(trgboxcent[1]-srcboxcent[1]) / segma[1];

	TabulateHermitPoly2D(2*Nterms-1,x,y,hermitX,hermitY);

	// generate factorial table
	fac[0]=1;0;
	for (i=1;i<Nterms;i++)
		fac[i]=-fac[i-1] / (double) i;


	tempexpptr=tempexp;
	for (j=0;j<Nterms;j++) {
		ptr=ffexp+j*Nterms;
		for (m=0; m <Nterms; m++) {
			sum=0.0;
			kptr=ptr;
			for (k=0;k<Nterms;k++) {
				// sum+=FFexp(j,k)*hermit(m+k)
				sum+=(*kptr++) * hermitY[m+k];
			}
			*tempexpptr++ =sum;		// tempexp(j,m)=sum
		}
	}



	locexpptr=localexp;
	for (l=0;l<Nterms;l++) {
		//ptr=tempexp+l*Nterms;
		for (m=0;m<Nterms;m++) {
			sum=0.0;
			
			jptr=tempexp+m;
			for (j=0;j<Nterms;j++) {
				sum+=(*jptr) * hermitX[l+j];
				jptr+=Nterms;
			}
			
			/*

			jptr=tempexp+Nterms*(Nterms-1)+m;
			for (j=Nterms-1;j>=0;j--) {
				sum+=(*jptr) * hermitX[l+j];
				jptr-=Nterms;
			}
			
			*/

			*locexpptr++ += sum * fac[l] * fac[m];
		}
	}
	

}

void ShiftExpansion3D(double * srcboxcent, 
											double * trgboxcent,
											double * ffexp,
											double * localexp,
											unsigned int Nterms,
											double * segma)
{
	double x,y,z;
	double hermitX[2*MAXTERMS+1];
	double hermitY[2*MAXTERMS+1];
	double hermitZ[2*MAXTERMS+1];

	double tempexp1[MAXTERMS*MAXTERMS*MAXTERMS];
	double tempexp2[MAXTERMS*MAXTERMS*MAXTERMS];

	double fac[MAXTERMS];


	int i,j,k,l,m,n;
	int nterms = (int) Nterms;
	int ntermsp2= nterms*nterms;

	double * ptr, *kptr, *jptr, *iptr, *tempexpptr , * locexpptr;
	double sum;

	// reset all coeff to 0;
	//for (i=0;i<Nterms*Nterms*Nterms;i++)
	//	localexp[i]=0.0;

	
	x=(trgboxcent[0]-srcboxcent[0]) / segma[0];
	y=(trgboxcent[1]-srcboxcent[1]) / segma[1];
	z=(trgboxcent[2]-srcboxcent[2]) / segma[2];

	TabulateHermitPoly3D(2*Nterms-1,x,y,z,hermitX,hermitY,hermitZ);

	// generate factorial table
	fac[0]=1;0;
	for (i=1;i<Nterms;i++)
		fac[i]=-fac[i-1] / (double) i;


	// tempexp1(i,j,n)= SUM_k FFexp(i,j,k) * hz(n+k)
	
	tempexpptr=tempexp1;
	for (i=0;i<nterms;i++) {
		for (j=0; j<nterms;j++) {
			// ptr = ffexp(i,j,0)
			ptr=ffexp+(i*nterms+j)*nterms;
			for( n=0; n<nterms; n++) {
				sum=0.0;
				kptr=ptr;
				for ( k=0; k<nterms;k++) {
					sum+=(*kptr++) * hermitZ[n+k];
				}

				*tempexpptr++=sum; // tempexp1(i,j,n)=sum
			}
		}
	}

	// tempexp2(i,m,n) = SUM_j tempexp1(i,j,n) * hy(m+j)
	tempexpptr=tempexp2;
	for (i=0;i<nterms;i++) {
		for (m=0; m<nterms;m++) {
			for( n=0; n<nterms; n++) {
				// ptr = tempexp1(i,0,n)
				jptr=tempexp1+i*ntermsp2+n;
				sum=0.0;
				for (j=0;j<nterms;j++) {
					sum+= (*jptr) * hermitY[m+j];
					jptr+=nterms;
				}
		
				*tempexpptr++ =sum;
			}
		}
	}

	// locexp(l,m,n) = SUM_i tempexp2(i,m,n) * hx(l+i)
	locexpptr=localexp;
	for (l=0;l<nterms;l++) {
		for (m=0; m<nterms;m++) {
			for( n=0; n<nterms; n++) {
				// ptr = tempexp2(0,m,n)
				iptr=tempexp2+m*nterms+n;
				sum=0.0;

				for (i=0;i<nterms;i++){
					sum+= (* iptr) * hermitX[l+i];
					iptr+=ntermsp2;
				}

				*locexpptr++ += sum * fac[l] * fac [m] * fac[n];
			}
		}
	}


}

/*
	EvalLocalExp_attarget2D :

	Evaluate local exp 'locexp' about 'cent[]' at target 'trg[]'

*/
double EvalLocalExp_attarget2D(double * locexp,
															 double * cent,
															 double * trg,
															 unsigned int Nterms,
															 double * segma)
{
	double x,y;
	double sum,parsum;
	int j,k;
	double * locexpptr;

	x=(trg[0]-cent[0]) / segma[0];
	y=(trg[1]-cent[1]) / segma[1];
	
	sum=0.0;

	locexpptr=locexp+Nterms*Nterms-1;  // locexpptr = locexp(Nterms,Nterms)

	for (j= Nterms-1;j>0;j--) {
		parsum=0.0;
		for (k= Nterms-1;k>0;k--) {
			parsum= ( parsum+ (*locexpptr--)) * y;	// localexp(j,k)
		}
		parsum += *locexpptr-- ; // localexp(j,0)
		sum = (sum+parsum)*x;
	}

	parsum=0.0;
	
	for (k=Nterms-1;k > 0 ; k--) {
		parsum = (parsum + locexp[k]) * y;		// locexp(0,k)
	}

	parsum += *locexp;	//locexp(0,0)

	sum+=parsum;

	return sum;



}

/*
	EvalLocalExp_attarget3D :

	Evaluate local exp 'locexp' about 'cent[]' at target 'trg[]'

*/
double EvalLocalExp_attarget3D(double * locexp,
															 double * cent,
															 double * trg,
															 unsigned int Nterms,
															 double * segma)
{
	double x,y,z;
	double sum,parsumy,parsumz;
	int i,j,k;
	double * locexpptr;

	x=(trg[0]-cent[0]) / segma[0];
	y=(trg[1]-cent[1]) / segma[1];
	z=(trg[2]-cent[2]) / segma[2];
	
	sum=0.0;

	locexpptr=locexp+Nterms*Nterms*Nterms-1;  // locexpptr = locexp(Nterms,Nterms,Nterms)

	for( i=Nterms-1;i>=0;i--) {
		parsumy=0;
		for (j= Nterms-1;j>=0;j--) {
			parsumz=0.0;
			for (k= Nterms-1;k>=0;k--) {
				parsumz=  (*locexpptr--) + z * parsumz;	// localexp(i,j,k)
			}

			parsumy=parsumz + y * parsumy;
		}


		sum = parsumy+ x * sum;
	}


	return sum;



}

/*
	Evaluate local expansions for targets in one box (trgbox).
*/
void FastGauss::EvalBoxLocalExpansion(int trgbox,double * locexp)
{
	double trg2D[2],trg3D[3];
	double * cent;
	int trg;
	
	if (Dim==2) {
		cent=BoxCenter+trgbox*Dim;						// pointer to box center

		trg=TrgsInBox[trgbox];
		while (trg >= 0) {

			// 2D code
			trg2D[0] = TRG(trg,0);
			trg2D[1] = TRG(trg,1);
			TargetField[trg] += EvalLocalExp_attarget2D(locexp,cent,trg2D,Nterms,fixedsegma);

			trg=TrgNext[trg];
		}		
	
	}
	else if (Dim==3) {
		cent=BoxCenter+trgbox*Dim;						// pointer to box center

		trg=TrgsInBox[trgbox];
		while (trg >= 0) {

			// 2D code
			trg3D[0] = TRG(trg,0);
			trg3D[1] = TRG(trg,1);
			trg3D[2] = TRG(trg,2);

			TargetField[trg] += EvalLocalExp_attarget3D(locexp,cent,trg3D,Nterms,fixedsegma);

			trg=TrgNext[trg];
		}		

	}

}

/*
	Evaluate local expansions for target boxes.
*/
void FastGauss::EvalLocalExpansions() 
{
	int trg;
	double * locexp;
	double trg2D[2],trg3D[3];
	double * cent;
	unsigned int i;	
	int M;

	if (Dim==2) {
		for (i=0;i<NBoxes;i++) {
			M=TrgsInBoxCount[i];
			if (M > NtMAX) {
				locexp=LocalExp+i*Nterms2pDim;	// pointer to local exp of box i
				cent=BoxCenter+i*Dim;						// pointer to box center

				trg=TrgsInBox[i];
				while (trg >= 0) {

					// 2D code
					trg2D[0] = TRG(trg,0);
					trg2D[1] = TRG(trg,1);
					TargetField[trg] += EvalLocalExp_attarget2D(locexp,cent,trg2D,Nterms,fixedsegma);

					trg=TrgNext[trg];
				}		
			}
		}
	} 
	else if (Dim==3) {
		for (i=0;i<NBoxes;i++) {
			M=TrgsInBoxCount[i];
			if (M > NtMAX) {
				locexp=LocalExp+i*Nterms2pDim;	// pointer to local exp of box i
				cent=BoxCenter+i*Dim;						// pointer to box center

				trg=TrgsInBox[i];
				while (trg >= 0) {

					// 2D code
					trg3D[0] = TRG(trg,0);
					trg3D[1] = TRG(trg,1);
					trg3D[2] = TRG(trg,2);
					TargetField[trg] += EvalLocalExp_attarget3D(locexp,cent,trg3D,Nterms,fixedsegma);

					trg=TrgNext[trg];
				}		
			}
		}

		
	}

}


int FastGauss::DetermineNterms(double tol,double r)
{

	double r2,z;
	int p;
	

	r2=r*r;

	z=(1.09)*(1.09) * r2 * r2  / ((1-r)*(1-r));
	p=1;
	do {
		p++;
		z*=r2 / (double) p;

	} while ( z >= tol);

	return (p);

}

int DetermineNtermsII(double tol,double r)
{

	double rp,z,e;
	int p;
	
	e=exp(1);

	p=(int) (r*r*e);
	do {
		p++;
		rp=sqrt(e/(p+1.0))*r;
		z=pow(rp,p)/(1-rp);

	} while ( z >= tol);

	return (p);

}


/*
int FastGauss::DetermineNterms(double tol)
{

	double r,r2,z;
	int p;
	unsigned int i;
	double maxr;

	maxr=0;

	for (i=0;i<Dim;i++) {
		r=1.0 / (Nboxes_eachside[i] * sqrt(2) * fixedsegma[i]);
		if (r > maxr)
			maxr=r;
	}
	r=maxr;

	r2=r*r;

	z=(1.09)*(1.09) * r2 * r2  / ((1-r)*(1-r));
	p=1;
	do {
		p++;
		z*=r2 / (double) p;

	} while ( z >= tol);

	return (p);

}
*/

int FastGauss::Evaluateboxes2D()
{

	unsigned int i,j;
	int N,M;
	int c1,c2,c3,c4;
	int trgboxidx;

	unsigned int NbrCount;
	//int * NbrList = new int [NBoxes];


	c1=0;
	c2=0;
	c3=0;
	c4=0;


	// process all boxes
	for (i=0;i<NBoxes;i++) {
		N=SrcsInBoxCount[i];
		
		if (N==0)
			continue;

		if (N <= NsMAX) {
			GenerateNbrs(i,NbrList,NbrCount);
			for (j=0;j<NbrCount;j++) {
				trgboxidx=NbrList[j];
				M=TrgsInBoxCount[trgboxidx];
				if (M > 0) {
					if (M <= NtMAX) {
						// direct evaluation
						DirectEvaluation(i,trgboxidx);
						c1++;
		
					}
					else // M > NLMAX
					{	
						// covert each source into taylor series
						BoxTaylorExpansion2D(i,trgboxidx);
						c2++;

					} 
				}
			}

		}
		else // N > NFMAX
		{
			double * srcboxcent;
			double * trgboxcent;

			srcboxcent=BoxCenter+i*Dim;  // pointer to center of box i (NDim vector)
			FarFieldExpansion2D(i,FFexp);
			GenerateNbrs(i,NbrList,NbrCount);
			for (j=0;j<NbrCount;j++) {
				trgboxidx=NbrList[j];
				M=TrgsInBoxCount[trgboxidx];
				if (M > 0) {
					if (M <= NtMAX) {
						// evaluate Far Field exp at all targets at box j
						EvalFarFieldExp(trgboxidx,srcboxcent,FFexp);					
						c3++;
					}
					else // M > NLMAX
					{
						// Convert FFexp to Localexp
						trgboxcent=BoxCenter+trgboxidx*Dim;

						ShiftExpansion2D(srcboxcent,trgboxcent,
							FFexp,LocalExp+trgboxidx*Nterms2pDim,Nterms,fixedsegma);
						
						c4++;
					} 
				}
			}

		}

	}	// for all boxes

	// evaluate local expansions
	//EvalLocalExpansions();

	
	double ct = (c1+c2+c3+c4) / 100.0;

	cout << endl;
	cout << " Statistics 2D:\n";
	cout << " Direct evaluation : " << c1 / ct << endl;
	cout << " taylor expansion  : " << c2 / ct << endl;
	cout << " direct Hermit expansion  : " << c3 / ct<< endl;
	cout << " Hermit expansion & taylor expansion  : " << c4  / ct<< endl;
	cout << endl;

	return 1;
}

int FastGauss::Evaluateboxes3D()
{
	unsigned int i,j;
	int N,M;
	int c1,c2,c3,c4;
	int trgboxidx;

	unsigned int NbrCount;
	int * NbrList = new int [NBoxes];


	c1=0;
	c2=0;
	c3=0;
	c4=0;


	// process all boxes
	for (i=0;i<NBoxes;i++) {
		N=SrcsInBoxCount[i];
		
		if (N==0)
			continue;

		if (N <= NsMAX) {
			GenerateNbrs(i,NbrList,NbrCount);
			for (j=0;j<NbrCount;j++) {
				trgboxidx=NbrList[j];
				M=TrgsInBoxCount[trgboxidx];
				if (M > 0) {
					if (M <= NtMAX) {
						// direct evaluation
						DirectEvaluation(i,trgboxidx);
						c1++;
		
					}
					else // M > NLMAX
					{	
						// covert each source into taylor series
						BoxTaylorExpansion3D(i,trgboxidx);
						c2++;

					} 
				}
			}

		}
		else // N > NFMAX
		{
			double * srcboxcent;
			double * trgboxcent;

			srcboxcent=BoxCenter+i*Dim;  // pointer to center of box i (NDim vector)
			FarFieldExpansion3D(i,FFexp);
			GenerateNbrs(i,NbrList,NbrCount);
			for (j=0;j<NbrCount;j++) {
				trgboxidx=NbrList[j];
				M=TrgsInBoxCount[trgboxidx];
				if (M > 0) {
					if (M <= NtMAX) {
						// evaluate Far Field exp at all targets at box j
						EvalFarFieldExp(trgboxidx,srcboxcent,FFexp);					
						c3++;
					}
					else // M > NLMAX
					{
						// Convert FFexp to Localexp
						trgboxcent=BoxCenter+trgboxidx*Dim;

						ShiftExpansion3D(srcboxcent,trgboxcent,
							FFexp,LocalExp+trgboxidx*Nterms2pDim,Nterms,fixedsegma);
						
						c4++;
					} 
				}
			}

		}

	}	// for all boxes

	// evaluate local expansions
	//EvalLocalExpansions();


	double ct = (c1+c2+c3+c4) / 100.0;

	cout << endl;
	cout << " Statistics 3D:\n";
	cout << " Direct evaluation : " << c1 / ct << endl;
	cout << " taylor expansion  : " << c2 / ct << endl;
	cout << " direct Hermit expansion  : " << c3 / ct<< endl;
	cout << " Hermit expansion & taylor expansion  : " << c4  / ct<< endl;
	cout << endl;


	return 1;
}


int FastGauss::FixedScaleEvaluation(double * source_locs,
																		double * source_strength,
																		double * source_segma,
																		double * target_locs,
																		double * & FFresults,
																		unsigned int NSources,
																		unsigned int NTargets,
																		double tolerance)
{

	// copy input to internal data structure
	Nsources=NSources;
	Ntargets=NTargets;

	
	sources	= source_locs;
	str			= source_strength;
	segma		= source_segma;
	targets	= target_locs;

	fixedsegma=segma;		// assume fixed segma for all sources.

	unsigned int i,j;
	double maxr,r,ro=0.5;



	double tol=tolerance;

	
	// determine number of boxes each side
	// initial r=0.5
	//Nboxes_eachside=(int) (1.0 / (sqrt(0.5) * fixedsegma) + 1 );

	Nboxes_eachside= new int[Dim];
	Nboxes_powers = new int [Dim*Dim];


	NBoxes=1;
	maxr=0;

	for (i=0;i<Dim;i++) {
		Nboxes_eachside[i]=(int) ceil(1.0 / (ro * sqrt(2) * fixedsegma[i])) ;
		r=1.0 / (Nboxes_eachside[i] * sqrt(2) * fixedsegma[i]);
		if (r > maxr)
			maxr=r;

		Nboxes_powers[i]=1;
		NBoxes *= Nboxes_eachside[i];
	}

	r=maxr;




	//Nboxes_powers[0]=1;
	for (i=0;i<Dim;i++) 
		for (j=1; j< Dim;j++){
		Nboxes_powers[j*Dim+i]=Nboxes_powers[(j-1)*Dim+i]*Nboxes_eachside[i];
	}

	//NBoxes=Nboxes_powers[Dim-1]*Nboxes_eachside;		


	// determine box Nbr size
	// assume r=0.5
	NbrSize=(int) floor(sqrt(-2*log(tol))+1);
	
	//r=1.0 / (Nboxes_eachside * sqrt(2) * fixedsegma);
	//NbrSize=(int) floor(sqrt(-log(tol)/2.0)/r +1);


	// determine number of hermit expansion terms
	//Nterms=DetermineNterms(tol,r)+1;
	Nterms=DetermineNtermsII(tol,r)+1;

	//NsMAX=80;
	//NsMAX=5;
	//NtMAX=Nterms+2;
	//NtMAX=1000;

	NsMAX=Nterms;
	NtMAX=Nterms;


	//MAX=0;
	//MAX=0;


	Ndisp=1;
	Nterms2pDim=1;
	for(i=0;i<Dim;i++) {
		Ndisp*=(2*NbrSize+1);
		Nterms2pDim	*= Nterms;
	}

	


	cout <<  " Fast Gauss setting : \n"; 
	cout << " tol = " << tol << endl;
	for (i=0;i<Dim;i++) 
		cout << " Number of boxes per dim " <<i  <<" = "<< Nboxes_eachside[i] << endl;

	cout << " Box Neighborhood per dim = " << (2*NbrSize+1) << endl;
	cout << " Number of terms = " << Nterms << endl;
	cout << " r= " << r << endl;

	// allocate buffers;
	if (! Allocate())
		return 0;

	TargetField=FFresults;

	if (! AllocateTargets())
			return 0;

	// reset result buffer
	for (i=0;i<Ntargets;i++)
		TargetField[i]=0.0;

	/*
	// reset local expansions
	for (i=0;i<NBoxes*Nterms2pDim;i++)
		LocalExp[i]=0;
	*/

	// generate Nbrs relative addresses
	GenerateDisplacements();

	// assign sources and targets into boxes
	GenerateBoxes();
	AssignSources();
	AssignTargets();


	if (Dim==2)
		Evaluateboxes2D();	
	else if (Dim==3)
		Evaluateboxes3D();	
	


	FFresults=TargetField;


	/*
	double ct = (c1+c2+c3+c4) / 100.0;

	cout << endl;
	cout << " Statistics :\n";
	cout << " Direct evaluation : " << c1 / ct << endl;
	cout << " taylor expansion  : " << c2 / ct << endl;
	cout << " direct Hermit expansion  : " << c3 / ct<< endl;
	cout << " Hermit expansion & taylor expansion  : " << c4  / ct<< endl;
	cout << endl;
	*/

	return 1;
}

int FastGauss::SourceExpansion2D()
{
	unsigned int i;
	int N;
	double * ffexpptr;

	ffexpptr=FFexp;
	Nffexp=0;

	// process all boxes
	for (i=0;i<NBoxes;i++) {
		N=SrcsInBoxCount[i];
		
		if ( N > NsMAX) {
			
			// allocate ffexp space
			SrcBoxFFexp[i]=ffexpptr;
			Nffexp++;
			assert(Nffexp < MAXFFEXP );

			FarFieldExpansion2D(i,ffexpptr);

			// increment ffexpptr
			ffexpptr+= Nterms2pDim;
		}
		else
			SrcBoxFFexp[i]=(double * ) NULL;

	}

	return 1;
}

int FastGauss::SourceExpansion3D()
{
	unsigned int i;
	int N;
	double * ffexpptr;

	ffexpptr=FFexp;
	Nffexp=0;

	// process all boxes
	for (i=0;i<NBoxes;i++) {
		N=SrcsInBoxCount[i];
		
		if ( N > NsMAX) {
			
			// allocate ffexp space
			SrcBoxFFexp[i]=ffexpptr;
			Nffexp++;

			FarFieldExpansion3D(i,ffexpptr);

			// increment ffexpptr
			ffexpptr+= Nterms2pDim;
		}
		else
			SrcBoxFFexp[i]=(double * ) NULL;

	}

	return 1;


}

int FastGauss::AllocateLocalExpansion()
{
		
	int utb;
	unsigned int i,j;
	int M;
	double * lptr;

	for (utb=0; utb < UsedTrgsBoxCount ; utb ++) {
		i=UsedTrgsBox[utb];
		M=TrgsInBoxCount[i];

		if (M > NtMAX && (TrgBoxLocalExp[i]==NULL)) {
			// allocate local exp space
			lptr= new double [Nterms2pDim];

			TrgBoxLocalExp[i]=lptr;

			// reset expansion
			for (j=0;j<Nterms2pDim;j++)
				lptr[j]=0.0;


			NLocalExp ++;		
		}
	}

	return 1;
}


int FastGauss::EvaluateTargetboxes2D()
{
	
	unsigned int i,j;
	int N,M;
	int c1,c2,c3,c4;
	int srcboxidx;
	double * srcboxcent;
	double * trgboxcent;


	unsigned int NbrCount;
	int * NbrList = new int [NBoxes];

	c1=0;
	c2=0;
	c3=0;
	c4=0;

	// process all boxes
	for (i=0;i<NBoxes;i++) {
		M=TrgsInBoxCount[i];
		if (M == 0) 
			continue;

		if (M <= NtMAX) {
			GenerateNbrs(i,NbrList,NbrCount);
			for (j=0;j<NbrCount;j++) {
				srcboxidx=NbrList[j];
				N=SrcsInBoxCount[srcboxidx];
				if (N > 0) {
					if ( N <=NsMAX) {
						// direct evaluation
						DirectEvaluation(srcboxidx,i);
						c1++;
					}
					else // N > NsMAX
					{
						srcboxcent=BoxCenter+srcboxidx*Dim;  // pointer to center of box i (NDim vector)
						// evaluate Far Field exp at all targets at box j
						assert(SrcBoxFFexp[srcboxidx]);
						EvalFarFieldExp(i,srcboxcent,SrcBoxFFexp[srcboxidx]);					
						c3++;
					}

				}
			}
		}
		else  // M > NtMAX
		{
			if ( ! TrgBoxLocalExpFlag[i]) {

				GenerateNbrs(i,NbrList,NbrCount);
				trgboxcent=BoxCenter+i*Dim;
				for (j=0;j<NbrCount;j++) {
					srcboxidx=NbrList[j];
					N=SrcsInBoxCount[srcboxidx];
					if (N>0) {
						if (N <=NsMAX) {
							// covert each source into taylor series
							BoxTaylorExpansion2D(srcboxidx,i);
							c2++;

						}
						else // N > NsMax
						{
							// Convert FFexp to Localexp
							srcboxcent=BoxCenter+srcboxidx*Dim;  // pointer to center of box srcboxidx (NDim vector)

							assert(SrcBoxFFexp[srcboxidx]);
							assert(TrgBoxLocalExp[i]);
							ShiftExpansion2D(srcboxcent,trgboxcent,
								SrcBoxFFexp[srcboxidx],TrgBoxLocalExp[i],Nterms,fixedsegma);
							
							c4++;

						}

					}
				}
				TrgBoxLocalExpFlag[i]=true;

			}
			EvalBoxLocalExpansion(i,TrgBoxLocalExp[i]);


		}
	}


	double ct = (c1+c2+c3+c4) / 100.0;

	cout << endl;
	cout << " Statistics 2D:\n";
	cout << " Direct evaluation : " << c1 / ct << endl;
	cout << " taylor expansion  : " << c2 / ct << endl;
	cout << " direct Hermit expansion  : " << c3 / ct<< endl;
	cout << " Hermit expansion & taylor expansion  : " << c4  / ct<< endl;
	cout << endl;


	return 1;


}

int FastGauss::EvaluateTargetboxes3D()
{

	unsigned int i,j;
	int N,M;
	int c1,c2,c3,c4;
	float t1,t2,t3,t4,t5;
	int srcboxidx;
	double * srcboxcent;
	double * trgboxcent;
	int utb;

	unsigned int NbrCount;
	// int * NbrList = new int [NBoxes];

	c1=0;
	c2=0;
	c3=0;
	c4=0;

	t1=0;
	t2=0;
	t3=0;
	t4=0;
	t5=0;

	// process all boxes
	//for (i=0;i<NBoxes;i++) {
	for (utb=0; utb < UsedTrgsBoxCount ; utb ++) {
		i=UsedTrgsBox[utb];
		M=TrgsInBoxCount[i];
		

		//if (M == 0) 
		//	continue;

		//assert (M > 0);

		if (M <= NtMAX) {
			GenerateNbrs(i,NbrList,NbrCount);
			for (j=0;j<NbrCount;j++) {
				srcboxidx=NbrList[j];
				N=SrcsInBoxCount[srcboxidx];
				if (N > 0) {
					if ( N <=NsMAX) {
						// direct evaluation
						//FG_TIMED_FUNCTION(DirectEvaluation(srcboxidx,i),t2);
						DirectEvaluation(srcboxidx,i);
						c1++;
					}
					else // N > NsMAX
					{
						srcboxcent=BoxCenter+srcboxidx*Dim;  // pointer to center of box i (NDim vector)
						// evaluate Far Field exp at all targets at box j
						assert(SrcBoxFFexp[srcboxidx]);
						//TIMED_FUNCTION(EvalFarFieldExp(i,srcboxcent,SrcBoxFFexp[srcboxidx]),t3);					
						EvalFarFieldExp(i,srcboxcent,SrcBoxFFexp[srcboxidx]);					
						c3++;
					}

				}
			}
		}
		else  // M > NtMAX
		{
			

			if ( ! TrgBoxLocalExpFlag[i]) {
				GenerateNbrs(i,NbrList,NbrCount);
				trgboxcent=BoxCenter+i*Dim;
				for (j=0;j<NbrCount;j++) {
					srcboxidx=NbrList[j];
					N=SrcsInBoxCount[srcboxidx];
					if (N>0) {
						
						if (N <=NsMAX) {
							// covert each source into taylor series
							//TIMED_FUNCTION(BoxTaylorExpansion3D(srcboxidx,i),t2);
							BoxTaylorExpansion3D(srcboxidx,i);
							c2++;

						}
						
						
						else // N > NsMax
						{
							// Convert FFexp to Localexp
							srcboxcent=BoxCenter+srcboxidx*Dim;  // pointer to center of box srcboxidx (NDim vector)

							assert(SrcBoxFFexp[srcboxidx]);
							assert(TrgBoxLocalExp[i]);

							//TIMED_FUNCTION(ShiftExpansion3D(srcboxcent,trgboxcent,SrcBoxFFexp[srcboxidx],LocalExp,Nterms,fixedsegma),t4);
							ShiftExpansion3D(srcboxcent,trgboxcent,SrcBoxFFexp[srcboxidx],TrgBoxLocalExp[i],Nterms,fixedsegma);
							
							c4++;

						}					
					}
				}
				TrgBoxLocalExpFlag[i]=true;
			}
			//TIMED_FUNCTION(EvalBoxLocalExpansion(i,LocalExp),t5);
			EvalBoxLocalExpansion(i,TrgBoxLocalExp[i]);


		}
	}

	
	double ct = (c1+c2+c3+c4) / 100.0;
	/*
	stat_1=c1 / ct;
	stat_2=c2 / ct;
	stat_3=c3 / ct;
	stat_4=c4 / ct;
  */
	
	timer_1=t1;
	timer_2=t2;

	

	
	cout << endl;
	cout << " Statistics 3D:\n";
	cout << " Direct evaluation : " << c1 / ct  <<endl;
	cout << " taylor expansion  : " << c2 / ct << endl;
	cout << " direct Hermit expansion  : " << c3 / ct<< endl;
	cout << " Hermit expansion & taylor expansion  : " << c4  / ct<<endl;
	cout << endl;
	



	return 1;
}

int FastGauss::BuildSourceStructure(unsigned int NSources,
																		double * source_locs,
																		double * source_strength,
																		double * source_segma,
																		double tolerance)
{

	// copy input to internal data structure
	Nsources=NSources;

	sources	= source_locs;
	str			= source_strength;
	segma		= source_segma;

	fixedsegma=segma;		// assume fixed segma for all sources.
	unsigned int i,j;
	double maxr,r,ro=2.0;

	double tol=tolerance;

	// determine number of boxes each side
	// initial r=0.5
	//Nboxes_eachside=(int) (1.0 / (sqrt(0.5) * fixedsegma) + 1 );

	Nboxes_eachside= new int[Dim];
	Nboxes_powers = new int [Dim*Dim];


	NBoxes=1;
	maxr=0;

	for (i=0;i<Dim;i++) {
		Nboxes_eachside[i]=(int) ceil(1.0 / (ro * sqrt(2) * fixedsegma[i])) ;
		r=1.0 / (Nboxes_eachside[i] * sqrt(2) * fixedsegma[i]);
		if (r > maxr)
			maxr=r;

		Nboxes_powers[i]=1;
		NBoxes *= Nboxes_eachside[i];
	}

	r=maxr;




	//Nboxes_powers[0]=1;
	for (i=0;i<Dim;i++) 
		for (j=1; j< Dim;j++){
		Nboxes_powers[j*Dim+i]=Nboxes_powers[(j-1)*Dim+i]*Nboxes_eachside[i];
	}

	//NBoxes=Nboxes_powers[Dim-1]*Nboxes_eachside;		


	// determine box Nbr size
	// assume r=0.5
	//NbrSize=(int) floor(sqrt(-2*log(tol))+1);
	
	//r=1.0 / (Nboxes_eachside * sqrt(2) * fixedsegma);
	NbrSize=(int) floor(sqrt(-log(tol)/2.0)/r +1);


	// determine number of hermit expansion terms
	//Nterms=DetermineNterms(tol,r)+1;
	Nterms=DetermineNtermsII(tol,r)+1;

	Ndisp=1;
	Nterms2pDim=1;
	for(i=0;i<Dim;i++) {
		Ndisp*=(2*NbrSize+1);
		Nterms2pDim	*= Nterms;
	}

	NsMAX=Nterms2pDim/Nterms;
	NtMAX=Nterms2pDim/Nterms;

	//NsMAX=Nterms-2;
	//NtMAX=Nterms-2;


	//NsMAX=80;
	//NtMAX=1000;

	//NsMAX=2000;
	//NtMAX=20;

	//MAX=0;
	//MAX=0;



	


	cout <<  " Fast Gauss setting : \n"; 
	cout << " tol = " << tol << endl;
	for (i=0;i<Dim;i++) 
		cout << " Number of boxes per dim " <<i  <<" = "<< Nboxes_eachside[i] << endl;

	cout << " Box Neighborhood per dim = " << (2*NbrSize+1) << endl;
	cout << " Number of terms = " << Nterms << endl;
	cout << " r= " << r << endl;

	// allocate buffers;
	if (! Allocate())
		return 0;


	/*
	// reset local expansions
	for (i=0;i<NBoxes*Nterms2pDim;i++)
		LocalExp[i]=0;
	*/

	// generate Nbrs relative addresses
	GenerateDisplacements();

	// assign sources and targets into boxes
	GenerateBoxes();
	AssignSources();

	if (Dim==2)
		SourceExpansion2D();	
	else if (Dim==3)
		SourceExpansion3D();	


	// reset local exp pointers and flags
	for(i=0; i<NBoxes; i++) {
		TrgBoxLocalExp[i]=NULL;
		TrgBoxLocalExpFlag[i]=false;		
	}

	NLocalExp=0;

	cout << " Number of FF exp =" << Nffexp << endl;

	return 1;
}



int FastGauss::EvaluateTargets(unsigned int NTargets,
															 double * target_locs,
															 double * & FFresults)
{
	unsigned int i;

	Ntargets=NTargets;
	targets	= target_locs;

	TargetField=FFresults;

	if ( ! PreAllocateTargetFlag || ! TargetField) {
		if (! AllocateTargets())
			return 0;
	}

	//timer_3=0.0;
	AssignTargets();

	// Allocate necessary Local Expansions space
	AllocateLocalExpansion();

	cout << " Number of Local exp =" << NLocalExp << endl;


	// reset result buffer
	for (i=0;i<Ntargets;i++)
		TargetField[i]=0.0;

	//timer_4=0.0;

	// evaluation
	if (Dim==2)
		EvaluateTargetboxes2D();	
	else if (Dim==3)
		EvaluateTargetboxes3D();	



	if (FFresults== NULL)
		FFresults=TargetField;

	return 1;
}




