#include "ObjectTracker.h"
//#include "VCAevent.h"
#include "ImageLib.h"
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


void FilterII(int * FGimage,
						unsigned char * object_label_image,
						int ObjectLabel,
						int P1label,
						int P2label,
						unsigned char * outimage,
						unsigned int rows,
						unsigned int cols)
{
	int Nbr[9];
	unsigned int r,c,i;
	int N1,N2;
	int idx,idx2;
	
	Nbr[0]=-cols-1;
	Nbr[1]=-cols;
	Nbr[2]=-cols+1;
	Nbr[3]=-1;
	Nbr[4]=0;
	Nbr[5]=1;
	Nbr[6]=cols-1;
	Nbr[7]=cols;
	Nbr[8]=cols+1;

	for (r=1;r<rows-1;r++)
		for (c=1;c<cols-1;c++) {
			idx=r*cols+c;
			if (FGimage[idx]==ObjectLabel)
				//if (object_label_image[idx] == 0) {
				if (object_label_image[idx] == P1label ||
					object_label_image[idx] == P2label ) {
					N1=0;
					N2=0;
					for (i=0;i<9;i++) {
						idx2=Nbr[i]+idx;
						if (object_label_image[idx2]==P1label)
							N1++;
						else if (object_label_image[idx2]==P2label)
							N2++;
					}
					if (N1-N2>1)
							outimage[idx]=P1label;
					else if (N2-N1>1)
							outimage[idx]=P2label;
				}
				else
					outimage[idx]=object_label_image[idx];
					
		}
}

void FilterI(int * FGimage,
						unsigned char * object_label_image,
						double * P1,
						double * P2,
						int ObjectLabel,
						int P1label,
						int P2label,
						unsigned char * outimage,
						unsigned int rows,
						unsigned int cols)
{
	int Nbr[9];
	unsigned int r,c,i;
	int N1,N2;
	int idx,idx2;
	
	Nbr[0]=-cols-1;
	Nbr[1]=-cols;
	Nbr[2]=-cols+1;
	Nbr[3]=-1;
	Nbr[4]=0;
	Nbr[5]=1;
	Nbr[6]=cols-1;
	Nbr[7]=cols;
	Nbr[8]=cols+1;

	for (r=1;r<rows-1;r++)
		for (c=1;c<cols-1;c++) {
			idx=r*cols+c;
			if (FGimage[idx]==ObjectLabel)
				//if (object_label_image[idx] == 0) {
				if (object_label_image[idx] != P1label &&
					object_label_image[idx] != P2label ) {
					N1=0;
					N2=0;
					for (i=0;i<9;i++) {
						idx2=Nbr[i]+idx;
						if (object_label_image[idx2]==P1label)
							N1++;
						else if (object_label_image[idx2]==P2label)
							N2++;
					}
					if (N1>N2)
							outimage[idx]=P1label;
					else if (N2>N1)
							outimage[idx]=P2label;
				}
				else
					outimage[idx]=object_label_image[idx];
					
		}
}

void GenerateOccModelOverlap(MovingObject ** ObjectsInMergeList,
											unsigned int ObjectsInMerge_cnt,
											unsigned int maxdepthindex,
											unsigned char * OccModelImage,
											int * OccList,
											int & OccListCount,
											unsigned int rows,
											unsigned int cols)
{

	Box tBox;
	unsigned int d,k;
	int l;
	int x,y;
	double a,b,_a2,_b2;
	int xo,yo;

	memset(OccModelImage,0,rows*cols);
	OccListCount=0;

	for (d=1;d<=maxdepthindex;d++)  {
		for (k=0;k<ObjectsInMerge_cnt;k++){
			if (ObjectsInMergeList[k]->OccDepth==d) {
				tBox=ObjectsInMergeList[k]->BBox;
				l=ObjectsInMergeList[k]->OccLabel;

				xo=(tBox.maxx+tBox.minx) /2;
				yo=(tBox.maxy+tBox.miny) /2;

				a=(tBox.maxx-tBox.minx) /2;
				b=(tBox.maxy-tBox.miny) /2;

				_a2= 1.0 / (a*a);
				_b2= 1.0 / (b*b);

				for(x=tBox.minx;x<tBox.maxx;x++)
					for(y=tBox.miny;y<tBox.maxy;y++) 
						if( ((x-xo)*(x-xo)*_a2+(y-yo)*(y-yo)*_b2) <=1.0) {
							if (OccModelImage[y*cols+x]) {
								OccList[OccListCount]=y*cols+x;
								OccListCount++;
							}
							OccModelImage[y*cols+x]=l;
						}
			}

		}

		
	}

}


void GenerateOccModel(MovingObject ** ObjectsInMergeList,
											unsigned int ObjectsInMerge_cnt,
											unsigned int maxdepthindex,
											unsigned char * OccModelImage,
											unsigned int rows,
											unsigned int cols)
{

	Box tBox;
	unsigned int d,k;
	int l;
	int x,y;
	double a,b,_a2,_b2;
	int xo,yo;

	memset(OccModelImage,0,rows*cols);

	for (d=1;d<=maxdepthindex;d++)  {
		for (k=0;k<ObjectsInMerge_cnt;k++){
			if (ObjectsInMergeList[k]->OccDepth==d) {
				tBox=ObjectsInMergeList[k]->BBox;
				l=ObjectsInMergeList[k]->OccLabel;

				xo=(tBox.maxx+tBox.minx) /2;
				yo=(tBox.maxy+tBox.miny) /2;

				a=(tBox.maxx-tBox.minx) /2;
				b=(tBox.maxy-tBox.miny) /2;

				_a2= 1.0 / (a*a);
				_b2= 1.0 / (b*b);

				for(x=tBox.minx;x<tBox.maxx;x++)
					for(y=tBox.miny;y<tBox.maxy;y++) 
						if( ((x-xo)*(x-xo)*_a2+(y-yo)*(y-yo)*_b2) <=1.0) {
							OccModelImage[y*cols+x]=l;
						}
			}

		}

		
	}

}


void UpdateVisibilityTable(unsigned char * label_image,
														unsigned char * OccModelImage,
													int Nobjects,
													double alpha,
													double ** PR,
													unsigned int rows,
													unsigned int cols)
{

	int C[3][3];
	int SM[3];

	//double PR[3][3];
	int N,i,j;
	unsigned int k;
	unsigned char L,M;

	double beta;

	N=3;  // Nobjects including background.

	memset(C,0,N*N*sizeof(int));
	memset(SM,0,N*sizeof(int));


	for (k=0;k<rows*cols;k++) {
		L=label_image[k];
		M=OccModelImage[k];

		if (L >0 ) {
			C[L][M] ++;
			SM[M] ++;		
		}

	}

	beta= (1-alpha) / (double) N;

	for (i=0;i<N;i++)
		for(j=0;j<N;j++) {
			PR[i][j] = alpha * ((double) C[i][j] / (double) SM[j]) + beta;
		}

	

}

void UpdateVisibility(unsigned char * OccModelImage,
											double ** V,
											double ** PR,
											unsigned int rows,
											unsigned int cols)
{

	unsigned int k;
	unsigned char M;
	int i,N;

	N=3;

	for (k=0;k<rows*cols;k++) {
		M=OccModelImage[k];

		for (i=0;i<N;i++)
			V[i][k]=PR[i][M];

		
	}


}


void UpdateVisibilityPr(unsigned char * label_image,
												unsigned char * OccModelImage,
												double ** V,
												int Nobjects,
												double alpha,
												unsigned int rows,
												unsigned int cols)


{

	int C[3][3];
	int SM[3];

	double PR[3][3];
	int N,i,j;
	unsigned int k;
	unsigned char L,M;

	double beta;

	N=3;  // Nobjects including background.

	memset(C,0,N*N*sizeof(int));
	memset(SM,0,N*sizeof(int));


	for (k=0;k<rows*cols;k++) {
		L=label_image[k];
		M=OccModelImage[k];

		if (L >0 ) {
			C[L][M] ++;
			SM[M] ++;		
		}

	}

	beta= (1-alpha) / (double) N;

	for (i=0;i<N;i++)
		for(j=0;j<N;j++) {
			PR[i][j] = alpha * ((double) C[i][j] / (double) SM[j]) + beta;
		}

	for (k=0;k<rows*cols;k++) {
		M=OccModelImage[k];

		for (i=0;i<N;i++)
			V[i][k]=PR[i][M];

		
	}


}

void EvaluateOccModel(int * object_label_image,
											int ObjectLabel,
											int * OccList,
											int OccListCount,
											unsigned char * seg_label_image,
											unsigned int rows,
											unsigned int cols,
											double & e1,
											double & e2) 
{

	int i;
	int n=0;
	int L1=0,L2=0,E=0;
	int idx;

	for (i=0;i<OccListCount;i++) {
		idx=OccList[i];
		assert(idx >0 && idx < rows*cols);
		if (object_label_image[idx]==ObjectLabel){
			n++;
			if (seg_label_image[idx]==1) {
				L1++;
			}
			else if (seg_label_image[idx]==2) {
				L2++;
			}
			else
				E++;
		}
		
	}

	e1=(n > 0) ? (100*L1) / n : 0;
	e2=(n > 0) ? (100*L2) / n : 0;


}



void ModifyProfiles(MovingObject * objectptr,
										int profilelength,
										double gama)
{
	double gama1=(1-gama)/3;

	for( int y=0;y<profilelength;y++) {
		objectptr->HeadModel->YProfile[y]= gama * objectptr->HeadModel->YProfile[y] + gama1;
		objectptr->TopModel->YProfile[y]= gama * objectptr->TopModel->YProfile[y] + gama1;
		objectptr->BottomModel->YProfile[y]= gama * objectptr->BottomModel->YProfile[y] + gama1;
																			

		
	}


}


double EvaluateHypothesis(int x1,int x2,int y1,int y2,
												Box evalbox,
												int * object_label_image,
												int ObjectLabel,
												MovingObject ** ObjectInMergeList,
												double * C1,double * C2, 
												double ** V,
												double * Px1,double * Px2,
												double * P1, double * P2,
												unsigned char * part_label_image,
												double OutlierTh,
												unsigned int rows,
												unsigned int cols)
{

	unsigned int miny,maxy,minx,maxx;

	miny=evalbox.miny;
	minx=evalbox.minx;
	maxy=evalbox.maxy;
	maxx=evalbox.maxx;

	double * Phead1 =ObjectInMergeList[0]->HeadModel->Prob_image ;
	double * Phead2 = ObjectInMergeList[1]->HeadModel->Prob_image;
	double * Ptop1 = ObjectInMergeList[0]->TopModel->Prob_image;
	double * Ptop2 = ObjectInMergeList[1]->TopModel->Prob_image;
	double * Pbottom1 = ObjectInMergeList[0]->BottomModel->Prob_image;
	double * Pbottom2 = ObjectInMergeList[1]->BottomModel->Prob_image;


	double sum1,sum2,N1,N2;
	double sume1,sume2;
	double eps=1e-40;
	double logeps=-92.1;
	double logoutlierth=log(OutlierTh);
	double px;

	int disp1,disp2;
	int x,y;
	int idx;



	sum1=0;
	sum2=0;
	sume1=0;
	sume2=0;
	N1=0;
	N2=0;
	disp1=rows-(ObjectInMergeList[0]->top_bottom_sep_disp+1);
	disp2=rows-(ObjectInMergeList[1]->top_bottom_sep_disp+1);

	memset(part_label_image,0,rows*cols);

	for (y=miny; y < maxy; y++)
		for (x=minx ; x<maxx ; x++) 
			if (object_label_image[y*cols+x]==ObjectLabel) {
				idx=y*cols+x;
				if ( C1[disp1+y1-y] == 0)
					P1[idx]=eps;
				else
					P1[idx]=((1*ObjectInMergeList[0]->HeadModel->YProfile[disp1+y1-y]*Phead1[idx]+
									1*ObjectInMergeList[0]->TopModel->YProfile[disp1+y1-y]*Ptop1[idx] +
									1*ObjectInMergeList[0]->BottomModel->YProfile[disp1+y1-y]*Pbottom1[idx]) 
									/ C1[disp1+y1-y]) *
									Px1[abs(x-x1)] ;
									 // * V[1][idx];
				
				if ( C2[disp2+y2-y] ==0)
					P2[idx]=eps;
				else
					P2[idx]=((1*ObjectInMergeList[1]->HeadModel->YProfile[disp2+y2-y]*Phead2[idx]+
									1*ObjectInMergeList[1]->TopModel->YProfile[disp2+y2-y]*Ptop2[idx] +
									1*ObjectInMergeList[1]->BottomModel->YProfile[disp2+y2-y]*Pbottom2[idx]) 
									/ C2[disp2+y2-y]) *
									Px2[abs(x-x2)] ;
								  // * V[2][idx];
				

				px=P1[idx]+P2[idx];

				if (P1[idx] > P2[idx]) {
					if (P1[idx] <= eps){
						P1[idx]=eps;
						sum1+=logeps;
						sume1+=1.0;
						part_label_image[y*cols+x]=4;
					}
					else if (P1[idx] <= OutlierTh) {
						part_label_image[y*cols+x]=4;
						//if (P2[idx] <= eps) {
						//	sum1+=logeps;
						//}
						//else sum1+=log(P2[idx]);
						sum1+=log(P1[idx]);
						sume1+=1-P1[idx]/px;
					}
					else {
						part_label_image[y*cols+x]=10;
						//blob_label_image[idx]=P1label;
						//seg_label_image[idx]=1;
						//if (P2[idx] <= eps) {
						//	sum1+=logeps;
						//}
						//else sum1+=log(P2[idx]);
						sum1+=log(P1[idx]);
						sume1+=1-P1[idx]/px;
					}
					//sum1+=(P1[idx]);
					N1++;
				}
				else {
					if (P2[idx] <= eps) {
						P2[idx]=eps;
						sum2+=logeps;
						sume2+=1.0;
						part_label_image[y*cols+x]=4;
					}
					else if (P2[idx] <= OutlierTh) {
						part_label_image[y*cols+x]=4;

						//if (P1[idx] <= eps) {
						//	sum2+=logeps;
						//}
						//else	sum2+=log(P1[idx]);
						sum2+=log(P2[idx]);
						sume2+=1-P2[idx]/px;

					}
					else {
						part_label_image[y*cols+x]=11;
						//blob_label_image[idx]=P2label;
						//seg_label_image[idx]=2;
						//if (P1[idx] <= eps) {
						//	sum2+=logeps;
						//}
						//else	sum2+=log(P1[idx]);
						sum2+=log(P2[idx]);
						sume2+=1-P2[idx]/px;

					}
					//sum2+=(P2[idx]);
					N2++;
				}

				

			}
	

			return -(sum1 / N1 + sum2 / N2);
			//return -(sum1+sum2)/(N1+N2);
			//return (sume1+sume2) / (N1+N2);

}


void ObjectTracker::MergeObjectColorSegmentII(unsigned int ObjectLabel,
							 unsigned char * image,
							 unsigned char * rgbimage,
							 int * object_label_image,
							 unsigned char * part_label_image,
							 KernelTableStruct * KernelTablePtr,
							 MovingObject * MergeObject,
							 MovingObject ** ObjectInMergeList,
							 unsigned int ObjectInMerge_cnt,
							 int OptFlag)
{
	int x,y;
	unsigned int i;
	unsigned int miny,maxy,minx,maxx;

	unsigned char * ptr;

	unsigned int FGPixelsCount;

	double OutlierTh = 1e-10; // 1e-8
	double OutlierTh2 = 1e-6;

	if (ObjectInMerge_cnt !=2) 
		return;

	double * Phead1 = &membuffer[0];
	double * Phead2 = &membuffer[1*rows*cols];
	double * Ptop1 = &membuffer[2*rows*cols];
	double * Ptop2 = &membuffer[3*rows*cols];
	double * Pbottom1 = &membuffer[4*rows*cols];
	double * Pbottom2 = &membuffer[5*rows*cols];

	double * P1 = &membuffer[6*rows*cols];
	double * P2 = &membuffer[7*rows*cols];

	double * V[3];

	V[0] = &VisibilityPr[0*rows*cols];
	V[1] = &VisibilityPr[1*rows*cols];
	V[2] = &VisibilityPr[2*rows*cols];

	double * C1 = new double [rows*2];
	double * C2 = new double [rows*2];

	double * Px1 = new double [cols];
	double * Px2 = new double [cols];
	

	int * blob_label_image=BlobLabelImage;
	unsigned char * part_label_imageI=PartLabelImageI;
	unsigned char * part_label_imageII=PartLabelImageII;
	unsigned char * seg_label_image=OccLabelImage;

	ColorBlob * P1Blob;
	ColorBlob * P2Blob;

	int P1label,P2label;
	int P1objectlabel,P2objectlabel;

	/*
	P1Blob=new ColorBlob(rows,cols);
	P2Blob=new ColorBlob(rows,cols);

	P1Blob->BlobID=ColorBlobCount;
	P1label=ColorBlobCount;
	ColorBlobList[ColorBlobCount++]=P1Blob;

	P2Blob->BlobID=ColorBlobCount;
	P2label=ColorBlobCount;
	ColorBlobList[ColorBlobCount++]=P2Blob;

	*/

	P1label=ObjectInMergeList[0]->BlobLabel;
	P2label=ObjectInMergeList[1]->BlobLabel;

	P1Blob=ColorBlobList[P1label];
	P2Blob=ColorBlobList[P2label];



	P1objectlabel=ObjectInMergeList[0]->ObjectLabel;
	P2objectlabel=ObjectInMergeList[1]->ObjectLabel;


	memset(blob_label_image,0,rows*cols*sizeof(int));

	miny=MergeObject->BBox.miny;
	minx=MergeObject->BBox.minx;
	maxy=MergeObject->BBox.maxy;
	maxx=MergeObject->BBox.maxx;




	for(i=0;i<rows*cols;i++) {
		Phead1[i]=0;
		Phead2[i]=0;
		Ptop1[i]=0;
		Ptop2[i]=0;
		Pbottom1[i]=0;
		Pbottom2[i]=0;
		P1[i]=0;
		P2[i]=0;
	}


	// evaluate part color probabilities in batches
	// fill in FG pixels data structure

	unsigned int Ntargets=0;
  int tmp;
	int c1_index,c2_index,c3_index;
	int rgbindex;

	FGPixelsCount=0;

	// reset rgbLUT
	for(i=0;i<256*256*256;i++) 
		rgbLUT[i]=-1;


	for (y=miny; y < maxy; y++) {
		for (x=minx ; x<maxx ; x++)
			if (object_label_image[y*cols+x]==ObjectLabel) 
			{
				FGPixelLocX[FGPixelsCount]=x;
				FGPixelLocY[FGPixelsCount]=y;


				tmp=(y*cols+x)*3;
				c1_index=rgbimage[tmp];
				c2_index=rgbimage[tmp+1];
				c3_index=rgbimage[tmp+2];

				rgbindex=c1_index*65536+c2_index*256+c3_index;

				if (rgbLUT[rgbindex]==-1)
				{
					FGtargets[Ntargets*3]=srgImage[tmp];
					FGtargets[Ntargets*3+1]=srgImage[tmp+1];
					FGtargets[Ntargets*3+2]=srgImage[tmp+2];

					FGPixelMap[FGPixelsCount]=Ntargets;
					rgbLUT[rgbindex]=Ntargets;
					Ntargets++;
				}
				else
				{
					FGPixelMap[FGPixelsCount]=rgbLUT[rgbindex];
				}

				//FGPixelValue[FGPixelsCount*3]=srgImage[3*(y*cols+x)];
				//FGPixelValue[FGPixelsCount*3+1]=srgImage[3*(y*cols+x)+1];
				//FGPixelValue[FGPixelsCount*3+2]=srgImage[3*(y*cols+x)+2];

				FGPixelsCount++;
			}
	}




		
	/*
	// write targets (foreground values) to file
	{
		FILE * targetsFP;
		char tfname [30];

		sprintf(tfname,"targets%d.txt",FrameNo);
		targetsFP=fopen(tfname,"wb");
		fwrite(FGPixelValue,sizeof(double),FGPixelsCount*3,targetsFP);
		fclose(targetsFP);

	}
	*/



	LARGE_INTEGER pht1,pht2,phf;
	float t1;
	float ct1,ct2,ct3,ct4,ct5,ct6;

	QueryPerformanceCounter(&pht1);


	ct1=ObjectInMergeList[0]->HeadModel->EstimateProbabilityFG(FGPixelsCount,FGPixelLocX,FGPixelLocY,FGPixelMap,FGtargets,NULL,Ntargets,Phead1);
	ct2=ObjectInMergeList[1]->HeadModel->EstimateProbabilityFG(FGPixelsCount,FGPixelLocX,FGPixelLocY,FGPixelMap,FGtargets,NULL,Ntargets,Phead2);
	ct3=ObjectInMergeList[0]->TopModel->EstimateProbabilityFG(FGPixelsCount,FGPixelLocX,FGPixelLocY,FGPixelMap,FGtargets,NULL,Ntargets,Ptop1);
	ct4=ObjectInMergeList[1]->TopModel->EstimateProbabilityFG(FGPixelsCount,FGPixelLocX,FGPixelLocY,FGPixelMap,FGtargets,NULL,Ntargets,Ptop2);
	ct5=ObjectInMergeList[0]->BottomModel->EstimateProbabilityFG(FGPixelsCount,FGPixelLocX,FGPixelLocY,FGPixelMap,FGtargets,NULL,Ntargets,Pbottom1);
	ct6=ObjectInMergeList[1]->BottomModel->EstimateProbabilityFG(FGPixelsCount,FGPixelLocX,FGPixelLocY,FGPixelMap,FGtargets,NULL,Ntargets,Pbottom2);


	QueryPerformanceCounter(&pht2);
	QueryPerformanceFrequency(&phf);
	t1=(float) (1000.0*(float)(pht2.QuadPart - pht1.QuadPart) / phf.QuadPart); 



	/*
	fprintf(FP,"%d : %d %d %f %d %d %d %d %d %d\n",FrameNo,FGPixelsCount,Ntargets,t1,
								ObjectInMergeList[0]->HeadModel->SubSampleCount,
								ObjectInMergeList[1]->HeadModel->SubSampleCount,
								ObjectInMergeList[0]->TopModel->SubSampleCount,
								ObjectInMergeList[1]->TopModel->SubSampleCount,
								ObjectInMergeList[0]->BottomModel->SubSampleCount,
								ObjectInMergeList[1]->BottomModel->SubSampleCount);

	*/
	fprintf(FP," %d %d %d %f    %f %f %f %f %f %f \n",FrameNo,FGPixelsCount,Ntargets,t1,ct1,ct2,ct3,ct4,ct5,ct6);

	/*
	// evaluate part color probabilities
	
	for (y=miny; y < maxy; y++)
		for (x=minx ; x<maxx ; x++)
			if (object_label_image[y*cols+x]==ObjectLabel) {
				ptr=image+3*(y*cols+x);
				Phead1[y*cols+x]=ObjectInMergeList[0]->HeadModel->EstimateProbability(ptr,KernelTablePtr);
				Phead2[y*cols+x]=ObjectInMergeList[1]->HeadModel->EstimateProbability(ptr,KernelTablePtr);
				Ptop1[y*cols+x]=ObjectInMergeList[0]->TopModel->EstimateProbability(ptr,KernelTablePtr);
				Ptop2[y*cols+x]=ObjectInMergeList[1]->TopModel->EstimateProbability(ptr,KernelTablePtr);
				Pbottom1[y*cols+x]=ObjectInMergeList[0]->BottomModel->EstimateProbability(ptr,KernelTablePtr);
				Pbottom2[y*cols+x]=ObjectInMergeList[1]->BottomModel->EstimateProbability(ptr,KernelTablePtr);
			}
	*/


	ObjectInMergeList[0]->HeadModel->Prob_image=Phead1;
	ObjectInMergeList[1]->HeadModel->Prob_image=Phead2;
	
	ObjectInMergeList[0]->TopModel->Prob_image=Ptop1;
	ObjectInMergeList[1]->TopModel->Prob_image=Ptop2;

	ObjectInMergeList[0]->BottomModel->Prob_image=Pbottom1;
	ObjectInMergeList[1]->BottomModel->Prob_image=Pbottom2;




  int h1,h2,w1,w2;
	int x1,x2,y1,y2;
	int idx;


	//h1=ObjectInMergeList[0]->BBox.maxy-ObjectInMergeList[0]->BBox.miny;
	//h2=ObjectInMergeList[1]->BBox.maxy-ObjectInMergeList[1]->BBox.miny;

	//w1=ObjectInMergeList[0]->BBox.maxx-ObjectInMergeList[0]->BBox.minx;
	//w2=ObjectInMergeList[1]->BBox.maxx-ObjectInMergeList[1]->BBox.minx;

	h1=ObjectInMergeList[0]->AveHeight;
	h2=ObjectInMergeList[1]->AveHeight;

	w1=ObjectInMergeList[0]->AveWidth;
	w2=ObjectInMergeList[0]->AveWidth;



	double segma1,segma2,S1,S2,t;
	double sum1,sum2,N1,N2;
	double eps=1e-40;
	double logeps=-92.1;

	segma1=(double) w1/(2.0*1.65);
	segma2=(double) w2/(2.0*1.65);
	S1=0;
	S2=0;
	for (x=0; x < cols;x++) {
		t=(double) x/segma1;
		Px1[x]=exp(-t*t);
		S1+=Px1[x];

		t=(double) x/segma2;
		Px2[x]=exp(-t*t);
		S2+=Px2[x];
	}

	for (x=0; x < cols;x++) {
		Px1[x] /=S1;
		Px2[x] /=S2;

	}

	// double gama=0.999;

	// add ncertainity to object Profiles
	//ModifyProfiles(ObjectInMergeList[0],rows*2,gama);
	//ModifyProfiles(ObjectInMergeList[1],rows*2,gama);

	// C1 , C2 are normalization constants for person 1 and person 2
	for (y=0; y < 2*rows; y++) {

		C1[y]=ObjectInMergeList[0]->HeadModel->YProfile[y]+
					ObjectInMergeList[0]->TopModel->YProfile[y]+
					ObjectInMergeList[0]->BottomModel->YProfile[y];

		C2[y]=ObjectInMergeList[1]->HeadModel->YProfile[y]+
					ObjectInMergeList[1]->TopModel->YProfile[y]+
					ObjectInMergeList[1]->BottomModel->YProfile[y];

	}

	// test 
	
	x1=ObjectInMergeList[0]->X;
	y1=ObjectInMergeList[0]->BBox.maxy;
	x2=ObjectInMergeList[1]->X;
	y2=ObjectInMergeList[1]->BBox.maxy;
	
	int x1min,x1max,x2min,x2max;
	int y1min,y1max,y2min,y2max;
	int xstep, int ystep;
	double minerror=1e60;
	double maxerror=0.0;
	double errorsum;
	int disp1,disp2;
	int w;
	//int x1opt,x2opt,y1opt,y2opt;


	xstep=2;
	ystep=2;
	w=4;

	x1min= x1-w;
	x1max= x1+w;
	y1min= y1-w;
	y1max= y1+w;

	x2min= x2-w;
	x2max= x2+w;
	y2min= y2-w;
	y2max= y2+w;


	memset(seg_label_image,0,rows*cols);

	x1opt=x1;
	x2opt=x2;
	y1opt=y1;
	y2opt=y2;
	
	

	
	// search

	if (OptFlag)	{

	
		for (x1=x1min;x1<=x1max;x1+=xstep)
			for (y1=y1min;y1<=y1max;y1+=ystep) 
				for (x2=x2min;x2<=x2max;x2+=xstep)
					for (y2=y2min;y2<=y2max;y2+=ystep)  
		{
				
		
			//GenerateOccModel(ObjectInMergeList,ObjectInMerge_cnt,2,OccModelImage,rows,cols);
			//UpdateVisibility(OccModelImage,V,VisPrTable,rows,cols);


			errorsum=EvaluateHypothesis(x1,x2,y1,y2,MergeObject->BBox,
							object_label_image,ObjectLabel,ObjectInMergeList,C1,C2,V,Px1,Px2,P1,P2,
							part_label_image,OutlierTh,rows,cols);

			if (errorsum < minerror) {
					x1opt=x1;
					x2opt=x2;
					y1opt=y1;
					y2opt=y2;
					minerror=errorsum;
					maxerror=errorsum;

			}
					
		}		
	
	}

	//GenerateOccModel(ObjectInMergeList,ObjectInMerge_cnt,2,OccModelImage,rows,cols);
	//UpdateVisibility(OccModelImage,V,VisPrTable,rows,cols);


	errorsum=EvaluateHypothesis(x1opt,x2opt,y1opt,y2opt,MergeObject->BBox,
						object_label_image,ObjectLabel,ObjectInMergeList,C1,C2,V,Px1,Px2,P1,P2,
						part_label_image,OutlierTh,rows,cols);

			 memset(part_label_imageI,0,rows*cols*sizeof(unsigned char));

			 FilterII(object_label_image,part_label_image,ObjectLabel,
				10,11,part_label_imageI,rows,cols);

			 memset(part_label_image,0,rows*cols*sizeof(unsigned char));
			 
			 FilterI(object_label_image,part_label_imageI,P1,P2,ObjectLabel,
				10,11,part_label_image,rows,cols);

		// label blob image
		for (i=0;i<rows*cols;i++) {
			if (part_label_image[i]==10) {
				blob_label_image[i]=P1label;
				seg_label_image[i]=1;

			}
			else if (part_label_image[i]==11) {
				blob_label_image[i]=P2label;
				seg_label_image[i]=2;
			}
			else {
				blob_label_image[i]=0;
				seg_label_image[i]=0;
			}
		}


			// extract blob profiles
			ExtractBlobProfiles(blob_label_image);
			P1Blob->UpdateAttributes(1.0,1);
			P2Blob->UpdateAttributes(1.0,1);

			//P1Blob->BBox.miny=P1Blob->BBox.maxy-h1;
			//P2Blob->BBox.miny=P2Blob->BBox.maxy-h2;

			//ObjectInMergeList[0]->BBox=P1Blob->BBox;
			//ObjectInMergeList[1]->BBox=P2Blob->BBox;

			ObjectInMergeList[0]->BBox.minx=P1Blob->medianX - w1 / 2;
			ObjectInMergeList[0]->BBox.maxx=P1Blob->medianX + w1 / 2;
			ObjectInMergeList[0]->BBox.miny=P1Blob->BBox.maxy - h1 ;
			ObjectInMergeList[0]->BBox.maxy=P1Blob->BBox.maxy;

			ObjectInMergeList[1]->BBox.minx=P2Blob->medianX - w2 / 2;
			ObjectInMergeList[1]->BBox.maxx=P2Blob->medianX + w2 / 2;
			ObjectInMergeList[1]->BBox.miny=P2Blob->BBox.maxy - h2 ;
			ObjectInMergeList[1]->BBox.maxy=P2Blob->BBox.maxy;



			ObjectInMergeList[0]->X=P1Blob->medianX;
			ObjectInMergeList[1]->X=P2Blob->medianX;


			if (OptFlag)	{
			
				ObjectInMergeList[0]->X=x1opt;
				ObjectInMergeList[1]->X=x2opt;

				ObjectInMergeList[0]->BBox.maxy=y1opt;
				ObjectInMergeList[1]->BBox.maxy=y2opt;
			
			}

			//ObjectInMergeList[0]->Y=y1opt;
			//ObjectInMergeList[1]->Y=y2opt;

			ObjectInMergeList[0]->nWPixels=P1Blob->NPixels;
			ObjectInMergeList[1]->nWPixels=P2Blob->NPixels;


			//determine Occ depth.
			double e1,e2;

			ObjectInMergeList[0]->VisibilityRatio=(100*ObjectInMergeList[0]->nWPixels) / (double) ObjectInMergeList[0]->AveNPixels;
			ObjectInMergeList[1]->VisibilityRatio=(100*ObjectInMergeList[1]->nWPixels) / (double) ObjectInMergeList[1]->AveNPixels;

			if (ObjectInMergeList[0]->VisibilityRatio >
					ObjectInMergeList[1]->VisibilityRatio )
			{
				ObjectInMergeList[0]->OccDepth=2;
				ObjectInMergeList[1]->OccDepth=1;			
				
			}
			else {
				ObjectInMergeList[0]->OccDepth=1;
				ObjectInMergeList[1]->OccDepth=2;
			}


			//	ObjectInMergeList[0]->OccDepth=1;
			//	ObjectInMergeList[1]->OccDepth=2;
			
			ObjectInMergeList[0]->OccLabel=1;
			ObjectInMergeList[1]->OccLabel=2;
			
			//GenerateOccModelOverlap(ObjectInMergeList,ObjectInMerge_cnt,2,OccModelImage,OccList,OccListCount,rows,cols);


			//EvaluateOccModel(object_label_image,ObjectLabel,OccList,OccListCount,
			//	seg_label_image,rows,cols,e1,e2);

			//ObjectInMergeList[0]->VisibilityRatio=e1;
			//ObjectInMergeList[1]->VisibilityRatio=e2;


			//x1opt=P1Blob->medianX;
			//y1opt=P1Blob->BBox.maxy;
			//x2opt=P2Blob->medianX;
			//y2opt=P2Blob->BBox.maxy;


			GenerateOccModel(ObjectInMergeList,ObjectInMerge_cnt,2,OccModelImage,rows,cols);
			UpdateVisibilityTable(seg_label_image,OccModelImage,3,0.9,VisPrTable,rows,cols);

			// update visibility probability
			//UpdateVisibilityPr(seg_label_image,OccModelImage,V,3,0,rows,cols);

			memset(part_label_imageII,0,rows*cols*sizeof(unsigned char));

	disp1=rows-(ObjectInMergeList[0]->top_bottom_sep_disp+1);
	disp2=rows-(ObjectInMergeList[1]->top_bottom_sep_disp+1);

	// segment parts
	double Ph,Pt,Pb;
	for (y=miny; y < maxy; y++)
		for (x=minx ; x<maxx ; x++) 
			if (object_label_image[y*cols+x]==ObjectLabel) {
				idx=y*cols+x;
			  if (blob_label_image[idx]==P1label){
				//if (part_label_imageI[idx]==10) {
					Ph=(ObjectInMergeList[0]->HeadModel->YProfile[disp1+y1opt-y]*Phead1[idx]) / C1[disp1+y1opt-y];
					Pt=(ObjectInMergeList[0]->TopModel->YProfile[disp1+y1opt-y]*Ptop1[idx]) / C1[disp1+y1opt-y];
					Pb=(ObjectInMergeList[0]->BottomModel->YProfile[disp1+y1opt-y]*Pbottom1[idx]) / C1[disp1+y1opt-y];

					if (Ph > Pt)
						if (Ph > Pb) {
							// head
							if (Ph > OutlierTh2) {
								part_label_imageII[idx]=3;
							}
							else {
							}
						} else {
							// bottom
							if (Pb > OutlierTh2) {
								part_label_imageII[idx]=2;
							}
							else {
							}


						}
					else 
						if (Pt > Pb) {
							// top
							if (Pt > OutlierTh2) {
								part_label_imageII[idx]=1;
							}
							else {
							}

						} else {
							// bottom
							if (Pb > OutlierTh2) {
								part_label_imageII[idx]=2;
							}
							else {
							}
						}		
				}
		  else if (blob_label_image[idx]==P2label){
			//else if (part_label_imageI[idx]==11) {

					Ph=(ObjectInMergeList[1]->HeadModel->YProfile[disp2+y2opt-y]*Phead2[idx]) / C2[disp2+y2opt-y];
					Pt=(ObjectInMergeList[1]->TopModel->YProfile[disp2+y2opt-y]*Ptop2[idx]) / C2[disp2+y2opt-y];
					Pb=(ObjectInMergeList[1]->BottomModel->YProfile[disp2+y2opt-y]*Pbottom2[idx]) / C2[disp2+y2opt-y];

					if (Ph > Pt)
						if (Ph > Pb) {
							// head
							if (Ph > OutlierTh2) {
								part_label_imageII[idx]=13;
							}
							else {
							}
						} else {
							// bottom
							if (Pb > OutlierTh2) {
								part_label_imageII[idx]=12;
							}
							else {
							}


						}
					else 
						if (Pt > Pb) {
							// top
							if (Pt > OutlierTh2) {
								part_label_imageII[idx]=11;
							}
							else {
							}

						} else {
							// bottom
							if (Pb > OutlierTh2) {
								part_label_imageII[idx]=12;
							}
							else {
							}
						}
						







				}
			}
		

	// relabel object image
	for (i=0;i<rows*cols;i++)
		if (part_label_imageI[i]==10)
			object_label_image[i]=P1objectlabel;
		else if (part_label_imageI[i]==11)
			object_label_image[i]=P2objectlabel;
		else object_label_image[i]=0;




	
			//sum=sum1 / (double) N1 + sum2 / (double) N2;
			/*
			sum=(sum1+sum2) / (double) (N1 + N2);
			if (sum > maxsum) {
				maxsum=sum;
				x1opt=x1;
				y1opt=y1;
				x2opt=x2;
				y2opt=y2;
			}
			*/
	//}
	
	
	// update Hbbox

	/*
	ObjectInMergeList[0]->BBox.minx=x1opt-w1/2;
	ObjectInMergeList[0]->BBox.maxx=x1opt+w1/2;
	ObjectInMergeList[0]->BBox.miny=y1opt-h1;
	ObjectInMergeList[0]->BBox.maxy=y1opt;
	ObjectInMergeList[0]->X=x1opt;

	ObjectInMergeList[1]->BBox.minx=x2opt-w2/2;
	ObjectInMergeList[1]->BBox.maxx=x2opt+w2/2;
	ObjectInMergeList[1]->BBox.miny=y2opt-h2;
	ObjectInMergeList[1]->BBox.maxy=y2opt;
	ObjectInMergeList[1]->X=x2opt;

	*/

	delete [] C1;
	delete [] C2;
	delete [] Px1;
	delete [] Px2;

	// reset ColorBlobCount
	//ColorBlobCount=1;
	//delete P1Blob;
	//delete P2Blob;
}


void ObjectTracker::MergeObjectColorSegment(unsigned int ObjectLabel,
							 unsigned char * image,
							 int * object_label_image,
							 unsigned char * part_label_image,
							 KernelTableStruct * KernelTablePtr,
							 MovingObject * MergeObject,
							 MovingObject ** ObjectInMergeList,
							 unsigned int ObjectInMerge_cnt)
{
	unsigned int x,y;
	unsigned int miny,maxy,minx,maxx;
	double p_top1,p_top2;

	unsigned char * ptr;

	double OutlierTh = 0.0001;

	miny=MergeObject->BBox.miny;
	minx=MergeObject->BBox.minx;
	maxy=MergeObject->BBox.maxy;
	maxx=MergeObject->BBox.maxx;


	if (ObjectInMerge_cnt==2) {

	for (y=miny; y < maxy; y++)
		for (x=minx ; x<maxx ; x++)
			if (object_label_image[y*cols+x]==ObjectLabel) {
				ptr=image+3*(y*cols+x);

				p_top1=ObjectInMergeList[0]->TopModel->EstimateProbability(ptr,KernelTablePtr);
				p_top2=ObjectInMergeList[1]->TopModel->EstimateProbability(ptr,KernelTablePtr);

				if (p_top1 > p_top2)
					if (p_top1 > OutlierTh)
						part_label_image[y*cols+x]=10;
					else 
						part_label_image[y*cols+x]=5;

				else if (p_top2 > OutlierTh)
						part_label_image[y*cols+x]=11;
					else
						part_label_image[y*cols+x]=5;


			}

	}


}


void ObjectTracker::ForegroundSegmentation(unsigned char * image,
										   unsigned char * FGimage,
										   int * blob_label_image)

{

	unsigned char * ptr;
	double p;

	ColorBlob * blob;

	unsigned int i,x,y;

	double OutlierTh=0.0001;

	// reset Pr Image. & label image
	memset(PrImage,0,rows*cols*sizeof(double));
	memset(blob_label_image,0,rows*cols*sizeof(int));


	for (i=1;i<ColorBlobCount;i++)
		if (ColorBlobList[i] && ColorBlobList[i]->BlobState==ACTIVE_BLOB) {

			blob=ColorBlobList[i];

			for (y=blob->PBBox.miny;y<blob->PBBox.maxy;y++)
				for (x=blob->PBBox.minx;x<blob->PBBox.maxx;x++)
				{

					if (FGimage[y*cols+x]) {
						ptr=image+3*(y*cols+x);
						p=blob->colormodel->EstimateProbability(ptr,&KernelTable);

						if (p > OutlierTh && p > PrImage[y*cols+x] ) {
							PrImage[y*cols+x]=p;
							blob_label_image[y*cols+x]=i;
						}
					}
				}

		}
}

void InitializeVisibilityPr(double * V,
														unsigned int Nobjects,
														unsigned int imagesize)
{
	double _N= 1.0 / (double) Nobjects;

	for(unsigned int k=0;k<imagesize*Nobjects;k++)
		V[k]=_N;

}

void ObjectTracker::InitializeSegmentation(Merge * mergePtr)
{
	int i;

	// reset Visibility Pr
	InitializeVisibilityPr(VisibilityPr,mergePtr->Merge_cnt+1,rows*cols);

	VisPrTable[0][0]=0;
	VisPrTable[0][1]=0;
	VisPrTable[0][2]=0;

	VisPrTable[1][0]=0.5;
	VisPrTable[1][1]=0.5;
	VisPrTable[1][2]=0.5;

	VisPrTable[2][0]=0.5;
	VisPrTable[2][1]=0.5;
	VisPrTable[2][2]=0.5;

	// smooth profiles
	for (i=0;i<mergePtr->Merge_cnt;i++) {

		mergePtr->ObjectsInMergeList[i]->HeadModel->SmoothProfile(1,&KernelTable);
		mergePtr->ObjectsInMergeList[i]->TopModel->SmoothProfile(1,&KernelTable);
		mergePtr->ObjectsInMergeList[i]->BottomModel->SmoothProfile(1,&KernelTable);


	}



	// initialize blobs

	ColorBlob * pBlob;
	
	for (i=0;i<mergePtr->Merge_cnt;i++) {
		
		// 
		pBlob=new ColorBlob(rows,cols,3);

		pBlob->BlobID=ColorBlobCount;
		mergePtr->ObjectsInMergeList[i]->BlobLabel=ColorBlobCount;
		ColorBlobList[ColorBlobCount++]=pBlob;

	}



		


}

