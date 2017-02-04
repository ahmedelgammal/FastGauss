// FGSegmentation.cpp: implementation of the FGSegmentation class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Tracker.h"
#include "FGSegmentation.h"
#include <assert.h>
#include <math.h>

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FGSegmentation::FGSegmentation(int Rows,
															 int Cols)
{

	rows=Rows;
	cols=Cols;
	imagesize=rows*cols;


	logFP=fopen("log.txt","w");

}

int FGSegmentation::Initialize()
{
	
	membuffer_cnt=MAX_PEOPLE_IN_OCCLUSION*MAX_BLOBS_PER_PERSON;
  membuffer= new double[rows*cols*membuffer_cnt]; 


	blob_label_image = new unsigned char [rows*cols];
	//blob_label_image = (unsigned char *) malloc(rows*cols);

	residual_image = new unsigned char [rows*cols];


	// Multi Object Representation
	TrackedObjects = new MultiObject(rows,cols);

	TrackedObjectList = TrackedObjects->TrackedObjectList; 
	ObjectCount=0;
	Target_prob_image=TrackedObjects->Target_prob_image;
	Max_prob_image=TrackedObjects->Max_prob_image;
	object_label_image=TrackedObjects->object_label_image;

	
	// Foreground  data atructure

	FGPixels= new Foreground(rows,cols);

	// assign ptrs
	FGColorValues= FGPixels->FGColorValues;
	FGColorValuesQuantized= FGPixels->FGColorValuesQuantized;

	FGPixelLocX= FGPixels->FGPixelLocX;
	FGPixelLocY= FGPixels->FGPixelLocY;
	FGPixelMap= FGPixels->FGPixelMap;
	
	//rgbLUT = new int [256*256*256];
	//FGPixel_rgbindex = new long int [rows*cols];


	/*
	if ((FGColorValues=(double *) malloc(rows*cols*sizeof(double)*3))==NULL)
		return 0;

	if ((FGPixelLocX=(unsigned int *) malloc(rows*cols*sizeof(unsigned int)))==NULL)
		return 0;

	if ((FGPixelLocY=(unsigned int *) malloc(rows*cols*sizeof(unsigned int)))==NULL)
		return 0;

	if ((FGPixelMap=(unsigned int *) malloc(rows*cols*sizeof(unsigned int)))==NULL)
		return 0;

	if ((rgbLUT=(int *) malloc(256*256*256*sizeof(int)))==NULL)
		return 0;
	*/

	// Hypothesis
	

	return 1;
}

FGSegmentation::~FGSegmentation()
{

}


int FGSegmentation::AddObject(MovingObject * ObjPtr)
{

	TrackedObjects->AddObject(ObjPtr);

	ObjectCount=TrackedObjects->ObjectCount;

	return 1;
	
}

int FGSegmentation::AllocateBuffers()
{
	int i,k;

	// re-assign  buffers for color probability images

	assert(ObjectCount<= MAX_PEOPLE_IN_OCCLUSION);

	k=0;
	for(i=0;i<ObjectCount;i++) {
		TrackedObjectList[i]->HeadBlob->colormodel->Prob_image= & membuffer[(k++)*imagesize];
		TrackedObjectList[i]->TopBlob->colormodel->Prob_image= & membuffer[(k++)*imagesize];
		TrackedObjectList[i]->BottomBlob->colormodel->Prob_image= & membuffer[(k++)*imagesize];
	}


	// reset all probability buffers to 0
	memset(membuffer,0,imagesize*(k)*sizeof(double));

	return 1;
}




float FGSegmentation::EvaluateColorProbabilities()
{
	int i;

	float t=(float) 0.0;

	
	for(i=0;i<ObjectCount;i++) {
		t+=TrackedObjectList[i]->HeadBlob->colormodel->EstimateProbabilityFG(FGPixels->FGPixelsCount,
			FGPixelLocX,FGPixelLocY,FGPixelMap,FGColorValues,FGColorValuesQuantized,FGPixels->FGColorValuesCount,
			TrackedObjectList[i]->HeadBlob->colormodel->Prob_image);

		t+=TrackedObjectList[i]->TopBlob->colormodel->EstimateProbabilityFG(FGPixels->FGPixelsCount,
			FGPixelLocX,FGPixelLocY,FGPixelMap,FGColorValues,FGColorValuesQuantized,FGPixels->FGColorValuesCount,
			TrackedObjectList[i]->TopBlob->colormodel->Prob_image);

		t+=TrackedObjectList[i]->BottomBlob->colormodel->EstimateProbabilityFG(FGPixels->FGPixelsCount,
			FGPixelLocX,FGPixelLocY,FGPixelMap,FGColorValues,FGColorValuesQuantized,FGPixels->FGColorValuesCount,
			TrackedObjectList[i]->BottomBlob->colormodel->Prob_image);

			
	}

	return t;

}



void FGSegmentation::BlobSegmentation(unsigned char * FG_label_image,
																			unsigned char * label_image)
{
	int objectlabel,i;
	double Ph,Pt,Pb;
	unsigned int x,y;
	int idx;


	memset(blob_label_image,0,rows*cols);

	for(i=0;i<ObjectCount;i++) {
		TrackedObjectList[i]->HeadBlob->ResetProfiles();
		TrackedObjectList[i]->TopBlob->ResetProfiles();
		TrackedObjectList[i]->BottomBlob->ResetProfiles();		
	}


	for (y=miny; y < maxy; y++) {
		for (x=minx ; x<maxx ; x++){
			idx=y*cols+x;
			if (FG_label_image[idx]==FGlabel) 
			{	
				objectlabel = label_image[idx];
				if (objectlabel > 0) {

					assert (objectlabel <= ObjectCount);
					i=objectlabel-1;

					// evaluate part probabilities.

					Ph=TrackedObjectList[i]->HeadBlob->Y_prob[y]*
						TrackedObjectList[i]->HeadBlob->X_prob[x]*
						TrackedObjectList[i]->HeadBlob->colormodel->Prob_image[idx];

					Pt=TrackedObjectList[i]->TopBlob->Y_prob[y]*
						TrackedObjectList[i]->TopBlob->X_prob[x]*
						TrackedObjectList[i]->TopBlob->colormodel->Prob_image[idx];
					
					Pb=TrackedObjectList[i]->BottomBlob->Y_prob[y]*
						TrackedObjectList[i]->BottomBlob->X_prob[x]*
						TrackedObjectList[i]->BottomBlob->colormodel->Prob_image[idx]; 
							
					// find maximum and label parts
					if (Ph > Pt)
						if (Ph > Pb) {
							// head
							if (Ph > BlobOutlierTh) {
								blob_label_image[idx]=HEAD_BLOB;
								TrackedObjectList[i]->HeadBlob->AddPixel(x,y);
							}
							else {
							}
						} else {
							// bottom
							if (Pb > BlobOutlierTh) {
								blob_label_image[idx]=BOTTOM_BLOB;
								TrackedObjectList[i]->BottomBlob->AddPixel(x,y);
							}
							else {
							}
						}
					else 
						if (Pt > Pb) {
							// top
							if (Pt > BlobOutlierTh) {
								blob_label_image[idx]=TOP_BLOB;
								TrackedObjectList[i]->TopBlob->AddPixel(x,y);
							}
							else {
							}

						} else {
							// bottom
							if (Pb > BlobOutlierTh) {
								blob_label_image[idx]=BOTTOM_BLOB;
								TrackedObjectList[i]->BottomBlob->AddPixel(x,y);
							}
							else {
							}
						}						
				}


			}
		}
	}

}


void FGSegmentation::FGSubtraction(unsigned char * FG_label_image,
																	 unsigned char * FG_Residual_image)
{

	unsigned int x,y;
	int idx;

	double FGOutlierTH=1e-8;


	memset(FG_Residual_image,0,rows*cols);

	for (y=miny; y < maxy; y++) {
		for (x=minx ; x<maxx ; x++){
			idx=y*cols+x;
			if (FG_label_image[idx]==FGlabel) 
			{	
				if (Max_prob_image[idx] < FGOutlierTH)
					FG_Residual_image[idx]=255;
			}
		}

	}





}

void FGSegmentation::Segment(int FG_label,
												unsigned char * FG_label_image,
												double * srgImage,
												unsigned char * srgQImage,
												unsigned char * rgbImage)
{

	int i;
	float ct; // time for color probability estimation

	ObjectHypothesis JointHypothesis[MAX_PEOPLE_IN_OCCLUSION];


	// assign foreground dimension

	minx=0;
	maxx=cols;
	miny=0;
	maxy=rows;


	// FGlabel
	FGlabel=FG_label;



	// allocate buffer for color probability images
	AllocateBuffers();

	FGPixels->SetBounds(minx,maxx,miny,maxy);
	FGPixels->ExtractForegroundValues(FG_label_image,FG_label,srgImage,srgQImage,rgbImage);

	FGPixelsCount=FGPixels->FGPixelsCount;
	FGColorValuesCount=FGPixels->FGColorValuesCount;

	ct=EvaluateColorProbabilities();

	TrackedObjects->InitializeSpatialProbabilities();

	// hypothesis generation and evaluation
	{


	}

	// pick up best hypothesis



	{
		LARGE_INTEGER pht1,pht2,phf;
		float t1,t2;
		double lp;

		for(i=0;i<ObjectCount;i++) {
			JointHypothesis[i].X=TrackedObjectList[i]->state_X;
			JointHypothesis[i].Y=TrackedObjectList[i]->state_Y;
			JointHypothesis[i].h=TrackedObjectList[i]->state_scale;		
		}


		QueryPerformanceCounter(&pht1);


		// given best hypothesis
		TrackedObjects->EvaluateHypothesis(JointHypothesis,FGPixels);


		
		/*
		QueryPerformanceCounter(&pht2);
		QueryPerformanceFrequency(&phf);
		t1=(float) (1000.0*(float)(pht2.QuadPart - pht1.QuadPart) / phf.QuadPart); 


		QueryPerformanceCounter(&pht1);
		
		lp=EvaluateHypothesis(FG_label_image);

		QueryPerformanceCounter(&pht2);		
		QueryPerformanceFrequency(&phf);
		t2=(float) (1000.0*(float)(pht2.QuadPart - pht1.QuadPart) / phf.QuadPart); 


		fprintf(logFP," %d %d  %d  %f  %f   %f    %f    %f    %f    %f\n",FGPixelsCount,FGColorValuesCount,ObjectCount,lp,ct,t1,t2,t1+t2,pt1,pt2);
		*/

	}


	// segment blobs
	
	BlobSegmentation(FG_label_image,object_label_image);
	
	float alpha;

	for(i=0;i<ObjectCount;i++) {
		//  Detect Separators
		TrackedObjectList[i]->DetectSeparators();
		TrackedObjectList[i]->UpdateBlobAttributes(1);
		// update object state

		alpha=(float) 0.95;

		TrackedObjectList[i]->state_X=(int) ((TrackedObjectList[i]->HeadBlob->medianX+
																	TrackedObjectList[i]->TopBlob->medianX+
																	TrackedObjectList[i]->BottomBlob->medianX) / (double) 3.0);

		//TrackedObjectList[i]->state_Y=TrackedObjectList[i]->HeadBlob->BBox.maxy;
		TrackedObjectList[i]->state_Y=(int) (0.5* (((float) TrackedObjectList[i]->head_top_sep+
			 (TrackedObjectList[i]->head_top_sep_disp * TrackedObjectList[i]->state_scale)) +
			 ((float) TrackedObjectList[i]->top_bottom_sep+
			 (TrackedObjectList[i]->top_bottom_sep_disp * TrackedObjectList[i]->state_scale)))); 


		TrackedObjectList[i]->state_scale=alpha*TrackedObjectList[i]->state_scale+
																(1-alpha)*TrackedObjectList[i]->scale;


	}


	FGSubtraction(FG_label_image,residual_image);


	FGPixels->Reset();


}
