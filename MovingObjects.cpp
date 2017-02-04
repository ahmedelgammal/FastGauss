/*
	Implementation of Class MovingObject

*/


#include "MovingObjects.h"
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

int BoxIntersection(Box A,Box B)
{
	int x1,x2,y1,y2;

	if (A.maxx <= B.minx || B.maxx <= A.minx)
		return 0;

	if (A.maxy <= B.miny || B.maxy <= A.miny)
		return 0;

	// x1 = MAX ( A.minx,B.minx )
	x1=(A.minx < B.minx) ? B.minx : A.minx;
	// x2 = MIN ( A.maxx,B.maxx )
	x2=(A.maxx > B.maxx) ? B.maxx : A.maxx;

	// y1 = MAX ( A.miny,B.miny )
	y1=(A.miny < B.miny) ? B.miny : A.miny;
	// y2 = MIN ( A.maxy,B.maxy )
	y2=(A.maxy > B.maxy) ? B.maxy : A.maxy;

	return (x2-x1+1)*(y2-y1+1);
}


int BoxXDistance(Box A,Box B)
{
	int dx;


	if (A.maxx < B.minx)
		dx=B.minx-A.maxx;
	else if (A.minx > B.maxx)
		dx=A.minx-B.maxx;
	else
		/* overlap */
		dx=0;

	return(dx);

}


MovingObject::MovingObject()
{

	CompTable = new MovingComponentPtr[MAX_COMP_PER_OBJECT];
	Comp_cnt = 0;

}

MovingObject::~MovingObject()
{
	
	delete CompTable;
}

void MovingObject::ResetFlags()
{

	StartSplitFlag=0;
	EndSplitFlag=0;
	ClassifiedFlag=0;
	DisappearedFlag = 0;

}

void MovingObject::InitializeObject(unsigned int r,
									unsigned int c,
									unsigned int t,
									char * ID)
{

	ObjectState=APPEARING;
	ObjectClass=UNCLASSIFIED;
	TrackingID=0;
	FirstAppearance=t;
	rows=r;
	cols=c;
	Comp_cnt=0;

	MatchFlag=1;
	ConflictFlag=0;
	MatchCount=0;	

	ResetFlags();

	// initialize history
	AveHeight=0;
	AveWidth=0;
	AveNPixels=0;

	X=0;
	Y=0;

	TrajCount=0;
	
	BBOffsetX=0;
	BBOffsetY=0;

	strcpy(FeatureFileName,ID);
}


void MovingObject::InitializeChildObject(unsigned int r,
										unsigned int c,
										unsigned int t,
										MovingComponent * comp)
{

	ObjectState=APPEARING;
	ObjectType=CHILD_OBJECT;
	ObjectClass=UNCLASSIFIED;
	TrackingID=0;
	FirstAppearance=t+1;
	rows=r;
	cols=c;
	Comp_cnt=0;

	MatchFlag=1;
	ConflictFlag=0;
	MatchCount=0;
	
	ResetFlags();

	// initialize history
	AveHeight=0;
	AveWidth=0;
	AveNPixels=0;

	X=0;
	Y=0;

	TrajCount=0;

	BBOffsetX=0;
	BBOffsetY=0;


	BBox=comp->BBox;

	area=(BBox.maxx-BBox.minx+1)*(BBox.maxy-BBox.miny+1);


}

void MovingObject::InitializeOcclusionObject(unsigned int r,
										unsigned int c,
										unsigned int t)
{
	//ObjectState=OCCLUSIONOBJECT;
	ObjectState=TRACKED;
	ObjectType=MERGE_OBJECT;

	TrackingID=0;
	FirstAppearance=t;
	rows=r;
	cols=c;
	Comp_cnt=0;

	MatchFlag=0;
	ConflictFlag=0;
	MatchCount=0;

	ResetFlags();

	// initialize history

	TrajCount=0;


	AveHeight=0;
	AveWidth=0;
	AveNPixels=0;

	BBox.minx=0;
	BBox.maxx=cols;
	BBox.miny=0;
	BBox.maxy=rows;

	//InitializeObjectProfiles();
}


////  Remove
/*
void MovingObject::InitializeMultipleChildTracking()
{

	// allocate memory the first time this funciton is called..

	if ( ChildList == (MovingObject **) NULL) {
		ChildList= new MovingObject * [MAX_COMP_PER_OBJECT];
	}



}

*/

void MovingObject::AddComp(MovingComponent * comp)
{
	CompTable[Comp_cnt]=comp;
	Comp_cnt++;
	MatchFlag++;
	assert (Comp_cnt < MAX_COMP_PER_OBJECT);
}


void MovingObject::RemoveComp(unsigned int CompNo)
{
	CompTable[CompNo]=NULL;
	MatchFlag--;
}


void MovingObject::UpdateObjectRBBox()
{
	int bth;
	int i,sum;

	bth=nWPixels*0.01;
	bth=(bth < 1 ? 1 : bth);

	i=0;
	sum=0;
	while (sum < bth) 
		sum+=HorizontalProfile[i++];

	RBBox.minx=(i>0 ? i-1 : 0);

	i=cols-1;
	sum=0;
	while (sum < bth) 
		sum+=HorizontalProfile[i--];

	RBBox.maxx=(i<cols-1 ? i+1 : cols-1);

	assert(RBBox.minx <= RBBox.maxx);


	i=0;
	sum=0;
	while (sum < bth) 
		sum+=VerticalProfile[i++];

	RBBox.miny=(i>0 ? i-1 : 0);


	i=rows-1;
	sum=0;
	while (sum < bth) 
		sum+=VerticalProfile[i--];
		
	RBBox.maxy=(i<rows-1 ? i+1 : rows-1);

	assert(RBBox.miny <= RBBox.maxy);

}


void MovingObject::UpdateObjectAttributes(int * label_image)
{
	// update object attributes using Overlaping Components
	int minx,maxx,miny,maxy;
	unsigned int i;

	MatchCount++;

	minx=cols;
	miny=rows;
	maxx=0;
	maxy=0;

	Cx=0;
	Cy=0;
	
	SB=0;
	SR=0;
	SG=0;

	nWPixels=0;
	for (i=0;i<Comp_cnt;i++)
		if (CompTable[i] !=NULL) {
			nWPixels+=CompTable[i]->nWPixels;
			Cx+=CompTable[i]->cx*(int) CompTable[i]->nWPixels;
			Cy+=CompTable[i]->cy*(int) CompTable[i]->nWPixels;

			SB+=CompTable[i]->SB;
			SR+=CompTable[i]->SR;
			SG+=CompTable[i]->SG;


			if (CompTable[i]->BBox.minx <minx)
				minx=CompTable[i]->BBox.minx;
			if (CompTable[i]->BBox.miny <miny)
				miny=CompTable[i]->BBox.miny;
			if (CompTable[i]->BBox.maxx > maxx)
				maxx=CompTable[i]->BBox.maxx;
			if (CompTable[i]->BBox.maxy > maxy)
				maxy=CompTable[i]->BBox.maxy;
			
		}

	Cx=Cx/(int) nWPixels;
	Cy=Cy/(int) nWPixels;

	SB=SB / (int) nWPixels;
	SR=SR / (int) nWPixels;
	SG=SG / (int) nWPixels;

	int h;

	h=maxy-miny;

	/*
	BodySegments[0]=miny+0.85 * h;
	BodySegments[1]=miny+0.75 * h;
	BodySegments[2]=miny+0.6*h;
	BodySegments[3]=miny+0.4*h;
	BodySegments[4]=miny+0.15*h;
	*/
	
	BodySegments[0]=miny+ HEAD_BAND_MIN * h;
	BodySegments[1]=miny+ TOP_BAND_MAX * h;
	BodySegments[2]=miny+ TOP_BAND_MIN *h;
	BodySegments[3]=miny+ BOTTOM_BAND_MAX *h;
	BodySegments[4]=miny+ BOTTOM_BAND_MIN *h;


	if (ObjectState != OCCLUSIONOBJECT) {
		BBox.minx=minx;
		BBox.maxx=maxx;
		BBox.miny=miny;
		BBox.maxy=maxy;
	} else {
		// occlusion object
		// never expand..!
		if (minx > BBox.minx)
			BBox.minx=minx;
		if (maxx < BBox.maxx)
			BBox.maxx=maxx;
		if (miny > BBox.miny)
			BBox.miny=miny;
		if (maxy < BBox.maxy)
			BBox.maxy=maxy;
	}

	/*
	if (ObjectState >=TRACKED) {
		assert(nWPixels==PnWPixels);
		UpdateObjectRBBox();

	}
	*/

	//UpdateObjectProfiles(label_image);
	
	// update object speed
	/*
	Xspeed=(BBox.maxx+BBox.minx)/2-X;
	Yspeed=(BBox.maxy+BBox.miny)/2-Y;
	*/

	Xspeed=Cx-X;
	Yspeed=Cy-Y;

	// update object location (X,Y)
	/*
	X=(BBox.maxx+BBox.minx)/2;
	Y=(BBox.maxy+BBox.miny)/2;
	*/

	X=Cx;
	Y=Cy;

	//TrajX[TrajCount]=X;
	//TrajY[TrajCount]=Y;
	//TrajCount++;
	
	if (TimeStamp==FirstAppearance){
		Xst=X;
		Yst=Y;
		Xspeed=0;
		Yspeed=0;

		AveR=SR;
		AveG=SG;
		AveB=SB;

	}

	// use current BBox as a Prediction..
	// UpdateObjectPrediction may overwrite this prediction.

	PredictedBBox=BBox;
	/*
	PredictedBBox.minx=minx;
	PredictedBBox.maxx=maxx;
	PredictedBBox.miny=miny;
	PredictedBBox.maxy=maxy;
	*/
	area=(maxx-minx+1)*(maxy-miny+1);

	LastAppearance=TimeStamp;

	// update Average Height and Width
	if (ObjectState > APPEARING && 
		ObjectState != OCCLUSIONOBJECT) {

		int n=((int) LastAppearance- (int) FirstAppearance- (int) APPEARING_TH-1);

		n=(n<1 ? 1 : n);

	/*
		AveHeight=(AveHeight*(n-1) + (BBox.maxy-BBox.miny+1))/(n);
		AveWidth=( AveWidth*(n-1) + (BBox.maxx-BBox.minx+1))/(n);
		AveNPixels= (AveNPixels*(n-1) + nWPixels) / (n);

	*/

		double alpha,beta;

		alpha=0.7;

		beta= (1.0/n > alpha ?  1.0 / n : alpha);

		AveHeight=(AveHeight*(1-beta) + (BBox.maxy-BBox.miny+1) * beta);
		AveWidth=( AveWidth*(1-beta) + (BBox.maxx-BBox.minx+1) * beta );
		AveNPixels= (unsigned int) (AveNPixels*(1-beta) + nWPixels * beta );

		AveR= (AveR*(n-1)+SR) / (n);
		AveB= (AveB*(n-1)+SB) / (n);
		AveG= (AveG*(n-1)+SG) / (n);

		/*
		AveXspeed= (AveXspeed*n+ Xspeed)/(n+1);
		AveYspeed= (AveYspeed*n+ Yspeed)/(n+1);
		*/

		/*
		BBOffsetX=(BBOffsetX*n+((int) Cx-BBox.minx))/(n+1);
		BBOffsetY=(BBOffsetY*n+((int) Cy-BBox.miny))/(n+1);
		*/


		int n1=((int) LastAppearance- (int) FirstAppearance+1);
		AveXspeed=(double) (X-Xst)/ n1;
		AveYspeed=(double) (Y-Yst)/ n1;


	}

	if(LastAppearance-FirstAppearance > CLASSIFICATION_TIME &&
		ObjectClass==UNCLASSIFIED)
		UpdateObjectClassification();

}


void MovingObject::InitializeObjectProfiles()
{

	HorizontalProfile= new unsigned int[cols];
	VerticalProfile= new unsigned int[rows];
	ResetObjectProfiles();

}

void MovingObject::ResetObjectProfiles()
{
	unsigned int i;

	for(i=0;i<rows;i++)
		VerticalProfile[i]=0;

	for(i=0;i<cols;i++)
		HorizontalProfile[i]=0;


	PnWPixels=0;
}

void MovingObject::UpdateObjectProfiles(int * label_image)
{


}


int MovingObject::IsStationary()
{
	double speed;
	double disp;

	speed=AveXspeed*AveXspeed+AveYspeed*AveYspeed;
	disp=(X-Xst)*(X-Xst)+(Y-Yst)*(Y-Yst);

	if ((speed < MOVING_OBJECT_MIN_SPEED &&
		disp  < 10) ) {

		return (1);
	}
	else return (0);


}

int MovingObject::IsSmallerThanPerson()
{

	int H,W;

	H=BBox.maxy-BBox.miny;
	W=BBox.maxx-BBox.minx;


	if (H< 20 || W < 10) {
		return (1);
	}
	else return (0);


}

void MovingObject::UpdateObjectClassification()
{
	double speed,freq;
	double AR,HW;
	double disp;

	speed=AveXspeed*AveXspeed+AveYspeed*AveYspeed;
	freq= (double) MatchCount/(double) (LastAppearance-FirstAppearance+1);
	AR=AveHeight/AveWidth;
	HW=AveHeight*AveWidth;
	disp=(X-Xst)*(X-Xst)+(Y-Yst)*(Y-Yst);

	if //(speed < MOVING_OBJECT_MIN_SPEED
		//&& disp  < 25) 
		//||
		(freq  < MOVING_OBJECT_MIN_MATCH_FREQUENCY )
	{
		ObjectClass=STATIONARY;
		ClassifiedFlag=1;
	}
	else if( AR > PERSON_MIN_AR &&
			 AR < PERSON_MAX_AR && 
			 (BBox.maxy-BBox.miny) > 60)
	{

		ObjectClass=PERSON;
		ClassifiedFlag=1;
	}
	else if(AR > VEICHLE_MIN_AR &&
			AR < VEICHLE_MAX_AR &&
			HW > VEICHLE_MIN_HW)
		ObjectClass=VEICHLE;
	else
		ObjectClass=UNCLASSIFIED;


}

void MovingObject::InitializeObjectPrediction()
{
	float StateVector[2];
	
	EstHeight=(float) (BBox.maxy-BBox.miny+1);
	EstWidth=(float) (BBox.maxx-BBox.minx+1);
	/*
	EstHeight=(float) AveHeight;
	EstWidth=(float) AveWidth;
	*/

	BBOffsetX=(float) Cx-BBox.minx;
	BBOffsetY=(float) Cy-BBox.miny;

	StateVector[0]=(float) X;
	StateVector[1]=(float) 0;
	X_estimate.initialize(StateVector);

	StateVector[0]=(float) Y;
	StateVector[1]=(float) 0;
	Y_estimate.initialize(StateVector);

/*
	Xst=X;
	Yst=Y;
*/
	// open Feature File
	//FeatureFile=fopen(FeatureFileName,"w");

}


void MovingObject::InitializeObjectTemplate(unsigned char * image,
											unsigned int * FGimage)
{
	TemplateMissCount=0;
	unsigned int tX,tY;

	tX=(unsigned int) X;
	tY=(unsigned int) Y;

	ObjectTemplate = new TemporalTemplate(150,150,rows,cols);
	ObjectTemplate->InitializeTemplate(BBox,tX,tY,image,(unsigned char *) FGimage);
	
}

void MovingObject::InitializeObjectPredictionWithTemplate(unsigned char * image,
														  unsigned int * FGimage)
{
	float StateVector[2];
	
	TemplateMissCount=0;

	template_X=(BBox.maxx+BBox.minx)/2;
	template_Y=(BBox.maxy+BBox.miny)/2;

	ObjectTemplate = new TemporalTemplate(150,150,rows,cols);
	ObjectTemplate->InitializeTemplate(BBox,template_X,template_Y,image,(unsigned char *) FGimage);
	
	StateVector[0]=(float) (template_X);
	StateVector[1]=(float) 0;
	template_X_estimate.initialize(StateVector);

	StateVector[0]=(float) (template_Y);
	StateVector[1]=(float) 0;
	template_Y_estimate.initialize(StateVector);


}

void MovingObject::UpdateObjectPrediction()
{
	float measurement[2];


	// update Height and width estimate
	EstHeight=EstHeight*ALPHA+(1-ALPHA)*(BBox.maxy-BBox.miny+1);
	EstWidth=EstWidth*ALPHA+(1-ALPHA)*(BBox.maxx-BBox.minx+1);

	//update offset estimate
	BBOffsetX=BBOffsetX*ALPHA+(1-ALPHA)*(Cx-BBox.minx);
	BBOffsetY=BBOffsetY*ALPHA+(1-ALPHA)*(Cy-BBox.miny);
										
	// update center estimate
	measurement[0]=X;
	measurement[1]=Xspeed;
	X_estimate.update(measurement,1);
	

	measurement[0]=Y;
	measurement[1]=Yspeed;
	Y_estimate.update(measurement,1);
	
	// generate BBox prediction
	// using Box center as the estimate..
	PredictedBBox.minx=X_estimate.x[0]-EstWidth/2;
	PredictedBBox.maxx=X_estimate.x[0]+EstWidth/2;
	PredictedBBox.miny=Y_estimate.x[0]-EstHeight/2;
	PredictedBBox.maxy=Y_estimate.x[0]+EstHeight/2;

	// generate BBox prediction
	// using center of mass as the estimate..
	PredictedBBox.minx=X_estimate.x[0]-BBOffsetX;
	PredictedBBox.miny=Y_estimate.x[0]-BBOffsetY;
	PredictedBBox.maxx=PredictedBBox.minx+EstWidth;
	PredictedBBox.maxy=PredictedBBox.miny+EstHeight;

}

void MovingObject::UpdateObjectTemplate(unsigned char * image,
										unsigned int * FGimage,
										FILE * FP)
{
	unsigned int estX,estY;
	unsigned int newX,newY;
	unsigned int Wx,Wy;
	double v;

	// estimated best correlation location
	estX=(unsigned int) X;
	estY=(unsigned int) Y;

	Wx=7;
	Wy=7;

	if (ObjectTemplate->FindBestCorrelation(image,estX,estY,Wx,Wy,newX,newY,v)) {
		//correlation success
		TemplateMissCount=0;
		ObjectTemplate->UpdateTemplate(BBox,newX,newY,image,(unsigned char *) FGimage);
	}
	else{
		// correlation miss
		TemplateMissCount++;
		if (TemplateMissCount >=10){
			//fprintf(FP,"initialize template");
			InitializeObjectTemplate(image,FGimage);
		}
	}

	CorrelationResult=v;
	//fprintf(FP,"f: %d  c: %4.0f",TimeStamp,v);
	

}

void MovingObject::UpdateObjectPredictionWithTemplate(unsigned char * image,
													  unsigned int * FGimage,
													  FILE * FP)
{
	float measurement[2];
	unsigned int estX,estY;
	unsigned int newX,newY;
	unsigned int dx,dy,Wx,Wy;
	double v;

	
	estX=(unsigned int) template_X_estimate.x[0];
	estY=(unsigned int) template_Y_estimate.x[0];

	dx=(unsigned int) template_X_estimate.P[0][0];
	dy=(unsigned int) template_Y_estimate.P[0][0];
	
	dx=(dx>5) ? 5 : dx;
	dy=(dy>5) ? 5 : dy;

	Wx=2*dx+1;
	Wy=2*dy+1;


	
	if (ObjectTemplate->FindBestCorrelation(image,estX,estY,Wx,Wy,newX,newY,v)) {
		//correlation success
		TemplateMissCount=0;
		// update center estimate

		measurement[0]=(float) newX;
		measurement[1]=(float) newX- (float) template_X;
		template_X_estimate.update(measurement,1);

		measurement[0]=(float) newY;
		measurement[1]=(float) newY- (float) template_Y;
		template_Y_estimate.update(measurement,1);

		template_X=newX;
		template_Y=newY;

		ObjectTemplate->UpdateTemplate(BBox,newX,newY,image,(unsigned char *) FGimage);
	}
	else {
		// correlation miss
		//TemplateMissCount++;
		if (TemplateMissCount >=20){
			//fprintf(FP,"initialize template");
			InitializeObjectPredictionWithTemplate(image,FGimage);
		}
		else {
			measurement[0]=(float) newX;
			measurement[1]=(float) newX- (float) template_X;
			template_X_estimate.update(measurement,0);

			measurement[0]=(float) newY;
			measurement[1]=(float) newY- (float) template_Y;
			template_Y_estimate.update(measurement,0);

			template_X=newX;
			template_Y=newY;
		}
	}
	/*
	fprintf(FP,"f: %d  c: %4.0f \n",TimeStamp,v,
		X_estimate.P[0][0],X_estimate.P[1][1],
		Y_estimate.P[0][0],Y_estimate.P[1][1]);
	*/

	estX=(unsigned int) template_X_estimate.x[0];
	estY=(unsigned int) template_Y_estimate.x[0];

	TemplateBBox=ObjectTemplate->GetBBox(estX,estY);
	//PredictedBBox=TemplateBBox;
}

void MovingObject::InitializeObjectSilhouette()
{
	pfc = new pfContour(rows,cols);

}

void MovingObject::UpdateObjectSilhouette(int * label_image)
{

	MPoint p;

	p.x=Cx;
	p.y=Cy;


	pfc->generate(label_image,&p,ObjectLabel);
}


// to be removed
void MovingObject::InitializeObjectColorModels()
{

	TopModel = new ColorModel(3,rows,cols);
	BottomModel = new ColorModel(3,rows,cols);
	HeadModel = new ColorModel(3,rows,cols);

	
}

void MovingObject::InitializeColorBlobs()
{

	HeadBlob = new ColorBlob(rows,cols,3);
	TopBlob = new ColorBlob(rows,cols,3);
	BottomBlob= new ColorBlob(rows,cols,3);

	Nseg=0;

	// to be removed
	TopModel = TopBlob->colormodel;
	BottomModel = BottomBlob->colormodel;
	HeadModel = HeadBlob->colormodel;

	
}

void MovingObject::UpdateColorBlobPrediction()
{

	HeadBlob->PBBox.minx=BBox.minx;
	HeadBlob->PBBox.maxx=BBox.maxx;

	HeadBlob->PBBox.miny=BodySegments[1];
	HeadBlob->PBBox.maxy=BBox.maxy;

	TopBlob->PBBox.minx=BBox.minx;
	TopBlob->PBBox.maxx=BBox.maxx;

	TopBlob->PBBox.miny=BodySegments[4];
	TopBlob->PBBox.maxy=BodySegments[0];


	BottomBlob->PBBox.minx=BBox.minx;
	BottomBlob->PBBox.maxx=BBox.maxx;

	BottomBlob->PBBox.miny=BBox.miny;
	BottomBlob->PBBox.maxy=BodySegments[1];
   

	
}

void MovingObject::CaptureObjectColorModels(unsigned char * image,
																						double * srgimage,
																						int * label_image)

{
	unsigned int x,y;

	BottomModel->Clear();
	TopModel->Clear();

	HeadModel->Clear();

	for(y = BodySegments[4]; y < BodySegments[3] ; y++)
		for(x = BBox.minx;x<BBox.maxx;x++)
			if (label_image[y*cols+x]==ObjectLabel)
				BottomModel->AddSample(image+3*(y*cols+x),srgimage+3*(y*cols+x));

	for(y = BodySegments[2]; y < BodySegments[1] ; y++)
		for(x = BBox.minx;x<BBox.maxx;x++)
			if (label_image[y*cols+x]==ObjectLabel)
				TopModel->AddSample(image+3*(y*cols+x),srgimage+3*(y*cols+x));

	for(y = BodySegments[0]; y < BBox.maxy ; y++)
		for(x = BBox.minx;x<BBox.maxx;x++)
			if (label_image[y*cols+x]==ObjectLabel)
				HeadModel->AddSample(image+3*(y*cols+x),srgimage+3*(y*cols+x));


 
	/*
	HeadModel->InitializeFastGaussComputation();
	TopModel->InitializeFastGaussComputation();
	BottomModel->InitializeFastGaussComputation();
	*/

	unsigned int comp=   FAST_GAUSS_3D; //DIRECT_EVALUATION; //  KERNEL_LUT; 


	HeadModel->InitializeComputation(comp);
	TopModel->InitializeComputation(comp);
	BottomModel->InitializeComputation(comp);
	
	/*
	// write sample to files
	char fname[30];


	sprintf(fname,"Obj%d_Head.txt",TrackingID);
	HeadModel->WriteSamples(fname);

	sprintf(fname,"Obj%d_Top.txt",TrackingID);
	TopModel->WriteSamples(fname);

	sprintf(fname,"Obj%d_Bottom.txt",TrackingID);
	BottomModel->WriteSamples(fname);
	*/
	
}

void MovingObject::DetectSeparators()
{

	int head_cnt_acc[240];
	int top_cnt_acc[240];
	int bottom_cnt_acc[240];

	unsigned int * head_cnt,* top_cnt, * bottom_cnt;
	unsigned int y;

	for (y=0;y <rows; y++) {
		head_cnt_acc[y]=0;
		top_cnt_acc[y]=0;
		bottom_cnt_acc[y]=0;
	}

	head_cnt=HeadBlob->VerticalProfile;
	top_cnt=TopBlob->VerticalProfile;
	bottom_cnt=BottomBlob->VerticalProfile;

	for (y=1;y<rows;y++) {
		head_cnt_acc[y]=head_cnt_acc[y-1]+head_cnt[y];
		top_cnt_acc[y]=top_cnt_acc[y-1]+top_cnt[y];
		bottom_cnt_acc[y]=bottom_cnt_acc[y-1]+bottom_cnt[y];
	}


	{
		int minmissclass_tb=bottom_cnt_acc[rows-1];
		int minmissclass_ht=top_cnt_acc[rows-1];
		int missclass_tb;
		int missclass_ht;

		//for (y=BBox.miny;y<BBox.maxy;y++) {

		for (y=0;y<rows;y++) {

			missclass_tb=top_cnt_acc[y]+(bottom_cnt_acc[rows-1]-bottom_cnt_acc[y]);
			missclass_ht=head_cnt_acc[y]+(top_cnt_acc[rows-1]-top_cnt_acc[y]);

			if (missclass_tb < minmissclass_tb) {
				minmissclass_tb=missclass_tb;
				top_bottom_sep=y;					
			}

			if (missclass_ht < minmissclass_ht) {
				minmissclass_ht=missclass_ht;
				head_top_sep=y;					
			}

			
		}
	}

}


void MovingObject::UpdateBlobAttributes(unsigned char AdaptationFlag)
{

	// separators are detected.

	// determine scale
	scale=(float) (head_top_sep-top_bottom_sep) / (float) initial_top_bottom_dist;

	HeadBlob->UpdateAttributes(scale,AdaptationFlag);
	TopBlob->UpdateAttributes(scale,AdaptationFlag);
	BottomBlob->UpdateAttributes(scale,AdaptationFlag);


	// object observed bounding box. maximal blob boxes
	BBox.maxy=HeadBlob->BBox.maxy;
	BBox.minx=BottomBlob->BBox.minx;
	
	BBox.minx=HeadBlob->BBox.minx;

	if (TopBlob->BBox.minx < BBox.minx)
		BBox.minx=TopBlob->BBox.minx;
	if (BottomBlob->BBox.minx < BBox.minx)
		BBox.minx=BottomBlob->BBox.minx;


	
	BBox.maxx=HeadBlob->BBox.maxx;

	if (TopBlob->BBox.maxx > BBox.maxx)
		BBox.maxx=TopBlob->BBox.maxx;
	if (BottomBlob->BBox.maxx > BBox.maxx)
		BBox.maxx=BottomBlob->BBox.maxx;




	if (AdaptationFlag) {
		Nseg++;
	
		head_top_sep_disp=((Nseg-1)*head_top_sep_disp+(BBox.maxy-head_top_sep)/scale)/(Nseg);

		top_bottom_sep_disp=((Nseg-1)*top_bottom_sep_disp+(BBox.maxy-top_bottom_sep)/scale)/(Nseg);

		HeadBlob->AveYdisp=0.5* head_top_sep_disp ;
		TopBlob->AveYdisp= 0.5 * (top_bottom_sep_disp + head_top_sep_disp);
		BottomBlob->AveYdisp = (BBox.maxy- BottomBlob->medianY) / scale;

	}



}

void MovingObject::ColorSegment(unsigned char * image, 
								int * object_label_image,
								unsigned char * part_label_image,
								KernelTableStruct * KernelTablePtr)
{
	unsigned char * ptr;
	double p_top,p_bottom,p_head;
	unsigned int x,y;
	
	int head_cnt[240];
	int top_cnt[240];
	int bottom_cnt[240];

	int head_cnt_acc[240];
	int top_cnt_acc[240];
	int bottom_cnt_acc[240];


	//double OutlierTh=0.0001;
	double OutlierTh=0.000001;

	// labels
	// 1	top
	// 2	bottom
	// 3	head
	// 4	bottom outlier
	// 5	outlier
	// 6  missclassification (relabeled)


	for (y=0;y <rows; y++) {
		head_cnt[y]=0;
		top_cnt[y]=0;
		bottom_cnt[y]=0;
		head_cnt_acc[y]=0;
		top_cnt_acc[y]=0;
		bottom_cnt_acc[y]=0;
	}


	HeadBlob->ResetProfiles();
	TopBlob->ResetProfiles();
	BottomBlob->ResetProfiles();


	for (y=BBox.miny; y < BodySegments[4]; y++)
		for (x=BBox.minx ; x<BBox.maxx ; x++)
			if (object_label_image[y*cols+x]==ObjectLabel) {
				ptr=image+3*(y*cols+x);
				p_bottom=BottomModel->EstimateProbability(ptr,KernelTablePtr);

				if (p_bottom > OutlierTh) {
						part_label_image[y*cols+x]=2;
						bottom_cnt[y]++;
						//BottomBlob->AddPixel(x,y);
				}
					else
						part_label_image[y*cols+x]=4; //5;

			}


	for (y=BodySegments[4]; y < BodySegments[1]; y++)
		for (x=BBox.minx ; x<BBox.maxx ; x++)
			if (object_label_image[y*cols+x]==ObjectLabel) {
				ptr=image+3*(y*cols+x);
				p_top=TopModel->EstimateProbability(ptr,KernelTablePtr);
				p_bottom=BottomModel->EstimateProbability(ptr,KernelTablePtr);

				if (p_top > p_bottom) {
					if (p_top > OutlierTh) {
						part_label_image[y*cols+x]=1;
						top_cnt[y]++;
						//TopBlob->AddPixel(x,y);

					}
					else 
						part_label_image[y*cols+x]=5;
				}
					
				else if (p_bottom > OutlierTh) {
						part_label_image[y*cols+x]=2;
						bottom_cnt[y]++;
						//BottomBlob->AddPixel(x,y);

				}
				else
						part_label_image[y*cols+x]=5;
			}

	for (y= BodySegments[1]; y < BBox.maxy; y++)
		for (x=BBox.minx ; x<BBox.maxx ; x++)
			if (object_label_image[y*cols+x]==ObjectLabel) {
				ptr=image+3*(y*cols+x);
				if (ptr[0]==255 && ptr[1]==255 && ptr[2]==255) {
					p_top=0.0;
					p_head=0.0;
				}
				else {
					p_top=TopModel->EstimateProbability(ptr,KernelTablePtr);
					p_head=HeadModel->EstimateProbability(ptr,KernelTablePtr);
				}
				

				if (p_top >= p_head) {
					if (p_top > OutlierTh) {
						part_label_image[y*cols+x]=1;
						top_cnt[y]++;
						//TopBlob->AddPixel(x,y);

					}
					else 
						part_label_image[y*cols+x]=5;
				}		
				else if (p_head > OutlierTh) {
						part_label_image[y*cols+x]=3;
						head_cnt[y]++;
						//HeadBlob->AddPixel(x,y);

				}
				else
						part_label_image[y*cols+x]=5;
			}





   // detect separators

	for (y=1;y<rows;y++) {
		head_cnt_acc[y]=head_cnt_acc[y-1]+head_cnt[y];
		top_cnt_acc[y]=top_cnt_acc[y-1]+top_cnt[y];
		bottom_cnt_acc[y]=bottom_cnt_acc[y-1]+bottom_cnt[y];
	}


	{
		int minmissclass_tb=bottom_cnt_acc[rows-1];
		int minmissclass_ht=top_cnt_acc[rows-1];
		int missclass_tb;
		int missclass_ht;

		for (y=BBox.miny;y<BBox.maxy;y++) {

			missclass_tb=top_cnt_acc[y]+(bottom_cnt_acc[rows-1]-bottom_cnt_acc[y]);
			missclass_ht=head_cnt_acc[y]+(top_cnt_acc[rows-1]-top_cnt_acc[y]);

			if (missclass_tb < minmissclass_tb) {
				minmissclass_tb=missclass_tb;
				top_bottom_sep=y;					
			}

			if (missclass_ht < minmissclass_ht) {
				minmissclass_ht=missclass_ht;
				head_top_sep=y;					
			}

			
		}
	}

	// relabel

	// bottom
	for (y=BBox.miny; y <= top_bottom_sep; y++)
		for (x=BBox.minx ; x<BBox.maxx ; x++)
			if (object_label_image[y*cols+x]==ObjectLabel) {
				if (part_label_image[y*cols+x]==BOTTOM_BLOB) {
						BottomBlob->AddPixel(x,y);				
				}
			}

	// top
	for (y=top_bottom_sep ; y <= head_top_sep; y++)
		for (x=BBox.minx ; x<BBox.maxx ; x++)
			if (object_label_image[y*cols+x]==ObjectLabel) {
				if (part_label_image[y*cols+x]==TOP_BLOB) {
						TopBlob->AddPixel(x,y);				
				}
			}


	// head
	for (y=head_top_sep ; y <BBox.maxy; y++)
		for (x=BBox.minx ; x<BBox.maxx ; x++)
			if (object_label_image[y*cols+x]==ObjectLabel) {
				if (part_label_image[y*cols+x]==HEAD_BLOB) {
						HeadBlob->AddPixel(x,y);				
				}
			}


	if (Nseg==0) {
		initial_top_bottom_dist=head_top_sep-top_bottom_sep;
	}

	scale=(float) (head_top_sep-top_bottom_sep) / (float) initial_top_bottom_dist;

	Nseg++;
	
	head_top_sep_disp=((Nseg-1)*head_top_sep_disp+(BBox.maxy-head_top_sep)/scale)/(Nseg);

	top_bottom_sep_disp=((Nseg-1)*top_bottom_sep_disp+(BBox.maxy-top_bottom_sep)/scale)/(Nseg);


	HeadBlob->UpdateAttributes(scale,1);
	TopBlob->UpdateAttributes(scale,1);
	BottomBlob->UpdateAttributes(scale,1);

	// update prob profiles


	HeadModel->UpdateProfile(head_cnt,BBox.miny,BBox.maxy,top_bottom_sep);
	TopModel->UpdateProfile(top_cnt,BBox.miny,BBox.maxy,top_bottom_sep);
	BottomModel->UpdateProfile(bottom_cnt,BBox.miny,BBox.maxy,top_bottom_sep);




	HeadBlob->AveYdisp=0.5* head_top_sep_disp ;
	TopBlob->AveYdisp= 0.5 * (top_bottom_sep_disp + head_top_sep_disp);
	BottomBlob->AveYdisp = (BBox.maxy- BottomBlob->medianY) / scale;

	 
}



void MovingObject::KillObject()
{
	/* Used for Objects being Tracked*/

	//delete HorizontalProfile;
	//delete VerticalProfile;

	ObjectState=DEAD;

	//fclose(FeatureFile);
}

void MovingObject::UpdateObjectState(unsigned char * image,
									 unsigned char * FGimage,
									 int * label_image,
									 FILE * FP)
{
	/* Change object state */
	switch (ObjectState){
	case APPEARING :	// appearing mode
		if (!MatchFlag)
			// KillObject();
			ObjectState=DEAD;
		else if(!ConflictFlag){

			if (TimeStamp-FirstAppearance >  APPEARING_TH)
			{
				ObjectState=INITIALIZING;
				//InitializeObjectProfiles();
			}
			UpdateObjectAttributes(label_image);
		}
		else // conflict matching
		{

			/*
			if (TimeStamp-FirstAppearance > APPEARING_TH)
				ObjectState=INITIALIZING;
				//InitializeObjectProfiles();

			UpdateObjectAttributes(label_image);

			*/
			ObjectState=DEAD;


		}
		break;

	case INITIALIZING :	// appearing mode
		if (!MatchFlag) {
			// no matching
			if (TimeStamp-LastAppearance > DISAPPEARING_TH)
				KillObject();
				//ObjectState=DEAD;
		}
		else if(!ConflictFlag){
			UpdateObjectAttributes(label_image);
			if (TimeStamp-FirstAppearance - APPEARING_TH > 
				(ObjectType==CHILD_OBJECT ? CHILD_INITIALIZING_TH : INITIALIZING_TH) )
			{
				InitializeObjectPrediction();
				//InitializeObjectSilhouette();
				//InitializeObjectColorModels(image,label_image);
				//InitializeObjectTemplate(image,FGimage);
				//InitializeObjectPredictionWithTemplate(image,FGimage);
				ObjectState=TRACKED;
			}
		}
		else // conflict matching
		{
			/*
			UpdateObjectAttributes(label_image);
			if (TimeStamp-FirstAppearance - APPEARING_TH > INITIALIZING_TH)
			{
				ObjectState=TRACKED;	
				InitializeObjectPrediction();
				//InitializeObjectTemplate(image,FGimage);
				//InitializeObjectPredictionWithTemplate(image,FGimage);
			}
			*/

			ObjectState=DEAD;

		}

		//ResetObjectProfiles();
		break;


	case TRACKED : // self mode
		if (!MatchFlag)	{
			// no matching
			//fprintf(FP,"no matching \n");
			if (TimeStamp-LastAppearance > DISAPPEARING_TH) 
			{
				KillObject();
				DisappearedFlag=1;
			}

		}
		else if (ConflictFlag)
		{
			// freeze object...
			ObjectState=DEAD;  //UNDER_MERGE;
				//ObjectState=UNDER_MERGE;
		}
		else if (MatchFlag == 1) 
		{	// 1-1 matching
			UpdateObjectAttributes(label_image);
			//UpdateObjectSilhouette(label_image);

			UpdateObjectPrediction();
			//UpdateObjectTemplate(image,FGimage,FP);
			//fprintf(FP," class :%d \n",ObjectClass);
			WriteFeatureFile();
			//UpdateObjectPredictionWithTemplate(image,FGimage,FP);

		}	
		else
		{	// 1-M matching

			UpdateObjectAttributes(label_image);
			UpdateObjectPrediction();

			ObjectState=TRACKED_MULTIPLE;
			StartSplitFlag=1;
			//InitializeMultipleChildTracking();



		}


		//ResetObjectProfiles();
		break;
	
	case TRACKED_MULTIPLE : // Multiple blob Tracking
		// Update Matching status using child info
		
		
		if (!MatchFlag)	{
			// no matching
			//fprintf(FP,"no matching \n");
			if (TimeStamp-LastAppearance > DISAPPEARING_TH) 
			{
				KillObject();
				DisappearedFlag=1;

			}
		}
		else if (MatchFlag == 1) 
		{	// 1-1 matching
			UpdateObjectAttributes(label_image);
			UpdateObjectPrediction();

			ObjectState=TRACKED;
			EndSplitFlag=1;

			//EndSplit();
			//EndMultipleBlobTracking();

			//UpdateObjectTemplate(image,FGimage,FP);
			//fprintf(FP," class :%d \n",ObjectClass);
			WriteFeatureFile();
			//UpdateObjectPredictionWithTemplate(image,FGimage,FP);

		}	
		else
		{	// 1-M matching

			UpdateObjectAttributes(label_image);
			UpdateObjectPrediction();

			//CheckSplit();
			//MultipleBlobTracking();


		}

		break;

	case OCCLUSIONOBJECT :	// OCCLUSION OBJECT
		if (!MatchFlag)	{
			// no matching
			//fprintf(FP,"no matching \n");
			if (TimeStamp-LastAppearance > 3)
				ObjectState=DEAD;
		}
		else{
			UpdateObjectAttributes(label_image);
		}
	
		//ResetObjectProfiles();

		break;

	case UNDER_MERGE :	// Object under Merge situation
		break;

	}

	
}

void MovingObject::WriteTrajectory()
{
	FILE * TrajFP;
	TrajFP=fopen("Traj","w");

	for (unsigned int i=0;i<TrajCount;i++)
		fprintf(TrajFP,"%d %d ",TrajX[i],TrajY[i]);

	fclose(TrajFP);

}

void MovingObject::WriteFeatureFile()
{
	// Write Feature File
	/*
	fprintf(FeatureFile,"%1.3f %3d %3d %1d %1d %4d %1.2f %3.1f %3.1f %1.2f %3.1f %3.1f %1.2f %2.3f %2.3f %f\n",
		(double) MatchCount/ (double) (LastAppearance-FirstAppearance+1),
		abs(X-Xst),
		abs(Y-Yst),
		abs(Xspeed),
		abs(Yspeed),
		nWPixels,
		(double) (BBox.maxy-BBox.miny+1)/(double) (BBox.maxx-BBox.minx+1),
		AveHeight,
		AveWidth,
		AveHeight/AveWidth,
		EstHeight,
		EstWidth,
		(double)EstHeight/ (double)EstWidth,
		X_estimate.P[0][0],
		Y_estimate.P[0][0],
		CorrelationResult
		);
	*/
	/*
	fprintf(FeatureFile,"%3d %1d %1d %3d %3d %1.2f %1.2f %4d %1.2f %d %d class : %d\n",
		TimeStamp-FirstAppearance,
		Xspeed,
		Yspeed,
		X-Xst,
		Y-Yst,
		(double) (X-Xst)*(X-Xst)+(Y-Yst)*(Y-Yst),
		AveXspeed*AveXspeed+AveYspeed*AveYspeed,
		nWPixels,
		AveHeight/AveWidth,
		(BBox.maxy-BBox.miny+1),
		(BBox.maxx-BBox.minx+1),
		ObjectClass
		);
	fflush(FeatureFile);
	*/
}

void MovingObject::UpdateObjectHistory()
{

	// delete Moving Components
	/*
	for (i=0;i<Comp_cnt;i++)
		delete CompTable[i];
	*/
	//initialize variables for next frame
	Comp_cnt=0;
	MatchFlag=0;
	ConflictFlag=0;


}

void MovingObject::DrawBB(unsigned char * image,unsigned char color)
{
#define IMAGE(i,j)  *(image+(i)*cols+(j))

  int x,y;
  for(x=BBox.minx;x<BBox.maxx;x++){
    IMAGE(BBox.miny,x)=color;
    IMAGE(BBox.maxy,x)=color;
  }

  for(y=BBox.miny;y<BBox.maxy;y++){
    IMAGE(y,BBox.minx)=color;
    IMAGE(y,BBox.maxx)=color;
  }
#undef IMAGE
}

void MovingObject::DrawPredictedBB(unsigned char * image,unsigned char color)
{
#define IMAGE(i,j)  *(image+(i)*cols+(j))

  int x,y;
  if (PredictedBBox.miny > 0)
	for(x=PredictedBBox.minx;x<PredictedBBox.maxx;x++)
	if (x>0 && x<cols){
		IMAGE(PredictedBBox.miny,x)=color;
	}


  if (PredictedBBox.maxy < rows)
	for(x=PredictedBBox.minx;x<PredictedBBox.maxx;x++)
	if (x>0 && x<cols){
		IMAGE(PredictedBBox.maxy,x)=color;
	}

	
  if (PredictedBBox.minx > 0)
	for(y=PredictedBBox.miny;y<PredictedBBox.maxy;y++)
	if (y>0 && y<rows){
		IMAGE(y,PredictedBBox.minx)=color;
	}

  if (PredictedBBox.maxx < cols)
	for(y=PredictedBBox.miny;y<PredictedBBox.maxy;y++)
	if (y>0 && y<rows){
		IMAGE(y,PredictedBBox.maxx)=color;
	}

#undef IMAGE
}

void MovingObject::DrawTemplateBB(unsigned char * image,unsigned char color)
{
#define IMAGE(i,j)  *(image+(i)*cols+(j))

  int x,y;
  if (TemplateBBox.miny > 0)
	for(x=TemplateBBox.minx;x<TemplateBBox.maxx;x++)
	if (x>0 && x<cols){
		IMAGE(TemplateBBox.miny,x)=color;
	}


  if (TemplateBBox.maxy < rows)
	for(x=TemplateBBox.minx;x<TemplateBBox.maxx;x++)
	if (x>0 && x<cols){
		IMAGE(TemplateBBox.maxy,x)=color;
	}

	
  if (TemplateBBox.minx > 0)
	for(y=TemplateBBox.miny;y<TemplateBBox.maxy;y++)
	if (y>0 && y<rows){
		IMAGE(y,TemplateBBox.minx)=color;
	}

  if (TemplateBBox.maxx < cols)
	for(y=TemplateBBox.miny;y<TemplateBBox.maxy;y++)
	if (y>0 && y<rows){
		IMAGE(y,TemplateBBox.maxx)=color;
	}

#undef IMAGE
}
