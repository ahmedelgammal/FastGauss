
#include <windows.h>

//#define NFMAX 0
//#define NLMAX 1000

#define MAXTERMS 30
#define MAXFFEXP 100
#define MAXDIM 5

class FastGauss {

private:
	double * sources;
	double * str;
	double * segma;
	double * targets;

	double * fixedsegma;	// N-Dim vector

	double * TargetField;

	unsigned int Nsources;	// number of sources
	unsigned int Ntargets;	// number of targets
	unsigned int Dim;				// 
	unsigned int Nterms;		// number of terms used in the expansions
	unsigned int Nterms2pDim;	// Nterms to the power (Dim)
	
	int * Nboxes_eachside;	// N-Dim number of source boxes along each dim.
	// int Nsource_boxes_eachside;  // number of source boxes along each dim.
	// int Ntarget_boxes_eachside;

	unsigned int NBoxes;							// total number of boxes;
	double * BoxCenter;		// N-Dim center of a source box.
	int * Box2Subscript; // Box2Subscript[i*Dim] : N-Dim subscript of box i

	// unsigned int NSrcBoxes;	// total number of source boxes
	// unsigned int NTrgBoxes;	// total number of target boxes

	
	int * SrcsInBoxCount;		// SrcsInBoxCount[i] : number of sources in box i
	int * SrcsInBox;				// SrcsInBox[i] : index of the first source in box i
	int * SrcNext;					// SrcNext[i] : index of next source in a linked list of sources in the same box

	// target boxes
	int * TrgsInBoxCount;		// TrgsInBoxCount[i] : number of targets in box i
	int * TrgsInBox;				// TrgsInBox[i] : index of the first target in box i
	int * TrgNext;					// TrgNext[i] : index of next target in a linked list of targets in the same box
	int * UsedTrgsBox;			// list of used target boxes; UsedTrgsBox[i]: index of ith used the target box
	int UsedTrgsBoxCount;		// Number of target boxes that actually contatin targets.
	// double * TrgBoxCenter;		// N-Dim center of a target box.

	int NbrSize;
	int Ndisp;
	int * NbrDisp;
	int * NbrRelativeIndex;
		
	int * NbrList;


	// expansions work space
	double * FFexp;		// used for Far Field Exp. size (Nterms^Dim * MAXFFEXP) 
	double ** SrcBoxFFexp;	// pointer to Far Field exp space for each source box.
	int	Nffexp;			// number of used ffexp

	double * LocalExp;	// used to accumalate  Local expansion for each target box. 
										// size: Nboxes*Nterms^Dim
	
	double ** TrgBoxLocalExp; // pointer to local exp space for each target box.
	int NLocalExp;	// number of used Local Exp
	unsigned char * TrgBoxLocalExpFlag; // flag for each target box to indicate whether the local expansion has 
																		 //	been previosly evaluated. This is useful for reruns with same sources


	int * Nboxes_powers ; // powers of Nboxes_eachside;

	// cut points
	int NsMAX;	// = NFMAX	cut point for number of sources in a box
	int NtMAX;	// = NLMAX	cut point for number of target in a box

	int * tmp;
	double * Dtmp;

	int PreAllocateTargetFlag;

	
public:
	// statistical counters
	double stat_1;
	double stat_2;
	double stat_3;
	double stat_4;

	double timer_1;
	double timer_2;
	double timer_3;
	double timer_4;

public:

	FastGauss(unsigned int D);
	~FastGauss();
	int FixedScaleEvaluation(double * source_locs,
														double * source_strength,
														double * source_segma,
														double * target_locs,
														double * & FFresults,
														unsigned int NSources,
														unsigned int NTargets,
														double tolerance);

	int BuildSourceStructure(unsigned int NSources,
														double * source_locs,
														double * source_strength,
														double * source_segma,
														double tolerance);

	int EvaluateTargets(unsigned int NTargets,
											double * target_locs,
											double * & FFresults);

	int PreAllocateTargets(int MaxTargets);

private:
	int Allocate();
	int AllocateTargets();
	void GenerateBoxes();
	void AssignSources();
	void AssignTargets();
	void GenerateDisplacements();
	void GenerateNbrs(int boxidx,int * NbrsList,unsigned int & NbrsCount);
	int  DetermineNterms(double tol,double r);
	void DirectEvaluation(int srcbox,int trgbox);
	void BoxTaylorExpansion2D(int srcbox,int trgbox);
	void BoxTaylorExpansion3D(int srcbox,int trgbox);
	void FarFieldExpansion2D(int srcbox, double * ffexp);
	void FarFieldExpansion3D(int srcbox, double * ffexp);
	double EvalFarFieldExp_atTarget2D(double * trg,double * cent,double * ffexp);
	double EvalFarFieldExp_atTarget3D(double * trg,double * cent,double * ffexp);
	void EvalFarFieldExp(int trgbox,double * cent,double * ffexp);	
	void EvalLocalExpansions(); 
	void EvalBoxLocalExpansion(int trgbox,double * locexp);
	int Evaluateboxes2D();
	int Evaluateboxes3D();
	int SourceExpansion2D();
	int SourceExpansion3D();
	int AllocateLocalExpansion();
	int EvaluateTargetboxes2D();
	int EvaluateTargetboxes3D();





};


/*--------------------macros-----------------------------*/


// TIMED_FUNCTION macro call the function "fn" and measure
// the time required to execute this function and add that
// to "x"

#define FG_TIMED_FUNCTION(fn,x) \
{			\
	LARGE_INTEGER __t1,__t2,__f;  \
	QueryPerformanceFrequency(&__f);\
	QueryPerformanceCounter(&__t1);\
	(fn);\
	QueryPerformanceCounter(&__t2);	\
	(x)+=(float) (1000.0*(float)(__t2.QuadPart - __t1.QuadPart) / __f.QuadPart); \
} 
