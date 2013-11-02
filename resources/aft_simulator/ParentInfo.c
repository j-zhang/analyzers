/*
* ParentInfo.c is a C program to calculate a vector of the parent times 
* and locations
*
*[AftParentTimes, xAftParentLocs, yAftParentLocs, zAftParentLocs, ...  
*     PointSourceAftParents] ...
*        = ParentInfo(TotalNumOfPointSourceAfts, ... 
*          NumOfPointSourceParents,StartAft,EndAft,ParentIDs, ... 
*          ParentTimes,xParentLocs,yParentLocs,zParentLocs);
*
*
*
* Matlab Function:
* [AddOn, lx, ly, lz, History] = mainshockInfo (AxxEnd, LengthX, Sp, ...
* Ep, IDold, tdist, xl, yl, zl)
*
*/

# include "mex.h"
#include "matrix.h"

/* C subroutine for calculating mainshockInfo */
void mainshockInfo (double AftParentTimes[], double xAftParentLocs[], 
double yAftParentLocs[], double zAftParentLocs[], 
double PointSourceAftParents[], int TotalNumOfPointSourceAfts, 
int NumOfPointSourceParents, double StartAft[], double EndAft[], 
double ParentIDs[], double ParentTimes[], double xParentLocs[], 
double yParentLocs[], double zParentLocs[])

{
	int j,k;
	for (j=0; j<=NumOfPointSourceParents-1; j++) {
		for (k=StartAft[j]-1; k<=EndAft[j]-1; k++) {
			AftParentTimes[k]=ParentTimes[j];
			xAftParentLocs[k]=xParentLocs[j];
			yAftParentLocs[k]=yParentLocs[j];
			zAftParentLocs[k]=zParentLocs[j];
			PointSourceAftParents[k]=ParentIDs[j];;
		}
	}
}


/* Interface to Matlab */

void mexFunction (int nlhs, mxArray *plhs[],  int nrhs, const mxArray *prhs[])
{

/* Declare vector pointers */
double *AftParentTimes, *xAftParentLocs, *yAftParentLocs, *zAftParentLocs;
double *PointSourceAftParents;
int TotalNumOfPointSourceAfts, NumOfPointSourceParents;
double *ParentIDs;
double *StartAft, *EndAft; 
double *ParentTimes, *xParentLocs, *yParentLocs, *zParentLocs;

/* Get scalars */
TotalNumOfPointSourceAfts = mxGetScalar (prhs[0]);
NumOfPointSourceParents = mxGetScalar (prhs[1]);

/* Create matrices for return variables */
plhs[0] = mxCreateDoubleMatrix (TotalNumOfPointSourceAfts, 1, mxREAL);
plhs[1] = mxCreateDoubleMatrix (TotalNumOfPointSourceAfts, 1, mxREAL);
plhs[2] = mxCreateDoubleMatrix (TotalNumOfPointSourceAfts, 1, mxREAL);
plhs[3] = mxCreateDoubleMatrix (TotalNumOfPointSourceAfts, 1, mxREAL);
plhs[4] = mxCreateDoubleMatrix (TotalNumOfPointSourceAfts, 1, mxREAL);


/* Assign pointers to the vectors */
AftParentTimes=mxGetPr(plhs[0]);
xAftParentLocs=mxGetPr(plhs[1]);
yAftParentLocs=mxGetPr(plhs[2]);
zAftParentLocs=mxGetPr(plhs[3]);
PointSourceAftParents=mxGetPr(plhs[4]);
StartAft=mxGetPr(prhs[2]);
EndAft=mxGetPr(prhs[3]);
ParentIDs=mxGetPr(prhs[4]);
ParentTimes=mxGetPr(prhs[5]);
xParentLocs=mxGetPr(prhs[6]);
yParentLocs=mxGetPr(prhs[7]);
zParentLocs=mxGetPr(prhs[8]);

/* Call the subroutine */
mainshockInfo(AftParentTimes, xAftParentLocs, yAftParentLocs, 
zAftParentLocs, PointSourceAftParents, TotalNumOfPointSourceAfts, 
NumOfPointSourceParents, StartAft, EndAft, ParentIDs, ParentTimes, 
xParentLocs, yParentLocs, zParentLocs);
}
