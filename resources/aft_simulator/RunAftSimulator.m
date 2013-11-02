%RunAftSimulator.m
%This code loads sample data, sets parameters, and runs the AftSimulator.m
%program.  Parameters for AftSimulator.m may be easily changed by changing them 
%within this code & saving.

%NOTE: The Earthquake Background Rate density grid for this sample file was
%prepared by Andy Michael.  Please reference:

%Hardebeck, Jeanne L., Karen R. Felzer, and Andrew J. Michael, Rigorous
%Observational Tests Contradict the Accelerating Moment Release Hypothesis,
%Journal of Geophysical Research, submitted, 2007. 


load SampleInput

%SampleInput Contains the following data:

%MeasCat: Southern California (south of 36.5 N) catalog from 1/1/1984 -
%12/31/2006, earthquakes M>=2.5.

%BkgndRateGrid_AllCal: Background earthquake rate grid for all of
%California

%BkgndRateGrid_SoCal: Background earthquake rate grid for Southern
%California only

%MeasCatPlanarSourceParams: Fault parameters for M>=6.5 earthquaeks in
%Southern California, 1932 - 2006.

%Below, setting parameters needed by the program

BkgndRateGrid = BkgndRateGrid_SoCal;  %Change to the _AllCal matrix for all of California

SimStartDate = 'January 1, 1984';
ReportStartDate = 'January 1, 2005';
ReportEndDate = 'July 15, 2005';

LowerLimitOnActiveQuakeMags = 2.5;
LowerLimitOnPlanarSourceMags = 6.5;
LowerLimitOnReportMags = 2.5;

UpperLimitOnSimMags = 7.9;

LowerLimitOnAftDistFromParent = 0.001;  %In km
UpperLimitOnAftDistFromParent = 500;   %In km

AftDecayWithDist=1.37;

CX=0.095;
P=1.3400;
A=0.0080;

[ActiveQuakes, ReportCat] ... 
      = AftSimulator(MeasCat, MeasCatPlanarSourceParams, BkgndRateGrid, ...
           SimStartDate, ReportStartDate, ReportEndDate, ... 
           LowerLimitOnActiveQuakeMags, LowerLimitOnPlanarSourceMags, ... 
           LowerLimitOnReportMags, UpperLimitOnSimMags,  ... 
           LowerLimitOnAftDistFromParent,UpperLimitOnAftDistFromParent, ... 
           AftDecayWithDist,CX, P, A);

