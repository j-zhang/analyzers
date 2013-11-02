function [ActiveQuakes, ReportCat]  ... 
      = AftSimulator(MeasCat, MeasCatPlanarSourceParams, BkgndRateGrid, ...
           SimStartDate, ReportStartDate, ReportEndDate, ... 
           LowerLimitOnActiveQuakeMags, LowerLimitOnPlanarSourceMags, ... 
           LowerLimitOnReportMags, UpperLimitOnSimMags,  ... 
           LowerLimitOnAftDistFromParent,UpperLimitOnAftDistFromParent, ... 
           AftDecayWithDist,CX, P, A); 
%This function simulates aftershocks of measured and background quakes and
%then simulates aftershocks of the aftershocks and so on to build complete
%simulated catalogs of background quakes and aftershocks.


%Written by Karen Felzer with help from Alan Felzer.  Latest revisions,
%October 31, 2007.

  
%REFERENCES

%Felzer, K. R., T. W. Becker, R. E. Abercrombie, G. Ekstršm, and J. R.
%Rice, Triggering of the 1999 Mw 7.1 Hector Mine earthquake by aftershocks
%of the 1992 Mw 7.3 Landers earthquake, J. Geophys. Res., 107, 2190,
%doi:10.1029/2001JB000911, 2002. 

%Hardebeck, Jeanne L., Karen R. Felzer, and Andrew J. Michael, Rigorous
%Observational Tests Contradict the Accelerating Moment Release Hypothesis,
%Journal of Geophysical Research, submitted, 2007. 

%Spatial aftershock decay with distance used in this program based on
%findings in:

%Felzer, K. R. and E. E. Brodsky, Decay of aftershock density with distance 
%indicates triggering by dynamic stress, Nature , 441, 735-738, 2006. 

%Statistically simulating aftershock sequences complete with full secondary
%aftershock triggering pioneered by Ogata (1988):

%Ogata, Y., Statistical models for earthquake occurrence and residual
%analysis for point processes", J. Am. Stat. Assoc., 83, 9-27, 1988.


%PROGRAM FILES

%There are four files:
%
% (1) The Matlab function AftSimulator.m (2) SampleInput.mat: File with
% sample data input, loaded by RunAftSimulator.m (see below). (3)
% RunAftSimulator.m which loads sample input data, sets parameters, and
% runs AftSimulator.m.  Parameters can be changed easily by changing them
% within AftSimulator.m. (4) A C program ParentInfo.c to speed up the
% execution of a for loop in
%     the subfunction AftsOfPointSourceQuakes
%
%If you wish to use the C program you need to uncomment the corresponding
%code in the subfunction AftsOfPointSourceQuakes and then compile
%ParentInfo.c for your computer by typing in the Matlab command window
%
%        >> mex ParentInfo.c
%
%The name of the compiled program will depend on your computer
%
%  ParentInfo.mexmac  for Macs ParentInfo.mexglx  for Linux ParentInfo.dll
%  for PCs
%
%See the function AftsOfPointSourceQuakes for more details 


%PROGRAM OVERVIEW

%This program uses empirical aftershock relationships and Monte Carlo
%methods to simulate background earthquakes and their aftershocks as well
%as the aftershocks of real earthquakes in a measured catalog. It then
%simulates aftershocks of the aftershocks and so on. Since all the
%generated quakes come from Monte Carlo simulations the program should be
%run multiple times if you wish to generate averages and error bars on
%earthquake rates.

%The final report contains all earthquakes that occur between the 
%ReportStartDate and ReportEndDate that have magnitudes larger than 
%LowerLimitOnReportMags.

%The magnitudes, times and location of each simulated aftershock are 
%randomly chosen according to the following empirical distributions:
%
%  (1) Gutenberg-Richter for magnitudes (2) Omori's Law for times (3) r^-n
%  for distances of aftershocks from hypocenters (Felzer-Brodsky).
%      The average value of n for Southern California is 1.37.  Note that
%      we refer to n in the program as AftDecayWithDist (one of the input
%      variables).

%Note that quakes referred to in this program as "active quakes" are those
%that are eligible to trigger simulated aftershocks, which means that they
%occur within a specified time range and are greater than or equal to
%LowerLimitOnActiveQuakeMags, a magnitude that should be chosen as small as
%possible so as to give reasonable results without taking unreasonable
%computation time.


%PROGRAM ORGANIZATION

%First we find the active quakes of the measured (real or pre-specified)
%catalog. These are all quakes from the measured catalog that occur after
%SimStartDate but before ReportStartDate and have magnitudes greater than
%LowerLimitOnActiveQuakesMags.

%The program then generates background quakes in the time interval from
%ReportStateDate to ReportEndDate. The locations of the background
%eathquakes are chosen randomly from densities in spatial boxes specified
%by the user in the BkgrndRateGrid input. Note that all background quakes
%have magnitudes large enough to be active and that all background quakes
%with magnitudes greater than greater than LowerLimitOnReportMags become
%quakes in the final Report. 

%The program then does the following:
%
%  (1) It simulates aftershocks of the active measured and background
%  quakes
%      and adds to the Report all aftershocks that occur between
%      ReportStartDate and ReportEndDate that have magnitudes greater than
%      LowerLimitOnReportMags
%
%  (2) It finds the active aftershocks - the aftershocks that occur between
%      ReportStartDate and ReportEndDate that have magnitudes greater than
%      LowerLimitOnActiveMags
%
%  (3) It simulates aftershocks of the active aftershocks
%
%  (4) It then simulates aftershocks of the new active aftershocks,
%  continuing this 
%      loop until there are no more active aftershocks.



%Note that quakes with magnitudes less than LowerLimitOnPlanarSourceQuakes 
%are treated as point sources and those with larger magnitudes as planar 
%sources. When the middle of the fault plane of an earthquake from the
%measured catalog is not well constrained we use the median of the first
%two days of the aftershock sequence to find it. For earthquakes larger
%than LowerLimitOnPlanarSourceQuakes generated within the program, we
%currently assign 75% 332 degree strike and 25% 242 degree strike, based on
%approximate assesement of trends in the Southern California catalog, and a
%90 degree dip.

%And finally the program not only outputs the final ReportCat but also a 
%catalog of all of the quakes that are active (allowed to generate their
%own aftershocks at any given time).  All earthquakes in the final report
%(ReportCat) have, in column 11, the number of their parent mainshock in
%the order listed in the ActiveQuakes output.  For example, if the value
%given in column 11 of ReportCat is 10, this means that that earthquake was
%triggered by the tenth earthquake listed in ActiveQuakes. A zero will be
%given in column 11 if the earthquake is generated as a background
%earthquake rather than an aftershock.

%------------------------ INPUT VARIABLES --------------------------% 

%MeasCat: Pre-existing measured earthquakes, presumably from a downloaded
%earthquake catalog .  Can also be just one specific earthquake, or
%scenario earthquake or earthquake(s) that you wish to simulate aftershocks
%for.  MeasCat needs to be in ten column format: year, month, day, hour,
%minute, second, lat, lon, depth, magnitude.  All columns should be
%separated by spaces, and there can be no non-numerical entries such as :
%or /.  If you do not want aftershocks of any pre-existing earthquakes in
%the simulation then just enter MeasCat as an empty matrix: [].

%MeasCatPlanarSourceParams: Seven column format: fault length, width,
%strike, dip, year, month, day.  At least one row for each earthquake in
%MeasCat that is >=LowerLimitOnPlanarSourceQuakes.  if there are no such
%earthquakes this can be input as an empty matrix, []. Listings should be
%in chronological order, not include earthquakes smaller than
%LowerLimitOnPlanarSourceQuakes, cover an area that is the same as the area
%covered by MeasCat, and cover at least as much time.


%BkgndRateGrid: Spatially varying background rate. Five column format:  
%minlon, maxlon, minlat, maxlat, earthquake rate of M>=4 per year in the
%box thus specified.  Give one row for each box that is to contain
%background seismicity.  If you do not want any background earthquakes in
%the simulation then enter BkgndRateGrid as an empty matrix: [].

%SimStartDate: Where we begin using data from the measured catalog, to be
%entered as a string, such as 'January 1, 2007', or a vector with year,
%month, day, hour, minute, second, such as [2007 1 1 0 30 29].

%ReportStartDate: First date of the ReportCat, to be entered as a string or
%vector, as above.

%ReportEndDate: Last date of the ReportCat, to be entered as a string or
%vector, as above.

%LowerLimitOnActiveQuakeMags: No active quake can have a smaller mag

%LowerLimitOnReportMags: No quake in the report can have a smaller mag
        
%UpperLimitOnSimMags: No quake generated in this program can have a larger 
%magnitude

%LowerLimitOnPlanarSourceMags: All quakes with at least this mag are
%modeled as planar source quakes
        
%UpperLimitOnAftDistFromParent: No aftershock can be farther away from
%their mainshock (their parent).

%LowerLimitOnAftDistFromParent: No aftershock can be closer to their
%mainshock (their parent).  Note that this value must be set at >0 to
%prevent the inverse power law from becoming unconstrained at r=0.

%AftDecayWithDist: This is the variable n in the equation r-n for 
%aftershock density as a function of distance from the mainshock.
        
%CX, P, A: Omori's law parameters, c, p, and a (or k). They should be 
%entered as direct triggering parameters, WHICH ARE DIFFERENT from
%parameters observed in the field for full (direct plus secondary)
%aftershock sequences.  See Felzer et al. (2002) and Felzer et al. (2003)
%for further details. The value of A should be appropriate for the rate of
%a mainshock of magnitude M triggering aftershocks of magnitude >= M.  The
%values of the parameters will also need to be different depending on the
%value of ActiveQuakeMinMag chosen.  

%RECOMMENDATIONS:  

% CX=0.095, P=1.34, and A=0.008 for ActiveQuakeMinMag=2.5 CX=0.09, P=1.34,
% and A=0.0068 for ActiveQuakeMinMag=0.5

%These values correspond to averages for Southern California. (from Felzer
%and Kilb, in preparation, 2007). Please note that if an innapropriate mix
%of parameters are chosen (A too large, and/or CX and/or P too small),
%aftershock production may become super-critical, meaning that on average
%(averaged accross all mainshock magnitudes), there is more than one
%aftershock produced per mainshock.  In this case aftershock production may
%accelerate and continue indefinitely, and the program will get stuck in an
%infinite loop.  As the program runs values are output to the screen
%indicating the current number of aftershocks being generated; if these
%numbers increase continuously and quickly then super-critical parameters
%have probably been input.



%---------------------- OUTPUT VARIABLES ---------------------------%

%ActiveQuakes: Six column format: time, xyz locations, mag, IDs

%ReportCat: Twelve column format: Columns 1-10: year, month, day, hour,
%min, sec, lat, lon, depth, mag Column 11: ID numbers of aftershock parents
%
%Column 12: Distance bewteen aftershocks and the nearest point on the fault
%plane of their parents


%-----------------OTHER CATALOGS IN THE PROGRAM---------------------%

%PlanarSourceParams: Ten column format: xyz location, fault length, width,
%strike, dip, year, month, day



%-----------------------BEGINNING OF CODE---------------------------%
    
%ACTIVE QUAKES IN THE MEASURED CATALOG
if isempty(MeasCat)
    MeasActiveQuakes = [];
    MeasActivePlanarSourceParams = [];
    CenterOfActiveQuakeLocs = [];
end

if ~isempty(MeasCat)  
    [MeasActiveQuakes,MeasActivePlanarSourceParams, ... 
     CenterOfActiveQuakeLocs] ... 
        = GetActiveQuakesInMeasCat(MeasCat, MeasCatPlanarSourceParams, ... 
          SimStartDate, ReportStartDate,  ... 
          LowerLimitOnActiveQuakeMags,LowerLimitOnPlanarSourceMags);
end


%Freeing up memory space
clear MeasCat
clear MeasCatPlanarSourceParams


%GENERATION OF BACKGROUND QUAKES FOR THE TIME OF THE REPORT. NOTE THAT ALL
%BACKGROUND QUAKES ARE ACTIVE.
if isempty(BkgndRateGrid)
    BkgndQuakes = [];
    BkgndPlanarSourceParams = [];
    BkgndReport = [];
end


if ~isempty(BkgndRateGrid)
[BkgndQuakes,BkgndPlanarSourceParams,BkgndReport] ...  
    = GenerateBkgndCat(BkgndRateGrid,ReportStartDate,ReportEndDate, ... 
        CenterOfActiveQuakeLocs,LowerLimitOnActiveQuakeMags, ... 
        UpperLimitOnSimMags,LowerLimitOnReportMags, ... 
        LowerLimitOnPlanarSourceMags);    
end


%Freeing up memory space
clear BkgndRateGrid


%COMBINING THE BACKGROUND AND ACTIVE MEASURED QUAKES
MeasAndBkgndActiveQuakes = [MeasActiveQuakes; BkgndQuakes];
MeasAndBkgndActivePlanarSourceParams ... 
    	= [MeasActivePlanarSourceParams; BkgndPlanarSourceParams];
    
SizeOfMeasAndBkgndActiveQuakes = size(MeasAndBkgndActiveQuakes);
NumOfMeasAndBkgndActiveQuakes = SizeOfMeasAndBkgndActiveQuakes(1);


%Freeing up memory space
    clear MeasActiveQuakes
    clear BkgndQuakes
    clear MeasActivePlanarSourceParams
    clear BkgndPlanarSourceParams


%WHEN THERE ARE NO ACTIVE QUAKES
if NumOfMeasAndBkgndActiveQuakes == 0
    disp('There are no active quakes to generate aftershocks')
end

%WHEN THERE ARE ACTIVE QUAKES WE ASSIGN THEM IDs AND THEN RUN ETAS TO GET
%THEIR AFTERSHOCKS AND THEN THE AFTERSHOCKS OF THEIR AFTERSHOCKS AND SO ON 
if NumOfMeasAndBkgndActiveQuakes > 0   
    MeasAndBkgndActiveQuakeIDs = [1:NumOfMeasAndBkgndActiveQuakes]';    
    
    MeasAndBkgndActiveQuakes ... 
        = [MeasAndBkgndActiveQuakes, MeasAndBkgndActiveQuakeIDs];  
    

    %NOW FOR ETAS
    [ActiveQuakes, ReportCat] ... 
        = RunETAS(MeasAndBkgndActiveQuakes, ... 
            MeasAndBkgndActivePlanarSourceParams,BkgndReport, ...
            SimStartDate,ReportStartDate,ReportEndDate, ... 
            CenterOfActiveQuakeLocs,LowerLimitOnActiveQuakeMags, ...
            UpperLimitOnSimMags,LowerLimitOnAftDistFromParent, ... 
            UpperLimitOnAftDistFromParent,LowerLimitOnReportMags, ...
            LowerLimitOnPlanarSourceMags,AftDecayWithDist,CX,P,A);
        
end

%PRINT
SizeOfReportCat = size(ReportCat)
MaxQuakeInReport = max(ReportCat(:,10))

end   %AftSimulator



%-----------------------------------------------------------------------%
%----------------Function GetActiveQuakesInMeasCat----------------------%
%-----------------------------------------------------------------------%

function [MeasActiveQuakes, MeasActivePlanarSourceParams, ... 
CenterOfActiveQuakeLocs] ... 
   = GetActiveQuakesInMeasCat(MeasCat, MeasCatPlanarSourceParams,  ... 
       SimStartDate,ReportStartDate, LowerLimitOnActiveQuakeMags, ... 
       LowerLimitOnPlanarSourceMags);
%This function finds the active quakes in the measured catalog - the quakes
%that occur after SimStartDate, before ReportStartDate and that have 
%magnitudes greater than or equal to LowerLimitOnActiveQuakeMag - the 
%lowest magnitude quake we ask ETAS to generate aftershocks for


%FINDING THE ACTIVE MEASURED QUAKES

%Start times in Julian days from 0/0/0000
SimStartTime=datenum(SimStartDate);
ReportStartTime=datenum(ReportStartDate);

%Times of meas cat quakes in Julian days from 0/0/0000
TimesOfMeasQuakes=datenum(MeasCat(:,1),MeasCat(:,2),MeasCat(:,3), ... 
    MeasCat(:,4),MeasCat(:,5),MeasCat(:,6));

%Magnitudes of meas cat quakes
MagsOfMeasQuakes=MeasCat(:,10);

%Active quakes in the measured catalog
IndicesOfActiveQuakes=find(TimesOfMeasQuakes>=SimStartTime  ... 
    & TimesOfMeasQuakes<=ReportStartTime  ... 
    & MagsOfMeasQuakes>=LowerLimitOnActiveQuakeMags);

%Number of measured catalog active quakes
SizeOfMeasCatActiveQuakes = size(IndicesOfActiveQuakes);
NumOfMeasCatActiveQuakes = SizeOfMeasCatActiveQuakes(1);


%When there are no active quakes
if NumOfMeasCatActiveQuakes==0
    MeasActiveQuakes = [];
    MeasActivePlanarSourceParams = [];
    CenterOfActiveQuakeLocs = [];
end


%When there are active quakes we put their parameters in a 5 column format
if NumOfMeasCatActiveQuakes > 0
    MeasCatActiveQuakes=MeasCat(IndicesOfActiveQuakes,:);

    %Times of active quakes in the measured catalog
    TimesOfMeasActiveQuakes=TimesOfMeasQuakes(IndicesOfActiveQuakes);

    %Magnitudes of active quakes in the measured catalog
    MagsOfMeasActiveQuakes=MagsOfMeasQuakes(IndicesOfActiveQuakes);

    %Changing lat, lon and depth to xyz coordinates
    CenterOfActiveQuakeLocs(1) = mean(MeasCatActiveQuakes(:,7));
    CenterOfActiveQuakeLocs(2) = mean(MeasCatActiveQuakes(:,8));
    xyLocsOfMeasActiveQuakes=latlon2xy2(MeasCatActiveQuakes(:,7:8), ... 
        CenterOfActiveQuakeLocs);
    zLocsOfMeasActiveQuakes = MeasCatActiveQuakes(:,9);
    xyzLocsOfMeasActiveQuakes=[xyLocsOfMeasActiveQuakes, ... 
        zLocsOfMeasActiveQuakes];

    %Measured active quakes in a 5 column format: 
    MeasActiveQuakes = [TimesOfMeasActiveQuakes,  ... 
        xyzLocsOfMeasActiveQuakes,MagsOfMeasActiveQuakes];

end


%FINDING THE ACTIVE MEASURED PLANAR SOURCE QUAKES

SizeOfMeasCatPlanarSourceQuakes = size(MeasCatPlanarSourceParams);
NumOfMeasCatPlanarSourceQuakes = SizeOfMeasCatPlanarSourceQuakes(1);


%When there are planar source quakes we look to see if any are active - if
%they fall in the right time interval
if NumOfMeasCatPlanarSourceQuakes > 0 
    %Times of planar source quakes in Julian days from 0/0/0000
    MeasCatPlanarSourceTimes=datenum(MeasCatPlanarSourceParams(:,5), ... 
        MeasCatPlanarSourceParams(:,6),MeasCatPlanarSourceParams(:,7));

    %Measured catalog active planar source quakes
    ActivePlanarSourceQuakes = ... 
        find(MeasCatPlanarSourceTimes>=SimStartTime ... 
            & MeasCatPlanarSourceTimes<=ReportStartTime);
        
    SizeOfActivePlanarSourceQuakes = size(ActivePlanarSourceQuakes);
    NumOfActivePlanarSourceQuakes = SizeOfActivePlanarSourceQuakes(1);
end


%ADDING XYZ LOCATIONS TO ACTIVE PLANAR SOURCE QUAKE SOURCE PARAMETERS

%When there are no planar source quakes 
if NumOfMeasCatPlanarSourceQuakes==0 | NumOfActivePlanarSourceQuakes==0
    MeasActivePlanarSourceParams = [];
    NumOfActivePlanarSourceQuakes = 0;
end

%When there are active planar source quakes
if NumOfActivePlanarSourceQuakes>0
    MeasCatActivePlanarSourceParams = ...
        MeasCatPlanarSourceParams(ActivePlanarSourceQuakes,:);

    %Adding the xyz locations of the planar source quakes - the medians of
    %the first two days of quakes located within 20 km of the epicenter -
    %to their fault parameters matrix

    PlanarSourceQuakes ...
        = find(MagsOfMeasActiveQuakes>=LowerLimitOnPlanarSourceMags);

    for j=1:NumOfActivePlanarSourceQuakes
        TimeDiffs = TimesOfMeasActiveQuakes  ...
            - TimesOfMeasActiveQuakes(PlanarSourceQuakes(j));

        LocDiffsX = xyzLocsOfMeasActiveQuakes(:,1)  ...
            - xyzLocsOfMeasActiveQuakes(PlanarSourceQuakes(j),1);
        LocDiffsY = xyzLocsOfMeasActiveQuakes(:,2) ...
            - xyzLocsOfMeasActiveQuakes(PlanarSourceQuakes(j),2);

        aftx = find(TimeDiffs>=0 & TimeDiffs<2 & LocDiffsX<=25 ...
            & LocDiffsY<=25);

        %Setting all of the fault center depths at 10 km -- especially for
        %the older earthquakes the depths are unreliable.

        if(length(aftx)>20)
            %Probably a decent aftershock sequence not dominated by
            %background quakes
            xyzActivePlanarSourceLocs(j,:) ...
                = [median(xyzLocsOfMeasActiveQuakes(aftx,1)), ...
                  median(xyzLocsOfMeasActiveQuakes(aftx,2)), 10];
        else
            xyzActivePlanarSourceLocs(j,:) ...
                = xyzLocsOfMeasActiveQuakes(PlanarSourceQuakes(j),:);
        end
    end

    %Fault parameters of active planar sources in a 10 column format: xyz
    %locations, length, width, strike, dip, year, month, day
    MeasActivePlanarSourceParams ...
        = [xyzActivePlanarSourceLocs, MeasCatActivePlanarSourceParams];

end

end  % GetMeasCatActiveQuakes



%-----------------------------------------------------------------------%
%--------------------Function GenerateBkgndCat--------------------------%
%-----------------------------------------------------------------------%

function [BkgndQuakes,BkgndPlanarSourceParams,BkgndReport] ... 
    = GenerateBkgndCat(BkgndRateGrid,ReportStartDate,ReportEndDate, ... 
        CenterOfActiveQuakeLocs,LowerLimitOnActiveQuakeMags, ... 
        UpperLimitOnSimMags,LowerLimitOnReportMags, ...
        LowerLimitOnPlanarSourceMags);
%This function generates background quakes and their fault parameters if
%they are large enough.


%Times in Julian days from 0/0/000
ReportStartTime = datenum(ReportStartDate);
ReportEndTime = datenum(ReportEndDate);
ReportTime=ReportEndTime-ReportStartTime;

%Changing background grid locations from latlon to xy coordinates
latlon1 = [BkgndRateGrid(:,3) BkgndRateGrid(:,1)];
latlon2 = [BkgndRateGrid(:,4) BkgndRateGrid(:,2)];


xycoords1 = latlon2xy2(latlon1,CenterOfActiveQuakeLocs); 
xycoords2 = latlon2xy2(latlon2,CenterOfActiveQuakeLocs); 

BkgndRateGrid(:,1) = xycoords1(:,1);
BkgndRateGrid(:,2) = xycoords2(:,1);
BkgndRateGrid(:,3) = xycoords1(:,2);
BkgndRateGrid(:,4) = xycoords2(:,2);


%Total number of background quakes in this run of the simulation
RateOfM4sInBkgndGrids=BkgndRateGrid(:,5);
AveNumOfBkgndQuakes=sum(RateOfM4sInBkgndGrids) ... 
    *(ReportTime/364.25)*10^(4 - LowerLimitOnActiveQuakeMags);  
NumOfBkgndQuakesInThisSim = poissinv(rand(1),AveNumOfBkgndQuakes); 

if NumOfBkgndQuakesInThisSim==0
    BkgndQuakes = [];
    BkgndPlanarSourceParams = [];
    BkgndReport = [];
end

if NumOfBkgndQuakesInThisSim > 0   
    %Probability a given background quake is in a given grid
    BkgndProbGrid = RateOfM4sInBkgndGrids/sum(RateOfM4sInBkgndGrids);

    %Cumulative sum of grid probabilities
    CumSumOfBkgndProbs = cumsum(BkgndProbGrid);

    %Making sure that CumSumOfBkgndProbs(end)=1 despite roundoff error
    CumSumOfBkgndProbs(end)=1;

    %Generating a random number between 0 and 1 for each background quake
    r=rand(NumOfBkgndQuakesInThisSim,1);

    %Randomly putting background quakes in the grids in proportion to the
    %probabilities that the grids have background quakes
    for k=1:NumOfBkgndQuakesInThisSim
        rk = find(CumSumOfBkgndProbs>=r(k));
        BkgndQuakeLoc(k) = rk(1);
    end


    %xyz locations of the background quakes in each grid cell.
    xdiff = BkgndRateGrid(BkgndQuakeLoc,2)-BkgndRateGrid(BkgndQuakeLoc,1);
    ydiff = BkgndRateGrid(BkgndQuakeLoc,4)-BkgndRateGrid(BkgndQuakeLoc,3);
   
    R = rand(NumOfBkgndQuakesInThisSim,3);
    xd = R(:,1).*xdiff + BkgndRateGrid(BkgndQuakeLoc,1);
    yd = R(:,2).*ydiff + BkgndRateGrid(BkgndQuakeLoc,3);
    zd = R(:,3)*20;  %random depth between 0 and 20

    xyzLocsOfBkgndQuakes = [xd yd zd];


    %Times of the background quakes
    TimesOfBkgndQuakes = rand(NumOfBkgndQuakesInThisSim,1) ... 
        *(ReportTime) + ReportStartTime;
    TimesOfBkgndQuakes = sort(TimesOfBkgndQuakes);


    %Magnitudes of the background quakes
    MagsOfBkgndQuakes ... 
       = GetMags(NumOfBkgndQuakesInThisSim,LowerLimitOnActiveQuakeMags, ... 
            UpperLimitOnSimMags);
     

    %Background active quakes in a 5 column format: 
    BkgndQuakes = [TimesOfBkgndQuakes, xyzLocsOfBkgndQuakes, ... 
        MagsOfBkgndQuakes];
end


%BACKGROUND QUAKE PLANAR SOURCE FAULT PARAMETERS

if NumOfBkgndQuakesInThisSim > 0 
    %Number of planar source background quakes
    PlanarSourceBkgndQuakes ... 
        = find(MagsOfBkgndQuakes>=LowerLimitOnPlanarSourceMags);
    
    SizeOfPlanarSourceBkgndQuakes = size(PlanarSourceBkgndQuakes);
    NumOfPlanarSourceBkgndQuakes = SizeOfPlanarSourceBkgndQuakes(1);

    if NumOfPlanarSourceBkgndQuakes==0
        BkgndPlanarSourceParams = [];
    end
    
    if NumOfPlanarSourceBkgndQuakes>0    
        xyzLocsOfPlanarSourceBkgndQuakes ... 
            = xyzLocsOfBkgndQuakes(PlanarSourceBkgndQuakes,:);
        MagsOfPlanarSourceBkgndQuakes ... 
            = MagsOfBkgndQuakes(PlanarSourceBkgndQuakes); 
                   
        [BkgndPlanarSourceParams] ... 
            = BuildPlanarSourceParams(xyzLocsOfPlanarSourceBkgndQuakes, ... 
                MagsOfPlanarSourceBkgndQuakes);  
    end
end


%INITIALIZING THE FINAL REPORT WITH QUAKES FROM THE BACKGROUND
if NumOfBkgndQuakesInThisSim > 0
    BkgndQuakesInReport ... 
        = find(MagsOfBkgndQuakes >= LowerLimitOnReportMags);
    NumOfBkgndQuakesInReport=length(BkgndQuakesInReport);

    if NumOfBkgndQuakesInReport==0
        BkgndReport = [];
    end
    
    if NumOfBkgndQuakesInReport>0
        ReportQuakeTimes = TimesOfBkgndQuakes(BkgndQuakesInReport);
        ReportQuakeLocs = xyzLocsOfBkgndQuakes(BkgndQuakesInReport,:);
        ReportQuakeMags = MagsOfBkgndQuakes(BkgndQuakesInReport);
        ParentsOfReportAfts = zeros(NumOfBkgndQuakesInReport,1);
        ReportAftDistancesFromParents=-ones(NumOfBkgndQuakesInReport,1);

        BkgndReport=[ReportQuakeTimes,ReportQuakeLocs,ReportQuakeMags, ...
            ParentsOfReportAfts,ReportAftDistancesFromParents];
    end
end
 
end   %GenerateBkgndCat



%-----------------------------------------------------------------------%
%---------------------------Function RunEtas----------------------------%
%-----------------------------------------------------------------------%

function [ActiveQuakes, ReportCat] ... 
    = RunETAS(MeasAndBkgndActiveQuakes, ... 
        MeasAndBkgndActivePlanarSourceParams,BkgndReport, ...
        SimStartDate,ReportStartDate,ReportEndDate, ... 
        CenterOfActiveQuakeLocs,LowerLimitOnActiveQuakeMags, ...
        UpperLimitOnSimMags,LowerLimitOnAftDistFromParent, ... 
        UpperLimitOnAftDistFromParent,LowerLimitOnReportMags, ...
        LowerLimitOnPlanarSourceMags,AftDecayWithDist,CX,P,A);
%This function simulates aftershocks of the active measured and background
%quakes and then simulates aftershocks of the aftershocks and so on.


% 1. Use the inverse Poisson function to generate aftershocks of the
%    active quakes 
% 2. Add to the simulation catalog all aftershocks that occur in the sim 
%    time with magnitudes >= LowerLimitOnReportMag
% 3. Repeat the process with all the newly generated aftershocks from 
%    Step 1 that are active quakes
% 4. Continue repeating until there are no more new aftershocks that are
%    active


%Times in Julian days from 0/0/0000
SimStartTime=datenum(SimStartDate);
ReportStartTime=datenum(ReportStartDate);
ReportEndTime=datenum(ReportEndDate);
SimTime=ReportEndTime-SimStartTime;

%Initializing the report and then freeing up memory space
Report = BkgndReport;
clear BkgndReport

%The average number of aftershocks expected during the simulation time with
%magnitudes greater than or equal to their mainshock
if(P ~= 1)
    AftProductivity = (A/(1-P))*((SimTime+CX)^(1-P) - CX^(1-P));
else
    AftProductivity = A*(log(SimTime+CX) - log(CX));
end  

%Number of active quakes from the MeasuredAndBkgndActiveQuakes
SizeOfActiveQuakes = size(MeasAndBkgndActiveQuakes);
NumOfActiveQuakes = SizeOfActiveQuakes(1);


%Initalizations for the first run through the while loop
NewActiveQuakes = MeasAndBkgndActiveQuakes;
NumOfNewActiveQuakes = NumOfActiveQuakes;
NewPlanarSourceParams = MeasAndBkgndActivePlanarSourceParams;

NumOfIterations = 0;


while NumOfNewActiveQuakes>0
    
    %Magnitudes of the new active quakes
    MagsOfNewActiveQuakes = NewActiveQuakes(:,5);
    
    %THE POINT SOURCE QUAKES
    IndicesOfPointSourceQuakes ... 
       = find(MagsOfNewActiveQuakes < LowerLimitOnPlanarSourceMags);
    PointSourceQuakes=NewActiveQuakes(IndicesOfPointSourceQuakes,:);
    NumOfPointSourceQuakes = length(IndicesOfPointSourceQuakes);
    
    %Print
    NumOfPointSourceQuakes
    
    
    %AFTERSHOCKS OF THE POINT SOURCE QUAKES
    if NumOfPointSourceQuakes==0
        PointSourceAfts = [];
    end
    
    if NumOfPointSourceQuakes>0    
       [PointSourceAfts] ... 
          = AftsOfPointSourceQuakes(PointSourceQuakes, ... 
              LowerLimitOnActiveQuakeMags,UpperLimitOnSimMags, ... 
              AftProductivity,SimTime,LowerLimitOnAftDistFromParent, ... 
              UpperLimitOnAftDistFromParent,AftDecayWithDist,CX,P,NumOfIterations);
    end
    
    
    %THE PLANAR SOURCE QUAKES
    IndicesOfPlanarSourceQuakes ... 
       = find(MagsOfNewActiveQuakes >= LowerLimitOnPlanarSourceMags);
    PlanarSourceActiveQuakes ... 
        = NewActiveQuakes(IndicesOfPlanarSourceQuakes,:);
    NumOfPlanarSourceQuakes = length(IndicesOfPlanarSourceQuakes);
       
    
    %Print
    NumOfPlanarSourceQuakes
    
    
    %AFTERSHOCKS OF THE PLANAR SOURCE QUAKES
    if NumOfPlanarSourceQuakes==0
        PlanarSourceAfts = [];
    end
    
    if NumOfPlanarSourceQuakes>0 
        %Generation of the planar source parameters
        if NumOfIterations>0
            xyzLocsOfPlanarSourceQuakes = PlanarSourceActiveQuakes(:,2:4);
            MagsOfPlanarSourceQuakes = PlanarSourceActiveQuakes(:,5); 
            
            [NewPlanarSourceParams] ... 
                = BuildPlanarSourceParams(xyzLocsOfPlanarSourceQuakes, ... 
                    MagsOfPlanarSourceQuakes);
        end

        %Generation of the aftershocks
        [PlanarSourceAfts] ... 
           = AftsOfPlanarSourceQuakes(PlanarSourceActiveQuakes, ... 
               NewPlanarSourceParams,LowerLimitOnActiveQuakeMags, ... 
               UpperLimitOnSimMags,AftProductivity,SimTime, ... 
               UpperLimitOnAftDistFromParent, ... 
               LowerLimitOnAftDistFromParent,AftDecayWithDist, CX,P);
    end


    %COMBINING THE POINT SOURCE AND PLANAR SOURCE AFTERSHOCKS
    NewAfts = [PointSourceAfts; PlanarSourceAfts];
    
    SizeOfNewAfts = size(NewAfts);
    NumOfNewAfts = SizeOfNewAfts(1);
    
    
    %PRINT
    NumOfNewAfts
    
    
    if NumOfNewAfts==0
        NumOfNewActiveQuakes = 0;
    end
            
    
    %THE NEW AFTERSHOCKS THAT ARE ACTIVE
    if NumOfNewAfts > 0
        TimesOfNewAfts = NewAfts(:,1);
        MagsOfNewAfts = NewAfts(:,5);
         
        IndicesOfNewActiveAfts = find(TimesOfNewAfts>=ReportStartTime ... 
            & TimesOfNewAfts<=ReportEndTime ...
            & MagsOfNewAfts>=LowerLimitOnActiveQuakeMags);
        NewActiveAfts = NewAfts(IndicesOfNewActiveAfts,:);   
        
        %The number of new active aftershocks 
        SizeOfNewActiveAfts = size(NewActiveAfts);
        NumOfNewActiveAfts = SizeOfNewActiveAfts(1);
        
        
        %PRINT
        NumOfNewActiveAfts
        
        
        %When there are no new active aftershocks then there are no new
        %active quakes
        if NumOfNewActiveAfts==0
            NumOfNewActiveQuakes = 0;
        end
        
    
        %When there are new active aftershocks they become the new active
        %quakes
        if NumOfNewActiveAfts > 0
            TimesOfNewActiveQuakes = NewActiveAfts(:,1);
            xyzLocsOfNewActiveQuakes = NewActiveAfts(:,2:4);
            MagsOfNewActiveQuakes = NewActiveAfts(:,5);
            NumOfNewActiveQuakes = NumOfNewActiveAfts;
            IDsOfNewActiveQuakes = [(NumOfActiveQuakes+1): ... 
                (NumOfActiveQuakes + NumOfNewActiveQuakes)]';
            
            %Combining the parameters of the new active quakes
            NewActiveQuakes = [TimesOfNewActiveQuakes, ...
                xyzLocsOfNewActiveQuakes,MagsOfNewActiveQuakes, ... 
                IDsOfNewActiveQuakes];
                        
                
            %COMBINING THE ACTIVE QUAKES WITH THE NEW ACTIVE QUAKES
            if NumOfIterations == 0
                ActiveQuakes = [MeasAndBkgndActiveQuakes; NewActiveQuakes];
            else
                ActiveQuakes = [ActiveQuakes; NewActiveQuakes];
            end
            
    
            %Total number of active quakes
            NumOfActiveQuakes = NumOfActiveQuakes + NumOfNewActiveQuakes;
    
    
            %THE NEW AFTERSHOCKS THAT GO IN THE REPORT
            ReportAfts = find(TimesOfNewAfts>=ReportStartTime ... 
                & TimesOfNewAfts<=ReportEndTime ... 
                & MagsOfNewAfts>=LowerLimitOnReportMags);
            NewReportAfts = NewAfts(ReportAfts,:);
    
            %Adding the new aftershocks to the report 
            Report = [Report; NewReportAfts];
    
        end 
    end      

    %Completed the iteration
    NumOfIterations = NumOfIterations + 1
    
    
    %print
    NumOfNewActiveQuakes

    
end  %while


%THE FINAL REPORT

%Report quake parameters
ReportQuakeTimes = Report(:,1);
xyzLocsOfReportQuakes = Report(:,2:4);
ReportQuakeMags = Report(:,5);
ParentsOfReportAfts = Report(:,6);
ReportAftDistFromParents = Report(:,7);

%Putting report quake times in traditional form
ReportQuakeTimes = datevec(ReportQuakeTimes);

%Putting report quake locations in traditional latlon form
LatLonReportLocs = xy2latlon(xyzLocsOfReportQuakes(:,1:2), ... 
    CenterOfActiveQuakeLocs);
Depth = xyzLocsOfReportQuakes(:,3);

%Putting the final report into a traditional earthquake catalog
ReportCat = [ReportQuakeTimes, LatLonReportLocs,Depth,ReportQuakeMags, ...
    ParentsOfReportAfts, ReportAftDistFromParents];

ReportCat = sortrows(ReportCat,[1 2 3 4 5 6]);

end  %RunEtas



%-----------------------------------------------------------------------%
%--------------------Function AftsOfPointSourceQuakes-------------------%
%-----------------------------------------------------------------------%

function [PointSourceAfts] ... 
        = AftsOfPointSourceQuakes(PointSourceQuakes, ... 
              LowerLimitOnActiveQuakeMags,UpperLimitOnSimMags, ... 
              AftProductivity,SimTime, LowerLimitOnAftDistFromParent, ... 
              UpperLimitOnAftDistFromParent,AftDecayWithDist,CX,P,NumOfIterations);
%This function generates aftershocks of the point source active quakes


PointSourceTimes = PointSourceQuakes(:,1);
xyzLocsOfPointSourceQuakes = PointSourceQuakes(:,2:4);
PointSourceMags = PointSourceQuakes(:,5);
PointSourceIDs = PointSourceQuakes(:,6);


%The average number of aftershocks produced by each point source quake  
AveNumOfAftsOfEachQuake ... 
   = AftProductivity.*10.^(PointSourceMags - LowerLimitOnActiveQuakeMags);

%Number of aftershocks produced by each active quake in this simulation
NumOfAftsOfEachQuake ... 
    = poissinv(rand(length(AveNumOfAftsOfEachQuake),1), ... 
        AveNumOfAftsOfEachQuake);
TotalNumOfPointSourceAfts=sum(NumOfAftsOfEachQuake); 


if TotalNumOfPointSourceAfts == 0
    PointSourceAfts = [];
end

if TotalNumOfPointSourceAfts>0   
    %Parents of the aftershocks
    PointSourceParents = find(NumOfAftsOfEachQuake>0);
    NumOfPointSourceParents = length(PointSourceParents);
    NumOfAftsOfEachPointSourceParent ... 
        = NumOfAftsOfEachQuake(PointSourceParents);

    
    %Print
    NumOfPointSourceParents       
        
        
    %Times, locations and IDs of point source parents - point source quakes
    %with aftershocks
    if NumOfPointSourceParents>0
        ParentTimes = PointSourceTimes(PointSourceParents);
        xParentLocs = xyzLocsOfPointSourceQuakes(PointSourceParents,1);
        yParentLocs = xyzLocsOfPointSourceQuakes(PointSourceParents,2);
        zParentLocs = xyzLocsOfPointSourceQuakes(PointSourceParents,3);
        ParentIDs = PointSourceIDs(PointSourceParents);
    end
    
    %Freeing up memory space
    clear PointSourceTimes
    clear xyzLocsOfPointSourceQuakes
    clear PointSourceIDs
    

    %A vector with the time, location and ID of each aftershock
    if NumOfPointSourceParents>0       
        CumSumOfAfts = cumsum(NumOfAftsOfEachPointSourceParent);
        StartAft = [1; CumSumOfAfts(1:end-1)+1];
        EndAft = CumSumOfAfts;

  
%-----------------------FOR LOOP IN MATLAB--------------------------%

        for j=1:NumOfPointSourceParents
           AftParentTimes(StartAft(j):EndAft(j),1) = ParentTimes(j);
           xAftParentLocs(StartAft(j):EndAft(j),1) = xParentLocs(j);
           yAftParentLocs(StartAft(j):EndAft(j),1) = yParentLocs(j);
           zAftParentLocs(StartAft(j):EndAft(j),1) = zParentLocs(j);
           PointSourceAftParents(StartAft(j):EndAft(j),1) = ParentIDs(j);
        end

%-------------------------------------------------------------------%
%
%To run the above for loop in C you need to
%  (1) Comment the above FOR LOOP IN MATLAB (2) Uncomment the below FOR
%  LOOP IN C (3) Compile the C function ParentInfo by typing in the command
%  window
%           
%                       >> mex ParentInfo.c
%
%-------------------------------------------------------------------%


%--------------------------FOR LOOP IN C----------------------------% 

   
%      [AftParentTimes, xAftParentLocs, yAftParentLocs, zAftParentLocs, ...  
%      PointSourceAftParents] ...
%        = ParentInfo(TotalNumOfPointSourceAfts, ... 
%          NumOfPointSourceParents,StartAft,EndAft,ParentIDs, ... 
%          ParentTimes,xParentLocs,yParentLocs,zParentLocs); 
      
%-------------------------------------------------------------------%
     

%print
if NumOfIterations > 0
    TotalNumOfPointSourceAfts, ... 
            NumOfPointSourceParents
        
        
end

        %Distances between aftershocks and their parents in km.
        AftDistFromPointSourceParents ... 
            = GetDistancesBetweenAftsAndParents(AftDecayWithDist, ... 
                TotalNumOfPointSourceAfts, ... 
                LowerLimitOnAftDistFromParent, ... 
                UpperLimitOnAftDistFromParent);
    
        %Converting the aftershock distances from their parents into x, y 
        %using a random theta.  Although not completely accurate we save 
        %computation time by making make the z distances of the aftershocks 
        %equal to those of their parents. As a result we don't need to  
        %worry about staying within seismogenic depth.  
        
        theta = rand(TotalNumOfPointSourceAfts,1)*2*pi;
     
        xAftDistFromParents = AftDistFromPointSourceParents.*cos(theta);
        yAftDistFromParents = AftDistFromPointSourceParents.*sin(theta);
        zAftDistFromParents = zeros(TotalNumOfPointSourceAfts,1);
    
    
        %xyz locations of the point source aftershocks
        xyzLocsOfPointSourceAfts ... 
            = [xAftParentLocs yAftParentLocs zAftParentLocs] ... 
                + [xAftDistFromParents, yAftDistFromParents, ... 
                zAftDistFromParents];
     
    
        %Times between aftershocks and their parents - GPowDistRc
        TimesBetweenAftsAndParents ... 
            = GetTimesBetweenAftsAndParents(P,CX, ... 
                TotalNumOfPointSourceAfts,0,SimTime);
            
    
        %Times of the point source aftershocks
        PointSourceAftTimes = AftParentTimes + TimesBetweenAftsAndParents;
        
    
        %Magnitudes of the point source aftershocks 
        PointSourceAftMags = GetMags(TotalNumOfPointSourceAfts, ... 
            LowerLimitOnActiveQuakeMags, UpperLimitOnSimMags);
    
    end

    PointSourceAfts =[PointSourceAftTimes,xyzLocsOfPointSourceAfts, ... 
        PointSourceAftMags,PointSourceAftParents, ... 
        AftDistFromPointSourceParents];
end   
    
end  %AftsOfPointSourceQuakes



%-----------------------------------------------------------------------%
%-------------------Function AftsOfPlanarSourceQuakes-------------------%
%-----------------------------------------------------------------------%

function [PlanarSourceAfts] ... 
           = AftsOfPlanarSourceQuakes(PlanarSourceActiveQuakes, ... 
               NewPlanarSourceParams,LowerLimitOnActiveQuakeMags, ... 
               UpperLimitOnSimMags,AftProductivity,SimTime, ... 
               UpperLimitOnAftDistFromParent, ... 
               LowerLimitOnAftDistFromParent, AftDecayWithDist, CX,P);
%This function generates aftershocks of the planar source active quakes

PlanarSourceTimes = PlanarSourceActiveQuakes(:,1); 
xyzLocsOfPlanarSourceQuakes = PlanarSourceActiveQuakes(:,2:4);
PlanarSourceMags = PlanarSourceActiveQuakes(:,5);
PlanarSourceIDs = PlanarSourceActiveQuakes(:,6);


%Number of planar source quakes
NumOfPlanarSourceQuakes = length(PlanarSourceTimes);

%Average number of aftershocks produced by each planar source quake  
AveNumOfAftsOfEachPlanarSourceQuake ... 
   =AftProductivity*(10.^(PlanarSourceMags - LowerLimitOnActiveQuakeMags));

%Number of aftershocks produced by each active quake in this simulation
NumOfAftsOfEachPlanarSourceQuake ... 
    = poissinv(rand(NumOfPlanarSourceQuakes,1), ... 
    AveNumOfAftsOfEachPlanarSourceQuake); 

%Total number of planar source aftershocks
TotalNumOfPlanarSourceAfts=sum(NumOfAftsOfEachPlanarSourceQuake)


%Aftershocks of the planar source quakes
if TotalNumOfPlanarSourceAfts==0
    PlanarSourceAfts = [];
end


if TotalNumOfPlanarSourceAfts>0
    %The planar source parents - the planar source quakes that have afts
    PlanarSourceParents = find(NumOfAftsOfEachPlanarSourceQuake>0);
    NumOfPlanarSourceParents = length(PlanarSourceParents);
    NumOfAftsOfEachPlanarSourceParent ... 
        = NumOfAftsOfEachPlanarSourceQuake(PlanarSourceParents);


    %XYZ locations of the parents of each aftershock.  These locations are
    %randomly chosen on the planar source faults and are the points from
    %from which the distances of the aftershocks are measured

    %Initialization
    xParentLocs = [];  yParentLocs = [];  zParentLocs = [];

    for j = 1:NumOfPlanarSourceParents
        FaultY = NewPlanarSourceParams(j,4);
        FaultZ = NewPlanarSourceParams(j,5);
    
        YLow = NewPlanarSourceParams(j,2) - FaultY/2;
        ZLow = NewPlanarSourceParams(j,3) - FaultZ/2;
    
        %Points on the planar source faults for each aftershock
        rpy = rand(NumOfAftsOfEachPlanarSourceParent(j),1).*FaultY;
        rpz = rand(NumOfAftsOfEachPlanarSourceParent(j),1)*FaultZ;
        rpx = zeros(NumOfAftsOfEachPlanarSourceParent(j),1);
    
        rpy = rpy - FaultY/2;
    
    
        %Matrix for rotation around strike
        z1 = NewPlanarSourceParams(PlanarSourceParents(j),6); 
        lambda1 = [cosd(z1) sind(z1); -sind(z1) cosd(z1)];

        %Matrix for rotation around dip
        if NewPlanarSourceParams(PlanarSourceParents(j),6) <= 180
            z2 = (NewPlanarSourceParams(PlanarSourceParents(j),7) - 90);
        else
            z2 = (90-NewPlanarSourceParams(PlanarSourceParents(j),7));
        end

        lambda2 = [cosd(z2) sind(z2); -sind(z2) cosd(z2)];
    

        %First rotating to a 0 degree strike and then to a 90 degree dip
        rot = lambda1*[rpx rpy]';
    
        rpx = rot(1,:)';
        rpy = rot(2,:)';  
        
        rot = lambda2*[rpz rpy]';
    
        rpz = rot(1,:)';
        rpy = rot(2,:)';
        
        Lx = rpx + NewPlanarSourceParams(PlanarSourceParents(j),1);
        Ly = rpy + NewPlanarSourceParams(PlanarSourceParents(j),2);
        Lz = rpz + NewPlanarSourceParams(PlanarSourceParents(j),3);
    
        %Adding Lx, Ly and Lz to the matrix of aftershock parent locations
        xParentLocs = [xParentLocs; Lx];
        yParentLocs = [yParentLocs; Ly];
        zParentLocs = [zParentLocs; Lz];

    end   %for loop


    %Times and IDs of planar source parents
    ParentTimes = PlanarSourceTimes(PlanarSourceParents);
    ParentIDs = PlanarSourceIDs(PlanarSourceParents);


    %Vectors for the time and ID of each planar source aftershock
    CumSumOfAfts = cumsum(NumOfAftsOfEachPlanarSourceParent);
    StartAft = [1; CumSumOfAfts(1:end-1)+1];
    EndAft = CumSumOfAfts;


    for j=1:NumOfPlanarSourceParents
        AftParentTimes(StartAft(j):EndAft(j),1) = ParentTimes(j);
        PlanarSourceAftParents(StartAft(j):EndAft(j),1) = ParentIDs(j);
    end
 

    %Distances between aftershocks and their parents. 
    AftDistFromPlanarSourceParents ... 
        = GetDistancesBetweenAftsAndParents(AftDecayWithDist, ... 
          TotalNumOfPlanarSourceAfts,LowerLimitOnAftDistFromParent, ... 
          UpperLimitOnAftDistFromParent);
    
    %Converting the aftershock distances from their parents into x,y using
    %a random theta.  Although not completely accurate we save computation
    %time by making make the z distances of the aftershocks equal to those
    %of the generating point on the mainshock fault. As a result we don't
    %need to worry about staying within seismogenic depth.  
        
    theta = rand(TotalNumOfPlanarSourceAfts,1)*2*pi;
     
    xAftDistFromParent = AftDistFromPlanarSourceParents.*cos(theta);
    yAftDistFromParent = AftDistFromPlanarSourceParents.*sin(theta);
    zAftDistFromParent = zeros(TotalNumOfPlanarSourceAfts,1);
    
    
    %xyz locations of the Planar source aftershocks
    xyzLocsOfPlanarSourceAfts ... 
       = [xParentLocs yParentLocs zParentLocs] ... 
         + [xAftDistFromParent yAftDistFromParent zAftDistFromParent];
     
    
    %Times between aftershocks and their parents 
    TimesBetweenAftsAndParents ... 
        = GetTimesBetweenAftsAndParents(P,CX, ... 
             TotalNumOfPlanarSourceAfts,0,SimTime);
    
    %Times of the planar source aftershocks
    PlanarSourceAftTimes = AftParentTimes + TimesBetweenAftsAndParents;
    
    
    %Magnitudes of the planar source aftershocks 
    PlanarSourceAftMags = GetMags(TotalNumOfPlanarSourceAfts, ... 
        LowerLimitOnActiveQuakeMags, UpperLimitOnSimMags);

    PlanarSourceAfts = [PlanarSourceAftTimes,xyzLocsOfPlanarSourceAfts, ... 
        PlanarSourceAftMags,PlanarSourceAftParents, ... 
        AftDistFromPlanarSourceParents];

end 

end  %AftsOfPlanarSourceQuakes



%-----------------------------------------------------------------------%
%----------------------------Function GetMags---------------------------%
%-----------------------------------------------------------------------%

function Mags = GetMags(NumOfAfts,MinMagOfActiveQuakes,MaxMagOfSimAft)

Mags = MinMagOfActiveQuakes - log10(rand(NumOfAfts,1));

MagsTooLarge = find(Mags > MaxMagOfSimAft);
NumOfMagsTooLarge = length(MagsTooLarge);

for i=1:NumOfMagsTooLarge
    
    while Mags(MagsTooLarge(i))>MaxMagOfSimAft
        Mags(MagsTooLarge(i)) = MinMagOfActiveQuakes - log10(rand(1));
    end
end
end   %GetMags



%-----------------------------------------------------------------------%
%-----------------Function BuildNewPlanarSourceParams-------------------%
%-----------------------------------------------------------------------%

function [NewPlanarSourceParams] ... 
    = BuildPlanarSourceParams(xyzLocs,Mags);
%This function generates fault parameters for simulated planar source
%quakes

NumOfPlanarSourceQuakes = length(Mags);

%Using WellsCopper to assign planar source lengths and widths with the
%following angles z1 and z2 in degrees 
z1 = 303;  z2 = 213;    
lambda1 = [cosd(z1) sind(z1); -sind(z1) cosd(z1)];
lambda2 = [cosd(z2) sind(z2); -sind(z2) cosd(z2)];

for k=1:NumOfPlanarSourceQuakes;
    [FaultY,FaultZ] = WellsCopper(Mags(k),1); 

    %Placing the hypocenter of each mainshock at a random point on the 
    %fault plane, and then determining the y, z limits of the fault plane.
    posY = rand(length(FaultY),1).*FaultY;
    posZ = rand(length(FaultY),1).*FaultZ;

    YLow = xyzLocs(k,2)-posY;
    ZLow = xyzLocs(k,3)-posZ;

    if ZLow<0  %Don't want any faults in the air!
       ZLow = 0;
    end
        
    midpoint = [0 FaultY/2 FaultZ/2];
     
    %Determining fault strike
    rs = rand(1);
    if (rs<=0.75)
        rot = lambda1*[midpoint(:,1) midpoint(:,2)]';
        s = 332;
    else 
        rot = lambda2*[midpoint(:,1) midpoint(:,2)]';
        s = 242;
    end
    
    rpx = rot(1,:)';
    rpy = rot(2,:)';
        
    midpoint(:,1:3) = [xyzLocs(k,1)+rpx xyzLocs(k,2)+rpy xyzLocs(k,3)];
        
    NewPlanarSourceParams(k,1:10) = [midpoint FaultY FaultZ s 90 0 0 0];
end

end   %BuildFaultParam



%-----------------------------------------------------------------------%
%------------Function GetDistancesBetweenAftsAndParents-----------------%
%-----------------------------------------------------------------------%

function DistancesBetweenAftsAndParents ... 
        = GetDistancesBetweenAftsAndParents(AftDecayWithDist,TotalNumOfAfts, ... 
          LowerLimitOnAftDistFromParent,UpperLimitOnAftDistFromParent);
%This function generates a distance distribution following an inverse power
%law with exponent m, contained within the LowerLimitOnAftDistFromParent
%and UpperLimitOnAftDistFromParent

p = rand(TotalNumOfAfts,1);

dmin = LowerLimitOnAftDistFromParent;   
dmax = UpperLimitOnAftDistFromParent;

a1 = p.*(dmax^(1-AftDecayWithDist) - dmin^(1-AftDecayWithDist));
DistancesBetweenAftsAndParents = (a1 + dmin^(1-AftDecayWithDist)).^(1/(1-AftDecayWithDist));  

end   %GetDistancesBetweenAftsAndParents



%-----------------------------------------------------------------------%
%---------------Function GetTimesBetweenAftsAndParents------------------%
%-----------------------------------------------------------------------%

function TimesBetweenAftsAndParents ... 
        = GetTimesBetweenAftsAndParents(P,CX,TotalNumOfAfts,tmin,tmax);
%This function generates random points from the function (t+c)^-p, between
%tmin and tmax.  0 may be entered for tmin.

r = rand(TotalNumOfAfts,1);

if(P~=1)
    
    a1 = (tmax + CX)^(1-P);
    a2 = (tmin + CX)^(1-P);
    a3 = r*a1 + (1-r)*a2;
    TimesBetweenAftsAndParents = a3.^(1/(1-P)) - CX;
    
else
    
    a1 = log(tmax+CX);
    a2 = log(tmin + CX);
    a3 = r*a1 + (1-r)*a2;
    TimesBetweenAftsAndParents = exp(a3) - CX;
     
end

end   %TimesBetweenAftsAndParents



%-----------------------------------------------------------------------%
%--------------------------Function xy2latlon---------------------------%
%-----------------------------------------------------------------------%

function cat = xy2latlon(cat,midpoint);

%This function takes a list of x,y and converts them to lat, lon given a
%lat lon centerpoint (given in midpoint).  cat is two columns, x and y. The
%new cat is the same thing in terms of lat and lon.

%NOTE: THIS PROGRAM IS ONLY MEANT TO BE USED FOR RELATIVELY SMALL AREAS
%(SAY 100 KM ACROSS!)


%km2deg
cat(:,2) = cat(:,2)/111.1949;

cat(:,2) = cat(:,2) + midpoint(1);


%km2deg
cat(:,1) = cat(:,1)/111.1949;

cat(:,1) = cat(:,1).*(1./cosd(cat(:,2)));

cat(:,1) = cat(:,1) + midpoint(2);

temp = cat(:,1);

cat(:,1)= cat(:,2);

cat(:,2)= temp;   

end   %xy2latlon



%-----------------------------------------------------------------------%
%--------------------------Function WellsCopper-------------------------%
%-----------------------------------------------------------------------%


function [l,w] = WellsCopper(M,FL)
%This short code gives a fault's length and width, according to its
%magnitude, in accordance with Wells & Coppersmith (1994) average equations
%(over all focal mechs).

logL = -2.44 + 0.59*M;

L = 10.^logL;

logW = -1.01 + 0.32*M;

W = 10.^logW;

l = FL*L;

w = FL*W;

end  %WellsCopper  



%-----------------------------------------------------------------------%
%--------------------------Function latlon2xy2--------------------------%
%-----------------------------------------------------------------------%

function xycoords = latlon2xy2(latlon,centerpoint)
%This function will take coordinates given in latitude and longitude and
%convert them to an x,y coordinate system with x pointing east, y pointing
%north.

%The zero point will be set at the point centerpoint entered in the format
%[lat lon]
lat = latlon(:,1);
lon = latlon(:,2);

%Number of latlon points being converted to xy coordinates
NumOfPoints = length(lat);


%Calculating the x and y coordinates of each point from its lat and lon
%coordinates
y = 111.1949*(lat - centerpoint(1));
x = 111.1949*(lon - centerpoint(2)).*cosd(lat);

%The xy coordinates
xycoords = [x y];

end  %latlon2xy2
   
