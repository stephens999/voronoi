#include "voronoi.hpp"
#include "readboundary.hpp"
#include <assert.h>
#include <sstream>

// VORONOI version 2.0.1

using namespace std; 

const string VORONOI_VERSION = "2.0.1";

// The following values are defaults and can change during the run
double YMIN = -999;  // note that in this code, X and Y are mixed up
double YMAX = -999; // X represents N-S and Y represent E-W. 
double XMIN = -999; // Except in the functions InElephantRange, InForest etc where they have the opposite correspondance
double XMAX = -999;  // [there is a switch in InRange to allow for this]
//double YMIN = -17;  // note that in this code, X and Y are mixed up
//double YMAX = 50; // X represents N-S and Y represent E-W. 
//double XMIN = -36; // Except in the functions InElephantRange, InForest etc where they have the opposite correspondance
//double XMAX = 31;  // [there is a switch in InRange to allow for this]
int GRIDSIZE = 67; 
int NITER = 100000;
int BURNIN = 5000;
int VLENGTH = 100; // length of vector of Voronoi points
int SAVANNAHONLY = -1;
int READGRID = 0;
int PRINTPROBS = 0;
bool SAMPLENAMEFILE = false;
bool KUHNERSIM = false;

int OFFSET = 0;

// Mary Kuhner addition:  permanent record of which grid squares are within (1)
// or outside (0) the species boundary.  Optimization plus improvement of
// meta-parameter handling.  9/2/2020

vector<vector<int> > MASK;

void error_and_exit(const string& msg) {
  cerr << msg << endl;
  exit(-1);
}

//------------------------------------------------------------------------------------------
// following functions cut & pasted from scat2.cpp(v2.0): jyamato 2020/09/03

void ReadInBoundary(ifstream & bfile, vector<double> & BoundaryX, vector<double> & BoundaryY)
{
   double x,y;
   do{
      bfile >> y;
      bfile >> x;
      cerr << y << "," << x << endl;
      x = PI*x/180;
      y = PI*y/180;
      BoundaryX.push_back(x);
      BoundaryY.push_back(y);
   } while(BoundaryX.size()<2 || x!=BoundaryX[0] || y!=BoundaryY[0]);

}

double isLeft(double x0, double y0, double x1, double y1,  double x2, double y2){
  return( (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0) );
}

//  winding number test for a point in a polygon
//  modelled on code from softsurfer.com, by Dan Sunday
//      Input:   x,y = a point,
//               BoundaryX and BoundaryY = points of a polygon with V[n]=V[0]
//      Return:  wn = the winding number (=0 only if (x,y) is outside polygon)
int IsInsideBoundary( double x, double y, const DoubleVec1d& BoundaryX, const DoubleVec1d& BoundaryY)
{
  if(BoundaryX.size() == 0) // if no boundary, just return 1
    return 1;

  int    wn = 0;    // the winding number counter

  // loop through all edges of the polygon
  for (int i=0; i < (BoundaryX.size()-1); i++) {   // edge from V[i] to V[i+1]
    if (BoundaryY[i] <= y) {         // start y <= P.y
      if (BoundaryY[i+1] > y)      // an upward crossing
        if (isLeft( BoundaryX[i], BoundaryY[i], BoundaryX[i+1], BoundaryY[i+1], x, y) > 0)  // P left of edge
          ++wn;            // have a valid up intersect
    }
    else {                       // start y > P.y (no test needed)
      if (BoundaryY[i+1] <= y)     // a downward crossing
        if (isLeft( BoundaryX[i], BoundaryY[i], BoundaryX[i+1] , BoundaryY[i+1], x, y) < 0)  // P right of edge
          --wn;            // have a valid down intersect
    }
  }
  return wn;
}


//------------------------------------------------------------------------------------------

bool InElephantRange(double x, double y){
  double x0,y0,x1,y1,a,b;

// check (crudely) for falling in sea or sahara
  if(y>0.29)
    return false; // Sahara region (N of Timbuktoo)
  
  x0 = 0.746; // exclude NE ethiopia
  y0 = -0.007;
  x1 = 0.648;
  y1 = 0.140;
  a = (y0-y1)/(x0-x1);
  b = y0 - a* x0;
  if(y > (a*x + b))
    return false;
  x0 = 0.572; 
  y0 = 0.182;
  a = (y0-y1)/(x0-x1);
  b = y0 - a* x0;
  if(y > (a*x + b))
    return false;
  x1 = 0.445;
  y1 = 0.135;
  a = (y0-y1)/(x0-x1);
  b = y0 - a* x0;
  if(y > (a*x + b)){
    x0 = 0.391;
    y0 = 0.253;
    a = (y0-y1)/(x0-x1);
    b = y0 - a* x0;
    if(y > (a*x + b))
      return false;
  }
  x0 = 0.391;
  y0 = 0.253;
  x1 = -0.022;
  y1 = 0.293;
  a = (y0-y1)/(x0-x1);
  b = y0 - a* x0;
  if(y > (a*x + b))
    return false;
  
  // exclude south of africa
  
  if(y < -0.368){
    x0 = 0.559;
    y0 = -0.501;
    x1 = 0.449;
    y1 = -0.440;
    a = (y0-y1)/(x0-x1);
    b = y0 - a* x0;
    if(y < (a*x + b))
      return false;
    x0 = 0.465;
    y0 = -0.368;
    a = (y0-y1)/(x0-x1);
    b = y0 - a* x0;
    if(y > (a*x + b))
      return false;
  }

  if(  (x < 0.157 ) 
       && ( y < (0.25*x+ PI*5.59/180) ) 
       && ( y < (-0.56*x+ PI*8.38/180) ) )
    return false; // bump in Lagos/Accra region
  
  if((x < 0.157) && (y<0.073) ) // SW coast
    return false;
  if((x < -0.2)) //W coast
    return false;
  if(y < -0.6) // S coast
    return false;
  if((x > 0.70) && (y<(-0.05)) )
    return false; // SE coast
  if((y < (-1.38 + 1.59 * x)))
    return false; // SEcoast, Durbin/Cidade de Nacala
  
  
  if((x > 0.70) && (y< -0.66 + (13.0/15)* x))
    return false; // Ecoast
 
  return true;
}

bool InForest(double x,double y){
  double x0,y0,x1,y1,a,b,a0,b0,a1,b1,a2,b2;
  bool forest = false;
  x0 = -0.003; // the part to the west of Mole
  y0 = 0.349;
  x1 = 0.038;
  y1 = 0.143;
  a = (y0-y1)/(x0-x1);
  b = y0 - a* x0;
  if( y < a*x + b){
    forest = true;
    //cout << "Is to the west of Mole" << endl;
  }
  
  // now check the rest of the forest areas
  x0 = 0.038;
  y0 = 0.143;
  x1 = 0.524;
  y1 = 0.087;
  a0 = (y0-y1)/(x0-x1);
  b0 = y0 - a0* x0;
  
  x0 = 0.489;
  y0 = -0.165;
  a1 = (y0-y1)/(x0-x1);
  b1 = y0 - a1* x0;
  
  x1 = 0.241;
  y1 = -0.157;
  a2 = (y0-y1)/(x0-x1);
  b2 = y0 - a2* x0;
  
  if( ( y < a0*x + b0) &&
      ( y > a1*x + b1) &&
      ( y > a2*x + b2))
    forest = true;

  if(y > 10*PI/180)
	  forest = false;
  
  return forest;  
}


double XCenterpoint(int j){
  return XMIN + (j+0.5) * (XMAX-XMIN)/GRIDSIZE;
}

double YCenterpoint(int k){
  return YMIN + (k+0.5) * (YMAX-YMIN)/GRIDSIZE;
}


bool InRange(double x,double y, bool READBOUNDARY, const DoubleVec1d& boundaryx, const DoubleVec1d& boundaryy,const Mapgrid& mymapgrid){
  // PI/180 because of degrees vs. radians.  x and y reversed because this program overall
  // and SCAT2, from which these three functions were taken, reverse the sense of x and y.
  // SORRY!
  if(READGRID) 
    return mymapgrid.in_range(x,y);
  if(READBOUNDARY)
     return (IsInsideBoundary(PI*y/180,PI*x/180,boundaryx,boundaryy));
  if(SAVANNAHONLY)
     return (InElephantRange(PI*y/180,PI*x/180) && !InForest(PI*y/180,PI*x/180));
  else
     return ((InElephantRange(PI*y/180,PI*x/180) && InForest(PI*y/180,PI*x/180)));
}


// test whether grid square (j,k) is in required range (Forest or Savannah)
// used only to create initial values in MASK
bool OldGridInRange(int j, int k, bool READBOUNDARY, const DoubleVec1d& boundaryx, const DoubleVec1d& boundaryy,const Mapgrid& mymapgrid){
  double x = XCenterpoint(j);
  double y = YCenterpoint(k);
  return InRange(x,y,READBOUNDARY,boundaryx,boundaryy,mymapgrid);
}

bool GridInRange(int j, int k) {
// relies on global MASK
  return (MASK[j][k] == 1);
}

// Mary Kuhner 9/2/2020
// Is there a legal location for every elephant in this Region configuration?
// used for both debugging and core logic
bool SumCountsLegal(const DoubleVec1d& sumcounts) {
  bool success = true;
  for (int i = 0; i < sumcounts.size(); ++i) {
    if (sumcounts[i] == 0) {
      success = false;
      cerr << "No sumcounts for elephant" << i << endl;
      break;
    }
  }
  return success;
}

bool SumCountsLegalWithCutoff(const DoubleVec1d& sumcounts, const DoubleVec1d& totalcounts, double cutoff) {
  bool success = true;
  for (int i = 0; i < sumcounts.size(); ++i) {
    if (sumcounts[i]/totalcounts[i] <= cutoff) {
      success = false;
      //cerr << "No cutoff-legal sumcounts for elephant" << i << endl;
      //cerr << "Sumcounts " << sumcounts[i] << " totalcounts " << totalcounts[i] << endl;
      break;
    }
  }
  return success;
}

// Debug function
bool ProgramStateValid(const DoubleVec1d& voronoiX, const DoubleVec1d& voronoiY, const
  IntVec1d& voronoiZ, const IntVec2d& bestpoint, const DoubleVec2d& dist, const IntVec2d& region, 
  const DoubleVec1d& sumcounts, const DoubleVec3d& counts, const DoubleVec2d& lcounts, int regionsize,
  const IntVec2d& mask, bool useboundaryfile, const DoubleVec1d& boundaryx, const DoubleVec1d& boundaryy,
  const Mapgrid& mymapgrid) {
  bool valid = true;
  // check dimensions of arrays
  if (voronoiX.size() != VLENGTH || voronoiY.size() != VLENGTH || voronoiZ.size() != VLENGTH) {
    cerr << "Array dimensions do not equal VLENGTH" << endl;
    valid = false;
  }

  // check legality of voronoiX and voronoiY coordinates
  for(int i=0; i < voronoiX.size(); i++) {
    if(voronoiX[i] < XMIN || voronoiX[i] > XMAX || voronoiY[i] < YMIN || voronoiY[i] > YMAX) {
    cerr << "Polygon center out of bounds at" << voronoiX[i] << ", " << voronoiY[i] << endl; 
    valid = false;
    }
  }

  // check that mask is correct
  for (int j=0; j<GRIDSIZE; j++) {
    for(int k=0; k<GRIDSIZE; k++) {
      if(OldGridInRange(j,k,useboundaryfile,boundaryx,boundaryy,mymapgrid) != mask[j][k]) {
        cerr << "MASK disagrees with range calculations for " << j << ", " << k << endl;
        valid = false;
      }
    }
  }

  // check that bestpoint and dist accurate to current voronoi X,Y
  IntVec2d newbestpoint(bestpoint);
  DoubleVec2d newdist(dist);
  ComputeBestPointAndDist(voronoiX,voronoiY,newbestpoint,newdist);
  if (newbestpoint != bestpoint) {
    cerr << "Bestpoint not consistent when recalculated" << endl;
    valid = false;
  }
  if (newdist != dist) {
    cerr << "Dist not consistent when recalculated" << endl;
    valid = false;
  }

  // check relationship between voronoiZ and region
  IntVec2d newregion(region);
  for (int j=0; j<GRIDSIZE; j++) {
    for(int k=0; k<GRIDSIZE; k++) {
      if(GridInRange(j,k) && voronoiZ[bestpoint[j][k]]) {
        newregion[j][k] == 1;
      } else {
        newregion[j][k] == 0;
      }
    }
  }
  if (newregion != region) {
    cerr << "VoronoiZ and REGION inconsistent" << endl;
    valid = false;
  }

// check correctness of regionsize
  int checkregionsize = 0;
  for(int j=0; j<GRIDSIZE; j++) {
    for(int k=0; k<GRIDSIZE; k++) {
      if (GridInRange(j,k)) {
        checkregionsize += region[j][k];
      }
    }
  }
  if(checkregionsize != regionsize) {
    cerr << "REGIONSIZE not consistent with contents of REGION" << endl;
    cerr << "Stored = " << regionsize << " calculated = " << checkregionsize << endl;
    valid = false;
  }

// check correctness of sumcounts
  int NIND = sumcounts.size();
  for(int i = 0; i < NIND; i++) {
    float check_sumcounts = 0.0;
    for(int j=0; j<GRIDSIZE; j++) {
      for(int k=0; k<GRIDSIZE; k++) {
        if(region[j][k]==1 && GridInRange(j,k)){
          check_sumcounts += counts[i][j][k];
        }
      }
    }
    if (check_sumcounts != sumcounts[i]) {
      cerr << "sumcounts inconsistent for sample " << i << " with computed " << check_sumcounts << " and stored " << sumcounts[i] << endl;
      valid = false;
    }

// check correctness of lcounts
    DoubleVec1d check_lcounts(VLENGTH,0.0);
    for(int j=0; j<GRIDSIZE; j++) {
      for(int k=0; k<GRIDSIZE; k++) {
        int mypoly = bestpoint[j][k];
        if(GridInRange(j,k)) {
          check_lcounts[mypoly] += counts[i][j][k];
        }
      }
    }
    for(int poly = 0; poly < VLENGTH; poly++) {
      if (check_lcounts[poly] != lcounts[i][poly]) {
        cerr << "lcounts inconsistent for sample " << i << " polygon " << poly << " with computed " << check_lcounts[poly] << " and stored " << lcounts[i][poly] << endl;
        valid = false;
      }
    }
  }

  return valid;
}


double dist_between(double x0, double y0, double x1, double y1) {
  double answer;
  answer = pow((x0-x1),2.0) + pow((y0-y1),2.0);
  answer = sqrt(answer);
  return answer;
}

// Function by Mary 8/27/2020
// take an x, y pair from SCAT2 that falls in an out of bounds
// VORONOI grid square, and move it to the nearest legal grid
// square; if no legal square is adjacent, abort the program.
// Return grid coordinates of legal square.

std::pair<int,int> make_legal(double x, double y) {

  // gridi and gridj are the grid coordinates (not lat/long!) of the grid
  // square where these SCAT2 coordinates would naturally fall
  int gridi = (int) trunc(GRIDSIZE * (x-XMIN)/(XMAX-XMIN));
  int gridj = (int) trunc(GRIDSIZE * (y-YMIN)/(YMAX-YMIN));

// find grid squares
  vector<pair<double,double> > candidates;
  for (int myi = gridi - 1; myi <= gridi + 1; ++myi) {
    if (myi >= 0 && myi < GRIDSIZE) {
      for (int myj = gridj - 1; myj <= gridj + 1; ++myj) {
        if (myj >= 0 && myj < GRIDSIZE) {
          // don't test the middle square as it's already known to be bad
          if (!(myi == gridi && myj == gridj)) {
            candidates.push_back(make_pair(myi,myj));
          }
        }
      }
    }
  }
 
  double min_dist = DBL_MAX;
  int mini = -1;
  int minj = -1;
 
  // test all candidates and find the closest
  for (int k = 0; k < candidates.size(); ++k) {
    int i = candidates[k].first;
    int j = candidates[k].second;
    if (!GridInRange(i,j)) continue;
    // measure distance between x,y and center of voronoi square
    double testx = XCenterpoint(i);
    double testy = YCenterpoint(j);
    double dist = dist_between(x,y,testx,testy);
    if (dist<min_dist) {
      min_dist = dist;
      mini = i;
      minj = j;
    }
  }

  if (mini == -1 || minj == -1) {
    cerr << "Location of elephant at " << x << "," << y << "far out of bounds" << endl;
    exit(-1);
  }
  return make_pair(mini,minj);
}

// number of (valid savannah/forest  counts in each voronoi cell
void ComputeLCount(const DoubleVec3d& counts, const IntVec2d& bestpoint, DoubleVec2d& lcounts){
  for(int i=0;i<counts.size();i++){
    for(int l = 0; l<VLENGTH; l++){
      lcounts[i][l] = 0;
    }
    for(int j=0;j<GRIDSIZE; j++){
      for(int k=0; k<GRIDSIZE; k++){
	if(GridInRange(j,k)){
	  lcounts[i][bestpoint[j][k]] += counts[i][j][k];
	} else {
          assert(counts[i][j][k] == 0);
        }
      }
    }
  }
}

void ComputeTotalCounts(const DoubleVec2d& lcounts, DoubleVec1d& totalcounts){
  for(int i=0; i<lcounts.size(); i++){
    totalcounts[i] = 0.0;
    for(int l=0; l<VLENGTH;l++){
      totalcounts[i] += lcounts[i][l];
    }
  }
}


void ComputeBestPointAndDist(const DoubleVec1d& vx, const DoubleVec1d& vy, IntVec2d& bestpoint, DoubleVec2d& dist){
  for(int j=0; j<GRIDSIZE; j++){
    for(int k=0; k<GRIDSIZE; k++){
      double x = XCenterpoint(j);
      double y = YCenterpoint(k);
      dist[j][k] = (x-vx[0])*(x-vx[0]) + (y-vy[0])*(y-vy[0]);
      bestpoint[j][k] = 0;
      for (int l=1; l<VLENGTH; l++){
	double newdist = (x-vx[l])*(x-vx[l]) + (y-vy[l])*(y-vy[l]);
	if(newdist < dist[j][k]){
	  bestpoint[j][k] = l;
	  dist[j][k] = newdist;
	}
      }
    }
  }    
}

//update the best point after index l of vx and vy has been changed
void UpdateBestPointAndDist(const DoubleVec1d& vx, const DoubleVec1d& vy, int l, IntVec2d& bestpoint, DoubleVec2d& dist) {
  for(int j=0; j<GRIDSIZE; j++){
    for(int k=0; k<GRIDSIZE; k++){
      double x = XCenterpoint(j);
      double y = YCenterpoint(k);
      if(bestpoint[j][k] == l){ // in this case need to do entire search of points
	dist[j][k] = (x-vx[0])*(x-vx[0]) + (y-vy[0])*(y-vy[0]);
	bestpoint[j][k] = 0;
	for (int l=1; l<VLENGTH; l++){
	  double newdist = (x-vx[l])*(x-vx[l]) + (y-vy[l])*(y-vy[l]);
	  if(newdist < dist[j][k]){
	    bestpoint[j][k] = l;
	    dist[j][k] = newdist;
	  }
	}
	
      } else { // just need to check if l is now closer
	double newdist = (x-vx[l])*(x-vx[l]) + (y-vy[l])*(y-vy[l]);
	if(newdist < dist[j][k]){
	  bestpoint[j][k] = l;
	  dist[j][k] = newdist;
	}
      }
    }
  }    
}


void ComputeSumCounts2(const DoubleVec2d& LCounts, const IntVec1d& vz, DoubleVec1d& SumCounts) {
  for(int i = 0;i<LCounts.size(); i++){
    SumCounts[i]=0;
    for(int l=0; l<VLENGTH; l++){
      SumCounts[i] += LCounts[i][l]*vz[l];
    }
  }
}
  


// compute grid of 1s and 0s corresponding to a voronoi tesselation
void ComputeRegion(const IntVec1d& vz, IntVec2d& Region, const IntVec2d& BestPoint) {
  for(int j=0; j<GRIDSIZE; j++){
    for(int k=0; k<GRIDSIZE; k++){
      Region[j][k] = vz[BestPoint[j][k]];
    }
  }
}

void ComputeCellSizeInRange(const IntVec2d& BestPoint, IntVec1d& CellSize) {
  for(int l =0;l<VLENGTH; l++)
    CellSize[l] = 0;
    
  for(int j=0; j<GRIDSIZE; j++){
    for(int k=0; k<GRIDSIZE; k++){
      if(GridInRange(j,k))
	CellSize[BestPoint[j][k]]++;
    }
  }
}

void OutputRegion(const IntVec2d& REGION, ofstream& output, bool OnlyRange) {
 
  for(int j = 0; j<GRIDSIZE; j++){
    for(int k=0; k<GRIDSIZE; k++){
      int inrange = GridInRange(j,k) || !OnlyRange;
      output << inrange*REGION[j][k] << " ";
    }
    output << endl;
  }
}


void ComputeRegionSizeInSav(const IntVec2d& REGION, int& REGIONSIZE) {
  REGIONSIZE = 0;
  for(int j = 0; j<(GRIDSIZE); j++){
    for(int k = 0; k<(GRIDSIZE); k++){
      if(GridInRange(j,k)){
	REGIONSIZE += REGION[j][k];
      }
    }
  }
}

void UpdateZs(double Vprob, const IntVec1d& CellSize, int& REGIONSIZE, DoubleVec1d& SUMCOUNTS, const DoubleVec2d& LCOUNTS, const DoubleVec1d& TOTALCOUNTS, IntVec2d& REGION, IntVec1d& VoronoiZ, const IntVec2d& BESTPOINT, double cutoff) {
  // update the Voronoi Zs
  int newregionsize;
  int NIND = LCOUNTS.size();
  double newloglik, currentloglik;
  double newsumcounts;

  for(int l=0; l<VLENGTH; l++){

    //compute relative prob of 0 and 1 v point
    newloglik = 0; currentloglik = 0;
    if(VoronoiZ[l] == 0){
      newregionsize = REGIONSIZE + CellSize[l];
      newloglik += log(Vprob);
      currentloglik += log(1-Vprob);
    }
    else{
      newregionsize = REGIONSIZE - CellSize[l]; 
      newloglik += log(1-Vprob);
      currentloglik += log(Vprob);
    }
    bool force_reject = false;
    if(newregionsize>0){
      for(int i = 0; i<NIND; i++){
	if(VoronoiZ[l] == 0)
	  newsumcounts = SUMCOUNTS[i] + LCOUNTS[i][l];
	else 
	  newsumcounts = SUMCOUNTS[i] - LCOUNTS[i][l];
        // Changed by Mary Kuhner 11/17/2020 to reject any solution which does not
        // sufficiently (i.e. "cutoff") explain every individual.
        if (newsumcounts/TOTALCOUNTS[i] < cutoff) {
          force_reject = true;
        }
	newloglik += log(newsumcounts/(TOTALCOUNTS[i] * newregionsize));
	currentloglik += log(SUMCOUNTS[i]/(TOTALCOUNTS[i] * REGIONSIZE));
      }
      
      if((ranf() < exp(newloglik - currentloglik) && !force_reject)){
        // acceptance logic
        VoronoiZ[l] = 1-VoronoiZ[l];
        ComputeRegion(VoronoiZ,REGION,BESTPOINT);
        REGIONSIZE = newregionsize;
	if(VoronoiZ[l] == 1){
	  for(int i = 0; i < NIND; i++){
	    SUMCOUNTS[i] += LCOUNTS[i][l];
	  }
	} else {
	  for(int i = 0; i < NIND; i++){
	    SUMCOUNTS[i] -= LCOUNTS[i][l];
	  }
        }
      }
    }
  }
}

// Currently, XMIN, XMAX, YMIN, YMAX, and GRIDSIZE are all program global constants
// COUNTS, despite the capitalization, is a local variable of the primary routine
void ReadScatInfile(ifstream& locatefile, const string& infile, int skip, int nmcmc, int ind, DoubleVec3d& COUNTS) {
  double x,y,dummy;
  // skipping SCAT burnin
  for(int index = 0; index<skip; index++){
    locatefile >> x;
    locatefile >> y;
    locatefile >> dummy;
  }
  
  for(int index = 0; index<nmcmc; index++){
    locatefile >> x;
    locatefile >> y;
    locatefile >> dummy;

    if(x>XMAX || x<XMIN || y>YMAX || y<YMIN){
      cerr << "Error: x or y out of bounds" << endl;
      cerr << "file = " << infile << endl;
      cerr << "x= " << x << endl;
      cerr << "y= " << y << endl;
      exit(1);
    }

    int j = (int) trunc(GRIDSIZE * (x-XMIN)/(XMAX-XMIN));
    int k = (int) trunc(GRIDSIZE * (y-YMIN)/(YMAX-YMIN));

    // Adjust mildly illegal input to the nearest legal grid square
    if(!GridInRange(j,k)){
      pair<int,int> newvals = make_legal(x,y);
      j = newvals.first;
      k = newvals.second;
    }
  	
    assert(GridInRange(j,k));
    COUNTS[ind][j][k] += 1;
  }
}

int main ( int argc, char** argv)
{
  cout << "VORONOI version " << VORONOI_VERSION << " ";
  cout << "Compiled " << __DATE__ << " " << __TIME__ << endl;
  bool READBOUNDARY = false;
  string bfilename;
  int SEED = 0;
  float REJECT_CUTOFF = 0.0;
  int nmcmc = 100;
  int skip = 100;
  map<string, string> filenames; 
  while( ( argc > 1 ) && ( argv[1][0] == '-' ) ) {
    switch(argv[1][1]) {
      
   case 'B':
	READBOUNDARY = true;
        bfilename = argv[1]+2;
        break;

   case 'C':
        ++argv;
        --argc;
        REJECT_CUTOFF = atof(&argv[1][0]);
        if (REJECT_CUTOFF < 0.0 || REJECT_CUTOFF >= 1.0) {
          cerr << "Invalid reject cutoff" << endl;
          exit(-1);
        }
        cout << "Reject cutoff = " << REJECT_CUTOFF << endl;
        break;

   case 'D':  // forest
	SAVANNAHONLY = 0;
	break;

   case 'd':  // savannah
	SAVANNAHONLY = 1;
	break;
     
    case 'g': // read in grid file, replacement for a boundary file
      ++argv; --argc;
      filenames["gridfile"] = argv[1];
      READGRID = 1;
      break;

   case 'k': // assume input data in kuhnersim directory tree
             // WARNING: ignores any input sample settings in the path file
             //  read by 'N' setting
      KUHNERSIM = true;
      SAMPLENAMEFILE = true;
      filenames["pathfile"] = "non_existant_pathfile_invoked";
      break;

   case 'm':  // number of steps in the VORONOI mcmc, including burnin
      ++argv; 
      --argc;
      NITER = atoi(&argv[1][0]);
      break; 

   case 'n':  // specify number of SCAT points per sample
      ++argv;
      --argc;
      nmcmc = atoi(&argv[1][0]);
      break;
     
   case 'N': // use namefile for input rather than standard "voronoiin.txt" style input
      SAMPLENAMEFILE = true;
      ++argv; --argc;
      filenames["pathfile"] = argv[1];
      break;

   case 'S': // seed
      ++argv;
      --argc;
      SEED = atoi(&argv[1][0]);
       cout << "Seed = " << SEED << endl;
    break;

    case 'v':  // print additional output format
      PRINTPROBS = 1;
      break;

    default: 
      cerr << "Error: option " << argv[1] << "unrecognized" << endl;
      return 1;
      
    }
    ++argv;
    --argc;
  }

  if(argc<3){
    cerr << "Usage is ./VORONOI samplefile outfile " << endl;
    exit(1);
  }

  bool BUILTIN = (SAVANNAHONLY != -1);

  if(!BUILTIN && !READBOUNDARY && !READGRID) {
    cerr << "Need to specify how boundaries are determined" << endl;
    exit(1);
  }

  if((BUILTIN && READBOUNDARY) || (BUILTIN && READGRID) || (READBOUNDARY && READGRID)) {
    cerr << "Two different boundary specifications given; please choose one" << endl;
    exit(1);
  }

  DoubleVec1d boundaryx,boundaryy;
  if(READBOUNDARY) {
    cout << "Reading in boundary data from boundary file " << bfilename.c_str() << endl;
    ifstream bfile(bfilename.c_str());
    if(!bfile.is_open()) {
      error_and_exit("Failed to open boundary file " + bfilename);
    }
    ReadInBoundary(bfile,boundaryx,boundaryy);
  }

  srandom(SEED);

  filenames["samples"] = argv[1];
  filenames["output"]  = argv[2];
   
  // Mary and Jon 10/13/2020
  Mapgrid mymapgrid;
  if(READGRID){
    cerr << "Reading in gridfile data from file " << filenames["gridfile"] << endl;
    ifstream gridfile (filenames["gridfile"].c_str());
    if (!gridfile.is_open()) {
      error_and_exit("Could not open grid file " + filenames["gridfile"]);
    }
    mymapgrid.Initialize(gridfile,MARGIN);
    GRIDSIZE = mymapgrid.get_gridsize();
    std::pair<int,int> latbounds = mymapgrid.get_latbounds();
    std::pair<int,int> longbounds = mymapgrid.get_longbounds();
    XMIN = latbounds.first;
    XMAX = latbounds.second;
    YMIN = longbounds.first;
    YMAX = longbounds.second; 
  } else {
    // WARNING warning here are the original program defaults, good for africa....
    // note again that X represent N-S and Y represent E-W
    //   except in the functions InElephantRange, InForest etc where the correspondance
    //   is switched.  The master function InRange handles this.
    YMIN = -17;
    YMAX = 50;
    XMIN = -36;
    XMAX = 31;
  }

  // Mary Kuhner 6/9/2020
  // Make a mask array which shows those grid cells that are not within the
  // species boundary.  This has three purposes:  (1)  Speed up calls to
  // GridInRange; (2)  Allow non-use of masked out grid cells in calculation
  // of the Vprob acceptance; (3) Avoid toggling Z array for masked out cells.
  // This is expected both to make the code faster and to remove a strong bias
  // in estimation of Vprob.

  int jcoord;
  int kcoord;
  vector<vector<int> > tempmask(GRIDSIZE,vector<int>(GRIDSIZE,0));
  MASK = tempmask;
  cout << endl;
  for (jcoord = 0; jcoord < GRIDSIZE; ++jcoord) {
    for (kcoord = 0; kcoord < GRIDSIZE ; ++kcoord) {
      if (OldGridInRange(jcoord,kcoord,READBOUNDARY,boundaryx,boundaryy,mymapgrid)){
        MASK[jcoord][kcoord] = 1;
      }
    }
  }

  string mapinfoname = filenames["output"] + "_mapinfo";
  ofstream mapinfofile (mapinfoname.c_str());
  string regionfilename = filenames["output"] + "_regions";
  ofstream regionfile (regionfilename.c_str());
  string indprobfilename = filenames["output"] + "_indprobs";
  ofstream indprobfile (indprobfilename.c_str());
  string samplefilename = filenames["samples"];
  ifstream samplefile(samplefilename.c_str());
  string tracefilename = filenames["output"] + "_trace";
  ofstream tracefile(tracefilename.c_str());
  
  // optional printprobs file
  string printprobsname = filenames["output"] + "_printprobs";
  ofstream printprobsfile; 
  if(PRINTPROBS) 
    printprobsfile.open(printprobsname.c_str(),ios::trunc);

  // print out info about map
  if (READGRID) {
    mymapgrid.WriteMapInfo(mapinfofile);
  } else {
    mapinfofile << "Using default map information for ";
    if (SAVANNAHONLY) {
      mapinfofile << "Savannah elephants " << endl;
    }
    if (READBOUNDARY) {
      mapinfofile << "Custom boundary file " << endl;
    }
    if (!SAVANNAHONLY && !READBOUNDARY) {
      mapinfofile << "Forest elephants " << endl;
    }
    mapinfofile << "Lower left corner at -36, -17" << endl;
    mapinfofile << "Upper right corner at 31, 50" << endl;
    mapinfofile << "North-South extent of grid:  67 squares" << endl;
    mapinfofile << "East-West extent of grid:  67 squares" << endl;
   }


  vector<double> VoronoiX(VLENGTH,0);
  vector<double> VoronoiY(VLENGTH,0);
  vector<int> VoronoiZ(VLENGTH,0);
  vector<int> CellSize(VLENGTH,0); // number of in-range grid squares in each polygon


  int NIND; // = lastsample - firstsample + 1;
 
  // read in samples to be located from samplefile
  vector<string> sampleids;
  vector<string> samplepaths;
  vector<int> samplevec;
  if (SAMPLENAMEFILE) {
    string sid;
    while (samplefile >> sid) {
      sampleids.push_back(sid);
    }
    NIND = int(sampleids.size());

    if (KUHNERSIM) {
      for(int i = 0; i < 9; i++) {
        samplepaths.push_back(to_string(i+1)+"/outputs/");
      }
    } else {
      ifstream pathfile (filenames["pathfile"].c_str());
      if (!pathfile.is_open()) {
        error_and_exit("Could not open path file " + filenames["pathfile"]);
      }
      string path_to_samples;
      while(pathfile >> path_to_samples) {
        if (path_to_samples.empty())
          continue;
        if (path_to_samples.back() != '/')
          path_to_samples += '/';
        samplepaths.push_back(path_to_samples);
      }
      pathfile.close();
      cout << "Samples will be pulled from " << to_string(samplepaths.size()) << " directories";
      cout << " found in " << filenames["pathfile"] << endl;
    }
  } else {
    samplefile >> NIND;
    for(int s=0;s<NIND;s++) {
      int foo;
      samplefile >> foo;
      samplevec.push_back(foo);
    }
  }

 //counts of number of times each ind is sampled in each grid square
  vector<vector<vector<double> > > COUNTS(NIND,vector<vector<double> >(GRIDSIZE,vector<double>(GRIDSIZE,0)));
  //posterior prob of each ind coming from each grid square
  vector<vector<double> > LCOUNTS(NIND, vector<double>(VLENGTH,0)); // number of times each individual is in each voronoi cell
  vector<double> TOTALCOUNTS(NIND,0);

  
  vector<vector<vector<double> > > INDPROBS(NIND,vector<vector<double> >(GRIDSIZE,vector<double>(GRIDSIZE,0)));
  vector<double> SUMCOUNTS(NIND,0); // sum of number of times each individual is in REGION


  //0/1 indicator of whether region is in or out
  vector<vector<int> > REGION(GRIDSIZE,vector<int>(GRIDSIZE,0));
  vector<vector<int> > BESTPOINT(GRIDSIZE,vector<int>(GRIDSIZE,0)); //which voronoi point each grid closest to
  vector<vector<double> >  DIST(GRIDSIZE,vector<double>(GRIDSIZE,0)); //distance of each grid point to bestpoint
  vector<vector<double> > PROB(GRIDSIZE,vector<double>(GRIDSIZE,0));
  int REGIONSIZE = 0;

  // Mary 9/2/2020:  moved region initialization to after data reading,
  // so it can be retried if the proposed region doesn't contain all
  // elephants

  if (SAMPLENAMEFILE) {
    for (int p = 0; p < samplepaths.size(); ++p) {
      for (int s = 0; s < sampleids.size(); ++s) {
        string filepath = samplepaths[p] + sampleids[s];
        ifstream locatefile(filepath.c_str());
        if(!locatefile.is_open()) {
          error_and_exit("Failed to open input file " + filepath);
        }
    
        ReadScatInfile(locatefile,filepath,skip,nmcmc,s,COUNTS);
      }
    }
  } else {
    for(char TAG = 'r'; TAG<='z'; TAG++){
       //cout << "Reading " << TAG << endl;
       for(int s = 0; s < NIND; s++){
         int sample = samplevec[s];
         // new code by Mary Kuhner 8/21/2020, fixing bug that
         // disallowed elephant numbers > 999
         string infile = std::to_string(sample);
         if (infile.length() == 1) {
           infile = string("00") + infile;
         }
         if (infile.length() == 2) {
           infile = string("0") + infile;
         }
         infile += TAG;
         //cout << "Reading from file " << infile << endl;

         ifstream locatefile(infile.c_str());
         if(!locatefile.is_open()) {
           error_and_exit("Failed to open input file " + infile);
         }
    
         ReadScatInfile(locatefile,infile,skip,nmcmc,s,COUNTS);
      }
    }
  }

  cerr << "Finished data initialization" << endl;
  double Vprob = 0.5; // prob of each voronoi point being a 1

  //initialise REGION
  // propose starting tesselations until a legal one is found
  // new code Mary Kuhner 9/2/2020
  bool force_solution = false;
  IntVec2d saveregion;
  for(int tries = 0; tries < 500; ++tries) {
    // cerr << "Trying tesselation " << tries << endl;
    for(int l = 0; l<VLENGTH; l++){
      VoronoiX[l] = XMIN + ranf() * (XMAX-XMIN);
      VoronoiY[l] = YMIN + ranf() * (YMAX-YMIN);
    }
    ComputeBestPointAndDist(VoronoiX,VoronoiY,BESTPOINT,DIST);
    ComputeCellSizeInRange(BESTPOINT,CellSize);
    // Initialize polygons outside species range to non-alive
    for(int l = 0; l<VLENGTH; l++) {
      if (force_solution) {
        VoronoiZ[l] = true;   // forcing everything else to be true
        continue;
      }
      VoronoiZ[l] = (ranf()<Vprob);  // picking at random
    }
    ComputeRegion(VoronoiZ,REGION,BESTPOINT);  
    saveregion = REGION;
    ComputeRegionSizeInSav(REGION,REGIONSIZE);
    // re-initialize SUMCOUNTS for this try!
    for (int i = 0; i < NIND; ++i) {
      for (int j = 0; j < GRIDSIZE; ++j) {
        for (int k = 0; k < GRIDSIZE; ++k) {
	  SUMCOUNTS[i] = 0.0;
        }
      }
    }
    // now compute it
    for (int i = 0; i < NIND; ++i) {
      for (int j = 0; j < GRIDSIZE; ++j) {
        for (int k = 0; k < GRIDSIZE; ++k) {
	  if(REGION[j][k]==1 && GridInRange(j,k)){
	    SUMCOUNTS[i] += COUNTS[i][j][k];
          }
        }
      }
    }
    // Test if this tesselation is legal
    ComputeLCount(COUNTS,BESTPOINT,LCOUNTS);
    ComputeTotalCounts(LCOUNTS,TOTALCOUNTS);

    if (SumCountsLegalWithCutoff(SUMCOUNTS,TOTALCOUNTS,REJECT_CUTOFF)) 
      break;  // found a legal tesselation

    // nope, wasn't a legal tesselation; can we try again?
    if (tries % 10 == 1) {
      cerr << "Retry" << tries << endl;
    }
    if (tries == 498) {  // too many tries
      force_solution = true;
      cerr << "Forcing an all-alive Voronoi tesselation as starting point" << endl;
    }
    continue;   // try another tesselation
  }
  if (saveregion != REGION) {
    cout << "Region has changed!" << endl;
  }
  assert(SumCountsLegalWithCutoff(SUMCOUNTS,TOTALCOUNTS,REJECT_CUTOFF));
  if(!ProgramStateValid(VoronoiX, VoronoiY, VoronoiZ, BESTPOINT, DIST,
    REGION, SUMCOUNTS, COUNTS, LCOUNTS, REGIONSIZE, MASK, READBOUNDARY,
    boundaryx, boundaryy, mymapgrid)) {
    cerr << "Failed right after inital tesselation" << endl;
    exit(-1);
  }

  //ComputeLCount(COUNTS,BESTPOINT,LCOUNTS);
  //ComputeTotalCounts(LCOUNTS,TOTALCOUNTS);

  int ACCEPT=0;
  
  double newloglik, currentloglik;

  int output_now = NITER / 5;
  int output_interval = output_now;
  cout << "Starting MCMC iterations" << flush;

  for(int iter = 0; iter<NITER; iter++){
    assert(SumCountsLegal(SUMCOUNTS));
    if(!ProgramStateValid(VoronoiX, VoronoiY, VoronoiZ, BESTPOINT, DIST,
      REGION, SUMCOUNTS, COUNTS, LCOUNTS, REGIONSIZE, MASK, READBOUNDARY,
      boundaryx, boundaryy, mymapgrid)) {
      cerr << "Failed at top of iter loop" << endl;
      exit(-1);
    }

    if (iter == output_now) {
      cout << "." << flush;
      output_now += output_interval;
    }

    UpdateZs(Vprob,CellSize,REGIONSIZE,SUMCOUNTS,LCOUNTS,TOTALCOUNTS,REGION,VoronoiZ,BESTPOINT,REJECT_CUTOFF);

    if(!ProgramStateValid(VoronoiX, VoronoiY, VoronoiZ, BESTPOINT, DIST,
      REGION, SUMCOUNTS, COUNTS, LCOUNTS, REGIONSIZE, MASK, READBOUNDARY,
      boundaryx, boundaryy, mymapgrid)) {
      cerr << "Failed after UpdateZs" << endl;
      exit(-1);
    }
    
    assert(SumCountsLegal(SUMCOUNTS));
    double newVprob = Vprob + rnorm(0,0.02);
    if(newVprob>0 && newVprob<1){
      double cells_on = 0;
      double total_cells = 0;
      for(int l=0; l<VLENGTH; l++) {
         // do not include polygons outside species range in the denominator
         if (CellSize[l] > 0) {
           total_cells++;
           cells_on += VoronoiZ[l];
         }
      }
      double cells_off = total_cells - cells_on;
      newloglik = cells_on * log(newVprob) + cells_off * log(1-newVprob);
      currentloglik = cells_on * log(Vprob) + cells_off * log(1-Vprob);
      if(ranf()< exp(newloglik - currentloglik))
	Vprob = newVprob;
      //cout << "Vprob = " << Vprob << endl;
      tracefile << iter << " " << Vprob << endl;
    }

    assert(SumCountsLegal(SUMCOUNTS));
    vector<double> oldVoronoiX(VoronoiX);
    vector<double> oldVoronoiY(VoronoiY);
    vector<vector<double> > oldLCOUNTS(LCOUNTS);
    vector<double> oldSUMCOUNTS(SUMCOUNTS); // sum of number of times each individual is in REGION
    vector<vector<int> > oldREGION(REGION);
    vector<vector<int> > oldBESTPOINT(BESTPOINT);
    vector<vector<double> >  oldDIST(DIST);
    int oldREGIONSIZE = REGIONSIZE;
    
    for(int l=0; l<VLENGTH; l++){
      VoronoiX[l] = VoronoiX[l] + rnorm(0,5);
      VoronoiY[l] = VoronoiY[l] + rnorm(0,5);
      if(VoronoiX[l] < XMIN)
	VoronoiX[l] = XMIN + (XMIN-VoronoiX[l]);
      if(VoronoiX[l] > XMAX)
	VoronoiX[l] = XMAX - (VoronoiX[l]-XMAX);
      if(VoronoiY[l] < YMIN)
	VoronoiY[l] = YMIN + (YMIN-VoronoiY[l]);
      if(VoronoiY[l] > YMAX)
	VoronoiY[l] = YMAX - (VoronoiY[l]-YMAX);	  
    }
   
    currentloglik = newloglik = 0;
    for(int i = 0; i<NIND; i++){
      currentloglik += log(SUMCOUNTS[i]/(TOTALCOUNTS[i] * REGIONSIZE));
    }
    ComputeBestPointAndDist(VoronoiX,VoronoiY,BESTPOINT,DIST);
    ComputeRegion(VoronoiZ,REGION,BESTPOINT);  
    ComputeRegionSizeInSav(REGION,REGIONSIZE);
    ComputeLCount(COUNTS,BESTPOINT,LCOUNTS);
    ComputeSumCounts2(LCOUNTS,VoronoiZ,SUMCOUNTS);
    
    bool reject=false;
    for(int i = 0; i<NIND; i++){
      if(SUMCOUNTS[i] / TOTALCOUNTS[i] <= REJECT_CUTOFF)
	reject = true;
      newloglik += log(SUMCOUNTS[i]/(TOTALCOUNTS[i] * REGIONSIZE));
    }
    
    if(reject || ranf() > exp(newloglik - currentloglik)){ // reject move
      VoronoiX = oldVoronoiX;
      VoronoiY = oldVoronoiY;
      LCOUNTS = oldLCOUNTS;
      SUMCOUNTS = oldSUMCOUNTS;
      REGION = oldREGION;
      BESTPOINT = oldBESTPOINT;
      DIST = oldDIST;
      REGIONSIZE = oldREGIONSIZE;
    } else{
      ComputeCellSizeInRange(BESTPOINT,CellSize);
      ACCEPT +=1;
    }
    assert(SumCountsLegal(SUMCOUNTS));
    if(!ProgramStateValid(VoronoiX, VoronoiY, VoronoiZ, BESTPOINT, DIST,
      REGION, SUMCOUNTS, COUNTS, LCOUNTS, REGIONSIZE, MASK, READBOUNDARY,
      boundaryx, boundaryy, mymapgrid)) {
      cerr << "Failed after main rearrangement logic" << endl;
      exit(-1);
    }

    if(iter>BURNIN){
      for(int j = 0; j<GRIDSIZE; j++){
 	for(int k=0; k<GRIDSIZE; k++){
	  int inrange = GridInRange(j,k);
	  regionfile << (inrange*REGION[j][k]) << " ";
	  PROB[j][k] += inrange * REGION[j][k];
	  for(int i = 0; i<NIND; i++){
	    if(SUMCOUNTS[i]==0){
	      cerr<< "Error: Sumcounts reached 0";
	      exit(1);
	    }
	    INDPROBS[i][j][k] += inrange*REGION[j][k]*COUNTS[i][j][k]/SUMCOUNTS[i];	  
	  }
	}
	regionfile << endl;
      }
    }
  }
  cout << "." << endl;
  
  string outputfilename = filenames["output"];
  ofstream output (outputfilename.c_str());
  
  for(int j = 0; j<GRIDSIZE; j++){
    for(int k=0; k<GRIDSIZE; k++){
      output << PROB[j][k]/(NITER - BURNIN) << " ";
    }
    output << endl;
  }

  for(int i =0; i<NIND; i++){
    // the implicit map between individual number and string identitifer is set up
    // when ReadScatInfile() is called; via the passed argument "s" in the calling routine.
    if(SAMPLENAMEFILE) {
      printprobsfile << "#" << sampleids[i] << endl;
    } else {
      printprobsfile << "#" << i << endl;
    }
    for(int j = 0; j<GRIDSIZE; j++){
      for(int k=0; k<GRIDSIZE; k++){
        assert(!(INDPROBS[i][j][k] != 0.0 && !GridInRange(j,k)));
        float probval = INDPROBS[i][j][k]/(NITER - BURNIN);
	indprobfile << probval << " ";
        if (PRINTPROBS) {
          float lati = XCenterpoint(j);
          float longi = YCenterpoint(k);
          printprobsfile << lati << " " << longi << " " << probval << endl;
        }
      }
      indprobfile << endl;
    }
  }

  
cout << endl << "Program done" << endl;
}
