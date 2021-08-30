#include "utility.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <numeric>
#include <cfloat>

const double PI = 3.141592;
const long MARGIN = 7;  // margin around grid

// VORONOI version 2.0.1

class Mapgrid;

typedef std::vector<double> DoubleVec1d;
typedef std::vector<std::vector<double> > DoubleVec2d;
typedef std::vector<std::vector<std::vector<double> > > DoubleVec3d;
typedef std::vector<int> IntVec1d;
typedef std::vector<std::vector<int> > IntVec2d;
typedef std::vector<std::vector<std::vector<int> > > IntVec3d;

bool InElephantRange(double x, double y);
bool InForest(double x, double y);
double XCenterpoint(int j);
double YCenterpoint(int k);
bool InRange(double x, double y, bool useboundaryfile, const DoubleVec1d& boundaryx,
  const DoubleVec1d& boundaryy, const Mapgrid& mymapgrid);
bool GridInRange(int j, int k);
bool OldGridInRange(int j, int k,bool useboundaryfile, const DoubleVec1d& boundaryx, 
  const DoubleVec1d& boundaryy,const Mapgrid& mymapgrid);
bool SumCountsLegal(const DoubleVec1d& sumcounts);
bool SumCountsLegalWithCutoff(const DoubleVec1d& sumcounts, double cutoff);

bool ProgramStateValid(const DoubleVec1d& voronoiX, const DoubleVec1d& voronoiY,
  const IntVec1d& voronoiZ, const IntVec2d& bestpoint, const DoubleVec2d& dist,
  const IntVec2d& region, const DoubleVec1d& sumcounts, const DoubleVec3d& 
  counts, const DoubleVec2d& lcounts, int regionsize, const IntVec2d& mask,
  bool useboundaryfile, const DoubleVec1d& boundaryx, const DoubleVec1d& 
  boundaryy, const Mapgrid& mymapgrid);

void error_and_exit(const string& msg);

double dist_between(double x0, double y0, double x1, double y1);
std::pair<int,int> make_legal(double x, double y);

void ComputeLCount(const DoubleVec3d& counts, const IntVec2d& bestpoint, DoubleVec2d& lcounts);
void ComputeTotalCounts(const DoubleVec2d& lcounts, DoubleVec1d& totalcounts);
void ComputeBestPointAndDist(const DoubleVec1d& vx, const DoubleVec1d& vy, IntVec2d& bestpoint, DoubleVec2d& dist);

void UpdateBestPointAndDist(const DoubleVec1d& vx, const DoubleVec1d& vy, int l, IntVec2d& bestpoint, DoubleVec2d& dist);

void ComputeSumCounts(const IntVec2d& Region, const DoubleVec3d& Counts, DoubleVec1d& SumCounts);
void ComputeSumCounts2(const DoubleVec2d& LCounts, const IntVec1d& vz, DoubleVec1d& SumCounts);

void ComputeRegion(const IntVec1d& vz, IntVec2d& Region, const IntVec2d& BestPoint);
void ComputeCellSizeInRange(const IntVec2d& BestPoint, IntVec1d& CellSize);

void OutputRegion(const IntVec2d& Region, ofstream& output, bool OnlyRange);
void ComputeRegionSizeInSav(const IntVec2d& Region, int& regionsize);

void UpdateZs(double Vprob, const IntVec1d& CellSize, int& regionsize, DoubleVec1d& sumcounts, const DoubleVec2d& lcounts, const DoubleVec1d& totalcounts, IntVec2d& region, IntVec1d& VoronoiZ, const IntVec2d& bestpoint, double cutoff);

void ReadScatInfile(ifstream& locatefile, const string& infile, int skip, int nmcmc, int ind, DoubleVec3d& COUNTS);
