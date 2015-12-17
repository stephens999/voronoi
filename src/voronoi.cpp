//#include "input.hpp"
#include "utility.hpp"

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <numeric>
#include <cmath>

using namespace std; 
const double PI = 3.141592; 

char TAG = 'r';  
double YMIN = -17;  // note that in this code, X and Y are mixed up
double YMAX = 50; // X represents N-S and Y represent E-W. 
double XMIN = -36; // Except in the functions InElephantRange, InForest etc where they have the opposite correspondance
double XMAX = 31;  // [there is a switch in InRange to allow for this]
int GRIDSIZE = 67;
int NITER = 10000;
int BURNIN = 5000;
int VLENGTH = 100; // length of vector of Voronoi points
int SAVANNAHONLY = 1;

int OFFSET = 0;

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


bool InRange(double x,double y){
  if(SAVANNAHONLY)
     return (InElephantRange(PI*y/180,PI*x/180) && !InForest(PI*y/180,PI*x/180));
  else
     return ((InElephantRange(PI*y/180,PI*x/180) && InForest(PI*y/180,PI*x/180)));
}


// test whether grid square (j,k) is in required range (Forest or Savannah)
bool GridInRange(int j, int k){
  double x = XCenterpoint(j);
  double y = YCenterpoint(k);
  return InRange(x,y);
}

// number of (valid savannah/forest  counts in each voronoi cell
void ComputeLCount(vector<vector< vector<double> > > & counts, vector<vector<int> > & bestpoint, vector<vector<double> > & lcounts){
  for(int i=0;i<counts.size();i++){
    for(int l = 0; l<VLENGTH; l++){
      lcounts[i][l] = 0;
    }
    for(int j=0;j<GRIDSIZE; j++){
      for(int k=0; k<GRIDSIZE; k++){
	if(GridInRange(j,k)){
	  lcounts[i][bestpoint[j][k]] += counts[i][j][k];
	}
      }
    }
  }
}

void ComputeTotalCounts( vector<vector<double> > & lcounts,vector<double> & totalcounts){
  for(int i=0; i<lcounts.size(); i++){
    for(int l=0; l<VLENGTH;l++){
      totalcounts[i] += lcounts[i][l];
    }
  }
}


void ComputeBestPointAndDist(vector<double> & vx,vector<double> & vy,vector<vector<int> > & bestpoint, vector<vector<double> > & dist){
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
void UpdateBestPointAndDist(vector<double> & vx,vector<double> & vy,int l, vector<vector<int> > & bestpoint, vector<vector<double> > & dist){
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

void ComputeSumCounts(vector<vector<int> > & Region, vector<vector<vector<double> > > & Counts,vector<double>  & SumCounts){
  for(int i = 0;i<Counts.size(); i++){
    SumCounts[i]=0;
    for(int j=0; j<GRIDSIZE; j++){
      for(int k=0; k<GRIDSIZE; k++){
	if(GridInRange(j,k))
	  SumCounts[i] += Region[j][k] * Counts[i][j][k];
      }
    }
  }
}
  
void ComputeSumCounts2(vector<vector<double> >  & LCounts,vector<int> & vz, vector<double>  & SumCounts){
  for(int i = 0;i<LCounts.size(); i++){
    SumCounts[i]=0;
    for(int l=0; l<VLENGTH; l++){
      SumCounts[i] += LCounts[i][l]*vz[l];
    }
  }
}
  


// compute grid of 1s and 0s corresponding to a voronoi tesselation
void ComputeRegion(vector<int> & vz,vector<vector<int> > & Region, vector<vector<int> > & BestPoint){
  for(int j=0; j<GRIDSIZE; j++){
    for(int k=0; k<GRIDSIZE; k++){
      Region[j][k] = vz[BestPoint[j][k]];
    }
  }
}

void ComputeCellSizeInRange(vector<vector<int> > & BestPoint, vector<int> & CellSize){
  for(int l =0;l<VLENGTH; l++)
    CellSize[l] = 0;
    
  for(int j=0; j<GRIDSIZE; j++){
    for(int k=0; k<GRIDSIZE; k++){
      if(GridInRange(j,k))
	CellSize[BestPoint[j][k]]++;
    }
  }
}

void OutputRegion(vector<vector<int> > & REGION,ofstream & output, bool OnlyRange){
 
  for(int j = 0; j<GRIDSIZE; j++){
    for(int k=0; k<GRIDSIZE; k++){
      int inrange = GridInRange(j,k) || !OnlyRange;
      output << inrange*REGION[j][k] << " ";
    }
    output << endl;
  }
}


void ComputeRegionSizeInSav(vector<vector<int> > & REGION, int & REGIONSIZE){
  REGIONSIZE = 0;
  for(int j = 0; j<(GRIDSIZE); j++){
    for(int k = 0; k<(GRIDSIZE); k++){
      if(GridInRange(j,k)){
	REGIONSIZE += REGION[j][k];
      }
    }
  }
}

void UpdateZs(double Vprob, vector<int> & CellSize, int REGIONSIZE, vector<double> & SUMCOUNTS, vector<vector<double> > & LCOUNTS, vector<double> & TOTALCOUNTS, vector<vector<int> > & REGION,vector<int> & VoronoiZ, vector<vector<int> > & BESTPOINT){
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
    if(newregionsize>0){
      for(int i = 0; i<NIND; i++){
	if(VoronoiZ[l] == 0)
	  newsumcounts = SUMCOUNTS[i] + LCOUNTS[i][l];
	else 
	  newsumcounts = SUMCOUNTS[i] - LCOUNTS[i][l];
	
	newloglik += log(newsumcounts/(TOTALCOUNTS[i] * newregionsize));
	currentloglik += log(SUMCOUNTS[i]/(TOTALCOUNTS[i] * REGIONSIZE));
      }
      
      if(ranf() < exp(newloglik - currentloglik)){
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

int main ( int argc, char** argv)
{
  int SEED = 0;
  map<string, string> filenames; 
  while( ( argc > 1 ) && ( argv[1][0] == '-' ) ) {
    switch(argv[1][1]) {
      
    case 't': //tag to use for location output files
      ++argv;
      --argc;
      TAG = argv[1][0];
      break;

   case 'S': // seed
      ++argv;
      --argc;
      SEED = atoi(&argv[1][0]);
       cout << "Seed = " << SEED << endl;
    break;

   case 'D':
	SAVANNAHONLY = 0;
	break;

   case 'd':
	SAVANNAHONLY = 1;
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

  srandom(SEED);

  int nmcmc = 100;
  int skip = 100;
  //int firstsample = atoi(argv[1]);
  //int lastsample = atoi(argv[2]);
  filenames["samples"] = argv[1];
  filenames["output"]  = argv[2];
  

  string regionfilename = filenames["output"] + "_regions";
  ofstream regionfile (regionfilename.c_str());
  string indprobfilename = filenames["output"] + "_indprobs";
  ofstream indprobfile (indprobfilename.c_str());
  string samplefilename = filenames["samples"];
  ifstream samplefile(samplefilename.c_str());

  vector<double> VoronoiX(VLENGTH,0);
  vector<double> VoronoiY(VLENGTH,0);
  vector<int> VoronoiZ(VLENGTH,0);
  vector<int> CellSize(VLENGTH,0); // number of grid squares in each voronoi cell


  int NIND; // = lastsample - firstsample + 1;
 
  // read in samples to be located from samplefile
  samplefile >> NIND;
  vector<int> samplevec(NIND,0);
  for(int s=0;s<NIND;s++)
    samplefile >> samplevec[s];


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


  //initialise REGION
  //for(int i = 0; i<NITER; i++){
  for(int l = 0; l<VLENGTH; l++){
   cout << "l=" << l << endl; 
   VoronoiX[l] = XMIN + ranf() * (XMAX-XMIN);
    VoronoiY[l] = YMIN + ranf() * (YMAX-YMIN);
    VoronoiZ[l] = (ranf()<0.1);
  }
  ComputeBestPointAndDist(VoronoiX,VoronoiY,BESTPOINT,DIST);
ComputeRegion(VoronoiZ,REGION,BESTPOINT);  
  ComputeRegionSizeInSav(REGION,REGIONSIZE);
 ComputeCellSizeInRange(BESTPOINT,CellSize);
  //OutputRegion(REGION,regionfile,false);
  //}

//   int SAVANNAHSIZE = 0;
//   double PSEUDOCOUNTS = 0.001;
//   for(int j = 1; j<(GRIDSIZE-1); j++){
//     for(int k = 1; k<(GRIDSIZE-1); k++){
//       if(GridInRange(j,k)){
// 	SAVANNAHSIZE+=1;
// 	for(int i = 0; i<NIND; i++){
// 	  COUNTS[i][j][k] = PSEUDOCOUNTS;
// 	  SUMCOUNTS[i] += REGION[j][k]*PSEUDOCOUNTS;
// 	}
//       }
//     }
//   }

//   cout << "Savannahsize = " << SAVANNAHSIZE << endl;  
//   double TOTALCOUNTS = nmcmc + SAVANNAHSIZE * PSEUDOCOUNTS;
  
  
  for(char TAG = 'r'; TAG<='z'; TAG++){
cout << "Reading " << TAG << endl;
     for(int s = 0; s < NIND; s++){
      int sample = samplevec[s];
      char firstdigit = (char) ('0' + (sample+OFFSET) / 100);
      char seconddigit = (char) ('0' +((sample+OFFSET) % 100) / 10);
      char thirddigit = (char)  ('0' +((sample+OFFSET) % 10));    
      
      char LOCATEFILE [5];
      LOCATEFILE[0] = firstdigit;
      LOCATEFILE[1] = seconddigit;
      LOCATEFILE[2] = thirddigit;
      LOCATEFILE[3] = TAG;
      LOCATEFILE[4] = 0;
      
      ifstream locatefile (LOCATEFILE);
      
      double x,y,dummy;
      
      for(int index = 0; index<skip; index++){
	locatefile >> x;
	locatefile >> y;
	locatefile >> dummy;
      }
      
      for(int index = 0; index<nmcmc; index++){
	locatefile >> x;
	locatefile >> y;
	locatefile >> dummy;
	int isinrange = InRange(x,y);
	if(!isinrange){
	  cerr << "Warning: Not in range" << endl;
	  cerr << "Sample =" << samplevec[s] << endl;
	  cerr << "Index = " << index << endl;
          cerr << "x= " << PI*x/180 << endl;
	  cerr << "y= " << PI*y/180 << endl;
	  //exit(1);
	}
	
	if(x>XMAX || x<XMIN || y>YMAX || y<YMIN){
	  cerr << "Error: x or y out of bounds" << endl;
	  cerr << "file = " << LOCATEFILE << endl;
	  cerr << "x= " << x << endl;
	  cerr << "y= " << y << endl;
	  exit(1);
	}
	
	int i = s; //sample - firstsample;
	int j = (int) trunc(GRIDSIZE * (x-XMIN)/(XMAX-XMIN));
	int k = (int) trunc(GRIDSIZE * (y-YMIN)/(YMAX-YMIN));
	
	//if(!GridInRange(j,k)){
	//	cerr << "Error: grid square out of bounds" << endl;
	//exit(1);
	//}
	
	COUNTS[i][j][k] += 1;
	if(REGION[j][k]==1 && GridInRange(j,k)){
	  SUMCOUNTS[i] +=1;
	}
	
      }
    }
  }

  ComputeLCount(COUNTS,BESTPOINT,LCOUNTS);
  ComputeTotalCounts(LCOUNTS,TOTALCOUNTS);
  
  int ACCEPT=0;
  
  double newloglik, currentloglik;

  double Vprob = 0.5; // prob of each voronoi point being a 1

  for(int iter = 0; iter<NITER; iter++){
    cout << "iter = " << iter << endl;

    UpdateZs(Vprob,CellSize,REGIONSIZE,SUMCOUNTS,LCOUNTS,TOTALCOUNTS,REGION,VoronoiZ,BESTPOINT);

    
    double newVprob = Vprob + rnorm(0,0.02);
    if(newVprob>0 && newVprob<1){
      double numZ = 0;
      for(int l=0; l<VLENGTH; l++)
	numZ += VoronoiZ[l];
      newloglik = numZ*log(newVprob) + (VLENGTH-numZ) * log(1-newVprob);
      currentloglik = numZ*log(Vprob) + (VLENGTH-numZ) * log(1-Vprob);
      if(ranf()< exp(newloglik - currentloglik))
	Vprob = newVprob;
      cout << "Vprob = " << Vprob << endl;
    }


    vector<double> oldVoronoiX(VoronoiX);
    vector<double> oldVoronoiY(VoronoiY);
    vector<vector<double> > oldLCOUNTS(LCOUNTS);
    vector<double> oldSUMCOUNTS(SUMCOUNTS); // sum of number of times each individual is in REGION
    vector<vector<int> > oldREGION(REGION);
    vector<vector<int> > oldBESTPOINT(BESTPOINT);
    vector<vector<double> >  oldDIST(DIST);
    
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
      if(SUMCOUNTS[i]==0)
	reject = true;
      newloglik += log(SUMCOUNTS[i]/(TOTALCOUNTS[i] * REGIONSIZE));
    }
    
   //  cout << "N: " << newloglik << endl;
//     cout << "C: " << currentloglik << endl; 
//     cout << "A: " << ACCEPT << endl;

    if(reject || ranf() > exp(newloglik - currentloglik)){ // reject move
      VoronoiX = oldVoronoiX;
      VoronoiY = oldVoronoiY;
      LCOUNTS = oldLCOUNTS;
      SUMCOUNTS = oldSUMCOUNTS;
      REGION = oldREGION;
      BESTPOINT = oldBESTPOINT;
      DIST = oldDIST;
    } else{
      ComputeCellSizeInRange(BESTPOINT,CellSize);
      ACCEPT +=1;
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
  
  string outputfilename = filenames["output"];
  ofstream output (outputfilename.c_str());
  
  for(int j = 0; j<GRIDSIZE; j++){
    for(int k=0; k<GRIDSIZE; k++){
      output << PROB[j][k]/(NITER - BURNIN) << " ";
    }
    output << endl;
  }

  for(int i =0; i<NIND; i++){
    for(int j = 0; j<GRIDSIZE; j++){
      for(int k=0; k<GRIDSIZE; k++){
	indprobfile << INDPROBS[i][j][k]/(NITER - BURNIN) << " ";
      }
      indprobfile << endl;
    }
  }

  
}




//     // if(iter>BURNIN){ 
// //       for(int j = 0; j<GRIDSIZE; j++){
// // 	for(int k=0; k<GRIDSIZE; k++){
// // 	  int insav = GridInRange(j,k);
// // 	  indprobfile << COUNTS[1][j][k] << " ";
// // 	}
// // 	indprobfile << endl;
// //       }
// //       for(int j = 0; j<GRIDSIZE; j++){
// // 	for(int k=0; k<GRIDSIZE; k++){
// // 	  int insav = GridInRange(j,k);
// // 	  indprobfile << insav << " ";
// // 	}
// // 	indprobfile << endl;
// //       }

// //       indprobfile << SUMCOUNTS[1] << endl;
// //     }   
    
//   }
  
  
//   string outputfilename = filenames["output"];
//   ofstream output (outputfilename.c_str());
  
//   for(int j = 0; j<GRIDSIZE; j++){
//     for(int k=0; k<GRIDSIZE; k++){
//       output << PROB[j][k]/(NITER - BURNIN) << " ";
//     }
//     output << endl;
//   }




//   cout << "Regionsize = " << REGIONSIZE << endl;

// }
  
