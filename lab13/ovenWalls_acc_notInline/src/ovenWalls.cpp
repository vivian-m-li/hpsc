//  ================================================================================
//  ||                                                                            ||
//  ||              ovenWalls                                                     ||
//  ||              ------------------------------------------------------        ||
//  ||              T H E R M A L   R A D I A T I O N                             ||
//  ||                                                                            ||
//  ||              D E M O N S T R A T I O N   C O D E                           ||
//  ||              ------------------------------------------------------        ||
//  ||                                                                            ||
//  ||       Developed by: Scott R. Runnels, Ph.D.                                ||
//  ||                     University of Colorado Boulder                         ||
//  ||                                                                            ||
//  ||                For: CU Boulder CSCI 4576/5576 and associated labs          ||
//  ||                                                                            ||
//  ||           Copyright 2020 Scott Runnels                                     ||
//  ||                                                                            ||
//  ||                     Not for distribution or use outside of the             ||
//  ||                     this course.                                           ||
//  ||                                                                            ||
//  ================================================================================

#include "ovenWalls.h"

class Cylinder
{
public:

  int node;
  int nCella;       // Number of cells in the axial direction
  int nCellc;       // Number of cells in the circular direction
  int nReala;       // Number of nodes in the axial direction
  int nRealc;       // Number of nodes in the circular direction
  int nField;       // Number of nodes
  double dtheta;    // Angular spacing of nodes, in degrees
  double dz;        // Axial spacing of nodes
  double radius;    // Cylinder radius
  double length;    // Cylinder length
  int    **face;    // Faces:  face[f][i] gives the nodes of face f
  int    nFaces;    // Number of faces
  double **coord;   // Vertices
  double **normal;  // Face normals
  double **center;  // Face centers

  Cylinder(){}
  
  Cylinder( int _nCella , int _nCellc , double _radius , double _length , int Out1In2 )
  {
    // Store user inputs
    
    nCella = _nCella;
    nCellc = _nCellc;
    
    nReala = nCella + 1 ;
    nRealc = nCellc     ;
  
    nField = nReala * nRealc; 
    nFaces = nCella * nCellc;

    radius = _radius;
    length = _length;
    
    dtheta = 360.   / nCellc;  
    dz     = length / nCella;

    // Form mesh

    coord  = Array2D_double(nField + 1 , 3);
    face   = Array2D_int   (nFaces + 1 , 5);
    normal = Array2D_double(nFaces + 1 , 3);
    center = Array2D_double(nFaces + 1 , 3);

    int faceCount = 0;
    
    for ( int j = 1 ; j <= nRealc ; ++j )
    for ( int i = 1 ; i <= nReala ; ++i )
      {
	int p = pid(i,j);

	double theta = dtheta * j;
	double zval  = dz     * (i-1);

	coord[p][0] = radius * cos(theta*3.1415/180.);
	coord[p][1] = radius * sin(theta*3.1415/180.);
	coord[p][2] = zval;

	// Form face, i.e., for this face number (which happens to be p)
	// collect the pids of the 4 nodes that comprise this face.

	int point[5];
	point[1] = p    ;
	point[2] = p + 1;
	point[3] = point[2] + nReala;
	point[4] = point[1] + nReala;

	// Correct point for when we have wrapped completely around the circle
	
	if ( j == nRealc )
	  {
	    point[3] = pid(i+1,1);
	    point[4] = pid(i,1);
	  }

	// Store the point values in this face, p (but not for the last point in an i-row

	if ( i < nReala)
	  {
	    ++faceCount;
	    for ( int i = 1 ; i <= 4 ; ++i ) face[faceCount][i] = point[i];
	  }
      }

    // Compute normals for each face
    
    for ( int f = 1 ; f <= nFaces ; ++f )
      {
	double xav = 0. , yav = 0. , zav = 0.;
	for ( int k = 1 ; k <= 4 ; ++k )
	  {
	    xav += coord[ face[f][k] ][0]  ;
	    yav += coord[ face[f][k] ][1]  ;
	    zav += coord[ face[f][k] ][2]  ;
	  }
	
	xav /= 4.;
	yav /= 4.;
	zav /= 4.;
	
	double mag = sqrt((xav*xav + yav*yav));
	
	normal[f][0] = xav/mag;
	normal[f][1] = yav/mag;
	normal[f][2] = 0.;
	center[f][0] = xav;
	center[f][1] = yav;
	center[f][2] = zav;

	if ( Out1In2 == 2 ) for ( int k = 0 ; k < 3 ; ++k ) normal[f][k] *= -1.;

      }
    
  }
    


  #include "plotter.h"
  int pid(int i , int j ) { return ( i + (j-1) * nReala); }

};


#include "geom.h"



//  ==
//  ||
//  ||
//  ||  Main Program
//  ||
//  ||
//  ==

int main(int argc, char *argv[])
{

  printf("\n");
  printf("\n");
  printf("Ray Tracing Demo Code\n");
  printf("\n");
  printf("\n");
  
  double length = 10.;
  double radius = 10.0;
  //             nCell   nCell                                  Facing In
  //             Axial   Angular     Radius         Length      Or Out
  //            ------  --------    ----------     ----------   -------
  Cylinder CylA(   20  ,  20     ,   radius     ,   length      ,   2      );
  Cylinder CylB(   20  ,  20     ,   radius/2.  ,   length      ,   1      );

  CylA.plot("CylA");
  CylB.plot("CylB");

  // Set up parameters for ray tracing

  VI facesBlocking;
  VI facesHit;

  // R A Y   T R A C E   F O R   S O U R C E   F A C E  1 to target face 100


  double *Ray0 = (double*) malloc ( sizeof(double)*3 );
  double *Ray  = (double*) malloc ( sizeof(double)*3 );
  int sourceFaceID = 5;
  int targetFaceID = 51;

  CylA.plotFace("source",sourceFaceID);
  
  // (1) Point on source face
    
  kLOOP Ray0[k] = CylA.center[ sourceFaceID ] [k];

  // (2) Ray to target
  
  int blockerFaceID = 0;
  double dotProduct;

  // Begin------------------------ Cylinder Mesh Data Structures

  int nFacesA = CylA.nFaces;
  int nFacesB = CylB.nFaces;
  int nFieldB = CylB.nField;
  
  double **centerA = Array2D_double( nFacesA + 1 , 3);
  double **normalA = Array2D_double( nFacesA + 1 , 3);
  
  double **normalB = Array2D_double( nFacesB + 1 , 3);
  int    **faceB   = Array2D_int   ( nFacesB + 1 , 5);
  
  double **coordB  = Array2D_double( nFieldB + 1 , 3);
  

  for ( int f = 1 ; f <= nFacesA ; ++f ) kLOOP centerA[f][k] = CylA.center[f][k];
  for ( int f = 1 ; f <= nFacesA ; ++f ) kLOOP normalA[f][k] = CylA.normal[f][k];
  for ( int f = 1 ; f <= nFacesB ; ++f ) kLOOP normalB[f][k] = CylB.normal[f][k];
  for ( int p = 1 ; p <= nFieldB ; ++p ) kLOOP coordB [p][k] = CylB.coord [p][k];
  for ( int f = 1 ; f <= nFacesB ; ++f ) for ( int i = 1 ; i <= 4 ; ++i) faceB[f][i] = CylB.face[f][i];


  // End-------------------------- Cylinder Mesh Data Structures

  // Begin------------------------ Block and Hit records

  int *A_wasHit  = new int [ nFacesA + 1 ];
  int *B_blocked = new int [ nFacesB + 1 ];

  for ( int f = 1 ; f <= nFacesA ; ++f ) A_wasHit [f] = 0;
  for ( int f = 1 ; f <= nFacesB ; ++f ) B_blocked[f] = 0;
    
  // End-------------------------- Block and Hit records
  
  
  // ====
  //  ||
  //  ||
  //  || BEGIN: Ray Trace from Cylinder A's sourceFaceID to all of Cylinder A's faces
  //  ||
  //  ||
  // ====

  double * Q0     = (double*) malloc ( sizeof(double)*3 );
  double * Q1     = (double*) malloc ( sizeof(double)*3 );
  double * Q2     = (double*) malloc ( sizeof(double)*3 );
  double * Q3     = (double*) malloc ( sizeof(double)*3 );
  double * normal = (double*) malloc ( sizeof(double)*3 );
  double * xyzInt = (double*) malloc ( sizeof(double)*3 );
  
  struct timespec t0, t1;

  int *blocked = new int [nFacesB+1];

  // =========================================================================================      
  // GPU: Timed
  StartTimer(t0);
  // =========================================================================================      


  for ( targetFaceID = 1 ; targetFaceID <= nFacesA ; ++targetFaceID )
    {
      kLOOP Ray [k] = centerA[ targetFaceID ] [k] - Ray0[k];

      // (4) Find blockers of ray

      for ( int f = 1 ; f <= nFacesB ; ++f ) blocked[f] = 0;

      dotProduct = 0.; kLOOP dotProduct += normalA[sourceFaceID][k]*normalA[targetFaceID][k];


      if ( dotProduct < .9999  )
	{
	  blockerFaceID = 0;

	  double dotProduct;
	  double dotProduct2;
	  
#pragma acc enter data copyin( coordB[ 0:nFieldB+1 ][0:3], faceB[ 0:nFacesB+1 ][0:5], normalB[ 0:nFacesB+1 ][0:3], blocked [ 0:nFacesB+1 ] , Ray[ 0:3], Ray0 [ 0:3] )
#pragma acc parallel loop private(  Q0[0:3], Q1[0:3], Q2[0:3], Q3[0:3], xyzInt[0:3], normal[0:3], dotProduct, dotProduct2 ) \
                          present( coordB [ 0:nFieldB+1 ][0:3], faceB  [ 0:nFacesB+1][0:5], normalB[ 0:nFacesB+1][0:3], blocked[0:nFacesB+1], Ray[0:3], Ray0[0:3] )
	  
	  for ( int BID = 1 ; BID <= nFacesB ; ++BID )
	    {
	      // Vertices of potential blocker
      
	      kLOOP Q0[k] = coordB[ faceB[ BID ][ 1 ] ][k];
	      kLOOP Q1[k] = coordB[ faceB[ BID ][ 2 ] ][k];
	      kLOOP Q2[k] = coordB[ faceB[ BID ][ 3 ] ][k];
	      kLOOP Q3[k] = coordB[ faceB[ BID ][ 4 ] ][k];

	       // Normal of potential blocker
      
	      kLOOP normal[k] = normalB[ BID   ][k];  

              if ( LineHitsFace( Ray0 , Ray , Q0 , Q1 , Q2 , Q3 , normal) > 0 ) blocked[BID] += 1;

	      
	    }
#pragma acc exit data copyout(blocked [0:nFacesB+1],Ray[0:3],Ray0[0:3])

	  int sumBlocked = 0;
	  for ( int f = 1 ; f <= nFacesB ; ++f ) sumBlocked += blocked[f];
	  
	  if ( sumBlocked == 0 ) A_wasHit [targetFaceID ] = 1;  
	}

    }

  
  // =========================================================================================      
  // GPU: Timed 
  EndTimer("GPU", t0, t1);
  // =========================================================================================      

  CylA.plotFacesIndicated("hit_gpu",A_wasHit);
  

  for ( int f = 1 ; f <= nFacesA ; ++f ) A_wasHit [f] = 0;
  for ( int f = 1 ; f <= nFacesB ; ++f ) B_blocked[f] = 0;


  // =========================================================================================      
  // CPU: Timed
  StartTimer(t0);
  // =========================================================================================      


  for ( targetFaceID = 1 ; targetFaceID <= nFacesA ; ++targetFaceID )
    {
      kLOOP Ray [k] = centerA[ targetFaceID ] [k] - Ray0[k];

      // (4) Find blockers of ray

      for ( int f = 1 ; f <= nFacesB ; ++f ) blocked[f] = 0;

      dotProduct = 0.; kLOOP dotProduct += normalA[sourceFaceID][k]*normalA[targetFaceID][k];


      if ( dotProduct < .9999  )
	{
	  blockerFaceID = 0;

	  int blockHasOccurred = 0;

	  while ( ! blockHasOccurred  && ++blockerFaceID <= nFacesB )
	    {
	      // Vertices of potential blocker
      
	      kLOOP Q0[k] = coordB[ faceB[ blockerFaceID ][ 1 ] ][k];
	      kLOOP Q1[k] = coordB[ faceB[ blockerFaceID ][ 2 ] ][k];
	      kLOOP Q2[k] = coordB[ faceB[ blockerFaceID ][ 3 ] ][k];
	      kLOOP Q3[k] = coordB[ faceB[ blockerFaceID ][ 4 ] ][k];

	      // Normal of potential blocker
      
	      kLOOP normal[k] = normalB[ blockerFaceID ][k];  

	      if ( LineHitsFace( Ray0 , Ray , Q0 , Q1 , Q2 , Q3 , normal) > 0 ) { blocked[blockerFaceID] = 1;  blockHasOccurred = 1; }

	    }
	  
	  if ( blocked > 0  ) B_blocked[blockerFaceID] = 1;

	  int sumBlocked = 0;
	  for ( int f = 1 ; f <= nFacesB ; ++f ) sumBlocked += blocked[f];
	  
	  if ( sumBlocked == 0 ) A_wasHit [targetFaceID ] = 1;  
	}

    }

  
  // =========================================================================================      
  // CPU: Timed 
  EndTimer("CPU", t0, t1);
  // =========================================================================================      










  // ====
  //  ||
  //  ||
  //  || END: Ray Trace
  //  ||
  //  ||
  // ====
  

  // Plot hits and blockers
  
  CylA.plotFacesIndicated("hit",A_wasHit);
  CylB.plotFacesIndicated("blockers",B_blocked);
  
  return 0;

}
