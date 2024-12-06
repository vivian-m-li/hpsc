//  ==
//  ||
//  ||
//  ||  Geometric Utilities
//  ||
//  ||
//  ==

#define _IC_MACRO_  kLOOP  { vAB[k] =  B[k] - A[k] ;   vCA[k] =  A[k] - C[k] ;   vAP[k] =  P[k] - A[k] ;  }  cB[0] = vAB[1]*vAP[2] - vAB[2]*vAP[1];  cB[1] = vAB[2]*vAP[0] - vAB[0]*vAP[2];  cB[2] = vAB[0]*vAP[1] - vAB[1]*vAP[0];  cC[0] = vCA[1]*vAP[2] - vCA[2]*vAP[1];  cC[1] = vCA[2]*vAP[0] - vCA[0]*vAP[2];  cC[2] = vCA[0]*vAP[1] - vCA[1]*vAP[0];  dot = 0. ; kLOOP dot += cB[k]*cC[k];

//  ==
//  ||
//  ||  dot: Computes the dot product of v and w
//  ||
//  ==

double Dot( double *v , double *w )
{
  return v[0]*w[0] + v[1]*w[1] + v[2]*w[2] ;     //  ~LabReplacePhrase~  ~n~    ~;~    ~   TO-DO  ~
}


//  ==
//  ||
//  ||  LineHitsFace: Returns 1 if the line specified by a point L0 and vector L
//  ||                intersect a face specified by four points Q0 - Q3.  
//  ||
//  ==

int LineHitsFace( double *Ray0 , double *Ray , double *Q0 , double *Q1 , double *Q2 , double *Q3 , double *normal)
{
  double xyzInt[3],  p0Ray0[3];
  double dotProduct;
  double numerator;
  double d;
  double A[3], B[3], C[3], vAB[3], vCA[3], vAP[3], cB[3], cC[3];
  double dot;
  double P[3];
  int    inside;

  // Compute intersection of the line Ray0/Ray with the plane with normal n
  
  kLOOP p0Ray0[k] = Q0[k] - Ray0[k];
  dotProduct = 0.;  kLOOP dotProduct += Ray[k]*normal[k];
  numerator  = 0.;  kLOOP numerator  += p0Ray0[k]*normal[k];
  if ( fabs(dotProduct) < 1.e-10 ) { kLOOP xyzInt[k] = 1.e+12; }  else { d = numerator / dotProduct; kLOOP xyzInt[k] = Ray0[k] + Ray[k] * d;}
  
  // Determine if xyzInt is inside the face (Q0, Q1, Q2, Q3)

  inside = 0;
  kLOOP P[k] = xyzInt[k];
  kLOOP    {      A[k] = Q0[k] ;  B[k] = Q1[k] ;   C[k] = Q3[k] ;     }  _IC_MACRO_  if ( dot >= 0. ) inside += 1;
  kLOOP    {      A[k] = Q1[k] ;  B[k] = Q0[k] ;   C[k] = Q2[k] ;     }  _IC_MACRO_  if ( dot >= 0. ) inside += 1;
  kLOOP    {      A[k] = Q2[k] ;  B[k] = Q1[k] ;   C[k] = Q3[k] ;     }  _IC_MACRO_  if ( dot >= 0. ) inside += 1;
  kLOOP    {      A[k] = Q3[k] ;  B[k] = Q2[k] ;   C[k] = Q0[k] ;     }  _IC_MACRO_  if ( dot >= 0. ) inside += 1;


  if (inside==4) return 1;
  
  return 0;


}
