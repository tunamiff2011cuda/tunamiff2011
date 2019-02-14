/*
 *      Tohoku University's Numerical-Analysis Model
 *          for Investigation of tsunami
 *           Far-field Tsunami version 
 */

#ifndef TUNAMI_FF_H
#define TUNAMI_FF_H

#define RADIEAR 6384.e+3          // Radius of the earth
#define GG 9.81         // Gravitational acceleration
#define Omega 7.29e-5        // Earth rotational period [1/sec]

#define MAX_VARS_PER_NODE 12

#define iD    0
#define iZ    1
#define iZmax 2
#define iM    3
#define iN    4
#define iR1   5
#define iR2   6
#define iR3   7
#define iR4   8
#define iR5   9
#define iTime 10
#define iTopo 11

struct TFPARAMS {
  char *oopsName;
  char *modelSubset;
  char *BathymetryInputfile;
  char *fileSource;
  char *fileGAGUEs;
  int DT;
  int time;
  int KEtotalTime;
  int KCgagueOUT;
  int GAGUERT;
  int MAXOUT;
  int outProgress;
  int KDspatialOut;
  int coriolis;
  float dmin;
  float GAGUEDIMX;
  float GAGUEDEMN;
  float GAGUEDEMX;
  float zout0TRel;
  float zout0TAbs;
  float zoutCT;
  float zoutZT;
  float zoutTT;
  float zoutAT;
  bool gpu;
  bool adjustZtop;
  bool verbose;
};

extern struct TFPARAMS MODELVARI;
extern int NLon,NLat;
extern double LonMin,LonMax,LatMin,LatMax;
extern double DLon,DLat;                 
extern double Dx,Dy;                     
extern float *R6;
extern float *C1;
extern float *C2;
extern float *C3;
extern float *C4;
extern int Imin;
extern int Imax;
extern int Jmin;
extern int Jmax;

#define idx(j,i) ((i-1)*NLat+j-1)
#define getLon(i) (LonMin+(i-1)*DLon)
#define getLat(j) (LatMin+(j-1)*DLat)

int tfRdepthParame();
int MODELPARAINITIA( int argc, char **argv );
void tfLogParams(void);
int RESETNODES();
int ICPROFPARA();
int MASSGBOUNDMOMENT();
int MASSGBOUNDMOMENTCor();

int PROPAOUTFL();
int FILEOTPROPA();
int PROPAMAX();
int tfReadGAGUEs();
int tfWriteGAGUEs();
int POINTGMAX();
int POINTGMAXCompact( int istage );

extern int NGAGUEs;
extern long *idxGAGUE;


#define printf_v( Args, ... )	if( MODELVARI.verbose ) printf( Args, ##__VA_ARGS__);

#include "NODEINI.h"

extern CNode *gNode;

#endif /* TUNAMI_FF_H */

