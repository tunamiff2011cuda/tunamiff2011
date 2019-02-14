/*
 *      Tohoku University's Numerical-Analysis Model
 *          for Investigation of tsunami
 *           Far-field Tsunami version 
 */

#ifndef OFAULTPARA_H
#define OFAULTPARA_H

#define FLT_POS_C  0      
#define FLT_POS_MT 1      
#define FLT_POS_BT 2      
#define FLT_POS_BB 3      
#define FLT_POS_MB 4      

#define FLT_ERR_DATA      1
#define FLT_ERR_ZTOP      3
#define FLT_ERR_INTERNAL  4
#define FLT_ERR_STRIKE    6

class OFAULTPARA
{

protected:

  int checked;
  int adjust_depth;
  double sind;
  double cosd;
  double sins;
  double coss;
  double tand;
  double coslat;
  double zbot;
  double wp;
  double dslip;
  double sslip;

  double mw2m0();

public:

  int refpos;                 
  double mw;
  double slip;
  double lon,lat,depth;
  double strike;
  double dip;
  double rake;
  double length,width;
  double mu;

  OFAULTPARA();
  ~OFAULTPARA();
  int read( char *faultparam );
  int check();
  double getMw();
  double getM0();
  double getZtop();
  int global2local( double glon, double glat, double& lx, double& ly );
  int local2global( double lx, double ly, double& glon, double& glat );
  int getDeformArea( double& lonmin, double& lonmax, double& latmin, double& latmax );
  int calculate( double lon, double lat, double& uz, double& ulon, double &ulat );
  int calculate( double lon, double lat, double& uz );
};

int DEFORMOKADA( double L,double W,double D,double sinD,double cosD,double U1,double U2,double x,double y,int flag_xy, double *Ux,double *Uy,double *Uz );

#endif // OFAULTPARA_H
