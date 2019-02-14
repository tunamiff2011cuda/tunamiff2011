/*
 *      Tohoku University's Numerical-Analysis Model
 *          for Investigation of tsunami
 *           Far-field Tsunami version 
 */

#ifndef DEFORMINITIA_H
#define DEFORMINITIA_H

#include "GRDUTIL.h"
#include "SPHERECAL.h"
#include "OFAULTPARA.h"


class DEFORMINITIA
{
protected:

  int finalized;
  int getDeformArea( int round, double& lonmin, double& lonmax, double& latmin, double& latmax );
  int setGrid( GRDUTIL& u );

public:

  int nfault;            
  double m0;             
  OFAULTPARA *fault;    

  DEFORMINITIA();
  ~DEFORMINITIA();
  int read( char *fname );
  int finalizeInput();
  double getM0();
  double getMw();
  int calculate( double lon, double lat, double& uz );
  int calculate( double lon, double lat, double& uz, double& ulon, double &ulat );
  int calculate( cObsArray& arr );
  int calculate( GRDUTIL& uZ );
  int calculate( GRDUTIL& uZ, GRDUTIL& uLon, GRDUTIL& uLat );
  char *sprint();

};

#endif // DEFORMINITIA_H
