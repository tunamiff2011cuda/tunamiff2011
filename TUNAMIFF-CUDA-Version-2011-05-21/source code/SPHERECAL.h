/*
 *      Tohoku University's Numerical-Analysis Model
 *          for Investigation of tsunami
 *           Far-field Tsunami version 
 */

#ifndef ONSPHERE_H
#define ONSPHERE_H

#define RADIEAR 6384.e+3          // Earth radius

class cObsArray
{

public:

  int nPos;
  int nObs;
  char **id;
  double *lon;
  double *lat;
  double **obs;

  cObsArray();
  ~cObsArray();
  int read( char *fname );
  int write( char *fname );
  int resetObs();
  int resetObs( int newnobs );
  int findById( char *id0 );
  long writeBin( FILE *fp );
  long readBin( FILE *fp );
  double residual( cObsArray& ref );
  double norm();

};


double GeoDistOnSphere( double lon1, double lat1, double lon2, double lat2 );
double GeoStrikeOnSphere( double lon1, double lat1, double lon2, double lat2 );

#endif  // ONSPHERE_H
