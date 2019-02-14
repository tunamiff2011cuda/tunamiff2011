/*
 *      Tohoku University's Numerical-Analysis Model
 *          for Investigation of tsunami
 *           Far-field Tsunami version 
 */

#ifndef GRDUTIL_H
#define GRDUTIL_H

class GRDUTIL
{

protected:


public:
  int nx,ny;
  int nnod;
  double noval;
  double xmin,xmax;
  double ymin,ymax;
  double dx,dy;
  double *val;

  GRDUTIL();
  GRDUTIL( const GRDUTIL& );
  GRDUTIL( double xmin0, double xmax0, int nx0, double ymin0, double ymax0, int ny0 );
  GRDUTIL( double xmin0, double xmax0, double dx0, double ymin0, double ymax0, double dy0 );
  ~GRDUTIL();

  void setNoval( double val );
  double getNoval();
  int initialize( double xmin0, double xmax0, int nx0, double ymin0, double ymax0, int ny0 );
  int initialize( double xmin0, double xmax0, double dx0, double ymin0, double ymax0, double dy0 );
  int readShape( const char *grdfile );
  int readHeader( const char *grdfile );
  int readGRD( const char *grdfile );
  int readXYZ( const char *xyzfile );
  int readRasterStream( FILE *fp, int ordering, int ydirection );
  GRDUTIL* extract( int i1, int i2, int j1, int j2 );
  GRDUTIL* extract( double x1, double x2, double y1, double y2 );
  int idx( int i, int j );
  double& operator() ( int i, int j );
  double& operator() ( int idx );
  GRDUTIL& operator= ( const GRDUTIL& grd );
  GRDUTIL& operator*= ( double coeff );
  GRDUTIL& operator+= ( GRDUTIL& grd );
  void getIJ( int idx, int& i, int& j );
  int getIJ( double x, double y, int& i, int& j );
  double getX( int i, int j );
  double getX( int idx );
  double getY( int i, int j );
  double getY( int idx );
  double getVal( int idx );
  double getVal( int i, int j );
  double getVal( double x, double y );
  double getMaxVal();
  double getMaxVal( int& i, int& j );
  double getMinVal();
  double getMinVal( int& i, int& j );
  double getMaxAbsVal();
  double getMaxAbsVal( int& i, int& j );
  double getMaxAbsValBnd();
  void setVal( double value, int i, int j );
  void setVal( double value, int idx );
  void reset();
  void resetVal();
  int getIntersectionRegion( const GRDUTIL &grd, int& imin, int& imax, int& jmin, int& jmax );
  int interpolateFrom( GRDUTIL &grd, int resetValues );
  int getNearestIdx( double x, double y );
  int getNearestIdx( double x, double y, double rangemin, double rangemax );
  int equalTo( GRDUTIL& grd );
  int isSameShapeTo( GRDUTIL& grd );
  void smooth( int radius );
  int writeGRD( const char *fname );
  int writeGRDbin( const char *fname );
  int writeXYZ( const char *fname );
};


#endif  // GRDUTIL_H
