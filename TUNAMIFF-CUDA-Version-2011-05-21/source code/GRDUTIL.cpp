/*
 *      Tohoku University's Numerical-Analysis Model
 *          for Investigation of tsunami
 *           Far-field Tsunami version
 *
 *Initialize with number of nodes
 *Initialize with grid step
 *only header read from GRD-file
 *shape read from GRD-file, reset values to zero
 *
 *Grid initialization from Golden Software GRD-file
 **check if bathymetry file is in ascii or binary format
 **Read Surfer GRD-file
 *
 *Grid initialization from XYZ-file
 **Open plain XYZ-file
 **Check xyz-file format, define number of nodes and min-max
 **Read first two lines and define if file is row- or column-wise
 **Define nx and ny
 **Other grid parameters dx, dy
 **Allocate memory for z-values
 **Read z-values
 *
 *Read grid from ASCII-raster file
 *Cut-off a subgrid
 *Total reset
 *Reset values to zero
 *
 *Get IJ-indices from the offset
 **J is the inner loop (increments faster than I)
 *
 *Get IJ-indices from coordinates
 *Get offset from IJ-indices
 *Get value at given node i, j, idx
 *Set value at given node i, j, idx
 *Get maximal value on a grid vmax
 *Get maximal value and its position on a grid l, lmax, vmax
 *Get minimal value on a grid l, vmin
 *Get minimal value and its position on a grid l, lmin, vmin
 *Get maximal absolute value on a grid l, vmax
 *Get maximal absolute value and its position on a grid l, lmax, vmax
 *Get maximal absolute value along grid boundaries i, j, vmax
 *Get value with interpolation (bi-linear) iO, jO, fx, fy, val_l, val_r, result
 *Get longitude from IJ-indeces
 *Get longitude from the offset
 *Get lattitude from IJ-indeces
 *Get lattitude from the offset
 *Check for the shapes are equal
 *Intersection region with another grid
 *Interpolate grid values from another O-grid (bi-linear) ierr,imin,imax,jmin,jmax, value
 *Get nearest node iO, jO, l, lmin, dist2, dist2min
 *Get nearest conditioned node i,j,i0,j0,rad,l,lmin,dist2,dist2min  
 *Smooth data on grid
 *write to GRD-file 
 *write to GRD-file: binary version
 *write in ascii 3 column-format 
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "UTILMDL.h"
#include "GRDUTIL.h"

#define DEFAULT_NOVAL 0.0
#define ROWWISE 1
#define COLWISE 2
#define TOP2BOT -1
#define BOT2TOP 1


int iabs( int i )
{
  return( (i<0) ? (-i) : i );
}

GRDUTIL::GRDUTIL( )
{
  xmin = xmax = dx = 0.;
  nx = 0;
  ymin = ymax = dy = 0.;
  ny = 0;
  nnod = 0;
  noval = DEFAULT_NOVAL;
  val = NULL;
}

GRDUTIL::GRDUTIL( const GRDUTIL& grd )
{
  xmin=grd.xmin; xmax=grd.xmax; dx=grd.dx; nx=grd.nx;
  ymin=grd.ymin; ymax=grd.ymax; dy=grd.dy; ny=grd.ny;
  nnod=grd.nnod;
  noval=grd.noval;

  val = new double[nnod];
  memcpy( val, grd.val, nnod*sizeof(double) );
}

GRDUTIL::GRDUTIL( double xmin0, double xmax0, int nx0, double ymin0, double ymax0, int ny0 )
{
  noval = DEFAULT_NOVAL;
  val = NULL;
  initialize( xmin0, xmax0, nx0, ymin0, ymax0, ny0 );
}
GRDUTIL::GRDUTIL( double xmin0, double xmax0, double dx0, double ymin0, double ymax0, double dy0 )
{
  noval = DEFAULT_NOVAL;
  val = NULL;
  initialize( xmin0, xmax0, dx0, ymin0, ymax0, dy0 );
}
GRDUTIL::~GRDUTIL()
{
  if(val) delete [] val;
}

int GRDUTIL::initialize( double xmin0, double xmax0, int nx0, double ymin0, double ymax0, int ny0 )
{

  xmin = xmin0;
  xmax = xmax0;
  nx = nx0;
  dx = (xmax-xmin)/(nx-1);

  ymin = ymin0;
  ymax = ymax0;
  ny = ny0;
  dy = (ymax-ymin)/(ny-1);

  nnod = nx*ny;

  if( val ) delete [] val;

  val = new double[nnod];
  for( int l=0; l<nnod; l++ ) val[l]=noval;

  return 0;
}
int GRDUTIL::initialize( double xmin0, double xmax0, double dx0, double ymin0, double ymax0, double dy0 )
{

  xmin = xmin0;
  xmax = xmax0;
  dx = dx0;
  nx = (int)((xmax-xmin)/dx + 1.e-10) + 1;

  ymin = ymin0;
  ymax = ymax0;
  dy = dy0;
  ny = (int)((ymax-ymin)/dy + 1.e-10) + 1;

  nnod = nx*ny;

  if( val ) delete [] val;

  val = new double[nnod];
  for( int l=0; l<nnod; l++ ) val[l]=noval;

  return 0;
}

void GRDUTIL::setNoval( double val )
{
  noval = val;
}


//=========================================================================
double GRDUTIL::getNoval()
{
  return noval;
}

int GRDUTIL::readHeader( const char *grdfile )
{
  FILE *fp;
  char dsaa_label[8];
  int ierr,isBin;
  short shval;
  float fval;
  double dval;

  if( (fp = fopen( grdfile, "rb" )) == NULL ) return Err.post(Err.msgOpenFile(grdfile));
  memset( dsaa_label, 0, 5 );
  ierr = fread( dsaa_label, 4, 1, fp );
  if( !strcmp( dsaa_label,"DSAA" ) )
    isBin = 0;
  else if( !strcmp( dsaa_label,"DSBB" ) )
    isBin = 1;
  else {
    fclose(fp);
    return Err.post(Err.msgReadFile( grdfile, 1, "not a GRD-file" ));
  }
  fclose(fp);

  if( isBin ) {
    fp = fopen( grdfile, "rb" );
    ierr = fread( dsaa_label, 4, 1, fp );
    ierr = fread( &shval, sizeof(short), 1, fp ); nx = shval;
    ierr = fread( &shval, sizeof(short), 1, fp ); ny = shval;
    ierr = fread( &xmin, sizeof(double), 1, fp ); ierr = fread( &xmax, sizeof(double), 1, fp );
    ierr = fread( &ymin, sizeof(double), 1, fp ); ierr = fread( &ymax, sizeof(double), 1, fp );
    ierr = fread( &dval, sizeof(double), 1, fp ); ierr = fread( &dval, sizeof(double), 1, fp ); 
  }
  else {
    fp = fopen( grdfile, "rt" );
    ierr = fscanf( fp, "%s", dsaa_label );
    ierr = fscanf( fp, " %d %d ", &nx, &ny );
    ierr = fscanf( fp, " %lf %lf ", &xmin, &xmax );
    ierr = fscanf( fp, " %lf %lf ", &ymin, &ymax );
    ierr = fscanf( fp, " %*s %*s " );   
  }

  fclose( fp );

  nnod = nx*ny;
  dx = (xmax-xmin)/(nx-1);
  dy = (ymax-ymin)/(ny-1);

  return 0;
}
int GRDUTIL::readShape( const char *grdfile )
{

  int ierr = readGRD( grdfile ); if(ierr) return ierr;

  for( int l=0; l<nnod; l++ ) val[l]=noval;

  return 0;
}

int GRDUTIL::readGRD( const char *grdfile )
{
  FILE *fp;
  char dsaa_label[8];
  int i,j,ierr,isBin;
  short shval;
  float fval;
  double dval;


  if( (fp = fopen( grdfile, "rb" )) == NULL ) return Err.post(Err.msgOpenFile(grdfile));
  memset( dsaa_label, 0, 5 );
  ierr = fread( dsaa_label, 4, 1, fp );
  if( !strcmp( dsaa_label,"DSAA" ) )
    isBin = 0;
  else if( !strcmp( dsaa_label,"DSBB" ) )
    isBin = 1;
  else {
    fclose(fp);
    return Err.post(Err.msgReadFile( grdfile, 1, "not GRD-file" ));
  }
  fclose(fp);


  if( isBin ) {
    fp = fopen( grdfile, "rb" );
    ierr = fread( dsaa_label, 4, 1, fp );
    ierr = fread( &shval, sizeof(short), 1, fp ); nx = shval;
    ierr = fread( &shval, sizeof(short), 1, fp ); ny = shval;
  }
  else {
    fp = fopen( grdfile, "rt" );
    ierr = fscanf( fp, "%s", dsaa_label );
    ierr = fscanf( fp, " %d %d ", &nx, &ny );
  }

  nnod = nx*ny;

  if( val ) delete [] val;
  val = new double[nnod];
  for( int l=0; l<nnod; l++ ) val[l]=noval;

  if( isBin ) {
    ierr = fread( &xmin, sizeof(double), 1, fp ); ierr = fread( &xmax, sizeof(double), 1, fp );
    ierr = fread( &ymin, sizeof(double), 1, fp ); ierr = fread( &ymax, sizeof(double), 1, fp );
    ierr = fread( &dval, sizeof(double), 1, fp ); ierr = fread( &dval, sizeof(double), 1, fp ); 
  }
  else {
    ierr = fscanf( fp, " %lf %lf ", &xmin, &xmax );
    ierr = fscanf( fp, " %lf %lf ", &ymin, &ymax );
    ierr = fscanf( fp, " %*s %*s " );  
  }

  dx = (xmax-xmin)/(nx-1);
  dy = (ymax-ymin)/(ny-1);

  for( j=0; j<ny; j++ ) {
    for( i=0; i<nx; i++ ) {

      if( isBin )
        ierr = fread( &fval, sizeof(float), 1, fp );
      else
        ierr = fscanf( fp, " %f ", &fval );

      val[idx(i,j)] = (double)fval;
    }
  }

  fclose( fp );

  return 0;
}

int GRDUTIL::readXYZ( const char *fname )
{
  #define MAXRECLEN 254
  FILE *fp;
  char record[254];
  int i,j,line;
  int ordering,ydirection;
  int nxny;
  double x0,x,y0,y,z;

  if( (fp = fopen( fname, "rt" )) == NULL )
    return utlPostError( utlErrMsgOpenFile(fname) );

  rewind( fp );
  for( line=0,nxny=0,xmin=ymin=RealMax,xmax=ymax=-RealMax; utlReadNextDataRecord( fp, record, &line ) != EOF; ) {
    if( sscanf( record, "%lf %lf %lf", &x, &y, &z ) != 3 ) return Err.post(Err.msgReadFile( fname, line, "X Y Z" ));
    nxny++;
    if( x < xmin ) xmin = x;
    if( x > xmax ) xmax = x;
    if( y < ymin ) ymin = y;
    if( y > ymax ) ymax = y;
  }

   rewind( fp );
  line = 0;
  utlReadNextDataRecord( fp, record, &line );
  sscanf( record, "%lf %lf %lf", &x0, &y0, &z );
  utlReadNextDataRecord( fp, record, &line );
  sscanf( record, "%lf %lf %lf", &x, &y, &z );
  if( x0==x && y0!=y )
    ordering = COLWISE;
  else if( x0!=x && y0==y )
    ordering = ROWWISE;
  else
    return Err.post(Err.msgReadFile( fname, line, "Cannot recognise data ordering" ));

   rewind( fp );
  line = 0;
  utlReadNextDataRecord( fp, record, &line );
  sscanf( record, "%lf %lf %lf", &x0, &y0, &z );
  if( ordering == ROWWISE ) {
    nx = 1;
    while( utlReadNextDataRecord( fp, record, &line ) != EOF ) {
      sscanf( record, "%lf %lf %lf", &x, &y, &z );
      if( y != y0 ) break;
      nx++;
    }
    ny = nxny / nx;
    if( y > y0 )
      ydirection = BOT2TOP;
    else
      ydirection = TOP2BOT;

  }
  else if( ordering == COLWISE ) {
    ny = 1;
    while( utlReadNextDataRecord( fp, record, &line ) != EOF ) {
      sscanf( record, "%lf %lf %lf", &x, &y, &z );
      if( x != x0 ) break;
      if( ny == 1 ) {
        if( y > y0 )
          ydirection = BOT2TOP;
        else
          ydirection = TOP2BOT;
      }
      ny++;
    }
    nx = nxny / ny;
  }

  if( nx*ny != nxny )
    return Err.post( "GRDUTIL::readXYZ -> nx*ny != nxny" );

   dx = (xmax-xmin)/(nx-1);
  dy = (ymax-ymin)/(ny-1);

  nnod = nx*ny;

  if( val ) delete [] val;
  val = new double[nnod];
  for( int l=0; l<nnod; l++ ) val[l]=noval;

  rewind( fp );
  line = 0;
  if( ordering == ROWWISE ) {
    if( ydirection == BOT2TOP ) {
      for( j=0; j<ny; j++ ) {
        for( i=0; i<nx; i++ ) {
          utlReadNextDataRecord( fp, record, &line );
          sscanf( record, "%*s %*s %lf", &val[idx(i,j)] );
        }
      }
    }
    else if( ydirection == TOP2BOT ) {
      for( j=ny-1; j>=0; j-- ) {
        for( i=0; i<nx; i++ ) {
          utlReadNextDataRecord( fp, record, &line );
          sscanf( record, "%*s %*s %lf", &val[idx(i,j)] );
        }
      }
    }
  }
  else if( ordering == COLWISE ) {
    if( ydirection == BOT2TOP ) {
      for( i=0; i<nx; i++ ) {
        for( j=0; j<ny; j++ ) {
          utlReadNextDataRecord( fp, record, &line );
          sscanf( record, "%*s %*s %lf", &val[idx(i,j)] );
        }
      }
    }
    else if( ydirection == TOP2BOT ) {
      for( i=0; i<nx; i++ ) {
        for( j=ny-1; j>=0; j-- ) {
          utlReadNextDataRecord( fp, record, &line );
          sscanf( record, "%*s %*s %lf", &val[idx(i,j)] );
        }
      }
    }
  }


  fclose( fp );

  return 0;
}
int GRDUTIL::readRasterStream( FILE *fp, int ordering, int ydirection )
{
  int i,j;
  double x0,x,y0,y,z;


  if( ordering == ROWWISE ) {
    if( ydirection == BOT2TOP ) {
      for( j=0; j<ny; j++ )
        for( i=0; i<nx; i++ )
          if( fscanf( fp, "%lf", &val[idx(i,j)] ) != 1 ) return Err.post("Unexpected EOF");
    }
    else if( ydirection == TOP2BOT ) {
      for( j=ny-1; j>=0; j-- )
        for( i=0; i<nx; i++ )
          if( fscanf( fp, "%lf", &val[idx(i,j)] ) != 1 ) return Err.post("Unexpected EOF");
    }
    else return Err.post("Unexpected direction");
  }
  else if( ordering == COLWISE ) {
    if( ydirection == BOT2TOP ) {
      for( i=0; i<nx; i++ )
        for( j=0; j<ny; j++ )
          if( fscanf( fp, "%lf", &val[idx(i,j)] ) != 1 ) return Err.post("Unexpected EOF");
    }
    else if( ydirection == TOP2BOT ) {
      for( i=0; i<nx; i++ )
        for( j=ny-1; j>=0; j-- )
          if( fscanf( fp, "%lf", &val[idx(i,j)] ) != 1 ) return Err.post("Unexpected EOF");
    }
    else return Err.post("Unexpected direction");
  }
  else return Err.post("Unexpected Ordering");

  return 0;
}

GRDUTIL* GRDUTIL::extract( int i1, int i2, int j1, int j2 )
{
  GRDUTIL *grd;

  if( i1 < 0 || i2 > (nx-1) || j1 < 0 || j2 > (ny-1) ) { Err.post("subGrid: bad ranges"); return NULL; }

  grd = new GRDUTIL( getX(i1,0), getX(i2,0), (i2-i1+1), getY(0,j1), getY(0,j2), (j2-j1+1) );

  for( int i=i1; i<=i2; i++ ) {
    for( int j=j1; j<=j2; j++ ) {
      (*grd).setVal( getVal(i,j), (i-i1), (j-j1) );
    }
  }

  return grd;
}

GRDUTIL* GRDUTIL::extract( double x1, double x2, double y1, double y2 )
{
  int i1,i2,j1,j2;

  i1 = (int)((x1-xmin)/dx); if( i1 < 0 ) i1 = 0;
  i2 = (int)((x2-xmin)/dx); if( i2 > (nx-1) ) i2 = nx-1;
  j1 = (int)((y1-ymin)/dy); if( j1 < 0 ) j1 = 0;
  j2 = (int)((y2-ymin)/dy); if( j2 > (ny-1) ) j2 = ny-1;

  return extract( i1, i2, j1, j2 );
}
void GRDUTIL::reset()
{
  xmin = xmax = dx = 0.;
  nx = 0;
  ymin = ymax = dy = 0.;
  ny = 0;

  nnod = 0;
  noval = DEFAULT_NOVAL;

  if( val != NULL ) delete [] val;
  val = NULL;
}

void GRDUTIL::resetVal()
{
  for( int l=0; l<nnod; l++ ) val[l]=noval;
}

GRDUTIL& GRDUTIL::operator= ( const GRDUTIL& grd )
{
  xmin=grd.xmin; xmax=grd.xmax; dx=grd.dx; nx=grd.nx;
  ymin=grd.ymin; ymax=grd.ymax; dy=grd.dy; ny=grd.ny;
  nnod=grd.nnod;

  if( val ) delete [] val;
  val = new double[nnod];
  memcpy( val, grd.val, nnod*sizeof(double) );

  return *this;
}
GRDUTIL& GRDUTIL::operator*= ( double coeff )
{

  for( int l=0; l<nnod; l++ )
    val[l] *= coeff;

  return *this;
}

GRDUTIL& GRDUTIL::operator+= ( GRDUTIL& grd )
{
  int l;
  int i,j,i1,i2,j1,j2;


  if( xmin==grd.xmin && xmax==grd.xmax && nx==grd.nx && ymin==grd.ymin && ymax==grd.ymax && ny==grd.ny ) {

    for( l=0; l<nnod; l++ )
      val[l] += grd.val[l];
  }
  else { 

    i1 = (int)((grd.xmin-xmin)/dx);
    if( i1 < 0 ) i1 = 0;
    i2 = (int)((grd.xmax-xmin)/dx) + 1;
    if( i2 > nx-1 ) i2 = nx-1;
    j1 = (int)((grd.ymin-ymin)/dy);
    if( j1 < 0 ) j1 = 0;
    j2 = (int)((grd.ymax-ymin)/dy) + 1;
    if( j2 > ny-1 ) j2 = ny-1;
    for( i=i1; i<=i2; i++ ) {
      for( j=j1; j<=j2; j++ ) {
        val[idx(i,j)] += grd.getVal( getX(i,j), getY(i,j) );
      }
    }
  }

  return *this;
}
double& GRDUTIL::operator() ( int i, int j )
{
  return( val[idx(i,j)] );
}

double& GRDUTIL::operator() ( int l )
{
  return( val[l] );
}

void GRDUTIL::getIJ( int idx, int& i, int& j )
{

  i = idx/ny;
  j = idx - i*ny;

}

int GRDUTIL::getIJ( double x, double y, int& i, int& j )
{
  i = (int)( (x-xmin)/(xmax-xmin)*nx );
  if( i<0 || i>(nx-1) ) return -1;
  j = (int)( (y-ymin)/(ymax-ymin)*ny );
  if( j<0 || j>(ny-1) ) return -1;

  return 0;
}
int GRDUTIL::idx( int i, int j )
{

 
  return( int( (int)ny*i + j ) );

}

double GRDUTIL::getVal( int i, int j )
{

  return( val[idx(i,j)] );

}

double GRDUTIL::getVal( int idx )
{

  return( val[idx] );

}

void GRDUTIL::setVal( double value, int i, int j )
{

  val[idx(i,j)] = value;

}
void GRDUTIL::setVal( double value, int idx )
{

  val[idx] = value;

}

double GRDUTIL::getMaxVal()
{
  int l;
  double vmax;


  for( vmax=-RealMax, l=0; l<nnod; l++ )
    if( val[l] > vmax )
      vmax = val[l];

  return vmax;
}

double GRDUTIL::getMaxVal( int& i, int& j )
{
  int l,lmax;
  double vmax;


  for( lmax=0,vmax=-RealMax, l=0; l<nnod; l++ ) {
    if( val[l] > vmax ) {
      vmax = val[l];
      lmax = l;
    }
  }

  getIJ( lmax, i, j );

  return vmax;
}

double GRDUTIL::getMinVal()
{
  int l;
  double vmin;


  for( vmin=RealMax, l=0; l<nnod; l++ )
    if( val[l] < vmin )
      vmin = val[l];

  return vmin;
}

double GRDUTIL::getMinVal( int& i, int& j )
{
  int l,lmin;
  double vmin;


  for( lmin=0,vmin=RealMax, l=0; l<nnod; l++ ) {
    if( val[l] < vmin ) {
      vmin = val[l];
      lmin = l;
    }
  }

  getIJ( lmin, i, j );

  return vmin;
}

double GRDUTIL::getMaxAbsVal()
{
  int l;
  double vmax;


  for( vmax=-RealMax, l=0; l<nnod; l++ )
    if( fabs(val[l]) > vmax )
      vmax = fabs(val[l]);

  return vmax;
}

double GRDUTIL::getMaxAbsVal( int& i, int& j )
{
  int l,lmax;
  double vmax;


  for( lmax=0,vmax=-RealMax, l=0; l<nnod; l++ ) {
    if( fabs(val[l]) > vmax ) {
      vmax = fabs(val[l]);
      lmax = l;
    }
  }

  getIJ( lmax, i, j );

  return vmax;
}

double GRDUTIL::getMaxAbsValBnd()
{
  int i,j;
  double vmax=-RealMax;

  for( i=0; i<nx; i++ ) {
    if( fabs(val[idx(i,0)]) > vmax ) vmax = fabs(val[idx(i,0)]);
    if( fabs(val[idx(i,ny-1)]) > vmax ) vmax = fabs(val[idx(i,ny-1)]);
  }

  for( j=0; j<ny; j++ ) {
    if( fabs(val[idx(0,j)]) > vmax ) vmax = fabs(val[idx(0,j)]);
    if( fabs(val[idx(nx-1,j)]) > vmax ) vmax = fabs(val[idx(nx-1,j)]);
  }

  return vmax;
}

double GRDUTIL::getVal( double x, double y )
{
  int i0,j0;
  double fx,fy,val_l,val_r,result;


  i0 = (int)((x-xmin)/dx);
  if( i0<0 || (i0+1)>=nx ) return(noval);
  j0 = (int)((y-ymin)/dy);
  if( j0<0 || (j0+1)>=ny ) return(noval);

  fx = (x - (xmin+dx*i0))/dx;
  fy = (y - (ymin+dy*j0))/dy;

  val_l = (1-fy)*getVal(i0,j0) + fy*getVal(i0,j0+1);
  val_r = (1-fy)*getVal(i0+1,j0) + fy*getVal(i0+1,j0+1);

  result = (1-fx)*val_l + fx*val_r;

  return( result );
}

double GRDUTIL::getX( int i, int j )
{
  return( xmin + i*dx );
}

double GRDUTIL::getX( int idx )
{
  int i,j;

  getIJ( idx, i, j );
  return( xmin + i*dx );
}
double GRDUTIL::getY( int i, int j )
{
  return( ymin + j*dy );
}

double GRDUTIL::getY( int idx )
{
  int i,j;

  getIJ( idx, i, j );
  return( ymin + j*dy );
}

int GRDUTIL::isSameShapeTo( GRDUTIL& grd )
{
  if( xmin != grd.xmin ) return 0;
  if( xmax != grd.xmax ) return 0;
  if( nx != grd.nx ) return 0;
  if( ymin != grd.ymin ) return 0;
  if( ymax != grd.ymax ) return 0;
  if( ny != grd.ny ) return 0;

  return 1;
}

int GRDUTIL::getIntersectionRegion( const GRDUTIL &grd, int& imin, int& imax, int& jmin, int& jmax )
{

  if( xmin < grd.xmin ) {
    if( xmax < grd.xmin ) return -1;
    imin = (int)((grd.xmin-xmin)/dx);
  }
  else if( xmin <= grd.xmax ) {
    imin = 0;
  }
  else
    return -1;

  if( xmax < grd.xmin )
    return -1;
  else if( xmax <= grd.xmax ) {
    imax = nx-1;
  }
  else {
    imax = (int)((grd.xmax-xmin)/dx);
  }

  if( ymin < grd.ymin ) {
    if( ymax < grd.ymin ) return -1;
    jmin = (int)((grd.ymin-ymin)/dy);
  }
  else if( ymin <= grd.ymax ) {
    jmin = 0;
  }
  else
    return -1;

  if( ymax < grd.ymin )
    return -1;
  else if( ymax <= grd.ymax ) {
    jmax = ny-1;
  }
  else {
    jmax = (int)((grd.ymax-ymin)/dy);
  }

  return 0;
}

int GRDUTIL::interpolateFrom( GRDUTIL &grd, int resetValues )
{
  int ierr,imin,imax,jmin,jmax;
  double value;

  if( resetValues ) resetVal();

  ierr = getIntersectionRegion( grd, imin, imax, jmin, jmax );
  if(ierr) return 0;

  for( int i=imin; i<=imax; i++ ) {
    for( int j=jmin; j<=jmax; j++ ) {
      value = grd.getVal( getX(i,j), getY(i,j) );
      if( value != grd.noval ) val[idx(i,j)] = value;
    }
  }

  return( 0 );
}

int GRDUTIL::getNearestIdx( double x0, double y0 )
{
  int i0,j0;
  int l,lmin;
  double dist2,dist2min;


  if( x0<xmin || x0>xmax || y0<ymin || y0>ymax )
    return -1;

  i0 = (int)((x0-xmin)/dx);
  j0 = (int)((y0-ymin)/dy);

  l = idx(i0,j0);
  dist2min = (x0-getX(l))*(x0-getX(l)) + (y0-getY(l))*(y0-getY(l));
  lmin = l;

  l = idx(i0+1,j0);
  dist2 = (x0-getX(l))*(x0-getX(l)) + (y0-getY(l))*(y0-getY(l));
  if( dist2 < dist2min ) { dist2min = dist2; lmin = l; }

  l = idx(i0,j0+1);
  dist2 = (x0-getX(l))*(x0-getX(l)) + (y0-getY(l))*(y0-getY(l));
  if( dist2 < dist2min ) { dist2min = dist2; lmin = l; }

  l = idx(i0+1,j0+1);
  dist2 = (x0-getX(l))*(x0-getX(l)) + (y0-getY(l))*(y0-getY(l));
  if( dist2 < dist2min ) { dist2min = dist2; lmin = l; }

  return lmin;
}

int GRDUTIL::getNearestIdx( double x0, double y0, double rangemin, double rangemax )
{
  int i,j,i0,j0,rad;
  int l,lmin;
  double dist2,dist2min;


  lmin = getNearestIdx( x0, y0 );
  if( lmin == -1 )
    return lmin;

  if( val[lmin] >= rangemin && val[lmin] <= rangemax )
    return lmin;

  getIJ( lmin, i0,j0 );

  for( lmin=-1,rad=1; rad<nx && rad<ny; rad++ ) {

    dist2min = RealMax;

    for( i=i0-rad; i<=i0+rad; i++ )
      for( j=j0-rad; j<=j0+rad; j++ ) {
        if( iabs(i-i0) != rad && iabs(j-j0) != rad ) continue;
        if( i<0 || i>nx-1 || j<0 || j>ny-1 ) continue;
        l = idx(i,j);
        if( val[l] < rangemin || val[l] > rangemax ) continue;
        dist2 = (x0-getX(l))*(x0-getX(l)) + (y0-getY(l))*(y0-getY(l));
        if( dist2 < dist2min ) { dist2min = dist2; lmin = l; }
      }

    if( lmin > 0 ) break;
  }

  return lmin;
}

void GRDUTIL::smooth( int radius )
{
  int i,j,k,ik,jk;
  int l;
  double sumwt;
  double wt[10] = { 1.0, 0.5, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25 };

  if( radius == 0 ) return;

  GRDUTIL tmp( *this );

  for( i=0; i<nx; i++ ) {
    for( j=0; j<ny; j++ ) {

      l = idx(i,j);

      for( sumwt=tmp(l)=0, k=0; k<=radius; k++ ) {
        for( ik=i-k; ik<=i+k; ik++ ) {
          if( ik<0 || ik>=nx ) continue;
          for( jk=j-k; jk<=j+k; jk++ ) {
            if( jk<0 || jk>=ny ) continue;
            if( iabs(ik-i)==k || iabs(jk-j)==k ) {
              tmp(l) += wt[k] * ((*this)(ik,jk));
              sumwt += wt[k];
            }
          }
        }
      }
      tmp(l) /= sumwt;
    }
  }

  *this = tmp;
}

int GRDUTIL::writeGRD( const char *fname )
{
  FILE *fp;
  int i,j,cnt;


  fp = fopen( fname, "wt" );

  fprintf( fp, "DSAA\n" );
  fprintf( fp, "%d %d\n", nx, ny );
  fprintf( fp, "%f %f\n", xmin, xmax );
  fprintf( fp, "%f %f\n", ymin, ymax );
  fprintf( fp, "%f %f\n", getMinVal(), getMaxVal() );

  for( cnt=0, j=0; j<ny; j++ ) {
    for( i=0; i<nx; i++ ) {
      cnt++;
      fprintf( fp, " %g", val[idx(i,j)] );
      if( cnt == 10 ) {
        fprintf( fp, "\n" );
        cnt = 0;
      }
    }
    fprintf( fp, "\n\n" );
    cnt = 0;
  }

  fclose( fp );

  return 0;
}

int GRDUTIL::writeGRDbin( const char *fname )
{
  FILE *fp;
  short i2buf;
  float r4buf;
  double r8buf;
  int i,j;


  fp = fopen( fname, "wb" );

  fwrite( "DSBB", 4,1, fp );
  i2buf = (short)nx; fwrite( &i2buf, sizeof(short), 1, fp );
  i2buf = (short)ny; fwrite( &i2buf, sizeof(short), 1, fp );
  fwrite( &xmin, sizeof(double), 1, fp );
  fwrite( &xmax, sizeof(double), 1, fp );
  fwrite( &ymin, sizeof(double), 1, fp );
  fwrite( &ymax, sizeof(double), 1, fp );
  r8buf = (double)getMinVal(); fwrite( &r8buf, sizeof(double), 1, fp );
  r8buf = (double)getMaxVal(); fwrite( &r8buf, sizeof(double), 1, fp );

  for( j=0; j<ny; j++ ) {
    for( i=0; i<nx; i++ ) {
      r4buf = (float)val[idx(i,j)];
      fwrite( &r4buf, sizeof(float), 1, fp );
    }
  }

  fclose( fp );

  return 0;
}
int GRDUTIL::writeXYZ( const char *xyzfile )
{
  FILE *fp;
  int i,j,l;

  fp = fopen( xyzfile, "wt" );

  for( j=0; j<ny; j++ ) {
    for( i=0; i<nx; i++ ) {
      l = idx(i,j);
      fprintf( fp, " %g %g %g\n", getX(l), getY(l), val[l] );
    }
  }

  fclose( fp );

  return 0;
}
