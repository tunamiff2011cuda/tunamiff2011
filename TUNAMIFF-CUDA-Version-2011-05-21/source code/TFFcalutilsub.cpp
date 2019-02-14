/*
 *      Tohoku University's Numerical-Analysis Model
 *          for Investigation of tsunami
 *           Far-field Tsunami version 
 *
 */
//====================================================

#include <stdio.h>         
#include <stdarg.h>        
#include <stdlib.h>        
#include <string.h>        
#include <math.h>          
#include <ctype.h>         
#include <time.h>          
#include "GRDUTIL.h"        
#include "DEFORMINITIA.h"   
#include "UTILMDL.h"       
#include "SPHERECAL.h"     
#include "OFAULTPARA.h"    
#include "TUNAMI_FF.h"     
//====================================================

int NLon,NLat;
double LonMin,LonMax,LatMin,LatMax;
double DLon,DLat;                 
double Dx,Dy;                     
float *R6;
float *C1;
float *C2;
float *C3;
float *C4;
struct TFPARAMS MODELVARI;
double fun_Chinnery( double (*fun)(double ksi, double eta), double x, double y );
double f_ssUx(double ksi, double eta);
double f_ssUy(double ksi, double eta);
double f_ssUz(double ksi, double eta);
double f_dsUx(double ksi, double eta);
double f_dsUy(double ksi, double eta);
double f_dsUz(double ksi, double eta);
double fun_R(double ksi, double eta);
double fun_X(double ksi, double eta);
double fun_yp(double ksi, double eta);
double fun_dp(double ksi, double eta);
double fun_I1(double ksi, double eta);
double fun_I2(double ksi, double eta);
double fun_I3(double ksi, double eta);
double fun_I4(double ksi, double eta);
double fun_I5(double ksi, double eta);

static double sdip;
static double cdip;
static double p;
static double q;
static double width;
static double length;
static double elast;
int Imin;
int Imax;
int Jmin;
int Jmax;
static int MaxGAGUEs;
int NGAGUEs;
static char **idGAGUE;
long *idxGAGUE;
static int *flagRunupGAGUE;
static int NtGAGUE;
static int *timeGAGUE;
static float **zoutGAGUE;
static char *IndexFile;
static int NrecSpatialWaveOutput;

#define SRC_GRD 1
#define SRC_FLT 2
#define ERRORFILE "errorlog.lg"
#define PI_RC 3.14159265358979
#define DISPLMAX 1000

/*
 *Implementation flow
 *
 **BathymetryInputfile - Bathymetry
 **fileSource - Source: Okada faults or Surfer grid
 **KEtotalTime - Time steps of Total simulation - end of computation, [sec]
 **oopsName - Simulation Run
 **coriolis - Use Coriolis force
 **MAXOUT - Periodic dumping of mariograms and cumulative SpatialWave-plots (wavemax, arrival times), [sec]
 **outProgress - Reporting simulation progress, [sec model time]
 **KDspatialOut - Time step lengths to output spatial wave profiles, [sec model time]
 **dmin - minimal calculation depth, [m]
 **DT - time step length in second
 **zout0TRel - Initial uplift: relative threshold
 **zout0TAbs - Initial uplift: absolute threshold, [m]
 **zoutAT - Threshold for SpatialWave arrival time (0 - do not calculate), [m]
 **zoutCT - Threshold for clipping of expanding computational area, [m]
 **zoutZT - Threshold for resetting the small zout (keep expanding area from unnesessary growing), [m]
 **zoutTT - Threshold for transparency, [m]
 **fileGAGUEs - Points Of Interest (GAGUEs) input file
 **GAGUEDIMX - GAGUE fitting: max search distance, [km]
 **GAGUEDEMN - GAGUE fitting: min depth, [m]
 **GAGUEDEMX - GAGUE fitting: max depth, [m]
 **GAGUERT - report of GAGUE loading
 **KCgagueOUT - Time step lengths in outputting the time history of water level, [sec]
 **gpu - GPU computation enable or disable
 *
 *
 */
int MODELPARAINITIA( int argc, char **argv )
{
  int argn,ierr;
    
	MODELVARI.BathymetryInputfile = strdup("Ogrid.grd"); 
 	MODELVARI.fileSource = strdup("faultdeform.grd"); 
 	MODELVARI.KEtotalTime = 600; 
 	MODELVARI.KEtotalTime *= 60; 
	MODELVARI.oopsName = strdup( "OOPS" );
	MODELVARI.coriolis = 1;
	MODELVARI.MAXOUT = 0;
	MODELVARI.outProgress = 60;
	MODELVARI.KDspatialOut = 60;
	MODELVARI.dmin = 10.;
	MODELVARI.DT = 0;  
	MODELVARI.zout0TRel = 0.01;
	MODELVARI.zout0TAbs = 0.0;
	MODELVARI.zoutAT = 0.001;
	MODELVARI.zoutCT = 1.e-4;
	MODELVARI.zoutZT = 1.e-5;
	MODELVARI.zoutTT = 0.0;
	MODELVARI.fileGAGUEs = NULL;
	MODELVARI.GAGUEDIMX = 10.0; 
	MODELVARI.GAGUEDIMX *= 1000.;
	MODELVARI.GAGUEDEMN = 1.0;
	MODELVARI.GAGUEDEMX = 10000.0;
	MODELVARI.GAGUERT = 0;
	MODELVARI.KCgagueOUT = 60;
	MODELVARI.gpu = true;
	//MODELVARI.gpu = false;
	MODELVARI.adjustZtop = false;
	MODELVARI.verbose = true;
   return 0;
}

/*
 *
 *Implementation flow
 *
 *Class Message
 *Class Error Message
 *
 */
int iscomment( int );

cMsg Msg;
cErrMsg Err;

cMsg::cMsg()
{
  enabled = 1;
  setchannel(MSG_OUTSCRN);
  setfilename( "defaultlog.lg" );
}

cMsg::~cMsg()
{

}

void cMsg::enable()
{
  enabled = 1;
}

void cMsg::disable()
{
  enabled = 0;
}

void cMsg::setfilename( const char* newfilename )
{
  sprintf( fname, "%s", newfilename );
}

void cMsg::setchannel( int newchannel )
{
  channel = newchannel;
}

int cMsg::print( const char* fmt, ... )
{
  if(!enabled) return 0;

  if( channel & MSG_OUTSCRN ) {
    va_list arglst;
    va_start(arglst, fmt);
    vprintf( fmt, arglst );
    va_end(arglst);
  }

  if( channel & MSG_OUTFILE ) {
    va_list arglst;
    va_start(arglst, fmt);
    FILE *fp = fopen( fname, "at" );
    vfprintf( fp, fmt, arglst );
    fclose( fp );
    va_end(arglst);
  }

  return 0;
}

cErrMsg::cErrMsg()
{
  enabled = 1;
  setfilename( "errorlog.lg" );
  setchannel(MSG_OUTSCRN|MSG_OUTFILE);
}

cErrMsg::~cErrMsg()
{

}

int cErrMsg::post( const char* fmt, ... )
{
  if(!enabled) return 0;

  if( channel & MSG_OUTSCRN ) {
    va_list arglst;
    va_start(arglst, fmt);
    vprintf( fmt, arglst );
    va_end(arglst);
    printf( "\n" );
  }

  if( channel & MSG_OUTFILE ) {
    va_list arglst;
    va_start(arglst, fmt);
    FILE *fp = fopen( fname, "at" );
    vfprintf( fp, fmt, arglst );
    fprintf( fp, "\n" );
    fclose( fp );
    va_end(arglst);
  }

  return -1;
}

char* cErrMsg::msgAllocateMem( void )
{
  char msg[256];

  sprintf( msg, "Error allocating memory" );

  return strdup(msg);
}

char* cErrMsg::msgOpenFile( const char *fname )
{
  char msg[256];

  sprintf( msg, "Cannot open file %s", fname );

  return strdup(msg);
}

char* cErrMsg::msgReadFile( const char *fname, int line, const char *expected )
{
  char msg[256];

  sprintf( msg, "Error reading file %s, line number %-d. Expected: %s", fname, line, expected );

  return strdup(msg);
}

char *utlCurrentTime()
{
  time_t timer;
  char *cp;

  timer = time(NULL);
  cp = asctime(localtime(&timer));
  cp[strlen(cp)-1] = '\0';

  return cp;
}


int utlPostError( const char *message )
{
  FILE *fp = fopen( ERRORFILE, "at" );
  fprintf( fp, "%s", utlCurrentTime() );
  fprintf( fp, " -> %s\n", message );
  fclose( fp );

  return -1;
}


char *utlErrMsgMemory()
{
  return strdup("Cannot allocate memory");
}


char *utlErrMsgOpenFile( const char *fname )
{
  char msg[256];

  sprintf( msg, "Cannot open file %s", fname );

  return strdup(msg);
}


char *utlErrMsgEndOfFile( const char *fname )
{
  char msg[256];

  sprintf( msg, "Unexpected end of file: %s", fname );

  return strdup(msg);
}


char *utlErrMsgReadFile( const char *fname, int line, const char *expected )
{
  char msg[256];

  sprintf( msg, "Error reading file %s, line number %-d. Expected: %s", fname, line, expected );

  return strdup(msg);
}

int utlCheckCommandLineOption( int argc, char **argv, const char *option, int letters_to_compare )
{
  int k;


  for( k=1; k<argc; k++ ) {
    if( argv[k][0] != '-' && argv[k][0] != '/' ) continue;
    if( !strncmp( argv[k]+1, option, letters_to_compare ) ) break;
  }

  if( k == argc )
    return 0;
  else
    return k;
}


int utlStringReadOption( char *record, char *option_name, char *contents )
{
  int found, length;
  char *cp, *cpe, buf[64];


  cp = record;
  found = 0;
  contents[0] = '\0';

  while( *cp != '\0' ) {

    cp = strchr( cp, '[' );
    if( cp == NULL ) break;

    cpe = strchr( cp+1, ']' );
    if( cpe == NULL ) break;

    length = cpe - cp - 1;
    strncpy( buf, cp+1, length );
    buf[length] = '\0';

    if( !strcmp( buf, option_name ) ) {   

      cp = strchr( cpe+1, '=' );
      if( cp == NULL ) break;
      while( isspace( *(++cp) ) ) ;

      if( *cp == '\0' )  
        ;
      else if( *cp == '[' )  
        ;
      else if( *cp == '"' ) {  
        cpe = strchr( cp+1, '"' );
        if( cpe == NULL ) break;
        length = cpe - cp - 1;
        strncpy( contents, cp+1, length );
        contents[length] = '\0';
      }
      else 
        sscanf( cp, "%s", contents );

      found = 1;
      break;

    }  

    cp++;
  }

  if( found )
    return 0;
  else
    return -1;
}


int utlReadSequenceOfInt( char *line, int *value )
{
  #define MAX_Ints 100
  int N=0;
  int itmp;
  char *cp=line;

  while( 1 )
    {
    while( !isdigit(*cp) && *cp != '\0' && *cp != ';' ) cp++;
    if( *cp == '\0' || *cp == ';' )  return N;
    if( sscanf( cp, "%d", &itmp ) )
      value[N++] = itmp;
    if( N == MAX_Ints ) return -N;
    while( isdigit(*cp) ) cp++;
    }
}


int utlReadSequenceOfDouble( char *line, double *value )
{
  #define MAX_Doubles 100
  int N=0;
  double dtmp;
  char *cp=line;

  while( 1 )
    {
    while( isspace( *cp ) ) cp++;
    if( *cp == '\0' || *cp == ';' )  return N;
    if( sscanf( cp, "%lf", &dtmp ) )
      value[N++] = dtmp;
    if( N == MAX_Doubles ) return -N;
    while( !isspace( *cp ) && *cp != '\0' ) cp++;
    }
}


int utlPickUpWordFromString( char *string, int n, char *word )
{
  char *cp;
  int i;

  for( cp = string, i = 1; 1; i++ )
    {
    while( isspace( *cp ) )
      cp++;
    if( *cp == '\0' )
      {
      word[0] = '\0';
      return -1;
      }
    if( i == n )
      break;
    while( !isspace( *cp ) && *cp != '\0' )
      cp++;
    }

  sscanf( cp, "%s", word );
  return 0;
}


char *utlPickUpWordFromString( char *string, int n )
{
  char *cp,buf[64];
  int i;

  for( cp = string, i = 1; 1; i++ )
    {
    while( isspace( *cp ) )
      cp++;
    if( *cp == '\0' )
      return NULL;
    if( i == n )
      break;
    while( !isspace( *cp ) && *cp != '\0' )
      cp++;
    }

  sscanf( cp, "%s", buf );
  return strdup(buf);
}


char *utlPickUpWordFromString( char *string, char *pos, char *word )
{
  char *cp;

  if( pos < string ) { word[0] = '\0'; return NULL; }
  if( pos > (string+strlen(string)) ) { word[0] = '\0'; return NULL; }

  for( cp = pos; isspace(*cp); cp++ ) ;
  if( *cp == '\0' ) { word[0] = '\0'; return NULL; }
  sscanf( cp, "%s", word );
  cp += strlen(word);

  return cp;
}


int utlCountWordsInString( char *line )
{
  int nwords=0;
  char *cp=line;

  while( 1 )
    {
    while( isspace( *cp ) ) cp++;
    if( *cp == '\0' || *cp == ';' )  return nwords;
    nwords++;
    while( !isspace( *cp ) && *cp != '\0' ) cp++;
    }

}

char *utlWords2String( int nwords, char **word )
{
  char *buf;
  int k,lenstr;


  for( lenstr=0,k=0; k < nwords; k++ )
    lenstr += strlen( word[k] );

  lenstr += (nwords + 1);  

  buf = new char[lenstr];

  memset( buf, 0, lenstr );

  for( k=0; k < nwords; k++ ) {
    if( k>0 ) strcat( buf, " " );
    strcat( buf, word[k] );
  }

  return buf;
}


int utlSubString( char *str, int p1, int p2, char *substr )
{
  char *cp;
  int k;


  if( p1<0 || p2<0 || p1>p2 || (unsigned)p2 > (strlen(str)-1) ) return -1;

  for( k=0,cp = &str[p1]; cp <= &str[p2]; cp++ )
    substr[k++] = *cp;
  substr[k] = '\0';

  return 0;
}


void utlPrintToString( char *string, int position, char *insert )
{
  char *cp;
  int i;

  for( i = 0; i < position; i++ )
    if( string[i] == '\0' )  string[i] = ' ';
  for( cp = insert, i = 0; *cp != '\0'; i++, cp++ )
    string[position+i] = *cp;

}


int iscomment( int ch )
{
  if( ch == ';' || ch == '!' )
    return 1;
  else
    return 0;
}


int utlFindNextInputField( FILE *fp )
{
  int ch;
  int lines = 0;

L1: while( isspace( (ch=fgetc(fp)) ) )
    if( ch == '\n' )  lines++;

  if( iscomment(ch) )
    {
    while( ( ch = fgetc( fp ) ) != '\n' && ch != EOF ) ;
    if( ch == '\n' )
      {
      lines++;
      goto L1;
      }
    else if( ch == EOF )
      return EOF;
    }
  else if( ch == EOF )
    return EOF;

  ungetc( ch, fp );

  return( lines );
}


int utlReadNextRecord( FILE *fp, char *record, int *line )
{
  int found = 0;
  char *cp, firstchar;

  while( !found && fgets( record, MaxFileRecordLength, fp ) )
    {
    for( cp = record; *cp == ' ' || *cp == '\t'; cp++ )
      ;
    if( *cp != '\n' && *cp != '\r' && *cp != ';' )
      found = 1;
    (*line)++;
    }

  if( !found )
    return EOF;
  else
    {
    firstchar = *cp;
    while( *cp != '\n' && *cp != '\r' && *cp != '\0' )
      cp++;
    *cp = '\0';
    return( firstchar );
    }
}


int utlReadNextDataRecord( FILE *fp, char *record, int *line )
{
  int found = 0;
  char *cp, firstchar;

  while( !found && fgets( record, MaxFileRecordLength, fp ) )
    {
    for( cp = record; *cp == ' ' || *cp == '\t'; cp++ )
      ;
    if( isdigit(*cp) || (*cp=='-' && isdigit(*(cp+1))) )
      found = 1;
    (*line)++;
    }

  if( !found )
    return EOF;
  else
    {
    firstchar = *cp;
    while( *cp != '\n' && *cp != '\r' && *cp != '\0' )
      cp++;
    *cp = '\0';
    return( firstchar );
    }
}


char *utlFileFindRecord( char *fname, char *pattern, char *record )
{
  char *cp;
  int line;
  FILE *fp;


  if( (fp = fopen(fname,"rt")) == NULL ) { utlPostError( utlErrMsgOpenFile(fname) ); return NULL; }

  line = 0;
  cp = NULL;
  while( utlReadNextRecord( fp, record, &line ) != EOF ) {

    if( (cp = strstr(record, pattern) ) == NULL ) continue;

    break;
  }

  fclose(fp);

  return cp;
}


int utlFileParceToString( FILE *fp, char *pattern )
{
  int ch,pos=0,ierr=1;


  while( (ch=fgetc(fp)) != EOF ) {

    if( ch != pattern[pos] ) {
      pos = 0;
      continue;
    }

    pos++;

    if( pattern[pos] == '\0' ) {
      ierr = 0;
      break;
    }

  }

  return ierr;
}


int utlGetNumberOfRecords( const char *fname )
{
  FILE *fp;
  int nrec,line=0;
  char record[1024];


  fp = fopen( fname, "rt" );
  if( fp == NULL ) return 0;

  nrec = 0;
  while( utlReadNextRecord( fp, record, &line ) != EOF )
    nrec++;

  fclose( fp );

  return nrec;
}


int utlGetNumberOfDataRecords( const char *fname )
{
  FILE *fp;
  int nrec,line=0;
  char record[1024];


  fp = fopen( fname, "rt" );
  if( fp == NULL ) return 0;

  nrec = 0;
  while( utlReadNextDataRecord( fp, record, &line ) != EOF )
    nrec++;

  fclose( fp );

  return nrec;
}


int utlFileReadOption( char *fname, char *option_name, char *contents )
{
  FILE *fp;
  int found,line;
  char record[1024];


  fp = fopen( fname, "rt" );
  if( !fp ) return 0;

  for( line=0, found=0; utlReadNextRecord( fp, record, &line ) != EOF; ) {
    if( utlStringReadOption( record, option_name, contents ) == 0 ) {
      found = 1;
      break;
    }
  }
  fclose( fp );

  if( found )
    return 0;
  else
    return -1;
}


char *utlFileChangeExtension( const char *fname, const char *newext )
{
  char buf[256], *cp;
  int len,extlen;

  len = strlen(fname);
  extlen = strlen(newext);
  memset( buf, 0, len+extlen+1 );
  strcpy( buf, fname );

  for( cp=&buf[len-1]; *cp!='.' && cp>=&buf[0]; cp-- ) ;

  if( *cp=='.' )
    sprintf( cp+1, "%s", newext );
  else
    sprintf( &buf[len], ".%s", newext );

  return strdup( buf );
}


char *utlFileAddExtension( const char *fname, const char *addext )
{
  char buf[256], *cp;
  int len,extlen;

  len = strlen(fname);
  extlen = strlen(addext);
  memset( buf, 0, len+extlen+1 );
  sprintf( buf, "%s.%s", fname, addext );

  return strdup( buf );
}


char *utlFileRemoveExtension( const char *fname )
{
  char buf[256], *cp;

  strcpy( buf, fname );
  for( cp=&buf[strlen(fname)-1]; *cp!='.' && cp>=&buf[0]; cp-- ) ;
  if( *cp=='.' ) *cp = '\0';

  return strdup( buf );
}


int utlFileRemoveExtension( char *newname, const char *fname )
{
  char *cp;

  strcpy( newname, fname );
  for( cp=&newname[strlen(newname)-1]; *cp!='.' && cp>=&newname[0]; cp-- ) ;
  if( *cp=='.' ) *cp = '\0';

  return 0;
}


int utlReadXYdata( char *fname, double **x, double **y )
{
  FILE *fp;
  char record[256];
  int n,line=0;
  double xv,yv;


  fp = fopen( fname, "rt" );
  if( fp == NULL )
    return utlPostError( utlErrMsgOpenFile(fname) );

  for( n=0; utlReadNextDataRecord( fp, record, &line ) != EOF; ) {
    if( sscanf( record, "%lf %lf", &xv, &yv ) != 2 )
      return utlPostError( utlErrMsgReadFile( fname, line, "X Y" ) );
    n++;
  }

  *x = new double [n];
  *y = new double [n];

  rewind( fp );
  line=0;
  for( int k=0; k<n; k++ ) {
    utlReadNextDataRecord( fp, record, &line );
    sscanf( record, "%lf %lf", &(*x)[k], &(*y)[k] );
  }

  return n;
}

double utlRandom( int& seed )
{

  seed = seed * 65539;
  if( seed < 0 ) seed = (seed + 2147483647) + 1;

  return( (double)seed * 0.4656613e-9 );
}


double utlRandom( double avg, double err )
{
  return( avg + ( (double)rand()/RAND_MAX-0.5 )*2*err );
}

double utlNormal( int& seed, double avg, double stdev )
{
  double r1, r2, w1, w2;

  r1 = utlRandom( seed );
  r2 = utlRandom( seed );
  w1 = log10(r1);
  w1 = sqrt(-2*w1);
  w2 = sin(2*PI_RC * r2);

  return ( avg + stdev*w1*w2 );
}


int utlPointInPolygon( int nvrtx, double *xvrtx, double *yvrtx, double x, double y )
{
  int l,ln,nsec;
  double yside;


  for( nsec=0, l=0; l<nvrtx; l++ ) {

    ln = l + 1;
    if( ln == nvrtx) ln = 0;

    if( (xvrtx[l]-x)*(xvrtx[ln]-x) >= 0 )
      continue;

    yside = yvrtx[l] + (yvrtx[ln]-yvrtx[l])/(xvrtx[ln]-xvrtx[l])*(x-xvrtx[l]);

    if( yside > y )
      nsec++;
  }

  if( (nsec/2)*2 != nsec )
    return 1;
  else
    return 0;

}


void utlTimeSplit( double ctime, int& nHour, int& nMin, int& nSec, int& nMsec )
{
  int fullSec;

  fullSec = (int)ctime;

  nMsec = (int)((ctime - fullSec)*1000 + 0.1);
  nHour = fullSec/3600;
  nMin = (fullSec - nHour*3600)/60;
  nSec = fullSec - nHour*3600 - nMin*60;

  return;
}


char *utlTimeSplitString( double ctime ) 
{
  int nHour,nMin,nSec,nMsec;
  static char buf[32];

  utlTimeSplit( ctime, nHour, nMin, nSec, nMsec );

  if( nMsec > 0 )
    sprintf( buf, "%2.2d:%2.2d:%2.2d.%3.3d", nHour,nMin,nSec,nMsec );
  else
    sprintf( buf, "%2.2d:%2.2d:%2.2d", nHour,nMin,nSec );

  return buf;
}
char *utlTimeSplitString( int ctime )
{
  int nHour,nMin,nSec,nMsec;
  static char buf[32];

  utlTimeSplit( (double)ctime, nHour, nMin, nSec, nMsec );

  sprintf( buf, "%2.2d:%2.2d:%2.2d", nHour,nMin,nSec );

  return buf;
}


int utlTimeHour( double ctime )
{
  int nHour,nMin,nSec,nMsec;

  utlTimeSplit( ctime, nHour, nMin, nSec, nMsec );

  return nHour;
}
int utlTimeHour( int ctime )
{
  int nHour,nMin,nSec,nMsec;

  utlTimeSplit( (double)ctime, nHour, nMin, nSec, nMsec );

  return nHour;
}


int utlTimeMin( double ctime )
{
  int nHour,nMin,nSec,nMsec;

  utlTimeSplit( ctime, nHour, nMin, nSec, nMsec );

  return nMin;
}
int utlTimeMin( int ctime )
{
  int nHour,nMin,nSec,nMsec;

  utlTimeSplit( (double)ctime, nHour, nMin, nSec, nMsec );

  return nMin;
}


int utlTimeSec( double ctime )
{
  int nHour,nMin,nSec,nMsec;

  utlTimeSplit( ctime, nHour, nMin, nSec, nMsec );

  return nSec;
}
int utlTimeSec( int ctime )
{
  int nHour,nMin,nSec,nMsec;

  utlTimeSplit( (double)ctime, nHour, nMin, nSec, nMsec );

  return nSec;
}
/*
 *Implementation flow
 *
 *Data input of water depth
 *Setting of parameters required in vectorized computation
 */
/*
FLOW of Implementation:
Check if bathymetry file is in ascii or binary format
 *allocate memory for GRIDNODE structure and for caching arrays
 *Make bathymetry from topography. Compute stable time step
 *Correct bathymetry for edge artefacts
 *Calculate caching grid parameters R1, R2, R3, R4, R5, R6 and C1, C2, C3, C4
*/

/*
VARIABLES
 *DLon, DLat, Dx, Dy
 *iTopo, iTime, iD, NLat, NLon, LatMin
 *MODELVARI
 */

/*
COEFFICIENTS
 *R1, R2, R3, R4, R5, R6 and
 *Coefficients used at open sea boundary  C1, C2, C3, C4 
*/

/*
CONSTANTS
 *gravity acceleration  GG 
 *Earth rotation period [1/sec]  Omega
 *Earth radius  RADIEAR 
 */
int tfRdepthParame()
{
  FILE *fp;
  char fileLabel[5];
  unsigned short shval;
  int ierr,isBin,i,j,m,k;
  float fval;
  double dval;

  CNode& Node = *gNode;
  if( (fp=fopen(MODELVARI.BathymetryInputfile,"rb")) == NULL ) return Err.post( Err.msgOpenFile(MODELVARI.BathymetryInputfile) );

  memset( fileLabel, 0, 5 );
  ierr = fread( fileLabel, 4, 1, fp );
  if( !strcmp( fileLabel,"DSAA" ) )
    isBin = 0;
  else if( !strcmp( fileLabel,"DSBB" ) )
    isBin = 1;
  else
    return Err.post( "%s: not GRD-file!", MODELVARI.BathymetryInputfile );

  fclose(fp);

  if( isBin ) {
    fp = fopen( MODELVARI.BathymetryInputfile, "rb" );
    ierr = fread( fileLabel, 4, 1, fp );
    ierr = fread( &shval, sizeof(unsigned short), 1, fp ); NLon = shval;
    ierr = fread( &shval, sizeof(unsigned short), 1, fp ); NLat = shval;
  }
  else {
    fp = fopen( MODELVARI.BathymetryInputfile, "rt" );
    ierr = fscanf( fp, "%s", fileLabel );
    ierr = fscanf( fp, " %d %d ", &NLon, &NLat );
  }

   if( Node.mallocMem() ) return Err.post( Err.msgAllocateMem() );

  if( isBin ) {
    ierr = fread( &LonMin, sizeof(double), 1, fp ); ierr = fread( &LonMax, sizeof(double), 1, fp );
    ierr = fread( &LatMin, sizeof(double), 1, fp ); ierr = fread( &LatMax, sizeof(double), 1, fp );
    ierr = fread( &dval, sizeof(double), 1, fp ); ierr = fread( &dval, sizeof(double), 1, fp ); 
  }
  else {
    ierr = fscanf( fp, " %lf %lf ", &LonMin, &LonMax );
    ierr = fscanf( fp, " %lf %lf ", &LatMin, &LatMax );
    ierr = fscanf( fp, " %*s %*s " );   
  }

  DLon = (LonMax - LonMin)/(NLon - 1);   
  DLat = (LatMax - LatMin)/(NLat - 1);

  Dx = RADIEAR * g2r( DLon );     
  Dy = RADIEAR * g2r( DLat );

  if( isBin ) {


      float *buf = new float[ NLat*NLon ];
	  ierr = fread( buf, sizeof(float), NLat*NLon, fp );

	  for( i=1; i<=NLon; i++ ) {
		for( j=1; j<=NLat; j++ ) {

		  m = idx(j,i);

		  if( isBin )
			fval = buf[ (j-1) * NLon + (i-1) ];

		  Node(m, iTopo) = fval;
		  Node(m, iTime) = -1;
		  Node(m, iD) = -fval;

		  if( Node(m, iD) < 0 ) {
			Node(m, iD) = 0.0f;
		  } else if( Node(m, iD) < MODELVARI.dmin ) {
			  Node(m, iD) = MODELVARI.dmin;
		  }

		}
	  }

	  delete[] buf;

  } else {

	  for( j=1; j<=NLat; j++ ) {
  		for( i=1; i<=NLon; i++ ) {

			m = idx(j,i);
			ierr = fscanf( fp, " %f ", &fval );

			Node(m, iTopo) = fval;
			Node(m, iTime) = -1;
			Node(m, iD) = -fval;

			if( Node(m, iD) < 0 ) {
			Node(m, iD) = 0.0f;
			} else if( Node(m, iD) < MODELVARI.dmin ) {
			  Node(m, iD) = MODELVARI.dmin;
			}

		}

	  }
  }

  for( k=1; k<MAX_VARS_PER_NODE-2; k++ ) {
	  Node.initMemory( k, 0 );
  }

  fclose( fp );

  if( !MODELVARI.DT ) { 

	double DTLOC=RealMax;

	for( i=1; i<=NLon; i++ ) {
	  for( j=1; j<=NLat; j++ ) {
		  m = idx(j,i);
		  if( Node(m, iD) == 0.0f ) continue;
		  DTLOC = My_min( DTLOC, 0.8 * (Dx*cosdeg(getLat(j))) / sqrt(GG*Node(m, iD)) );
	  }
	}

    if( DTLOC > 15 ) MODELVARI.DT = 15;
    else if( DTLOC > 10 ) MODELVARI.DT = 10;
    else if( DTLOC > 5 ) MODELVARI.DT = 5;
    else if( DTLOC > 2 ) MODELVARI.DT = 2;
    else if( DTLOC > 1 ) MODELVARI.DT = 1;
    else return Err.post("Bathymetry requires too small time step (<1sec)");
  }

   for( i=1; i<=NLon; i++ ) {
    if( Node(idx(1,i), iD) != 0 && Node(idx(2,i), iD) == 0 ) Node(idx(1,i), iD) = 0.;
    if( Node(idx(NLat,i), iD) != 0 && Node(idx(NLat-1,i), iD) == 0 ) Node(idx(NLat,i), iD) = 0.;
  }
  for( j=1; j<=NLat; j++ ) {
    if( Node(idx(j,1), iD) != 0 && Node(idx(j,2), iD) == 0 ) Node(idx(j,1), iD) = 0.;
    if( Node(idx(j,NLon), iD) != 0 && Node(idx(j,NLon-1), iD) == 0 ) Node(idx(j,NLon), iD) = 0.;
  }


  for( j=1; j<=NLat; j++ ) {
      // Setting up R6
    R6[j] = cosdeg( LatMin + (j-0.5)*DLat );
  }

  for( i=1; i<=NLon; i++ ) {
    for( j=1; j<=NLat; j++ ) {

      m = idx(j,i);

      if( Node(m, iD) == 0 ) continue;

      // Setting up R1
      Node(m, iR1) = MODELVARI.DT/Dy/R6[j];

      // Setting up R2 and R3
      if( i != NLon ) {
        if( Node(m+NLat, iD) != 0 ) {
          Node(m, iR2) = 0.5*GG*MODELVARI.DT/Dy/R6[j]*(Node(m, iD)+Node(m+NLat, iD));
          Node(m, iR3) = 0.5*MODELVARI.DT*Omega*sindeg( LatMin + (j-0.5)*DLat );
        }
      }
      else {
    	Node(m, iR2) = 0.5*GG*MODELVARI.DT/Dy/R6[j]*Node(m, iD)*2;
    	Node(m, iR3) = 0.5*MODELVARI.DT*Omega*sindeg( LatMin + (j-0.5)*DLat );
      }

      if( j != NLat ) {
        if( Node(m+1, iD) != 0 ) {
       // Setting up R4 and R5
         Node(m, iR4) = 0.5*GG*MODELVARI.DT/Dy*(Node(m, iD)+Node(m+1, iD));
          Node(m, iR5) = 0.5*MODELVARI.DT*Omega*sindeg( LatMin + j*DLat );
        }
      }
 
      else {
    	Node(m, iR2) = 0.5*GG*MODELVARI.DT/Dy*Node(m, iD)*2;
    	Node(m, iR3) = 0.5*MODELVARI.DT*Omega*sindeg( LatMin + j*DLat );
      }

    }
  }

  for( i=1; i<=NLon; i++ ) {
      // Setting up C1 and C3
    C1[i] = 0;
    if( Node(idx(1,i), iD) != 0 ) C1[i] = 1./sqrt(GG*Node(idx(1,i), iD));
    C3[i] = 0;
    if( Node(idx(NLat,i), iD) != 0 ) C3[i] = 1./sqrt(GG*Node(idx(NLat,i), iD));
  }

  for( j=1; j<=NLat; j++ ) {
      // Setting up C2 and C4
    C2[j] = 0;
    if( Node(idx(j,1), iD) != 0 ) C2[j] = 1./sqrt(GG*Node(idx(j,1), iD));
    C4[j] = 0;
    if( Node(idx(j,NLon), iD) != 0 ) C4[j] = 1./sqrt(GG*Node(idx(j,NLon), iD));
  }

  return 0;
}
/*
 *Implementation flow
 *
 *Read Input File of Deformation Grid
 **check input file type: GRD 
 **
 **load GRD file
 **
 **read grid resolution, grid dimensions 
 **integrate for tsunami energy
 **effective source
 **remove noise in the source
 **calculated (if needed) arrival threshold (negative value means it is relative)
 **transfer uplift onto tsunami grid and define deformed area for acceleration
 *
 */
int ICPROFPARA()
{
  char dsaa_label[8];
  int i,j,ierr,srcType;
  double lon,lat,dz,absuzmax,absuzmin;
  FILE *fp;
  GRDUTIL uZ;

  CNode& Node = *gNode;

  if( (fp = fopen( MODELVARI.fileSource, "rb" )) == NULL ) return Err.post( Err.msgOpenFile(MODELVARI.fileSource) );
  memset( dsaa_label, 0, 5 );
  ierr = fread( dsaa_label, 4, 1, fp );
  if( !strcmp( dsaa_label,"DSAA" ) || !strcmp( dsaa_label,"DSBB" ) )
    srcType = SRC_GRD;
  else
    srcType = SRC_GRD;

 fclose(fp);


  if( srcType == SRC_GRD) {
    ierr = uZ.readGRD( MODELVARI.fileSource ); if(ierr) return ierr;
  }

  absuzmax = uZ.getMaxAbsVal();

  if( (MODELVARI.zout0TRel + MODELVARI.zout0TAbs) != 0 ) {

    absuzmin = RealMax;
    if( MODELVARI.zout0TRel != 0 ) absuzmin = MODELVARI.zout0TRel*absuzmax;
    if( MODELVARI.zout0TAbs != 0 && MODELVARI.zout0TAbs < absuzmin ) absuzmin = MODELVARI.zout0TAbs;

    for( i=0; i<uZ.nx; i++ ) {
      for( j=0; j<uZ.ny; j++ ) {
        if( fabs(uZ(i,j)) < absuzmin ) uZ(i,j) = 0;
      }
    }

  }

  if( MODELVARI.zoutAT < 0 ) MODELVARI.zoutAT = absuzmax * fabs(MODELVARI.zoutAT);

  Imin = NLon; Imax = 1; Jmin = NLat; Jmax = 1;

  for( i=1; i<=NLon; i++ ) {
    for( j=1; j<=NLat; j++ ) {

      lon = getLon(i);
      lat = getLat(j);

      if( Node(idx(j,i), iD) != 0. )
        dz = Node(idx(j,i), iZ) = uZ.getVal( lon,lat );
      else
        dz = Node(idx(j,i), iZ) = 0.;

      if( fabs(dz) > MODELVARI.zoutCT ) {
        Imin = My_min( Imin, i );
        Imax = My_max( Imax, i );
        Jmin = My_min( Jmin, j );
        Jmax = My_max( Jmax, j );
      }

    }
  }

  if( Imin == NLon ) return Err.post( "Zero initial displacement" );

  Imin = My_max( Imin - 2, 2 );
  Imax = My_min( Imax + 2, NLon-1 );
  Jmin = My_max( Jmin - 2, 2 );
  Jmax = My_min( Jmax + 2, NLat-1 );

  return 0;
}
/*
 *Implementation flow
 *
 *
 */

int PROPAOUTFL()
{
  FILE *fp;
  char buf[64];

  NrecSpatialWaveOutput = 0;

  return 0;
}



int FILEOTPROPA()
{
  FILE *fp;
  short nOutI,nOutJ;
  int i,j,m;
  float ftmp;
  double dtmp,lonOutMin,lonOutMax,latOutMin,latOutMax;
  char record[128];

  CNode& Node = *gNode;

  NrecSpatialWaveOutput++;

  nOutI = Imax-Imin+1;
  lonOutMin = getLon(Imin); lonOutMax = getLon(Imax);
  nOutJ = Jmax-Jmin+1;
  latOutMin = getLat(Jmin); latOutMax = getLat(Jmax);

  sprintf( record, "zout_%5.5d.grd",MODELVARI.time );
 fp = fopen( record, "wb" );
  fwrite( "DSBB", 4, 1, fp );
  fwrite( &nOutI, sizeof(short), 1, fp );
  fwrite( &nOutJ, sizeof(short), 1, fp );
  fwrite( &lonOutMin, sizeof(double), 1, fp );
  fwrite( &lonOutMax, sizeof(double), 1, fp );
  fwrite( &latOutMin, sizeof(double), 1, fp );
  fwrite( &latOutMax, sizeof(double), 1, fp );
  dtmp = -1.; fwrite( &dtmp, sizeof(double), 1, fp );
  dtmp = +1.; fwrite( &dtmp, sizeof(double), 1, fp );
  for( j=Jmin; j<=Jmax; j++ ) {
    for( i=Imin; i<=Imax; i++ ) {
      m = idx(j,i);
      if( fabs(Node(m, iZ)) < MODELVARI.zoutTT )
        ftmp = (float)9999;
      else
        ftmp = (float)Node(m, iZ);
      fwrite( &ftmp, sizeof(float), 1, fp );
    }
  }
  fclose( fp );
  return 0;
}


int PROPAMAX()
{
  FILE *fp;
  short nOutI,nOutJ;
  int i,j,m;
  float ftmp;
  double dtmp,lonOutMin,lonOutMax,latOutMin,latOutMax;
  char record[128];

  CNode& Node = *gNode;

  nOutI = Imax-Imin+1;
  lonOutMin = getLon(Imin); lonOutMax = getLon(Imax);
  nOutJ = Jmax-Jmin+1;
  latOutMin = getLat(Jmin); latOutMax = getLat(Jmax);

  sprintf( record, "zmax_a.grd");
 fp = fopen( record, "wb" );
  fwrite( "DSBB", 4, 1, fp );
  fwrite( &nOutI, sizeof(short), 1, fp );
  fwrite( &nOutJ, sizeof(short), 1, fp );
  fwrite( &lonOutMin, sizeof(double), 1, fp );
  fwrite( &lonOutMax, sizeof(double), 1, fp );
  fwrite( &latOutMin, sizeof(double), 1, fp );
  fwrite( &latOutMax, sizeof(double), 1, fp );
  dtmp = 0.; fwrite( &dtmp, sizeof(double), 1, fp );
  dtmp = 1.; fwrite( &dtmp, sizeof(double), 1, fp );
  for( j=Jmin; j<=Jmax; j++ ) {
    for( i=Imin; i<=Imax; i++ ) {
      ftmp = (float)Node(idx(j,i), iZmax);
      fwrite( &ftmp, sizeof(float), 1, fp );
    }
  }
  fclose( fp );

  return 0;
}
/*
 *
 *Implementation flow 
 *
 *
 */

int RESETNODES()
{
  int i,j,m;

  CNode& Node = *gNode;

  for( i=1; i<=NLon; i++ ) {
    for( j=1; j<=NLat; j++ ) {

      m = idx(j,i);
      Node(m, iZ) = Node(m, iZmax) = Node(m, iM) = Node(m, iN) = 0;
      Node(m, iTime) = -1;
    }
  }

  Imin = 1; Imax = NLon; Jmin = 1; Jmax = NLat;

  return 0;
}
/*
 *
 *Implementation flow
 *
 *Constructor
 *Destructor
 *Reset observation array
 *Fully reset observation array to a new number of observations
 *Find tgpoint
 *
 *Read observations from text file
 **Get number of values per location
 **check if site ID's are available in the first column. Criterium: first character is not digit
 **nObs = number of columns minus 2 (lon lat) minus 1 (if idPresent)
 **Get number of positions
 **Allocate memory
 **Read data
 *
 *Write to simple text file
 *Write to binary stream
 *Read from binary stream
 *Calculate observation residual
 *Calculate norm of observations
 *
 *Haversine formula for distance between any two points on the Earth surface
 **multiplyer to convert from degrees into radians
 **Earth radius in km along equator
 *
 *GeoStrike Calculation on sphere
 **multiplyer to convert from degrees into radians
 **multiplyer to convert from radians into degrees
 **Earth radius in km along equator
 */
cObsArray::cObsArray()
{
  nPos = 0;
  nObs = 0;
  id = NULL;
  lon = NULL;
  lat = NULL;
  obs = NULL;
}
cObsArray::~cObsArray()
{
  if( lon ) delete [] lon;
  if( lat ) delete [] lat;
  if( obs ) {
    for( int n=0; n<nPos; n++ )
      delete [] obs[n];
    delete [] obs;
  }
  if( id ) {
    for( int n=0; n<nPos; n++ )
      delete [] id[n];
    delete [] id;
  }

}
int cObsArray::resetObs()
{
  if( !obs ) return 0;

  for( int n=0; n<nPos; n++ ) {
    memset( obs[n], 0, nObs*sizeof(double) );
  }

  return 0;
}
int cObsArray::resetObs( int newnobs )
{

  if( obs ) {
    for( int n=0; n<nPos; n++ )
      delete [] obs[n];
    delete [] obs;
  }

  nObs = newnobs;
  obs = new double* [nPos];
  for( int n=0; n<nPos; n++ ) {
    obs[n] = new double [nObs];
    memset( obs[n], 0, nObs*sizeof(double) );
  }

  return 0;
}
int cObsArray::findById( char *id0 )
{
  int n;

  for( n=0; n<nPos; n++ )
    if( !strcmp( id0, id[n] ) ) return n;

  return -1;
}
int cObsArray::read( char *fname )
{
  FILE *fp;
  char record[256], buf[64];
  int n,k,idPresent,line=0;
  double lonr,latr;


  fp = fopen(fname,"rt");
  if( fp == NULL ) return Err.post(Err.msgOpenFile(fname));

  line = 0;
  if( utlReadNextRecord( fp, record,  &line ) == EOF ) return Err.post("Unexpected EOF: %s", fname);
  fclose(fp);
  sscanf( record, "%s", buf );
  if( isalpha( buf[0] ) )
    idPresent = 1;
  else
    idPresent = 0;
  nObs = utlCountWordsInString( record ) - 2 - idPresent;
  if( nObs < 0 ) return Err.post(Err.msgReadFile( fname, line, "expected: [id lon lat [obs...]" ));

  nPos = utlGetNumberOfRecords( fname );

  lon = new double [nPos];
  lat = new double [nPos];
  if( idPresent ) {
    id = new char* [nPos];
    for( n=0; n<nPos; n++ )
      id[n] = NULL;
  }
  if( nObs > 0 ) {
    obs = new double* [nPos];
    for( n=0; n<nPos; n++ ) {
      obs[n] = new double [nObs];
      memset( obs[n], 0, nObs*sizeof(double) );
    }
  }

  fp = fopen(fname,"rt");
  line = 0;
  for( n=0; n<nPos; n++ ) {
    if( utlReadNextRecord( fp, record,  &line ) == EOF ) return Err.post("Unexpected EOF: %s", fname);

    if( utlCountWordsInString( record ) != (2+idPresent+nObs) ) return Err.post(Err.msgReadFile( fname, line, "invalid number of values" ));

    if( idPresent ) {
      if( sscanf( record, "%s %lf %lf", buf, &lon[n], &lat[n] ) != 3 ) return Err.post(Err.msgReadFile( fname, line, "expected: id lon lat obs..." ));
      id[n] = strdup(buf);
    }
    else {
      if( sscanf( record, "%lf %lf", &lon[n], &lat[n] ) != 2 ) return Err.post(Err.msgReadFile( fname, line, "expected: lon lat obs..." ));
    }

    for( k=0; k<nObs; k++ ) {
      if( utlPickUpWordFromString( record, 3+idPresent+k, buf ) != 0 ) return Err.post(Err.msgReadFile( fname, line, "expected: id lon lat obs..." ));
      if( sscanf( buf, "%lf", &obs[n][k] ) != 1 ) return Err.post(Err.msgReadFile( fname, line, "expected: id lon lat obs..." ));
    }
  }
  fclose( fp );

  return 0;
}
int cObsArray::write( char *fname )
{
  FILE *fp;
  int n,k;


  fp = fopen(fname,"wt");

  for( n=0; n<nPos; n++ ) {

    if( id )
      fprintf( fp, "%s %g %g", id[n], lon[n],lat[n] );
    else
      fprintf( fp, "%g %g", lon[n],lat[n] );

    for( k=0; k<nObs; k++ )
      fprintf( fp, " %g", obs[n][k] );

    fprintf( fp, "\n" );
  }

  fclose( fp );

  return 0;
}
long cObsArray::writeBin( FILE *fp )
{
  long bytes_written;
  float fbuf;


  bytes_written = 0;
  for( int n=0; n<nPos; n++ ) {
    for( int k=0; k<nObs; k++ ) {
      fbuf = (float)obs[n][k];
      fwrite( &fbuf, sizeof(float), 1, fp );
      bytes_written += sizeof(float);
    }
  }

  return bytes_written;
}

long cObsArray::readBin( FILE *fp )
{
  long bytes_read;
  float fbuf;


  bytes_read = 0;
  for( int n=0; n<nPos; n++ ) {
    for( int k=0; k<nObs; k++ ) {
      if( fread( &fbuf, sizeof(float), 1, fp ) != 1 )
        return utlPostError("Unexpected EOF");
      obs[n][k] = (double)fbuf;
      bytes_read += sizeof(float);
    }
  }

  return bytes_read;
}

double cObsArray::residual( cObsArray& ref )
{
  double resid=0.;

  for( int n=0; n<nPos; n++ ) {
    for( int k=0; k<nObs; k++ ) {
      resid += pow( (obs[n][k] - ref.obs[n][k]), 2. );
    }
  }

  resid = sqrt( resid )/nPos/nObs;

  return resid;
}

double cObsArray::norm()
{
  double norm=0.;

  for( int n=0; n<nPos; n++ ) {
    for( int k=0; k<nObs; k++ ) {
      norm += obs[n][k]*obs[n][k];
    }
  }

  norm = sqrt( norm )/nPos/nObs;

  return norm;
}

double GeoDistOnSphere( const double lon1, const double lat1, const double lon2, const double lat2 )
{
  const double G2R = 3.14159265358979/180.;  
  const double REARTH = 6378.137; 
  double a,c,dist,rad;

  a = pow( sin(G2R*(lat2-lat1)/2), 2. ) + cos(G2R*lat1)*cos(G2R*lat2)*pow( sin(G2R*(lon2-lon1)/2), 2. );
  rad = sqrt(a);
  if( rad > 1 ) rad = 1;
  c = 2 * asin(rad);
  dist = REARTH*c;

  return dist;
}

double GeoStrikeOnSphere( double lon1, double lat1, double lon2, double lat2 )
{
  const double G2R = 3.14159265358979/180.;  
  const double R2G = 180./3.14159265358979;  
  const double REARTH = 6378.137; 
  double strike, distRad;


  if( (lat1 == lat2) && (lon1 == lon2) ) {
    strike = 0.;
  }
  else if( lon1 == lon2 ) {
    if( lat1 > lat2 )
      strike = 180.;
    else
      strike = 0.;
  }
  else {
    distRad = GeoDistOnSphere( lon1,lat1,lon2,lat2 ) / REARTH;
    strike = R2G * asin( cos(G2R*lat2)*sin(G2R*(lon2-lon1)) / sin(distRad) );

    if( (lat2 > lat1) && (lon2 > lon1) ) {
    }
    else if( (lat2 < lat1) && (lon2 < lon1) ) {
      strike = 180.0 - strike;
    }
    else if( (lat2 < lat1) && (lon2 > lon1) ) {
      strike = 180.0 - strike;
    }
    else if( (lat2 > lat1) && (lon2 < lon1) ) {
      strike += 360.0;
    }

  }

  return strike;
}

