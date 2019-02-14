/*
 *      Tohoku University's Numerical-Analysis Model
 *          for Investigation of tsunami
 *           Far-field Tsunami version 
 */

#ifndef UTILMDL_H
#define UTILMDL_H

#include <math.h>

#define MaxFileRecordLength 16384
#define RealMax 1.e+30
#define RealMin 1.e-30

#define PI_RC  3.14159265358979
#define g2r(x)  (((double)(x))*PI_RC/180)
#define r2g(x)  (((double)(x))/PI_RC*180)
#define cosdeg(x) cos(g2r(x))
#define sindeg(x) sin(g2r(x))
#define tandeg(x) tan(g2r(x))

#define My_max(a,b)  (((a) > (b)) ? (a) : (b))
#define My_min(a,b)  (((a) < (b)) ? (a) : (b))

#define MSG_OUTSCRN 1
#define MSG_OUTFILE 2
class cMsg
{
protected:
  int enabled;
  int channel;
  char fname[64];

public:
  cMsg();
  ~cMsg();
  void enable();
  void disable();
  void setchannel( int newchannel );
  void setfilename( const char* newfilename );
  int print( const char* fmt, ... );
  };
extern cMsg Msg;


class cErrMsg : public cMsg
{

public:
  cErrMsg();
  ~cErrMsg();
  int post( const char* fmt, ... );
  char* msgAllocateMem();
  char* msgOpenFile( const char *fname );
  char* msgReadFile( const char *fname, int line, const char *expected );
  };
extern cErrMsg Err;


class cLog : public cMsg
{
private:
  int timestamp_enabled;

public:
  cLog();
  ~cLog();
  void start( const char* filename );
  void timestamp_enable();
  void timestamp_disable();
  int print( const char* fmt, ... );
  };
extern cLog Log;


char *utlCurrentTime();

int utlPostError( const char* message );
char *utlErrMsgMemory();
char *utlErrMsgOpenFile( const char *fname );
char *utlErrMsgEndOfFile( const char *fname );
char *utlErrMsgReadFile( const char *fname, int line, const char *expected );
int utlCheckCommandLineOption( int argc, char **argv, const char *option, int letters_to_compare );
int utlStringReadOption( char *record, char *option_name, char *contents );
int utlReadSequenceOfInt( char *line, int *value );
int utlReadSequenceOfDouble( char *line, double *value );
int utlPickUpWordFromString( char *string, int n, char *word );
char *utlPickUpWordFromString( char *string, int n );
char *utlPickUpWordFromString( char *string, char *pos, char *word );
int utlCountWordsInString( char *record );
char *utlWords2String( int nwords, char **word );
int utlSubString( char *logstr, int p1, int p2, char *substr );
void utlPrintToString( char *string, int position, char *insert );
int utlFindNextInputField( FILE *fp );
int utlReadNextRecord( FILE *fp, char *record, int *line );
int utlReadNextDataRecord( FILE *fp, char *record, int *line );
char *utlFileFindRecord( char *fname, char *pattern, char *record );
int utlFileParceToString( FILE *fp, char *pattern );
int utlGetNumberOfRecords( const char *fname );
int utlGetNumberOfDataRecords( const char *fname );
int utlFileReadOption( char *fname, char *option_name, char *contents );
char *utlFileChangeExtension( const char *fname, const char *newext );
char *utlFileRemoveExtension( const char *fname );
int utlFileRemoveExtension( char *newname, const char *fname );
char *utlFileAddExtension( const char *fname, const char *addext );
int utlReadXYdata( char *fname, double **x, double **y );
double utlRandom( int& seed );
double utlRandom( double avg, double err );
double utlNormal( int& seed, double avg, double stdev );
int utlPointInPolygon( int nvrtx, double *xvrtx, double *yvrtx, double x, double y );
void utlTimeSplit( double ctime, int& nHour, int& nMin, int& nSec, int& nMsec );
char *utlTimeSplitString( double ctime );
char *utlTimeSplitString( int ctime );
int utlTimeHour( double timesec );
int utlTimeHour( int timesec );
int utlTimeMin( double timesec );
int utlTimeMin( int timesec );
int utlTimeSec( double timesec );
int utlTimeSec( int timesec );

#endif  // UTILMDL_H
