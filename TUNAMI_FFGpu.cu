/*
 *      Tohoku University's Numerical-Analysis Model
 *          for Investigation of tsunami
 *           Far-field Tsunami version 
 *
 *TUNAMIFF main - Implementation flow 
 *
 *Read parameters and use default (as necessary)
 *Read bathymetry
 *Intial Condition - Init tsunami-grid
 *Copy to GPU
 **cuda memory copy for 2-dim(d, z, zMax, fM, fN, cR1, cR2, cR4, tArr)
 **as well for 1-dim (cR6, cB1, cB2, cB3, cB4, g_MinMax)
 **and move global variables into datastructure
 *Execute GnodeCudaTUNAMIFF Kernel
 *Execute MASSGBOUNDMOMENT Kernel
 *output on commnd line
 *Copy from GPU
 **cuda memory copy for 2-dim (zMax, tArr)
 **intermediate copy - ignore copy requests if data already present on CPU side
 **And copy as necessary with data change 
 *Write out files
 *Free up CUDA node memory 
 */

#define HEADER "\n \n \
 *      Tohoku University's Numerical-Analysis Model \n \
 *          for Investigation of tsunami \n \
 *           Far-field Tsunami version\n \
 *                      2011-05-21\n \
\n \n"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "UTILMDL.h"
#include "TUNAMI_FF.h"

#ifdef __CUDACC__
#include "GnodeCudaTUNAMIFF.cuh"
#endif

CNode *gNode;

double diff(timespec start, timespec end) {

	timespec temp;
	if ((end.tv_nsec-start.tv_nsec)<0) {
		temp.tv_sec = end.tv_sec-start.tv_sec-1;
		temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
	} else {
		temp.tv_sec = end.tv_sec-start.tv_sec;
		temp.tv_nsec = end.tv_nsec-start.tv_nsec;
	}

	return (double)((double)temp.tv_nsec / 1000000000.0 + (double)temp.tv_sec);
}


int main( int argc, char **argv )
{
  char buf[1024];
  int ierr,argn,elapsed;
  int lastProgress,lastPropagation,lastDump;
  int loop;

  printf(HEADER);
  Err.setchannel(MSG_OUTFILE);

  ierr = MODELPARAINITIA( argc, argv ); if(ierr) return ierr;


  for( argn=1; argn<argc; argn++ ) {
    strcat( buf, " " );
    strcat( buf, argv[argn] );
  }

  if( MODELVARI.gpu ) {
#ifdef __CUDACC__
	  gNode = new CGpuNode();
#endif
  } else {
	  gNode = new CStructNode();
	  //gNode = new CArrayNode();
  }

  CNode& Node = *gNode;

  ierr = tfRdepthParame(); if(ierr) return ierr;

  ierr = ICPROFPARA(); if(ierr) return ierr;

  if( MODELVARI.KDspatialOut ) PROPAOUTFL();

  Node.copyToGPU();

  // Main loop

  timespec start, end;
  clock_gettime(CLOCK_MONOTONIC, &start);

  for( MODELVARI.time=0,loop=1,lastProgress=MODELVARI.outProgress,lastPropagation=MODELVARI.KDspatialOut,lastDump=0;
    MODELVARI.time<=MODELVARI.KEtotalTime; loop++,MODELVARI.time+=MODELVARI.DT,lastProgress+=MODELVARI.DT,lastPropagation+=MODELVARI.DT ) {

      if( MODELVARI.fileGAGUEs && MODELVARI.KCgagueOUT && ((MODELVARI.time/MODELVARI.KCgagueOUT)*MODELVARI.KCgagueOUT == MODELVARI.time) ) {
    }

    Node.run();

    elapsed = ((int)clock())/CLOCKS_PER_SEC;

    if( MODELVARI.outProgress ) {
      if( lastProgress >= MODELVARI.outProgress ) {
        printf( "%s\n", utlTimeSplitString(MODELVARI.time) );
        lastProgress = 0;
      }
    }

    fflush(stdout);

    if( MODELVARI.KDspatialOut ) {
      if( lastPropagation >= MODELVARI.KDspatialOut ) {
    	Node.copyIntermediate();
        FILEOTPROPA();
        lastPropagation = 0;
      }
    }

    if( MODELVARI.MAXOUT ) {
      if( (elapsed-lastDump) >= MODELVARI.MAXOUT ) {
    	Node.copyIntermediate();
        PROPAMAX();
        lastDump = elapsed;
      }
    }

  } // main loop
  clock_gettime(CLOCK_MONOTONIC, &end);
  Node.copyIntermediate();
  Node.copyFromGPU();

  PROPAMAX();

  Node.freeMem();

  printf_v("Executiontime (sec): %.3lf\n", diff(start, end));

  delete gNode;

  return 0;
}


