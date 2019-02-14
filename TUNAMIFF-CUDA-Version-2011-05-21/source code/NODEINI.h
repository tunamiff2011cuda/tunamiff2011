/*
 *      Tohoku University's Numerical-Analysis Model
 *          for Investigation of tsunami
 *           Far-field Tsunami version 
 */

#ifndef NODEINI_H
#define NODEINI_H

#include <stdlib.h>
#include <string.h>
#include "TUNAMI_FF.h"

#define CHKRET( x ) if( (x) == NULL ) return 1;

typedef float Float[MAX_VARS_PER_NODE];

class CNode {

public:
	virtual ~CNode() {};
	virtual float& operator()( const int idx1, const int idx2 ) = 0;
	virtual int mallocMem() = 0;
	virtual int copyToGPU() = 0;
	virtual int copyFromGPU() = 0;
	virtual int copyIntermediate() = 0;
	virtual int copyGAGUEs() = 0;
	virtual int freeMem() = 0;
	virtual int run() = 0;

	virtual void initMemory( int index, int val ) = 0;
};

class CStructNode : public CNode {

public:
	Float *node;

public:
	inline float& operator()( const int idx1, const int idx2 ) {

		return node[idx1][idx2];
	}

	void initMemory( int index, int val ) {

		int m;
		for( int i=1; i<=NLon; i++ ) {
		  for( int j=1; j<=NLat; j++ ) {
				m = idx(j,i);
				this->operator ()(m, index) = val;
		  }
		}
	}

	int mallocMem() {

		CHKRET( this->node = (Float*) malloc( sizeof(Float) * NLon * NLat) );

		CHKRET( R6 = (float*) malloc( sizeof(float) * (NLat+1) ) );
		CHKRET( C1 = (float*) malloc( sizeof(float) * (NLon+1) ) );
		CHKRET( C3 = (float*) malloc( sizeof(float) * (NLon+1) ) );
		CHKRET( C2 = (float*) malloc( sizeof(float) * (NLat+1) ) );
		CHKRET( C4 = (float*) malloc( sizeof(float) * (NLat+1) ) );

		return 0;
	}

	int freeMem() {

		free( this->node );
		free( R6 );
		free( C1 );
		free( C2 );
		free( C3 );
		free( C4 );

		return 0;
	}

	int run() {

		if( MODELVARI.coriolis )
		  return MASSGBOUNDMOMENTCor();

		return MASSGBOUNDMOMENT();
	}

	int copyToGPU() { return 0; }
	int copyFromGPU() {	return 0; }
	int copyIntermediate() { return 0; }
	int copyGAGUEs() { return 0; }

};

#pragma pack(push, 1)
class CArrayNode : public CNode {

protected:
	float *d;
	float *h;
	float *zMax;
	float *fM;
	float *fN;
	float *cR1;
	float *cR2;
	float *cR3;
	float *cR4;
	float *cR5;
	float *tArr;
	float *topo;

public:
	virtual float& operator()( const int idx1, const int idx2 ) {

		return ((float**)&d)[idx2][idx1];
	}

	void *getBuf( int idx ) { return ((float**)&d)[idx]; }

	virtual void initMemory( int index, int val ) {

		memset( getBuf(index), 0, NLat * NLon * sizeof(float) );

	}

	virtual int mallocMem() {

		CHKRET( this->d = (float*) malloc( sizeof(float) * NLon * NLat ) );
		CHKRET( this->h = (float*) malloc( sizeof(float) * NLon * NLat ) );
		CHKRET( this->zMax = (float*) malloc( sizeof(float) * NLon * NLat ) );
		CHKRET( this->fM = (float*) malloc( sizeof(float) * NLon * NLat ) );
		CHKRET( this->fN = (float*) malloc( sizeof(float) * NLon * NLat ) );
		CHKRET( this->cR1 = (float*) malloc( sizeof(float) * NLon * NLat ) );
		CHKRET( this->cR2 = (float*) malloc( sizeof(float) * NLon * NLat ) );
		CHKRET( this->cR3 = (float*) malloc( sizeof(float) * NLon * NLat ) );
		CHKRET( this->cR4 = (float*) malloc( sizeof(float) * NLon * NLat ) );
		CHKRET( this->cR5 = (float*) malloc( sizeof(float) * NLon * NLat ) );
		CHKRET( this->tArr = (float*) malloc( sizeof(float) * NLon * NLat ) );
		CHKRET( this->topo = (float*) malloc( sizeof(float) * NLon * NLat ) );


        CHKRET( R6 = (float*) malloc( sizeof(float) * (NLat+1) ) );
		CHKRET( C1 = (float*) malloc( sizeof(float) * (NLon+1) ) );
		CHKRET( C3 = (float*) malloc( sizeof(float) * (NLon+1) ) );
		CHKRET( C2 = (float*) malloc( sizeof(float) * (NLat+1) ) );
		CHKRET( C4 = (float*) malloc( sizeof(float) * (NLat+1) ) );

		return 0;
	}

	virtual int freeMem() {

		free( this->d );
		free( this->h );
		free( this->zMax );
		free( this->fM );
		free( this->fN );
		free( this->cR1 );
		free( this->cR2 );
		free( this->cR3 );
		free( this->cR4 );
		free( this->cR5 );
		free( this->tArr );
		free( this->topo );

		free( R6 );
		free( C1 );
		free( C2 );
		free( C3 );
		free( C4 );

		return 0;
	}

	virtual int run() {

		if( MODELVARI.coriolis )
		  return MASSGBOUNDMOMENTCor();

		return MASSGBOUNDMOMENT();
	}

	virtual int copyToGPU() { return 0; }
	virtual int copyFromGPU() {	return 0; }
	virtual int copyIntermediate() { return 0; }
	virtual int copyGAGUEs() { return 0; }

};
#pragma pack(pop)

#endif /* NODEINI_H */

