/*
 *      Tohoku University's Numerical-Analysis Model
 *          for Investigation of tsunami
 *           Far-field Tsunami version 
 */

#ifndef TF_GPUNODE_H
#define TF_GPUNODE_H


#include "TUNAMI_FF.h"
#include "NODEINI.h"
#include <stdio.h>

#define CUDA_CALL(x) if( x != cudaSuccess ) { fprintf( stderr, "Error in file %s on line %u: %s\n", __FILE__, __LINE__, cudaGetErrorString( cudaGetLastError() ) ); return 1; }

#undef idx

class Params {

public:
	int mTime;
	int nI;
	int nJ;
	int iMin;
	int iMax;
	int jMin;
	int jMax;
	float zoutAT;
	float zoutZT;
	float zoutCT;


	size_t pI;
	size_t lpad;
};

class KernelData {

public:

	float *d;
	float *z;
	float *zMax;
	float *fM;
	float *fN;
	float *cR1;
	float *cR2;
	float *cR4;
	float *tArr;

	float *cR6;
	float *cB1;
	float *cB2;
	float *cB3;
	float *cB4;

	Params params;

	int4 *g_MinMax;

	__device__ int le( int ij ) { return ij - params.pI; }
	__device__ int ri( int ij ) { return ij + params.pI; }
	__device__ int up( int ij ) { return ij + 1; }
	__device__ int dn( int ij ) { return ij - 1; }
	__host__ __device__ int idx( int i, int j ) { return (j-1) + (i-1) * params.pI + params.lpad; }
};



class CGpuNode : public CArrayNode {

protected:
	KernelData data;


	size_t pitch;

	bool copied;

	cudaEvent_t evtStart[5];
	cudaEvent_t evtEnd[5];
	float dur[5];

public:
	CGpuNode();
	int mallocMem();
	int copyToGPU();
	int copyFromGPU();
	int copyIntermediate();
	int copyGAGUEs();
	int freeMem();
	int run();
};

#endif /* TF_GPUNODE_H */

