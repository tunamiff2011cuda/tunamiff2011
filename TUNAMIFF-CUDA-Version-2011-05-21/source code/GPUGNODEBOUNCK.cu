/*
 *      Tohoku University's Numerical-Analysis Model
 *          for Investigation of tsunami
 *           Far-field Tsunami version 
 *
 */

#include "GnodeCudaTUNAMIFF.cuh"
#include "MASSGBOUNDMOMENT.cuh"
#include <cstring>

/*
 *
 *Implementation flow
 *
 *Computation of the equation of continuity
 *Open sea, boundary condition, open boundaries
 *Computation of the equation of motion
 *Making area of computation be with in the area under consideration and
 *Enlargement of the area of computation as the tsunami propagates
 */

__global__ void runMASSKernel( KernelData data ) {
//Computation of the equation of continuity	MASS
/*
VARIABLES
Water Level 			z
Discharge In I-direction	fM
Discharge In J-direction	fN
Area where the tsunami exists 	iMin, jMin, iMax, jMax
and the computation is carried out
//
Highest water level	zMax
zoutZT
COEFFICIENTS:
Coefficients given 	cR1 and cR6=COS (THETA fM+1/2)
*/

  Params& dp = data.params;

  int i = blockIdx.y * blockDim.y + threadIdx.y + dp.iMin;
  int j = blockIdx.x * blockDim.x + threadIdx.x + dp.jMin;
  int ij = data.idx(i,j);
  float absH;


  if( i <= dp.iMax && j <= dp.jMax && data.d[ij] != 0 ) {

	  float zz = data.z[ij] - data.cR1[ij] * ( data.fM[ij] - data.fM[data.le(ij)] + data.fN[ij] * data.cR6[j] - data.fN[data.dn(ij)]*data.cR6[j-1] );

	  absH = fabs(zz);

	  if( absH < dp.zoutZT ) {
		  zz = 0.f;
	  } else if( zz > data.zMax[ij] ) {
		  data.zMax[ij] = zz;
		  
	  }

	  if( dp.zoutAT && data.tArr[ij] < 0 && absH > dp.zoutAT )
	  	  data.tArr[ij] = dp.mTime;

	  data.z[ij] = zz;
  }

}

__global__ void runMOMENTKernel( KernelData data ) {
//Computation of the equation of motion MOMENT
/*
VARIABLES
Water Level 			z
Discharge In I-direction	fM
Discharge In J-direction	fN
Area where the tsunami exists 	iMin, jMin, iMax, jMax
and the computation is carried out
//
COEFFICIENTS:
Coefficients given 	cR2 and cR4
*/

	Params& dp = data.params;

	int i = blockIdx.y * blockDim.y + threadIdx.y + dp.iMin;
	int j = blockIdx.x * blockDim.x + threadIdx.x + dp.jMin;
	int ij = data.idx(i,j);

	if( i <= dp.iMax && j <= dp.jMax && data.d[ij] != 0 ) {

	  float zz = data.z[ij];

	  if( data.d[data.ri(ij)] != 0 ) {
		  data.fM[ij] = data.fM[ij] - data.cR2[ij]*(data.z[data.ri(ij)] - zz);
	  }

	  if( data.d[data.up(ij)] != 0 )
		  data.fN[ij] = data.fN[ij] - data.cR4[ij]*(data.z[data.up(ij)] - zz);

	}

}

__global__ void runGBOUNDKernel( KernelData data ) {
//Open sea, boundary condition GBOUND
/*
VARIABLES
Water Level 			z
Discharge In I-direction	fM
Discharge In J-direction	fN
//
COEFFICIENTS:
Coefficients given 	
(THETA M+1/2)in radian  = cB1
(THETA M)in radian 	= cB2
(THETA M-1/2)in radian	= cB3
coefficent		= cB4
*/

	KernelData& DT = data;
	Params& dp = data.params;

	int id = blockIdx.x * blockDim.x + threadIdx.x + 2;
	int ij;

	if( id <= dp.nI-1 ) {
	  ij = DT.idx(id,1);
	  DT.z[ij] = sqrtf( powf(DT.fN[ij],2.0f) + 0.25f*powf((DT.fM[ij] + DT.fM[DT.le(ij)]),2.0f) )*DT.cB1[id-1];
	  if( DT.fN[ij] > 0 ) DT.z[ij] = -DT.z[ij];
	}

	if( id <= dp.nI-1 ) {
	  ij = DT.idx(id,dp.nJ);
	  DT.z[ij] = sqrtf( powf(DT.fN[DT.dn(ij)],2.0f) + 0.25f*powf((DT.fM[ij] + DT.fM[DT.dn(ij)]),2.0f) )*DT.cB3[id-1];
	  if( DT.fN[DT.dn(ij)] < 0 ) DT.z[ij] = -DT.z[ij];
	}

	if( id <= dp.nJ-1 ) {
	  ij = DT.idx(1,id);
	  DT.z[ij] = sqrtf( powf(DT.fM[ij],2.0f) + 0.25f*powf((DT.fN[ij] + DT.fN[DT.dn(ij)]),2.0f) )*DT.cB2[id-1];
	  if( DT.fM[ij] > 0 ) DT.z[ij] = -DT.z[ij];
	}

	if( id <= dp.nJ-1 ) {
	  ij = DT.idx(dp.nI,id);
	  DT.z[ij] = sqrtf( powf(DT.fM[DT.le(ij)],2.0f) + 0.25f*powf((DT.fN[ij] + DT.fN[DT.dn(ij)]),2.0f) )*DT.cB4[id-1];
	  if( DT.fM[DT.le(ij)] < 0 ) DT.z[ij] = -DT.z[ij];
	}

	if( id == 2 ) {
	  ij = DT.idx(1,1);
	  DT.z[ij] = sqrtf( powf(DT.fM[ij],2.0f) + powf(DT.fN[ij],2.0f) )*DT.cB1[0];
	  if( DT.fN[ij] > 0 ) DT.z[ij] = -DT.z[ij];

	  ij = DT.idx(dp.nI,1);
	  DT.z[ij] = sqrtf( powf(DT.fM[DT.le(ij)],2.0f) + powf(DT.fN[ij],2.0f) )*DT.cB1[dp.nI-1];
	  if( DT.fN[ij] > 0 ) DT.z[ij] = -DT.z[ij];

	  ij = DT.idx(1,dp.nJ);
	  DT.z[ij] = sqrtf( powf(DT.fM[ij],2.0f) + powf(DT.fN[DT.dn(ij)],2.0f) )*DT.cB3[0];
	  if( DT.fN[DT.dn(ij)] < 0 ) DT.z[ij] = -DT.z[ij];

	  ij = DT.idx(dp.nI,dp.nJ);
	  DT.z[ij] = sqrtf( powf(DT.fM[DT.le(ij)],2.0f) + powf(DT.fN[DT.dn(ij)],2.0f) )*DT.cB3[dp.nI-1];
	  if( DT.fN[DT.dn(ij)] < 0 ) DT.z[ij] = -DT.z[ij];
	}
}

__global__ void runMOMENTBOUNDKernel( KernelData data ) {
/*
VARIABLES
Water Level 			z
Discharge In I-direction	fM
Discharge In J-direction	fN
//
COEFFICIENTS:
Coefficients given 	cR2 and cR4
*/

	KernelData& DT = data;
	Params& dp = data.params;

	int id = blockIdx.x * blockDim.x + threadIdx.x + 1;
	int ij;

	if( id <= dp.nI-1 ) {
	  ij = DT.idx(id,1);
	  DT.fM[ij] = DT.fM[ij] - DT.cR2[ij]*(DT.z[DT.ri(ij)] - DT.z[ij]);
	}

	if( id <= dp.nJ ) {
	  ij = DT.idx(1,id);
	  DT.fM[ij] = DT.fM[ij] - DT.cR2[ij]*(DT.z[DT.ri(ij)] - DT.z[ij]);
	}

	if( id <= dp.nI-1 ) {
	  ij = DT.idx(id,dp.nJ);
	  DT.fM[ij] = DT.fM[ij] - DT.cR2[ij]*(DT.z[DT.ri(ij)] - DT.z[ij]);
	}

	if( id <= dp.nJ-1 ) {
	  ij = DT.idx(1,id);
	  DT.fN[ij] = DT.fN[ij] - DT.cR4[ij]*(DT.z[DT.up(ij)] - DT.z[ij]);
	}

	if( id <= dp.nI ) {
	  ij = DT.idx(id,1);
	  DT.fN[ij] = DT.fN[ij] - DT.cR4[ij]*(DT.z[DT.up(ij)] - DT.z[ij]);
	}

	if( id <= dp.nJ-1 ) {
	  ij = DT.idx(dp.nI,id);
	  DT.fN[ij] = DT.fN[ij] - DT.cR4[ij]*(DT.z[DT.up(ij)] - DT.z[ij]);
	}

}

__global__ void runAlimitBlimitGridKernel( KernelData data ) {

/*
//Making area of computation be with in the area under consideration as wellas
//Enlargement of the area of computation as the tsunami propagates
*/
/*
VARIABLES
Water Level 			z
Area where the tsunami exists 	iMin, jMin, iMax, jMax
and the computation is carried out
zoutCT, x, y, w
*/

	Params& dp = data.params;

	int id = blockIdx.x * blockDim.x + threadIdx.x + 1;

#if (__CUDA_ARCH__ >= 130)

	if( id >= dp.jMin && id <= dp.jMax ) {

	  if( fabsf(data.z[data.idx(dp.iMin+2,id)]) > dp.zoutCT )
		  atomicAdd( &(data.g_MinMax->x), 1 );

	  if( fabsf(data.z[data.idx(dp.iMax-2,id)]) > dp.zoutCT )
		  atomicAdd( &(data.g_MinMax->y), 1 );
	}

	if( id >= dp.iMin && id <= dp.iMax ) {

	  if( fabsf(data.z[data.idx(id,dp.jMin+2)]) > dp.zoutCT )
		  atomicAdd( &(data.g_MinMax->z), 1 );

	  if( fabsf(data.z[data.idx(id,dp.jMax-2)]) > dp.zoutCT )
		  atomicAdd( &(data.g_MinMax->w), 1 );
	}

#else

        if( id == 1 ) {

          for( int j = dp.jMin; j <= dp.jMax; j++ ) {
            
            if( fabsf(data.z[data.idx(dp.iMin+2,j)]) > dp.zoutCT ) {
                data.g_MinMax->x = 1;
                break;
            }

          }

          for( int j = dp.jMin; j <= dp.jMax; j++ ) {

            if( fabsf(data.z[data.idx(dp.iMax-2,j)]) > dp.zoutCT ) {
               data.g_MinMax->y = 1;
               break;
            }

          }

          for( int i = dp.iMin; i <= dp.iMax; i++ ) {
        
            if( fabsf(data.z[data.idx(i,dp.jMin+2)]) > dp.zoutCT ) {
              data.g_MinMax->z = 1;
              break;
            }

          }

          for( int i = dp.iMin; i <= dp.iMax; i++ ) {

            if( fabsf(data.z[data.idx(i,dp.jMax-2)]) > dp.zoutCT ) {
              data.g_MinMax->w = 1;
              break;
            }

          }

        }

#endif

}
/*
 *
 *Implementation flow
 *
 *Initilisation of Parameters and Fields
 *CUDA memory allocation for 2-dim (d, z, zMax, fM, fN, cR1, cR2, cR4, tArr)
 *cR3, cR5 for coriolis
 * and 1-dim (cR6, cB1, cB2, cB3, cB4, g_MinMax)
 * copy to gpu
 **cuda memory copy for 2-dim(d, z, zMax, fM, fN, cR1, cR2, cR4, tArr)
 **as well for 1-dim (cR6, cB1, cB2, cB3, cB4, g_MinMax)
 ** and move global variables into datastructure
 * copy from gpu
 **cuda memory copy for s-dim (zMax, tArr)
 **intermediate copy - ignore copy requests if data already present on CPU side
 *Run kernal on GPU nodes; And copy as necessary with data change 
 *CUDA memory free for all 2-dim and 1-dim variables
 */

CGpuNode::CGpuNode() {

	pitch = 0;
	copied = true;

	for( int i = 0; i < 5; i++ ) {
		cudaEventCreate( &(evtStart[i]) );
		cudaEventCreate( &(evtEnd[i]) );
		dur[i] = 0.0;
	}

}

int CGpuNode::mallocMem() {

	CArrayNode::mallocMem();

	Params& dp = data.params;


	dp.nI = NLon;
	dp.nJ = NLat;
	dp.zoutAT = MODELVARI.zoutAT;
	dp.zoutCT = MODELVARI.zoutCT;
	dp.zoutZT = MODELVARI.zoutZT;
	dp.lpad = 31;

	size_t nJ_aligned = dp.nJ + dp.lpad;


	CUDA_CALL( cudaMallocPitch( &(data.d), &pitch, nJ_aligned * sizeof(float), dp.nI ) );
	CUDA_CALL( cudaMallocPitch( &(data.z), &pitch, nJ_aligned * sizeof(float), dp.nI ) );
	CUDA_CALL( cudaMallocPitch( &(data.zMax), &pitch, nJ_aligned * sizeof(float), dp.nI ) );
	CUDA_CALL( cudaMallocPitch( &(data.fM), &pitch, nJ_aligned * sizeof(float), dp.nI ) );
	CUDA_CALL( cudaMallocPitch( &(data.fN), &pitch, nJ_aligned * sizeof(float), dp.nI ) );
	CUDA_CALL( cudaMallocPitch( &(data.cR1), &pitch, nJ_aligned * sizeof(float), dp.nI ) );
	CUDA_CALL( cudaMallocPitch( &(data.cR2), &pitch, nJ_aligned * sizeof(float), dp.nI ) );
	CUDA_CALL( cudaMallocPitch( &(data.cR4), &pitch, nJ_aligned * sizeof(float), dp.nI ) );
	CUDA_CALL( cudaMallocPitch( &(data.tArr), &pitch, nJ_aligned * sizeof(float), dp.nI ) );



	CUDA_CALL( cudaMalloc( &(data.cR6), dp.nJ * sizeof(float) ) );
	CUDA_CALL( cudaMalloc( &(data.cB1), dp.nI * sizeof(float) ) );
	CUDA_CALL( cudaMalloc( &(data.cB2), dp.nJ * sizeof(float) ) );
	CUDA_CALL( cudaMalloc( &(data.cB3), dp.nI * sizeof(float) ) );
	CUDA_CALL( cudaMalloc( &(data.cB4), dp.nJ * sizeof(float) ) );

	CUDA_CALL( cudaMalloc( &(data.g_MinMax), sizeof(int4) ) );


	dp.pI = pitch / sizeof(float);

	return 0;
}

int CGpuNode::copyToGPU() {

	Params& dp = data.params;


        Jmin -= (Jmin-2) % 32;
       
  
        dp.iMin = Imin;
	dp.iMax = Imax;
        dp.jMin = Jmin;
	dp.jMax = Jmax;


	CUDA_CALL( cudaMemcpy2D( data.d + dp.lpad, pitch, d, dp.nJ * sizeof(float), dp.nJ * sizeof(float), dp.nI, cudaMemcpyHostToDevice ) );
	CUDA_CALL( cudaMemcpy2D( data.z + dp.lpad, pitch, h, dp.nJ * sizeof(float), dp.nJ * sizeof(float), dp.nI, cudaMemcpyHostToDevice ) );
	CUDA_CALL( cudaMemcpy2D( data.zMax + dp.lpad, pitch, zMax, dp.nJ * sizeof(float), dp.nJ * sizeof(float), dp.nI, cudaMemcpyHostToDevice ) );
	CUDA_CALL( cudaMemcpy2D( data.fM + dp.lpad, pitch, fM, dp.nJ * sizeof(float), dp.nJ * sizeof(float), dp.nI, cudaMemcpyHostToDevice ) );
	CUDA_CALL( cudaMemcpy2D( data.fN + dp.lpad, pitch, fN, dp.nJ * sizeof(float), dp.nJ * sizeof(float), dp.nI, cudaMemcpyHostToDevice ) );
	CUDA_CALL( cudaMemcpy2D( data.cR1 + dp.lpad, pitch, cR1, dp.nJ * sizeof(float), dp.nJ * sizeof(float), dp.nI, cudaMemcpyHostToDevice ) );
	CUDA_CALL( cudaMemcpy2D( data.cR2 + dp.lpad, pitch, cR2, dp.nJ * sizeof(float), dp.nJ * sizeof(float), dp.nI, cudaMemcpyHostToDevice ) );
	CUDA_CALL( cudaMemcpy2D( data.cR4 + dp.lpad, pitch, cR4, dp.nJ * sizeof(float), dp.nJ * sizeof(float), dp.nI, cudaMemcpyHostToDevice ) );
	CUDA_CALL( cudaMemcpy2D( data.tArr + dp.lpad, pitch, tArr, dp.nJ * sizeof(float), dp.nJ * sizeof(float), dp.nI, cudaMemcpyHostToDevice ) );

	CUDA_CALL( cudaMemcpy( data.cR6, R6, dp.nJ * sizeof(float), cudaMemcpyHostToDevice ) );
	CUDA_CALL( cudaMemcpy( data.cB1, C1, dp.nI * sizeof(float), cudaMemcpyHostToDevice ) );
	CUDA_CALL( cudaMemcpy( data.cB2, C2, dp.nJ * sizeof(float), cudaMemcpyHostToDevice ) );
	CUDA_CALL( cudaMemcpy( data.cB3, C3, dp.nI * sizeof(float), cudaMemcpyHostToDevice ) );
	CUDA_CALL( cudaMemcpy( data.cB4, C4, dp.nJ * sizeof(float), cudaMemcpyHostToDevice ) );

	return 0;
}
int CGpuNode::copyFromGPU() {

	Params& dp = data.params;

	CUDA_CALL( cudaMemcpy2D( zMax, dp.nJ * sizeof(float), data.zMax + dp.lpad, pitch, dp.nJ * sizeof(float), dp.nI, cudaMemcpyDeviceToHost ) );
	CUDA_CALL( cudaMemcpy2D( tArr, dp.nJ * sizeof(float), data.tArr + dp.lpad, pitch, dp.nJ * sizeof(float), dp.nI, cudaMemcpyDeviceToHost ) );

	return 0;
}

int CGpuNode::copyIntermediate() {


	if( copied )
		return 0;

	Params& dp = data.params;

	CUDA_CALL( cudaMemcpy2D( h, dp.nJ * sizeof(float), data.z + dp.lpad, pitch, dp.nJ * sizeof(float), dp.nI, cudaMemcpyDeviceToHost ) );


	copied = true;

	return 0;
}

int CGpuNode::copyGAGUEs() {

	Params& dp = data.params;

	if( copied )
		return 0;

	for( int n = 0; n < NGAGUEs; n++ ) {

		int i = idxGAGUE[n] / dp.nJ + 1;
		int j = idxGAGUE[n] % dp.nJ + 1;

		int id = data.idx( i, j );

		CUDA_CALL( cudaMemcpy( h + idxGAGUE[n], data.z + dp.lpad + id, sizeof(float), cudaMemcpyDeviceToHost ) );
	}

	return 0;
}

int CGpuNode::freeMem() {


	CUDA_CALL( cudaFree( data.d ) );
	CUDA_CALL( cudaFree( data.z ) );
	CUDA_CALL( cudaFree( data.zMax ) );
	CUDA_CALL( cudaFree( data.fM ) );
	CUDA_CALL( cudaFree( data.fN ) );
	CUDA_CALL( cudaFree( data.cR1 ) );
	CUDA_CALL( cudaFree( data.cR2 ) );
	CUDA_CALL( cudaFree( data.cR4 ) );
	CUDA_CALL( cudaFree( data.tArr ) );

	CUDA_CALL( cudaFree( data.cR6 ) );
	CUDA_CALL( cudaFree( data.cB1 ) );
	CUDA_CALL( cudaFree( data.cB2 ) );
	CUDA_CALL( cudaFree( data.cB3 ) );
	CUDA_CALL( cudaFree( data.cB4 ) );

	CUDA_CALL( cudaFree( data.g_MinMax ) );

	float total_dur = 0.f;
	for( int j = 0; j < 5; j++ ) {
		printf_v("Duration %u (msec): %.3f\n", j, dur[j]);
		total_dur += dur[j];
	}
	printf_v("Duration total (msec): %.3f\n",total_dur);

	CArrayNode::freeMem();

	return 0;
}

int CGpuNode::run() {

	Params& dp = data.params;

	int nThreads = 256;
	int xThreads = 32;
	int yThreads = nThreads / xThreads;

	int NJ = dp.jMax - dp.jMin + 1;
	int NI = dp.iMax - dp.iMin + 1;
	int xBlocks = ceil( (float)NJ / (float)xThreads );
	int yBlocks = ceil( (float)NI / (float)yThreads );

	dim3 threads( xThreads, yThreads );
	dim3 blocks( xBlocks, yBlocks );

	int nBlocks = ceil( (float)max(dp.nI,dp.nJ) / (float)nThreads );

	dp.mTime = MODELVARI.time;

	CUDA_CALL( cudaEventRecord( evtStart[0], 0 ) );
	runMASSKernel<<<blocks,threads>>>( data );
	CUDA_CALL( cudaEventRecord( evtEnd[0], 0 ) );
	CUDA_CALL( cudaEventRecord( evtStart[1], 0 ) );
	runGBOUNDKernel<<<nBlocks,nThreads>>>( data );
	CUDA_CALL( cudaEventRecord( evtEnd[1], 0 ) );
	CUDA_CALL( cudaEventRecord( evtStart[2], 0 ) );
	runMOMENTKernel<<<blocks,threads>>>( data );
	CUDA_CALL( cudaEventRecord( evtEnd[2], 0 ) );
	CUDA_CALL( cudaEventRecord( evtStart[3], 0 ) );
	runMOMENTBOUNDKernel<<<nBlocks,nThreads>>>( data );
	CUDA_CALL( cudaEventRecord( evtEnd[3], 0 ) );
	CUDA_CALL( cudaEventRecord( evtStart[4], 0 ) );
	CUDA_CALL( cudaMemset( data.g_MinMax, 0, sizeof(int4) ) );
	runAlimitBlimitGridKernel<<<nBlocks,nThreads>>>( data );
	CUDA_CALL( cudaEventRecord( evtEnd[4], 0 ) );

	int4 MinMax;
	CUDA_CALL( cudaMemcpy( &MinMax, data.g_MinMax, sizeof(int4), cudaMemcpyDeviceToHost ) );
	cudaDeviceSynchronize();

	if( MinMax.x )
	    Imin = dp.iMin = max( dp.iMin-1, 2 );

	if( MinMax.y )
	    Imax = dp.iMax = min( dp.iMax+1, dp.nI-1 );

	if( MinMax.z )
	    Jmin = dp.jMin = max( dp.jMin-32, 2 );

	if( MinMax.w )
	    Jmax = dp.jMax = min( dp.jMax+1, dp.nJ-1 );

	float _dur;
	for( int j = 0; j < 5; j++ ) {
		cudaEventElapsedTime( &_dur, evtStart[j], evtEnd[j]);
		dur[j] += _dur;
	}


	copied = false;

	return 0;
}


