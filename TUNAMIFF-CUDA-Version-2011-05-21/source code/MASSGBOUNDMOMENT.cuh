/*
 *      Tohoku University's Numerical-Analysis Model
 *          for Investigation of tsunami
 *           Far-field Tsunami version 
 */

#ifndef TF_KERNELS_H
#define TF_KERNELS_H

/*
Flow of TUNAMI FF simulation main program

Input of Water Depth and Initial profile
Initial condition: Still water level
Check of the area of computation
Equation of continuity
Open Sea boundary condition
Equation of Motion
Check of the area of computation
*/

__global__ void runMASSKernel( KernelData data );
__global__ void runGBOUNDKernel( KernelData data );
__global__ void runMOMENTKernel( KernelData data );
__global__ void runMOMENTBOUNDKernel( KernelData data );
__global__ void runAlimitBlimitGridKernel( KernelData data );

#endif /* TF_KERNELS_H */
/*
 * TUNAMI_FF
 */
