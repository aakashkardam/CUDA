#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//#include "timer.h"
//include "check.h"
#include <chrono>

#define nthreads 32
#define SOFTENING 1e-9f

/*
 * Each body contains x, y, and z coordinate positions,
 * as well as velocities in the x, y, and z directions.
 */

typedef struct { float x, y, z, vx, vy, vz; } Body;

void randomizeBodies(float *data, int n) { // this should remain a host function
  for (int i = 0; i < n; i++) {
    data[i] = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
  }
}

/*
 * This function calculates the gravitational impact of all bodies in the system
 * on all others, but does not update their positions.
 */

__global__
void bodyForce(Body *p, float dt, int n) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i<n) 
  {
    float Fx = 0.0f; float Fy = 0.0f; float Fz = 0.0f;
    for(int tl = 0; tl <gridDim.x; tl++)
    {
      __shared__ float3 position_in_shared_mem[nthreads]; // usign shared memory on gpu
      float position_in_x = p[tl * blockDim.x + threadIdx.x].x;
      float position_in_y = p[tl * blockDim.x + threadIdx.x].y;
      float position_in_z = p[tl * blockDim.x + threadIdx.x].z;
      position_in_shared_mem[threadIdx.x] = make_float3(position_in_x, position_in_y, position_in_z);
      __syncthreads();
      for(int j=0; j<nthreads; j++)
      {
        float dx = position_in_shared_mem[j].x - p[i].x;
	float dy = position_in_shared_mem[j].y - p[i].y;
	float dz = position_in_shared_mem[j].z - p[i].z;
	float distSqr = dx*dx + dy*dy + dz*dz + SOFTENING;
        float invDist = rsqrtf(distSqr);
        float invDist3 = invDist * invDist * invDist;

	Fx += dx * invDist3; Fy += dy * invDist3; Fz += dz * invDist3;
      }
      __syncthreads();
     }
     p[i].vx += dt*Fx; p[i].vy += dt*Fy; p[i].vz += dt*Fz;
  }
}

int main(const int argc, const char** argv) {


  int nBodies = 2<<11;
  int salt = 0;
  if (argc > 1) nBodies = 2<<atoi(argv[1]);

  /*
   * This salt is for assessment reasons. Tampering with it will result in automatic failure.
   */

  if (argc > 2) salt = atoi(argv[2]);

  const float dt = 0.01f; // time step
  const int nIters = 10;  // simulation iterations

  int nblock = (nBodies + nthreads - 1)/nthreads;

  cudaDeviceProp props;
  int deviceId;
  cudaGetDevice(&deviceId); 
  cudaGetDeviceProperties(&props, deviceId);
  int computeCapabilityMajor = props.major;
  int computeCapabilityMinor = props.minor;
  int multiProcessorCount = props.multiProcessorCount;
  int warpSize = props.warpSize;
  int bytes = nBodies * sizeof(Body);
  float *buf;

  //buf = (float *)malloc(bytes);
  cudaMallocManaged(&buf,bytes);
  cudaMemPrefetchAsync(buf,bytes,deviceId);

  Body *p = (Body*)buf;


  randomizeBodies(buf, 6 * nBodies); // Init pos / vel data

  double totalTime = 0.0;


  printf("Device ID: %d\nNumber of SMs: %d\nCompute Capability Major: %d\nCompute Capability Minor: %d\nWarp Size: %d\n", deviceId, multiProcessorCount, computeCapabilityMajor, computeCapabilityMinor, warpSize);
  /*
   * This simulation will run for 10 cycles of time, calculating gravitational
   * interaction amongst bodies, and adjusting their positions to reflect.
   */
 
  /*******************************************************************/
  for (int iter = 0; iter < nIters; iter++) {
 auto start=std::chrono::high_resolution_clock::now();	  
    
 bodyForce<<<nblock,nthreads>>>(p,dt,nBodies);
  /*
   * This position integration cannot occur until this round of `bodyForce` has completed.
   * Also, the next round of `bodyForce` cannot begin until the integration is complete.
   */

    cudaDeviceSynchronize();
    for (int i = 0 ; i < nBodies; i++) { // integrate position
      p[i].x += p[i].vx*dt;
      p[i].y += p[i].vy*dt;
      p[i].z += p[i].vz*dt;
    }

  /*******************************************************************/
    //const double tElapsed = GetTimer() / 1000.0;
    //totalTime += tElapsed;
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> tElapsed = finish - start;
    totalTime += tElapsed.count();
  }

  cudaDeviceSynchronize();
  double avgTime = totalTime / (double)(nIters);
  float billionsOfOpsPerSecond = 1e-9 * nBodies * nBodies / avgTime;

#ifdef ASSESS
  checkPerformance(buf, billionsOfOpsPerSecond, salt);
#else
  //checkAccuracy(buf, nBodies);
  printf("%d Bodies: average %0.3f Billion Interactions / second\n", nBodies, billionsOfOpsPerSecond);
  salt += 1;
#endif
  /*******************************************************************/


  //free(buf);
  cudaFree(buf);
}

