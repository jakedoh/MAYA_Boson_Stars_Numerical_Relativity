#include <cctk.h>
#include <cctk_Arguments.h>

#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
using namespace std;



namespace HelloWorldCUDA {
  
  // Check a return value, and if there is an error, output a
  // human-readable error message
  void check_error(cudaError_t cerr, char const *msg = "", ...)
#ifdef __GNUC__
    __attribute__((format (printf, 2, 3)))
#endif
    ;
  void check_error(cudaError_t cerr, char const *msg, ...)
  {
    if (cerr) {
      if (strcmp(msg, "")) {
        va_list ap;
        va_start(ap, msg);
        char *usermsg;
        vasprintf(&usermsg, msg, ap);
        va_end(ap);
        CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "CUDA Error %d: %s:\n%s",
                   int(cerr), cudaGetErrorString(cerr), usermsg);
        free(usermsg);
      } else {
        CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "CUDA Error %d: %s\n",
                   int(cerr), cudaGetErrorString(cerr));
      }
    }
  }



  extern "C"
  void HelloWorldCUDA_initial(CCTK_ARGUMENTS)
  {
    // Output device properties
    cudaDeviceProp prop;
    cudaError_t cerr = cudaGetDeviceProperties(&prop, 0);
    check_error(cerr, "Could not get device properties");
    printf("CUDA device properties (device %d):\n", 0);
    printf("   name:                        %s\n", prop.name);
    printf("   totalGlobalMem:              %zu\n", prop.totalGlobalMem);
    printf("   sharedMemPerBlock:           %zu\n", prop.sharedMemPerBlock);
    printf("   regsPerBlock:                %d\n", prop.regsPerBlock);
    printf("   warpSize:                    %d\n", prop.warpSize);
    printf("   memPitch:                    %zu\n", prop.memPitch);
    printf("   maxThreadsPerBlock:          %d\n", prop. maxThreadsPerBlock);
    printf("   maxThreadsDim:               %d %d %d\n", prop.maxThreadsDim[0], prop.maxThreadsDim[1], prop.maxThreadsDim[2]);
    printf("   maxGridSize:                 %d %d %d\n", prop.maxGridSize[0], prop.maxGridSize[1], prop.maxGridSize[2]);
    printf("   clockRate:                   %d\n", prop.clockRate);
    printf("   totalConstMem                %zu\n", prop.totalConstMem);
    printf("   major:                       %d\n", prop.major);
    printf("   minor:                       %d\n", prop.minor);
    printf("   textureAlignment:            %zu\n", prop.textureAlignment);
    printf("   deviceOverlap:               %d\n", prop.deviceOverlap);
    printf("   multiProcessorCount:         %d\n", prop.multiProcessorCount);
    printf("   kernelExecTimeoutEnabled:    %d\n", prop.kernelExecTimeoutEnabled);
    printf("   integrated:                  %d\n", prop.integrated);
    printf("   canMapHostMemory:            %d\n", prop.canMapHostMemory);
    printf("   computeMode:                 %d\n", prop.computeMode);
    printf("   maxTexture1D:                %d\n", prop.maxTexture1D);
    printf("   maxTexture2D:                %d %d\n", prop.maxTexture2D[0], prop.maxTexture2D[1]);
    printf("   maxTexture3D:                %d %d %d\n", prop.maxTexture3D[0], prop.maxTexture3D[1], prop.maxTexture3D[2]);
    printf("   maxTexture1DLayered:         %d %d\n", prop.maxTexture1DLayered[0], prop.maxTexture1DLayered[1]);
    printf("   maxTexture2DLayered:         %d %d %d\n", prop.maxTexture2DLayered[0], prop.maxTexture2DLayered[1], prop.maxTexture2DLayered[2]);
    printf("   surfaceAlignment:            %zu\n", prop.surfaceAlignment);
    printf("   concurrentKernels:           %d\n", prop.concurrentKernels);
    printf("   ECCEnabled:                  %d\n", prop.ECCEnabled);
    printf("   pciBusID:                    %d\n", prop.pciBusID);
    printf("   pciDeviceID:                 %d\n", prop.pciDeviceID);
    printf("   pciDomainID:                 %d\n", prop.pciDomainID);
    printf("   tccDriver:                   %d\n", prop.tccDriver);
    printf("   asyncEngineCount:            %d\n", prop.asyncEngineCount);
    printf("   unifiedAddressing:           %d\n", prop.unifiedAddressing);
    printf("   memoryClockRate:             %d\n", prop.memoryClockRate);
    printf("   memoryBusWidth:              %d\n", prop.memoryBusWidth);
    printf("   l2CacheSize:                 %d\n", prop.l2CacheSize);
    printf("   maxThreadsPerMultiProcessor: %d\n", prop.maxThreadsPerMultiProcessor);
  }
  
  
  
  __global__ void add(const int val1, const int val2, int* res)
  {
    const size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t j = blockIdx.y * blockDim.y + threadIdx.y;
    const size_t k = blockIdx.z * blockDim.z + threadIdx.z;
    if (i==0 && j==0 && k==0) {
      *res = val1 + val2;
    }
  }
  
  extern "C"
  void HelloWorldCUDA_evol(CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS;
    
    const int val1 = cctk_iteration;
    const int val2 = 3;
    int res = 42;                 // poison
    
    const dim3 blockDim(4, 4, 4);
    const dim3 gridDim((cctk_lsh[0] + blockDim.x - 1) / blockDim.x,
                       (cctk_lsh[1] + blockDim.y - 1) / blockDim.y,
                       (cctk_lsh[2] + blockDim.z - 1) / blockDim.z);
    add<<<gridDim, blockDim>>>(val1, val2, &res);
    
    CCTK_VInfo(CCTK_THORNSTRING, "CUDA says: %d + %d = %d", val1, val2, res);
  }
  
} // namespace HelloWorldCUDA
