#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

// Two utility functions from "common"
extern int output_device_info(cl_device_id);
extern char *OpenCL_err_code(cl_int);

#define OpenCL_VWarn(...)                                                      \
  CCTK_VWarn(opencl_error_code, __LINE__, __FILE__, CCTK_THORNSTRING,          \
             __VA_ARGS__)
#define OpenCL_VInfo(...) CCTK_VInfo(CCTK_THORNSTRING, __VA_ARGS__)

//-------------------------------------------------------------------------------
void OpenCL_PrintInfo(CCTK_ARGUMENTS) {
  int err;                           // error code returned from OpenCL calls
  cl_uint numPlatforms;              // number of platforms
  cl_uint numDevices;                // number of devices for a given platform
  cl_char vendor_name[1024] = {0};   // string to hold vendor name
  cl_char platform_name[1024] = {0}; // string to hold name of platform

  CCTK_INT opencl_error_code = 0;

  // How many platforms are visible to this program?
  err = clGetPlatformIDs(0, NULL, &numPlatforms);
  if (err != CL_SUCCESS || numPlatforms <= 0) {
    OpenCL_VWarn("Error: Failed to find the platform! %s",
                 OpenCL_err_code(err));
    return;
  }

  OpenCL_VInfo("%d Platforms found", numPlatforms);

  // find the IDs of the platforms
  cl_platform_id platformId[numPlatforms]; // array of platform IDs
  err = clGetPlatformIDs(numPlatforms, platformId, NULL);
  if (err != CL_SUCCESS) {
    OpenCL_VWarn("Error: Failed to find  platform IDs! %s",
                 OpenCL_err_code(err));
    return;
  }

  for (unsigned int i = 0; i < numPlatforms; i++) {
    err = clGetPlatformInfo(*(platformId + i), CL_PLATFORM_NAME,
                            sizeof(platform_name), &platform_name, NULL);
    if (err != CL_SUCCESS) {
      OpenCL_VWarn("Error: Failed to access platform name!");
      return;
    }
    OpenCL_VInfo("Platform %d is %s", i + 1, platform_name);

    err = clGetPlatformInfo(*(platformId + i), CL_PLATFORM_VENDOR,
                            sizeof(vendor_name), &vendor_name, NULL);
    if (err != CL_SUCCESS) {
      OpenCL_VWarn("Error: Failed to access platform vendor!");
      return;
    }
    OpenCL_VInfo("  from %s", vendor_name);

    err = clGetDeviceIDs(*(platformId + i), CL_DEVICE_TYPE_ALL, 0, NULL,
                         &numDevices);
    if (err != CL_SUCCESS) {
      OpenCL_VWarn("Error: Failed to find number of devices! %s",
                   OpenCL_err_code(err));
      return;
    }

    OpenCL_VInfo("  %d devices found", numDevices);

    cl_device_id device_id[numDevices]; // compute device id
    err = clGetDeviceIDs(*(platformId + i), CL_DEVICE_TYPE_ALL, numDevices,
                         device_id, NULL);
    if (err != CL_SUCCESS) {
      OpenCL_VWarn("Error: Failed to find defice IDs %s", OpenCL_err_code(err));
      return;
    }

    for (unsigned int j = 0; j < numDevices; j++) {
      err = output_device_info(*(device_id + j));
    }
  }
  return;
}
