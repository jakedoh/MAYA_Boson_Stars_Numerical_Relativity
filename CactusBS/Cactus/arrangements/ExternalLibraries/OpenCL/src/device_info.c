//------------------------------------------------------------------------------
//
// Name:     device_info()
//
// Purpose:  Function to output key parameters about the input OpenCL device.
// License:  public domain
//
//
// RETURN:   The OCL_SUCESS or the error code from one of the OCL function
//           calls internal to this function
//
// HISTORY:  Written by Tim Mattson, June 2010
//
//------------------------------------------------------------------------------
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <stdlib.h>
#include <string.h>

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

//  define VERBOSE if you want to print info about work groups sizes
//#define  VERBOSE 1
#ifdef VERBOSE
extern char *OpenCL_err_code(cl_int);
#endif

#define OpenCL_VWarn(...)                                                      \
  CCTK_VWarn(opencl_error_code, __LINE__, __FILE__, CCTK_THORNSTRING,          \
             __VA_ARGS__)
#define OpenCL_VInfo(...) CCTK_VInfo(CCTK_THORNSTRING, __VA_ARGS__)

int output_device_info(cl_device_id device_id);
int output_device_info(cl_device_id device_id) {
  int err; // error code returned from OpenCL calls
  cl_device_type
      device_type;    // Parameter defining the type of the compute device
  cl_uint comp_units; // the max number of compute units on a device
  cl_char vendor_name[1024] = {
      0}; // string to hold vendor name for compute device
  cl_char device_name[1024] = {0}; // string to hold name of compute device
#ifdef VERBOSE
  cl_uint max_work_itm_dims;
  size_t max_wrkgrp_size;
#endif

  CCTK_INT opencl_error_code = 0;

  err = clGetDeviceInfo(device_id, CL_DEVICE_NAME, sizeof(device_name),
                        &device_name, NULL);
  if (err != CL_SUCCESS) {
    OpenCL_VWarn("Error: Failed to access device name!");
    return EXIT_FAILURE;
  }
  OpenCL_VInfo("  Device: %s", device_name);

  err = clGetDeviceInfo(device_id, CL_DEVICE_TYPE, sizeof(device_type),
                        &device_type, NULL);
  if (err != CL_SUCCESS) {
    OpenCL_VWarn("Error: Failed to access device type information!");
    return EXIT_FAILURE;
  }

  err = clGetDeviceInfo(device_id, CL_DEVICE_VENDOR, sizeof(vendor_name),
                        &vendor_name, NULL);
  if (err != CL_SUCCESS) {
    OpenCL_VWarn("Error: Failed to access device vendor name!");
    return EXIT_FAILURE;
  }
  if (device_type == CL_DEVICE_TYPE_GPU)
    OpenCL_VInfo("    GPU from %s", vendor_name);
  else if (device_type == CL_DEVICE_TYPE_CPU)
    OpenCL_VInfo("    CPU from %s", vendor_name);
  else
    OpenCL_VInfo("    non CPU or GPU processor from %s", vendor_name);

  err = clGetDeviceInfo(device_id, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint),
                        &comp_units, NULL);
  if (err != CL_SUCCESS) {
    OpenCL_VWarn("Error: Failed to access device number of compute units !");
    return EXIT_FAILURE;
  }
  OpenCL_VInfo("      with a max of %d compute units", comp_units);

#ifdef VERBOSE
  //
  // Optionally print information about work group sizes
  //
  err = clGetDeviceInfo(device_id, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS,
                        sizeof(cl_uint), &max_work_itm_dims, NULL);
  if (err != CL_SUCCESS) {
    OpenCL_VWarn("Error: Failed to get device Info "
                 "(CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS)!",
                 OpenCL_err_code(err));
    return EXIT_FAILURE;
  }

  size_t max_loc_size[max_work_itm_dims];
  err = clGetDeviceInfo(device_id, CL_DEVICE_MAX_WORK_ITEM_SIZES,
                        max_work_itm_dims * sizeof(size_t), max_loc_size, NULL);
  if (err != CL_SUCCESS) {
    OpenCL_VWarn(
        "Error: Failed to get device Info (CL_DEVICE_MAX_WORK_ITEM_SIZES)!",
        OpenCL_err_code(err));
    return EXIT_FAILURE;
  }
  err = clGetDeviceInfo(device_id, CL_DEVICE_MAX_WORK_GROUP_SIZE,
                        sizeof(size_t), &max_wrkgrp_size, NULL);
  if (err != CL_SUCCESS) {
    OpenCL_VWarn(
        "Error: Failed to get device Info (CL_DEVICE_MAX_WORK_GROUP_SIZE)!",
        OpenCL_err_code(err));
    return EXIT_FAILURE;
  }
  OpenCL_VInfo("work group, work item information");
  OpenCL_VInfo("\n max loc dim ");
  for (int i = 0; i < max_work_itm_dims; i++)
    OpenCL_VInfo(" %d ", (int)(*(max_loc_size + i)));
  OpenCL_VInfo("\n");
  OpenCL_VInfo(" Max work group size = %d\n", (int)max_wrkgrp_size);

#endif

  return CL_SUCCESS;
}
