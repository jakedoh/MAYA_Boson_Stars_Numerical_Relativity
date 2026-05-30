//------------------------------------------------------------------------------
//
// Name:     OpenCL_err_code()
//
// Purpose:  Function to output descriptions of errors for an input error code
// License:  public domain
//
//
// RETURN:   echoes the input error code
//
// HISTORY:  Written by Tim Mattson, June 2010
//
//------------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

char *OpenCL_err_code(cl_int err_in);
char *OpenCL_err_code(cl_int err_in) {
  switch (err_in) {
  case CL_INVALID_PLATFORM:
    return ("CL_INVALID_PLATFORM");
    break;
  case CL_INVALID_DEVICE_TYPE:
    return ("CL_INVALID_DEVICE_TYPE");
    break;
  case CL_INVALID_CONTEXT:
    return ("CL_INVALID_CONTEXT");
    break;
  case CL_INVALID_DEVICE:
    return ("CL_INVALID_DEVICE");
    break;
  case CL_INVALID_VALUE:
    return ("CL_INVALID_VALUE");
    break;
  case CL_INVALID_QUEUE_PROPERTIES:
    return ("CL_INVALID_QUEUE_PROPERTIES");
    break;
  case CL_OUT_OF_RESOURCES:
    return ("CL_OUT_OF_RESOURCES");
    break;
  case CL_INVALID_PROGRAM_EXECUTABLE:
    return ("CL_INVALID_PROGRAM_EXECUTABLE");
    break;
  case CL_INVALID_KERNEL:
    return ("CL_INVALID_KERNEL");
    break;
  case CL_INVALID_KERNEL_ARGS:
    return ("CL_INVALID_KERNEL_ARGS");
    break;
  case CL_INVALID_WORK_DIMENSION:
    return ("CL_INVALID_WORK_DIMENSION");
    break;
  case CL_INVALID_GLOBAL_OFFSET:
    return ("CL_INVALID_GLOBAL_OFFSET");
    break;
  case CL_INVALID_WORK_GROUP_SIZE:
    return ("CL_INVALID_WORK_GROUP_SIZE");
    break;
  case CL_INVALID_WORK_ITEM_SIZE:
    return ("CL_INVALID_WORK_ITEM_SIZE");
    break;
  case CL_INVALID_IMAGE_SIZE:
    return ("CL_INVALID_IMAGE_SIZE");
    break;
  case CL_INVALID_EVENT_WAIT_LIST:
    return ("CL_INVALID_EVENT_WAIT_LIST");
    break;
  case CL_INVALID_MEM_OBJECT:
    return ("CL_INVALID_MEM_OBJECT");
    break;
  case CL_MEM_COPY_OVERLAP:
    return ("CL_MEM_COPY_OVERLAP");
    break;
  case CL_MEM_OBJECT_ALLOCATION_FAILURE:
    return ("CL_MEM_OBJECT_ALLOCATION_FAILURE");
    break;
  case CL_OUT_OF_HOST_MEMORY:
    return ("CL_OUT_OF_HOST_MEMORY");
    break;
  case CL_PROFILING_INFO_NOT_AVAILABLE:
    return ("CL_PROFILING_INFO_NOT_AVAILABLE");
    break;
  case CL_INVALID_EVENT:
    return ("CL_INVALID_EVENT");
    break;
  default:
    return ("unknown error.");
    break;
  }
  return "";
}
