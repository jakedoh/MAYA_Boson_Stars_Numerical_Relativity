 /*@@
   @header    ioHDF5UtilGH.h
   @date      Fri Oct 6 2000
   @author    Thomas Radke
   @desc
              The GH extensions structure for thorn IOHDF5Util.
   @enddesc
   @version   $Header$
 @@*/

#ifndef _IOUTILHDF5_IOUTILHDF5GH_H_
#define _IOUTILHDF5_IOUTILHDF5GH_H_ 1

#define H5_USE_16_API 1
#include <hdf5.h>
#include "CactusBase/IOUtil/src/ioutil_Utils.h"


/********************************************************************
 ********************    Macro Definitions   ************************
 ********************************************************************/
/* check return code of HDF5 call and print a warning in case of an error */
#define HDF5_ERROR(fn_call)                                                   \
        {                                                                     \
          hid_t _error_code = fn_call;                                        \
                                                                              \
                                                                              \
          if (_error_code < 0)                                                \
          {                                                                   \
            CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,              \
                        "HDF5 call '%s' returned error code %d",              \
                        #fn_call, (int)_error_code);                          \
          }                                                                   \
        }

/* append an attribute to a given file object */
#define WRITE_ATTRIBUTE(name, value, object, ioHDF5UtilGH, dim, datatype)     \
        {                                                                     \
          int _len;                                                           \
          hid_t _attr, _dataspace;                                            \
          /* this union is only there to fix a bug in the HDF5 API */         \
          union                                                               \
          {                                                                   \
            const void *const_val;                                            \
            void *non_const_val;                                              \
          } _v;                                                               \
          hsize_t _dim = dim;                                                 \
                                                                              \
                                                                              \
          _v.const_val = value;                                               \
          if (H5Tget_class (datatype) == H5T_STRING)                          \
          {                                                                   \
            _len = strlen ((const char *) _v.const_val);                      \
            if (_len == 0)      /* HDF5 doesn't like zero-len strings */      \
            {                                                                 \
              _len++;                                                         \
            }                                                                 \
            HDF5_ERROR (H5Tset_size (datatype, _len));                        \
          }                                                                   \
          if (_dim > 0)                                                       \
          {                                                                   \
            _dataspace = ioHDF5UtilGH->array_dataspace;                       \
            HDF5_ERROR (H5Sset_extent_simple (_dataspace, 1, &_dim, NULL));   \
          }                                                                   \
          else                                                                \
          {                                                                   \
            _dataspace = ioHDF5UtilGH->scalar_dataspace;                      \
          }                                                                   \
          HDF5_ERROR (_attr = H5Acreate (object, name, datatype, _dataspace,  \
                                         H5P_DEFAULT));                       \
          HDF5_ERROR (H5Awrite (_attr, datatype, _v.non_const_val));          \
          HDF5_ERROR (H5Aclose (_attr));                                      \
        }

/* read an attribute in a given file object */
#define READ_ATTRIBUTE(object, attrname, requested_type, buffer)              \
        {                                                                     \
          hid_t _attr, _attrtype;                                             \
          hsize_t _asize = 0;                                                 \
                                                                              \
                                                                              \
          _attr = H5Aopen_name (object, attrname);                            \
          if (_attr >= 0)                                                     \
          {                                                                   \
            if (requested_type == H5T_C_S1)                                   \
            {                                                                 \
              HDF5_ERROR (_attrtype = H5Aget_type (_attr));                   \
              HDF5_ERROR (_asize = H5Tget_size (_attrtype));                  \
              if (_asize + 1 >= sizeof (buffer))                              \
              {                                                               \
                CCTK_WARN (1, "Can't read " attrname " attribute (too long)");\
              }                                                               \
            }                                                                 \
            else                                                              \
            {                                                                 \
              _attrtype = requested_type;                                     \
            }                                                                 \
            if (H5Aread (_attr, _attrtype, buffer) < 0)                       \
            {                                                                 \
              CCTK_WARN (1, "Can't read '" attrname "' attribute");           \
            }                                                                 \
            if (requested_type == H5T_C_S1)                                   \
            {                                                                 \
              ((char *) buffer)[_asize] = 0;                                  \
              HDF5_ERROR (H5Tclose (_attrtype));                              \
            }                                                                 \
            HDF5_ERROR (H5Aclose (_attr));                                    \
          }                                                                   \
          else                                                                \
          {                                                                   \
            CCTK_WARN (1, "Can't find '" attrname "' attribute");             \
          }                                                                   \
        }


/*** Define the different datatypes used for HDF5 I/O
     NOTE: the complex datatype is defined dynamically at runtime in Startup.c
 ***/
/* byte type is easy */
#define HDF5_BYTE   H5T_NATIVE_UCHAR

/* floating point types are architecture-independent,
   ie. a float has always 4 bytes, and a double has 8 bytes
     HDF5_REAL  is used for storing reals of the generic CCTK_REAL type
     HDF5_REALn is used to explicitely store n-byte reals */
#ifdef  HAVE_CCTK_REAL4
#define HDF5_REAL4  H5T_NATIVE_FLOAT
#endif
#ifdef  HAVE_CCTK_REAL8
#define HDF5_REAL8  H5T_NATIVE_DOUBLE
#endif
#ifdef  HAVE_CCTK_REAL16
#define HDF5_REAL16 (sizeof (CCTK_REAL16) == sizeof (long double) ?           \
                     H5T_NATIVE_LDOUBLE : -1)
#endif


#ifdef  CCTK_REAL_PRECISION_16
#define HDF5_REAL   HDF5_REAL16
#elif   CCTK_REAL_PRECISION_8
#define HDF5_REAL   HDF5_REAL8
#elif   CCTK_REAL_PRECISION_4
#define HDF5_REAL   HDF5_REAL4
#endif


/* integer types are architecture-dependent:
     HDF5_INT  is used for communicating integers of the generic CCTK_INT type
     HDF5_INTn is used to explicitely communicate n-byte integers */
#ifdef  HAVE_CCTK_INT8
#define HDF5_INT8   (sizeof (CCTK_INT8) == sizeof (int) ? H5T_NATIVE_INT :    \
                     sizeof (CCTK_INT8) == sizeof (long) ? H5T_NATIVE_LONG :  \
                     sizeof (CCTK_INT8) == sizeof (long long) ?               \
                     H5T_NATIVE_LLONG : -1)
#endif

#ifdef  HAVE_CCTK_INT4
#define HDF5_INT4   (sizeof (CCTK_INT4) == sizeof (int) ? H5T_NATIVE_INT :    \
                     sizeof (CCTK_INT4) == sizeof (short) ?                   \
                     H5T_NATIVE_SHORT : -1)
#endif

#ifdef  HAVE_CCTK_INT2
#define HDF5_INT2   (sizeof (CCTK_INT2) == sizeof (short) ?                   \
                     H5T_NATIVE_SHORT : -1)
#endif

#ifdef  HAVE_CCTK_INT1
#define HDF5_INT1   H5T_NATIVE_SCHAR
#endif

#ifdef  CCTK_INTEGER_PRECISION_8
#define HDF5_INT    HDF5_INT8
#elif   CCTK_INTEGER_PRECISION_4
#define HDF5_INT    HDF5_INT4
#elif   CCTK_INTEGER_PRECISION_2
#define HDF5_INT    HDF5_INT2
#elif   CCTK_INTEGER_PRECISION_1
#define HDF5_INT    HDF5_INT1
#endif


/* constants to index timers within the timers array */
#define CP_PARAMETERS_TIMER   0
#define CP_VARIABLES_TIMER    1
#define CP_TOTAL_TIMER        2
#define RECOVERY_TIMER        3
#define IOHDF5_NUM_TIMERS     4

/* names of the groups that hold global attributes and parameters */
#define GLOBAL_ATTRIBUTES_GROUP "Global Attributes"
#define CACTUS_PARAMETERS_GROUP "Cactus Parameters"
#define ALL_PARAMETERS          "All Parameters"


#ifdef __cplusplus
extern "C"
{
#endif


/* structure describing a given recovery file */
typedef struct
{
  /* flag indicating valid file info */
  int is_HDF5_file;

  /* HDF5 file handle */
  hid_t file;

  /* complete file name for recovery */
  char *filename;

  /* number of total and I/O processors, the associated IO processor,
     flag telling whether data was written chunked or unchunked */
  int nprocs, ioproc_every, ioproc, unchunked;

  /* whether file contains a Cactus version ID (used to distinguish checkpoint
     files with old/new timelevel naming scheme) */
  int has_version;

} fileinfo_t;


/* IOHDF5Util's GH extension structure */
typedef struct
{
  /* predefined dataspaces for writing scalar and array attributes */
  hid_t scalar_dataspace, array_dataspace;

  /* predefined datatype for writing CCTK_COMPLEX types */
  hid_t HDF5_COMPLEX, HDF5_COMPLEX8, HDF5_COMPLEX16, HDF5_COMPLEX32;

  /* predefined datatype for writing C string string attributes */
  hid_t HDF5_STRING;

} ioHDF5UtilGH;


/* exported functions */
hid_t IOHDF5Util_DataType (const ioHDF5UtilGH *myGH, int cctk_type);
void IOHDF5Util_DumpParameters (const cGH *GH, int all, hid_t group);
void IOHDF5Util_DumpGHExtensions (const cGH *GH, hid_t group);
int  IOHDF5Util_DumpGH (const cGH *GH, const int *timers, hid_t file);
void IOHDF5Util_DumpCommonAttributes (const cGH *GH, const ioRequest *request,
                                      hid_t dataset);
int IOHDF5Util_DumpVar (const cGH *GH, const ioRequest *request, hid_t file);

int IOHDF5Util_RecoverParameters (const fileinfo_t *filenfo);
int IOHDF5Util_RecoverGHextensions (cGH *GH, const fileinfo_t *filenfo);
int IOHDF5Util_RecoverVariables (cGH *GH, const fileinfo_t *filenfo);

#ifdef __cplusplus
} // extern "C"
#endif

#endif  /* _IOUTILHDF5_IOUTILHDF5GH_H_ */
