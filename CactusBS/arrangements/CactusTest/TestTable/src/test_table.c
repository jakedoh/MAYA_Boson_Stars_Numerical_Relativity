/* test_table.c -- test driver for key/value tables API */
/* $Header: */

/*
** prototypes for functions local to this file
 * TestTable_driver
**
** CHECK_SET_GET_{INT,GENERIC_INT,REAL,COMPLEX}
** CHECK_SET_GET_{INT,REAL,GENERIC_REAL,COMPLEX}_ARRAY
**
** test_nonexistent_tables
** test_table_create_destroy
** test_set_get
** test_set_get_array
** test_iterators
** test_delete_table_entry
** test_set_create_from_string
** test_set_get_string
** test_set_get_pointers
** test_set_get_fpointers
** test_clone
** check_table_contents
** check_table_contents_ij
** check_table_contents_real1
** check_table_contents_real_e
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "util_Table.h"
#include "util_ErrorCodes.h"

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

/* FIXME: C99 defines <stdbool.h>, we should include that or a fake version */
typedef int bool;
#define true    1
#define false   0

/******************************************************************************/

/*
 * prototypes for functions in this file
 */

/* top-level driver */
int TestTable_driver(void);

/* functions private to this file */
static void test_nonexistent_tables(void);
static void test_table_create_destroy(void);
static void test_set_get(int handle, bool case_insensitive);
static void test_set_get_array(int handle);
static void test_iterators(int handle);
static void test_delete_table_entry(int handle, bool case_insensitive);
static int test_set_create_from_string(void);
static void test_set_get_string(int handle, bool case_insensitive);
static void test_set_get_pointers(int handle);
static void test_set_get_fpointers(int handle);
static void test_clone(int handle);
static void check_table_contents(int handle, bool order_up_flag);
static void check_table_contents_ij(int handle, int ihandle);
static void check_table_contents_real1(int handle, int ihandle);
static void check_table_contents_real_e(int handle, int ihandle);

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/*
 * This function is a top-level driver to test the key/value API.
 */
int TestTable_driver(void)
{
  printf("testing key/value table routines\n");
  fflush(stdout);

  test_nonexistent_tables();
  test_table_create_destroy();

    {
  const int handle = Util_TableCreate(UTIL_TABLE_FLAGS_DEFAULT);
  assert( handle >= 0 );
  assert( Util_TableSetInt(handle, 42, "foo/") == UTIL_ERROR_TABLE_BAD_KEY );

    {
  const int HANDLE = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);

  #ifdef UTIL_TABLE_DEBUG
  printf("--- printing handle=%d table (should be empty)\n", handle);
  print_table(handle);
  printf("--- about to test set/get on handle=%d table\n", handle);
  #endif
  test_set_get(handle, false);
  test_set_get(HANDLE, true);

  test_iterators(handle);
  test_delete_table_entry(handle, false);
  test_iterators(HANDLE);
  test_delete_table_entry(HANDLE, true);

  test_set_get_array(handle);
  test_set_get_array(HANDLE);

  test_set_get_string(handle, false);
  test_set_get_pointers(handle);
  test_set_get_fpointers(handle);

    {
  const int HANDLE2 = test_set_create_from_string();
  test_clone(HANDLE2);
  test_set_get_string(HANDLE2, true);

  printf("all ok!\n" );
    }
    }
    }

  return 0;
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/*
 * This macro tests set/get of explicitly-typed integral scalars.
 */
#define CHECK_SET_GET_INT(handle, type,                                 \
                          key_already_exists, case_insensitive,         \
                          set_fn, get_fn)                               \
    {                                                                   \
  type x = 42;                                                          \
                                                                        \
  printf(                                                               \
    "CHECK_GET_SET_INT(handle=%d, type=%s,\n"                           \
    "                  key_already_exists=%d, case_insensitive=%d,\n"   \
    "                  set_fn=%s, get_fn=%s)\n",                        \
         handle, #type,                                                 \
         key_already_exists, case_insensitive,                          \
         #set_fn, #get_fn);                                             \
  fflush(stdout);                                                       \
                                                                        \
  assert( set_fn(handle, x, "int_x") == key_already_exists );           \
                                                                        \
  x = 1;                                                                \
  assert( get_fn(handle, &x, "int_x") == 1 );                           \
  assert( x == 42 );                                                    \
                                                                        \
  if (case_insensitive)                                                 \
  {                                                                     \
    x = 2;                                                              \
    assert( get_fn(handle, &x, "Int_X") == 1 );                         \
    assert( x == 42 );                                                  \
  }                                                                     \
  else                                                                  \
  {                                                                     \
    assert( get_fn(handle, &x, "Int_X")                                 \
            == UTIL_ERROR_TABLE_NO_SUCH_KEY );                          \
  }                                                                     \
    }                                                           /* end macro */

/******************************************************************************/

/*
 * This macro tests set/get of generically-typed integral scalars.
 */
#define CHECK_SET_GET_GENERIC_INT(handle, type_code, type,              \
                                  key_already_exists, case_insensitive) \
    {                                                                   \
  type gx = 42;                                                         \
                                                                        \
  printf(                                                               \
    "CHECK_SET_GET_GENERIC_INT(handle=%d, type_code=%d, type=%s,\n"     \
    "                          key_already_exists=%d, case_insensitive=%d)\n", \
         handle, type_code, #type,                                      \
         key_already_exists, case_insensitive);                         \
  fflush(stdout);                                                       \
                                                                        \
  assert( Util_TableSetGeneric(handle,                                  \
                               type_code, (const void *) &gx,           \
                               "gint_x") == key_already_exists );       \
                                                                        \
  gx = 1;                                                               \
  assert( Util_TableGetGeneric(handle,                                  \
                               type_code, (void *) &gx,                 \
                               "gint_x") == 1 );                        \
  assert( gx == 42 );                                                   \
                                                                        \
  if (case_insensitive)                                                 \
  {                                                                     \
    gx = 2;                                                             \
    assert( Util_TableGetGeneric(handle,                                \
                                 type_code, (void *) &gx,               \
                                 "GInt_X") == 1 );                      \
    assert( gx == 42 );                                                 \
  }                                                                     \
  else                                                                  \
  {                                                                     \
    assert( Util_TableGetGeneric(handle,                                \
                                 type_code, (void *) &gx,               \
                                 "GInt_X")                              \
            == UTIL_ERROR_TABLE_NO_SUCH_KEY );                          \
  }                                                                     \
    }                                                           /* end macro */

/******************************************************************************/

/*
 * This macro tests set/get of explicitly-typed real scalars.
 */
#define CHECK_SET_GET_REAL(handle, type,                                \
                           key_already_exists, case_insensitive,        \
                           set_fn, get_fn)                              \
    {                                                                   \
  type y = 42.25;                                                       \
                                                                        \
  printf(                                                               \
    "CHECK_SET_GET_REAL(handle=%d, type=%s,\n"                          \
    "                   key_already_exists=%d, case_insensitive=%d,\n"  \
    "                   set_fn=%s, get_fn=%s)\n",                       \
         handle, #type,                                                 \
         key_already_exists, case_insensitive,                          \
         #set_fn, #get_fn);                                             \
  fflush(stdout);                                                       \
                                                                        \
  assert( set_fn(handle, y, "REAL_y") == key_already_exists );          \
                                                                        \
  y = 1.25;                                                             \
  assert( get_fn(handle, &y, "REAL_y") == 1 );                          \
  assert( y == 42.25 );                                                 \
                                                                        \
  if (case_insensitive)                                                 \
  {                                                                     \
    y = 1.5;                                                            \
    assert( get_fn(handle, &y, "real_y") == 1 );                        \
    assert( y == 42.25 );                                               \
  }                                                                     \
  else                                                                  \
  {                                                                     \
    assert( get_fn(handle, &y, "real_y")                                \
            == UTIL_ERROR_TABLE_NO_SUCH_KEY );                          \
  }                                                                     \
    }                                                           /* end macro */

/******************************************************************************/

/*
 * This macro tests set/get of explicitly-typed complex scalars.
 */
#define CHECK_SET_GET_COMPLEX(handle, type,                             \
                              key_already_exists, case_insensitive,     \
                              set_fn, get_fn)                           \
    {                                                                   \
  type z = CCTK_Cmplx(42.25, 105.5);                                    \
                                                                        \
  printf(                                                               \
    "CHECK_SET_GET_COMPLEX(handle=%d, type=%s,\n"                       \
    "                      key_already_exists=%d, case_insensitive=%d,\n" \
    "                      set_fn=%s, get_fn=%s)\n",                    \
         handle, #type,                                                 \
         key_already_exists, case_insensitive,                          \
         #set_fn, #get_fn);                                             \
  fflush(stdout);                                                       \
                                                                        \
  assert( set_fn(handle, z, "COMPlex_Z") == key_already_exists );       \
                                                                        \
  z = CCTK_Cmplx(1.25, -2.78);                                          \
  assert( get_fn(handle, &z, "COMPlex_Z") == 1 );                       \
  assert( CCTK_CmplxReal(z) == 42.25 );                                 \
  assert( CCTK_CmplxImag(z) == 105.5 );                                 \
                                                                        \
  if (case_insensitive)                                                 \
  {                                                                     \
    z = CCTK_Cmplx(1.5, -2.83);                                         \
    assert( get_fn(handle, &z, "COMPLEX_Z") == 1 );                     \
    assert( CCTK_CmplxReal(z) == 42.25 );                               \
    assert( CCTK_CmplxImag(z) == 105.5 );                               \
  }                                                                     \
  else                                                                  \
  {                                                                     \
    assert( get_fn(handle, &z, "COMPLEX_Z")                             \
            == UTIL_ERROR_TABLE_NO_SUCH_KEY );                          \
  }                                                                     \
    }                                                           /* end macro */

/******************************************************************************/

/*
 * This macro tests set/get of explicitly-typed integral arrays.
 */
#define CHECK_SET_GET_INT_ARRAY(handle, type,                           \
                                key_already_exists,                     \
                                set_fn, get_fn)                         \
    {                                                                   \
  static type xx[5] = { 41, 42, 48, 45, 47 };                           \
                                                                        \
  printf(                                                               \
    "CHECK_SET_GET_INT_ARRAY(handle=%d, type=%s,\n"                     \
    "                        key_already_exists=%d,\n"                  \
    "                        set_fn=%s, get_fn=%s)\n",                  \
         handle, #type,                                                 \
         key_already_exists,                                            \
         #set_fn, #get_fn);                                             \
  fflush(stdout);                                                       \
                                                                        \
  assert( set_fn(handle, 3, xx, "xx") == key_already_exists );          \
                                                                        \
  xx[0] = 14;  xx[1] = 15;  xx[2] = 16;  xx[3] = 17;  xx[4] = 19;       \
                                                                        \
  /* try to get 4 values, but only 3 were stored ==> should only get 3 */ \
  assert( get_fn(handle, 4, xx, "xx") == 3 );                           \
  assert( xx[0] == 41 );                                                \
  assert( xx[1] == 42 );                                                \
  assert( xx[2] == 48 );                                                \
  assert( xx[3] == 17 );                                                \
  assert( xx[4] == 19 );                                                \
    }                                                           /* end macro */

/******************************************************************************/

/*
 * This macro tests set/get of explicitly-typed real arrays.
 */
#define CHECK_SET_GET_REAL_ARRAY(handle, type,                          \
                                 key_already_exists,                    \
                                 set_fn, get_fn)                        \
    {                                                                   \
  static type yy[5] = { 41.25, 42.5, 48.0, 45.75, 47.125 };             \
                                                                        \
  printf(                                                               \
    "CHECK_SET_GET_REAL_ARRAY(handle=%d, type=%s,\n"                    \
    "                         key_already_exists=%d,\n"                 \
    "                         set_fn=%s, get_fn=%s)\n",                 \
         handle, #type,                                                 \
         key_already_exists,                                            \
         #set_fn, #get_fn);                                             \
  fflush(stdout);                                                       \
                                                                        \
  assert( set_fn(handle, 4, yy, "yy") == key_already_exists );          \
                                                                        \
  yy[0] = 14.0;  yy[1] = 15.5;  yy[2] = 16.0;                           \
  yy[3] = 17.5;  yy[4] = 19.5;                                          \
                                                                        \
  /* only get 3 of 4 stored values */                                   \
  assert( get_fn(handle, 3, yy, "yy") == 4 );                           \
  assert( yy[0] == 41.25 );                                             \
  assert( yy[1] == 42.5 );                                              \
  assert( yy[2] == 48.0 );                                              \
  assert( yy[3] == 17.5 );                                              \
  assert( yy[4] == 19.5 );                                              \
    }                                                           /* end macro */

/******************************************************************************/

/*
 * This macro tests set/get of generically-typed real arrays.
 */
#define CHECK_SET_GET_GENERIC_REAL_ARRAY(handle, type_code, type,       \
                                         key_already_exists)            \
    {                                                                   \
  static type gyy[5] = { 41.25, 42.5, 48.0, 45.75, 47.125 };            \
                                                                        \
  printf(                                                               \
    "CHECK_SET_GET_GENERIC_REAL_ARRAY(handle=%d, type_code=%d, type=%s,\n" \
    "                                 key_already_exists=%d)\n",        \
         handle, type_code, #type,                                      \
         key_already_exists);                                           \
  fflush(stdout);                                                       \
                                                                        \
  assert( Util_TableSetGenericArray(handle,                             \
                                    type_code, 4, (const void *) gyy,   \
                                    "gyy") == key_already_exists );     \
                                                                        \
  gyy[0] = 14.0;  gyy[1] = 15.5;  gyy[2] = 16.0;                        \
  gyy[3] = 17.5;  gyy[4] = 19.5;                                        \
                                                                        \
  /* only get 3 of 4 stored values */                                   \
  assert( Util_TableGetGenericArray(handle,                             \
                                    type_code, 3, gyy,                  \
                                    "gyy") == 4 );                      \
  assert( gyy[0] == 41.25 );                                            \
  assert( gyy[1] == 42.5 );                                             \
  assert( gyy[2] == 48.0 );                                             \
  assert( gyy[3] == 17.5 );                                             \
  assert( gyy[4] == 19.5 );                                             \
    }                                                           /* end macro */

/******************************************************************************/

/*
 * This macro tests set/get of explicitly-typed complex arrays.
 */
#define CHECK_SET_GET_COMPLEX_ARRAY(handle, type,                       \
                                    key_already_exists,                 \
                                    set_fn, get_fn)                     \
    {                                                                   \
  type zz[5]                                                            \
    = { CCTK_Cmplx(3.5,1.25), CCTK_Cmplx(9.5,4.5), CCTK_Cmplx(0.5,8.0), \
        CCTK_Cmplx(5.0,5.5), CCTK_Cmplx(.5,7.25) };                     \
                                                                        \
  printf(                                                               \
    "CHECK_SET_GET_COMPLEX_ARRAY(handle=%d, type=%s,\n"                 \
    "                            key_already_exists=%d,\n"              \
    "                            set_fn=%s, get_fn=%s)\n",              \
         handle, #type,                                                 \
         key_already_exists,                                            \
         #set_fn, #get_fn);                                             \
  fflush(stdout);                                                       \
                                                                        \
  assert( set_fn(handle, 4, zz, "zz") == key_already_exists );          \
                                                                        \
  zz[0] = CCTK_Cmplx(10.25, 11.75);                                     \
  zz[1] = CCTK_Cmplx(-2.5, 3.5);                                        \
  zz[2] = CCTK_Cmplx(14.0, -8.5);                                       \
  zz[3] = CCTK_Cmplx(0.25, 8.875);                                      \
  zz[4] = CCTK_Cmplx(-0.25, -0.75);                                     \
                                                                        \
  /* only get 3 of 4 stored values */                                   \
  assert( get_fn(handle, 3, zz, "zz") == 4 );                           \
  assert( CCTK_CmplxReal(zz[0]) == 3.5 );      assert( CCTK_CmplxImag(zz[0]) == 1.25 );           \
  assert( CCTK_CmplxReal(zz[1]) == 9.5 );      assert( CCTK_CmplxImag(zz[1]) == 4.5 );            \
  assert( CCTK_CmplxReal(zz[2]) == 0.5 );      assert( CCTK_CmplxImag(zz[2]) == 8.0 );            \
  assert( CCTK_CmplxReal(zz[3]) == 0.25 );     assert( CCTK_CmplxImag(zz[3]) == 8.875 );          \
  assert( CCTK_CmplxReal(zz[4]) == -0.25 );    assert( CCTK_CmplxImag(zz[4]) == -0.75 );          \
    }                                                           /* end macro */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/*
 * This function tests that various operations on nonexistent tables
 * and iterators give error returns.
 */
static
  void test_nonexistent_tables(void)
{
  printf("test_nonexistent_tables()\n");
  fflush(stdout);

  assert( Util_TableDestroy(42) == UTIL_ERROR_BAD_HANDLE );
  assert( Util_TableQueryFlags(-42) == UTIL_ERROR_BAD_HANDLE );
  assert( Util_TableQueryNKeys(42) == UTIL_ERROR_BAD_HANDLE );
  assert( Util_TableDeleteKey(-1, "pickle") == UTIL_ERROR_BAD_HANDLE );
  assert( Util_TableSetInt(-1, 42, "fourty-two") == UTIL_ERROR_BAD_HANDLE );
  assert( Util_TableGetReal(-1, NULL, "something wierd")
          == UTIL_ERROR_BAD_HANDLE );

  assert( Util_TableItCreate(42) == UTIL_ERROR_BAD_HANDLE );
  assert( Util_TableItDestroy(42) == UTIL_ERROR_BAD_HANDLE );
  assert( Util_TableItQueryIsNull(42) == UTIL_ERROR_BAD_HANDLE );
  assert( Util_TableItQueryIsNonNull(42) == UTIL_ERROR_BAD_HANDLE );
  assert( Util_TableItQueryTableHandle(42) == UTIL_ERROR_BAD_HANDLE );
  assert( Util_TableItQueryKeyValueInfo(42,
                                        0, NULL,
                                        NULL, NULL) == UTIL_ERROR_BAD_HANDLE );
  assert( Util_TableItAdvance(42) == UTIL_ERROR_BAD_HANDLE );
  assert( Util_TableItResetToStart(42) == UTIL_ERROR_BAD_HANDLE );
  assert( Util_TableItSetToNull(42) == UTIL_ERROR_BAD_HANDLE );
}

/******************************************************************************/

/*
 * This function tests creation and destruction of tables.
 * It also tests
 * - querying flags words
 * - querying NKeys and MaxKeyLength for empty tables
 * - deleting keys from empty tables
 *
 * It assumes that no tables exist when the function is called,
 * and it eventually destroys all the tables it creates.
 *
 * Bugs:
 * Parts of this test are tied to the present implementation -- it
 * assumes a specific strategy for allocating handles.
 */
static
  void test_table_create_destroy(void)
{
  printf("test_table_create_destroy()\n");
  fflush(stdout);

  int handles[4];
  assert( (handles[0] = Util_TableCreate(UTIL_TABLE_FLAGS_DEFAULT)         ) >= 0 );
  assert( (handles[1] = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE)) > handles[0] );
  assert( (handles[2] = Util_TableCreate(UTIL_TABLE_FLAGS_DEFAULT)         ) > handles[1] );
  assert( (handles[3] = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE)) > handles[2] );
  assert( Util_TableQueryFlags(handles[0]) == UTIL_TABLE_FLAGS_DEFAULT );
  assert( Util_TableQueryFlags(handles[1]) == UTIL_TABLE_FLAGS_CASE_INSENSITIVE );
  assert( Util_TableQueryFlags(handles[2]) == UTIL_TABLE_FLAGS_DEFAULT );
  assert( Util_TableQueryFlags(handles[3]) == UTIL_TABLE_FLAGS_CASE_INSENSITIVE );

  assert( Util_TableDeleteKey(handles[3], "pickle") == UTIL_ERROR_TABLE_NO_SUCH_KEY );
  assert( Util_TableDeleteKey(handles[3], "Pickle") == UTIL_ERROR_TABLE_NO_SUCH_KEY );
  assert( Util_TableDeleteKey(handles[3], "PICKLE") == UTIL_ERROR_TABLE_NO_SUCH_KEY );

  assert( Util_TableDestroy(handles[2]) == 0 );

  assert( Util_TableCreate(0x43) == handles[2] );
  assert( Util_TableQueryFlags(handles[2]) == 0x43);

  assert( Util_TableDestroy(handles[1]) == 0 );
  assert( Util_TableQueryNKeys(handles[0]) == 0 );
  assert( Util_TableQueryMaxKeyLength(handles[1]) == UTIL_ERROR_BAD_HANDLE );
  assert( Util_TableQueryMaxKeyLength(handles[2]) == 0 );
  assert( Util_TableDeleteKey(handles[3], "pickle") == UTIL_ERROR_TABLE_NO_SUCH_KEY );

  assert( Util_TableDestroy(handles[1]) == UTIL_ERROR_BAD_HANDLE );

  assert( Util_TableDestroy(handles[0]) == 0 );
  assert( Util_TableDestroy(handles[2]) == 0 );
  assert( Util_TableDestroy(handles[3]) == 0 );
}

/******************************************************************************/

/*
 * This function tests set/get of various-typed scalars.
 * It also tests querying NKeys, MaxKeyLength, and some keys.
 *
 * It assumes the table is empty when this function is called;
 * it leaves entries in the table.
 */
static
  void test_set_get(int handle, bool case_insensitive)
{
  printf("test_set_get(handle=%d, case_insensitive=%d)\n",
         handle, (int)case_insensitive);
  fflush(stdout);

  /*
   * Note we put a test of a type that's guaranteed to be defined...
   * - at the *beginning* of each group of tests, so we can properly
   *   assert whether or not the key was already in table beforehand.
   * - at the *end* of each group of tests, so the final table contents
   *   are known independently of which types are and aren't defined.
   */

  /* integers */
  CHECK_SET_GET_INT(handle, CCTK_INT, 0, case_insensitive,
                    Util_TableSetInt, Util_TableGetInt);
  CHECK_SET_GET_INT(handle, CCTK_BYTE, 1, case_insensitive,
                    Util_TableSetByte, Util_TableGetByte);
  CHECK_SET_GET_INT(handle, CCTK_CHAR, 1, case_insensitive,
                    Util_TableSetChar, Util_TableGetChar);
  #ifdef HAVE_CCTK_INT1
  CHECK_SET_GET_INT(handle, CCTK_INT1, 1, case_insensitive,
                    Util_TableSetInt1, Util_TableGetInt1);
  #endif
  #ifdef HAVE_CCTK_INT2
  CHECK_SET_GET_INT(handle, CCTK_INT2, 1, case_insensitive,
                    Util_TableSetInt2, Util_TableGetInt2);
  #endif
  #ifdef HAVE_CCTK_INT4
  CHECK_SET_GET_INT(handle, CCTK_INT4, 1, case_insensitive,
                    Util_TableSetInt4, Util_TableGetInt4);
  #endif
  #ifdef HAVE_CCTK_INT8
  CHECK_SET_GET_INT(handle, CCTK_INT8, 1, case_insensitive,
                    Util_TableSetInt8, Util_TableGetInt8);
  #endif
  #ifdef HAVE_CCTK_INT16
  CHECK_SET_GET_INT(handle, CCTK_INT16, 1, case_insensitive,
                    Util_TableSetInt16, Util_TableGetInt16);
  #endif
  CHECK_SET_GET_INT(handle, CCTK_INT, 1, case_insensitive,
                    Util_TableSetInt, Util_TableGetInt);
  assert( Util_TableQueryNKeys(handle) == 1 );
  assert( Util_TableQueryMaxKeyLength(handle) == (int)strlen("int_x") );

  /* generic scalars which are actually integers */
  CHECK_SET_GET_GENERIC_INT(handle, CCTK_VARIABLE_INT, CCTK_INT,
                            0, case_insensitive);
  CHECK_SET_GET_GENERIC_INT(handle, CCTK_VARIABLE_CHAR, CCTK_CHAR,
                            1, case_insensitive);
  CHECK_SET_GET_GENERIC_INT(handle, CCTK_VARIABLE_BYTE, CCTK_BYTE,
                            1, case_insensitive);
  #ifdef HAVE_CCTK_INT1
  CHECK_SET_GET_GENERIC_INT(handle, CCTK_VARIABLE_INT1, CCTK_INT1,
                            1, case_insensitive);
  #endif
  #ifdef HAVE_CCTK_INT2
  CHECK_SET_GET_GENERIC_INT(handle, CCTK_VARIABLE_INT2, CCTK_INT2,
                            1, case_insensitive);
  #endif
  #ifdef HAVE_CCTK_INT4
  CHECK_SET_GET_GENERIC_INT(handle, CCTK_VARIABLE_INT4, CCTK_INT4,
                            1, case_insensitive);
  #endif
  #ifdef HAVE_CCTK_INT8
  CHECK_SET_GET_GENERIC_INT(handle, CCTK_VARIABLE_INT8, CCTK_INT8,
                            1, case_insensitive);
  #endif
  #ifdef HAVE_CCTK_INT16
  CHECK_SET_GET_GENERIC_INT(handle, CCTK_VARIABLE_INT16, CCTK_INT16,
                            1, case_insensitive);
  #endif
  CHECK_SET_GET_GENERIC_INT(handle, CCTK_VARIABLE_INT, CCTK_INT,
                            1, case_insensitive);
  assert( Util_TableQueryNKeys(handle) == 2 );
  assert( Util_TableQueryMaxKeyLength(handle) == (int)strlen("gint_x") );
  assert( Util_TableDeleteKey(handle, "gint_x") == 0 );
  assert( Util_TableQueryNKeys(handle) == 1 );

  /* complex numbers */
  CHECK_SET_GET_COMPLEX(handle, CCTK_COMPLEX, 0, case_insensitive,
                        Util_TableSetComplex, Util_TableGetComplex);
  #ifdef HAVE_CCTK_REAL4
  CHECK_SET_GET_COMPLEX(handle, CCTK_COMPLEX8, 1, case_insensitive,
                        Util_TableSetComplex8, Util_TableGetComplex8);
  #endif
  #ifdef HAVE_CCTK_REAL8
  CHECK_SET_GET_COMPLEX(handle, CCTK_COMPLEX16, 1, case_insensitive,
                        Util_TableSetComplex16, Util_TableGetComplex16);
  #endif
  #ifdef HAVE_CCTK_REAL16
  CHECK_SET_GET_COMPLEX(handle, CCTK_COMPLEX32, 1, case_insensitive,
                        Util_TableSetComplex32, Util_TableGetComplex32);
  #endif
  CHECK_SET_GET_COMPLEX(handle, CCTK_COMPLEX, 1, case_insensitive,
                        Util_TableSetComplex, Util_TableGetComplex);
  assert( Util_TableQueryNKeys(handle) == 2 );
  assert( Util_TableQueryMaxKeyLength(handle) == (int)strlen("COMPlex_Z") );

  /* reals */
  CHECK_SET_GET_REAL(handle, CCTK_REAL, 0, case_insensitive,
                     Util_TableSetReal, Util_TableGetReal);
  #ifdef HAVE_CCTK_REAL4
  CHECK_SET_GET_REAL(handle, CCTK_REAL4, 1, case_insensitive,
                     Util_TableSetReal4, Util_TableGetReal4);
  #endif
  #ifdef HAVE_CCTK_REAL8
  CHECK_SET_GET_REAL(handle, CCTK_REAL8, 1, case_insensitive,
                     Util_TableSetReal8, Util_TableGetReal8);
  #endif
  #ifdef HAVE_CCTK_REAL16
  CHECK_SET_GET_REAL(handle, CCTK_REAL16, 1, case_insensitive,
                     Util_TableSetReal16, Util_TableGetReal16);
  #endif
  CHECK_SET_GET_REAL(handle, CCTK_REAL, 1, case_insensitive,
                     Util_TableSetReal, Util_TableGetReal);
  assert( Util_TableQueryNKeys(handle) == 3 );
  assert( Util_TableQueryMaxKeyLength(handle) == (int)strlen("COMPlex_Z") );

    {
  CCTK_INT type_code, N_elements;
  assert( Util_TableQueryValueInfo(handle, &type_code, &N_elements, "COMPlex_Z")
          == 1 );
  assert( type_code == CCTK_VARIABLE_COMPLEX );
  assert( N_elements == 1 );

  assert( Util_TableQueryValueInfo(handle, &type_code, &N_elements, "pickle")
          == 0 );

  assert( Util_TableQueryValueInfo(handle, NULL, NULL, "int_x") == 1 );
  assert( Util_TableQueryValueInfo(handle, NULL, NULL, "Int_x")
          == (case_insensitive ? 1 : 0) );
  assert( Util_TableQueryValueInfo(handle, NULL, NULL, "real_y")
          == (case_insensitive ? 1 : 0) );
  assert( Util_TableQueryValueInfo(handle, NULL, NULL, "COMPLEX_Z")
          == (case_insensitive ? 1 : 0) );
    }
}

/******************************************************************************/

/*
 * This function tests set/get of various-typed arrays.
 *
 * It assumes the table is empty when this function is called;
 * it leaves entries in the table.
 */
static
  void test_set_get_array(int handle)
{
  printf("test_set_get_array(handle=%d)\n", handle);
  fflush(stdout);

  /* the comments of  test_set_get()  about test ordering, also apply here */

  /* integers */
  CHECK_SET_GET_INT_ARRAY(handle, CCTK_INT, 0,
                          Util_TableSetIntArray, Util_TableGetIntArray);
  CHECK_SET_GET_INT_ARRAY(handle, CCTK_CHAR, 1,
                          Util_TableSetCharArray, Util_TableGetCharArray);
  CHECK_SET_GET_INT_ARRAY(handle, CCTK_BYTE, 1,
                          Util_TableSetByteArray, Util_TableGetByteArray);
  #ifdef HAVE_CCTK_INT1
  CHECK_SET_GET_INT_ARRAY(handle, CCTK_INT1, 1,
                          Util_TableSetInt1Array, Util_TableGetInt1Array);
  #endif
  #ifdef HAVE_CCTK_INT2
  CHECK_SET_GET_INT_ARRAY(handle, CCTK_INT2, 1,
                          Util_TableSetInt2Array, Util_TableGetInt2Array);
  #endif
  #ifdef HAVE_CCTK_INT4
  CHECK_SET_GET_INT_ARRAY(handle, CCTK_INT4, 1,
                          Util_TableSetInt4Array, Util_TableGetInt4Array);
  #endif
  #ifdef HAVE_CCTK_INT8
  CHECK_SET_GET_INT_ARRAY(handle, CCTK_INT8, 1,
                          Util_TableSetInt8Array, Util_TableGetInt8Array);
  #endif
  CHECK_SET_GET_INT_ARRAY(handle, CCTK_INT, 1,
                          Util_TableSetIntArray, Util_TableGetIntArray);

  /* reals */
  CHECK_SET_GET_REAL_ARRAY(handle, CCTK_REAL, 0,
                           Util_TableSetRealArray, Util_TableGetRealArray);
  #ifdef HAVE_CCTK_REAL4
  CHECK_SET_GET_REAL_ARRAY(handle, CCTK_REAL4, 1,
                           Util_TableSetReal4Array, Util_TableGetReal4Array);
  #endif
  #ifdef HAVE_CCTK_REAL8
  CHECK_SET_GET_REAL_ARRAY(handle, CCTK_REAL8, 1,
                           Util_TableSetReal8Array, Util_TableGetReal8Array);
  #endif
  #ifdef HAVE_CCTK_REAL16
  CHECK_SET_GET_REAL_ARRAY(handle, CCTK_REAL16, 1,
                           Util_TableSetReal16Array, Util_TableGetReal16Array);
  #endif
  CHECK_SET_GET_REAL_ARRAY(handle, CCTK_REAL, 1,
                           Util_TableSetRealArray, Util_TableGetRealArray);

  /* generic arrays which are actually reals */
  CHECK_SET_GET_GENERIC_REAL_ARRAY(handle, CCTK_VARIABLE_REAL, CCTK_REAL, 0);
  #ifdef HAVE_CCTK_REAL4
  CHECK_SET_GET_GENERIC_REAL_ARRAY(handle, CCTK_VARIABLE_REAL4, CCTK_REAL4, 1);
  #endif
  #ifdef HAVE_CCTK_REAL8
  CHECK_SET_GET_GENERIC_REAL_ARRAY(handle, CCTK_VARIABLE_REAL8, CCTK_REAL8, 1);
  #endif
  #ifdef HAVE_CCTK_REAL16
  CHECK_SET_GET_GENERIC_REAL_ARRAY(handle,
                                   CCTK_VARIABLE_REAL16, CCTK_REAL16, 1);
  #endif
  CHECK_SET_GET_GENERIC_REAL_ARRAY(handle, CCTK_VARIABLE_REAL, CCTK_REAL, 1);

  /* complex numbers */
  CHECK_SET_GET_COMPLEX_ARRAY(handle, CCTK_COMPLEX, 0,
                              Util_TableSetComplexArray,
                              Util_TableGetComplexArray);
  #ifdef HAVE_CCTK_REAL4
  CHECK_SET_GET_COMPLEX_ARRAY(handle, CCTK_COMPLEX8, 1,
                              Util_TableSetComplex8Array,
                              Util_TableGetComplex8Array);
  #endif
  #ifdef HAVE_CCTK_REAL8
  CHECK_SET_GET_COMPLEX_ARRAY(handle, CCTK_COMPLEX16, 1,
                              Util_TableSetComplex16Array,
                              Util_TableGetComplex16Array);
  #endif
  #ifdef HAVE_CCTK_REAL16
  CHECK_SET_GET_COMPLEX_ARRAY(handle, CCTK_COMPLEX32, 1,
                              Util_TableSetComplex32Array,
                              Util_TableGetComplex32Array);
  #endif
  CHECK_SET_GET_COMPLEX_ARRAY(handle, CCTK_COMPLEX, 1,
                              Util_TableSetComplexArray,
                              Util_TableGetComplexArray);
}

/******************************************************************************/

/*
 * This function tests iterating through a table and resetting an iterator.
 * It assumes the initial table contents are those generated by
 * test_get_set() , namely 3 keys "REAL_y", "COMPlex_Z", "int_x".
 *
 * Bugs:
 * This test is tied to the present implementation -- it assumes a
 * specific ordering of table elements returned by an iterator.
 */
static
  void test_iterators(int handle)
{
  printf("test_iterators(handle=%d)\n", handle);
  fflush(stdout);

    {
  const int ihandle = Util_TableItCreate(handle);
  assert( ihandle >= 0 );
  assert( Util_TableItQueryTableHandle(ihandle) == handle );
  assert( Util_TableItQueryIsNonNull(ihandle) == 1);
  assert( Util_TableItQueryIsNull(ihandle) == 0);

  /* set up the key buffer */
    {
  const int max_key_length = Util_TableQueryMaxKeyLength(handle);
  assert( max_key_length == (int)strlen("COMPlex_Z") );
    {
  const int N_key_buffer = max_key_length + 1;
  char *const key_buffer = malloc(N_key_buffer);
  assert( key_buffer != NULL );

  /* walk the table to verify iterator traversal */
    {
  CCTK_INT type_code, N_elements;

  /* REAL_y */
  assert( Util_TableItQueryKeyValueInfo(ihandle,
                                        N_key_buffer, key_buffer,
                                        &type_code, &N_elements)
          == (int)strlen("REAL_y") );
  assert( strcmp(key_buffer, "REAL_y") == 0 );
  assert( type_code == CCTK_VARIABLE_REAL );
  assert( N_elements == 1 );

  assert( Util_TableItAdvance(ihandle) == 1 );

  /* COMPlex_Z */
  type_code = 123456;
  N_elements = 54321;
  assert( Util_TableItQueryKeyValueInfo(ihandle,
                                        N_key_buffer, key_buffer,
                                        &type_code, &N_elements)
          == (int)strlen("COMPlex_Z") );
  assert( strcmp(key_buffer, "COMPlex_Z") == 0 );
  assert( type_code == CCTK_VARIABLE_COMPLEX );
  assert( N_elements == 1 );

  assert( Util_TableItAdvance(ihandle) == 1 );

  /* clone the iterator and check the clone */
    {
  const int clone_ihandle = Util_TableItClone(ihandle);
  type_code = 123456;
  N_elements = 54321;
  assert( Util_TableItQueryKeyValueInfo(clone_ihandle,
                                        N_key_buffer, key_buffer,
                                        &type_code, &N_elements)
          == (int)strlen("int_x") );
  assert( strcmp(key_buffer, "int_x") == 0 );
  assert( type_code == CCTK_VARIABLE_INT );
  assert( N_elements == 1 );

  /* int_x */
  type_code = 123456;
  N_elements = 54321;
  assert( Util_TableItQueryKeyValueInfo(ihandle,
                                        N_key_buffer, key_buffer,
                                        &type_code, &N_elements)
          == (int)strlen("int_x") );
  assert( strcmp(key_buffer, "int_x") == 0 );
  assert( type_code == CCTK_VARIABLE_INT );
  assert( N_elements == 1 );

  /* advance past last table entry ==> "null-pointer" state */
  assert( Util_TableItAdvance(ihandle) == 0 );
  assert( Util_TableItQueryIsNull(ihandle) == 1);
  assert( Util_TableItQueryIsNonNull(ihandle) == 0);

  /* advance again ==> stays in "null-pointer" state */
  assert( Util_TableItAdvance(ihandle) == 0 );
  assert( Util_TableItQueryIsNull(ihandle) == 1);
  assert( Util_TableItQueryIsNonNull(ihandle) == 0);
  assert( Util_TableItQueryKeyValueInfo(ihandle,
                                        0, NULL,
                                        NULL, NULL)
          == UTIL_ERROR_TABLE_ITERATOR_IS_NULL );

  /* test reset to starting point */
  assert( Util_TableItResetToStart(ihandle) == 1 );
  assert( Util_TableItQueryIsNonNull(ihandle) == 1 );
  assert( Util_TableItQueryIsNull(ihandle) == 0 );

  /* COMPlex_Z */
  type_code = 123456;
  N_elements = 54321;
  assert( Util_TableItAdvance(ihandle) == 1 );
  assert( Util_TableItQueryKeyValueInfo(ihandle,
                                        N_key_buffer, key_buffer,
                                        &type_code, &N_elements)
          == (int)strlen("COMPlex_Z") );
  assert( strcmp(key_buffer, "COMPlex_Z") == 0 );
  assert( type_code == CCTK_VARIABLE_COMPLEX );
  assert( N_elements == 1 );

  /* test reset to "null-pointer" state */
  assert( Util_TableItSetToNull(ihandle) == 0 );
  assert( Util_TableItQueryIsNull(ihandle) == 1);
  assert( Util_TableItQueryIsNonNull(ihandle) == 0);

  /* test set to key "REAL_y" */
  assert( Util_TableItSetToKey(ihandle, "REAL_y") == 0 );
  assert( Util_TableItQueryIsNonNull(ihandle) == 1);
  assert( Util_TableItQueryIsNull(ihandle) == 0);
  assert( Util_TableItQueryKeyValueInfo(ihandle,
                                        N_key_buffer, key_buffer,
                                        &type_code, &N_elements)
          == (int)strlen("REAL_y") );
  assert( strcmp(key_buffer, "REAL_y") == 0 );
  assert( type_code == CCTK_VARIABLE_REAL );
  assert( N_elements == 1 );

  assert( Util_TableItDestroy(ihandle) == 0 );
  assert( Util_TableItDestroy(clone_ihandle) == 0 );

  free(key_buffer);
    }
    }
    }
    }
    }
}

/******************************************************************************/

/*
 * This function tests deleting table entries.
 * It assumes the initial table contents are those generated by
 * test_get_set() , namely 3 keys {"REAL_y", "COMPlex_Z", "int_x"}.
 *
 * Bugs:
 * This test is tied to the present implementation -- it assumes a
 * specific ordering of table elements returned by an iterator.
 */
static
  void test_delete_table_entry(int handle, bool case_insensitive)
{
  printf("test_delete_table_entry(handle=%d, case_insensitive=%d)\n",
         handle, (int)case_insensitive);
  fflush(stdout);

    {
  /* set up the key buffer */
  const int max_key_length = Util_TableQueryMaxKeyLength(handle);
  assert( max_key_length == (int)strlen("COMPlex_Z") );
    {
  const int N_key_buffer = max_key_length + 1;
  char *const key_buffer = malloc(N_key_buffer);
  assert( key_buffer != NULL );

  /*
   * delete the starting table entry "REAL_y"
   * (this is a special case in the implementation)
   */

  assert( Util_TableQueryNKeys(handle) == 3 );
  assert( Util_TableDeleteKey(handle,
                              case_insensitive ? "rEAL_y" : "REAL_y")
          == 0 );
  assert( Util_TableQueryNKeys(handle) == 2 );

  /* walk the table again to verify remaining keys {"COMPlex_Z", "int_x"} */
  assert( Util_TableQueryNKeys(handle) == 2 );
    {
  int ihandle = Util_TableItCreate(handle);
  assert( ihandle >= 0 );
  assert( Util_TableItQueryTableHandle(ihandle) == handle );
  assert( Util_TableItQueryIsNonNull(ihandle) == 1);
  assert( Util_TableItQueryIsNull(ihandle) == 0);

  /* COMPlex_Z */
    {
  CCTK_INT type_code = 123456;
  CCTK_INT N_elements = 54321;
  assert( Util_TableItQueryKeyValueInfo(ihandle,
                                        N_key_buffer, key_buffer,
                                        &type_code, &N_elements)
          == (int)strlen("COMPlex_Z") );
  assert( strcmp(key_buffer, "COMPlex_Z") == 0 );
  assert( type_code == CCTK_VARIABLE_COMPLEX );
  assert( N_elements == 1 );

  /* int_x */
  type_code = 123456;
  N_elements = 54321;
  assert( Util_TableItAdvance(ihandle) == 1 );
  assert( Util_TableItQueryKeyValueInfo(ihandle,
                                        N_key_buffer, key_buffer,
                                        &type_code, &N_elements)
          == (int)strlen("int_x") );
  assert( strcmp(key_buffer, "int_x") == 0 );
  assert( type_code == CCTK_VARIABLE_INT );
  assert( N_elements == 1 );

  /* advance past last table entry ==> "null-pointer" state */
  assert( Util_TableItAdvance(ihandle) == 0 );
  assert( Util_TableItQueryIsNull(ihandle) == 1);
  assert( Util_TableItQueryIsNonNull(ihandle) == 0);

  /* delete the last key "int_x" */
  assert( Util_TableDeleteKey(handle,
                              case_insensitive ? "INT_X" : "int_x")
          == 0 );

  /* walk the table again to verify remaining key {"COMPlex_Z"} */
  assert( Util_TableQueryNKeys(handle) == 1 );
    {
  int ihandle2 = Util_TableItCreate(handle);
  assert( ihandle2 >= 0 );
  assert( Util_TableItQueryTableHandle(ihandle2) == handle );
  assert( Util_TableItQueryIsNonNull(ihandle2) == 1);
  assert( Util_TableItQueryIsNull(ihandle2) == 0);

  /* COMPlex_Z */
  type_code = 123456;
  N_elements = 54321;
  assert( Util_TableItQueryKeyValueInfo(ihandle2,
                                        N_key_buffer, key_buffer,
                                        &type_code, &N_elements)
          == (int)strlen("COMPlex_Z") );
  assert( strcmp(key_buffer, "COMPlex_Z") == 0 );
  assert( type_code == CCTK_VARIABLE_COMPLEX );
  assert( N_elements == 1 );

  /* advance past last table entry ==> "null-pointer" state */
  assert( Util_TableItAdvance(ihandle2) == 0 );
  assert( Util_TableItQueryIsNull(ihandle2) == 1);
  assert( Util_TableItQueryIsNonNull(ihandle2) == 0);

  /* delete the last key "COMPlex_Z" */
  assert( Util_TableQueryNKeys(handle) == 1 );
  assert( Util_TableDeleteKey(handle,
                              case_insensitive ? "INT_X" : "int_x")
          == UTIL_ERROR_TABLE_NO_SUCH_KEY );
  assert( Util_TableQueryNKeys(handle) == 1 );
  assert( Util_TableDeleteKey(handle,
                              case_insensitive ? "compLEX_z" : "COMPlex_Z")
          == 0 );
  assert( Util_TableQueryNKeys(handle) == 0 );

    {
  /* check that table is indeed now empty */
  int ihandle3 = Util_TableItCreate(handle);
  assert( ihandle3 >= 0 );
  assert( Util_TableItQueryIsNull(ihandle3) == 1);
  assert( Util_TableItQueryIsNonNull(ihandle3) == 0);

  /* clean up our iterators */
  assert( Util_TableItDestroy(ihandle2) == 0 );
  assert( Util_TableItDestroy(42) == UTIL_ERROR_BAD_HANDLE );
  assert( Util_TableItDestroy(ihandle3) == 0 );
  assert( Util_TableItDestroy(ihandle) == 0 );
  free(key_buffer);
    }
    }
    }
    }
    }
    }
}

/******************************************************************************/

/*
 * This function tests
 *      Util_TableSetFromString()
 *      Util_TableCreateFromString()
 * It returns the handle of one of the newly-created tables.
 *
 * Bugs:
 * This test is tied to the present implementation -- it assumes a
 * specific ordering of table elements returned by an iterator.
 */
static
  int test_set_create_from_string(void)
{
  printf("test_set_create_from_string()\n");
  fflush(stdout);

    {
  /*
   * Test an empty string
   */
  const int handle = Util_TableCreateFromString("");
  assert( Util_TableQueryNKeys(handle) == 0 );

  /*
   * Test some error cases
   */
  assert( Util_TableSetFromString(handle, "foo" ) == UTIL_ERROR_BAD_INPUT );
  assert( Util_TableSetFromString(handle, "foo/" ) == UTIL_ERROR_BAD_INPUT );

  assert( Util_TableCreateFromString("foo" ) == UTIL_ERROR_BAD_INPUT );
  assert( Util_TableCreateFromString("foo/" ) == UTIL_ERROR_BAD_INPUT );
  assert( Util_TableCreateFromString("foo=" ) == UTIL_ERROR_BAD_INPUT );
  assert( Util_TableCreateFromString("foo/=12" ) == UTIL_ERROR_TABLE_BAD_KEY );
  assert( Util_TableCreateFromString("foo=12;" ) == UTIL_ERROR_BAD_INPUT );
  assert( Util_TableCreateFromString("foo=12,0" ) == UTIL_ERROR_BAD_INPUT );
  assert( Util_TableCreateFromString("foo='" ) == UTIL_ERROR_BAD_INPUT );
  assert( Util_TableCreateFromString("foo=\"" ) == UTIL_ERROR_BAD_INPUT );
  assert( Util_TableCreateFromString("foo=\"'" ) == UTIL_ERROR_BAD_INPUT );
  assert( Util_TableCreateFromString("foo={" ) == UTIL_ERROR_BAD_INPUT );
  assert( Util_TableCreateFromString("foo=}" ) == UTIL_ERROR_BAD_INPUT );
  assert( Util_TableCreateFromString("foo={0" ) == UTIL_ERROR_BAD_INPUT );
  assert( Util_TableCreateFromString("foo=0}" ) == UTIL_ERROR_BAD_INPUT );
  assert( Util_TableCreateFromString("foo={bar}" ) == UTIL_ERROR_BAD_INPUT );
  assert( Util_TableCreateFromString(" foo = { \"\r\t\n\v\" } " ) ==
          UTIL_ERROR_BAD_INPUT );
  assert( Util_TableCreateFromString("foo={0 1.0}" ) ==
          UTIL_ERROR_TABLE_NO_MIXED_TYPE_ARRAY );

  /*
   * Test some "good" strings with single values
   */
    {
  CCTK_INT  int_value = 0;
  CCTK_REAL real_value = 0.0;
  CCTK_CHAR string_value[10] = "";
  CCTK_INT  int_array[4] = {0, 0, 0, 0};
  CCTK_REAL real_array[4] = {0.0, 0.0, 0.0, 0.0};
  CCTK_INT  type, nelems;

  assert( Util_TableSetFromString(handle, "foo=12" ) == 1 );
  assert( Util_TableSetFromString(handle, " foo=12" ) == 1 );
  assert( Util_TableSetFromString(handle, "foo =12" ) == 1 );
  assert( Util_TableSetFromString(handle, "foo= 12" ) == 1 );
  assert( Util_TableSetFromString(handle, "foo=12 " ) == 1 );
  assert( Util_TableSetFromString(handle, " foo = 12 " ) == 1 );
  assert( Util_TableSetFromString(handle, " foo = .0 " ) == 1 );
  assert( Util_TableSetFromString(handle, " foo = 12. " ) == 1 );

  assert( Util_TableSetFromString(handle, " foo = +12 " ) == 1 );
  assert( Util_TableGetInt (handle, &int_value, "foo") == 1 );
  assert( int_value == 12 );

  assert( Util_TableSetFromString(handle, " foo = 012 " ) == 1 );
  assert( Util_TableGetInt (handle, &int_value, "foo") == 1 );
  assert( int_value == 12 );

  assert( Util_TableSetFromString(handle, " foo = -12.0 " ) == 1 );
  assert( Util_TableGetReal (handle, &real_value, "foo") == 1 );
  assert( real_value == -12.0 );

  assert( Util_TableSetFromString(handle, " foo = '\\nbar\\r' " ) == 1 );
  assert( Util_TableGetCharArray (handle, sizeof (string_value), string_value,
                                  "foo") == (int) strlen ("\\nbar\\r") );
  string_value[strlen ("\\nbar\\r")] = 0;
  assert( strcmp ((char *) string_value, "\\nbar\\r") == 0 );

  assert( Util_TableSetFromString(handle, " foo = \"\tbar\v\" " ) == 1 );
  assert( Util_TableGetCharArray (handle, sizeof (string_value), string_value,
                                  "foo") == (int) strlen ("\tbar\v") );
  string_value[strlen ("\tbar\v")] = 0;
  assert( strcmp ((char *) string_value, "\tbar\v") == 0 );

  assert( Util_TableSetFromString(handle, " foo = {} " ) == 1 );
  assert( Util_TableQueryValueInfo (handle, &type, &nelems, "foo") == 1 );
  assert( type == CCTK_VARIABLE_INT && nelems == 0 );

  assert( Util_TableSetFromString(handle, " foo = {\t-1 +2 -3 +4\t} " ) == 1 );
  assert( Util_TableGetIntArray (handle, 5, int_array, "foo") == 4 );
  assert( int_array[0] == -1 && int_array[1] == 2 &&
          int_array[2] == -3 && int_array[3] == 4 );

  assert( Util_TableSetFromString(handle,
                                  " foo = {\t-1.1\t+2.2\t-3.3\t+4.4\t} ") == 1);
  assert( Util_TableGetRealArray (handle, 100, real_array, "foo") == 4 );
  assert( real_array[0] == -1.1 && real_array[1] == 2.2 &&
          real_array[2] == -3.3 && real_array[3] == 4.4 );

  assert( Util_TableDeleteKey (handle, "foo") == 0 );
  assert( Util_TableQueryNKeys (handle) == 0 );
    }

  /*
   * Test some "good" strings with multiple values of different types
   */
    {
  CCTK_CHAR string_value[10] = "";
  CCTK_INT  int_array[3] = {0, 0, 0};
  CCTK_REAL real_array[3] = {0.0, 0.0, 0.0};

  const int handle2 = Util_TableCreate(UTIL_TABLE_FLAGS_DEFAULT);
  assert( Util_TableSetFromString(handle2, "ij = 42\t"
                                           "real1 = 3.5\r"
                                           "string = 'string'\v"
                                           "int_array = {0 1 2}\n"
                                           "real_array = {0.0 1.1 2.2}") == 5);
  assert( Util_TableQueryNKeys(handle2) == 5 );
  assert( Util_TableSetFromString(handle2, "real_e=2.75") == 1);
  assert( Util_TableQueryNKeys(handle2) == 6 );
  assert( Util_TableGetIntArray (handle2, 5, int_array, "int_array") == 3 );
  assert( Util_TableGetRealArray (handle2, 100, real_array, "real_array") == 3);
  assert( Util_TableGetCharArray (handle2, sizeof (string_value), string_value,
                                  "string") == (int) strlen ("string") );
  string_value[strlen ("string")] = 0;

  assert( int_array[0] == 0 && int_array[1] == 1 && int_array[2] == 2 );
  assert( real_array[0] == 0.0 && real_array[1] == 1.1 && real_array[2] == 2.2);
  assert( strcmp ((char *) string_value, "string") == 0 );

  assert( Util_TableDeleteKey (handle2, "real_array") == 0 );
  assert( Util_TableDeleteKey (handle2, "int_array") == 0 );
  assert( Util_TableDeleteKey (handle2, "string") == 0 );
  check_table_contents(handle2, true);

  assert( Util_TableDestroy(handle2) == 0 );
    }
    {
  const int handle3
          = Util_TableCreateFromString("  ij=42 real1=3.5 real_e=2.75  ");
  assert( handle3 >= 0 );

  assert( Util_TableQueryFlags(handle3) == UTIL_TABLE_FLAGS_CASE_INSENSITIVE );
  assert( Util_TableQueryNKeys(handle3) == 3 );
  check_table_contents(handle3, true);

  return handle3;
    }
    }
}

/******************************************************************************/

/*
 * This function tests  Util_Table{Set,Get}String()
 */
static
  void test_set_get_string(int handle, bool case_insensitive)
{
  printf("test_set_get_string(handle=%d, case_insensitive=%d)\n",
         handle, (int)case_insensitive);
  fflush(stdout);

  assert( Util_TableSetString(handle, "Germany", "AEI") == 0 );
  assert( Util_TableSetString(handle, "Golm", "AEI") == 1 );

    {
  CCTK_INT type_code, N_elements;
  assert( Util_TableQueryValueInfo(handle,
                                   &type_code, &N_elements,
                                   case_insensitive ? "aei" : "AEI") == 1 );
  assert( type_code == CCTK_VARIABLE_CHAR );
  assert( N_elements == (int)strlen("Golm") );

    {
  const int N_buffer = N_elements+1;
  char *const buffer = (char *) malloc(N_buffer);
  assert( buffer != NULL );
  assert( Util_TableGetCharArray(handle,
                                 N_buffer, (CCTK_CHAR *) buffer,
                                 "AEI") == (int)strlen("Golm") );
  assert( Util_TableGetString(handle,
                              0, NULL,
                              "AEI") == (int)strlen("Golm") );
  assert( Util_TableGetString(handle,
                              N_buffer, buffer,
                              case_insensitive ? "aEI" : "AEI")
          == (int)strlen("Golm") );
  assert( strcmp(buffer, "Golm") == 0 );

  /* check getting string longer than buffer */
  assert( Util_TableSetString(handle, "Max-Planck", "famous") == 0 );
  assert( Util_TableGetString(handle,
                              N_buffer, buffer,
                              case_insensitive ? "FAMouS" : "famous")
          == UTIL_ERROR_TABLE_STRING_TRUNCATED );
  assert( strcmp(buffer, "Max-") == 0 );
    }
    }
}

/******************************************************************************/

/*
 * This function tests
 *      Util_Table{Set,Get}Pointer()
 *      Util_Table{Set,Get}PointerArray()
 */
static
  void test_set_get_pointers(int handle)
{
  printf("test_set_get_pointers(handle=%d)\n", handle);
  fflush(stdout);

    {
  CCTK_INT i, j;
  assert( Util_TableSetPointer(handle, (CCTK_POINTER) &i, "i_ptr") == 0);
  assert( Util_TableSetPointer(handle, (CCTK_POINTER) &j, "j_ptr") == 0);

    {
  CCTK_POINTER iptr, jptr;
  assert( Util_TableGetPointer(handle, &iptr, "i_ptr") == 1 );
  assert( (CCTK_INT*) iptr == &i );
  assert( Util_TableGetPointer(handle, &jptr, "j_ptr") == 1 );
  assert( (CCTK_INT*) jptr == &j );

    {
  CCTK_POINTER ijptr[2];
  ijptr[0] = &j;
  ijptr[1] = &i;
  assert( Util_TableSetPointerArray(handle, 2, ijptr, "ijptr") == 0 );

    {
  CCTK_POINTER ijptr_copy[5];
  assert( Util_TableGetPointerArray(handle, 5, ijptr_copy, "ijptr") == 2 );
  assert( ijptr_copy[0] = ijptr[0] );
  assert( ijptr_copy[1] = ijptr[1] );
    }
    }
    }
    }
}

/******************************************************************************/

/*
 * This function tests
 *      Util_Table{Set,Get}FPointer()
 *      Util_Table{Set,Get}FPointerArray()
 */
static
  void test_set_get_fpointers(int handle)
{
  printf("test_set_get_fpointers(handle=%d)\n", handle);
  fflush(stdout);

    {
  CCTK_FPOINTER dptr = (CCTK_FPOINTER) & test_set_get_pointers;
  CCTK_FPOINTER fptr = (CCTK_FPOINTER) & test_set_get_fpointers;
  assert( Util_TableSetFPointer(handle, dptr, "dptr") == 0);
  assert( Util_TableSetFPointer(handle, fptr, "fptr") == 0);

    {
  CCTK_FPOINTER dptr_copy, fptr_copy;
  assert( Util_TableGetFPointer(handle, &dptr_copy, "dptr") == 1 );
  assert( dptr_copy == dptr );
  assert( Util_TableGetFPointer(handle, &fptr_copy, "fptr") == 1 );
  assert( fptr_copy == fptr );

    {
  CCTK_FPOINTER fdptrs[2];
  fdptrs[0] = fptr;
  fdptrs[1] = dptr;
  assert( Util_TableSetFPointerArray(handle, 2, fdptrs, "fdptrs") == 0 );
  
    {
  CCTK_FPOINTER fdptrs_copy[5];
  assert( Util_TableGetFPointerArray(handle, 2, fdptrs_copy, "fdptrs") == 2 );
  assert( fdptrs_copy[0] = fdptrs[0] );
  assert( fdptrs_copy[1] = fdptrs[1] );
    }
    }
    }
    }
}

/******************************************************************************/

/*
 * This function tests cloning a table.  We assume that on entry the
 * table contains the 3 keys
 *      real_e = 2.75
 *      real1 = 3.5
 *      ij = 42
 * in that order.
 */
void test_clone(int handle)
{
  printf("test_clone(handle=%d)\n", handle);
  fflush(stdout);

  check_table_contents(handle, true);     /* entry assumption */
    {
  const int clone_handle = Util_TableClone(handle);

  /* make sure we didn't modify the table that we cloned */
  check_table_contents(handle, true);

  /* check the clone */
  check_table_contents(clone_handle, false);

  /*
   * check that the tables are now distinct by
   * changing each of them and checking that the other is unchanged
   */
  assert( Util_TableSetInt(handle, 105, "universal number") == 0 );
  check_table_contents(clone_handle, false);
  assert( Util_TableDeleteKey(handle, "universal number") == 0 );

  assert( Util_TableSetInt(clone_handle, 105, "universal number #2") == 0 );
  check_table_contents(handle, true);
  assert( Util_TableDeleteKey(clone_handle, "universal number #2") == 0 );

  assert( Util_TableDestroy(clone_handle) == 0 );
    }
}

/******************************************************************************/

/*
 * This function does a sequence of assert() calls to verify that
 * a table contains the 3 keys
 *      real_e = 2.75
 *      real1 = 3.5
 *      ij = 42
 * in a specified order.
 *
 * Arguments:
 * order_up_flag = true --> check for the order {real_e,real1,ij}
 *                 false --> check for the order {ij,real1,real_e}
 *
 * Bugs:
 * Having to test for a specific order is a kludge.
 */
static
  void check_table_contents(int handle, bool order_up_flag)
{
  printf("check_table_contents(handle=%d, order_up_flag=%d)\n",
         handle, (int) order_up_flag);
  fflush(stdout);

  assert( Util_TableQueryNKeys(handle) == 3 );

  /* set up the key buffer */
    {
  const int max_key_length = Util_TableQueryMaxKeyLength(handle);
  assert( max_key_length == (int)strlen("real_e") );
    {
  const int N_key_buffer = max_key_length + 1;
  char *const key_buffer = malloc(N_key_buffer);
  assert( key_buffer != NULL );

  /* walk through the table to verify contents in the right order */
    {
  const int ihandle = Util_TableItCreate(handle);
  if (order_up_flag)
  {
    check_table_contents_real_e(handle,ihandle);
    assert( Util_TableItAdvance(ihandle) == 1 );
    check_table_contents_real1(handle,ihandle);
    assert( Util_TableItAdvance(ihandle) == 1 );
    check_table_contents_ij(handle, ihandle);
  }
  else
  {
    check_table_contents_ij(handle, ihandle);
    assert( Util_TableItAdvance(ihandle) == 1 );
    check_table_contents_real1(handle,ihandle);
    assert( Util_TableItAdvance(ihandle) == 1 );
    check_table_contents_real_e(handle,ihandle);
  }
  assert( Util_TableItAdvance(ihandle) == 0 );
  assert(Util_TableItDestroy(ihandle) == 0 );
    }
    }
    }
}

/******************************************************************************/

/*
 * This function does a sequence of assert() calls to verify that
 * the longest key in a table has the length of "real_e", and that
 * a table iterator points to the key  ij = 42 .
 */
static
  void check_table_contents_ij(int handle, int ihandle)
{
  printf("check_table_contents_ij(handle=%d, ihandle=%d)\n", handle, ihandle);
  fflush(stdout);

  /* set up the key buffer */
    {
  const int max_key_length = Util_TableQueryMaxKeyLength(handle);
  assert( max_key_length == (int)strlen("real_e") );
    {
  const int N_key_buffer = max_key_length + 1;
  char *const key_buffer = malloc(N_key_buffer);
  assert( key_buffer != NULL );

  /* check ij = 42 */
    {
  CCTK_INT type_code = 123456;
  CCTK_INT N_elements = 54321;
  assert( Util_TableItQueryKeyValueInfo(ihandle,
                                        N_key_buffer, key_buffer,
                                        &type_code, &N_elements)
          == (int)strlen("ij") );
  assert( strcmp(key_buffer, "ij") == 0 );
  assert( type_code == CCTK_VARIABLE_INT );
  assert( N_elements == 1 );
    {
  CCTK_INT value_int;
  assert( Util_TableGetInt(handle, &value_int, "ij") == 1 );
  assert( value_int == 42 );
    }
    }
    }
    }
}

/******************************************************************************/

/*
 * This function does a sequence of assert() calls to verify that
 * the longest key in a table has the length of "real_e", and that
 * a table iterator points to the key  real1 = 3.5 .  (N.b. 3.5 is
 * assumed to be exactly representable as a floating point number.
 * This should be ok for any reasonable floating-point format.)
 */
static
  void check_table_contents_real1(int handle, int ihandle)
{
  printf("check_table_contents_real1(handle=%d, ihandle=%d)\n",
         handle, ihandle);
  fflush(stdout);

  /* set up the key buffer */
    {
  const int max_key_length = Util_TableQueryMaxKeyLength(handle);
  assert( max_key_length == (int)strlen("real_e") );
    {
  const int N_key_buffer = max_key_length + 1;
  char *const key_buffer = malloc(N_key_buffer);
  assert( key_buffer != NULL );

  /* check real1 = 3.5 */
    {
  CCTK_INT type_code = 123456;
  CCTK_INT N_elements = 54321;
  assert( Util_TableItQueryKeyValueInfo(ihandle,
                                        N_key_buffer, key_buffer,
                                        &type_code, &N_elements)
          == (int)strlen("real1") );
  assert( strcmp(key_buffer, "real1") == 0 );
  assert( type_code == CCTK_VARIABLE_REAL );
  assert( N_elements == 1 );
    {
  CCTK_REAL value_real;
  assert( Util_TableGetReal(handle, &value_real, key_buffer) == 1 );
  assert( value_real == 3.5 );
    }
    }
    }
    }
}

/******************************************************************************/

/*
 * This function does a sequence of assert() calls to verify that
 * the longest key in a table has the length of "real_e", and that
 * a table iterator points to the key  real_e = 2.75 .  (N.b. 2.75 is
 * assumed to be exactly representable as a floating point number.
 * This should be ok for any reasonable floating-point format.)
 */
static
  void check_table_contents_real_e(int handle, int ihandle)
{
  printf("check_table_contents_real_e(handle=%d, ihandle=%d)\n",
         handle, ihandle);
  fflush(stdout);

  /* set up the key buffer */
    {
  const int max_key_length = Util_TableQueryMaxKeyLength(handle);
  assert( max_key_length == (int)strlen("real_e") );
    {
  const int N_key_buffer = max_key_length + 1;
  char *const key_buffer = malloc(N_key_buffer);
  assert( key_buffer != NULL );

  /* check real_e = 2.75 */
    {
  CCTK_INT type_code = 123456;
  CCTK_INT N_elements = 54321;
  assert( Util_TableItQueryKeyValueInfo(ihandle,
                                        N_key_buffer, key_buffer,
                                        &type_code, &N_elements)
          == (int)strlen("real_e") );
  assert( strcmp(key_buffer, "real_e") == 0 );
  assert( type_code == CCTK_VARIABLE_REAL );
  assert( N_elements == 1 );
    {
  CCTK_REAL value_real;
  assert( Util_TableGetReal(handle, &value_real, key_buffer) == 1 );
  assert( value_real == 2.75 );

  value_real = 314159.271828;
  assert( Util_TableGetGeneric(handle,
                               CCTK_VARIABLE_REAL, (void *) &value_real,
                               key_buffer)
          == 1 );
  assert( value_real == 2.75 );

  value_real = 314159.271828;
  assert( Util_TableGetGenericArray(handle,
                                    CCTK_VARIABLE_REAL, 1, (void *) &value_real,
                                    key_buffer)
          == 1 );
  assert( value_real == 2.75 );
    }
    }
    }
    }
}
