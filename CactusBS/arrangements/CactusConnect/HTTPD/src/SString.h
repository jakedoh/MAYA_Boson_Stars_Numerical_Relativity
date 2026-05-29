 /*@@
  @file      SString.h
  @date      02.04.2004
  @author    Steve White
  @desc      C Module for generic operations on strings
  @enddesc
  @version $Header$
  @@*/
#ifndef _SSTRING_H
#define _SSTRING_H

#include <stddef.h>

typedef struct String_tag String;       /* The abstract type for a String */

typedef char SSCHAR;                    /* Could be Unicode... */

typedef enum { SSFALSE, SSTRUE } SSBOOL;

#ifdef __cplusplus
extern "C" 
{
#endif
                                /* String creation and deletion */
String          *String_New( void );
String          *String_Copy( const String *other );
String          *String_Make( const SSCHAR *other );
void            String_Delete( String * );
                                /* Number of characters in string */
size_t          StringLength( const String * );
                                /* Accessors */
SSCHAR          StringNthChar( const String *, size_t n );
String *        StringSetNthChar( String *, size_t n, SSCHAR c );
                                /* Copying */
String *        StringSet( String *dest, const String *source );
String *        StringSetRange( String *dest, const String *source,
                                size_t start, size_t length );
String *        StringSetToPartAfter( String *dest, const String *source, 
                                size_t position );
                                /* Conversion to and from C string */
String *        StringSetToCString( String *s, const SSCHAR *c_str );
String *        StringConcatCString( String *s, const SSCHAR *c_str );
int             StringCompareCString( const String *s, const SSCHAR *c_str );
String *        StringInsertCString( String * s, const char * c_str,
                                                size_t position );
String *        StringSetToBuffer( String *s, const SSCHAR *buf,
                                                size_t len );
const SSCHAR *  StringGetBuffer( const String * s );
void            StringCopyBuffer( const String * s, SSCHAR * buf,
                                                size_t bufsize );
                                /* Searching */
SSBOOL          StringFindStringFrom( const String *s, const String *toFind,
                                size_t *position );
SSBOOL          StringFindCharFrom( const String *s, SSCHAR theChar,
                                size_t *position );
size_t          StringSetNextToken( const String *s, String *token,
                                const SSCHAR *delim, size_t index );
                                /* Comparison */
SSBOOL          StringEquals( const String *a, const String *b );
int             StringCompare( const String *a, const String *b );
                                /* Insertion and Concatenation */
String *        StringInsert( String * s, const String * other, size_t pos );
String *        StringConcat( String * s, const String * other);
String *        StringInsertChar( String * s, SSCHAR c, size_t pos );
                                /* Deleting and Trimming */
String *        StringDeleteChar( String * s, size_t pos );
String *        StringDeleteRange( String * s, size_t begin, size_t end );
String *        StringTruncate( String *, size_t n );
String *        StringTrimLeading( String *, size_t n );
                                /* For line reading and writing */
void            StringSetLineEndChars( String *str, const SSCHAR * );
                                /* Printing to stdout */
void            StringPrint( const String *str );
void            StringFormatPrint( const String *str, const SSCHAR *format );
                                /* Numeric conversions */
double          StringAsDouble( const String *str );
long            StringAsLong( const String *str );
unsigned long   StringAsUnsignedLong( const String *str );

typedef enum SSFORMAT_TAG
{
        SFMT_DEFAULT            = 0,
        SFMT_LEFT_ALIGN         = 1,
        SFMT_PAD_ZERO           = 1<<1,
        SFMT_LOWERCASE          = 1<<2,
        SFMT_ADD_SIGN_SPACE     = 1<<3,
        SFMT_PRINT_SIGN         = 1<<4,
        SFMT_ALT                = 1<<5
} SSFORMAT;

typedef enum SSDOUBLE_FORMAT_TAG 
{
        SFMT_EXPONENTIAL        = 1<<6,
        SFMT_DISCRETIONARY_EXP  = 1<<7
} SSDOUBLE_FORMAT;

typedef enum SSINT_FORMAT_TAG
{
        SFMT_HEX                = 1<<6,
        SFMT_OCTAL              = 1<<7
} SSINT_FORMAT;
#define SFMT_ALWAYS_POINT       SFMT_ALT
#define SFMT_DONT_TRIM_ZEROS    SFMT_ALT

String *        StringConcatFormattedDecimal( String *s, long int d,
                        int width, int precision, SSINT_FORMAT f );
String *        StringConcatFormattedUnsigned( String *s, unsigned long int n,
                        int width, int precision, SSINT_FORMAT f );
String *        StringConcatFormattedDouble( String *s, double d,
                        int width, int precision, SSDOUBLE_FORMAT f  );
String *        StringConcatDecimal( String *str, long int d );
String *        StringConcatDouble( String *str, double d );

#ifdef __cplusplus
}
#endif

#endif
