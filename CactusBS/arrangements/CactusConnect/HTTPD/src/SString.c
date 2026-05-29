/*@@
  @file      SString.c
  @date      02.04.2004
  @author    Steve White
  @desc      Module for generic operations on strings
  @enddesc
  @version $Header$
  @@*/

#include "SString.h"
#include "SStringIO.h"

#include <string.h>
#include <stdlib.h>

#include "util_String.h"

#ifndef MIN
#define MIN( a, b ) ( (a) < (b) ? (a) : (b) )
#endif
#ifndef MAX
#define MAX( a, b ) ( (a) > (b) ? (a) : (b) )
#endif

typedef enum { MAKE_NULL_ALLOC, NEW_NULL_ALLOC, BUFFER_NULL_ALLOC }
SSTRING_ERROR;
static void String_HandleSeriousError( SSTRING_ERROR );

#define LINE_END_BUFSIZE 3

typedef struct String_tag
{
  size_t  length;
  SSCHAR  *chars;
  SSCHAR  line_end[LINE_END_BUFSIZE];
} String_placeholder;

static SSCHAR   kLineEndChar = '\n';

static SSBOOL   BufferForThisManyChars( String *s, size_t size );

size_t
StringLength( const String * s )
{
  return s->length;
}

const SSCHAR *
StringGetBuffer( const String * s )
{
  return s->chars;
}

void
StringCopyBuffer( const String * s, SSCHAR * buf, size_t bufSize )
{
  strncpy( buf, s->chars, bufSize );
  buf[bufSize-1] = '\0';
}

/*
 * Returns 0 if index is not in the string
 */
SSCHAR
StringNthChar( const String * s, size_t n )
{
  if( n < s->length )
  {
    return s->chars[n];
  }
  else
  {
    return '\0';
  }
}

String *
String_New()
{
  String *s = (String *)calloc( 1, sizeof( String ) );
  if( s != NULL )
  {
    BufferForThisManyChars( s, 0 );
  }
  else
  {
    String_HandleSeriousError( NEW_NULL_ALLOC );
  }             
  return s;
}

String *
String_Make( const SSCHAR * c_string )
{
  String *s = (String *)calloc( 1, sizeof( String ) );
  if( s != NULL )
  {
    if( c_string != NULL )
    {
      size_t len = strlen( c_string );
      s->chars = strdup( c_string );
      s->length = len;
    }
    else
    {
      BufferForThisManyChars( s, 0 );
    }
    s->line_end[0] = kLineEndChar;
  }
  else
  {
    String_HandleSeriousError( MAKE_NULL_ALLOC );
  }
  return s;
}

/* note this is a little ambiguous.
 * needs to be made clear this makes a new string
 */
String *
String_Copy( const String * other )
{
  String *s = StringSetToCString( String_New(),
                                  StringGetBuffer( other ) );
  memcpy( s->line_end, other->line_end, LINE_END_BUFSIZE
          * sizeof( SSCHAR ) );
  return s;
}

void
String_Delete( String * s )
{
  if( s != NULL )
  {
    free( s->chars );
  }
  free( s );
}

String *
StringTruncate( String * s, size_t n )
{
  if( n < s->length )
  {
    s->chars[n] = '\0';
    s->length = n;
  }
  return s;
}

String *
StringTrimLeading( String * s, size_t n )
{
  const size_t orig_len = StringLength( s );

  if( orig_len > 0 && n > 0 )
  {
    const size_t position = MIN( orig_len, n ) - 1;
    const size_t new_len = orig_len - position;
                
    memmove( s->chars, s->chars + position, new_len );
    s->chars[new_len] = '\0';
    s->length = new_len;
  }
  return s;
}

String *
StringSetNthChar( String * s, size_t n, SSCHAR c )
{
  if( n < s->length )
  {
    if( c == '\0' )
    {
      StringTruncate( s, n );
    }
    else
    {
      s->chars[n] = c;
    }
  }
  return s;
}
/* Allocates if there is no buffer
 * Re-allocates (without changing remaining string) if there is a buffer
 * Puts a null char at the location of the corresponding size
 * DOES NOT pad with zeros, or initialize the buffer beyond this
 */
SSBOOL
BufferForThisManyChars( String *s, size_t size )
{
  const size_t    blocksize = sizeof( SSCHAR ) * ( size + 1 );

  if( s->chars == NULL )
  {
    s->chars = (SSCHAR *)malloc( blocksize );
  }
  else if( s->length < size )
  {
    s->chars = (SSCHAR *)realloc( s->chars, blocksize );
  }
  else if( s->length > ( size << 1 ) ) /* Don't resize if larger
                                          but in ballpark */
  {
    s->length = 0;
    s->chars = (SSCHAR *)malloc( blocksize );
  }

  if( s->chars != NULL )
  {
    s->length = size;
    if( s->chars != NULL )
      s->chars[size] = '\0';
  }
  else
  {
    String_HandleSeriousError( BUFFER_NULL_ALLOC );
  }

  return s->chars != NULL;
}

String *
StringSetToCString( String * s, const SSCHAR *c_string )
{
  return StringSetToBuffer( s, c_string, strlen( c_string ) );
}

size_t
StringSetNextToken( const String *s, String *token,
                    const SSCHAR *delim, size_t start )
{
  if( start < s->length )
  {
    const size_t    ndelims = strlen( delim );
    size_t          next = s->length;
    size_t          d;

    for( d = 0; d < ndelims; d++ )
    {
      size_t nextDelim = start;
      if( StringFindCharFrom( s, delim[d], &nextDelim ) )
        next = MIN( nextDelim, next ); 
    }

    StringSetRange( token, s, start, next - start );

    return next < s->length ? next + 1 : next;
  }
  else
  {
    return 0;
  }
}

String *
StringConcatCString( String * s, const SSCHAR *c_string )
{
  const size_t orig_len = s->length;
  if( BufferForThisManyChars( s, orig_len + strlen( c_string ) ) )
  {
    strcpy( s->chars + orig_len, c_string );
  }
  return s;
}

String *
StringSetToBuffer( String * s, const SSCHAR *buf, size_t len )
{
  if( s != NULL && buf != NULL )
    if( BufferForThisManyChars( s, len ) )
    {
      strncpy( s->chars, buf, len );
      s->chars[ len ] = '\0';
    }
  return s;
}

String *
StringSet( String * s, const String * other )
{
  return StringSetToCString( s, StringGetBuffer( other ) );
}
/*
 * Take care! What if this String's buffer is the same as the other's?
 * For now, this can happen only if the other String is the same as this
 * String.
 *
 * This needs thought for every String that takes a String argument,
 * and for the ones that take C strings as well.
 */
String *
StringSetRange( String * s, const String * other,
                size_t first, size_t length )
{
  const size_t    other_length = other->length;
  const size_t    veryFirst = MIN( first + 1, other_length ) - 1,
      minLength = MIN( first + length, other_length ),
      newLength = minLength - veryFirst;

  return StringSetToBuffer( s, other->chars + veryFirst, newLength );
}

String *
StringSetToPartAfter( String * s, const String * other, size_t position )
{
  return StringSetRange( s, other,
                         position + 1, StringLength( other ) - position );
}

/*
 * Last argument 'position' is both input and output.
 * Specifies where to begin search: pass 0 to search whole string
 * On output, is position where char was found.
 */
SSBOOL
StringFindCharFrom( const String *s, SSCHAR theChar, size_t * position )
{
  SSCHAR  * charPtr;

  if( s->length > 0
      && position
      && *position < s->length
      && ( charPtr = (SSCHAR *)strchr( s->chars + *position, theChar ) )
      != NULL )
  {
    *position = ( charPtr - s->chars ) / sizeof( SSCHAR );
    return SSTRUE;
  }
  return SSFALSE;
}

SSBOOL
StringFindStringFrom( const String * s, const String * other,
                      size_t * position )
{
  SSCHAR  * charPtr;

  if( s->length >= other->length
      && position
      && *position < s->length
      && s->length > 0
      && other->length > 0
      && ( charPtr = (SSCHAR *)strstr( s->chars + *position, other->chars ) )
      != NULL )
  {
    *position = ( charPtr - s->chars ) / sizeof( SSCHAR );
    return SSTRUE;
  }
  return SSFALSE;
}

SSBOOL
StringEquals( const String * a, const String * b )
{
  if( a->length > 0 )
  {
    return strncmp( a->chars, b->chars, a->length ) == 0;
  }
  else if( b->length == 0 )
  {
    return SSTRUE;
  }
  return SSFALSE;
}

int
StringCompare( const String * a, const String * b )
{
  return StringCompareCString( a, b->chars );
}

int
StringCompareCString( const String * a, const char * b )
{
  const size_t    length = MIN( a->length, strlen( b ) );

  if( length > 0 )
  {
    return strncmp( a->chars, b, length );
  }
  else if( strlen( b ) == 0 )
  {
    return 0;
  }
  else            /* b is empty but a isn't */
  {
    return 1;
  }
}

String *
StringInsert( String * s, const String * other, size_t position )
{
  const size_t    oldLength = s->length, otherLength = other->length;
  String          *old = String_Copy( s );

  if( otherLength > 0
      && position <= oldLength
      && BufferForThisManyChars( s, oldLength + otherLength ) )
  {
    if( position != 0 )
    {
      strncpy( s->chars, old->chars, position );
    }

    strncpy( s->chars + position, other->chars, otherLength );

    if( position < oldLength )
    {
      strncpy( s->chars + position + otherLength,
               old->chars + position, oldLength - position );
    }
    s->chars[s->length] = '\0';
  }
  String_Delete( old );

  return s;
}

String *
StringInsertCString( String * s, const char * c_string, size_t position )
{
  String *other = String_Make( c_string );
  StringInsert( s, other, position );
  String_Delete( other );
  return s;
}
/*
 * What to do if position is off end of string?  We do nothing
 */
String *
StringInsertChar( String * s, SSCHAR c, size_t position )
{
  if( position <= s->length )
  {
    BufferForThisManyChars( s, s->length + 1 );
    if( position + 2 < s->length )
    {
      s->chars[s->length] = '\0';
      memmove( s->chars + position + 1, s->chars + position,
               s->length - 1 - position );
      s->chars[position] = c;
    }
  }
  return s;
}

String *
StringDeleteChar( String * s, size_t position )
{
  return StringDeleteRange( s, position, position );
}

String *
StringDeleteRange( String * s, size_t begin, size_t end )
{
  if( begin <= end && end < s->length )
  {
    if( end + 1 < s->length )
    {
      memmove( s->chars + begin, s->chars + end + 1,
               ( s->length - end ) - 1 );
    }
    s->length -= end - begin + 1;
    s->chars[s->length] = '\0';
  }
  return s;
}

String *
StringConcat( String * s, const String * other )
{
  return StringInsert( s, other, s->length );
}

/*
 * On unix the default \n works; on the Mac, you might want \r
 */
void
StringSetLineEndCharacter( String * s, const SSCHAR *end )
{
  strncpy( s->line_end, end, LINE_END_BUFSIZE );
  s->line_end[LINE_END_BUFSIZE - 1] = '\0';
}

String *
StringReadToEndOfLine( String * s, FILE * file )
{
  return StringReadToDelimiter( s, file, kLineEndChar );
}

String *
StringReadToDelimiter( String * s, FILE * file, SSCHAR delim )
{
  int     next;

  while( ( next = fgetc( file ) ) != EOF 
         && (SSCHAR)next != delim )
  {
    StringInsertChar( s, (SSCHAR)next, s->length );
  }

  return s;
}

void
StringPrint( const String * s )
{
  fprintf( stdout, "%s", s->chars );
}

void
StringPrintToFile( const String * s, FILE * file )
{
  fprintf( file, "%s", s->chars );
}

void
StringFormatPrint( const String * s, const SSCHAR *format )
{
  fprintf( stdout, format, s->chars );
}

void
StringFormatPrintToFile( const String * s, const SSCHAR *format, FILE * file )
{
  fprintf( file, format, s->chars );
}

#define DECIMALBUFSIZE 64
String *
StringConcatDecimal( String * s, long d )
{
  char buf[DECIMALBUFSIZE] = { '\0' };
  Util_snprintf( buf, sizeof( buf ), "%ld", d );
  return StringConcatCString( s, buf );
}

String *
StringConcatDouble( String * s, double d )
{
  char buf[DECIMALBUFSIZE] = { '\0' };
  Util_snprintf( buf, sizeof( buf ), "%f", d );
  return StringConcatCString( s, buf );
}
/*
 *
 * 0) It's inappropritate to handle string and character convesions.
 *    And there are some very specific-use conversions, such as 't'
 *
 * 1) I'm not sure I completely understand all the conversions.  
 *    specifically, I don't get G and g.
 *
 * 2) There are several standards, including C99 and SUSv2
 *
 * 3) I haven't done long long or long double
 *
 * 4) Mixed up notions of unsigned with hex and octal...wrong?
 *
 * 
 * 
 */
static void
addNumericMods( char * format, SSFORMAT f )
{
  strcat( format, "%" );
  if( f & SFMT_LEFT_ALIGN )
  {
    strcat( format, "-" );
  }
  if( f & SFMT_PRINT_SIGN )
  {
    strcat( format, "+" );
  }
  if( f & SFMT_ADD_SIGN_SPACE )
  {
    strcat( format, " " );
  }
  if( f & SFMT_PAD_ZERO )
  {
    strcat( format, "0" );
  }
  if( f & SFMT_ALT )
  {
    strcat( format, "#" );
  }
}
#define EMPTYSTRING { '\0' }
String *
StringConcatFormattedDecimal( String *s, long int d,
                              int width, int precision, SSINT_FORMAT f )
{
  char buf[DECIMALBUFSIZE] = EMPTYSTRING;
  char format[16] = EMPTYSTRING;
  addNumericMods( format, f );
  strcat( format, "*.*ld" );
  Util_snprintf( buf, sizeof( buf ), format, width, precision, d );
  return StringConcatCString( s, buf );
}

String *
StringConcatFormattedUnsigned( String *s, unsigned long int n,
                               int width, int precision, SSINT_FORMAT f )
{
  char buf[DECIMALBUFSIZE] = EMPTYSTRING;
  char format[16] = EMPTYSTRING;
  addNumericMods( format, f );
  strcat( format, "*.*l" );
  if( f & SFMT_HEX )
  {
    if( f & SFMT_LOWERCASE )
    {
      strcat( format, "x" );
    }
    else
    {
      strcat( format, "X" );
    }
  }
  else if( f & SFMT_OCTAL )
  {
    strcat( format, "o" );
  }
  else
  {
    strcat( format, "u" );
  }
  Util_snprintf( buf, sizeof( buf ), format, width, precision, n );
  return StringConcatCString( s, buf );
}

String *
StringConcatFormattedDouble( String *s, double d,
                             int width, int precision, SSDOUBLE_FORMAT f  )
{
  char buf[DECIMALBUFSIZE] = EMPTYSTRING;
  char format[16] = EMPTYSTRING;
  addNumericMods( format, f );
  if( f & SFMT_ALWAYS_POINT )
  {
    strcat( format, "#" );
  }
  strcat( format, "*.*" );
  if( f & SFMT_EXPONENTIAL )
  {
    if( f & SFMT_LOWERCASE )
    {
      strcat( format, "e" );
    }
    else
    {
      strcat( format, "E" );
    }
  }
  if( f & SFMT_DISCRETIONARY_EXP )
  {
    if( f & SFMT_LOWERCASE )
    {
      strcat( format, "g" );
    }
    else
    {
      strcat( format, "G" );
    }
  }
  else
  {
    strcat( format, "f" );
  }
  Util_snprintf( buf, sizeof( buf ), format, width, precision, d );
  return StringConcatCString( s, buf );
}

void
String_HandleSeriousError( SSTRING_ERROR e )
{
  /* to be filled in on a per-implementation basis */
}

