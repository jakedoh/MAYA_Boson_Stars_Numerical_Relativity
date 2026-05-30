 /*@@
  @file      SString_Namespace.h
  @date      02.04.2004
  @author    Steve White
  @desc      Module for generic operations on strings
  @enddesc
  @version $Header$
  @@*/
#ifndef _SSTRING_NAMESPACE_H
#define _SSTRING_NAMESPACE_H

#include "SString.h"

/* A poor man's namespace for the String module */

#define Set( a, b ) \
		StringSet( a, b )
#define SetSubString( a, b, p, l ) \
		StringSetSubString( a, b, p, l )
#define SetToCString( a, b ) \
		StringSetToCString( a, b )
#define SetToPartAfter( a, b, p ) \
		StringSetToPartAfter( a, b, p )
#define InsertCString( a, b, c ) \
		StringInsertCString( a, b, c )
#define ConcatCString( a, b ) \
		StringConcatCString( a, b )
#define CompareCString( a, b ) \
		StringCompareCString( a, b )
#define SetToBuffer( a, b, l ) \
		StringSetToBuffer( a, b, l )
#define GetBuffer( a ) \
		StringGetBuffer( a )
#define Length( p ) \
		StringLength( p )
#define	NthChar( s, n ) \
		StringNthChar( s, n )
#define SetNthChar( s, n, c ) \
		StringSetNthChar( s, n, c )
#define Truncate( s, n ) \
		StringTruncate( s, n )
#define TrimLeading( s, n ) \
		StringTrimLeading( s, n )
#define	FindStringFrom( s, c, p ) \
		StringFindStringFrom( s, c, p )
#define	FindCharFrom( s, c, p ) \
		StringFindCharFrom( s, c, p )
#define	SetNextToken( s, c, p, r ) \
		StringSetNextToken( s, c, p, r )
#define Compare( a, b ) \
		StringCompare( a, b )
#define Equals( a, b ) \
		StringEquals( a, b )
#define Insert( a, b, p ) \
		StringInsert( a, b, p )
#define InsertChar( a, b, p ) \
		StringInsertChar( a, b, p )
#define Concat( a, b ) \
		StringConcat( a, b )
#define Print( a ) \
		StringPrint( a )
#define FormatPrint( a, s ) \
		StringFormatPrint( a, s )
#define ConcatDecimal( a, s ) \
		StringConcatDecimal( a, s )
#define ConcatHex( a, s ) \
		StringConcatHex( a, s )
#define ConcatOctal( a, s ) \
		StringConcatOctal( a, s )
#define ConcatDouble( a, s ) \
		StringConcatDouble( a, s )
#define ConcatFormattedDecimal( a, d, s1, s2, f ) \
		StringConcatFormattedDecimal( a, d, s1, s2, f )
#define ConcatFormattedUnsigned( a, d, s1, s2, f ) \
		StringConcatFormattedUnsigned( a, d, s1, s2, f )
#define ConcatFormattedDouble( a, d, s1, s2, f ) \
		StringConcatFormattedDouble( a, d, s1, s2, f )

#endif
