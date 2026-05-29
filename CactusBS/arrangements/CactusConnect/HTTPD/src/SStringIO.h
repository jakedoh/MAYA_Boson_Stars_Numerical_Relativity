 /*@@
  @file      SStringIO.h
  @date      02.04.2004
  @author    Steve White
  @desc      Extensions to Strings module involvint file IO
  @enddesc
  @version $Header$
  @@*/
#ifndef _SSTRINGIO_H
#define _SSTRINGIO_H

#include "SString.h"
#include <stdio.h>

#ifdef __cplusplus
extern "C"
{
#endif
				/* File utilities */
String *	StringReadToDelimiter( String *str, FILE * file, char delim );
String *	StringReadToEndOfLine( String *str, FILE * is );
String *	StringReadLine( String *str, FILE * is );

void		StringPrintToFile( const String *str, FILE * is );
void		StringFormatPrintToFile( const String *str, const char *format,
					FILE * is );

#ifdef __cplusplus
}
#endif

#endif
