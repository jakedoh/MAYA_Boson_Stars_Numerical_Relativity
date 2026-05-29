 /*@@
  @file      SStringIO_Namespace.h
  @date      02.04.2004
  @author    Steve White
  @desc      Extensions to Strings module involvint file IO
  @enddesc
  @version $Header$
  @@*/
#ifndef _SSTRINGIO_NAMESPACE_H
#define _SSTRINGIO_NAMESPACE_H

#include "SStringIO.h"

#define ReadToDelimiter( a, f, d )	StringReadToDelimiter( a, f, d )
#define ReadToEndOfLine( a, f )		StringReadToEndOfLine( a, f )
#define ReadLine( a, f )		StringReadLine( a, f )
#define PrintToFile( a, f )		StringPrintToFile( a, f )
#define FormatPrintToFile( a, s, f )	StringFormatPrintToFile( a, s, f )

#endif
