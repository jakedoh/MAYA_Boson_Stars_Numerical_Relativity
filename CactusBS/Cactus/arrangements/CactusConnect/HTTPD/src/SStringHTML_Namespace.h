 /*@@
  @file      SStringHTMLNamespace.h
  @date      02.04.2004
  @author    Steve White
  @desc      Extension to Strings module with function specific to HTML
  @enddesc
  @version $Header$
  @@*/
#ifndef _SSTRINGHTML_NAMESPACE_H
#define _SSTRINGHTML_NAMESPACE_H

#include "SStringHTML.h"


#define EncodeHTML( s ) \
		StringEncodeHTML( s )
#define SetToEncodedHTMLCString( s, c ) \
		StringSetToEncodedHTMLCString( s, c )

#endif
