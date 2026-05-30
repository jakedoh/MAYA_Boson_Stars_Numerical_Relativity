 /*@@
  @file      SStringHTML.h
  @date      02.04.2004
  @author    Steve White
  @desc      Extension to Strings module with function specific to HTML
  @enddesc
  @version $Header$
  @@*/
#ifndef _SSTRINGHTML_H
#define _SSTRINGHTML_H

#include "SString.h"

/*
const String * StringHTMLBodyTag = StringHTMLBeginTag( BODY, NULL );
const String * StringHTMLBodyEndTag = StringHTMLEndTag( BODY, NULL );
String *	StringSGMLEntity( SGMLEntityNum n );
*/

/* At least encode unsafe ASCII according to rfc1738.html
 * < > & # % " ' space tab { } | \ ^ ~ ] [ ` ? ; : , @ = 
 * as well as 80-FF 00-1F and 7F
 * should be changed to %nn where nn is the hex representation of
 * the ASCII value.
 *
 * If move to bigger characters will have to encode more.
 * */
#ifdef __cplusplus
extern "C"
{
#endif

String * StringSetToEncodedURLCString( String * str, const char *c );
String * StringEncodeURL( String * str );

/* Replace < > & with character entities
 * */
String * StringSetToEncodedHTMLCString( String * str, const char *c );
String * StringEncodeHTML( String * str );

#ifdef __cplusplus
}
#endif

#endif
