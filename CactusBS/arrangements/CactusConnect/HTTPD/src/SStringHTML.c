 /*@@
  @file      SStringHTML.c
  @date      02.04.2004
  @author    Steve White
  @desc      Extension to Strings module with function specific to HTML
  @enddesc
  @version $Header$
  @@*/
#include "SStringHTML.h"

static String *
StringReplaceCharWithCString( String * str, SSCHAR c, const SSCHAR *cstr );

String * 
StringSetToEncodedHTMLCString( String * str, const SSCHAR *c )
{
        return StringEncodeHTML( StringSetToCString( str, c ) );
}

String * 
StringEncodeHTML( String * str )
{
        StringReplaceCharWithCString( str, '&', "&amp;" );
        StringReplaceCharWithCString( str, '<', "&lt;" );
        StringReplaceCharWithCString( str, '>', "&gt;" );
        return str;
}

String *
StringReplaceCharWithCString( String * str, SSCHAR c, const SSCHAR *cstr )
{
        size_t position = 0;
        while( StringFindCharFrom( str, c, &position ) )
        {
                StringDeleteChar( str, position );
                StringInsertCString( str, cstr, position );
                position ++;
        }
        return str;
}
