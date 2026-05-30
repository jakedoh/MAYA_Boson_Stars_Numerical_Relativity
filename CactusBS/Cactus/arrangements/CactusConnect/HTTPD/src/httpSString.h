 /*@@
   @header    http_SString.h
   @date      April 7 14:19:23 2004
   @author    Steve White
   @desc 
   Routines exported by the Content stuff.
   @enddesc
   @version $Header$
 @@*/

#ifndef __HTTP_SSTRING_H__
#define __HTTP_SSTRING_H__ 1

#include "SString.h"

#define EMPTYSTRING {'\0'}

#ifdef __cplusplus
extern "C" 
{
#endif

int HTTP_SendString(httpRequest *request, const String* message);
SSBOOL HTTP_GetHeaderValueString(const httpRequest *request,
                const String *header, String *value);
#ifdef __cplusplus
}
#endif

#endif 
