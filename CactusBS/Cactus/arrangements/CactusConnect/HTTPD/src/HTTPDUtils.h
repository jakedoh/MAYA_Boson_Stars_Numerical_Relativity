 /*@@
   @header    HTTPDUtils.h
   @date      April 7 2004
   @author    Steve White
   @desc 
   Routines exported by HTTPD to other thorns.
   @enddesc
   @version $Header$
 @@*/

#ifndef __HTTP_UTILS_H__
#define __HTTP_UTILS_H__ 1

#define HTTP_QUICKLINK 1

#include "http_Request.h"
#include "http_SString.h"
#include "http_Content.h"
#include "SString_Namespace.h"

#ifdef __cplusplus
extern "C" 
{
#endif

int SetHTML_ContentHeader(const cGH *cctkGH, int choice, String *mess,
                const String *menu);
int SetHTML_ContentFooter(const cGH *cctkGH, int choice, String *mess);

#ifdef __cplusplus
}
#endif

#endif 
