 /*@@
   @header    Content.h
   @date      Sun Sep 17 14:19:23 2000
   @author    Tom Goodale
   @desc 
   Routines exported by the Content stuff.  Was httpd_Content.h; still
   exported by that name.
   @enddesc
   @version $Header$
 @@*/

#ifndef __HTTP_CONTENT_H__
#define __HTTP_CONTENT_H__ 1

#define HTTP_QUICKLINK 1

#include "SString.h"
#include "httpRequest.h"

#ifdef __cplusplus
extern "C" 
{
#endif

int  HTTP_ContentLink(const char *URL, 
                      const char *name,
                      const char *description,
                      int flags);

int  HTTP_ContentSendFromFile(httpRequest *request, int filedes);

int  HTTP_ContentHeader(const cGH *cctkGH, int choice, int len, char *mess,
                      const char *menu);
int  HTTP_ContentFooter(const cGH *cctkGH, int choice, int len, char *mess);

int  HTTP_SetContentHeaderString(const cGH *cctkGH, int choice, String *mess,
                      const String *menu);
int  HTTP_SetContentFooterString(const cGH *cctkGH, int choice, String *mess);

void HTTP_SetHeadInfo( String *header);

void HTTP_SetDoctype( String *header);

#ifdef __cplusplus
}
#endif

#endif /* __HTTP_CONENT_H__ */
