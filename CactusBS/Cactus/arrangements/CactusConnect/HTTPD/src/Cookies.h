 /*@@
   @header    Cookies.h
   @date      Mon Sep 18 22:40:06 2000
   @author    Tom Goodale
   @desc 
   Functions to manipulate cookies.  Was http_Cookies.h; still exported by
   that name
   @enddesc
   @version $Header$
 @@*/

#ifndef __HTTP_COOKIES_H__
#define __HTTP_COOKIES_H__ 1

#include "httpRequest.h"

#ifdef __cplusplus
extern "C" 
{
#endif

int HTTP_CookieSend(httpRequest *request,
                    const char *name, 
                    const char *value, 
                    const char *path,
                    const char *domain,
                    const char *expires,
                    int secure);

int HTTP_CookieCancel(httpRequest *request,
                      const char *name, 
                      const char *path);

char *HTTP_CookieGet(httpRequest *request,
                     const char *name);


#ifdef __cplusplus
}
#endif

#endif /* __HTTP_COOKIES_H__ */

