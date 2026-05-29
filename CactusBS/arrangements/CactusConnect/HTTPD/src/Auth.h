 /*@@
   @header    Auth.h
   @date      Fri Sep 15 13:20:01 2000
   @author    Tom Goodale
   @desc 
   Was httpd_Auth.h; still exported by that name.
   @enddesc
   @version $Header$
 @@*/

#ifndef __HTTP_AUTH_H__
#define __HTTP_AUTH_H__ 1

#include "httpRequest.h"

#ifdef __cplusplus
extern "C" 
{
#endif

int HTTP_AuthAddUser(const char *database, 
                     const char *name,
                     const char *password,
                     const char *encryption_scheme);

int HTTP_AuthenticateBasic(httpRequest *request,
                           const char *database,
                           char *user,
                           int length);

#ifdef __cplusplus
}
#endif

#endif /* __HTTP_AUTH_H__ */
