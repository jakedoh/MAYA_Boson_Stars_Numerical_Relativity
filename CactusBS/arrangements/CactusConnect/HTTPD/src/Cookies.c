 /*@@
   @file      Cookies.c
   @date      Mon Sep 18 21:08:37 2000
   @author    Tom Goodale
   @desc 
   Cookie stuff.
   @enddesc 
   @version $Header$
 @@*/

#include "cctk.h"

#include <stdlib.h>

#include "util_String.h"

#include "httpRequest.h"
#include "Cookies.h"
#include "SString_Namespace.h"
static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusConnect_HTTPD_Cookies_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    HTTP_CookieSend
   @date       Mon Sep 18 22:42:35 2000
   @author     Tom Goodale
   @desc 
   Sends a cookie to a browser.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
int HTTP_CookieSend(httpRequest *request,
                    const char *name, 
                    const char *value, 
                    const char *path,
                    const char *domain,
                    const char *expires,
                    int secure)
{
  String *message = String_New();

  SetToCString(message, "Set-Cookie: ");
  ConcatCString(message, name);
  ConcatCString(message, "=");
  ConcatCString(message, value);

  if(path)
  {
    ConcatCString(message, "; path=");
    ConcatCString(message, path);
  }

  if(domain)
  {
    ConcatCString(message, "; domain=");
    ConcatCString(message, domain);
  }

  if(expires)
  {
    ConcatCString(message, "; expires=");
    ConcatCString(message, expires);
  }

  if(secure)
  {
    ConcatCString(message, "; secure");
  }    

  ConcatCString(message, "\r\n");

  HTTP_SendString(request, message);

  String_Delete( message );
  return 0;
}

 /*@@
   @routine    HTTP_CookieCancel
   @date       Mon Sep 18 22:43:04 2000
   @author     Tom Goodale
   @desc 
   Cancels a cookie.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
int HTTP_CookieCancel(httpRequest *request,
                      const char *name, 
                      const char *path)
{
  String *message = String_New();

  /* Clear the value */
  SetToCString(message, "Set-Cookie: ");
  ConcatCString(message, name);
  ConcatCString(message, "=");

  if(path)
  {
    ConcatCString(message, "; path=");
    ConcatCString(message, path);
  }

  /* Pick a date in the past */
  ConcatCString(message, "; expires Sun Sep 17 21:57:45 CEST 2000");
  ConcatCString(message, "\r\n");

  HTTP_SendString(request, message);

  String_Delete( message );
  return 0;
}

 /*@@
   @routine    HTTP_CookieGet
   @date       Mon Sep 18 22:43:20 2000
   @author     Tom Goodale
   @desc
               Gets the value of a cookie from a request.
               The allocated return string must be freed by the caller.
   @enddesc
@@*/

char *HTTP_CookieGet(httpRequest *request, const char *cookie_name)
{
  char *cookie_value = NULL;
  String *attribute = String_Make("Cookie");
  String *header = String_New();

  /* Get the cookie header */
  if( HTTP_GetHeaderValueString( request, attribute, header ) )
  {
    size_t value_start = 0;
    String *name = String_Make( cookie_name );
    ConcatCString (name, "=");

    /* Search for "<name>=" */
    if (FindStringFrom (header, name, &value_start))
    {
      /* truncate "<value>" at the next ';' char */
      size_t value_end = value_start + StringLength (name);
      value_start = value_end;
      if (FindCharFrom (header, ';', &value_end))
      {
        StringTruncate (header, value_end);
      }
      cookie_value = Util_Strdup (GetBuffer (header) + value_start);
    }
    String_Delete(name);
  }
  String_Delete(header);
  String_Delete(attribute);

  return cookie_value;
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/

