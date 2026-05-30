 /*@@
   @file      Groups.c
   @date      Wed Sep 24 23:47:43 2000
   @author    Gabrielle Allen
   @desc 
   Pages about thorns
   @enddesc
   @version $Header$
 @@*/

#include <stdio.h>

#include "cctk.h"

#include "util_String.h"

#include "httpRequest.h"
#include "Content.h"
#include "SString_Namespace.h"
#include "SStringHTML_Namespace.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusConnect_HTTPD_Thorns_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

static int ThornMainPage(const cGH *cctkGH, httpRequest *request, void *data);
static int ThornPage(const cGH *cctkGH, httpRequest *request, void *data);


/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

int HTTPi_RegisterThornPages(void);

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    HTTPi_RegisterThornsPages
   @date       Wed Sep 14 11:29:43 2000
   @author     Gabrielle Allen
   @desc 
   Httpd utils registration routine.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
#define THORN_NAME_MAXLENGTH (27+20)
int HTTPi_RegisterThornPages(void)
{
  int i;
  char pagename[THORN_NAME_MAXLENGTH] = { '\0' };

  /* Register the group info page. */
  HTTP_RegisterPage("/Thorns", ThornMainPage, NULL);

  HTTP_ContentLink("/Thorns/index.html", "Thorns",
                   "Information from Flesh and individual thorns",
                   HTTP_QUICKLINK);

  for (i = 0; i < CCTK_NumCompiledThorns (); i++)
  {
    const char *thorn = CCTK_CompiledThorn(i);
    char *namecopy;

    sprintf(pagename,"/Thorns/%s", thorn);

    namecopy = Util_Strdup(thorn); /*SW isn't this a memory leak?*/

    HTTP_RegisterPage(pagename, ThornPage, namecopy);
  }

  return 0;
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/


/******************************************************************************
 *************************** Thorn Page **************************************
 ******************************************************************************/

 /*@@
   @routine    ThornMainPage
   @date       Thu Sep 14 23:47:43 2000
   @author     Gabrielle Allen
   @desc 
   Displays the thorn main page.
   @enddesc 
   @calls     
   @calledby   
@@*/
static int ThornMainPage(const cGH *cctkGH, httpRequest *request, void *data)
{
  int i;
  int retval = -1;
  int foundone = 0;
  const char *thorn;
  String *message = String_New();

  /* avoid compiler warning about unused parameter */
  data = data;

  HTTP_SendOKHeader( request );

  HTTP_SetDoctype( message );
  HTTP_SendString(request, message);
  /* Start the page */
  HTTP_Send(request, "<html><head><title>Cactus Thorns</title>\n");

  HTTP_SetHeadInfo( message);
  HTTP_SendString(request, message );

  HTTP_Send(request,"<style type=\"text/css\">\n"
                  "\t.thorns td { text-align: left; } \n"
                  "</style>\n");
  HTTP_Send(request,"</head>\n<body>\n");

  /* HTTP_SendString out the header part. */

  HTTP_SetContentHeaderString(cctkGH,0,message,NULL);

  retval = HTTP_SendString(request, message);

  retval = HTTP_Send(request, "<h1>Thorns</h1>\n"
         "<p>These pages describe the thorns used in this simulation.</p>\n");

  retval = HTTP_Send(request, "<table><tr><td>\n");

  for (i = 0; i < CCTK_NumCompiledThorns (); i++)
  {
    thorn = CCTK_CompiledThorn (i);
    if (CCTK_IsThornActive (thorn))
    {
      if (!foundone)
      {
        HTTP_Send(request,
               "<h2>Active Thorns</h2>\n"
               "<div class=\"center\">\n"
               "<table class=\"thorns\" cellspacing=\"0\" cellpadding=\"5\">\n"
               "<tr>\n"
               "<th>Thorn Name</th>\n"
               "<th>Implementation</th>\n"
               "</tr>\n");
        foundone++;
      }
      SetToCString(message, "<tr>\n<td><a href=\"/Thorns/");
      ConcatCString(message, thorn);
      ConcatCString(message, "/\">");
      ConcatCString(message, thorn);
      ConcatCString(message, "</a></td>\n<td>");
      ConcatCString(message, CCTK_ThornImplementation(thorn));
      ConcatCString(message, "</td>\n</tr>\n");
      HTTP_SendString(request, message);
    }
  }

  if (foundone)
  {
    HTTP_Send(request,"</table></div>\n");
  }

  retval = HTTP_Send(request,"</td><td>");

  foundone = 0;
  for (i = 0; i < CCTK_NumCompiledThorns (); i++)
  {
    thorn = CCTK_CompiledThorn (i);
    if (!CCTK_IsThornActive (thorn))
    {

      if (!foundone)
      {
        HTTP_Send(request,
               "<h2>Dormant Thorns</h2>\n"
               "<div class=\"centered\">\n"
               "<table class=\"thorns\" cellspacing=\"0\" cellpadding=\"5\">\n"
               "<tr>\n"
               "<th>Thorn Name</th>\n"
               "<th>Implementation</th>\n"
               "</tr>\n");
        foundone++;
      }

      SetToCString(message, "<tr>\n<td>\n");
      ConcatCString(message, thorn);
      ConcatCString(message, "</td>\n<td>");
      ConcatCString(message, CCTK_ThornImplementation(thorn));
      ConcatCString(message, "</td>\n</tr>\n");
      HTTP_SendString(request, message);
      
    }
  }

  if (foundone)
  {
    HTTP_Send(request,"</table>\n</div>\n");
  }

  retval = HTTP_Send(request,"</td>\n</tr>\n</table>\n");

  /* Write out the footer part. */

  HTTP_SetContentFooterString(cctkGH,0,message);
  retval = HTTP_SendString(request, message);

  String_Delete( message );
  return retval;
}



 /*@@
   @routine    ThornPage
   @date       Sun Sep 24 15:13:55 2000
   @author     Gabrielle Allen
   @desc 
   Individual thorn information
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
static int ThornPage(const cGH *cctkGH, httpRequest *request, void *data)
{
  int retval=0;
  String *message = String_New();
  const char *thorn = (const char *)data;
  
  HTTP_SendOKHeader(request);

  HTTP_SetDoctype( message );
  HTTP_SendString(request, message);

  /* Start the page */
  SetToCString(message, "<html>\n<head>\n<title>Thorn Page : ");
  ConcatCString(message, thorn);
  ConcatCString(message, "</title>\n");

  HTTP_SendString(request, message);
  HTTP_SetHeadInfo( message);
  HTTP_SendString(request, message );

  HTTP_Send(request,"</head>\n<body>\n");

  HTTP_SetContentHeaderString(cctkGH,0,message,NULL);

  HTTP_SendString(request, message);

  SetToCString(message, "<h1>Thorn ");
  ConcatCString(message, thorn);
  ConcatCString(message, "</h1>\n");

  HTTP_SendString(request, message);
  
  SetToCString(message,"<p>This page will include all the information about "
                       "thorn \'");
  ConcatCString(message, thorn);
  ConcatCString(message, "\'.\n For now, only information about the parameters"
                         " is given.</p>\n");

  ConcatCString(message, "<ul>\n"
                         "<li><a href=\"/Parameters/");
  ConcatCString(message, thorn );
  ConcatCString(message, "\">Parameters</a></li>\n"
                         "</ul>\n");
  HTTP_SendString(request, message);

  /* Write out the footer part. */
    
  HTTP_SetContentFooterString(cctkGH, 0, message);
  retval = HTTP_SendString(request, message);

  String_Delete( message );
  return retval;
}

