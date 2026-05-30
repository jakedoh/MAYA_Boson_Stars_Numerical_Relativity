 /*@@
   @file      Content.c
   @date      Wed Sep 13 23:47:43 2000
   @author    Tom Goodale
   @desc 
              The actual content of the web pages.
              Note that this need not be in this thorn, since
              it should just use the interface defined in http_Request.h
              and should never touch the internals.
   @enddesc
   @version   $Id$
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"

#include <stdlib.h>
#include <string.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include "cctk_Version.h"

#include "util_String.h"
#include "util_Network.h"

#include "httpRequest.h"

#include "Auth.h"
#include "Steer.h"
#include "Cookies.h"

#include "cctk_Parameters.h"

#include "SString_Namespace.h"

#include "Content.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusConnect_HTTPD_Content_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

struct httpStaticPage
{
  char *page;
  unsigned int length;
  char *mime_type;
};

struct httpLink
{
  struct httpLink *next;
  char *URL;
  char *name;
  char *description;
  int flags;
};

#define HOSTLENGTH 255

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

static int RegisterImages(void);

static int MainPage(const cGH *cctkGH, httpRequest *request, void *data);
static int AboutPage(const cGH *cctkGH, httpRequest *request, void *data);

static int ShowStaticPage(const cGH *cctkGH, httpRequest *request, void *data);

static int CompareStrings(const void *string1, const void *string2);

static int ControlPage(const cGH *cctkGH, httpRequest *request, void *data);
static int ControlSet(const cGH *cctkGH, httpRequest *request);
static int ControlTerminationPage(const cGH *cctkGH, httpRequest *request);

static int CookieTestPage(const cGH *cctkGH, httpRequest *request, void *data);

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

void HTTP_ContentWork(CCTK_ARGUMENTS);
int HTTP_RegisterPages(void);
int HTTPi_RegisterGroupsPages(void);
int HTTPi_RegisterThornPages(void);
int HTTPi_RegisterParameterPages(void);


/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

struct httpLink *ContentLinks = NULL;

static const char *notauthorized_page =
"<html>\n<title>Error 401: Not Authorized</title>\n</head>\n"
"<body>\nYou are not authorized to access this page\n</body>\n<html>\n";


/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/
int
HTTP_SendString( httpRequest *request, const String * message )
{
  return HTTP_Send(request, GetBuffer( message ) );
}

void
HTTP_SendOKHeader( httpRequest *request )
{
  HTTP_SendOKRefreshHeader( request, -1 );
}

void
HTTP_SendOKRefreshHeader( httpRequest *request, int secs )
{
  /* Status message */
  HTTP_Send(request, "HTTP/1.0 200 OK\r\n");
  /* Content-Type */
  HTTP_Send(request, "Content-Type: text/html\r\n");

  if( secs >= 0 )
  {
    String * timeString = String_Make( "Refresh: " );
    ConcatDecimal( timeString, secs );
    ConcatCString( timeString, "\r\n" );
    HTTP_SendString( request, timeString );
    String_Delete( timeString );
  }
  HTTP_Send( request, "\r\n" );
}

void
SSUtil_SplitFilename(const char **dir,const char **file, const String *message)
{
  Util_SplitFilename((char **)dir,(char **)file,GetBuffer(message));
}

void
CCTK_GetRunTitleString( String *s )
{
  char buf[1024] = EMPTYSTRING;
  CCTK_RunTitle(sizeof(buf)-1,buf);
  SetToCString( s, buf );
}

void
CCTK_GetParameterFilenameString( String *name )
{
  char buf[1024] = EMPTYSTRING;
  CCTK_ParameterFilename(sizeof(buf)-1,buf);
  SetToCString( name, buf );
}
 /*@@
   @routine    HTTP_ContentWork
   @date       Sat Sep 16 15:22:59 2000
   @author     Tom Goodale
   @desc 
   Does any work which needs to be done by the content pages.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
void HTTP_ContentWork(CCTK_ARGUMENTS)
{
#if 0
  HTTP_SteerDispatch();
#endif
}


 /*@@
   @routine    HTTP_RegisterPages
   @date       Wed Sep 13 23:47:43 2000
   @author     Tom Goodale
   @desc 
   Main page registration routine.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
int HTTP_RegisterPages(void)
{
  DECLARE_CCTK_PARAMETERS
#ifdef HTTP_DEBUG
  printf("Registering index.html\n");
#endif

  /* Register the master page. */
  HTTP_RegisterPage("/index.html", MainPage, NULL);

  /* Register parameter control stuff. */
 
  HTTP_RegisterPage("/control.html", ControlPage, NULL);

  HTTP_ContentLink("/control.html", "Cactus Control",
                   "Control Panel for this run",
                   HTTP_QUICKLINK);

  /* Register thorn pages */

  HTTPi_RegisterThornPages();

  /* Register the server description page */
  HTTP_RegisterPage("/About.html", AboutPage, NULL);

  /* Register parameter stuff */
  
  HTTPi_RegisterParameterPages();

  /* Register Groups Pages */
  HTTPi_RegisterGroupsPages();

  HTTP_AuthAddUser("user",user,password,encryption_scheme);

  /* Register images */
  RegisterImages();

  HTTP_RegisterPage("/cookies.html", CookieTestPage, NULL);

  return 0;
}

 /*@@
   @routine    HTTP_ContentLink
   @date       Sun Sep 17 13:17:56 2000
   @author     Tom Goodale
   @desc 
   Register content links.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
int HTTP_ContentLink(const char *URL, 
                     const char *name,
                     const char *description,
                     int flags)
{
  int retval = -1;
  struct httpLink *last;
  struct httpLink *current;
  struct httpLink *hlink = (struct httpLink *)malloc(sizeof(struct httpLink));

  if(hlink)
  {
    hlink->URL         = Util_Strdup(URL);
    hlink->name        = Util_Strdup(name);
    hlink->description = Util_Strdup(description);
    hlink->flags       = flags;
    hlink->next        = NULL;

    /* Put on list */
    if(ContentLinks)
    {
      for(current=ContentLinks; current; current = current->next)
      {
        last = current;
      }
      last->next = hlink;
    }
    else
    {
      ContentLinks = hlink;
    }

    retval = 0;
  }

  return retval;
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/

 /*@@
   @routine    CompareStrings
   @date       Thu Sep 14 18:57:52 2000
   @author     Tom Goodale
   @desc 
   Case independent string comparison to pass to qsort.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
static int CompareStrings(const void *string1, const void *string2)
{
  return Util_StrCmpi(*(const char * const *)string1, *(const char * const *)string2);
}

 /*@@
   @routine    TimeListItem
   @date       10.04.2004
   @author     Steve White
   @desc 
   Reduce repetition in MainPage below
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
static void
TimeListItem( String *message, int time, const char *units )
{
  ConcatCString(message, " <li><span class=\"hilite\">" );
  ConcatDecimal(message, time );
  ConcatCString(message, "</span> " );
  ConcatCString(message, units );
  ConcatCString(message, "</li>\n" );
}

/******************************************************************************
 ***************************** Main Page **************************************
 ******************************************************************************/


 /*@@
   @routine    MainPage
   @date       Wed Sep 13 23:47:43 2000
   @author     Tom Goodale
   @desc 
   Displays the main page.
   @enddesc 
   @calls     
   @calledby   
   @history 
   @hdate Thu Sep 14 10:54:22 2000 @hauthor Tom Goodale
   @hdesc  Copied content format from original http thorn
           developed by Werner Benger, with the aesthetic
           enhancements of Gabrielle Allen, John Shalf and
           Ed Seidel.
   @endhistory 

@@*/

static int MainPage(const cGH *cctkGH, httpRequest *request, void *data)
{

  DECLARE_CCTK_PARAMETERS

  int retval = -1;
  String *message = String_New();
  String *title = String_New();
  String *menu = String_New();
  const char *dir = NULL;
  const char *file = NULL;
  char *this_user = NULL;
  struct httpLink *hlink = NULL;
  int seconds,minutes,hours,days,weeks,months,years,millenia;
  char host[HOSTLENGTH+1] = EMPTYSTRING;

  /* avoid compiler warning about unused parameter */
  data = data;

  HTTP_SendOKRefreshHeader( request, refresh_seconds );

  HTTP_SetDoctype( message );
  HTTP_SendString(request, message);
  /* Start the page */

  HTTP_Send(request, "<html>\n<head>\n");
  HTTP_Send(request, "<title>Running CACTUS Status Information</title>\n");
  HTTP_SetHeadInfo( message);
  HTTP_SendString(request, message );
  HTTP_Send(request, "<style type=\"text/css\">\n"
       "\ttd.authenticate { text-align: center; background-color: #E5FFA2; }\n"
       "</style>\n");

  HTTP_Send(request, "</head>\n<body>\n");
  /* Write out the main header part */

  /* LIST COMPILED THORNS */
  {
    int i;
    int nthorns = CCTK_NumCompiledThorns();
    if( nthorns > 0 )
    {
      const char **thorns = (const char **)malloc(nthorns * sizeof(char *));
      if( thorns != NULL )
      {
        for(i=0; i < nthorns; i++)
          thorns[i] = CCTK_CompiledThorn (i);

        qsort(thorns, nthorns, sizeof(char *), CompareStrings);
  
        SetToCString( menu,"<h3>Active Thorns:</h3>\n");
        for(i=0; i < nthorns; i++)
          if (CCTK_IsThornActive(thorns[i]))
          {
            ConcatCString( menu,"<a href=\"/Thorns/");
            ConcatCString( menu, thorns[i]);
            ConcatCString( menu,"/\">");
            ConcatCString( menu, thorns[i]);
            ConcatCString( menu,"</a><br />\n");
          }
        free(thorns);
      }
    }
  }

  HTTP_SetContentHeaderString(cctkGH, 0, message, menu);
  HTTP_SendString(request, message );

  HTTP_Send(request, 
         "<div class=\"banner\">\n"
         "<img src=\"/Images/wwwcactuscodeorg.jpg\""
         " alt=\"Cactus\" /></div>\n");

  CCTK_GetRunTitleString(title);

  /* Some blurb */
  SetToCString(message, 
         "<div class=\"centered\">\n"
         "<table width=\"60%\">\n"
         "<tr>\n"
         "<td>\n"
         "<h2>");
  Concat(message, title );
  ConcatCString(message, 
         "</h2>\n"
         "<p>This browser is connected to a Cactus simulation which \n"
         "contains a web server thorn. This thorn provides information \n"
         " and control for the simulation.</p>\n"
         "<table cellpadding=\"15\">\n"
         "<tr>\n"
         "<td class=\"authenticate\">\n"
         "<p><strong>Before controlling any features of the simulation, users \n"
         "must <a href=\"/control.html\">authenticate</a>.</strong></p>\n"
         "</td></tr></table>\n"
         "</td>\n"
         "</tr>\n"
         "</table>\n"
         "</div>\n");

  HTTP_SendString(request, message );

  HTTP_Send(request, 
         "<table cellpadding=\"10\">\n"
         "<tr>\n"
         "<td>\n");

  /* AVAILABLE OPTIONS */

  if(ContentLinks)
  {
    SetToCString(message, 
         "<h3>Available options:</h3>\n"
         "<dl>\n");

    for(hlink = ContentLinks; hlink; hlink=hlink->next)
    {
      ConcatCString(message, "<dt> <a href=\"");
      ConcatCString(message, hlink->URL );
      ConcatCString(message, "\">" );
      ConcatCString(message, hlink->name );
      ConcatCString(message, "</a> </dt>\n<dd>" );
      ConcatCString(message, hlink->description );
      ConcatCString(message, "</dd>\n" );
    }
    
    ConcatCString(message, "</dl>\n" );
    HTTP_SendString(request, message );
  }

  HTTP_Send(request, 
         "</td>\n"
         "<td>\n");


  /* CONFIGURATION DETAILS */

  SetToCString(message, 
          "<h3>Simulation:</h3>\n"
          "<ul>\n"
          "<li>Flesh version <span class=\"hilite\"> ");
  ConcatCString(message, CCTK_FullVersion() );
  ConcatCString(message, 
          "</span></li>\n"
          "<li>Flesh compiled on <span class=\"hilite\">");
  ConcatCString(message, CCTK_CompileDate() );
  ConcatCString(message, 
          "</span>\n"
          " at <span class=\"hilite\">");
  ConcatCString(message, CCTK_CompileTime() );
  ConcatCString(message, 
          "</span></li>\n");

  HTTP_SendString(request, message );

  seconds = CCTK_RunTime();
  minutes = seconds/60;
  seconds = seconds-minutes*60;
  hours = minutes/60;
  minutes = minutes-hours*60;
  days = hours/24;
  hours = hours-days*24;
  weeks = days/7;
  days=days-weeks*7;
  months = weeks/4;
  weeks=weeks-months*4;
  years = months/12;
  months=months-years*12;
  millenia=years/1000;
  years=years-millenia*1000;

  SetToCString(message, "<li>Time since start up\n<ul>");
  if (millenia)
  {
    TimeListItem( message, millenia, "millenia" );
  }
  if (years)
  {
    TimeListItem( message, years, "years" );
  }
  if (months)
  {
    TimeListItem( message, months, "months" );
  }
  if (weeks)
  {
    TimeListItem( message, weeks, "weeks" );
  }
  if (days)
  {
    TimeListItem( message, days, "days" );
  }
  if (hours)
  {
    TimeListItem( message, hours, "hours" );
  }
  if (minutes)
  {
    TimeListItem( message, minutes, "minutes" );
  }
  if (seconds)
  {
    TimeListItem( message, seconds, "seconds" );
  }

  ConcatCString(message, "</ul></li>\n<li> Parameter filename <span class=\"parfile\">" );
  HTTP_SendString(request, message );

  CCTK_GetParameterFilenameString(message);
  SSUtil_SplitFilename(&dir,&file,message);
  HTTP_Send(request, file);

  HTTP_Send(request, "</span></li>\n");

  if (cctkGH && cctkGH->cctk_iteration)
  {
    SetToCString(message, "<li>Estimated time per iteration:<ul>\n");
    ConcatCString(message, " <li><span class=\"hilite\">");
    ConcatDouble(message, CCTK_RunTime()/(double)cctkGH->cctk_iteration );
    ConcatCString(message, "</span> seconds</li></ul></li>\n");

    HTTP_SendString(request, message );

    SetToCString(message, "<li>Estimated time to completion:\n<ul>");

    if (cctk_final_time<cctk_initial_time)
    {
      seconds = (cctk_itlast-cctkGH->cctk_iteration)*CCTK_RunTime()/
        cctkGH->cctk_iteration;
    }
    else
    {
      seconds = (cctk_final_time-cctkGH->cctk_time)/cctkGH->cctk_delta_time*
        CCTK_RunTime()/cctkGH->cctk_iteration;
    }

    minutes = seconds/60;
    seconds = seconds-minutes*60;
    hours = minutes/60;
    minutes = minutes-hours*60;
    days = hours/24;
    hours = hours-days*24;
    weeks = days/7;
    days=days-weeks*7;
    months = weeks/4;
    weeks=weeks-months*4;
    years = months/12;
    months=months-years*12;
    millenia=years/1000;
    years=years-millenia*1000;
  

    if (millenia)
    {
      TimeListItem( message, millenia, "millenia" );
    }
    if (years)
    {
      TimeListItem( message, years, "years" );
    }
    if (months)
    {
      TimeListItem( message, months, "months" );
    }
    if (weeks)
    {
      TimeListItem( message, months, "months" );
    }
    if (days)
    {
      TimeListItem( message, days, "days" );
    }
    if (hours)
    {
      TimeListItem( message, hours, "hours" );
    }
    if (minutes)
    {
      TimeListItem( message, minutes, "minutes" );
    }
    if (seconds)
    {
      TimeListItem( message, seconds, "seconds" );
    }
    
    ConcatCString(message, "</ul></li>\n");

    HTTP_SendString(request, message );
  }

  Util_GetHostName(host,HOSTLENGTH);
  host[HOSTLENGTH] = 0;

  if (CCTK_nProcs(cctkGH) == 1)
  {
    SetToCString(message, "<li>Single processor run</li>\n"
                            "<li>Running on <span class=\"hilite\">" );
    ConcatCString(message, host);
    ConcatCString(message, "</span></li>\n");
  }
  else
  {
    SetToCString(message, " <li>Multiprocessor run on " );
    ConcatDecimal(message, CCTK_nProcs(cctkGH));
    ConcatCString(message, " CPUs</li>\n"
                        " <li>Processor 0 running on <span class=\"hilite\">");
    ConcatCString(message, host);
    ConcatCString(message, "</span></li>\n");
  }
  HTTP_SendString(request, message );

  this_user = getenv ("USER");
  if (this_user)
  {
    SetToCString(message, " <li> Started by <span class=\"hilite\">");
    ConcatCString(message, this_user);
    ConcatCString(message, "</span></li>\n");
    HTTP_SendString(request, message );
  }

  HTTP_Send(request, "</ul>\n");

  /* Finish table started by blurb */
  HTTP_Send(request, "</td>\n</tr>\n</table>\n\n");
 
  /* Write out the footer part. */

  HTTP_SetContentFooterString(cctkGH, 0, message);
  retval = HTTP_SendString(request, message );

  String_Delete( message );
  String_Delete( title );
  String_Delete( menu );

  return retval;
}

/******************************************************************************
 ***************************** Images *****************************************
 ******************************************************************************/

#include "images/wwwcactuscodeorg.h"
#include "images/www.h"

 /*@@
   @routine    RegisterImages
   @date       Wed Sep 13 23:47:43 2000
   @author     Tom Goodale
   @desc 
   Registers images.  Images are static pages, so it
   registers a function which takes a data structure and
   displays it.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
static int RegisterImages(void)
{
  struct httpStaticPage *image2 = NULL;
  struct httpStaticPage *image
           = (struct httpStaticPage *)malloc(sizeof(struct httpStaticPage));

  if(image)
  {
    image->page = (char *)malloc(sizeof(wwwcactuscodeorg));
    memcpy(image->page, wwwcactuscodeorg, sizeof(wwwcactuscodeorg));
    image->length = sizeof(wwwcactuscodeorg);
    image->mime_type = (char *)malloc(strlen("image/jpeg")+1);
    strcpy(image->mime_type, "image/jpeg");

#ifdef HTTP_DEBUG
    printf("Registering /Images/wwwcactuscodeorg.jpg\n");
#endif

    HTTP_RegisterPage("/Images/wwwcactuscodeorg.jpg", ShowStaticPage, (void *)image);
  }

  image2  = (struct httpStaticPage *)malloc(sizeof(struct httpStaticPage));

  if(image2)
  {
    image2->page = (char *)malloc(sizeof(www));
    memcpy(image2->page, www, sizeof(www));
    image2->length = sizeof(www);
    image2->mime_type = (char *)malloc(strlen("image/gif")+1);
    strcpy(image2->mime_type, "image/gif");

    HTTP_RegisterPage("/Images/www.gif", ShowStaticPage, (void *)image2);
  }

  return 0;
}

 /*@@
   @routine    ShowStaticPages
   @date       Wed Sep 13 23:47:43 2000
   @author     Tom Goodale
   @desc 
   Displays a static page.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
static int ShowStaticPage(const cGH *cctkGH, httpRequest *request, void *data)
{
  int retval = -1;

  /* avoid compiler warning about unused parameter */
  cctkGH = cctkGH;

  if(data)
  {
    struct httpStaticPage *page = (struct httpStaticPage *)data;
    String *message = String_New();

    HTTP_Send(request, "HTTP/1.0 200 OK\r\n"); 

    SetToCString(message, "Content-Length: ");
    ConcatDecimal(message, page->length );
    ConcatCString(message, "\r\nContent-Type: ");
    ConcatCString(message, page->mime_type );
    ConcatCString(message, "\r\n\r\n");

    HTTP_SendString(request, message);

    retval = HTTP_Write(request, page->page, page->length);

    String_Delete( message );
  }
  
  return retval;
}

#define USER_LENGTH 255
 /*@@
   @routine    ControlPage
   @date       Sun Sep 17 14:37:59 2000
   @author     Tom Goodale
   @desc 
   Cactus Control Panel
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
static int ControlPage(const cGH *cctkGH, httpRequest *request, void *data)
{
  DECLARE_CCTK_PARAMETERS

  int notauthorised = 0;
  char thisuser[USER_LENGTH+1] = EMPTYSTRING;

  data = data; /* avoid compiler warning about unused parameter */

  notauthorised = HTTP_AuthenticateBasic(request, "user", thisuser,
                                         sizeof(thisuser));

  if(!notauthorised)
  {
    /* Ok the person is authorised. */
    if(HTTP_NumArguments( request ) == 0)
    {
      String *message = String_New();
      /* No arguments, so just display the page */
      HTTP_Send(request,"HTTP/1.0 200 Ok\r\n"); 

      HTTP_Send(request,"WWW-Authenticate: Basic realm=\"Cactus Control\"\r\n"); 

      HTTP_CookieSend(request,"user", thisuser, "/", NULL, NULL, 0);
 
      HTTP_Send(request,"Content-Type: text/html\r\n\r\n");

      HTTP_SetDoctype( message );
      HTTP_SendString(request, message);
      /* Start the page */
      HTTP_Send(request, "<html><head>\n");
      HTTP_Send(request, "<title>Cactus Control and Status Page</title>\n");

      HTTP_SetHeadInfo( message);
      HTTP_SendString(request, message );
      HTTP_Send(request, "<style type=\"text/css\">\n"
              "\t.controls td { text-align: left; vertical-align: middle; } \n"
              "</style>\n" );

      HTTP_Send(request, "</head>\n<body>\n");
      
      HTTP_SetContentHeaderString(cctkGH, 0,message,NULL);

      ConcatCString(message, "<h1>Control and Status Page</h1>\n");

      HTTP_SendString(request, message);

      HTTP_Send(request,
             "<p>This page is the control center for interacting with\n"
             " the current simulation. It is possible to steer certain\n"
             " parameters, as well as pause, restart, or terminate the"
             " simulation.</p>\n");

      HTTP_Send(request, 
             "<div class=\"centered\">\n"
             "<form action=\"/control.html\" method=\"get\">\n");

      HTTP_Send(request, 
             "<h4> Run Control </h4>\n"
             "<p> Select if the run should be paused, running normally, "
             "or terminated.\n"
             " You may also single step to the next iteration.</p>\n");

      HTTP_Send(request,
             "<table class=\"controls\" cellspacing=\"5\" cellpadding=\"5\">\n"
             "<tr>\n");
      
      SetToCString(message, 
              "<td><input type=\"radio\" name=\"runstate\" ");
      ConcatCString(message, pause ? "checked=\"checked\"" : "");
      ConcatCString(message, 
              " value=\"PAUSE\" /> PAUSE</td>\n");
      HTTP_SendString(request, message);
      
      SetToCString(message, 
              "<td><input type=\"radio\" name=\"runstate\" ");
      ConcatCString(message, pause ? "" : "checked=\"checked\"");
      ConcatCString(message, 
              " value=\"RUN\" /> RUN</td>\n");
      HTTP_SendString(request, message);
      
      HTTP_Send(request,
             "<td><input type=\"radio\" name=\"runstate\" "
             "value=\"TERMINATE\" /> TERMINATE</td>\n");

      HTTP_Send(request,
             "<td><input type=\"submit\" name=\"step\" "
             "value=\"STEP\" /></td>\n");
      
      HTTP_Send(request,
             "</tr></table>\n"
             "<table>\n"
             "<tr><td><input type=\"submit\" value=\"OK\" /></td>\n"
             "<td> <input type=\"reset\" /></td>\n"
             "</tr>\n"
             "</table>\n");

      HTTP_Send(request, 
           "<h4> Run Until </h4>\n"
           "<p> The following parameters allow you to select an iteration\n"
           "number or physical time at which the code will pause.\n"
           " You may also choose to pause if a particular expression made up\n"
           " of grid scalars, simulation time and iteration is true. \n"
           " Note that even if 'run' is selected above, the settings here have"
           " precedence. </p>\n");

      SetToCString(message,
             "<table>\n");

      ConcatCString(message,"<tr><td>Iteration</td>"
             "<td><input type=\"text\" value=\"");
      ConcatDecimal(message, until_it );
      ConcatCString(message, 
             "\" name=\"iteration\" /></td>\n"
             "<td><input type=\"checkbox\" value=\"yes\" ");
      ConcatCString(message, until_it_active ? "checked=\"checked\"" : "" );
      ConcatCString(message, 
            " name=\"until_it_active\" /></td>\n"
              "</tr>\n");
      ConcatCString(message, 
             "<tr><td>Time</td>"
             "<td><input type=\"text\" value=\"");
      ConcatDouble(message, until_time );
      ConcatCString(message, 
             "\" name=\"time\" /></td>\n"
             "<td><input type=\"checkbox\" value=\"yes\" ");
      ConcatCString(message, until_time_active ? "checked=\"checked\"" : "" );
      ConcatCString(message, 
             " name=\"until_time_active\" /></td>\n"
              "</tr>\n");
      ConcatCString(message, 
             "<tr><td>Expression</td>\n"
             "<td><input type=\"text\" value=\"");
      ConcatCString(message, until_expression );
      ConcatCString(message, 
             "\" name=\"expression\" /></td>\n"
             "<td><input type=\"checkbox\" value=\"yes\" ");
      ConcatCString(message, until_expression_active ? " checked=\"checked\"" : "" );
      ConcatCString(message, 
             " name=\"until_expression_active\" /></td>\n"
              "</tr>\n");

      ConcatCString(message,"</table>\n");
      
      HTTP_SendString(request, message);

      HTTP_Send(request,
             "</form>\n"
             "</div>\n");
      
      /* Write out the footer part. */
      
      HTTP_SetContentFooterString(cctkGH, 0, message);
      HTTP_SendString(request, message);

      String_Delete( message );
    }
    else
    {
      /* Arguments, so control simulation */
      ControlSet(cctkGH, request);
    }
  }
  else
  {
    /* Not authorised */
    HTTP_Send(request,"HTTP/1.0 401 Unauthorized\r\n"); 

    HTTP_Send(request,"WWW-Authenticate: Basic realm=\"Cactus Control\"\r\n"); 

    HTTP_CookieCancel(request,"user", "/");

    HTTP_Send(request,"Content-Type: text/html\r\n\r\n");
  
    HTTP_Send(request, notauthorized_page);
  }

  return 0;
}

 /*@@
   @routine    ControlSet
   @date       Sun Sep 17 14:39:50 2000
   @author     Tom Goodale
   @desc 
   Set the status of the simulation based on the controls.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
static int ControlSet(const cGH *cctkGH, httpRequest *request)
{
  DECLARE_CCTK_PARAMETERS;

  String *message = String_New();
  const char *value = HTTP_ArgumentValue(request,"runstate");

  if(value)
  {
    switch(*value)
    {
      case 'T' : HTTP_SteerQueue(CCTK_THORNSTRING, "terminate", "yes");
                 if(pause)
                 {
                   HTTP_SteerQueue(CCTK_THORNSTRING, "pause", "no");
                 }
                 ControlTerminationPage(cctkGH, request);
                 break;
      case 'P' : HTTP_SteerQueue(CCTK_THORNSTRING, "pause", "yes");
                 break;
      case 'R' : HTTP_SteerQueue(CCTK_THORNSTRING, "pause", "no");
                 break;
      default  :
        CCTK_VWarn(1, __LINE__,__FILE__,CCTK_THORNSTRING,
                  "Unknown runstate \"%s\"", value );
    }
  }

  /* Is single stepping switched on ? */
  value = HTTP_ArgumentValue(request,"step");
  
  if(value)
  {
    HTTP_SteerQueue(CCTK_THORNSTRING, "single_step", "yes");
  }

  /****************************************************************************
   *************** Is running until an iteration switched on ? ****************
   */
  value = HTTP_ArgumentValue(request,"until_it_active");
  
  if(value)
  {
    /* Value exists, so must be yes */
    if(!until_it_active)
    {
      HTTP_SteerQueue(CCTK_THORNSTRING, "until_it_active", "yes");
    }
  }
  else
  {
    /* Value doesn't exist, so must be no */
    if(until_it_active)
    {
      HTTP_SteerQueue(CCTK_THORNSTRING, "until_it_active", "no");
    }
  }    

  /* Now look at the value and steer it if it's changed */
  value = HTTP_ArgumentValue(request,"iteration");

  if(value)
  {
    if(atoi(value) != until_it)
    {
      HTTP_SteerQueue(CCTK_THORNSTRING, "until_it", value);
    }
  }

  /****************************************************************************
   *************** Is running until a particular time switched on ? ***********
   */
  value = HTTP_ArgumentValue(request,"until_time_active");
  
  if(value)
  {
    /* Value exists, so must be yes */
    if(!until_time_active)
    {
      HTTP_SteerQueue(CCTK_THORNSTRING, "until_time_active", "yes");
    }
  }
  else
  {
    /* Value doesn't exist, so must be no */
    if(until_time_active)
    {
      HTTP_SteerQueue(CCTK_THORNSTRING, "until_time_active", "no");
    }
  }    

  /* Now look at the value and steer it if it's changed */
  value = HTTP_ArgumentValue(request,"time");

  if(value)
  {
    if(atof(value) != until_time)
    {
      HTTP_SteerQueue(CCTK_THORNSTRING, "until_time", value);
    }
  }
  
  /****************************************************************************
   *************** Is running until an expression is true switched on ? *******
   */
  value = HTTP_ArgumentValue(request,"until_expression_active");
  
  if(value)
  {
    /* Value exists, so must be yes */
    if(!until_expression_active)
    {
      HTTP_SteerQueue(CCTK_THORNSTRING, "until_expression_active", "yes");
    }
  }
  else
  {
    /* Value doesn't exist, so must be no */
    if(until_expression_active)
    {
      HTTP_SteerQueue(CCTK_THORNSTRING, "until_expression_active", "no");
    }
  }    

  /* Now look at the value and steer it if it's changed */
  value = HTTP_ArgumentValue(request,"expression");

  if(value)
  {
    if(!CCTK_Equals(value,until_expression))
    {
      HTTP_SteerQueue(CCTK_THORNSTRING, "until_expression", value);
    }
  }

  /****************************************************************************
   **************** Now redirect the browser to the normal page ***************
   */
  /* SW: I don't think this is working right.  On my browsers, it just
   * displays the ControlTerminationPage, with the following lines
   * of text rendered at the bottom. */

  /* Status message */
  if(HTTP_MajorVersion( request ) < 1 || 
     (HTTP_MajorVersion( request ) == 1 && HTTP_MinorVersion( request ) < 1))
  {
    /* Older browsers don't understand 303 */
    SetToCString(message,"HTTP/1.0 302 Found\r\n");
  }
  else
  {
    SetToCString(message,"HTTP/1.0 303 See Other\r\n");
  }
  
  ConcatCString(message, "Location: /control.html\r\n\r\n");
  
  HTTP_SendString(request, message);

  String_Delete( message );

  return 0;
}

 /*@@
   @routine    ControlTerminationPage
   @date       Sun Sep 17 14:40:16 2000
   @author     Tom Goodale
   @desc 
   Page to be shown on termination of Cactus.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
static int ControlTerminationPage(const cGH *cctkGH, httpRequest *request)
{
  int retval = -1;
  String *message = String_New();
  HTTP_SendOKHeader( request );

  HTTP_SetDoctype( message );
  HTTP_SendString(request, message);
  /* Start the page */
  HTTP_Send(request, "<html><head>\n");
  HTTP_Send(request, "<title>Running CACTUS Status Information : Terminated</title>\n");

  HTTP_SetHeadInfo( message);
  HTTP_SendString(request, message );

  HTTP_Send(request, "</head>\n<body>\n");

  /* Write out the main header part */
  HTTP_SetContentHeaderString(cctkGH,1,message,NULL);
  HTTP_SendString(request, message);

  HTTP_Send(request, "<h1>Simulation Home Page</h1>\n");

  /* Some blurb */
  HTTP_Send(request,
       "<div class=\"centered\">\n"
       "<table cellspacing=\"5\" cellpadding=\"5\" border=\"0\"><tr><td>\n"
       "<h3>Simulation web server:</h3>\n" 
       "<p>This browser is connected to a Cactus simulation which \n"
       "contains a web server thorn. This thorn allows you to monitor \n"
       "the simulation,\n"
       "and view and change parameters</p>\n"
       "<p>Depending on which other thorns are active, there may be \n"
       "additional features available, such as the viewing and \n"
       "downloading of output files</p>\n");

  /* CONFIGURATION DETAILS */

  SetToCString(message,
         "<h3>Simulation:</h3>\n"
         "<ul> <li>Flesh version <span class=\"hilite\"> ");
  ConcatCString(message, CCTK_FullVersion() );
  ConcatCString(message,
          "</span></li>\n"
          "<li>Flesh compiled on <span class=\"hilite\"> " __DATE__" </span>\n"
          "at <span class=\"hilite\"> "__TIME__" </span></li>\n");

  HTTP_SendString(request, message);

  if (cctkGH)
  {
    if (CCTK_nProcs(cctkGH) == 1)
    {
      SetToCString(message,"<li>Single processor run</li>\n");
    }
    else
    {
      SetToCString(message," <li>Multiprocessor run on ");
      ConcatDecimal( message, CCTK_nProcs(cctkGH) );
      ConcatCString( message, " CPUs</li>\n");
    }
    HTTP_SendString(request, message);
  }

  HTTP_Send(request,"</ul>\n");

  /************************************************************************/

  /* NEW COLUMN */

  HTTP_Send(request, "</td><td>");

  /************************************************************************/

  /* CURRENT STATE OF SIMULATION */

  if (cctkGH)
  {

    SetToCString(message, 
                     "<h3> Current state:</h3> \n"
                     "<ul><li>Physical time <span class=\"hilite\"> ");
    ConcatDouble(message, cctkGH->cctk_time );
    ConcatCString(message, 
                     "</span></li> \n"
                     "<li>Iteration number <span class=\"hilite\"> ");
    ConcatDecimal(message, cctkGH->cctk_iteration );
    ConcatCString(message, "</li>\n");
  }
  else
  {
    SetToCString(message, "<li>Current cactus state is unknown</li>\n");
  }

  HTTP_SendString(request, message);
  HTTP_Send(request, "<li>This Cactus run is over.</li>\n");

  HTTP_Send(request, "</ul>");

  /* LIST COMPILED THORNS */
  {
    int i;
    int nthorns = CCTK_NumCompiledThorns();
    const char **thorns = (const char **)malloc(nthorns * sizeof(char *));

    for(i=0; i < nthorns; i++)
    {
      thorns[i] = CCTK_CompiledThorn (i);
    }

    /* Sort the thorns */
    qsort(thorns, nthorns, sizeof(char *), CompareStrings);

    SetToCString(message, 
       "<h3>Compiled thorns:</h3>\n"
       "<p>This list shows all thorns compiled into the executable; those \n"
       "thorns which have been activated in the parameter file for this \n"
       "simulation are shown in <span class=\"hilite\">red</span></p>\n"
       "<ul>\n");

    for(i=0; i < nthorns; i++)
    {
      ConcatCString(message, "<li>");
      if (CCTK_IsThornActive(thorns[i]))
      {
        ConcatCString(message, "  <span class=\"hilite\">");
        ConcatCString(message, thorns[i]);
        ConcatCString(message, " </span>");
      }
      else
      {
        ConcatCString(message, thorns[i]);
      }
      ConcatCString(message, "</li>\n");
    }
    ConcatCString(message, "</ul>\n");

    free(thorns);
  }

  HTTP_SendString(request, message);
  

  /* Finish table started by blurb */
  HTTP_Send(request, "</td></tr></table>\n");
 
  /* Write out the footer part. */

  HTTP_SetContentFooterString(cctkGH,0,message);
  retval = HTTP_SendString(request, message);

  String_Delete( message );
  return retval;
}

 /*@@
   @routine    HTTP_ContentSendFromFile
   @date       Sun Sep 17 17:35:57 2000
   @author     Tom Goodale
   @desc 
   Reads data from the filedescriptor and sends it to the HTTP request.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
int HTTP_ContentSendFromFile(httpRequest *request, int filedes)
{
  int bytes_sent = 0;
  int n_bytes = 0;
  char buffer[4098] = EMPTYSTRING;

  bytes_sent = 0;
  while((n_bytes = read(filedes, buffer,sizeof(buffer))) > 0)
  {
    HTTP_Write(request, buffer, n_bytes);
    bytes_sent += n_bytes;
  }

  if(n_bytes == -1)
  {
    bytes_sent = -bytes_sent;
  }

  return bytes_sent;
}


 /*@@
   @routine    AboutPage
   @date       Sun Sep 17 2000
   @author     Gabrielle Allen
   @desc 
   Displays a page about the web server
   @enddesc 
   @calls     
   @calledby   
@@*/
static int AboutPage(const cGH *cctkGH, httpRequest *request, void *data)
{
  int retval = -1;
  String *message = String_New();

  /* avoid compiler warning about unused parameter */
  data = data;

  HTTP_SendOKHeader( request );

  HTTP_SetDoctype( message );
  HTTP_SendString(request, message);
  /* Start the page */
  SetToCString(message, "<html><head>\n<title>About Cactus Server</title>\n");
  HTTP_SetHeadInfo( message);
  HTTP_SendString(request, message );

  HTTP_Send(request, "</head>\n<body>\n");
  HTTP_SetContentHeaderString(cctkGH,0,message,NULL);
  HTTP_SendString(request, message);

  HTTP_Send(request, "<h1>About this Web Server</h1>\n");

  SetToCString(message, "<p>These web pages are served by a simulation \n"
         "which is using the Cactus Code and Computational ToolKit, \n"
         "a freely available, parallel, collaborative,"
         " portable and \n"
         "modular programming environment for HPC.</p>\n"
         "<p>The HTTPD module, or <dfn>thorn</dfn> "
         "which is serving these pages\n"
         " can be added to any Cactus application to provide on-line \n"
         "monitoring and control of simulations from any web browser.</p>\n"
         "<p>This HTTPD server and thorn interface has been designed and \n"
         "and implemented by Tom Goodale, based on the original idea and \n"
         "implementation by Werner Benger.</p>\n");

  retval = HTTP_SendString(request, message);

  SetToCString(message, "<p>For more information about Cactus, visit our "
         "permanent home page at \n"
         "<a href=\"http://www.cactuscode.org\">www.cactuscode.org</a></p>\n");

  retval = HTTP_SendString(request, message);

  /* Write out the footer part. */
  
  HTTP_SetContentFooterString(cctkGH,0,message);
  retval = HTTP_SendString(request, message);

  String_Delete( message);
  return retval;
}


 /*@@
   @routine    CookieTestPage
   @date       Mon Sep 18 23:28:19 2000
   @author     Tom Goodale
   @desc 
   Test and example page for cookies.  This will disappear
   soon.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
static int CookieTestPage(const cGH *cctkGH, httpRequest *request, void *data)
{
  int retval = -1;
  String *message = String_New();
  const char *value = NULL;
  char *value2 = NULL;

  /* avoid compiler warning about unused parameter */
  data = data;

  /* Status message */
  HTTP_Send(request,"HTTP/1.0 200 OK\r\n");

  /* Cookie */
  HTTP_CookieSend(request, "user1", "foobar4", NULL,NULL,NULL,0);

  HTTP_CookieSend(request, "user2", "foobar3", NULL,NULL,NULL,0);

  HTTP_CookieSend(request, "user3", "foobar2", NULL,NULL,NULL,0);

  HTTP_CookieSend(request, "user4", "foobar1", NULL,NULL,NULL,0);

  HTTP_Send(request,"Content-Type: text/html\r\n\r\n");

  HTTP_SetDoctype( message );
  HTTP_SendString(request, message);
  /* Start the page */
  HTTP_Send(request, "<html><head><title>Cookie Test</title>\n");
  HTTP_SetHeadInfo( message);
  HTTP_SendString(request, message );

  HTTP_Send(request, "</head>\n<body>\n");

  HTTP_SetContentHeaderString(cctkGH,0,message,NULL);
  HTTP_SendString(request, message);

  HTTP_Send(request, "<h1>Cookie Test</h1>\n");

  HTTP_Send(request, "<div class=\"centered\">");

  value = HTTP_HeaderValue(request, "Cookie");

  SetToCString(message, "<p>Cookie was '");
  ConcatCString(message, value);
  ConcatCString(message, "'</p>\n");

  HTTP_SendString(request, message);

  value2 = HTTP_CookieGet(request,"user3");

  SetToCString(message, "<p>Cookie from decoder was '");
  ConcatCString(message, value2);
  ConcatCString(message, "'</p>\n");

  free(value2);

  HTTP_SendString(request, message);

  HTTP_Send(request,"</div>\n");
  
  /* Write out the footer part. */

  HTTP_SetContentFooterString(cctkGH,0,message);
  retval = HTTP_SendString(request, message);

  String_Delete( message );
  return retval;
}
