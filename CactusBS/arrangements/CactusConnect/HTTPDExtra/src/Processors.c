 /*@@
   @file      Processors.c
   @date      Wed Nov 8 2000
   @author    Gabrielle Allen
   @desc 
   Pages about processors
   @enddesc
   @version $Header$ 
 @@*/

#include "cctk.h"

#include "util_String.h"

#include "httpextra_HostNames.h"
#include "http_Content.h"

/* SW Temporary, while testing the SString module*/
#include "CactusConnect/HTTPD/src/SString.h"
#include "CactusConnect/HTTPD/src/SString_Namespace.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusConnect_HTTPDExtra_Processors_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

static int ProcessorsPage(const cGH *cctkGH, httpRequest *request, void *data);


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
   @routine    HTTPDExtra_RegisterProcessorsPages
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
int HTTPDExtra_RegisterProcessorsPages(void)
{
  /* Register the group info page. */
  HTTP_RegisterPage("/Processors", ProcessorsPage, NULL);

  HTTP_ContentLink("/Processors/index.html", "Processor Information",
                   "Processor layout and properties",
                   HTTP_QUICKLINK);
  return 0;
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/


/******************************************************************************
 ***************************** Groups Page **************************************
 ******************************************************************************/

 /*@@
   @routine    ProcessorsPage
   @date       Thu Sep 14 23:47:43 2000
   @author     Gabrielle Allen
   @desc 
   Displays the processor description page.
   @enddesc 
   @calls     
   @calledby   
@@*/
static int ProcessorsPage(const cGH *cctkGH, httpRequest *request, void *data)
{
  int retval = 0;
  int nprocs = 0,np = 0;
  String *message = String_New();

  HTTP_SendOKHeader( request );

  HTTP_SetDoctype( message );
  HTTP_SendString(request, message);

  /* Start the page */
  HTTP_Send(request, "<html>\n<head>\n");
  HTTP_Send(request, "<title>Cactus Simulation Processor Information</title>\n");

  HTTP_SetHeadInfo( message);
  HTTP_SendString(request, message );
  
  HTTP_Send(request, "<style type=\"text/css\">\n");
  HTTP_Send(request, "  th, td { padding-left: 1em; padding-right: 1em; }\n");
  HTTP_Send(request, "  td.name { text-align: left; }\n");
  HTTP_Send(request, "  td.number { text-align: center; }\n");
  HTTP_Send(request, "</style>\n");
  HTTP_Send(request, "</head>\n<body>\n");

  /* HTTP_Write out the header part. */
  HTTP_SetContentHeaderString(cctkGH, 0, message, NULL);
  retval = HTTP_SendString(request, message);
  
  HTTP_Send(request, "<h1>Processor Information</h1>\n");

  HTTP_Send(request, 
         "<div class=\"centered\">\n<table rules=\"cols\">\n"
         "<tr>"
         "<th>Number</th>\n"
         "<th>Machine Name</th>\n"
         "</tr>");
         
  nprocs = CCTK_nProcs(cctkGH);
  for (np=0;np<nprocs;np++)
  {
    SetToCString( message, "<tr><td class=\"number\">" );
    ConcatDecimal( message, np );
    ConcatCString( message, "</td><td class=\"name\">" );
    ConcatCString( message, HTTPDExtra_RemoteHostName(np) );
    ConcatCString( message, "</td></tr>\n");
    HTTP_SendString(request, message );
  }

  retval = HTTP_Send(request, "</table>\n</div>\n");
    
  HTTP_SetContentFooterString(cctkGH, 0, message);
  retval = HTTP_SendString(request, message);

  String_Delete( message );
  return retval;
}

