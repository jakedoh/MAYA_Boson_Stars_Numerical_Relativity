 /*@@
   @file      Groups.c
   @date      Wed Sep 14 23:47:43 2000
   @author    Gabrielle Allen
   @desc 
   Pages about groups
   @enddesc
   @version $Header$
 @@*/

#include <stdlib.h>

#include "cctk.h"

#include "httpRequest.h"
#include "Content.h"

#include "SString_Namespace.h"
static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusConnect_HTTPD_Groups_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

static int GroupsPage(const cGH *cctkGH, httpRequest *request, void *data);

/*static int watch[1024];*/

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

int HTTPi_RegisterGroupsPages(void);

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    HTTPi_RegisterGroupPages
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
int HTTPi_RegisterGroupsPages(void)
{
  /* Register the group info page. */
  HTTP_RegisterPage("/Groups", GroupsPage, NULL);

  HTTP_ContentLink("/Groups/index.html", "Groups and Variables",
                   "Information about grid variables and groups",
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
   @routine    GroupsPage
   @date       Thu Sep 14 23:47:43 2000
   @author     Gabrielle Allen
   @desc 
   Displays the group description page.
   @enddesc 
   @calls     
   @calledby   
@@*/
static int GroupsPage(const cGH *cctkGH, httpRequest *request, void *data)
{
  int retval;
  String *message = String_New();
  int i,j;
  int ngroups,nvars;
  cGroup gdata;
  char *groupname;

  /* avoid compiler warning about unused parameter */
  data = data;

  HTTP_SendOKHeader( request );

  /* Start the page */
  HTTP_SetDoctype( message );
  HTTP_SendString(request, message);

  HTTP_Send(request,"<html><head>\n");
  HTTP_Send(request,"<title>Cactus Simulation Group Information</title>\n");

  HTTP_SetHeadInfo( message);
  HTTP_SendString(request, message );

  HTTP_Send(request,"<style type=\"text/css\">\n"
                  "\ttable.groups { background-color: #E9F4D3; } \n"
                  "\t\t.groups td { text-align: left; } \n"
                  "</style>\n");

  HTTP_Send(request,"</head>\n<body>\n");

  /* HTTP_Write out the header part. */

  HTTP_SetContentHeaderString(cctkGH,0,message,NULL);

  retval = HTTP_SendString(request, message);

  ngroups = CCTK_NumGroups();
  
  retval = HTTP_Send(request, "<h1>Groups and Grid Variables</h1>\n"
         "<p>These pages describe the grid variables and groups \n"
         "active in this simulation.</p>\n");

  SetToCString(message,
          "<p>This simulation contains ");
  ConcatDecimal(message, CCTK_NumGroups());
  ConcatCString(message,
          " groups, and ");
  ConcatDecimal(message, CCTK_NumVars());
  ConcatCString(message,
          " variables, "
          "set in ");
  ConcatDecimal(message, CCTK_MaxDim());
  ConcatCString(message,
          "-space dimensions. \nGroups for which storage is currently\n"
          "assigned are written in <span class=\"hilite\">red</span>. \n"
          "The numbers in square brackets are the group and variable indices."
          "</p>\n");
  retval = HTTP_SendString(request, message);

  retval = HTTP_Send(request,"<form action=\"/Groups.html\" method=\"get\">\n");

  retval = HTTP_Send(request,"<div class=\"centered\">\n"
         "<table class=\"groups\" width=\"100%\" cellpadding=\"5\" "
         "cellspacing=\"5\">\n"
         "<tr><th>Groups</th><th>Group Properties</th>"
         "<th>Variables</th></tr>\n");

  for(i=0; i < ngroups; i++)
  {
    SetToCString(message,"<tr>");

    groupname = CCTK_GroupName(i);
    if (CCTK_QueryGroupStorageI(cctkGH,i) > 0)
    {
      ConcatCString(message, " <td>[");
      ConcatDecimal(message, i);
      ConcatCString(message, "] <span class=\"hilite\">");
      ConcatCString(message, groupname);
      ConcatCString(message, "</span></td>\n"); 
    }
    else
    {
      ConcatCString(message, " <td>[");
      ConcatDecimal(message, i);
      ConcatCString(message, "] ");
      ConcatCString(message, groupname);
      ConcatCString(message, "</td>\n");
    }

    free(groupname);

    HTTP_SendString(request, message);

    /* Group Description */

    SetToCString(message, "<td>");

    if (CCTK_GroupData(i,&gdata) > -1)
    {
      switch (CCTK_GroupTypeI(i))
      {
      case CCTK_SCALAR:
        ConcatCString(message,"Grid scalar");
        break;
      case CCTK_ARRAY:
        ConcatCString(message,"Grid array");
        break;
      case CCTK_GF:
        ConcatCString(message,"Grid function");
        break;
      }        
          
      ConcatCString(message, " " );
      ConcatCString(message, CCTK_VarTypeName(gdata.vartype) );
      ConcatCString(message, " <br />\n(");
      ConcatDecimal(message, CCTK_VarTypeSize(gdata.vartype));
      ConcatCString(message, " bytes)");
      if (!(CCTK_GroupTypeI(i) == CCTK_SCALAR))
      {
        ConcatCString(message, " <br />\nDimension ");
        ConcatDecimal(message, gdata.dim);
        ConcatCString(message, " <br />\nTimelevels ");
        ConcatDecimal(message, gdata.numtimelevels);
      }
      HTTP_SendString(request, message);
    }

    HTTP_Send(request,"</td>");

    nvars = CCTK_NumVarsInGroupI(i);
    SetToCString(message,"<td>");
    for(j=CCTK_FirstVarIndexI(i); j < CCTK_FirstVarIndexI(i)+nvars; j++)
    {
      ConcatCString(message, "[" );
      ConcatDecimal(message, j);
      ConcatCString(message, "] " );
      ConcatCString(message, CCTK_VarName(j) );
      ConcatCString(message, "<br />\n" );
    }
    ConcatCString(message,"</td></tr>");
    HTTP_SendString(request, message);

  }

  HTTP_Send(request,"</table></div></form>\n");
    
  /* Write out the footer part. */

  HTTP_SetContentFooterString(cctkGH,0,message);
  retval = HTTP_SendString(request, message);

  String_Delete( message );

  return retval;
}

