 /*@@
   @file      Parameters.c
   @date      Wed Sep 13 23:47:43 2000
   @author    Tom Goodale
   @desc 
   Files for displaying and steering parameters
   @enddesc
   @version $Header$
 @@*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include "cctk.h"

#include "util_String.h"

#include "httpRequest.h"

#include "Auth.h"
#include "Steer.h"

#include "Content.h"

#include "cctk_Parameters.h"

#include "SString_Namespace.h"
#include "SStringHTML_Namespace.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusConnect_HTTPD_Parameters_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

static int MainParameterPage(const cGH *cctkGH, httpRequest *request, void *data);
static int ThornParameterPage(const cGH *cctkGH, httpRequest *request, void *data);
static int ParameterPage(const cGH *cctkGH, httpRequest *request, void *data);


/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

int HTTPi_RegisterParameterPages(void);

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

static const char *notauthorized_page =
"<HTML>\n<HEAD><TITLE>Error 401: Not Authorized</TITLE></HEAD>\
<BODY>You are not authorized to access this page</BODY>\n<HTML>\n";


#define USER_LENGTH 255
/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    HTTPi_RegisterParameterPages
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
int HTTPi_RegisterParameterPages(void)
{
  DECLARE_CCTK_PARAMETERS
  
  int i;
  char pagename[27+20+100]; /* Thorns have maximum length
                             then added 100 for parameters */
  union
  {
    const cParamData *pData;
    void *non_const_pData;
  } u;

  /* Two ways to do this - can either just have one function
   * registered as /Parameters which then checks request->residual,
   * to see what the actual accessed page is, or can register a 
   * main page and one for each thorn.
   * Choose the latter for the mo to keep functions smaller, and also
   * since the compiled thorn list is static at the moment.
   */

  HTTP_RegisterPage("/Parameters/index.html", MainParameterPage, NULL);
  
  HTTP_ContentLink("/Parameters/index.html", "Parameters",
                   "Parameter Information and Control",
                   HTTP_QUICKLINK);

  for (i = 0; i < CCTK_NumCompiledThorns (); i++)
  {
    const char *thorn = CCTK_CompiledThorn(i);
    char *namecopy;
    int first = 1;

    sprintf(pagename,"/Parameters/%s", thorn);

    namecopy = Util_Strdup(thorn);

    HTTP_RegisterPage(pagename, ThornParameterPage, namecopy);

    /* Walk through all parameters of given implementation. */
      while(CCTK_ParameterWalk(first, thorn, NULL, &u.pData) == 0)
      {
        first = 0;
        sprintf(pagename,"/Parameters/%s/%s",thorn,u.pData->name);
        HTTP_RegisterPage(pagename, ParameterPage, u.non_const_pData);
      }
  }

  return 0;
}


/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/


static void
SendHTTP_Redirect_Header( httpRequest *request )
{
  /* Now redirect the browser to the normal page */
  /* Status message */
  if(HTTP_MajorVersion( request ) < 1 || 
     (HTTP_MajorVersion( request ) == 1 && HTTP_MinorVersion( request ) < 1))
  {
    /* Older browsers don't understand 303 */
    HTTP_Send(request,"HTTP/1.0 302 Found\r\n");
  }
  else
  {
    HTTP_Send(request,"HTTP/1.0 303 See Other\r\n");
  }
}

static void
SendHTTP_Uauthorized_Header( httpRequest *request )
{
  HTTP_Send(request,"HTTP/1.0 401 Unauthorized\r\n"); 
    
  HTTP_Send(request,"WWW-Authenticate: Basic realm=\"foo\"\r\n"); 
      
  HTTP_Send(request,"Content-Type: text/html\r\n\r\n");
      
}

/******************************************************************************
 ***************************** Parameter Pages ********************************
 ******************************************************************************/


 /*@@
   @routine    MainParameterPage
   @date       Wed Sep 13 23:47:43 2000
   @author     Tom Goodale
   @desc 
   Displays the parameter page.
   @enddesc 
   @calls     
   @calledby   
   @history 
   @hdate Sat Sep 16 15:13:16 2000 @hauthor Tom Goodale
   @hdesc  Copied content format from original http thorn
           developed by Werner Benger, with the aesthetic
           enhancements of Gabrielle Allen, John Shalf and
           Ed Seidel.
   @endhistory 
@@*/
static int MainParameterPage(const cGH *cctkGH, httpRequest *request, void *data)
{
  int retval = -1;
  String *message = String_New();
  int i;

  /* avoid compiler warning about unused parameter */
  data = data;

  HTTP_SendOKHeader( request );

  HTTP_SetDoctype( message );
  HTTP_SendString(request, message);
  /* Start the page */
  HTTP_Send(request, "<html><head><title>Cactus Parameters Request</title>\n");
  HTTP_SetHeadInfo( message);
  HTTP_SendString(request, message );
  Truncate( message, 0 );

  HTTP_Send(request, "\n</head>\n<body>\n");

  HTTP_SetContentHeaderString(cctkGH,0,message,NULL);
  HTTP_SendString(request, message);

  HTTP_Send(request,
         "<h1>Check/Modify Parameters</h1>\n"
         "<p>From this page you can check the values of all parameters for \n"
         "the simulation, and modify any parameters which have been \n"
         "designated as <i>steerable</i></p>\n"
         "<p>Parameters can be viewed for all <i>Active Thorns</i>, that is,\n"
         "for thorns which have been activated in the parameter file for the\n"
         " simulation. \n"
         "Select one of the active thorns for this simulation from the list \n"
         "below to view all of its parameters</p>\n"
         "<p>Steerable parameters can be identified by the presence of a \n"
         "form input box, to change the value of a parameter, simply edit \n"
         "the value in the box and press the submit button to register the \n"
         "new values.</p>\n"
         "<div class=\"centered\">\n"
         "<table cellspacing=\"5\" cellpadding=\"5\">\n"
         "<tr><th>Thorn Name</th><th>Implementation</th></tr>\n");

  for (i = 0; i < CCTK_NumCompiledThorns (); i++)
  {
    const char *thorn = CCTK_CompiledThorn (i);
    
    if (CCTK_IsThornActive (thorn))
    {
      SetToCString(message, "<tr>\n<td><a href=\"/Parameters/");
      ConcatCString(message, thorn);
      ConcatCString(message, "/\">");
      ConcatCString(message, thorn);
      ConcatCString(message, "</a></td>\n<td>");
      ConcatCString(message, CCTK_ThornImplementation(thorn));
      ConcatCString(message, "</td>\n</tr>\n");
      HTTP_SendString(request, message);
      
    }
  }

  HTTP_Send(request,"</table></div>");

  /* Write out the footer part. */
  
  HTTP_SetContentFooterString(cctkGH, 0, message);
  retval = HTTP_SendString(request, message);

  String_Delete( message );
  return retval;
}

 /*@@
   @routine    ThornParameterPage
   @date       Sat Sep 16 15:13:55 2000
   @author     Tom Goodale
   @desc
   Parameter page for a thorn.
   Checks to see if it is a setting request or 
   an info request and behaves accordingly.
 
   Copied content format from original http thorn
   developed by Werner Benger, with the aesthetic
   enhancements of Gabrielle Allen, John Shalf and
   Ed Seidel.

   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
static int ThornParameterPage(const cGH *cctkGH, httpRequest *request, void *data)
{
  int retval=0;
  int i;
  String * message = String_New();
  String * menu = String_New();
  String * temp = String_New();
  const char *thorn = (const char *)data;
  const cParamData *pData;
  char user[USER_LENGTH+1] = EMPTYSTRING;


  int notauthorised = HTTP_AuthenticateBasic(request, "user", user,
                                            USER_LENGTH);
  int readonly = notauthorised;
  
  if(HTTP_NumArguments( request ) > 0)
  {
    /* This is a parameter set request */

    if(!notauthorised)
    {
      if(!readonly)
      {
        const httpArg *argument;
        /* Queue parameters for steering */
        int first = 1;
        while((argument = HTTP_ArgumentWalk(request, first)) != NULL)
        {
          first = 0;
          CCTK_VWarn(5, __LINE__,__FILE__,CCTK_THORNSTRING,
                     "Setting %s::%s%s",
                     thorn,
                     HTTP_ArgName( argument),
                     HTTP_ArgValue( argument ) );
          HTTP_SteerQueue(thorn, HTTP_ArgName( argument),
                          HTTP_ArgValue( argument ));
        }
      }
      SendHTTP_Redirect_Header( request );

      ConcatCString(message, "Location: /Parameters/");
      ConcatCString(message, thorn);
      ConcatCString(message, "/\r\n\r\n");

      HTTP_SendString(request, message);
    }
    else
    {
      SendHTTP_Uauthorized_Header( request );
  
      HTTP_Send(request, notauthorized_page);
    }

  } 
  else
  {    
    /* Display the page. */
    /* Status message */
    HTTP_SendOKHeader( request );

    HTTP_SetDoctype( message );
    HTTP_SendString(request, message);

    /* Start the page */
    SetToCString(message, "<html><head>\n<title>Cactus Parameters Request : ");
    ConcatCString(message, thorn);
    ConcatCString(message, "</title>\n");
    HTTP_SendString(request, message );

    HTTP_SetHeadInfo( message);
    HTTP_SendString(request, message );

    HTTP_Send(request, "<style type=\"text/css\">\n"
            "\t.paramtable td { text-align: left; vertical-align: middle; } \n"
            "\t.paramtable td.description { font-size: small; } \n"
            "</style> \n");

    HTTP_Send(request, "</head>\n<body>\n");

    Truncate( message, 0 );

    if (CCTK_NumCompiledThorns()>0)
    {
      SetToCString(menu, "<h3>Parameters:</h3>\n");
    }
    /* Menu for this page */
    for (i = 0; i < CCTK_NumCompiledThorns (); i++)
    {
      const char *menuthorn = CCTK_CompiledThorn (i);      
      if (CCTK_IsThornActive (menuthorn))
      {
        ConcatCString(menu, " <a href=\"/Parameters/");
        ConcatCString(menu, menuthorn);
        ConcatCString(menu, "/\">");
        ConcatCString(menu, menuthorn);
        ConcatCString(menu, "</a><br />\n");
      }
    }
    HTTP_SetContentHeaderString(cctkGH,0,message,menu);
    HTTP_SendString(request, message);

    if (!CCTK_IsThornActive(thorn))
    {
      SetToCString(message, "<strong> Thorn ");
      ConcatCString(message, thorn);
      ConcatCString(message, " is not active !!!</strong><br />\n");
      HTTP_SendString(request, message);
    }
    else
    {
      int nfixed=0;
      int nsteerable=0;
      int first = 1;
      /* Send table of available parameters for given thorn */
      /* Steerable parameters can be edited in a FORM. */
      SetToCString(message,
              "<h1>Check/Modify Parameters</h1>\n"
              "<h2>Thorn ");
      ConcatCString(message, thorn);
      ConcatCString(message,
              "</h2>\n");
      ConcatCString(message,
        "<p>Parameters in Cactus can be either <dfn>fixed</dfn> or <dfn>steerable</dfn>. \n"
        "Steerable parameters are those which can be modified during a simulation.\n"
        "To change steerable parameters, edit the values and press "
           " the submit button.\n"
        "Before applying the new parameter value, Cactus will check that the\n"
        "parameter lies in the allowed range. \nNote that there"
        " is currently no log of parameters which have been modified.</p>\n"
        "<p>The tables below show the parameter name, current value and description.\n"
        "The default value of each parameter is shown in brackets at the end of the description.</p>"
             "\n");

      HTTP_SendString(request, message);
      HTTP_Send(request, "<div class=\"centered\">\n");
      
      if(!readonly )
      {
        SetToCString(message,"<form action=\"/Parameters/");
        ConcatCString(message, thorn);
        ConcatCString(message, "/\">\n");
        HTTP_SendString(request, message);
      }

      /* Walk through all steerable parameters of given implementation. */
      
      while(CCTK_ParameterWalk(first, thorn, NULL, &pData) == 0)
      {
        char *value = CCTK_ParameterValString (pData->name, pData->thorn);
        Truncate(message,0);
        first = 0;

        if(value)
        {          
          if (pData->steerable == CCTK_STEERABLE_ALWAYS)
          {

            if (nsteerable == 0)
            {
              ConcatCString(message,"<h2>Steerable Parameters</h2>\n");
              if (readonly)
              {
                ConcatCString(message,
                       "<p><em>The following parameters are steerable, "
                       "but you do not have authorisation to change them. \n"
                       "To change parameters you must first register on the "
                       "<a href=\"/control.html\">Simulation Control Page</a>."
                       "</em></p>\n");
              }
              ConcatCString(message,
                        "<table class=\"paramtable\" cellpadding=\"5\" cellspacing=\"5\">\n");
            }
            nsteerable++;

            if (!readonly)
            {
              int param_type;
              if (pData->type == PARAMETER_BOOLEAN)
              {
                /* Steerable boolean */
                int param_bool = 
                  *((const CCTK_INT *)CCTK_ParameterGet(pData->name,thorn,&param_type));
                ConcatCString(message, "<tr>\n<td><a href=\"/Parameters/");
                ConcatCString(message, pData->thorn);
                ConcatCString(message, "/");
                ConcatCString(message, pData->name);
                ConcatCString(message, "\">");
                ConcatCString(message, pData->name);
                ConcatCString(message, "</a></td>\n"
                        "<td>"
                        "Yes \n<input type=\"radio\" name=\"");
                ConcatCString(message, pData->name);
                ConcatCString(message, "\" ");
                ConcatCString(message, param_bool ? "checked=\"checked\"" : "");
                ConcatCString(message, " value=\"1\" />"
                        "&nbsp;"
                        "No \n<input type=\"radio\" name=\"");
                ConcatCString(message, pData->name);
                ConcatCString(message, "\" ");
                ConcatCString(message, param_bool ? "" : "checked=\"checked\"");
                ConcatCString(message, " value=\"0\" />"
                        "</td>\n"
                        "<td class=\"description\">");
                SetToEncodedHTMLCString( temp, pData->description );
                Concat(message, temp );
                ConcatCString(message, " (");
                ConcatCString(message, pData->defval);
                ConcatCString(message, ")</td>\n"
                        "</tr>\n");
              }
              else if (pData->type == PARAMETER_KEYWORD)
              {
                t_range *range;
                /* Steerable keyword */
                CCTK_ParameterGet(pData->name,thorn,&param_type);

                ConcatCString(message,"<tr>\n<td><a href=\"/Parameters/");
                ConcatCString(message, pData->thorn);
                ConcatCString(message, "/");
                ConcatCString(message, pData->name);
                ConcatCString(message, "\">");
                ConcatCString(message, pData->name);
                ConcatCString(message, "</a></td>\n");
                ConcatCString(message,"<td> <select name=\"");
                ConcatCString(message, pData->name);
                ConcatCString(message, "\" size=\"1\">\n");
                for(range = pData->range; range ; range = range->next)
                {
                  if (CCTK_Equals(value,range->range))
                  {
                    ConcatCString(message, "<option selected>");
                    ConcatCString(message, range->range);
                    ConcatCString(message, "\n");
                  }
                  else
                  {
                    ConcatCString(message, "<option>");
                    ConcatCString(message, range->range);
                    ConcatCString(message, "\n");
                  }
                }                
                ConcatCString(message,"</select>\n");
                ConcatCString(message,"</td>\n<td class=\"description\">");
                SetToEncodedHTMLCString( temp, pData->description );
                Concat(message, temp);
                ConcatCString(message, " (");
                SetToEncodedHTMLCString( temp, pData->defval );
                Concat(message, temp);
                ConcatCString(message, ")</td>\n</tr>\n");
              }
              else
              {
                /* Steerable nonboolean */
                ConcatCString(message, "<tr>\n<td><a href=\"/Parameters/");
                ConcatCString(message, pData->thorn);
                ConcatCString(message, "/");
                ConcatCString(message, pData->name);
                ConcatCString(message, "\">");
                ConcatCString(message, pData->name);
                ConcatCString(message, "</a></td>\n"
                                       "<td><input type=\"text\" name=\"");
                ConcatCString(message, pData->name);
                ConcatCString(message, "\" value=\"");
                ConcatCString(message, value);
                ConcatCString(message, "\" /></td>\n"
                                       "<td class=\"description\">");
                SetToEncodedHTMLCString( temp, pData->description );
                Concat(message, temp);
                ConcatCString(message, " (");
                SetToEncodedHTMLCString( temp, pData->defval );
                Concat(message, temp);
                ConcatCString(message, ")</td>\n</tr>\n");
              }
            }
            else
            {
              /* Steerable but no authority */
              ConcatCString(message, "<tr><td><a href=\"/Parameters/");
              ConcatCString(message, pData->thorn);
              ConcatCString(message, "/");
              ConcatCString(message, pData->name);
              ConcatCString(message, "\">");
              ConcatCString(message, pData->name);
              ConcatCString(message, "</a></td>\n<td>");
              ConcatCString(message, value);
              ConcatCString(message, "</td><td class=\"description\">");
              SetToEncodedHTMLCString( temp, pData->description );
              Concat(message, temp);
              ConcatCString(message, " (");
              SetToEncodedHTMLCString( temp, pData->defval );
              Concat(message, temp);
              ConcatCString(message, ")</td></tr>\n");
            }
          }
          free (value);
        }
        HTTP_SendString(request, message);

      } 

      if (nsteerable>0)
      {
        HTTP_Send(request, "</table>\n");
      }
      else
      {
        HTTP_Send(request, "<p>This thorn has no steerable parameters.</p>\n");
      }

      if(!readonly && nsteerable>0)
      {
        HTTP_Send(request,
               "<input type=\"submit\" value=\"Update all parameters\" />\n"
               "</form>\n");
      }


      /* Walk through non-all steerable parameters of given implementation. */
      first = 1;
      
      while(CCTK_ParameterWalk(first, thorn, NULL, &pData) == 0)
      {
        char *value;
        first = 0;
        Truncate(message,0);
        
        value = CCTK_ParameterValString (pData->name, pData->thorn);
        if(value)
        {          
          if (!(pData->steerable == CCTK_STEERABLE_ALWAYS))
          {
            
            if (nfixed == 0)
            {
              ConcatCString(message,"<h2>Fixed Parameters</h2>\n"
               "<table class=\"paramtable\" cellpadding=\"5\" cellspacing=\"5\">");
            }
            nfixed++;

            /* FIXME:  This is a hack - should put in parameter tags. */
            if(strcmp(thorn,CCTK_THORNSTRING) ||
               (strcmp(pData->name,"user")     && 
                strcmp(pData->name,"password") &&
                strcmp(pData->name,"encryption_scheme")))
            {
        
                ConcatCString(message, "<tr>\n<td><a href=\"/Parameters/");
                ConcatCString(message, pData->thorn);
                ConcatCString(message, "/");
                ConcatCString(message, pData->name);
                ConcatCString(message, "\">");
                ConcatCString(message, pData->name);
                ConcatCString(message, "</a></td>\n<td>");
                ConcatCString(message, value);
                ConcatCString(message, "</td>\n<td class=\"description\">");
                SetToEncodedHTMLCString( temp, pData->description );
                Concat(message, temp);
                ConcatCString(message, " (");
                SetToEncodedHTMLCString( temp, pData->defval );
                Concat(message, temp);
                ConcatCString(message, ")</td>\n</tr>\n");
            }
          }
          free(value);
        }
        HTTP_SendString(request, message);

      } 

      if (nfixed>0)
      {
        HTTP_Send(request, "</table>\n");
      }
      else
      {
        HTTP_Send(request, "<p>This thorn has no fixed parameters.</p>\n");
      }

      HTTP_Send(request, "</div>\n");
    }

    /* Write out the footer part. */
    
    HTTP_SetContentFooterString(cctkGH, 0, message);
    retval = HTTP_SendString(request, message);
  } /* n_arguments > 0 */
  String_Delete( message );
  String_Delete( menu );
  String_Delete( temp );
  return retval;
}

 /*@@
   @routine    ParameterPage
   @date       Thu Sep 21 15:13:55 2000
   @author     Gabrielle Allen
   @desc
   Page for an individual parameter
   Lots copied from ThornParameterPage
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
static int ParameterPage(const cGH *cctkGH, httpRequest *request, void *data)
{
  int retval=0;
  String * message = String_New();
  String * menu = String_New();
  String * temp = String_New();
  int first = 1;
  const cParamData *pData = (cParamData *)data;
  const cParamData *pDataWalk=NULL;
  char *value;
  char user[USER_LENGTH+1] = EMPTYSTRING;
  int notauthorised = HTTP_AuthenticateBasic(request, "user", user,
                  USER_LENGTH);
  int readonly = notauthorised;
  
  if(HTTP_NumArguments( request ) > 0)
  {
    /* This is a parameter set request */

    if(!notauthorised)
    {
      if(!readonly)
      {
        const httpArg *argument;
        /* Queue parameters for steering */
        first = 1;
        while((argument = HTTP_ArgumentWalk(request, first)) != NULL)
        {
          first = 0;
          fprintf(stderr, "Setting %s::%s to %s\n", pData->thorn,
                          HTTP_ArgName( argument),
                            HTTP_ArgValue( argument ));
          HTTP_SteerQueue(pData->thorn, HTTP_ArgName( argument),
                            HTTP_ArgValue( argument ));
        }
      }
      SendHTTP_Redirect_Header( request );
      ConcatCString(message, "Location: /Parameters/");
      ConcatCString(message, pData->thorn);
      ConcatCString(message, "/\r\n\r\n");

      HTTP_SendString(request, message);
    }
    else
    {
      SendHTTP_Uauthorized_Header( request );
  
      HTTP_Send(request, notauthorized_page);
    }
  } 
  else
  {    
    t_range *range;
    /* Display the page. */
    HTTP_SendOKHeader( request );

    HTTP_SetDoctype( message );
    HTTP_SendString(request, message);
    /* Start the page */
    HTTP_Send(request, "<html><head>\n");

    SetToCString(message, "<title>Cactus Parameter Request : ");
    ConcatCString(message, pData->name);
    ConcatCString(message, "</title>\n");

    HTTP_SendString(request, message);

    HTTP_SetHeadInfo( message);
    HTTP_SendString(request, message );
    HTTP_Send(request, "<style type=\"text/css\">\n"
                       ".paramtable td.description { font-size: small; } \n"
                       "</style>\n");

    HTTP_Send(request,"</head>\n<body>\n");

    Truncate( message, 0 );

    /* Menu for this page */
    first = 1;
    while(CCTK_ParameterWalk(first, pData->thorn, NULL, &pDataWalk) == 0)
    {
      if (first==1)
      {
        ConcatCString(menu,"<h3>");
        ConcatCString(menu, pDataWalk->thorn);
        ConcatCString(menu, ":</h3>\n");
      }
      first = 0;
      ConcatCString(menu, " <a href=\"/Parameters/");
      ConcatCString(menu, pDataWalk->thorn);
      ConcatCString(menu, "/");
      ConcatCString(menu, pDataWalk->name);
      ConcatCString(menu, "\">");
      ConcatCString(menu, pDataWalk->name);
      ConcatCString(menu, "</a><br />\n");
    }
    
    HTTP_SetContentHeaderString(cctkGH,0,message,menu);
    HTTP_SendString(request, message);

    SetToCString(message,"<h1>");
    ConcatCString(message, pData->thorn);
    ConcatCString(message, ": ");
    ConcatCString(message, pData->name);
    ConcatCString(message, "</h1>\n");
    HTTP_SendString(request, message);

    SetToCString(message,
            "<div class=\"centered\">Return to all parameters for this \n"
            "<a href=\"Parameters/");
    ConcatCString(message,
            CCTK_ThornImplementation(pData->thorn));
    ConcatCString(message,
            "\">thorn</a>.</div> ");

    HTTP_SendString(request, message);
            
    value = CCTK_ParameterValString (pData->name, pData->thorn);

    /* FIXME: HACK need a tag on parameter definition */
    if(strcmp(pData->thorn,CCTK_THORNSTRING) ||
       (strcmp(pData->name,"user")     && 
        strcmp(pData->name,"password") &&
        strcmp(pData->name,"encryption_scheme")))
    {

    SetToCString(message,"<div class=\"centered\">\n"
            "<table class=\"thornparams\" cellpadding=\"5\" cellspacing=\"5\" "
            " border=\"1\">\n" );
    ConcatCString(message, "<tr>\n<th>Name:</th>\n<td>");
    ConcatCString(message, pData->name);
    ConcatCString(message, "</td>\n</tr>\n" );
    ConcatCString(message, "<tr>\n<th>Thorn:</th>\n<td>");
    ConcatCString(message, pData->thorn);
    ConcatCString(message, "</td>\n</tr>\n" );
    ConcatCString(message, "<tr>\n<th>Implementation:</th>\n<td>");
    ConcatCString(message, CCTK_ThornImplementation(pData->thorn));
    ConcatCString(message, "</td>\n</tr>\n" );
    ConcatCString(message, "<tr>\n<th>Current value:</th>\n<td>");
    ConcatCString(message, value);
    ConcatCString(message, "</td>\n</tr>\n" );
    ConcatCString(message, "<tr>\n<th>Description::</th>\n<td>");
    SetToEncodedHTMLCString( temp, pData->description );
    Concat(message, temp);
    ConcatCString(message, "</td>\n</tr>\n" );
    ConcatCString(message, "<tr>\n<th>Default:::</th>\n<td>");
    SetToEncodedHTMLCString( temp, pData->defval );
    Concat(message, temp);
    ConcatCString(message, "</td>\n</tr>");
    HTTP_SendString(request, message);
    SetToCString(message,"<tr>\n<th>Steerable:</th>\n<td>");

    switch(pData->steerable)
    {
      case CCTK_STEERABLE_ALWAYS : 
        ConcatCString(message,"Always");
        break;
      case CCTK_STEERABLE_NEVER : 
        ConcatCString(message,"Never");
        break;
      case CCTK_STEERABLE_RECOVER : 
        ConcatCString(message,"Recovery");
        break;
      default  :
        ConcatCString(message,"Not matched");
    }
    ConcatCString(message,"</td>\n</tr>\n");
    HTTP_SendString(request, message);

    SetToCString(message,"<tr>\n"
           "<th>Type:</th>\n<td>");
    switch(pData->type)
    {
      case PARAMETER_BOOLEAN : 
        ConcatCString(message,"Boolean");
        break;
      case PARAMETER_REAL : 
        ConcatCString(message,"Real");
        break;
      case PARAMETER_INTEGER : 
        ConcatCString(message,"Integer");
        break;
      case PARAMETER_SENTENCE : 
        ConcatCString(message,"Sentence");
        break;
      case PARAMETER_STRING : 
        ConcatCString(message,"String");
        break;
      case PARAMETER_KEYWORD : 
        ConcatCString(message,"Keyword");
        break;
      default :
        ConcatCString(message,"Not matched");
    }
    ConcatCString(message,"</td>\n</tr>\n");

    HTTP_SendString(request, message);

    
    SetToCString(message,"<tr><th>Scope:</th><td>");
    switch(pData->scope)
    {
      case SCOPE_GLOBAL : 
        ConcatCString(message,"Global</td></tr>\n");
        break;
      case SCOPE_RESTRICTED : 
        ConcatCString(message,"Restricted</td></tr>\n");
        break;
      case SCOPE_PRIVATE : 
        ConcatCString(message,"Private</td></tr>\n");
        break;
      default :
        ConcatCString(message,"Not matched</td></tr>\n");
    }

    HTTP_SendString(request, message);

    Truncate(message,0);
    first = 1;
    for(range = pData->range; range ; range = range->next)
    {
      if (first==1)
      {
        ConcatCString(message,
                "<tr>\n"
                "<th>\n"
                "Range:"
                "</th>\n<td><dl>\n"
                );
      }
      first = 0;
      ConcatCString(message, "<dt>");
      ConcatCString(message, range->range);
      ConcatCString(message, "</dt>\n<dd>");
      SetToEncodedHTMLCString( temp, range->description );
      Concat(message, temp);
      ConcatCString(message, "</dd>\n");
      if (!CCTK_Equals(range->origin,pData->thorn))
      {
        ConcatCString(message, "<br />[Extended by thorn ");
        ConcatCString(message, range->origin);
        ConcatCString(message, "]\n");
      }
    }                
    if(first==0)
    {
      ConcatCString(message, "</dl></td>\n</tr>\n");
    }

    HTTP_SendString(request, message);
      
    SetToCString(message, "<tr>\n<th>Times Set:</th>\n<td>");
    ConcatDecimal(message, pData->n_set);
    ConcatCString(message, "</td>\n</tr>\n</table>\n</div>\n");
    HTTP_SendString(request, message);

    }
    else
    {
      HTTP_Send(request,"<p>Hidden parameter, information is not available</p>\n");
    }


    /* Write out the footer part. */

    HTTP_SetContentFooterString(cctkGH, 0, message);
    retval = HTTP_SendString(request, message);
  }
  String_Delete( message );
  String_Delete( menu );
  String_Delete( temp );
  return retval;
}
