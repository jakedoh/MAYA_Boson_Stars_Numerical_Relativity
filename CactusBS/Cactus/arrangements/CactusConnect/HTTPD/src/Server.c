 /*@@
   @file      Server.c
   @date      Wed Sep 13 20:10:24 2000
   @author    Tom Goodale
   @desc 
   Web serving and utility routines.   
   @enddesc
   @version  $Header$
 @@*/

#include <stdlib.h>
#include <string.h>

#include "cctk.h"

#include "cctk_Parameters.h"

#include "httpd_Map.h"
#include "util_String.h"

#include "httpd.h"

#include "Steer.h"
#include "Expression.h"

#include "SString_Namespace.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusConnect_HTTPD_Server_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

typedef struct
{
  int (*function)(const cGH *, httpRequest *, void *);
  void *data;
} httpPage;

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

static httpPage *CreatePageData(int (*function)(const cGH *,httpRequest *, void *), void *data);
static httpPage *FindPage(const char *path, const char **residual);

static int StatusUntilIt (const cGH *cctkGH);
static int StatusUntilTime (const cGH *cctkGH);
static int StatusUntilExpression (cGH *cctkGH);

static double evaluator(const char *name, void *data);


/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

static uMap pages = NULL;

static const char *notfound_page =
"<html>\n<head>\n<title>Error 404: Not Found</title>\n</head>\n"
"<body>\nThe URI you requested could not be found\n</body>\n</html>\n";

static const char *notimplemented_page =
"<html>\n<head>\n<title>Error 501: Not Implemented</title>\n</head>\n"
"<body>\nThe requested method is not implemented\n</body>\n<html>\n";

#define INITIAL_SIZE 16

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    HTTP_RequestGet
   @date       Wed Sep 13 20:10:24 2000
   @author     Tom Goodale
   @desc 
   Routine to deal with a GET request.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

int HTTP_RequestGET(cGH *cctkGH, httpRequest *request)
{
  int retval = -1;
  const char *residual = NULL;
  httpPage *pagedata = FindPage(HTTP_URI( request ), &residual );
  HTTP_SetResidual( request, residual );

  if( pagedata )
  {
    retval = pagedata->function(cctkGH, request, pagedata->data);
  }

  if(retval < 0)
  {
    HTTP_Send(request,"HTTP/1.0 404 Not Found\r\n");
    
    HTTP_Send(request,"Content-Type: text/html\r\n\r\n");

    HTTP_Send(request, notfound_page);
  }

  return retval;
}

 /*@@
   @routine    HTTP_RequestUnsupported
   @date       Thu Sep 14 12:18:56 2000
   @author     Tom Goodale
   @desc 
   
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
int HTTP_RequestUnsupported(cGH *cctkGH, httpRequest *request)
{
  cctkGH = cctkGH; /* avoid compiler warning about unused parameter */

  HTTP_Send(request,"HTTP/1.0 501 Not Implemented\r\n");

  HTTP_Send(request,"Content-Type: text/html\r\n\r\n");
  
  HTTP_Send(request, notimplemented_page );

  return 0;
}

 /*@@
   @routine    HTTP_RegisterPage
   @date       Wed Sep 13 20:10:24 2000
   @author     Tom Goodale
   @desc 
   Routine to register a web page.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
int HTTP_RegisterPage(const char *path, int (*function)(const cGH *, httpRequest *, void *), void *data)
{
  int retval = -1;

  /* Create the hash table if it's not already been created */
  if(! pages)
  {
    pages = Httpd_MapCreate();
  }

  if(Httpd_MapData(pages, strlen(path), path))
  {
    CCTK_VWarn(1, __LINE__,__FILE__,CCTK_THORNSTRING,
               "Page exists already:\n\"%s\"", path);
  }
  else
  {
    httpPage *pagedata = CreatePageData(function, data);

    if(pagedata)
    {
      retval = Httpd_MapStore(pages, strlen(path), path, (void *)pagedata);
    }
    else
    {
      retval = -2;
    }
  }

  return retval;
}


 /*@@
   @routine    HTTP_UpdateState
   @date       Sat Sep 16 21:12:23 2000
   @author     Tom Goodale
   @desc 
   Updates the state of the server from the latest parameter values.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
int HTTP_UpdateState(cGH *cctkGH, httpState *state)
{
  DECLARE_CCTK_PARAMETERS
  int stepping = single_step;
  static int steering = 1; 
  static int last_steered = -1;


  state->paused = (! stepping) && 
                  (pause || (cctkGH &&
                  ((until_it_active         && StatusUntilIt(cctkGH)) ||
                  (until_time_active   && StatusUntilTime(cctkGH)) ||
                  (until_expression_active && StatusUntilExpression(cctkGH)))));

  /* If single step is selected, need to reset the parameter to force a pause
   * next time.
   */
  if(stepping)
  {
    HTTP_SteerQueue(CCTK_THORNSTRING, "single_step", "no");
  }

  if(cctkGH && ! state->paused)
  {
    if(steering)
    {
      if(steering_frequency > 0 &&
         steering_frequency+last_steered <= cctkGH->cctk_iteration)
      {
        /* Time to steer */
        state->steer = 1;
        last_steered = cctkGH->cctk_iteration;
      }
      else
      {
        state->steer = 0;
      }
      
      if(steering_frequency <= 0)
      {
        /* Once steering is switched off, can't restart it ! */
        steering = 0;
      }
    }
  }
  else
  {
    state->steer = 1;
  }

  state->terminate = terminate;
  if(!state->paused)
  {  
    state->timeout_seconds = timeout_seconds;
    state->timeout_useconds = timeout_useconds;
  }
  else
  {
    state->timeout_seconds = -1;
    state->timeout_useconds = -1;
  }

  return 0;
}

 /*@@
   @routine    HTTP_Terminate
   @date       Sat Sep 16 21:40:27 2000
   @author     Tom Goodale
   @desc 
   Terminate the simulation.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
int HTTP_Terminate(cGH *cctkGH)
{
  HTTP_SteerQueue("Cactus","terminate_next","yes");

  return 0;
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/

 /*@@
   @routine    CreatePageData
   @date       Wed Sep 13 20:10:24 2000
   @author     Tom Goodale
   @desc 
   Creates a httpPage structure.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
static httpPage *CreatePageData(int (*function)(const cGH *, httpRequest *, void *), void *data)
{
  httpPage *pagedata = (httpPage *)malloc(sizeof(httpPage));

  if(pagedata)
  {
    pagedata->function = function;
    pagedata->data = data;
  }

  return pagedata;
}

 /*@@
   @routine    FindPage
   @date       Wed Sep 13 20:10:24 2000
   @author     Tom Goodale
   @desc 
   Finds a page, if it exists.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
static httpPage *FindPage(const char *path, const char **residual)
{
  httpPage *pagedata = NULL;

#ifdef HTTP_DEBUG
  printf ("Searching for '%s'\n", path);
#endif

  if(pages && path)
  {
    String *temp = String_New();
    /* Check for index.html */
    if(strlen(path) == 0 || strcmp( path, "/") == 0 )
    {
#ifdef HTTP_DEBUG
      printf("Looking for '%sindex.html'\n", path);
#endif
      SetToCString( temp, path);
      ConcatCString( temp,"index.html");

      pagedata = Httpd_MapData(pages, Length(temp), GetBuffer(temp));
    }
    else if((pagedata = Httpd_MapData(pages, strlen(path), path)))
    {
      /* Or exact path */
    }
    else
    {
#ifdef HTTP_DEBUG
      printf("Looking for '%s/index.html'\n", path);
#endif
      SetToCString( temp, path);
      ConcatCString( temp,"/index.html");

      pagedata = Httpd_MapData(pages, Length(temp), GetBuffer(temp));

    }
    *residual = NULL;

    if(!pagedata && strlen( path ) > 0)
    {
      const char *position;
      /* Ok, now cycle through.  Know it doesn't end with a slash */
      for(position = path+strlen(path)-1; position >= path; position--)
      {
        if(*position == '/')
        {
#ifdef HTTP_DEBUG
          printf("Looking for '%s' less '%s' \n", path, position);
#endif
          if((pagedata = Httpd_MapData(pages, position-path, path)))
          {
            *residual = position+1;
            break;
          }
        }
      }
    }
    String_Delete( temp );
    temp = NULL;
  }

  return pagedata;
}

 /*@@
   @routine    StatusUntilIt
   @date       Tue Sep 19 22:59:34 2000
   @author     Tom Goodale
   @desc 
   Is the iteration greater than the desired stop value ?
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
static int StatusUntilIt (const cGH *cctkGH)
{
  DECLARE_CCTK_PARAMETERS


  return (cctkGH && until_it <= cctkGH->cctk_iteration);
}


 /*@@
   @routine    StatusUntilTime
   @date       Tue Sep 19 23:00:16 2000
   @author     Tom Goodale
   @desc 
   Is the simulation time greater than the stop value ?
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
static int StatusUntilTime (const cGH *cctkGH)
{
  DECLARE_CCTK_PARAMETERS


  return (cctkGH && until_time <= cctkGH->cctk_time);
}


 /*@@
   @routine    StatusUntilExpression
   @date       Tue Sep 19 23:00:47 2000
   @author     Tom Goodale
   @desc 
   Is the stop expression true ?
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
static int StatusUntilExpression (cGH *cctkGH)
{
  DECLARE_CCTK_PARAMETERS
  static char *parsed_expression = NULL;
  static int times_set = -1;
  int retval = 0;

  /* See if we need to parse the expression again. */
  int new_times_set = CCTK_ParameterQueryTimesSet("until_expression",
                                              CCTK_THORNSTRING);
  
  if(new_times_set > times_set)
  {
    times_set = new_times_set;
    if(parsed_expression)
    {
      free(parsed_expression);
    }
    parsed_expression = HTTP_ExpressionParse(until_expression);
  }

  if(parsed_expression && strlen(parsed_expression) > 0 && cctkGH)
  {
    /* Make a copy */
    char *copy = Util_Strdup(parsed_expression);
    
    /* Evaluate the expression */
    retval = HTTP_ExpressionEvaluate(copy, evaluator, cctkGH);

    /* Free the copy */
    free(copy);
  }

  return retval;
}


 /*@@
   @routine    evaluator
   @date       Tue Sep 19 23:07:20 2000
   @author     Tom Goodale
   @desc 
   Takes the name of a gridscalar and returns its value.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
static double evaluator(const char *expression, void *data)
{
  double retval = 0.0;
  cGH *cctkGH = (cGH *)data;
  int varindex = CCTK_VarIndex(expression);

  if(varindex > -1)
  {
    int vartype  = CCTK_VarTypeI(varindex);
    int *pointer = CCTK_VarDataPtrI(cctkGH, 0, varindex);
    
    switch(vartype)
    {
      case CCTK_VARIABLE_BYTE :
        retval = *((CCTK_BYTE *)pointer);
        break;
      case CCTK_VARIABLE_INT :
        retval = *((CCTK_INT *)pointer);
        break;
      case CCTK_VARIABLE_REAL :
        retval = *((CCTK_REAL *)pointer);
        break;
      default :
        CCTK_WARN( 0, "Unsupported variable type" );
        retval = 0.0;
    }
  }
  else if(CCTK_Equals(expression,"time"))
  {
    retval = cctkGH->cctk_time;
  }
  else if(CCTK_Equals(expression,"iteration"))
  {
    retval = cctkGH->cctk_iteration;
  }    
  else
  {
    retval = strtod(expression,NULL);
  }

  return retval;
}
