 /*@@
   @file      Startup.c
   @date      Wed Sep 13 21:26:56 2000
   @author    Tom Goodale
   @desc 
   Scheduled routines for the HTTPD thorn.
   @enddesc 
   @version $Header$
 @@*/

#include "cctk.h"
#include <stdio.h>
#include <stdlib.h>

#ifdef HAVE_CAPABILITY_PTHREADS
#include <pthread.h>
#endif

#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "httpd.h"
#include "Steer.h"
#include "Redirect.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusConnect_HTTPD_Startup_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

#ifdef HAVE_CAPABILITY_PTHREADS
static void * HTTP_Thread(void *cctkGH);
static void HTTP_SetupPollingThread(cGH *cctkGH);
#endif

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

int HTTP_StartServer(void);
int HTTP_FirstServ(void);
void HTTP_Work(CCTK_ARGUMENTS);
int HTTP_Shutdown(void);

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

#ifdef HAVE_CAPABILITY_PTHREADS
static pthread_t polling_thread;
static int thread_started = 0; 
#endif

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/



 /*@@
   @routine    HTTP_StartServer
   @date       Wed Sep 13 21:26:56 2000
   @author     Tom Goodale
   @desc 
   Starts the webserver.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

int HTTP_StartServer(void)
{
  int use_port;
  char *environ_port;
  DECLARE_CCTK_PARAMETERS


#ifndef HAVE_CAPABILITY_PTHREADS
  if (use_pthreads)
  {
    CCTK_WARN (1, "Parameter 'HTTPD::use_pthreads' is set to \"yes\" but you "
                  "didn't compile with PTHREADS. Setting will be ignored.");
  }
#endif

  /* get HTTPD port to use from the 'port' parameter or the shell environment */
  use_port = port;
  environ_port = getenv ("HTTPD_PORT");
  if (environ_port)
  {
    use_port = atoi (environ_port);
    if (use_port >= 1 && use_port <= 65535)
    {
      CCTK_VInfo (CCTK_THORNSTRING, "Environment variable HTTPD_PORT is set to "
                  "%d and will override parameter setting from 'HTTPD::port'",
                  use_port);
    }
    else
    {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Environment variable HTTPD_PORT is set to '%s' which cannot "
                  "be parsed as a valid port number. Using parameter setting "
                  "'HTTPD::port = %d' instead.", environ_port, (int) port);
      use_port = port;
    }
  }

  if(CCTK_MyProc(NULL) == 0)
  {
    /* Does the server provide any pages by default ? */
    if(provide_pages)
    {
      HTTP_RegisterPages();
    }
    
    HTTP_SetupServer(use_port, queue_length, hunt);
  }

  /* Do we need to redirect to a master server ? */
  HTTP_SetupRedirect(use_port,queue_length,hunt);

  return 0;
}

 /*@@
   @routine    HTTP_FirstServ
   @date       Tue Nov 12 20:32:36 2002
   @author     Tom Goodale
   @desc 
   Listen for connection requests.
   @enddesc 
   @calls     
   @calledby   
   @history 
   @hdate Tue Nov 12 20:33:07 2002
   @hauthor Tom Goodale
   @hdesc This was originally in HTTP_StartServer but
          that made it impossible to spawn a task telling
          it where the webserver was if httpd::pause was "true".
   @endhistory 
 
 @@*/
int HTTP_FirstServ(void)
{

  httpState state;

  HTTP_UpdateState(NULL, &state);

  do
  {
    if(HTTP_IsServer())
    {
      HTTP_Poll(NULL, state.timeout_seconds, state.timeout_useconds);
    }
    
    if(state.steer == 1)
    {
      HTTP_SteerDispatch();

      HTTP_UpdateState(NULL, &state);
    }

    if(state.terminate)
    {
      HTTP_Terminate(NULL);
    }

  } while(state.paused);

  return 0;
}

 /*@@
   @routine    HTTP_Work
   @date       Wed Sep 13 21:26:56 2000
   @author     Tom Goodale
   @desc 
   Main working routine for the webserver.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void HTTP_Work(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  httpState state;

#ifdef HAVE_CAPABILITY_PTHREADS
  if(HTTP_IsServer() && use_pthreads && ! thread_started)
  {
    HTTP_SetupPollingThread(cctkGH);
  }
#endif

  HTTP_UpdateState(cctkGH, &state);
  
  do
  {
#ifdef HAVE_CAPABILITY_PTHREADS
    if(!use_pthreads)
    {
#endif
      if(HTTP_IsServer())
      {
        HTTP_Poll(cctkGH, state.timeout_seconds, state.timeout_useconds);
      }
#ifdef HAVE_CAPABILITY_PTHREADS
    }
#endif

    if(state.steer == 1)
    {
      HTTP_SteerDispatch();

      HTTP_UpdateState(cctkGH, &state);
    }

    if(state.terminate)
    {
      HTTP_Terminate(cctkGH);
    }

  } while(state.paused);

}

 /*@@
   @routine    HTTP_Shutdown
   @date       Wed Sep 13 21:26:56 2000
   @author     Tom Goodale
   @desc 
   Shutdown routine for the webserver.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
int HTTP_Shutdown(void)
{
  DECLARE_CCTK_PARAMETERS
#ifdef HAVE_CAPABILITY_PTHREADS
  /* Wait for the polling thread to exit */
  if(HTTP_IsServer() && use_pthreads && thread_started)
  {
    CCTK_ParameterSet("terminate", CCTK_THORNSTRING, "yes");
    pthread_join(polling_thread, NULL);
  }
#endif

  if(HTTP_IsServer())
  {
    HTTP_ShutdownServer();
  }

  return 0;
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/

 /*@@
   @routine    HTTP_Thread
   @date       Wed Oct 25 17:04:59 2000
   @author     Tom Goodale
   @desc 
   The http thread routine.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
   @var     cctkGH
   @vdesc   The cGH
   @vtype   cGH *
   @vio     inout
   @vcomment 
 
   @endvar 

@@*/
#ifdef HAVE_CAPABILITY_PTHREADS
static void * HTTP_Thread(void *cctkGH)
{
  CCTK_INT terminate;
  int type;

  do
  {
    /* Only poll for a few seconds in case we are terminated */
    HTTP_Poll((cGH *)cctkGH, 10, 0);

    terminate = *(const CCTK_INT *)CCTK_ParameterGet("terminate",CCTK_THORNSTRING, &type);

  } while(! terminate);

  return NULL;
}

 /*@@
   @routine    HTTP_SetupPollingThread
   @date       Wed Oct 25 17:05:38 2000
   @author     Tom Goodale
   @desc 
   Sets up the http polling thread.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
   @var     cctkGH
   @vdesc   The cGH
   @vtype   cGH *
   @vio     inout
   @vcomment 
 
   @endvar 

@@*/
static void HTTP_SetupPollingThread(cGH *cctkGH)
{
  if(pthread_create(&polling_thread, NULL, HTTP_Thread, (void *)cctkGH))
  { 
    perror("pthread_create: ");
    CCTK_Exit(cctkGH,99);
  }
  else
  { 
    thread_started = 1;
  }
}
#endif /* HAVE_CAPABILITY_PTHREADS */
