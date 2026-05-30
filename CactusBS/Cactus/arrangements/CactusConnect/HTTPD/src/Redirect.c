 /*@@
   @file      Redirect.c
   @date      Fri May 18 10:03:15 2001
   @author    Tom Goodale
   @desc
              File to contain John Shalf's HTTP redirect stuff.
   @enddesc
   @version   $Id$
 @@*/

#include <stdlib.h>
#include <string.h>

#include "cctk.h"

#ifdef CCTK_MPI
#include <mpi.h>
#endif /* CCTK_MPI */

#include "util_Network.h"    /* Util_GetHostName() */

#include "httpd.h"

#include "httpRequest.h"
#include "Content.h"

#include "Redirect.h"
#include "SString_Namespace.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusConnect_HTTPD_Redirect_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

static void RegisterRedirect(void);
static int RedirectPage(const cGH *cctkGH, httpRequest *request, void *data);

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

static char httpredirect=0;
static char httpmaster[1024];

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    HTTP_SetupRedirect
   @date       Wed April 20 20:39:15 2001
   @author     John Shalf
   @desc
   Creates a socket to listen on purely for server redirect to the root node.
   @enddesc
   @calls
   @calledby
   @history

   @endhistory

@@*/

#define ADDRLEN 1024

int HTTP_SetupRedirect(int port,
                       int queue_size,
                       int hunt)
{
  int retval = 0;
  /* should not hunt */
  int i;
  char hostnm[ADDRLEN] = EMPTYSTRING;
  char *alladdr;

  int nprocs   = CCTK_nProcs(NULL);
  int proc     = CCTK_MyProc(NULL);

#ifdef HTTPD_DEBUG
  CCTK_INFO("enter setup redirect------------");
#endif
  memset(hostnm,0,sizeof(hostnm));
  Util_GetHostName(hostnm, sizeof(hostnm));
  alladdr = (char *) calloc (sizeof(hostnm), nprocs);

#ifdef CCTK_MPI
  { /* push stack */
#ifdef HTTPD_DEBUG
    CCTK_INFO("all gather");
#endif

    MPI_Allgather(hostnm,sizeof(hostnm),MPI_CHAR,
                  alladdr,sizeof(hostnm),MPI_CHAR,
                  MPI_COMM_WORLD);
#ifdef HTTPD_DEBUG
    CCTK_INFO("collected");
#endif
  }
#endif

  /* set our master */
#ifdef HTTPD_DEBUG
  printf("alladdr is %s\n", alladdr);
#endif
  strncpy(httpmaster,alladdr, sizeof(httpmaster)); /* the 0th element is the master hostname */
  httpmaster[sizeof(httpmaster) - 1] = '\0';

  /* so compare my addr to the list sequentially */
  for(i=0;i<nprocs;i++)
  {
#ifdef HTTPD_DEBUG
    printf("Cycle through %u\n",i);
#endif
    if(!strcmp(alladdr+i*sizeof(hostnm),hostnm))
    {
      /* we matched addresses */
      if(i<proc || i==0)
      {
        break; /* I'm not the lowest order */
      }
      else
      {
      /* otherwise, I am the lowest ranked processor with this
         address that isn't the primary webserver.  So
         I should do a server redirect */
        httpredirect=1; /* this has a redirect socket */
        HTTP_SetupServer(port,queue_size,hunt);

        CCTK_VInfo(CCTK_THORNSTRING, "HTTP requests will be redirected to "
                   "primary webserver node on\n                http://%s:%d/",
                   HTTP_Master(), (unsigned int) HTTP_Port());

        RegisterRedirect();
        retval = 1;
        break; /* one redirect is enough */
      }
    }
  }
  if(! retval)
  {
#ifdef HTTPD_DEBUG
    CCTK_INFO("No Redirect**************\n");
#endif
    httpredirect=0;
  }

  free (alladdr);

  return retval;
}

 /*@@
   @routine    HTTP_Master
   @date       Friday April 20 11:24:40 2001
   @author     John Shalf
   @desc
   Reports the hostname of the node the HTTP server is listening on.
   The returns the hostname of the master as a string.
   @enddesc
   @calls
   @calledby
   @history

   @endhistory
@@*/
const char *HTTP_Master(void)
{
  return httpmaster;
}

 /*@@
   @routine    HTTP_IsServer
   @date       Friday April 20 11:24:40 2001
   @author     John Shalf
   @desc
           True if this node has a server socket open.
   @enddesc
   @calls
   @calledby
   @history

   @endhistory
@@*/
int HTTP_IsServer(void)
{
  return (CCTK_MyProc (NULL) == 0 || httpredirect);
}

 /*@@
   @routine    HTTP_IsRedirect
   @date       Friday April 20 11:24:40 2001
   @author     John Shalf
   @desc
           True only if this node has a server open and that server
        is exclusively used to redirect web traffic back to the root
    node.
   @enddesc
   @calls
   @calledby
   @history

   @endhistory
@@*/
int HTTP_IsRedirect(void)
{
  return httpredirect;
}


#if 0 /* this would be more efficient, but PUGH doesn't define
         a datatype for 8-byte integers.  So instead,
         setupRedirect will use the actual hostname. */
CCTK_INT8 Util_GetHostAddr(char *hostname)
{
  hostent *he;
  char hostname[1025];
  int addrlen;
  CCTK_INT8 addr;

  Util_GetHostName(hostname, 1024);
  he=gethostbyname(hostname);
  addrlen=he->h_length;
  if(addrlen>sizeof(CCTK_INT8)) addrlen=CCTK_INT8;
  addr=0;
  memcpy(&addr,he->h_addr,addrlen);
  return addr;
}
#endif



/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/

 /*@@
   @routine    RegisterRedirect
   @date       Fri April 20 23:47:43 2001
   @author     John Shalf
   @desc
   Redirection Page Registration.
   @enddesc
   @calls
   @calledby
   @history

   @endhistory

@@*/
static void RegisterRedirect(void)
{
#ifdef HTTP_DEBUG
    printf("Registering Redirection Page");
#endif
  HTTP_RegisterPage("/index.html",RedirectPage,NULL);
}

 /*@@
   @routine    RedirectPage
   @date       Fri April 20 23:47:43 2001
   @author     John Shalf
   @desc
   Redirects the webbrowser to the main server page
   @enddesc
   @calls
   @calledby
   @history
   @hdate Thu Sep 14 10:54:22 2000 @hauthor John Shalf
   @hdesc For clusters where its difficult to know which node is going to
   have the webserver process, this opens a port on *all* distinct nodes
   and redirects all web connection attempts to the root node for processing.
   @endhistory

@@*/

static int RedirectPage(const cGH *cctkGH, httpRequest *request, void *data)
{
  int retval = -1;
  String *message = String_New();

  /* avoid compiler warning about unused parameter */
  data = data;

  HTTP_SendOKHeader(request);

  /* Start the page */
  SetToCString(message,"<html><head><title>Server Redirect</title>\n");
  /* Write out the main header part */
  HTTP_SendString(request, message);
  SetToCString(message,"<meta http-equiv=\"refresh\" content=\"1; URL=http://");
  ConcatCString(message,HTTP_Master());
  ConcatCString(message, ":");
  ConcatDecimal(message,(unsigned int) HTTP_Port());
  ConcatCString(message, "\" />\n</head>\n<body>\n");
  HTTP_SendString(request,message);
  /* ********** Server Redirect To Master ************* */
  SetToCString(message,"<h1>Redirect to master host=");
  ConcatCString(message,HTTP_Master());
  ConcatCString(message, ":");
  ConcatDecimal(message,(unsigned int) HTTP_Port());
  ConcatCString(message, "</h1>\n");
  HTTP_SendString(request,message);

  HTTP_SetContentFooterString(cctkGH, 0, message);
  retval = HTTP_SendString(request, message);

  String_Delete( message );
  return retval;
}
