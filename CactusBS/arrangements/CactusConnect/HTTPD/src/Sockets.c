 /*@@
   @file      Sockets.c
   @date      Wed Sep 13 20:39:15 2000
   @author    Tom Goodale
   @desc
   Routines which deal with sockets.
   These should probably move into thorn Socket at some point if
   they aren't already there.
   @enddesc
   @version $Header$
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "util_Network.h"
#include "SocketUtils.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_SIGNAL_H
#include <signal.h>
#endif

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif /* HAVE_UNISTD_H */
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif /* HAVE_SYS_TIME_H */
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif /* HAVE_SYS_TYPES_H */
#ifdef HAVE_SYS_SOCKET_H
#include <sys/socket.h>
#endif /* HAVE_SYS_SOCKET_H */
#ifdef HAVE_NETINET_IN_H
#include <netinet/in.h>
#endif /* HAVE_NETINET_IN_H */
#ifdef HAVE_NETDB_H
#include <netdb.h>
#endif /* HAVE_NETDB_H */
#ifdef HAVE_ARPA_INET_H
#include <arpa/inet.h>
#endif /* ARPA_INET_H */
#ifdef HAVE_WINSOCK2_H
#include <winsock2.h>
#endif /* HAVE_WINSOCK2_H */
#include <errno.h>

#include "httpd.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusConnect_HTTPD_Sockets_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/
#ifndef MSG_NOSIGNAL
#define MSG_NOSIGNAL 0
#endif

typedef enum {closed, open} httpSocketState;

typedef struct HTTPSocketTag
{
  httpSocket *prev;
  httpSocket *next;
  unsigned long filedes;
  httpSocketState state;
} httpSocket_PLACEHOLDER;

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

static httpSocket *SocketCreate(unsigned long int filedes);
static void SocketDestroy(httpSocket *this);
static void SocketClose(httpSocket *this);

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/* Active file descriptors */
static fd_set active_fd_set;

/* Main server socket */
static SOCKET sock;
static SOCKET minsock;
static SOCKET maxsock;

static httpSocket *socklist = NULL;

static unsigned long int httpport = 0;

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    HTTP_Send
   @date       10.04.2004
   @author     Steve White
   @desc
   Writes all of an HTTP reply.
   @enddesc
   @calls
   @calledby
   @history

   @endhistory

@@*/
int
HTTP_Send( httpRequest * request, const char * message )
{
  return HTTP_Write(request, message, strlen(message)); 
}

/*@@
   @routine    HTTP_SetupServer
   @date       Wed Sep 13 20:39:15 2000
   @author     Tom Goodale
   @desc
   Creates a socket to listen on.
   @enddesc
   @calls
   @calledby
   @history

   @endhistory

@@*/
int HTTP_SetupServer(int port, int queue_size, int hunt)
{
  char hostname[1025] = {'\0'};
  unsigned int realport = 0;

  /* Some systems need special logic for starting up TCP. */
  Socket_InitializeSocketLayer ();

  /* Create the socket and set it up to accept connections. */
  sock = Socket_TCPOpenServerSocket (port, hunt ? &realport : NULL, queue_size);
  if (sock == INVALID_SOCKET)
  {
    CCTK_WARN (0, "HTTPD Failed to create server socket");
  }

  Util_GetHostName(hostname, sizeof(hostname));

  httpport = hunt ? realport : (unsigned long int) port;

  printf("Server started on http://%s:%lu/\n", hostname, httpport);
  if (CCTK_IsFunctionAliased("Send_Twitter_Msg"))
  {
    char msg[141];
    char *this_user = getenv("USER");
    if (!this_user) this_user = "unknown";
    /* We do not want 'hostname' here because that can be wrong and not be overwritten */
    char *myhostname = getenv("HOSTNAME");
    if (!myhostname) myhostname = hostname;
    /* If there is an environment variable HOSTPORT, announce this instead (for ssh forwarding) */
    if (getenv("HOSTPORT"))
      snprintf(msg, 140, "HTTPD by %s on http://%s:%s/\n", this_user, myhostname, getenv("HOSTPORT"));
    else
      snprintf(msg, 140, "HTTPD by %s on http://%s:%lu/\n", this_user, myhostname, httpport);
    Send_Twitter_Msg(msg);
  }
  else printf("Not announcing location via Twitter.\n");

  minsock = sock;
  maxsock = sock;

  /* Initialize the set of active sockets. */
  FD_ZERO (&active_fd_set);
  FD_SET (sock, &active_fd_set);

  /* Not all architectures support the MSG_NOSIGNAL flag to send and recv,
   * so setup to ignore the SIGPIPE via normal signal hadnling as well.
   */

#ifdef SIGPIPE
  signal(SIGPIPE,SIG_IGN);
#endif

  return httpport; /* return the actual port here */
}

 /*@@
   @routine    HTTP_ShutdownServer
   @date       Wed Sep 13 20:39:15 2000
   @author     Tom Goodale
   @desc
   Closes all sockets we are interested in.
   @enddesc
   @calls
   @calledby
   @history

   @endhistory

@@*/
int HTTP_ShutdownServer(void)
{
  SOCKET i;

  /* Close all sockets in our active set */
  for(i = maxsock; i >= minsock; i--)
  {
    if(FD_ISSET(i, &active_fd_set))
    {
      Socket_CloseSocket (i);
    }
  }

  return 0;
}

 /*@@
   @routine    HTTP_Poll
   @date       Wed Sep 13 20:39:15 2000
   @author     Tom Goodale
   @desc
   Main workhorse routine.
   Looks for activity on any of the sockets which are open
   and dispatches work.
   @enddesc
   @calls
   @calledby
   @history

   @endhistory

@@*/
int HTTP_Poll(cGH *cctkGH, long sec, long usec)
{
  fd_set read_fd_set = active_fd_set;
  struct sockaddr_in clientname;
  struct timeval timeout;
  httpSocket *this;


  timeout.tv_sec = sec;
  timeout.tv_usec = usec;

  /* Check if any input is available on one or more active sockets. */

  if (select (FD_SETSIZE, &read_fd_set, NULL, NULL, sec >= 0 ? &timeout : NULL)
      == SOCKET_ERROR)
  {
    perror ("select");
    CCTK_Abort(cctkGH, EXIT_FAILURE);
  }

  /* Service all the sockets with input pending. */
  if (FD_ISSET (sock, &read_fd_set))
  {
    /* Connection request on original socket. */
#ifdef HAVE_SOCKLEN_T
    socklen_t size = sizeof (clientname);
#else
    int size = sizeof (clientname);
#endif
    SOCKET new = accept (sock, (struct sockaddr *) &clientname, &size);
    if (new == INVALID_SOCKET)
    {
      perror ("accept");
      /*CCTK_Abort(cctkGH, EXIT_FAILURE);*/
    }
    else
    {
      if (*(const CCTK_INT *) CCTK_ParameterGet ("verbose", CCTK_THORNSTRING,
                                                 NULL))
      {
        CCTK_VInfo (CCTK_THORNSTRING,
                    "Server: connect from host %s, port %hd",
                    inet_ntoa (clientname.sin_addr),
                    ntohs (clientname.sin_port));
      }
      SocketCreate(new);
    }
  }

  this = socklist;
  while( this )
  {
    httpSocket *next = this->next;
    /* Data arriving on an already-connected socket. */
    if (FD_ISSET (this->filedes, &read_fd_set) &&
        HTTP_ReadFromClient (cctkGH, this) < 0)
    {
      SocketDestroy(this);
    }
    this = next;
  }

  return (0);
}


 /*@@
   @routine    HTTP_Write
   @date       Fri Sep 15 18:47:41 2000
   @author     Tom Goodale
   @desc
   Writes part or all of an HTTP reply.
   @enddesc
   @calls
   @calledby
   @history

   @endhistory

@@*/
int HTTP_Write(httpRequest *request, const char *buffer, size_t count)
{
  int retval = -1;
  size_t bytes_sent = 0;
  size_t halfcount;
  httpSocket *connection = HTTP_Connection( request);
  int done = 0;
  int tmp;
  fd_set this_set;

  if(connection->state == open)
  {
    FD_ZERO (&this_set);
    while(!done)
    {
      /* Try and send the data.  Make sure we don't get a SIGPIPE if the
       * other end is closed.
       */
      retval = send(connection->filedes,
                    buffer+bytes_sent,
                    count-bytes_sent,
                    MSG_NOSIGNAL);

      /* Check for errors */
      if(retval < 0)
      {
        switch(errno)
        {
#ifdef EMSGSIZE
          case EMSGSIZE : /* The socket requires that message be sent atomically,
                           *  and the size of the message to be sent made this
                           * impossible.
                           */
            /* Split message in two and try again. */
            halfcount = count/2;
            retval = HTTP_Write(request, buffer, halfcount);
            if(retval > -1)
            {
              tmp = HTTP_Write(request, buffer+halfcount, count-halfcount);
              if(tmp > -1)
              {
                retval += tmp;
              }
              else
              {
                retval = tmp;
              }
            }
            break;
#endif

#ifdef ENOBUFS
          case ENOBUFS :
            /* The output queue for a network interface was full.  This generally indicates that the  interface  has
             * stopped sending, but may be caused by transient congestion.  (This cannot occur in Linux, packets are
             * just silently dropped when a device queue overflows.)
             */
            retval = 0;
            break;
#endif
          case EBADF    : /* An invalid descriptor was specified. */
#ifdef ENOTSOCK
          case ENOTSOCK : /* The filedescriptor is not a socket. */
#endif
#ifdef ECONNRESET
          case ECONNRESET : /* Connection reset by peer */
#endif
          case EPIPE : /* The local end has been shut down on a connection oriented socket.  In this case the process will also
                        *  receive a SIGPIPE unless MSG_NOSIGNAL is set.
                        */
            SocketClose(connection);
            break;

          case EINTR : /* A signal occurred.*/

          case ENOMEM :  /* No memory available. */

          case EINVAL : /* Invalid argument passed. */

          case EFAULT   : /* An invalid user space address was specified for a parameter. */
            break;

#if 0
          case EAGAIN      :
          case EWOULDBLOCK :
            /* The socket is marked non-blocking and the requested operation would block. */
#endif
          default :
            CCTK_WARN (1, "Error:  unknown errno ");
        }
      }
      else
      {
        bytes_sent += retval;
      }

      if(retval >= 0 && bytes_sent < count)
      {
        /* Wait until can write */
        FD_SET (connection->filedes, &this_set);
        if (select (FD_SETSIZE, NULL, &this_set, NULL, NULL) < 0)
        {
          perror ("select");
          CCTK_Abort(NULL, EXIT_FAILURE);
        }
      }
      else
      {
        done = 1;
      }
    }
  }

#ifdef HTTP_DEBUG
  if(retval != -1)
  {
    fprintf (stderr, "Wrote: `%s'\n", buffer);
  }
  else
  {
    fprintf(stderr, "Couldn't write to the socket 8-(\n");
  }
#endif

  return retval;
}

 /*@@
   @routine    HTTP_Read
   @date       Mon Sep 18 10:14:03 2000
   @author     Tom Goodale
   @desc
   Reads part or all of an HTTP request.
   @enddesc
   @calls
   @calledby
   @history

   @endhistory

@@*/
int HTTP_Read(httpRequest *request, char *buffer, size_t count)
{
  int retval = -1;

  httpSocket *connection = HTTP_Connection( request );

  if(connection->state == open)
  {
    /* Currently don't do anything fancy. */
    retval = recv(connection->filedes, buffer, count,MSG_NOSIGNAL);

#ifdef HTTP_DEBUG
    fprintf (stderr, "Read: `%.*s'\n", retval, buffer);
#endif
  }

  return retval;
}

 /*@@
   @routine    HTTP_Port
   @date       Tue Oct  3 11:24:40 2000
   @author     Tom Goodale
   @desc
   Reports the port the HTTP server is listening on.
   @enddesc
   @calls
   @calledby
   @history

   @endhistory

@@*/
unsigned long int HTTP_Port(void)
{
  return httpport;
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
 /*@@
   @routine    SocketCreate
   @date       Thu Sep 21 15:03:32 2000
   @author     Tom Goodale
   @desc
   Creates an httpSocket structure and links it to the list.
   @enddesc
   @calls
   @calledby
   @history

   @endhistory

@@*/
static httpSocket *SocketCreate(unsigned long int filedes)
{
  httpSocket *this = (httpSocket *)malloc(sizeof (httpSocket));

  if(this)
  {
    this->filedes = filedes;
    this->state = open;
    this->prev = NULL;
    this->next = socklist;

    if(socklist)
    {
      socklist->prev = this;
    }

    socklist = this;

    FD_SET (this->filedes, &active_fd_set);
  }

  return this;
}

 /*@@
   @routine    SocketDestroy
   @date       Thu Sep 21 15:04:03 2000
   @author     Tom Goodale
   @desc
   Destroys an httpSocket structure and removes it from the list.
   @enddesc
   @calls
   @calledby
   @history

   @endhistory

@@*/
static void SocketDestroy(httpSocket *this)
{
  if(this)
  {
    if(this->state == open)
    {
      SocketClose(this);
    }

    if(this->next)
    {
      this->next->prev = this->prev;
    }

    if(this->prev)
    {
      this->prev->next=this->next;
    }
    else
    {
      socklist = this->next;
    }

    free(this);
  }
}

 /*@@
   @routine    SocketClose
   @date       Thu Sep 21 15:04:27 2000
   @author     Tom Goodale
   @desc
   Closes the socket associated with an httpSocket structure.
   @enddesc
   @calls
   @calledby
   @history

   @endhistory

@@*/
static void SocketClose(httpSocket *this)
{
  if(this)
  {
    if(this->state == open)
    {
      Socket_CloseSocket (this->filedes);
      this->state=closed;
      FD_CLR (this->filedes, &active_fd_set);
    }
  }
}
