/*@@
   @file      Utils.c
   @date      1991
   @author    John Shalf
   @desc
              Routines which deal with sockets.
   @enddesc
   @history
   @hdate     Thu May 25 13:45:29 2000
   @hauthor   Tom Goodale
   @hdesc     Moved this file into Cactus
   @hdate     Thu 14 May 2002
   @hauthor   Thomas Radke
   @hdesc     Merged with HTTPD socket code
   @endhistory
   @version   $Id$
 @@*/

#include "SocketUtils.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#ifdef HAVE_NETDB_H
#include <netdb.h>
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_WINSOCK2_H
#include <windows.h>
#include <winsock2.h>
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif
#ifdef HAVE_SYS_FILIO_H
#include <sys/filio.h>
#endif
#ifdef HAVE_SYS_IOCTL_H
#include <sys/ioctl.h>
#endif
#ifdef HAVE_SYS_SOCKET_H
#include <sys/socket.h>
#endif
#ifdef HAVE_NETINET_IN_H
#include <netinet/in.h>
#endif
#ifdef HAVE_ARPA_INET_H
#include <arpa/inet.h>
#endif

static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusConnect_Socket_Utils_c)

/********************************************************************
 *********************     Macro Definitions   **********************
 ********************************************************************/
/* SunOS doesn't know INADDR_NONE */
#ifndef INADDR_NONE
#define INADDR_NONE  (-1)
#endif

#ifndef MSG_NOSIGNAL
#define MSG_NOSIGNAL 0
#endif

#ifdef  SOCKET_HAVE_UNIX_SOCKETS
#define IOCTL_SOCKET(a, b, c)    ioctl(a, b, c)
#else
#define errno                    WSAGetLastError()
#define IOCTL_SOCKET(a, b, c)    ioctlsocket (a, b, (u_long *) (c))
#endif

#define GET_INADDR(hostname, sin_addr)                                        \
        {                                                                     \
          struct hostent *host_entry;                                         \
                                                                              \
                                                                              \
          host_entry = gethostbyname (hostname);                              \
          if (host_entry)                                                     \
          {                                                                   \
            memcpy (&sin_addr, host_entry->h_addr_list[0], host_entry->h_length); \
          }                                                                   \
          else                                                                \
          {                                                                   \
            sin_addr.s_addr = inet_addr (hostname);                           \
            if (sin_addr.s_addr == INADDR_NONE)                               \
            {                                                                 \
              CCTK_VWarn (3, __LINE__, __FILE__, CCTK_THORNSTRING,            \
                          "Can't find network address for host '%s'",         \
                          hostname);                                          \
              return (INVALID_SOCKET);                                        \
            }                                                                 \
          }                                                                   \
        }


/********************************************************************
 *********************     Internal Functions   *********************
 ********************************************************************/
static SOCKET TCPOpenSocket (void);
static SOCKET UDPOpenSocket (struct in_addr sin_addr, unsigned int port);


 /*@@
   @routine    Socket_InitializeSocketLayer
   @date       Wed Sep 13 20:39:15 2000
   @author     Tom Goodale
   @desc
               Special code for starting up the socket layer.
               This is only neccessary for Windows.
   @enddesc

   @returntype int
   @returndesc
               0 for success, negative otherwise
   @endreturndesc
@@*/
int Socket_InitializeSocketLayer (void)
{
  int retval;
#ifdef HAVE_WINSOCK2_H
  WSADATA wsaData;


  retval = WSAStartup (MAKEWORD (2, 0), &wsaData);
#else
  retval = 0;
#endif /* HAVE_WINSOCK2_H */

  return (retval == 0 ? 0 : -1);
}


 /*@@
   @routine    Socket_TCPOpenClientSocket
   @date       1991
   @author     John Shalf
   @desc
               Opens a TCP client socket by connecting to a server host with
               given port number.
   @enddesc

   @var        hostname
   @vdesc      hostname of the server to connect to
   @vtype      const char *
   @vio        in
   @endvar
   @var        port
   @vdesc      port number on the server to connect to
   @vtype      unsigned int
   @vio        in
   @endvar

   @returntype SOCKET
   @returndesc
               a valid socket descriptor, or INVALID_SOCKET for failure
   @endreturndesc
@@*/
SOCKET Socket_TCPOpenClientSocket (const char *hostname, unsigned int port)
{
  SOCKET s;
  struct sockaddr_in server;


  memset (&server, 0, sizeof (server));
  server.sin_family = AF_INET;
  server.sin_port = htons (port);
  GET_INADDR (hostname, server.sin_addr);

  s = TCPOpenSocket ();
  if (s != INVALID_SOCKET &&
      connect (s, (struct sockaddr *) &server, sizeof (server)) == SOCKET_ERROR)
  {
    CCTK_VWarn (3, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Couldn't connect to host '%s' port %u: %s",
                hostname, port, strerror (errno));
    Socket_CloseSocket (s);
    s = INVALID_SOCKET;
  }

  return (s);
}


/*@@
   @routine    Socket_TCPOpenServerSocket
   @date       1991
   @author     John Shalf
   @desc
               Opens a TCP server socket on the given port.
               If the port is already taken, it will hunt for the next
               available port. The port is set up to be listened on.
   @enddesc

   @var        port
   @vdesc      port number to bind server socket
   @vtype      unsigned int
   @vio        in
   @endvar

   @returntype SOCKET
   @returndesc
               a valid socket descriptor, or INVALID_SOCKET for failure
   @endreturndesc
@@*/
SOCKET Socket_TCPOpenServerSocket (unsigned int port, unsigned int *hunt,
                                   int backlog)
{
  SOCKET s;
  int is_bound;
  const int on = 1;
  struct sockaddr_in addr;


  s = TCPOpenSocket ();

#ifdef SO_REUSEADDR
  /* try to reuse the port if possible */
  if (s != INVALID_SOCKET)
  {
    if (setsockopt (s, SOL_SOCKET, SO_REUSEADDR, &on, sizeof (on)) ==
        SOCKET_ERROR)
    {
      CCTK_VWarn (3, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Couldn't set SO_REUSEADDR to port %u: %s",
                  port, strerror (errno));
      Socket_CloseSocket (s);
      s = INVALID_SOCKET;
    }
  }
#endif

  if (s != INVALID_SOCKET)
  {
    do
    {
      /* give the socket a name */
      addr.sin_family = AF_INET;
      addr.sin_port = htons (port);
      addr.sin_addr.s_addr = htonl (INADDR_ANY);
      is_bound = bind (s, (struct sockaddr *) &addr, sizeof (addr)) !=
                 SOCKET_ERROR;
      if (! is_bound && hunt)
      {
        /* hunt for a new port */
        CCTK_VWarn (3, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Port %u taken, trying next port", port++);
      }
    } while (! is_bound && hunt);

    if (is_bound)
    {
      if (listen (s, backlog) == SOCKET_ERROR)
      {
        CCTK_VWarn (3, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Couldn't listen on port %u: %s", port, strerror (errno));
        Socket_CloseSocket (s);
        s = INVALID_SOCKET;
      }
      else if (hunt)
      {
        *hunt = port;
      }
    }
    else
    {
      CCTK_VWarn (3, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Couldn't bind socket on port %u: %s", port, strerror(errno));
      Socket_CloseSocket (s);
      s = INVALID_SOCKET;
    }
  }

  return (s);
}


/*@@
   @routine    Socket_CloseSocket
   @date       Tue 15 May 2002
   @author     Thomas Radke
   @desc
               Closes the given socket descriptor.
   @enddesc

   @var        s
   @vdesc      socket to close
   @vtype      SOCKET
   @vio        in
   @endvar
@@*/
void Socket_CloseSocket (SOCKET s)
{
#ifdef SOCKET_HAVE_UNIX_SOCKETS
  close (s);
#else
  closesocket (s);
#endif
}


/*@@
   @routine    Socket_TCPBlockingWrite
   @date       1991
   @author     John Shalf
   @desc
** purpose: Handles error conditions when writing to a file descriptor.
**      If a buffer is too large to write in one block, this routine will
**      divide the buffer into two segments and recursively call itself to
**      send each of the halves.
   @enddesc
@@*/
int Socket_TCPBlockingWrite (SOCKET s, const char *buffer, int buflen)
{
  int n;
  int nstore;


  n = send (s, buffer, buflen, MSG_NOSIGNAL);

  if(n != SOCKET_ERROR)
  {
    return n;
  }
  else
  {
    switch(errno) /* note use of global variable errno */
    {
      case EBADF:
#ifdef PERRORS
        perror("invalid file descriptor");
#endif
        break;
      case EPIPE:
#ifdef PERRORS
        perror("attemped to write to an unconnected socket");
#endif
        break;
      case EFBIG:
#ifdef PERRORS
        perror("datasize too large to write.");
        perror("will attempt recovery by sending in smaller pieces");
#endif
        /* subdivide buffer and call TCPBlockingWrite() recursively */
        nstore=n; /* preserve variable */
        Socket_TCPBlockingWrite(s,buffer,buflen>>1);
        Socket_TCPBlockingWrite(s,buffer+(buflen>>1),buflen-(buflen>>1));
        n=nstore; /* restore variable */
        break;
      case EFAULT:
#ifdef PERRORS
        perror("invalid buffer address");
#endif
        break;
      case EINVAL:
#ifdef PERRORS
        perror("file descriptor points to unusable device.");
#endif
        break;
      case EIO:
#ifdef PERRORS
        perror("an io error occured");
#endif
        break;
#ifdef EWOULDBLOCK
      case EWOULDBLOCK:
#ifdef PERRORS
        perror("Non-blocking I/O is specified by ioctl");
        perror("but socket has no data (would have blocked).");
#endif
#endif
        break;
    }
  }
  return n; /* default, don't know what happened */
}

/*@@
   @routine    Socket_TCPBlockingRead
   @date       1991
   @author     John Shalf
   @desc
** purpose: Handles error conditions when reading from a file descriptor.
**      With TCP stream connections, the record size of the recieved stream
**      may not coincide with the record written.  This routine assembles
**      a fragmented record into a single array by blocking until buflen
**      characters are recieved, or write returns a premature EOF.
   @enddesc
@@*/
int Socket_TCPBlockingRead(SOCKET s,char *buffer, int buflen)
{
  int n,accum=0;


  while ((n = recv (s, buffer, buflen, MSG_NOSIGNAL)) != SOCKET_ERROR)
  {
    buffer = buffer + n;
    buflen -= n;
    accum+=n;
  }
  if (n == SOCKET_ERROR)
  {
    switch(errno) /* note use of global variable errno */
    {
      case EBADF:
#ifdef PERRORS
        perror("invalid file descriptor");
#endif
        break;
      case EFAULT:
#ifdef PERRORS
        perror("invalid buffer address");
#endif
        break;
      case EINTR:
#ifdef PERRORS
        perror("operation interrupted by a signal");
        perror("disable signals and try again");
#endif
        break;
#ifdef EWOULDBLOCK
      case EWOULDBLOCK:
#ifdef PERRORS
        perror("Non-blocking I/O is specified by ioctl");
        perror("but socket has no data (would have blocked).");
#endif
#endif
        break;
    }
  }

  if(n<0)
  {
    return n;
  }
  else
  {
    return accum;
  }
}

int Socket_SetNonBlocking (SOCKET s)
{
  int retval;
  unsigned int on = 1;


  retval = IOCTL_SOCKET (s, FIONBIO, &on) != SOCKET_ERROR ? 0 : -1;

  return(retval);
}


SOCKET Socket_UDPOpenClientSocket (const char *hostname, unsigned int port)
{
  SOCKET s;
  struct in_addr sin_addr;


  GET_INADDR (hostname, sin_addr);

  s = UDPOpenSocket (sin_addr, port);

  return (s);
}


SOCKET Socket_UDPOpenServerSocket (unsigned int port)
{
  SOCKET s;
  struct in_addr sin_addr;


  sin_addr.s_addr = INADDR_ANY;

  s = UDPOpenSocket (sin_addr, port);

  return (s);
}


/********************************************************************
 *********************     Internal Routines   **********************
 ********************************************************************/
static SOCKET TCPOpenSocket (void)
{
  SOCKET s;
  struct protoent *protocol;


  protocol = getprotobyname ("tcp");
  if (! protocol)
  {
    CCTK_WARN (3, "Can't find TCP protocol");
    return (INVALID_SOCKET);
  }

  s = socket (PF_INET, SOCK_STREAM, protocol->p_proto);
  if (s == INVALID_SOCKET)
  {
    CCTK_VWarn (3, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Couldn't create socket: %s", strerror (errno));
  }

  return (s);
}


static SOCKET UDPOpenSocket (struct in_addr sin_addr, unsigned int port)
{
  SOCKET s;
  struct protoent *protocol;
  struct sockaddr_in addr;


  memset (&addr, 0, sizeof (addr));
  addr.sin_family = AF_INET;
  addr.sin_port = htons (port);
  addr.sin_addr = sin_addr;

  protocol = getprotobyname ("udp");
  if (! protocol)
  {
    CCTK_WARN (3, "Can't find UDP protocol");
    return (INVALID_SOCKET);
  }

  s = socket (PF_INET, SOCK_DGRAM, protocol->p_proto);
  if (s == INVALID_SOCKET)
  {
    CCTK_VWarn (3, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Couldn't create socket on port %u: %s",
                port, strerror (errno));
  }
  else if (bind (s, (struct sockaddr *) &addr, sizeof (addr)) == SOCKET_ERROR)
  {
    CCTK_VWarn (3, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Couldn't bind to port %u: %s", port, strerror (errno));
    Socket_CloseSocket (s);
    s = INVALID_SOCKET;
  }

  return (s);
}
