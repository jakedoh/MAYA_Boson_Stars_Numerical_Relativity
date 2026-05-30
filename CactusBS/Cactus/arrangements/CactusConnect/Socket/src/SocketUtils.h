 /*@@
   @header    SocketUtils.h
   @date      1991
   @author    John Shalf
   @desc
              Header file for thorn Socket which defines routines
              to deal with sockets.
   @enddesc
   @version   $Header$
 @@*/

#ifndef _SOCKET_SOCKETUTILS_H_
#define _SOCKET_SOCKETUTILS_H_ 1

#ifdef __cplusplus
extern "C"
{
#endif

#include "cctk.h"     /* HAVE_WINSOCK2_H */

/* check what sockets type we have (Unix or Windows sockets)
   Note that only MS compilers require to use Windows sockets
   but gcc under Windows does not. */
#if ! defined(HAVE_WINSOCK2_H) || defined(__GNUC__)
#define SOCKET_HAVE_UNIX_SOCKETS  1
#endif

/* define the socket data type and constants denoting an invalid socket
   and an error code returned from a socket call */
#ifdef  SOCKET_HAVE_UNIX_SOCKETS
#define SOCKET          int
#define INVALID_SOCKET  -1
#define SOCKET_ERROR    -1
#else
/* the SOCKET type and the constants are defined in this header */
#include <winsock2.h>
#endif


int Socket_InitializeSocketLayer (void);
SOCKET Socket_TCPOpenServerSocket (unsigned int port, unsigned int *hunt,
                                   int backlog);
SOCKET Socket_TCPOpenClientSocket (const char *hostname, unsigned int port);
SOCKET Socket_UDPOpenServerSocket (unsigned int port);
SOCKET Socket_UDPOpenClientSocket (const char *hostname, unsigned int port);
void Socket_CloseSocket (SOCKET s);

int Socket_TCPBlockingWrite (SOCKET s, const char *buffer, int buflen);
int Socket_TCPBlockingRead (SOCKET s, char *buffer, int buflen);

int Socket_SetNonBlocking (SOCKET s);

#ifdef __cplusplus
}
#endif

#endif /* _SOCKET_SOCKETUTILS_H */
