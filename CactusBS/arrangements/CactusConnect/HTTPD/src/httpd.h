 /*@@
   @header    httpd.h
   @date      Wed Sep 13 20:15:23 2000
   @author    Tom Goodale
   @desc 
   
   @enddesc
   @version $Header$
 @@*/

#include "httpRequest.h"

#ifndef __HTTPD_H__
#define __HTTPD_H__ 1

typedef struct
{
  int steer;
  int paused;
  int terminate;
  int timeout_seconds;
  int timeout_useconds;
} httpState;

#ifdef __cplusplus
extern "C" 
{
#endif

int HTTP_ReadFromClient (cGH *cctkGH, void *connection);

int HTTP_SetupServer(int port, int queue_size, int hunt);
int HTTP_ShutdownServer(void);

int HTTP_Poll(cGH *cctkGH, long sec, long usec);

/* Request stuff */
int HTTP_RequestGET(cGH *cctkGH, httpRequest *request);
int HTTP_RequestUnsupported(cGH *cctkGH, httpRequest *request);

/* State stuff */ 
int HTTP_UpdateState(cGH *cctkGH, httpState *state);
int HTTP_Terminate(cGH *cctkGH);

/* Content stuff */
int HTTP_RegisterPages(void);

#ifdef __cplusplus
}
#endif

#endif /* __HTTP_H__ */
