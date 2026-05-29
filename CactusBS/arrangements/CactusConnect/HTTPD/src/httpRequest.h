 /*@@
   @header    httpRequest.h
   @date      Wed Sep 13 23:49:30 2000
   @author    Tom Goodale
   @desc 
   Was http_Request.h
   @enddesc
   @version $Header$
 @@*/


#ifndef __HTTP_REQUEST_H__
#define __HTTP_REQUEST_H__ 1

typedef struct HTTPSocketTag httpSocket;
typedef struct HTTPArg httpArg;
/* This is the main structure for storing data about a request. */
typedef struct httpRequestTag httpRequest;


#ifdef __cplusplus
extern "C" 
{
#endif

int HTTP_RegisterPage(const char *path,
                   int (*function)(const cGH *, httpRequest *, void *),
                   void *data);

const char *HTTP_ArgumentValue(const httpRequest *request, const char *arg);
const httpArg *HTTP_ArgumentWalk(httpRequest *request, int first);

const char *HTTP_HeaderValue(const httpRequest *request, const char *header);

int HTTP_Write(httpRequest *request, const char *buffer, size_t count);
int HTTP_Read(httpRequest *request, char *buffer, size_t count);

int HTTP_Send( httpRequest *request, const char * message );
void HTTP_SendOKHeader( httpRequest *request);
void HTTP_SendOKRefreshHeader( httpRequest *request, int secs );

unsigned long int HTTP_Port(void);

unsigned int HTTP_MajorVersion( const httpRequest *request );
unsigned int HTTP_MinorVersion( const httpRequest *request );
const char * HTTP_URI( const httpRequest *request );
void HTTP_SetResidual( httpRequest *request, const char *residual );
const char * HTTP_Residual( const httpRequest *request );
unsigned int HTTP_NumArguments( const httpRequest *request );
const char * HTTP_GetArgument( const httpRequest *request, unsigned int index );
httpSocket * HTTP_Connection( const httpRequest *request );
const char * HTTP_ArgName( const httpArg *arg );
const char * HTTP_ArgValue( const httpArg *arg );
#ifdef __cplusplus
}
#endif

#include "httpSString.h"

#endif /* __HTTP_REQUEST_H__ */
