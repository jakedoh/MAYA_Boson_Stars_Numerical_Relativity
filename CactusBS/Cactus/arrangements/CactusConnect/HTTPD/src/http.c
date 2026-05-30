 /*@@
   @file      http.c
   @date      Wed Sep 13 20:44:31 2000
   @author    Tom Goodale
   @desc 
   Stuff to parse and deal with HTTP requests.
   @enddesc
   @version $Header$
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"

#include <stdio.h>
#include <stdlib.h>

#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif
#ifdef HAVE_WINSOCK2_H
#include <winsock2.h>
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <string.h>

#include "util_String.h"

#include "httpd.h"
#include "httpd_Map.h"
#include "SString_Namespace.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusConnect_HTTPD_http_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

typedef struct HTTPArg
{
  /* Public data */
  char *arg;
  char *value;
  
  /* Private data */ 
  struct HTTPArg *next;
  struct HTTPArg *all_next;
} httpArg_PLACEHOLDER;

/* This is the main structure for storing data about a request. */
typedef struct httpRequestTag
{
  /* Public members of the structure. */
  char *body;        /* The body of the request */
  int body_length;   /* How long the body is */

  const char *method;      /* The HTTP method */
  char *uri;         /* The URI */
  const char *residual;    /* What's left of the URI after subtracting the found page */
  
  /* HTTP version numbers */
  int http_major_version;  
  int http_minor_version;

  /* How many arguments there are */
  int n_arguments;

  /* These are all private members of this structure */

  /* The connection data */
  void *connection;

  /* The request header lines */
  uMap headers;

  /* Stuff for arguments */

  /* First a hash table to look the data up quickly */
  uMap arguments;

  /* Now a linked list to allow walking. */
  httpArg *firstarg;
  httpArg *lastarg;
  httpArg *currentarg;

} httpRequest_PLACEHOLDER;

#define BLANK_HTTPREQUEST {NULL, 0, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL } 
  /* Accessors for the httpRequest structure */
unsigned int HTTP_MajorVersion( const httpRequest *request )
{
  return request->http_major_version;
}

unsigned int HTTP_MinorVersion( const httpRequest *request )
{
  return request->http_minor_version;
}

unsigned int HTTP_NumArguments( const httpRequest *request )
{
  return request->n_arguments;
}

const char * HTTP_URI( const httpRequest *request )
{
  return request->uri;
}

const char * HTTP_Residual( const httpRequest *request )
{
  return request->residual;
}

void HTTP_SetResidual( httpRequest *request, const char *residual )
{
  request->residual = residual;
}

const char * HTTP_ArgValue( const httpArg *arg )
{
  return arg->value;
}

const char * HTTP_ArgName( const httpArg *arg )
{
  return arg->arg;
}

httpSocket * HTTP_Connection( const httpRequest *request )
{
  return request->connection;
}

SSBOOL HTTP_GetHeaderValueString(const httpRequest *request,
                const String *header, String *value)
{
  const char * v = HTTP_HeaderValue(request, GetBuffer( header ) );
  if( v )
  {
    SetToCString( value, v );
    return SSTRUE;
  }
  else
    return SSFALSE;
}

struct httpHeader
{
  char *line;
};

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

static char *NextLine(char *buffer, int *start, int size);
static int  DealWithRequest(cGH *cctkGH, httpRequest *request, char *buffer, 
                            int bufsize);
static int  AddHeader(httpRequest *request, const char *line);

static int  StripArgs(httpRequest *request, char *request_uri);
static void Decode(char *string);
static int  InitialiseRequest(httpRequest *request);
static int  ClearRequest(httpRequest *request);

static void DestroyHeader(struct httpHeader *header);
static void DestroyArgument(httpArg *argument);

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

#define INITIAL_SIZE 16
#define MAXMSG  2048
#define EMPTY_STRING  {'\0'}

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    HTTP_ReadFromClient
   @date       Wed Sep 13 20:44:31 2000
   @author     Tom Goodale
   @desc 
   Reads a request from a client socket.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
int HTTP_ReadFromClient(cGH *cctkGH, void *connection)
{
  int retval = -1;
  int keepalive = 0;
  char inbuffer[MAXMSG] = EMPTY_STRING;
  httpRequest request = BLANK_HTTPREQUEST;
  int request_buffer_length = 0;
  char *request_buffer = NULL;
  struct timeval timeout = {1, 0};
  int found_eoh = 0;
  int failure_count = 0;

  InitialiseRequest(&request);

  request.connection = connection;


  while(! found_eoh)
  {
    int nbytes = HTTP_Read(&request, inbuffer, MAXMSG);
   
    if (nbytes < 0)
    {
      /* Read error. */
      perror ("HTTP_Read");
      CCTK_WARN(1, "Error reading from socket");
    }
    else if (nbytes == 0)
    {      
      CCTK_WARN(1, "Got 0 bytes when waiting for end of HTTP request header");
      if(failure_count < 3)
      {
        select (0, NULL, NULL, NULL, &timeout);
        failure_count++;
      }
      else
      {
        /* Give up waiting for this socket */
        CCTK_WARN(1, "Too many failures reading from socket.  Giving up.");
        request_buffer_length = 0;
        break;
      }
    }
    else
    {
      /* Reallocate space for buffer.
       * Add space for one extra char so can null-terminate to allow string search.
       */
      char *tmp = (char *)realloc(request_buffer, (request_buffer_length+nbytes+1)*sizeof(char));
      if(tmp)
      {
        request_buffer = tmp;
      }
      else
      {
        CCTK_WARN(0, "Memory allocation failure.");
      }

      /* Copy the new data into the buffer */
      memcpy(request_buffer+request_buffer_length, inbuffer, nbytes);
      request_buffer_length += nbytes;
      /* Add a null byte in the extra allocated space to allow strstr */
      request_buffer[request_buffer_length] = 0;

      /* Look for end of HTTP header. 
       * Comment out the check for \n\n for strict compliance, but this is
       * useful for debugging with telnet.
       */
      if(strstr(request_buffer, "\r\n\r\n") || strstr(request_buffer, "\n\n"))
      {
        found_eoh = 1;
#ifdef HTTP_DEBUG
        fprintf (stderr, "Found end of HTTP header.\n");
#endif
      }
    }
  }

  if(request_buffer_length > 0)
  {
    /* Data read. */
#ifdef HTTP_DEBUG
    fprintf (stderr, "Server: got message: `%s'\n", request_buffer);
#endif

    keepalive = DealWithRequest(cctkGH, &request, request_buffer, request_buffer_length);
  }

  if(keepalive)
  {
    retval = 0;
  }
  else
  {
    retval = -1;
  }

  if(request_buffer)
  {
    free(request_buffer);
  }

  return retval;
}

 /*@@
   @routine    HTTP_ArgumentValue
   @date       Thu Sep 14 16:40:11 2000
   @author     Tom Goodale
   @desc 
   Returns the value (if any) of an argument.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
const char *HTTP_ArgumentValue(const httpRequest *request, const char *arg)
{
  if(request->arguments)
  {
    const httpArg *value = (httpArg *)Httpd_MapData(request->arguments, strlen(arg), arg);

    if(value)
    {
      return value->value;
    }
  }

  return NULL;
}

 /*@@
   @routine    HTTP_ArgumentWalk
   @date       Sat Sep 16 14:41:01 2000
   @author     Tom Goodale
   @desc 
   Walks the argument list.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
const httpArg *HTTP_ArgumentWalk(httpRequest *request, int first)
{
  if(first || ! (request->currentarg))
  {
    request->currentarg = request->firstarg;
  }
  else
  {
    request->currentarg = request->currentarg->all_next;
  }

  return request->currentarg;
}

 /*@@
   @routine    HTTP_HeaderValue
   @date       Thu Sep 14 19:55:04 2000
   @author     Tom Goodale
   @desc 
   Gets the value of an HTTP header line.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
const char *HTTP_HeaderValue(const httpRequest *request, const char *header)
{
  struct httpHeader *header_line;

  if(request->headers)
  {
    header_line = Httpd_MapData(request->headers, strlen(header), header);

    if(header_line)
    {
      return header_line->line;
    }
  }

  return NULL;
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/

 /*@@
   @routine    NextLine
   @date       Wed Sep 13 20:44:31 2000
   @author     Tom Goodale
   @desc 
   Gets the next line from an HTTP header.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
static char *NextLine(char *buffer, int *start, int size)
{
  char *line = NULL;
  int found = 0;
  char *pos;

  for(pos = buffer + *start; pos-buffer < size-1; pos++)
  {
    if(*pos == '\r' && *(pos+1) == '\n')
    {
      *pos = 0;
      found = 1;
      break;
    }
  }
   
  if(found)
  {
    line = buffer + *start;
    *start = pos-buffer+2;
  }
  else
  {
    *start = size+1;
  }

  return line;
}

 /*@@
   @routine    DealWithRequest
   @date       Wed Sep 13 20:44:31 2000
   @author     Tom Goodale
   @desc 
   Takes an http request buffer and deals with it.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
static int DealWithRequest(cGH *cctkGH, httpRequest *request, char *buffer, int bufsize)
{
  char *line;
  char *tmp = buffer;
  int start = 0;
  char *method = NULL;
  char *request_uri = NULL;
  char *http_version = NULL;
  DECLARE_CCTK_PARAMETERS

#ifdef HTTP_DEBUG
  printf("Dealing with request\n");
#endif

  line = NextLine(tmp, &start, bufsize);

  if(line)
  {
    method = line;

    for(; *line; line++)
    {
      if(*line == ' ')
      {
        *line = 0;
        line++;
        if(!request_uri)
        {
          request_uri = line;
        }
        else
        {
          http_version = line;
          break;
        }
      }
    }
  }

  if(method)
  {
#ifdef HTTP_DEBUG
    printf("Method: %s\n", method);
#endif

    request->method = method; 
  }

  if(request_uri)
  {
    if (verbose)
    {
      CCTK_VInfo (CCTK_THORNSTRING, "URI:    %s", request_uri);
    }

    request->uri = request_uri; 
    StripArgs(request, request->uri);
  }

  if(http_version)
  {
#ifdef HTTP_DEBUG
    printf("HTTP:   %s\n", http_version);
#endif
    sscanf(http_version, "%*5s%d%*1[.]%d", &(request->http_major_version),
                                          &(request->http_minor_version));
  }
  else
  {
    request->http_major_version = 0;
    request->http_minor_version = 90;
  }

#ifdef HTTP_DEBUG
  printf("HTTP Major version: %d\n", request->http_major_version);
  printf("HTTP Minor version: %d\n", request->http_minor_version);
#endif

#ifdef HTTP_DEBUG
  printf("Start of request header\n");
#endif

  while((line = NextLine(tmp, &start, bufsize)) != NULL)
  {
    if(! request->body && *line != 0)
    {
      AddHeader(request, line);
    }

#ifdef HTTP_DEBUG
    if(line)
    {
      printf("Line is '%s'\n", line);
    }
#endif

    if(*line == 0)
    {
#ifdef HTTP_DEBUG
      printf("End of request header\n");

      printf("Start of request body\n");
#endif
      request->body = line + 2;
    }
  }

#ifdef HTTP_DEBUG
  printf("End of request body\n");
#endif

  request->body_length = buffer + bufsize - request->body;

#ifdef HTTP_DEBUG
  printf("Replying...\n");
#endif

  if(!strcmp(request->method, "GET"))
  {
    HTTP_RequestGET(cctkGH, request);
  }
  else
  {
    HTTP_RequestUnsupported(cctkGH, request);
  }
  
#ifdef HTTP_DEBUG
  printf("Dealt with.\n");
#endif

  ClearRequest(request);

  return 0;
}



 /*@@
   @routine    AddHeader
   @date       Wed Sep 13 20:44:31 2000
   @author     Tom Goodale
   @desc 
   Stores an HTTP header line on a request.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
static int AddHeader(httpRequest *request, const char *line)
{
  int keylength = -1;
  struct httpHeader *header_line;
  /* Split the string */
  char *value = strchr(line, ':');

  if(value)
  {
    keylength = value - line;
  
    value++;

    /* Strip off leading whitespace */
    for(; *value && (*value == ' ' || *value == '\t') ; value++);

  }
  else
  {
    CCTK_VWarn(1, __LINE__,__FILE__,CCTK_THORNSTRING,
               "Invalid header line:\n\"%s\"", line );
  }
  
  if(value)
  {
    if(!request->headers)
    {
      /* Need to create the hash table */
      request->headers = Httpd_MapCreate();
    }

    if(request->headers)
    {
      /* Does the line already exist ? */
      header_line = (struct httpHeader *)Httpd_MapData(request->headers, keylength, line);

      if(header_line)
      {
        char *temp = (char *)realloc(header_line->line, strlen(header_line->line)+strlen(value)+2);

        if(temp)
        {
          header_line->line = temp;
          sprintf(header_line->line,"%s,%s", header_line->line,value);
        }
        else
        {
          CCTK_WARN( 0, "Out of memory!\n" );
        }
      }
      else
      {
        /* Ok, must create it */
        header_line = (struct httpHeader *)malloc(sizeof(struct httpHeader));

        if(header_line)
        {
          header_line->line = (char *)malloc(strlen(value)+1);

          if(header_line->line)
          {
            strcpy(header_line->line, value);

            Httpd_MapStore(request->headers, keylength, line, (void *)header_line);
          }
        }
      }
    }
  }

  return 0;
}

 /*@@
   @routine    StripArgs
   @date       Thu Sep 14 16:36:56 2000
   @author     Tom Goodale
   @desc 
   Strips the arguments for a submitted URI from the
   URI and decodes any HTTP encodings.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
static int StripArgs(httpRequest *request, char *request_uri)
{
  char *position;
  /* Do we have arguments ? */
  if((position = strchr(request_uri, '?')))
  {
    char *token;
    /* Mark the end of the URI */
    *position = 0;

    /* Create the hash table */
    request->arguments = Httpd_MapCreate();

    /* Parse the argument list */
    position++;

    token = strtok(position, "&");

    while(token)
    {
      char *value = strchr(token, '=');
      
      if(value)
      {
        httpArg *argument;

        *value = 0;

        value++;

        /* Remove any encoding in the token and value */
        Decode(token);
        Decode(value);

        argument = (httpArg *)malloc(sizeof(httpArg));

        if(argument)
        {
          httpArg *last = NULL;
          httpArg *stored;

          argument->next      = NULL;
          argument->all_next  = NULL;
          argument->arg   = Util_Strdup(token);
          argument->value = Util_Strdup(value);

          if((stored = (httpArg *)Httpd_MapData(request->arguments, strlen(token), token)))
          {
            /* Argument already exists */

            /* Find the last one on the list */
            for(; stored; stored = stored->next)
            {
              last = stored; 
            }

            /* Append to local list */
              /* SW it isn't perfecty obvious this non NULL */
            last->next = argument;
          }
          else
          {
            Httpd_MapStore(request->arguments, strlen(token), token, (void *)argument);
            request->n_arguments++;
          }
          /* Append to global list. */
          if(request->lastarg)
          {
            request->lastarg->all_next = argument;
          }
          request->lastarg = argument;
          if(!request->firstarg)
          {
            request->firstarg = argument;
          }
        }
        else
        {
          CCTK_WARN( 0, "Out of memory!\n" );
        }
      }
      else
      {
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Argument \"%s\" has no value!", token );
      }     
      token=strtok(NULL, "&");  
    }
  }
  else
  {
    request->arguments = NULL;
  }

  /* Finally remove any encoding in the URI */
  Decode(request_uri);

  return 0;
}

 /*@@
   @routine    InitialiseRequest
   @date       Sat Sep 16 14:29:17 2000
   @author     Tom Goodale
   @desc 
   Initialises a request data structure.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
static int InitialiseRequest(httpRequest *request)
{

  /* Initialise the public data */
  request->body        = NULL;
  request->body_length = 0;
  request->method      = NULL;
  request->uri         = NULL;
  request->residual    = NULL;

  request->http_major_version = -1;
  request->http_minor_version = -1;

  request->n_arguments = 0;

  /* Initialise the private data */

  request->connection  = NULL;
  request->headers     = NULL;
  request->arguments   = NULL;
  
  request->firstarg    = NULL;
  request->lastarg     = NULL;
  request->currentarg  = NULL;

  return 0;
}

 /*@@
   @routine    ClearRequest
   @date       Wed Sep 13 20:44:31 2000
   @author     Tom Goodale
   @desc 
   Frees all memory associated with a request.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
static int ClearRequest(httpRequest *request)
{
  if(request->arguments)
  {
    Httpd_MapDestroy(request->arguments, (void (*)(void *))DestroyArgument);
  }

  if(request->headers)
  {
    Httpd_MapDestroy(request->headers, (void (*)(void *))DestroyHeader);
  }

  return 0;
}

 /*@@
   @routine    Decode
   @date       Thu Sep 14 16:25:25 2000
   @author     Tom Goodale
   @desc 
   Decode HTTP encodings -
   + -> ' '
   %2-bytehexadecimal -> ascii
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
static void Decode(char *string)
{
  char *position;
  char *to;
  char hexadecimal[3] = {0, 0, 0};

  for(position=string, to=position; *position; position++, to++)
  {
    switch(*position)
    {
      case '+' : *to = ' '; break;
      case '%' : hexadecimal[0] = *(position+1);
                 hexadecimal[1] = *(position+2);
                 position += 2;
                 *to = strtol(hexadecimal, NULL, 16);
                 break;
      default :  *to = *position;
    }
  }

  *to = 0;
}

 /*@@
   @routine    DestroyHeader
   @date       Wed Sep 13 20:44:31 2000
   @author     Tom Goodale
   @desc 
   Frees all memory associated with a header structure.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
static void DestroyHeader(struct httpHeader *header)
{
  free(header->line);
  free(header);
}

 /*@@
   @routine    DestroyArgument
   @date       Sat Sep 16 14:10:57 2000
   @author     Tom Goodale
   @desc 
   Frees all memory associated with an argument.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
static void DestroyArgument(httpArg *argument)
{
  httpArg *next;

  for(; argument; argument = next)
  {
    next = argument->next;
    free(argument->arg);
    free(argument->value);
    free(argument);
  }
}

