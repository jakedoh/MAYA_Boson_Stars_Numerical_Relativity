 /*@@
   @file      Steer.c
   @date      Fri Sep 15 09:51:01 2000
   @author    Tom Goodale
   @desc 
   Parameter steering interface for the webserver.
   Perhaps should make these required driver functions ?
   @enddesc
   @version $Header$
 @@*/

#ifndef TEST_HTTP_STEER
#include "cctk.h"
#endif

#ifdef HAVE_CAPABILITY_PTHREADS
#include <pthread.h>
#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#ifdef CCTK_MPI
#include "mpi.h"
#endif

#ifndef TEST_HTTP_STEER
#include "cctk_Parameter.h"
#endif

#include "Steer.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusConnect_HTTPD_Steer_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

static int CommunicateBuffer(void);
static int SteerParameters(void);
#ifdef CCTK_MPI
static void ByteSwap(void *buf,int nelements,int elementsize);
#endif /* CCTK_MPI */

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

static char *queuebuffer    = NULL;
static int queuebuffer_size = 0;

#ifdef HAVE_CAPABILITY_PTHREADS
static pthread_mutex_t steer_mutex = PTHREAD_MUTEX_INITIALIZER;
#endif

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    HTTP_SteerQueue
   @date       Fri Sep 15 09:52:13 2000
   @author     Tom Goodale
   @desc 
   Adds a paremeter and its new value to the queue to be steered.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
int HTTP_SteerQueue(const char *thorn, const char *parameter, const char *value)
{
  int retval = 0;
  int buffer_length = 0;
  int parameter_length = 0;

#ifdef HAVE_CAPABILITY_PTHREADS
  pthread_mutex_lock(&steer_mutex);
#endif

  buffer_length = queuebuffer ? strlen(queuebuffer) : 0;

  /* Added length will be lengthof(thorn::par=val\n) */
  parameter_length=strlen(thorn) + 2 + strlen(parameter)+ 1 + strlen(value)+1;

  if(buffer_length+parameter_length+1 > queuebuffer_size)
  {
    char *tmp = (char *)realloc(queuebuffer, buffer_length+parameter_length+1);

    if(tmp)
    {
      queuebuffer = tmp;
      queuebuffer_size = buffer_length+parameter_length+1;
    }
    else
      retval = -1;
  }

  if(!retval)
  {
    sprintf(queuebuffer,"%s%s::%s=%s\n", 
            buffer_length ? queuebuffer : "",
            thorn,parameter, value);
  }

#ifdef HAVE_CAPABILITY_PTHREADS
  pthread_mutex_unlock(&steer_mutex);
#endif

  return retval;
}

 /*@@
   @routine    HTTP_SteerDispatch
   @date       Fri Sep 15 09:58:07 2000
   @author     Tom Goodale
   @desc 
   
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
int HTTP_SteerDispatch(void)
{
#ifdef HAVE_CAPABILITY_PTHREADS
  pthread_mutex_lock(&steer_mutex);
#endif

  CommunicateBuffer();

  SteerParameters();

#ifdef HAVE_CAPABILITY_PTHREADS
  pthread_mutex_unlock(&steer_mutex);
#endif

  return 0;
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/

 /*@@
   @routine    CommunicateBuffer
   @date       Fri Sep 15 10:22:11 2000
   @author     Tom Goodale
   @desc 
   Communicates the contents of the queue buffer to all procs. 
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
static int CommunicateBuffer(void)
{
#ifdef CCTK_MPI
  int rank = 0;
  int nprocs = 1;
  int buffer_size = 0;
  CCTK_INT8 transmit_buffer_size;

  /* Work out how many processes there are. */
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  
  if(nprocs > 1)
  {
    /* Work out if this is proc 0 or not. */
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank == 0)
    {
      /* Sending process */

      /* How much data must we send. */
      buffer_size = queuebuffer ? strlen(queuebuffer) : 0;

      /* Now account for final \0 */
      if(buffer_size)
      {
        buffer_size++;
      }

      transmit_buffer_size = buffer_size;
      
      ByteSwap(&transmit_buffer_size, 1, 8);

      /* Tell other procs how much data to receive. */
      MPI_Bcast(&transmit_buffer_size, 8, MPI_BYTE, 0, MPI_COMM_WORLD);

      if(buffer_size > 0)
      {
        /* Send data to other procs */
        MPI_Bcast(queuebuffer, buffer_size, MPI_BYTE, 0, MPI_COMM_WORLD);
      }
    }
    else
    {
      /* Receiving process */

      /* Find out how much data to receive. */
      MPI_Bcast(&transmit_buffer_size, 8, MPI_BYTE, 0, MPI_COMM_WORLD);

      ByteSwap(&transmit_buffer_size, 1, 8);

      buffer_size = transmit_buffer_size;
      
#ifdef TEST_HTTP_STEER
      fprintf(stderr, "%d - buffer size %d queuebuffer size %d\n", rank, buffer_size, queuebuffer_size);
#endif

      if(buffer_size > queuebuffer_size)
      {
        /* Don't use realloc, as don't need to copy data */
        if(queuebuffer)
        {
          free(queuebuffer);
        }
        
        queuebuffer = (char *)malloc(buffer_size);

        if(!queuebuffer)
        {
          /* Have no choice except to abort as other processes would
             block. 
          */
          MPI_Abort(MPI_COMM_WORLD, 99);
        }

        queuebuffer_size = buffer_size;
      }

      /* Receive data from processor 0 */

      if(buffer_size > 0)
      {
        MPI_Bcast(queuebuffer, buffer_size, MPI_BYTE, 0, MPI_COMM_WORLD);
      }
    }
  }

#endif /* CCTK_MPI */

  return 0;
}

 /*@@
   @routine    SteerParameters
   @date       Fri Sep 15 10:38:40 2000
   @author     Tom Goodale
   @desc 
   Steers parameters from queuebuffer.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
static int SteerParameters(void)
{
  int retval = 0;
  char *value = NULL;
 
  if(queuebuffer)
  {
    char *token = strtok(queuebuffer, "\n");

    while(token)
    {
      char *thorn = token;

      char *parameter = strchr(token, ':');

      if(parameter)
      {
        /* Get rid of two colons */
        *parameter = 0;
        parameter++;
        *parameter = 0;
        parameter++;

        value = strchr(parameter,'=');

        if(value)
        {
          *value = 0;
          value++;

          CCTK_ParameterSet(parameter, /* The name of the parameter  */
                            thorn,     /* The originating thorn      */
                            value);    /* The value of the parameter */

          retval++;
        }
      }
      token = strtok(NULL,"\n");
    }
  }

  if(queuebuffer)
  {
    *queuebuffer = 0;
  }

  return retval;
}

 /*@@
   @routine    ByteSwap
   @date       Mon Oct 09 2000
   @author     John Shalf
   @desc 
   
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
   @var     buf
   @vdesc   buffer to byteswap
   @vtype   void *
   @vio     inout
   @vcomment 
 
   @endvar 
   @var     nelements
   @vdesc   number of data elements in array
   @vtype   int
   @vio     in
   @vcomment 
 
   @endvar 
   @var     elementsize
   @vdesc   size of each element in buffer
   @vtype   int
   @vio     in
   @vcomment 
 
   @endvar 

@@*/
#ifdef CCTK_MPI
static void ByteSwap(void *buf,int nelements,int elementsize)
{
#ifndef WORDS_BIGENDIAN
  char *buffer=(char *)buf;
  int i;
  int s,d;
  
  for(i=0;i<nelements;i++,buffer+=elementsize)
  {
    /* do the swap thing on each element */
    for(s=0, d=elementsize-1; s < d; s++,d--)
    {
      char      c=buffer[s]; 
      buffer[s]=buffer[d]; 
      buffer[d]=c;
    }
  } 
#endif
}
#endif /* CCTK_MPI */

#ifdef TEST_HTTP_STEER

int CCTK_ParameterSet(char *parameter, char *thorn, char *value)
{
  int rank = 0;

#ifdef CCTK_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  fprintf(stderr, "%d - %s::%s='%s'\n", rank, thorn, parameter, value);

  return 0;
}

int main(int argc, char *argv[])
{
  int i;
  int rank = 0;
  char value[20] = {'\0'};

#ifdef CCTK_MPI
  MPI_Init(&argc, &argv);
#endif

#ifdef CCTK_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  for(i = 0; i < 6; i++)
  {
    if(rank == 0)
    {
      sprintf(value,"iter %d\n", i);

      switch(i)
      {
        case 3 :
        case 5 :
          HTTP_SteerQueue("a","b",value);
          HTTP_SteerQueue("b","c",value);
          HTTP_SteerQueue("c","d",value);
        case 1 :
          HTTP_SteerQueue("d","e",value);
          HTTP_SteerQueue("e","f",value);
          HTTP_SteerQueue("f","g",value);
        default :
      }
    }

    HTTP_SteerDispatch();
  }

#ifdef CCTK_MPI
  MPI_Finalize();
#endif
  
  return 0;
}
      
#endif /* TEST_HTTP_STEER */    
