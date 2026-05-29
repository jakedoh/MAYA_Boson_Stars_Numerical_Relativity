 /*@@
   @file      HostNames.c
   @date      Tue Nov  7 17:36:35 2000
   @author    Tom Goodale
   @desc 
   Routines to collect data about all hosts in a Ccatus job.
   @enddesc 
   @version $Header$
 @@*/

#ifndef TEST_HOSTNAMES
#include "cctk.h"
#include "util_Network.h"
#endif

#include <stdlib.h>
#include <string.h>

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
#ifdef HAVE_WINSOCK2_H
#include <winsock2.h>
#endif /* HAVE_WINSOCK2_H */
#include <errno.h>

#ifdef CCTK_MPI
#include "mpi.h"
#endif /* CCTK_MPI */

#include "httpextra_HostNames.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusConnect_HTTPDExtra_HostNames_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

#define HOSTDATALENGTH 255

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

char *hostdata = NULL;

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    HTTPDExtra_CollateHostData
   @date       Tue Nov  7 18:10:01 2000
   @author     Tom Goodale
   @desc 
   Gets data about all hosts in the parallel job.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
void HTTPDExtra_CollateHostData(void)
{
  int rank = 0;
  int nprocs = 1;
  char thisdata[HOSTDATALENGTH+1];

  Util_GetHostName(thisdata, HOSTDATALENGTH);

  thisdata[HOSTDATALENGTH] = 0;

#ifdef CCTK_MPI
  /* Work out how many processes there are. */
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  /* Work out if this is proc 0 or not. */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  if(rank == 0)
  {
    hostdata=(char *)malloc((HOSTDATALENGTH+1)*nprocs);
    
    if(!hostdata)
    {
      CCTK_WARN(0, "Could not allocate memory");
    }
  }
  
#ifdef CCTK_MPI
  MPI_Gather(thisdata, HOSTDATALENGTH+1, MPI_BYTE, 
             hostdata, HOSTDATALENGTH+1, MPI_BYTE,
             0, MPI_COMM_WORLD);
#else
  strncpy(hostdata, thisdata, HOSTDATALENGTH);
#endif

}
 
 /*@@
   @routine    HTTPDExtra_RemoteHostData
   @date       Tue Nov  7 18:10:38 2000
   @author     Tom Goodale
   @desc 
   Gets name of a remote host indexed by host process number.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
const char *HTTPDExtra_RemoteHostName(int host)
{
  return hostdata+host*(HOSTDATALENGTH+1);
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/

#ifdef TEST_HOSTNAMES

int main(int argc, char *argv[])
{
  int rank = 0;
  int nprocs = 1;

#ifdef CCTK_MPI
  MPI_Init(&argc, &argv);

  /* Work out how many processes there are. */
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  
  /* Work out if this is proc 0 or not. */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  HTTPDExtra_CollateHostData();

  if(rank == 0)
  {
    for(rank = 0; rank < nprocs; rank++)
    {
      printf("Host %d is %s\n", rank, HTTPDExtra_RemoteHostName(rank));
    }
  }

#ifdef CCTK_MPI
  MPI_Finalize();
#endif

  return 0;
}

#endif /* TESTHOSTNAMES */
