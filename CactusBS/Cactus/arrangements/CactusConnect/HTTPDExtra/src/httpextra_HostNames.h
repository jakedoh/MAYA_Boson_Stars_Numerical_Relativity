 /*@@
   @header    httpextra_HostNames.h
   @date      Tue Nov  7 23:57:31 2000
   @author    Tom Goodale
   @desc 
   Prototypes for routines which give data about hosts in the machine.
   @enddesc
   @version $Header$
 @@*/

#ifndef _HOSTNAMES_H_
#define _HOSTNAMES_H_ 1

#ifdef __cplusplus
extern "C" 
{
#endif

void HTTPDExtra_CollateHostData(void);

const char *HTTPDExtra_RemoteHostName(int host);

#ifdef __cplusplus
}
#endif

#endif /* _HOSTNAMES_H_ */
