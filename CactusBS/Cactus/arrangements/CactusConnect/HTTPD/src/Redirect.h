 /*@@
   @header    Redirect.h
   @date      Fri May 18 11:08:25 2001
   @author    Tom Goodale
   @desc 
   Rdirection stuff.  Was http_Redirect.h.
   @enddesc
   @version
 @@*/
#ifndef _HTTP_REDIRECT_H_
#define _HTTP_REDIRECT_H_ 1

#ifdef __cplusplus
extern "C" 
{
#endif

int HTTP_SetupRedirect(int port, 
                       int queue_size,
                       int hunt);

int HTTP_RegisterRedirect(void);

const char *HTTP_Master(void);
int HTTP_IsServer(void);
int HTTP_IsRedirect(void);



#ifdef __cplusplus
}
#endif

#endif /* _HTTP_REDIRECT_H_ */
