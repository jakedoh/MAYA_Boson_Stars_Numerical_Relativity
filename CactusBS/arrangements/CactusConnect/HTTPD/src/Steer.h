 /*@@
   @header    Steer.h
   @date      Fri Sep 15 11:31:04 2000
   @author    Tom Goodale
   @desc 
   Webserver parameter steering routines.  Was http_Steer.h; still exported
   by that name
   @enddesc 
   @version $Header$
 @@*/

#ifndef __HTTP_STEER_H__
#define __HTTP_STEER_H__ 1

#ifdef __cplusplus
extern "C" 
{
#endif

int HTTP_SteerQueue(const char *thorn, const char *parameter,
                    const char *value);

int HTTP_SteerDispatch(void);

#ifdef __cplusplus
}
#endif

#endif /* __HTTP_STEER_H__ */
