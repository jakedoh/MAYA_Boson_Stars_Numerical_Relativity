 /*@@
   @header    base64.h
   @date      Fri Sep 15 12:03:21 2000
   @author    Tom Goodale
   @desc 
   Header for base64 encoding stuff.
   @enddesc 
   @version $Header$
 @@*/

#ifndef __BASE64_H__
#define __BASE64_H__ 1

#ifdef __cplusplus
extern "C" 
{
#endif

int HTTP_b64_ntop(unsigned char *src, size_t srclength, 
                  char *target, size_t targsize);

int HTTP_b64_pton(const char *src, unsigned char *target, size_t targsize);

#ifdef __cplusplus
}
#endif

#endif /* __BASE64_H__ */

