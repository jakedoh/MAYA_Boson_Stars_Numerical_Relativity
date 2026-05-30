 /*@@
   @header    Expression.h
   @date      Tue Sep 19 22:02:45 2000
   @author    Tom Goodale
   @desc 
   Header for expression stuff.  Was http_Expression.h
   @enddesc
   @version $Header$
 @@*/

#ifndef __HTTP_EXPRESSION_H__
#define __HTTP_EXPRESSION_H__ 1

#ifdef __cplusplus
extern "C" 
{
#endif

char *HTTP_ExpressionParse(const char *expression);

double HTTP_ExpressionEvaluate(char *buffer, 
                               double (eval)(const char *, void *),
                               void *data);

#ifdef __cplusplus
}
#endif

#endif /* __HTTP_EXPRESSION_H__ */
