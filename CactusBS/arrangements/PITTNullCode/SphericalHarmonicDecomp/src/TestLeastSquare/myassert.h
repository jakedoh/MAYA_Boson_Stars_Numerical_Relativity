#ifndef MYASSERT_H
#define MYASSERT_H

#include<stdio.h>
#include<stdlib.h>

#define myassert(x)                                        \
  do {                                                     \
       const long int  _err = (long int) (x);              \
       if (!_err)                                          \
       {                                                   \
         fprintf(stderr,                                   \
            "Test '%s' failed at line '%d' in file '%s'\n",\
            #x, __LINE__, __FILE__);                       \
         fflush(stderr);                                   \
         exit(-1);                                         \
       }                                                   \
    } while(0)

#endif
