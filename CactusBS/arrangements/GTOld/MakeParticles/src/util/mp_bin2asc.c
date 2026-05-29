#if 0
 gcc -std=c99 -Wall mp_bin2asc.c -o mp_bin2asc
exit
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#define DIM(x) (sizeof(x)/sizeof((x)[0]))

/* information type encodings */
#include "../makeparticles_types.h"

int main(int argc, char **argv)
{
  // buf has to match CCTK_REAL
  double *buf;
  int max_nvals;
  enum info_type_t info_type;

  max_nvals = 0;
  for(int t = 0 ; t < DIM(info) ; t++)
    if(info[t].nvals > max_nvals)
      max_nvals = info[t].nvals;
  buf = malloc(sizeof(*buf)*max_nvals);
  assert(NULL != buf);
 
  if(argc != 2)
  {
    fprintf(stderr, "usage: %s [-fallback|-luminosity|-sph|-everything] <infile.bin >outfile.asc\n", argv[0]);
    exit(1);
  }

  // decode type argument
  info_type = info_invalid;
  for(int t = 0 ; t < DIM(info) ; t++)
  {
    if(argv[1][0] == '-' && strcmp(argv[1]+1, info[t].name) == 0)
    {
      info_type = info[t].type;
      break;
    }
  }
  if(info_invalid == info_type)
  {
    fprintf(stderr, "unknown type '%s'\n", argv[1]);
    exit(1);
  }

  // output header and data
  fputs(info[info_type].header, stdout);
  while(fread(buf, sizeof(buf[0]), info[info_type].nvals, stdin) == info[info_type].nvals)
  {
    if(info_type == info_fallback)
    {
      fprintf(stdout, info[info_fallback].fmt, buf[0], buf[1], buf[2], buf[3],
          buf[4], buf[5], buf[6], buf[7], buf[8], buf[9], buf[10], buf[11],
          buf[12], buf[13], buf[14], buf[15], buf[16], buf[17], buf[18], buf[19]);
    }
    else if(info_type == info_luminosity)
    {
      fprintf(stdout, info[info_luminosity].fmt, buf[0], buf[1], buf[2], buf[3],
          buf[4], buf[5], buf[6], buf[7], buf[8], buf[9], buf[10], buf[11]);
    }
    else if(info_type == info_sph)
    {
      fprintf(stdout, info[info_luminosity].fmt, buf[0], buf[1], buf[2], buf[3],
          buf[4], buf[5], buf[6], buf[7], buf[8]);
    }
    else if(info_type == info_everything)
    {
      fprintf(stdout, info[info_everything].fmt, buf[0], buf[1], buf[2], buf[3],
          buf[4], buf[5], buf[6], buf[7], buf[8], buf[9], buf[10], buf[11],
          buf[12], buf[13], buf[14], buf[15], buf[16], buf[17], buf[18], buf[19],
          buf[20], buf[21], buf[22]);
    }
    else
      assert(0 == "Internal error: cannot handle type.");

  } 
  if(ferror(stdin))
  {
    fprintf(stderr, "error while reading: %s", strerror(errno));
    exit(1);
  }

  return 0;
}
