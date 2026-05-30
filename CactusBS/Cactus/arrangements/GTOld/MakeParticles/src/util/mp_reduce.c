#if 0
 gcc -lm -std=c99 -Wall mp_bin2roman.c -o mp_bin2roman
exit
#endif
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <assert.h>

#define DIM(x) (sizeof(x)/sizeof((x)[0]))
#define SQR(x) ((x)*(x))

/* information type encodings */
#include "../makeparticles_types.h"

// find column of data named "columns"
int strwcmp(const char *a, const char *b)
{
  while(*a && *b && !isspace(*a) && !isspace(*b) && *a == *b) {
    a += 1;
    b += 1;
  }
  return (!*a || isspace(*a)) && (!*b || isspace(*b)) ? 0 : *a-*b;
}
int find_column(const char *header, const char *column)
{
  int col = 0;

  while((header = strchr(header, ':')) != NULL && strwcmp(header+1, column) != 0) {
    col += 1;
    header += 1;
  }

  return header != NULL ? col : -1;
}

// very simple minded converter to produce files maybe more suitable for Roman Shcherbakov
int main(int argc, char **argv)
{
  // buf has to match CCTK_REAL
  double *in_buf;
  int max_nvals;

  enum out_ids {o_x, o_y, o_z, o_rho, o_eps};
  float out_buf[5];
  enum info_type_t info_type;
  // indices we care about
  int i_x, i_y, i_z, i_rho, i_eps;
 
  if(argc != 2)
  {
    fprintf(stderr, "usage: %s [-fallback|-luminosity|-sph|-everything] <infile.dat >outfile.dat\n", argv[0]);
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
  
  i_x = find_column(info[info_type].header, "x");
  i_y = find_column(info[info_type].header, "y");
  i_z = find_column(info[info_type].header, "z");
  i_rho = find_column(info[info_type].header, "rho");
  i_eps = find_column(info[info_type].header, "eps");
  assert(i_x >= 0);
  assert(i_y >= 0);
  assert(i_z >= 0);
  assert(i_rho >= 0);
  assert(i_eps >= 0);

  max_nvals = 0;
  for(int t = 0 ; t < DIM(info) ; t++)
    if(info[t].nvals > max_nvals)
      max_nvals = info[t].nvals;
  in_buf = malloc(sizeof(*in_buf)*max_nvals);
  assert(NULL != in_buf);

  while(fread(in_buf, sizeof(*in_buf), info[info_type].nvals, stdin) == info[info_type].nvals)
  {
    out_buf[o_x] = (float)in_buf[i_x];
    out_buf[o_y] = (float)in_buf[i_y];
    out_buf[o_z] = (float)in_buf[i_z];

    out_buf[o_rho] = (float)in_buf[i_rho];
    out_buf[o_eps] = (float)in_buf[i_eps];

    fwrite(out_buf, sizeof(out_buf), 1, stdout);
  } 
  if(!feof(stdin))
  {
    fprintf(stderr, "error while reading: %s\n", strerror(errno));
    exit(1);
  }

  return 0;
}
