#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "util.h"

void util_verify_hdf5_file(const char *infile)
{
  assert(infile);
  const size_t ilen = strlen(infile);

  if(strcmp(infile+ilen-3, ".h5"))
  {
    fprintf(stderr, 
          "Incorrect file extension for input file: %s\n",
          infile);
    exit(-1);
  }

  if (!util_fexist(infile))
  {
    fprintf(stderr, "Input file %s does not exist\n", infile);
    exit(-1);
  }
}

/* make an outfile by replacing FILE.h5 with FILE_MOD.h5
 * check existence of FILE_MOD.h5. return 0 on success,
 */
int util_gen_out_hdf5_filename(const char *infile,
            char **outfile, const char *mod)
{
  assert(infile);
  const size_t ilen = strlen(infile);
  assert(outfile);

  assert(*outfile==NULL);

  char *dup = strdup(infile);

  assert(dup);

  const size_t outsize = ilen + sizeof(mod)+1;
  *outfile = calloc(outsize, 1);
  assert(*outfile);

  dup[ilen-3] = '\0';
  snprintf(*outfile, outsize, "%s_%s.h5", dup, mod);
  free(dup);


  if (util_fexist(*outfile))
  {
    fprintf(stderr, "Output file %s already exists\n", *outfile);
    return -1;
  };
  return 0;
}

int util_fexist(const char *name)
{
  FILE *check = NULL;

  check = fopen(name, "r");

  if (check)
  {
    fclose(check);
  }

  return check ? 1 : 0;
}
