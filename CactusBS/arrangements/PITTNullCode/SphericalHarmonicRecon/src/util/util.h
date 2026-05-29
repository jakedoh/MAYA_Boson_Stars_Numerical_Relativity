#ifndef UTIL_H
#define UTIL_H

extern void util_verify_hdf5_file(const char *infile);

extern int util_gen_out_hdf5_filename(const char *infile,
            char **outfile, const char *mod);

extern int util_fexist(const char *name);

#endif
