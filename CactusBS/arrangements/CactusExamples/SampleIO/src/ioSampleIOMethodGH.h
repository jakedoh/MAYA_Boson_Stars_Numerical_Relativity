 /*@@
   @header    ioSampleIOMethodGH.h
   @date      Tue 14 May 2002
   @author    Thomas Radke
   @desc
              The extensions to the GH structure from thorn SampleIO
   @enddesc
   @version   $Header$
 @@*/

#ifndef IOSAMPLEIOMETHOD_IOSAMPLEIOMETHODGH_H_
#define IOSAMPLEIOMETHOD_IOSAMPLEIOMETHODGH_H_ 1

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SampleIOGH
{
  /* default number of times to output */
  int out_every_default;

  /* number of times to output each individual variable */
  CCTK_INT *out_every;

  /* the last iteration output for each variable */
  int *out_last;

  /* current list of variables to output */
  char *out_vars;

  /* stop on I/O parameter parsing errors ? */
  int stop_on_parse_errors;

} ioSampleIOMethodGH;


/* functions which are registered as "SampleIO" I/O method */
int SampleIO_OutputGH (const cGH *GH);
int SampleIO_OutputVarAs (const cGH *GH, const char *fullname,
                          const char *alias);
int SampleIO_TimeFor (const cGH *GH, int vindex);
int SampleIO_TriggerOutput (const cGH *GH, int vindex);

/* other function prototypes for thorn-internal routines */
int SampleIO_Write (const cGH *GH, int vindex);
void SampleIO_CheckSteerableParameters (const cGH *GH, ioSampleIOMethodGH *myGH);


#ifdef __cplusplus
} // extern "C"
#endif

#endif /* IOSAMPLEIOMETHOD_IOSAMPLEIOMETHODGH_H_ */
