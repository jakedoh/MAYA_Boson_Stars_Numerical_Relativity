 /*@@
   @file      SampleBnd.h
   @date      6 May 2003
   @author    David Rideout
   @desc
              Prototypes for boundary routines
   @enddesc
 @@*/


#ifndef _SAMPLEBND_H_
#define _SAMPLEBND_H_


#ifdef __cplusplus
extern "C"
{
#endif

/* data type for pointer to function which implements a physical boundary 
   condition: */
typedef CCTK_INT (*const phys_bc_fn_ptr)(const CCTK_POINTER_TO_CONST, const CCTK_INT, 
                                   const CCTK_INT *, const CCTK_INT *, 
                                   const CCTK_INT *, const CCTK_INT *);

/* prototype for routine registed as providing 'LinExtrap' boundary condition */
int BndLinExtrap (const CCTK_POINTER_TO_CONST GH, const CCTK_INT num_vars, 
                  const CCTK_INT *var_indices,  const CCTK_INT *faces, 
                  const CCTK_INT *widths, const CCTK_INT *table_handles);

#ifdef __cplusplus
}
#endif

#endif  /* _SAMPLEBND_H_ */
