/* Macro definitions to speed up the code for thorn ReconstructPPM */

/* Spatial Determinant */
#define SPATIAL_DET(gxx_,gxy_,gxz_,gyy_,gyz_,gzz_) \
                   (-(gxz_)**2*(gyy_) + 2*(gxy_)*(gxz_)*(gyz_) - (gxx_)*(gyz_)**2 - (gxy_)**2*(gzz_) \
                   + (gxx_)*(gyy_)*(gzz_))

