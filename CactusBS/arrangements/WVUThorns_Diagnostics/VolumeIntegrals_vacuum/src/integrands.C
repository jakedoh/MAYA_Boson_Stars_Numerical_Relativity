#ifndef WVU_INTEGRAND_C
#define WVU_INTEGRAND_C

/* Integrand for L2 norms */
inline void L2_integrand(double *VolIntegrand1, int index,double *f,double *x,double *y,double *z) {
  double fL = f[index];
  VolIntegrand1[index] = fL*fL;
}


/* Center of Lapse: */
inline void CoL_integrand(double *VolIntegrand1,double *VolIntegrand2,double *VolIntegrand3,double *VolIntegrand4, int index,double *lapse,double *x,double *y,double *z) {
  double one_minus_lapseL = pow(1.0 - lapse[index],80); // <- Yields pretty consistent results with CoM integrand.
  VolIntegrand1[index] = one_minus_lapseL*x[index];
  VolIntegrand2[index] = one_minus_lapseL*y[index];
  VolIntegrand3[index] = one_minus_lapseL*z[index];
  VolIntegrand4[index] = one_minus_lapseL;
}
#endif
