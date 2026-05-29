enum info_type_t {info_invalid=-1, info_fallback=0, info_luminosity, info_sph, info_everything};
static const struct {
  const enum info_type_t type;
  const char *name;
  const char *header;
  const char *fmt;
  const int nvals;
} info[] = {
  {info_fallback, "fallback time", 
   "# 1:x 2:y 3:z 4:v^x 5:v^y 6:v^z 7:gxx 8:gxy 9:gxz 10:gyy 11:gyz 12:gzz 13:beta^x 14:beta^y 15:beta^z 16:alpha 17:mass 18:energy_per_unit_mass 19:rho 20:eps\n",
   "%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
   20 /* number of output columns */
  },
  {info_luminosity, "luminosity", 
   "% 1:x 2:y 3:z 4:rho 5:eps 6:mass 7:detg 8:w_lorentz 9:Df4Y 10:Df4z 11:dx 12:weight*dV 13:velx 14:vely 15:velz\n",
   "%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
   15 /* number of output columns */
  },
  {info_sph, "sph particles", 
   "# 1:x 2:y 3:z 4:v^x 5:v^y 6:v^z 7:rho 8:eps 9:mass\n",
   "%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
   9 /* number of output columns */
  },
  {info_everything, "everything", 
   "# 1:x 2:y 3:z 4:v^x 5:v^y 6:v^z 7:gxx 8:gxy 9:gxz 10:gyy 11:gyz 12:gzz 13:beta^x 14:beta^y 15:beta^z 16:alpha 17:rho 18:eps 19:press 20:edge 21:weight 22:dx 23:reflevel\n",
   "%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %g\n",
   23
  },
};
