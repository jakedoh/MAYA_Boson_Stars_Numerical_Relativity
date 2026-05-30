#include "Geodesics.hh"
#include "carpet.hh"
//#include <gsl/gsl_sf_synchrotron.h>

extern "C" void Geodesics_CalcRHSEM(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(start_t > cctk_time) return; // Save some cpu time if we aren't evolving yet...

  int lsh[1];
  int ierr = CCTK_GrouplshGN (cctkGH, 1, lsh, "Geodesics::geo_pos");
  int lnum_geo = lsh[0];

  // If we arent on the finest refinement level, set the RHS to zero.
  // this keeps the geodesics synced up with cctk_time, otherwise
  // they will be updated once for each refinement level using that reflevel's
  // timestep.
  //
  // See http://cactuscode.org/pipermail/developers/2007-September/005433.html
  if (Carpet::reflevel != Carpet::reflevels - 1) {
    for( int i=0; i < lnum_geo; i++ ) {
      geo_x_rhs[i]=0.0;
      geo_y_rhs[i]=0.0;
      geo_z_rhs[i]=0.0;

      geo_ux_rhs[i]=0.0;
      geo_uy_rhs[i]=0.0;
      geo_uz_rhs[i]=0.0;

      if(CCTK_Equals(geo_moltest,"yes"))  geo_t_rhs[i] = 0.0;
    }
  } 
  else 
  {

    if( verbose > 2)
      CCTK_VInfo(CCTK_THORNSTRING,
          "Geodesics_CalcRHS(): starting");

    const CCTK_INT num_input_arrays = 17;
    const CCTK_INT num_output_arrays = 50;

    /*
    // ADM vars
    CCTK_STRING slapse = "admbase::alp";
    CCTK_STRING sgxx = "admbase::gxx";
    CCTK_STRING sgyy = "admbase::gyy";
    CCTK_STRING sgzz = "admbase::gzz";
    CCTK_STRING sgxy = "admbase::gxy";
    CCTK_STRING sgxz = "admbase::gxz";
    CCTK_STRING sgyz = "admbase::gyz";
    CCTK_STRING sbetax = "admbase::betax";
    CCTK_STRING sbetay = "admbase::betay";
    CCTK_STRING sbetaz = "admbase::betaz";
     */

    // BSSN vars
CCTK_STRING slapse, sgxx, sgyy, sgzz, sgxy, sgxz, sgyz, sbetax, sbetay, sbetaz, schi;
if(CCTK_IsThornActive("Kranc2BSSNChiMatter") != 0) {
     slapse = "Kranc2BSSNChiMatter::alpha";
     sgxx = "Kranc2BSSNChiMatter::h11";
     sgyy = "Kranc2BSSNChiMatter::h22";
     sgzz = "Kranc2BSSNChiMatter::h33";
     sgxy = "Kranc2BSSNChiMatter::h21";
     sgxz = "Kranc2BSSNChiMatter::h31";
     sgyz = "Kranc2BSSNChiMatter::h32";
     sbetax = "Kranc2BSSNChiMatter::beta1";
     sbetay = "Kranc2BSSNChiMatter::beta2";
     sbetaz = "Kranc2BSSNChiMatter::beta3";
     schi = "Kranc2BSSNChiMatter::chi";
}
if(CCTK_IsThornActive("Kranc2BSSNChi") != 0) {
     slapse = "Kranc2BSSNChi::alpha";
     sgxx = "Kranc2BSSNChi::h11";
     sgyy = "Kranc2BSSNChi::h22";
     sgzz = "Kranc2BSSNChi::h33";
     sgxy = "Kranc2BSSNChi::h21";
     sgxz = "Kranc2BSSNChi::h31";
     sgyz = "Kranc2BSSNChi::h32";
     sbetax = "Kranc2BSSNChi::beta1";
     sbetay = "Kranc2BSSNChi::beta2";
     sbetaz = "Kranc2BSSNChi::beta3";
     schi = "Kranc2BSSNChi::chi";
}

    // E&M vars
    CCTK_STRING sbx = "EMevoPlasmaChi::Bff1";
    CCTK_STRING sby = "EMevoPlasmaChi::Bff2";
    CCTK_STRING sbz = "EMevoPlasmaChi::Bff3";
    CCTK_STRING sex = "EMevoPlasmaChi::Eff1";
    CCTK_STRING sey = "EMevoPlasmaChi::Eff2";
    CCTK_STRING sez = "EMevoPlasmaChi::Eff3";


    //    CCTK_STRING sjx = "EMevoPlasmaChi::current1";
    //    CCTK_STRING sjy = "EMevoPlasmaChi::current2";
    //    CCTK_STRING sjz = "EMevoPlasmaChi::current3";


    const CCTK_INT input_array_indices[num_input_arrays]
      = { CCTK_VarIndex(slapse), CCTK_VarIndex(sgxx), CCTK_VarIndex(sgyy),
        CCTK_VarIndex(sgzz), CCTK_VarIndex(sgxy), CCTK_VarIndex(sgxz),
        CCTK_VarIndex(sgyz), CCTK_VarIndex(sbetax), CCTK_VarIndex(sbetay),
        CCTK_VarIndex(sbetaz), CCTK_VarIndex(schi), CCTK_VarIndex(sbx),
        CCTK_VarIndex(sby), CCTK_VarIndex(sbz), CCTK_VarIndex(sex),
        CCTK_VarIndex(sey), CCTK_VarIndex(sez)//, CCTK_VarIndex(sjx),
        //CCTK_VarIndex(sjy), CCTK_VarIndex(sjz), 
      };

    /*    for(int i = 0 ; i < num_input_arrays ; i++)
          {
          if (input_array_indices[i] < 0)
          {
          const CCTK_STRING varnames[num_input_arrays]
          = { slapse,
          sgxx, sgyy, sgzz, sgxy, sgxz, sgyz,
          sbetax, sbetay, sbetaz, schi,
          sbx, sby, sbz,
          sex, sey, sez };

          CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
          "Geodesics: couldn't find variable \"%s\"!",
          varnames[i]);         
          return;
          }
          }
     */
    const void* interp_coords[3]
      = { (const void *) geo_x,
        (const void *) geo_y,
        (const void *) geo_z };

    // Output arrays for the interpolator
    CCTK_REAL a[lnum_geo], a_dx[lnum_geo], a_dy[lnum_geo], a_dz[lnum_geo];
    CCTK_REAL g_xx[lnum_geo], g_xx_dx[lnum_geo], g_xx_dy[lnum_geo], g_xx_dz[lnum_geo];
    CCTK_REAL g_yy[lnum_geo], g_yy_dx[lnum_geo], g_yy_dy[lnum_geo], g_yy_dz[lnum_geo];
    CCTK_REAL g_zz[lnum_geo], g_zz_dx[lnum_geo], g_zz_dy[lnum_geo], g_zz_dz[lnum_geo];
    CCTK_REAL g_xy[lnum_geo], g_xy_dx[lnum_geo], g_xy_dy[lnum_geo], g_xy_dz[lnum_geo];
    CCTK_REAL g_xz[lnum_geo], g_xz_dx[lnum_geo], g_xz_dy[lnum_geo], g_xz_dz[lnum_geo];
    CCTK_REAL g_yz[lnum_geo], g_yz_dx[lnum_geo], g_yz_dy[lnum_geo], g_yz_dz[lnum_geo];
    CCTK_REAL beta_x[lnum_geo], beta_x_dx[lnum_geo], beta_x_dy[lnum_geo], beta_x_dz[lnum_geo];
    CCTK_REAL beta_y[lnum_geo], beta_y_dx[lnum_geo], beta_y_dy[lnum_geo], beta_y_dz[lnum_geo];
    CCTK_REAL beta_z[lnum_geo], beta_z_dx[lnum_geo], beta_z_dy[lnum_geo], beta_z_dz[lnum_geo];
    CCTK_REAL chiL[lnum_geo], chi_dx[lnum_geo], chi_dy[lnum_geo], chi_dz[lnum_geo];
    CCTK_REAL bx[lnum_geo], by[lnum_geo], bz[lnum_geo];
    CCTK_REAL ex[lnum_geo], ey[lnum_geo], ez[lnum_geo];
    //    CCTK_REAL jx[lnum_geo], jy[lnum_geo], jz[lnum_geo];
    //    CCTK_REAL pex[lnum_geo], pey[lnum_geo], pez[lnum_geo];


    void * output_arrays[num_output_arrays]
      = { (void *) a, (void *) a_dx, (void *) a_dy, (void *) a_dz,
        (void *) g_xx, (void *) g_xx_dx, (void *) g_xx_dy, (void *) g_xx_dz,
        (void *) g_yy, (void *) g_yy_dx, (void *) g_yy_dy, (void *) g_yy_dz,
        (void *) g_zz, (void *) g_zz_dx, (void *) g_zz_dy, (void *) g_zz_dz,
        (void *) g_xy, (void *) g_xy_dx, (void *) g_xy_dy, (void *) g_xy_dz,
        (void *) g_xz, (void *) g_xz_dx, (void *) g_xz_dy, (void *) g_xz_dz,
        (void *) g_yz, (void *) g_yz_dx, (void *) g_yz_dy, (void *) g_yz_dz,
        (void *) beta_x, (void *) beta_x_dx, (void *) beta_x_dy, (void *) beta_x_dz,
        (void *) beta_y, (void *) beta_y_dx, (void *) beta_y_dy, (void *) beta_y_dz,
        (void *) beta_z, (void *) beta_z_dx, (void *) beta_z_dy, (void *) beta_z_dz,
        (void *) chiL, (void *) chi_dx, (void *) chi_dy, (void *) chi_dz,
        (void *) bx, (void *) by, (void *) bz,
        (void *) ex, (void *) ey, (void *) ez,
        //        (void *) jx, (void *) jy, (void *) jz 
        // (void *) pex, (void *) pey, (void *) pez
      };

    const CCTK_INT output_array_types[num_output_arrays]
      = { CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        //        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        //        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL };

    const CCTK_INT operand_indices[num_output_arrays]
      = { 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3,
        4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7,
        8, 8, 8, 8, 9, 9, 9, 9, 10,10,10,10,
        11,12,13,14,15,16//,17,18,19 
      };

    const CCTK_INT operation_codes[num_output_arrays]
      = { 0, 1, 2, 3,
        0, 1, 2, 3,
        0, 1, 2, 3,
        0, 1, 2, 3,
        0, 1, 2, 3,
        0, 1, 2, 3,
        0, 1, 2, 3,
        0, 1, 2, 3,
        0, 1, 2, 3,
        0, 1, 2, 3,
        0, 1, 2, 3,
        0,0,0,0,0,0//,0,0,0 
      };

    const int operator_handle = CCTK_InterpHandle(interpolator_name);

    if (operator_handle < 0)
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
          "Geodesics_Interp(): couldn't find interpolator \"%s\"!",
          interpolator_name);         /*NOTREACHED*/

    // This tells the interpolater what arguments it should use, e.g.
    // interpolation order.
    int param_table_handle = Util_TableCreate(UTIL_TABLE_FLAGS_DEFAULT);

    // User parameters can be set in the parameter file
    ierr = Util_TableSetFromString(param_table_handle, interpolator_pars);

    if (ierr < 0)
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
          "Geodesics_Interp(): bad interpolator parameter(s) \"%s\"!",
          interpolator_pars);         /*NOTREACHED*/

    if (interp_num_timelevels != -1)
    {
      // Don't interpolate in time.  I do this because I think that carpet's
      // dertermination of when it should interpolate in time is bad.  It
      // should never interpolate this in time.
      ierr = Util_TableSetInt(param_table_handle, interp_num_timelevels,
          "InterpNumTimelevels");
    }

    // Perform the interpolation from the global mesh (this is pretty important)
    ierr += Util_TableSetInt(param_table_handle, 1,
        "want_global_mode");

    // Dont show warnings from the interpolator if we arent verbose (they can output ALOT to stdout)
    if(verbose == 0) ierr += Util_TableSetInt(param_table_handle, 1, "suppress_warnings");

    // TODO: This error checking prolly shouldnt be +='ed liked this
    ierr += Util_TableSetIntArray(param_table_handle,  
        num_output_arrays, operand_indices,  
        "operand_indices");

    ierr += Util_TableSetIntArray(param_table_handle,  
        num_output_arrays, operation_codes,  
        "operation_codes");

    if (ierr < 0)
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
          "Geodesics_Interp(): bad interpolator parameter(s) \"%d\"!", 1);        

    const int coord_system_handle = CCTK_CoordSystemHandle(coord_system);
    if (coord_system_handle < 0)
      CCTK_VWarn(-1, __LINE__, __FILE__, CCTK_THORNSTRING,
          "Geodesics_Interp(): can't get coordinate system handle for\n"
          "                           for coordinate system \"%s\"!",
          coord_system);   

    // Call the interpolator
    ierr = CCTK_InterpGridArrays(cctkGH,
        3,
        operator_handle,
        param_table_handle,
        coord_system_handle,
        lnum_geo,
        CCTK_VARIABLE_REAL,
        interp_coords,
        num_input_arrays,
        input_array_indices,
        num_output_arrays,
        output_array_types,
        output_arrays);

    // Some vars that dont need to be stored as arrays
    CCTK_REAL e4phi;  
    CCTK_REAL iDetG; 
    CCTK_REAL iGxx, iGyy, iGzz, iGxy, iGxz, iGyz;
    CCTK_REAL b_x, b_y, b_z;
    CCTK_REAL e_x, e_y, e_z;

    for( int i=0; i < lnum_geo; i++ )  {
      if(active_geo[i]==-1) { //|| start_t > cctk_time){ //we shouldnt get this far if start_t > cctk_time...
        geo_x_rhs[i] = 0.0;
        geo_y_rhs[i] = 0.0;
        geo_z_rhs[i] = 0.0;

        geo_ux_rhs[i] = 0.0;
        geo_uy_rhs[i] = 0.0;
        geo_uz_rhs[i] = 0.0;

        geo_iut[i] = 0.0;
        geo_iux[i] = 0.0;
        geo_iuy[i] = 0.0;
        geo_iuz[i] = 0.0;

      } else {

        // Switch metric (and derivatives) to ADM for calculation
        // put a lower bound on chi (otherwise we may divide by zero)
        CCTK_REAL chi = fmax(chiL[i],geo_chieps); 
        e4phi = 1.0/chi;

        g_xx[i] *= e4phi;
        g_yy[i] *= e4phi;
        g_zz[i] *= e4phi;
        g_xy[i] *= e4phi;
        g_xz[i] *= e4phi;
        g_yz[i] *= e4phi;

        g_xx_dx[i] = e4phi*(g_xx_dx[i] - g_xx[i]*chi_dx[i]);
        g_yy_dx[i] = e4phi*(g_yy_dx[i] - g_yy[i]*chi_dx[i]);
        g_zz_dx[i] = e4phi*(g_zz_dx[i] - g_zz[i]*chi_dx[i]);
        g_xy_dx[i] = e4phi*(g_xy_dx[i] - g_xy[i]*chi_dx[i]);
        g_xz_dx[i] = e4phi*(g_xz_dx[i] - g_xz[i]*chi_dx[i]);
        g_yz_dx[i] = e4phi*(g_yz_dx[i] - g_yz[i]*chi_dx[i]);

        g_xx_dy[i] = e4phi*(g_xx_dy[i] - g_xx[i]*chi_dy[i]);
        g_yy_dy[i] = e4phi*(g_yy_dy[i] - g_yy[i]*chi_dy[i]);
        g_zz_dy[i] = e4phi*(g_zz_dy[i] - g_zz[i]*chi_dy[i]);
        g_xy_dy[i] = e4phi*(g_xy_dy[i] - g_xy[i]*chi_dy[i]);
        g_xz_dy[i] = e4phi*(g_xz_dy[i] - g_xz[i]*chi_dy[i]);
        g_yz_dy[i] = e4phi*(g_yz_dy[i] - g_yz[i]*chi_dy[i]);

        g_xx_dz[i] = e4phi*(g_xx_dz[i] - g_xx[i]*chi_dz[i]);
        g_yy_dz[i] = e4phi*(g_yy_dz[i] - g_yy[i]*chi_dz[i]);
        g_zz_dz[i] = e4phi*(g_zz_dz[i] - g_zz[i]*chi_dz[i]);
        g_xy_dz[i] = e4phi*(g_xy_dz[i] - g_xy[i]*chi_dz[i]);
        g_xz_dz[i] = e4phi*(g_xz_dz[i] - g_xz[i]*chi_dz[i]);
        g_yz_dz[i] = e4phi*(g_yz_dz[i] - g_yz[i]*chi_dz[i]);

        // Calc inverse metric
        iDetG = 1.0/(g_xx[i]*g_yy[i]*g_zz[i] + 2.0*g_xy[i]*g_xz[i]*g_yz[i] 
            - (pow(g_xx[i]*g_yz[i],2) + pow(g_yy[i]*g_xz[i],2) + pow(g_zz[i]*g_xy[i],2)));

        iGxx = iDetG*(g_yy[i]*g_zz[i] - g_yz[i]*g_yz[i]);
        iGyy = iDetG*(g_xx[i]*g_zz[i] - g_xz[i]*g_xz[i]);
        iGzz = iDetG*(g_xx[i]*g_yy[i] - g_xy[i]*g_xy[i]);
        iGxy = iDetG*(g_xz[i]*g_yz[i] - g_zz[i]*g_xy[i]);
        iGxz = iDetG*(g_xy[i]*g_yz[i] - g_yy[i]*g_xz[i]);
        iGyz = iDetG*(g_xy[i]*g_xz[i] - g_xx[i]*g_yz[i]);

        // raise U_i to U^i
        geo_iux[i] = iGxx*geo_ux[i] + iGxy*geo_uy[i] + iGxz*geo_uz[i];
        geo_iuy[i] = iGxy*geo_ux[i] + iGyy*geo_uy[i] + iGyz*geo_uz[i];
        geo_iuz[i] = iGxz*geo_ux[i] + iGyz*geo_uy[i] + iGzz*geo_uz[i];

        // Calc U^0 (aka gamma) from the normalization condition
        if(CCTK_Equals(geo_type,"time")) { //U_mu*U^mu = -1
          geo_iut[i] = sqrt(geo_m[i]*geo_m[i] + iGxx*geo_ux[i]*geo_ux[i]
              + iGyy*geo_uy[i]*geo_uy[i] 
              + iGzz*geo_uz[i]*geo_uz[i] 
              + 2.0 *(iGxy*geo_ux[i]*geo_uy[i]
                + iGxz*geo_ux[i]*geo_uz[i] 
                + iGyz*geo_uy[i]*geo_uz[i]))
            /a[i];

        }
        if(CCTK_Equals(geo_type,"null")) { //U_mu*U^mu = 0
          geo_iut[i] = sqrt( iGxx*geo_ux[i]*geo_ux[i]
              + iGyy*geo_uy[i]*geo_uy[i]
              + iGzz*geo_uz[i]*geo_uz[i]
              + 2.0 *(iGxy*geo_ux[i]*geo_uy[i]
                + iGxz*geo_ux[i]*geo_uz[i]
                + iGyz*geo_uy[i]*geo_uz[i]))
            /a[i];
        } 

        // Calc velocity
        geo_x_rhs[i] = (1.0/geo_iut[i])*(iGxx*geo_ux[i] + iGxy*geo_uy[i] + iGxz*geo_uz[i]) - beta_x[i];
        geo_y_rhs[i] = (1.0/geo_iut[i])*(iGxy*geo_ux[i] + iGyy*geo_uy[i] + iGyz*geo_uz[i]) - beta_y[i];
        geo_z_rhs[i] = (1.0/geo_iut[i])*(iGxz*geo_ux[i] + iGyz*geo_uy[i] + iGzz*geo_uz[i]) - beta_z[i];

        // Calc acceleration
        // lapse contrib
        geo_ux_rhs[i] = -a[i] * geo_iut[i] * a_dx[i];
        geo_uy_rhs[i] = -a[i] * geo_iut[i] * a_dy[i];
        geo_uz_rhs[i] = -a[i] * geo_iut[i] * a_dz[i];

        // shift contrib 
        geo_ux_rhs[i] += beta_x_dx[i]*geo_ux[i] 
          + beta_y_dx[i]*geo_uy[i]
          + beta_z_dx[i]*geo_uz[i];
        geo_uy_rhs[i] += beta_x_dy[i]*geo_ux[i]
          + beta_y_dy[i]*geo_uy[i]
          + beta_z_dy[i]*geo_uz[i];
        geo_uz_rhs[i] += beta_x_dz[i]*geo_ux[i]
          + beta_y_dz[i]*geo_uy[i]
          + beta_z_dz[i]*geo_uz[i];

        // metric contrib
        geo_ux_rhs[i] += 0.5*(1.0/geo_iut[i]) 
          *(g_xx_dx[i]*geo_iux[i]*geo_iux[i]
              + g_yy_dx[i]*geo_iuy[i]*geo_iuy[i]
              + g_zz_dx[i]*geo_iuz[i]*geo_iuz[i]
              + 2.0 *(g_xy_dx[i]*geo_iux[i]*geo_iuy[i]
                + g_xz_dx[i]*geo_iux[i]*geo_iuz[i]
                + g_yz_dx[i]*geo_iuy[i]*geo_iuz[i]));

        geo_uy_rhs[i] += 0.5*(1.0/geo_iut[i])
          *(g_xx_dy[i]*geo_iux[i]*geo_iux[i]
              + g_yy_dy[i]*geo_iuy[i]*geo_iuy[i]
              + g_zz_dy[i]*geo_iuz[i]*geo_iuz[i]
              + 2.0 *(g_xy_dy[i]*geo_iux[i]*geo_iuy[i]
                + g_xz_dy[i]*geo_iux[i]*geo_iuz[i]
                + g_yz_dy[i]*geo_iuy[i]*geo_iuz[i]));

        geo_uz_rhs[i] += 0.5*(1.0/geo_iut[i])
          *(g_xx_dz[i]*geo_iux[i]*geo_iux[i]
              + g_yy_dz[i]*geo_iuy[i]*geo_iuy[i]
              + g_zz_dz[i]*geo_iuz[i]*geo_iuz[i]
              + 2.0 *(g_xy_dz[i]*geo_iux[i]*geo_iuy[i]
                + g_xz_dz[i]*geo_iux[i]*geo_iuz[i]
                + g_yz_dz[i]*geo_iuy[i]*geo_iuz[i]));


        // E&M contrib
        // Lower E^i, B^i to E_i, B_i
        b_x = g_xx[i]*bx[i] + g_xy[i]*by[i] + g_xz[i]*bz[i];
        b_y = g_xy[i]*bx[i] + g_yy[i]*by[i] + g_yz[i]*bz[i];
        b_z = g_xz[i]*bx[i] + g_yz[i]*by[i] + g_zz[i]*bz[i];
        e_x = g_xx[i]*ex[i] + g_xy[i]*ey[i] + g_xz[i]*ez[i];
        e_y = g_xy[i]*ex[i] + g_yy[i]*ey[i] + g_yz[i]*ez[i];
        e_z = g_xz[i]*ex[i] + g_yz[i]*ey[i] + g_zz[i]*ez[i];

        // Calc lorentz force 
        geo_ux_rhs[i] += geo_q[i]/geo_m[i]/geo_iut[i] 
          * (geo_iut[i]*e_x + geo_iuy[i]*b_z - geo_iuz[i]*b_y);
        geo_uy_rhs[i] += geo_q[i]/geo_m[i]/geo_iut[i] 
          * (geo_iut[i]*e_y + geo_iuz[i]*b_x - geo_iux[i]*b_z);
        geo_uz_rhs[i] += geo_q[i]/geo_m[i]/geo_iut[i] 
          * (geo_iut[i]*e_z + geo_iux[i]*b_y - geo_iuy[i]*b_x);

//        if(CCTK_IsThornActive("Pic")) {
//          geo_ux_rhs[i] += pic_ux_rhs[i];
//          geo_uy_rhs[i] += pic_uy_rhs[i];
//          geo_uz_rhs[i] += pic_uz_rhs[i];
//          //cout << pic_ux_rhs[i] << endl;
//        }

        if(CCTK_Equals(calc_sync,"yes")) {

          CCTK_REAL b_mag = sqrt(b_x*b_x +b_y*b_y + b_z*b_z);

          CCTK_REAL bhat_x = b_x/b_mag;
          CCTK_REAL bhat_y = b_y/b_mag;
          CCTK_REAL bhat_z = b_z/b_mag;

          CCTK_REAL vperp_x = geo_iuy[i]*bhat_z - geo_iuz[i]*bhat_y;
          CCTK_REAL vperp_y = geo_iuz[i]*bhat_x - geo_iux[i]*bhat_z;
          CCTK_REAL vperp_z = geo_iux[i]*bhat_y - geo_iuy[i]*bhat_x;

          CCTK_REAL vperp_mag = sqrt(vperp_x*vperp_x + vperp_y*vperp_y + vperp_z*vperp_z);

          CCTK_REAL r_0 = 2.81794*pow(10.0,-13); //in cgs
          geo_sync[i] = (2.0/3.0) * r_0*r_0 * vperp_mag*vperp_mag * geo_iut[i]*geo_iut[i] * b_mag*b_mag;
//          geo_sync[i] *= 3*pow(10.0,8) * 1.21 * pow(10.0,19);
          geo_sync[i] = log10(geo_sync[i]);

          CCTK_REAL v_mag = sqrt(geo_iux[i]*geo_iux[i] + geo_iuy[i]*geo_iuy[i] + geo_iuz[i]*geo_iuz[i]);
          CCTK_REAL sinpitchangle = vperp_mag/v_mag;
          if(isnan(sinpitchangle)) sinpitchangle = 1.0;

          geo_nu_c[i] = 2.0*3.141592653589793* 3 * geo_iut[i] * geo_iut[i] * geo_q[i] * b_mag * sinpitchangle / (2*geo_m[i]);
          geo_nu_c[i] = log10(abs(geo_nu_c[i]));
        }

        // Reverse time if requested
        if(CCTK_Equals(time_direction,"backward")) {
          geo_x_rhs[i] *= -1.0;
          geo_y_rhs[i] *= -1.0;
          geo_z_rhs[i] *= -1.0;

          geo_ux_rhs[i] *= -1.0;
          geo_uy_rhs[i] *= -1.0;
          geo_uz_rhs[i] *= -1.0;
        } //I assume this is all time reversal will take (TODO check that)

        // Store E_i and B_i if they are being output
        if(CCTK_Equals(em_output,"yes")) {
          geo_bx[i] = b_x;
          geo_by[i] = b_y;
          geo_bz[i] = b_z;
          geo_ex[i] = e_x;
          geo_ey[i] = e_y;
          geo_ez[i] = e_z;
        }

        /*        // Store current if its being output
                  if(CCTK_Equals(current_output,"yes")) {
                  geo_jx[i] = jx[i];
                  geo_jy[i] = jy[i];
                  geo_jz[i] = jz[i];

        // And b/c the current isnt valid near boundries...
        if(isnan(jx[i])) geo_jx[i] = 0.0;
        if(isnan(jy[i])) geo_jy[i] = 0.0;
        if(isnan(jz[i])) geo_jz[i] = 0.0;

        }
         */

        // Set dt/dt (its a test...)
        if(CCTK_Equals(geo_moltest,"yes"))  geo_t_rhs[i] = 1.0;

// TODO: geo_test isn't right
/*        if(CCTK_Equals(geo_test,"yes")) {

          //Test that p_i= G_ij * p^j 
          test1_x[i] = geo_ux[i] - (g_xx[i]*geo_iux[i] + g_xy[i]*geo_iuy[i] + g_xz[i]*geo_iuz[i]);
          test1_y[i] = geo_uy[i] - (g_xy[i]*geo_iux[i] + g_yy[i]*geo_iuy[i] + g_yz[i]*geo_iuz[i]);
          test1_z[i] = geo_uz[i] - (g_xz[i]*geo_iux[i] + g_yz[i]*geo_iuy[i] + g_zz[i]*geo_iuz[i]);

          //Test that p_i= G_ij * p^j 
          test2_x[i] = geo_iux[i] - (iGxx*geo_ux[i] + iGxy*geo_uy[i] + iGxz*geo_uz[i]);
          test2_y[i] = geo_iuy[i] - (iGxy*geo_ux[i] + iGyy*geo_uy[i] + iGyz*geo_uz[i]);
          test2_z[i] = geo_iuz[i] - (iGxz*geo_ux[i] + iGyz*geo_uy[i] + iGzz*geo_uz[i]);

          // Test that U_mu * U^mu = -1
          test3[i] = 1.0 + geo_ux[i]*geo_iux[i] + geo_uy[i]*geo_iuy[i] + geo_uz[i]*geo_iuz[i] // u_i*u^i 
            + geo_iut[i]*geo_iut[i] * (-a[i]*a[i] +    // (u^0)^2*g_00
                + beta_x[i] * (g_xx[i]*beta_x[i] + g_xy[i]*beta_y[i] + g_xz[i]*beta_z[i])
                + beta_y[i] * (g_xy[i]*beta_x[i] + g_yy[i]*beta_y[i] + g_yz[i]*beta_z[i])
                + beta_z[i] * (g_xz[i]*beta_x[i] + g_yz[i]*beta_y[i] + g_zz[i]*beta_z[i]))
            + 2.0*geo_iut[i] // 2*u^0*(u_i*beta^i)
            * (geo_ux[i]*beta_x[i] + geo_uy[i]*beta_y[i] +geo_uz[i]*beta_z[i]);

          // Output to stdout
          if( verbose > 1 ) {
            CCTK_VInfo(CCTK_THORNSTRING,
                "Geodesics_Test: dtau=%f r=%f test1 = %f,%f,%f   test2 = %f",
                (double)1.0/geo_iut[i],
                sqrt(geo_x[i]*geo_x[i]+geo_y[i]*geo_y[i]+geo_z[i]*geo_z[i]),
                (double)test1_x[i],
                (double)test1_y[i],
                (double)test1_z[i],
                (double)test3[i]); }
        }
*/
        // Some outputs (in case you want to flood stdout)
        if( verbose > 1 ) {
          CCTK_VInfo(CCTK_THORNSTRING,
              "Geodesics_CalcRHS(): Particle %d @ x=%f,y=%f,z=%f",
              i,
              (double)geo_x[i],
              (double)geo_y[i],
              (double)geo_z[i]); }
        if( verbose > 1 ) {
          CCTK_VInfo(CCTK_THORNSTRING,
              "Geodesics_CalcRHS(): Particle %d @ px=%f,py=%f,pz=%f",
              i,
              (double)geo_ux[i],
              (double)geo_uy[i],
              (double)geo_uz[i]); }
        if( verbose > 1 ) {
          CCTK_VInfo(CCTK_THORNSTRING,
              "Geodesics_CalcRHS(): Particle %d @ ax=%f,ay=%f,az=%f",
              i,
              (double)geo_ux_rhs[i],
              (double)geo_uy_rhs[i],
              (double)geo_uz_rhs[i]); }
      }
    }

    if( verbose > 2)
      CCTK_VInfo(CCTK_THORNSTRING,
          "Geodesics_CalcRHS(): finished");
  }
}
