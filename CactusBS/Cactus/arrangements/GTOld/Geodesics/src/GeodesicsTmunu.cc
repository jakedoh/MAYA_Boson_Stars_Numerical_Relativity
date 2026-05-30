#include "Geodesics.hh"
#include "carpet.hh"

extern "C" void Geodesics_Tmunu(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(start_t > cctk_time) return; // Save some cpu time if we aren't evolving yet...

  int lsh[1];
  int ierr = CCTK_GrouplshGN (cctkGH, 1, lsh, "Geodesics::geo_pos");
  int lnum_geo = lsh[0];

  if( verbose > 2)
    CCTK_VInfo(CCTK_THORNSTRING,
        "Pic_FindCell(): starting");

  int i,j,k;
  for (k=0; k<cctk_lsh[2]; k++) {
    for (j=0; j<cctk_lsh[1]; j++) {
      for (i=0; i<cctk_lsh[0]; i++) {
        int const index = CCTK_GFINDEX3D(cctkGH,i,j,k);
        geo_gi[index] = i;
        geo_gj[index] = j;
        geo_gk[index] = k;
        //geopercell[index]=0;
        //geo_rho[index] = 0.0;
        //if(cctk_iteration<=1)pic_phi[index] = 0.0;
        //          Mcoeff[index] = 0.0;
      }
    }
  }

  const CCTK_INT num_input_arrays = 4;
  const CCTK_INT num_output_arrays = 4;

  CCTK_STRING spi = "Geodesics::geo_gi";
  CCTK_STRING spj = "Geodesics::geo_gj";
  CCTK_STRING spk = "Geodesics::geo_gk";
  CCTK_STRING sdg = "ADMAnalysis::detofg";


  const CCTK_INT input_array_indices[num_input_arrays]
    = { CCTK_VarIndex(spi) , CCTK_VarIndex(spj), CCTK_VarIndex(spk), CCTK_VarIndex(sdg)};
  /*
     for(int i = 0 ; i < num_input_arrays ; i++)
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
  CCTK_REAL pi[lnum_geo], pj[lnum_geo], pk[lnum_geo], det_g[lnum_geo];

  void * output_arrays[num_output_arrays]
    = { (void *) pi, (void *) pj, (void *) pk, (void *) det_g };

  const CCTK_INT output_array_types[num_output_arrays]
    = { CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL };

  const CCTK_INT operand_indices[num_output_arrays]
    = { 0,1,2,3 };

  const CCTK_INT operation_codes[num_output_arrays]
    = { 0, 0, 0, 0 }; 

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
    ierr = Util_TableSetInt(param_table_handle, 2 /*interp_num_timelevels*/,
        "InterpNumTimelevels");
  }

  // Perform the interpolation from the global mesh (this is pretty important)
//      ierr += Util_TableSetInt(param_table_handle, 1,
//          "want_global_mode");

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

  for( int i=0; i < lnum_geo; i++ )  {
    //cout << pi[i] << endl;

    CCTK_INT ii = (int)pi[i];
    CCTK_INT ij = (int)pj[i];
    CCTK_INT ik = (int)pk[i];

    int gfind = CCTK_GFINDEX3D(cctkGH,ii,ij,ik);
    //geopercell[gfind]++;
    //cout << i << '\t' << pic_i[i] << '\t' << pic_j[i] << '\t' << pic_k[i] << endl;
//    cout << ii << '\t' << ij << '\t' << ik << endl;
//    cout << gfind << endl;


//    cout << geo_m[i] << endl;
//    cout << geo_iut[i] << endl;

if (gfind >= 0){
//    cout << geo_x[i] << endl;
//    cout << x[gfind] << endl;
//   cout << geo_m[i] << endl;
//    cout << geo_iut[i] << endl;

    CCTK_REAL dx = geo_x[i] - x[gfind];
    CCTK_REAL dy = geo_y[i] - y[gfind];
    CCTK_REAL dz = geo_z[i] - z[gfind];

    CCTK_REAL tx = CCTK_DELTA_SPACE(0) - dx;//1.0 - dx;
    CCTK_REAL ty = CCTK_DELTA_SPACE(1) - dy;//1.0 - dy;
    CCTK_REAL tz = CCTK_DELTA_SPACE(2) - dz;//1.0 - dz;
    // TODO this dont seem to be centered right...

    CCTK_REAL geo_nn = 1.0 / ( geo_iut[i] * sqrt((double)det_g[i])
                       *CCTK_DELTA_SPACE(0)*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(2));

    if (isnan(geo_nn)) geo_nn = 0.0;
 
    //cout << geo_nn  << endl;
    //        tt = geo_ut[i] * geo_ut[i];

    CCTK_INT ind000 = CCTK_GFINDEX3D(cctkGH,ii,ij,ik);
    CCTK_INT ind100 = CCTK_GFINDEX3D(cctkGH,ii+1,ij,ik);
    CCTK_INT ind010 = CCTK_GFINDEX3D(cctkGH,ii,ij+1,ik);
    CCTK_INT ind110 = CCTK_GFINDEX3D(cctkGH,ii+1,ij+1,ik);
    CCTK_INT ind001 = CCTK_GFINDEX3D(cctkGH,ii,ij,ik+1);
    CCTK_INT ind101 = CCTK_GFINDEX3D(cctkGH,ii+1,ij,ik+1);
    CCTK_INT ind011 = CCTK_GFINDEX3D(cctkGH,ii,ij+1,ik+1);
    CCTK_INT ind111 = CCTK_GFINDEX3D(cctkGH,ii+1,ij+1,ik+1);

    eTtt[ind000] += geo_m[i] * geo_nn * geo_iut[i] * geo_iut[i] * tx * ty * tz;
    eTtt[ind100] += geo_m[i] * geo_nn * geo_iut[i] * geo_iut[i] * dx * ty * tz;
    eTtt[ind010] += geo_m[i] * geo_nn * geo_iut[i] * geo_iut[i] * tx * dy * tz;
    eTtt[ind110] += geo_m[i] * geo_nn * geo_iut[i] * geo_iut[i] * dx * dy * tz;
    eTtt[ind001] += geo_m[i] * geo_nn * geo_iut[i] * geo_iut[i] * tx * ty * dz;
    eTtt[ind101] += geo_m[i] * geo_nn * geo_iut[i] * geo_iut[i] * dx * ty * dz;
    eTtt[ind011] += geo_m[i] * geo_nn * geo_iut[i] * geo_iut[i] * tx * dy * dz;
    eTtt[ind111] += geo_m[i] * geo_nn * geo_iut[i] * geo_iut[i] * dx * dy * dz;

    eTtx[ind000] += geo_m[i] * geo_nn * geo_iut[i] * geo_iux[i] * tx * ty * tz;
    eTtx[ind100] += geo_m[i] * geo_nn * geo_iut[i] * geo_iux[i] * dx * ty * tz;
    eTtx[ind010] += geo_m[i] * geo_nn * geo_iut[i] * geo_iux[i] * tx * dy * tz;
    eTtx[ind110] += geo_m[i] * geo_nn * geo_iut[i] * geo_iux[i] * dx * dy * tz;
    eTtx[ind001] += geo_m[i] * geo_nn * geo_iut[i] * geo_iux[i] * tx * ty * dz;
    eTtx[ind101] += geo_m[i] * geo_nn * geo_iut[i] * geo_iux[i] * dx * ty * dz;
    eTtx[ind011] += geo_m[i] * geo_nn * geo_iut[i] * geo_iux[i] * tx * dy * dz;
    eTtx[ind111] += geo_m[i] * geo_nn * geo_iut[i] * geo_iux[i] * dx * dy * dz;

    eTty[ind000] += geo_m[i] * geo_nn * geo_iut[i] * geo_iuy[i] * tx * ty * tz;
    eTty[ind100] += geo_m[i] * geo_nn * geo_iut[i] * geo_iuy[i] * dx * ty * tz;
    eTty[ind010] += geo_m[i] * geo_nn * geo_iut[i] * geo_iuy[i] * tx * dy * tz;
    eTty[ind110] += geo_m[i] * geo_nn * geo_iut[i] * geo_iuy[i] * dx * dy * tz;
    eTty[ind001] += geo_m[i] * geo_nn * geo_iut[i] * geo_iuy[i] * tx * ty * dz;
    eTty[ind101] += geo_m[i] * geo_nn * geo_iut[i] * geo_iuy[i] * dx * ty * dz;
    eTty[ind011] += geo_m[i] * geo_nn * geo_iut[i] * geo_iuy[i] * tx * dy * dz;
    eTty[ind111] += geo_m[i] * geo_nn * geo_iut[i] * geo_iuy[i] * dx * dy * dz;

    eTtz[ind000] += geo_m[i] * geo_nn * geo_iut[i] * geo_iuz[i] * tx * ty * tz;
    eTtz[ind100] += geo_m[i] * geo_nn * geo_iut[i] * geo_iuz[i] * dx * ty * tz;
    eTtz[ind010] += geo_m[i] * geo_nn * geo_iut[i] * geo_iuz[i] * tx * dy * tz;
    eTtz[ind110] += geo_m[i] * geo_nn * geo_iut[i] * geo_iuz[i] * dx * dy * tz;
    eTtz[ind001] += geo_m[i] * geo_nn * geo_iut[i] * geo_iuz[i] * tx * ty * dz;
    eTtz[ind101] += geo_m[i] * geo_nn * geo_iut[i] * geo_iuz[i] * dx * ty * dz;
    eTtz[ind011] += geo_m[i] * geo_nn * geo_iut[i] * geo_iuz[i] * tx * dy * dz;
    eTtz[ind111] += geo_m[i] * geo_nn * geo_iut[i] * geo_iuz[i] * dx * dy * dz;

    eTxx[ind000] += geo_m[i] * geo_nn * geo_iux[i] * geo_iux[i] * tx * ty * tz;
    eTxx[ind100] += geo_m[i] * geo_nn * geo_iux[i] * geo_iux[i] * dx * ty * tz;
    eTxx[ind010] += geo_m[i] * geo_nn * geo_iux[i] * geo_iux[i] * tx * dy * tz;
    eTxx[ind110] += geo_m[i] * geo_nn * geo_iux[i] * geo_iux[i] * dx * dy * tz;
    eTxx[ind001] += geo_m[i] * geo_nn * geo_iux[i] * geo_iux[i] * tx * ty * dz;
    eTxx[ind101] += geo_m[i] * geo_nn * geo_iux[i] * geo_iux[i] * dx * ty * dz;
    eTxx[ind011] += geo_m[i] * geo_nn * geo_iux[i] * geo_iux[i] * tx * dy * dz;
    eTxx[ind111] += geo_m[i] * geo_nn * geo_iux[i] * geo_iux[i] * dx * dy * dz;

    eTxy[ind000] += geo_m[i] * geo_nn * geo_iux[i] * geo_iuy[i] * tx * ty * tz;
    eTxy[ind100] += geo_m[i] * geo_nn * geo_iux[i] * geo_iuy[i] * dx * ty * tz;
    eTxy[ind010] += geo_m[i] * geo_nn * geo_iux[i] * geo_iuy[i] * tx * dy * tz;
    eTxy[ind110] += geo_m[i] * geo_nn * geo_iux[i] * geo_iuy[i] * dx * dy * tz;
    eTxy[ind001] += geo_m[i] * geo_nn * geo_iux[i] * geo_iuy[i] * tx * ty * dz;
    eTxy[ind101] += geo_m[i] * geo_nn * geo_iux[i] * geo_iuy[i] * dx * ty * dz;
    eTxy[ind011] += geo_m[i] * geo_nn * geo_iux[i] * geo_iuy[i] * tx * dy * dz;
    eTxy[ind111] += geo_m[i] * geo_nn * geo_iux[i] * geo_iuy[i] * dx * dy * dz;

    eTxz[ind000] += geo_m[i] * geo_nn * geo_iux[i] * geo_iuz[i] * tx * ty * tz;
    eTxz[ind100] += geo_m[i] * geo_nn * geo_iux[i] * geo_iuz[i] * dx * ty * tz;
    eTxz[ind010] += geo_m[i] * geo_nn * geo_iux[i] * geo_iuz[i] * tx * dy * tz;
    eTxz[ind110] += geo_m[i] * geo_nn * geo_iux[i] * geo_iuz[i] * dx * dy * tz;
    eTxz[ind001] += geo_m[i] * geo_nn * geo_iux[i] * geo_iuz[i] * tx * ty * dz;
    eTxz[ind101] += geo_m[i] * geo_nn * geo_iux[i] * geo_iuz[i] * dx * ty * dz;
    eTxz[ind011] += geo_m[i] * geo_nn * geo_iux[i] * geo_iuz[i] * tx * dy * dz;
    eTxz[ind111] += geo_m[i] * geo_nn * geo_iux[i] * geo_iuz[i] * dx * dy * dz;

    eTyy[ind000] += geo_m[i] * geo_nn * geo_iuy[i] * geo_iuy[i] * tx * ty * tz;
    eTyy[ind100] += geo_m[i] * geo_nn * geo_iuy[i] * geo_iuy[i] * dx * ty * tz;
    eTyy[ind010] += geo_m[i] * geo_nn * geo_iuy[i] * geo_iuy[i] * tx * dy * tz;
    eTyy[ind110] += geo_m[i] * geo_nn * geo_iuy[i] * geo_iuy[i] * dx * dy * tz;
    eTyy[ind001] += geo_m[i] * geo_nn * geo_iuy[i] * geo_iuy[i] * tx * ty * dz;
    eTyy[ind101] += geo_m[i] * geo_nn * geo_iuy[i] * geo_iuy[i] * dx * ty * dz;
    eTyy[ind011] += geo_m[i] * geo_nn * geo_iuy[i] * geo_iuy[i] * tx * dy * dz;
    eTyy[ind111] += geo_m[i] * geo_nn * geo_iuy[i] * geo_iuy[i] * dx * dy * dz;

    eTyz[ind000] += geo_m[i] * geo_nn * geo_iuy[i] * geo_iuz[i] * tx * ty * tz;
    eTyz[ind100] += geo_m[i] * geo_nn * geo_iuy[i] * geo_iuz[i] * dx * ty * tz;
    eTyz[ind010] += geo_m[i] * geo_nn * geo_iuy[i] * geo_iuz[i] * tx * dy * tz;
    eTyz[ind110] += geo_m[i] * geo_nn * geo_iuy[i] * geo_iuz[i] * dx * dy * tz;
    eTyz[ind001] += geo_m[i] * geo_nn * geo_iuy[i] * geo_iuz[i] * tx * ty * dz;
    eTyz[ind101] += geo_m[i] * geo_nn * geo_iuy[i] * geo_iuz[i] * dx * ty * dz;
    eTyz[ind011] += geo_m[i] * geo_nn * geo_iuy[i] * geo_iuz[i] * tx * dy * dz;
    eTyz[ind111] += geo_m[i] * geo_nn * geo_iuy[i] * geo_iuz[i] * dx * dy * dz;

    eTzz[ind000] += geo_m[i] * geo_nn * geo_iuz[i] * geo_iuz[i] * tx * ty * tz;
    eTzz[ind100] += geo_m[i] * geo_nn * geo_iuz[i] * geo_iuz[i] * dx * ty * tz;
    eTzz[ind010] += geo_m[i] * geo_nn * geo_iuz[i] * geo_iuz[i] * tx * dy * tz;
    eTzz[ind110] += geo_m[i] * geo_nn * geo_iuz[i] * geo_iuz[i] * dx * dy * tz;
    eTzz[ind001] += geo_m[i] * geo_nn * geo_iuz[i] * geo_iuz[i] * tx * ty * dz;
    eTzz[ind101] += geo_m[i] * geo_nn * geo_iuz[i] * geo_iuz[i] * dx * ty * dz;
    eTzz[ind011] += geo_m[i] * geo_nn * geo_iuz[i] * geo_iuz[i] * tx * dy * dz;
    eTzz[ind111] += geo_m[i] * geo_nn * geo_iuz[i] * geo_iuz[i] * dx * dy * dz;

    /*      if( verbose > 1 ) {
            CCTK_VInfo(CCTK_THORNSTRING,
            "Pic_CalcRHS(): Particle %d @ ax=%f,ay=%f,az=%f",
            i,
            (double)geo_ux_rhs[i],
            (double)geo_uy_rhs[i],
            (double)geo_uz_rhs[i]); }
     */      //}
  }
}
if( verbose > 2)
  CCTK_VInfo(CCTK_THORNSTRING,
      "Pic_FindCell(): finished");


}
