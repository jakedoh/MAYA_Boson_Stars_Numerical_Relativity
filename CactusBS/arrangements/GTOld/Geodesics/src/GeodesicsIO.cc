#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include "Geodesics.hh"

extern "C" void Geodesics_InitOutput(CCTK_ARGUMENTS);
extern "C" void Geodesics_Output(CCTK_ARGUMENTS);

static string output_filename(int i)
{
  // get the cactus output directory so that we don't spill data
  // everywhere.
  CCTK_STRING *out_dir  = (CCTK_STRING *) CCTK_ParameterGet("out_dir",
      "IOUtil", NULL);
  stringstream file_name;

  // This is the name of the output file
  file_name << "geosnap_" << i << ".okc" << ends;
//  file_name << "geo_time=" << i << ".okc" << ends;

  // Check to make sure you have a valid output directory
  // If not, just return the bare name and no path.
  if (*out_dir == 0) {
    CCTK_WARN(1,"Parameter IOUtil::out_dir not found");
    return file_name.str();
  }

  // If you do have a valid dir, the prepend the path to the
  // bare filename.
  return string(*out_dir) + string("/") + file_name.str();
}

extern "C" void Geodesics_InitOutput(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(start_t > cctk_time || cctk_iteration == 0) return; // Save some cpu time if we aren't evolving yet...

  if (CCTK_MyProc(cctkGH) == 0 && cctk_iteration%output_every == 0) { // Only write header from proc 0
   
    CCTK_REAL xmn, xmx, ymn, ymx, zmn, zmx;

    CCTK_CoordRange(cctkGH, &xmn, &xmx, -1, "x", "cart3d");
    CCTK_CoordRange(cctkGH, &ymn, &ymx, -1, "y", "cart3d");
    CCTK_CoordRange(cctkGH, &zmn, &zmx, -1, "z", "cart3d");

    ofstream outfile;
    outfile.open(output_filename((int)cctk_iteration).c_str(),
        ofstream::out | ofstream::app);      

    // Find number of vars for output on the point mesh
    int numout = 11;

    //if(CCTK_Equals(calc_gamma,"yes")) numout++; 
    if(CCTK_Equals(calc_sync,"yes")) numout+= 2;
    if(CCTK_Equals(em_output,"yes")) numout += 6;
    if(CCTK_Equals(acc_output,"yes")) numout += 3;
//    if(CCTK_Equals(current_output,"yes")) numout += 3;
//    if(CCTK_Equals(geo_test,"yes")) numout += 7;
    if(CCTK_Equals(geo_moltest,"yes")) numout++;

    // Output the first half of the XMDV header
    outfile << setprecision(io_precision);
    outfile << numout << " " << (int)num_geo << " 12" << endl
      << "x" << endl
      << "y" << endl
      << "z" << endl
      << "id" << endl
      << "u^t" << endl
      << "u^x" << endl
      << "u^y" << endl
      << "u^z" << endl
      << "u_x" << endl
      << "u_y" << endl
      << "u_z" << endl;
    if(CCTK_Equals(calc_sync,"yes")) { outfile 
      << "sync" << endl
      << "nu_c" << endl; }
    if(CCTK_Equals(em_output,"yes")) { outfile
      << "Bx" << endl
      << "By" << endl
      << "Bz" << endl
      << "Ex" << endl
      << "Ey" << endl
      << "Ez" << endl; }
    if(CCTK_Equals(acc_output,"yes")) { outfile
      << "acc_x" << endl
      << "acc_y" << endl
      << "acc_z" << endl; }
//    if(CCTK_Equals(current_output,"yes")) { outfile
//      << "Jx" << endl
//      << "Jy" << endl
//      << "Jz" << endl; }
//    if(CCTK_Equals(geo_test,"yes")) { outfile 
//      << "test1x" << endl
//      << "test1y" << endl
//      << "test1z" << endl
//      << "test2x" << endl
//      << "test2y" << endl
//      << "test2z" << endl
//      << "test3" << endl; }
    if(CCTK_Equals(geo_moltest,"yes")) outfile << "time" << endl;

    // Output the second half of the header
    // (none of this matters to visit, but is required to make the file a valid xmdv)
    //TODO: fix min (they are all y atm)
    outfile 
      << ymn << '\t' << xmx << '\t' << "10" << endl
      << ymn << '\t' << ymx << '\t' << "10" << endl
      << ymn << '\t' << zmx << '\t' << "10" << endl
      << 0 << '\t' << num_geo << '\t' << "10" << endl
      << -1.0 << '\t' << 1.0 << '\t' << "10" << endl
      << -1.0 << '\t' << 1.0 << '\t' << "10" << endl
      << -1.0 << '\t' << 1.0 << '\t' << "10" << endl
      << -1.0 << '\t' << 1.0 << '\t' << "10" << endl
      << -1.0 << '\t' << 1.0 << '\t' << "10" << endl
      << -1.0 << '\t' << 1.0 << '\t' << "10" << endl
      << -1.0 << '\t' << 1.0 << '\t' << "10" << endl;
    if(CCTK_Equals(calc_sync,"yes")) { outfile 
      << 0.0 << '\t' << 1.0 << '\t' << "10" << endl
      << 0.0 << '\t' << 1.0 << '\t' << "10" << endl; }
    if(CCTK_Equals(em_output,"yes")) { outfile
      << -1.0 << '\t' << 1.0 << '\t' << "10" << endl
      << -1.0 << '\t' << 1.0 << '\t' << "10" << endl
      << -1.0 << '\t' << 1.0 << '\t' << "10" << endl
      << -1.0 << '\t' << 1.0 << '\t' << "10" << endl
      << -1.0 << '\t' << 1.0 << '\t' << "10" << endl
      << -1.0 << '\t' << 1.0 << '\t' << "10" << endl; }
    if(CCTK_Equals(acc_output,"yes")) { outfile
      << -1.0 << '\t' << 1.0 << '\t' << "10" << endl
      << -1.0 << '\t' << 1.0 << '\t' << "10" << endl
      << -1.0 << '\t' << 1.0 << '\t' << "10" << endl; }
//    if(CCTK_Equals(current_output,"yes")) { outfile
//      << -1.0 << '\t' << 1.0 << '\t' << "10" << endl
//      << -1.0 << '\t' << 1.0 << '\t' << "10" << endl
//      << -1.0 << '\t' << 1.0 << '\t' << "10" << endl; }
//    if(CCTK_Equals(geo_test,"yes")) { outfile 
//      << -10.0 << '\t' << 10.0 << '\t' << "10" << endl
//      << -10.0 << '\t' << 10.0 << '\t' << "10" << endl
//      << -10.0 << '\t' << 10.0 << '\t' << "10" << endl
//      << -10.0 << '\t' << 10.0 << '\t' << "10" << endl
//      << -10.0 << '\t' << 10.0 << '\t' << "10" << endl
//      << -10.0 << '\t' << 10.0 << '\t' << "10" << endl
//      << -10.0 << '\t' << 10.0 << '\t' << "10" << endl; }
    if(CCTK_Equals(geo_moltest,"yes"))  outfile << 0.0 << '\t' << 1.0 << '\t' << "1000000" << endl;

    outfile.close();
  }
  CCTK_Barrier(cctkGH); // Wait till the other procs finish (to preserve order in the file)
}

extern "C" void Geodesics_Output(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(start_t > cctk_time || cctk_iteration == 0) return; // Save some cpu time if we aren't evolving yet...

  int nprocs = CCTK_nProcs(cctkGH);

  for(int proc=0; proc < nprocs; proc++) // Loop over procs in a known way to preserve order in the file
  {
    if(proc==CCTK_MyProc(cctkGH) && cctk_iteration%output_every == 0){

      int lsh[1];
      int const ierr = CCTK_GrouplshGN (cctkGH, 1, lsh, "Geodesics::geo_pos");
      int lnum_geo = lsh[0];

      ofstream outfile;
      outfile.open(output_filename((int)cctk_iteration).c_str(),
          ofstream::out | ofstream::app);

      outfile << setprecision(io_precision);

      for( int i=0; i < lnum_geo; i++ ) {
        //Dont output inactive geodesics (visit will stop reading at a nan)
        if( active_geo[i] != -1){ outfile 
            << geo_x[i] << '\t'
            << geo_y[i] << '\t'
            << geo_z[i] << '\t'
            << i+lnum_geo*proc << '\t' //this should give the particles constant id's in [0,num_geo]
            << geo_iut[i] << '\t'
            << geo_iux[i] << '\t'
            << geo_iuy[i] << '\t'
            << geo_iuz[i] << '\t'
            << geo_ux[i] << '\t'
            << geo_uy[i] << '\t'
            << geo_uz[i] << '\t';
          if(CCTK_Equals(calc_sync,"yes")) { outfile 
            << geo_sync[i] << '\t'
            << geo_nu_c[i] << '\t'; }
          if(CCTK_Equals(em_output,"yes")) { outfile
            << geo_bx[i] << '\t'
              << geo_by[i] << '\t'
              << geo_bz[i] << '\t'
              << geo_ex[i] << '\t'
              << geo_ey[i] << '\t'
              << geo_ez[i] << '\t'; }
          if(CCTK_Equals(acc_output,"yes")) { outfile
            << geo_ux_rhs[i] << '\t'
              << geo_uy_rhs[i] << '\t'
              << geo_uz_rhs[i] << '\t'; }
 //         if(CCTK_Equals(current_output,"yes")) { outfile
 //           << geo_jx[i] << '\t'
//              << geo_jy[i] << '\t'
//              << geo_jz[i] << '\t'; }
//          if(CCTK_Equals(geo_test,"yes")) { outfile 
//            << test1_x[i] << '\t'
//              << test1_y[i] << '\t'
//              << test1_z[i] << '\t'
//              << test2_x[i] << '\t'
//              << test2_y[i] << '\t'
//              << test2_z[i] << '\t'
//              << test3[i] << '\t'; }
          if(CCTK_Equals(geo_moltest,"yes")) outfile << geo_t[i] << '\t';

          outfile << endl;
        } else
        { outfile //if the geodesic isnt active we still want its id in the output (for visit)
            << 0.0 << '\t'
            << 0.0 << '\t'
            << 0.0 << '\t'
            << i+lnum_geo*proc << '\t' //this should give the particles constant id's in [0,num_geo]
            << 1.0 << '\t'
            << 0.0 << '\t'
            << 0.0 << '\t'
            << 0.0 << '\t'
            << 0.0 << '\t'
            << 0.0 << '\t'
            << 0.0 << '\t';
          if(CCTK_Equals(calc_sync,"yes")) { outfile 
            << 1.0 << '\t'
            << 1.0 << '\t'; }
          if(CCTK_Equals(em_output,"yes")) { outfile
            << 0.0 << '\t'
              << 0.0 << '\t'
              << 0.0 << '\t'
              << 0.0 << '\t'
              << 0.0 << '\t'
              << 0.0 << '\t'; }
          if(CCTK_Equals(em_output,"yes")) { outfile
            << 0.0 << '\t'
              << 0.0 << '\t'
              << 0.0 << '\t'; }
//          if(CCTK_Equals(current_output,"yes")) { outfile
//            << 0.0 << '\t'
//              << 0.0 << '\t'
//              << 0.0 << '\t'; }
//          if(CCTK_Equals(geo_test,"yes")) { outfile
//            << 0.0 << '\t'
//              << 0.0 << '\t'
//              << 0.0 << '\t'
//              << 0.0 << '\t'
//              << 0.0 << '\t'
//              << 0.0 << '\t'
//              << 0.0 << '\t'; }
          if(CCTK_Equals(geo_moltest,"yes")) outfile << 0.0 << '\t';

          outfile << endl;
        }
      }
      outfile.close();
    }
    CCTK_Barrier(cctkGH); // Wait till the other procs finish (to preserve order in the file)
  }
}
