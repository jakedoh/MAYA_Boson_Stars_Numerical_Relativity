#include <stdio.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "mpi.h"

#include "vars.hh"


using namespace std;

extern "C" void SphericalHarmonicReconGeneric_Startup(CCTK_ARGUMENTS)
{
   DECLARE_CCTK_ARGUMENTS;
   DECLARE_CCTK_PARAMETERS;
   
   assert(SHR::initialized == 0);
   
   SHR::initialized = 1;
   
   if (verbose) {
      stringstream str;
      str << "char. dt = " << cctk_delta_time << endl;
      CCTK_INFO(str.str().c_str());
   }
   
   if (CCTK_IsFunctionAliased("GetMPICommWorld"))
    {
      SHR::comm_world = *(MPI_Comm *)GetMPICommWorld(cctkGH);
    }
    else
    {
      SHR::comm_world = MPI_COMM_WORLD;
    }
   
   SHR::time_derivative_in_file = time_derivative_in_file;
   if (SHR::time_derivative_in_file)
      CCTK_INFO("Time derivatives are set to be contained in the file(s).");
   else
      CCTK_INFO("Taking time derivatives via finite differences from coefficients.");
   
   if (time_interpolate) {
      ostringstream s;
      s << "Using 5th-order time interpolation to interpolate to constant du=" << cctk_delta_time << " if necessary.";
      CCTK_INFO(s.str().c_str());
   }
   
   if (CCTK_Equals(sYlm_convention, "Condon-Shortley-Phase"))
   {
      SHR::use_Condon_Shortley_phase_factor = true;
      CCTK_INFO("Using Condon-Shortley phase convention for sYlm (factor of (-1)^m) for input worldtube modes.");
   }
   else
      CCTK_INFO("Using Goldberg67 convention of sYlm for input worldtube modes.");
   
   // open database, i.e. scan input file containing the harmonic coefficients
   CCTK_INFO("Opening and scanning file containing harmonic coefficients for the Cauchy variables.");
   CCTK_VInfo(CCTK_THORNSTRING,"Trying to open ASCII files in \"%s\".", path);
   
   // get all different databases-names
   vector<string> fnames;
   if (!CCTK_EQUALS(file, "")) {
      // we have a single file!
      fnames.resize(30, file);
   } else {
      // we have multiple files!
      fnames.push_back(file_gxx[0]);
      fnames.push_back(file_gxy[0]);
      fnames.push_back(file_gxz[0]);
      fnames.push_back(file_gyy[0]);
      fnames.push_back(file_gyz[0]);
      fnames.push_back(file_gzz[0]);
      fnames.push_back(file_shiftx[0]);
      fnames.push_back(file_shifty[0]);
      fnames.push_back(file_shiftz[0]);
      fnames.push_back(file_lapse[0]);
      fnames.push_back(file_gxx[1]);
      fnames.push_back(file_gxy[1]);
      fnames.push_back(file_gxz[1]);
      fnames.push_back(file_gyy[1]);
      fnames.push_back(file_gyz[1]);
      fnames.push_back(file_gzz[1]);
      fnames.push_back(file_shiftx[1]);
      fnames.push_back(file_shifty[1]);
      fnames.push_back(file_shiftz[1]);
      fnames.push_back(file_lapse[1]);
      fnames.push_back(file_gxx[2]);
      fnames.push_back(file_gxy[2]);
      fnames.push_back(file_gxz[2]);
      fnames.push_back(file_gyy[2]);
      fnames.push_back(file_gyz[2]);
      fnames.push_back(file_gzz[2]);
      fnames.push_back(file_shiftx[2]);
      fnames.push_back(file_shifty[2]);
      fnames.push_back(file_shiftz[2]);
      fnames.push_back(file_lapse[2]);
   }
   
   
   for (int i=0; i < 3*NUM_METRIC_COMPONENTS; ++i)
   {
      int found = -1;
      for (int j=0; j < i; ++j)
      {
         if (fnames[i] == fnames[j]) 
         {
            found = j;
            break;
         }
      }
      if (found < 0)
      {
         string fullname = path;
         if (!fullname.empty() && fullname[fullname.size()-1] != '/') {
            // an empty string means that the filenames contain full paths
            fullname.append("/");
         }
         fullname.append(fnames[i]);

         CCTK_VInfo(CCTK_THORNSTRING, "Parsing %s", fullname.c_str());
         if (CCTK_EQUALS(format, "ASCII"))
            SHR::db[i] = new SHR::SPH_db_ASCII(fullname, verbose, cached_timesteps+1, column_time, column_iteration, column_radius);
         if (CCTK_EQUALS(format, "DAT"))
            SHR::db[i] = new SHR::SPH_db_DAT(fullname, verbose, cached_timesteps+1, lmaxInFile, column_time, column_iteration, column_radius, column_lmax, column_n_variables, column_data, false);
         if (CCTK_EQUALS(format, "DAT-v2"))
            SHR::db[i] = new SHR::SPH_db_DAT(fullname, verbose, cached_timesteps+1, lmaxInFile, column_time, column_iteration, column_radius, column_lmax, column_n_variables, column_data, true);
         if (CCTK_EQUALS(format, "SpEC-H5"))
            SHR::db[i] = new SHR::SPH_db_SpEC_H5(fullname, verbose, cached_timesteps+1);
         if (CCTK_EQUALS(format, "SpEC-H5-v2"))
            SHR::db[i] = new SHR::SPH_db_SpEC_H5_v2(fullname, verbose, cached_timesteps+1);
         
         // check if timestep-size found in file matches characteristic timestep!
         if (fabs(SHR::db[i]->delta_t()-CCTK_DELTA_TIME) > 1e-8)
         {
            if (CCTK_DELTA_TIME > SHR::db[i]->delta_t())
	       CCTK_INFO("Time step-size of harmonic amplitudes of Cauchy variables is greater than characteristic time step-size!");
	    else
	       CCTK_INFO("Time step-size of harmonic amplitudes of Cauchy variables is smaller than characteristic time step-size!");
	    
	    {
	       stringstream str;
	       str << "  (intitial) Cauchy dt  = " << SHR::db[i]->delta_t() << endl;
	       CCTK_INFO(str.str().c_str());
            }
            {
               stringstream str;
               str << "             char.  dt  = " << CCTK_DELTA_TIME << endl;
               CCTK_INFO(str.str().c_str());
	    }
	    
	    if (!time_interpolate)
	       CCTK_WARN(0, "Time interpolation has not been requested! Stopping.");
         }
         else
         {
            CCTK_INFO("Initial Cauchy dt found in file and char. dt are identical.");
         }
         
         // check if specified radius matches selected radius in file
         if (column_radius >= 0)
            if (fabs(cr - SHR::db[i]->radius(sphere_number)) > 1e-8)
               CCTK_WARN(0, "Radius of selected sphere and sphere number do not match!");
            
         // check if lmax matches lmax found in file (DAT format only)
         
         
         // check if n_variables matches n_variables found in file (DAT format only)
         
         
      }
      else
         SHR::db[i] = SHR::db[found];
    
   }

   // set up time derivative business.
   int n_timelevels = 1; 
   int timelevel = 0; 
   int t_fd_order = 0;
   if (!time_derivative_in_file) {
      if (time_fd_order == 2) {
	 n_timelevels = 3; // second-order time-derivative needs 3 timelevels
	 timelevel = 1;    // and we operate on the second timelevel (counting from 0)
	 t_fd_order = 2;
      }
      if (time_fd_order == 4) {
	 n_timelevels = 5; // fourth-order time-derivative needs 5 timelevels
	 timelevel = 2;    // and we operate on the third timelevel (counting from 0)
	 t_fd_order = 4;
      }
   }
   
   // set up interpolation business
   int n_timelevels_for_coeffs_in_file = n_timelevels;
   
   // get initial "nominal" time in files.
   // If we are not using time interpolation and don't need to take time derivatives
   // then this is the very first time step found in the boundary data file
   if (time_interpolate)
   {
      n_timelevels_for_coeffs_in_file = 6;  // 6 for 5th-order interpolation
      if(disable_auto_time_offset==0)
        SHR::initial_time_in_file = SHR::db[0]->time(timelevel+2);
      ostringstream s;
      s << "Initial nominal time found in file is t_0 = " << SHR::initial_time_in_file << " corresponding to timestep n = " << timelevel+2 << " (counting from 0).";
      CCTK_INFO(s.str().c_str());
   }
   else
   {
      if(disable_auto_time_offset==0)
        SHR::initial_time_in_file = SHR::db[0]->time(timelevel);
      ostringstream s;
      s << "Initial nominal time found in file is t_0 = " << SHR::initial_time_in_file << " corresponding to timestep n = " << timelevel << " (counting from 0).";
      CCTK_INFO(s.str().c_str());
   }
   
   // initialize decomposed variables
   SHR::C[0] = new SHR::spherical_decomposed_variable<CCTK_REAL>(cctkGH, n_timelevels_for_coeffs_in_file, n_timelevels, timelevel, t_fd_order, CCTK_DELTA_TIME, column_gxx[0], sphere_number, *SHR::db[0], SHR::comm_world);
   SHR::C[1] = new SHR::spherical_decomposed_variable<CCTK_REAL>(cctkGH, n_timelevels_for_coeffs_in_file, n_timelevels, timelevel, t_fd_order, CCTK_DELTA_TIME, column_gxy[0], sphere_number, *SHR::db[1], SHR::comm_world);
   SHR::C[2] = new SHR::spherical_decomposed_variable<CCTK_REAL>(cctkGH, n_timelevels_for_coeffs_in_file, n_timelevels, timelevel, t_fd_order, CCTK_DELTA_TIME, column_gxz[0], sphere_number, *SHR::db[2], SHR::comm_world);
   SHR::C[3] = new SHR::spherical_decomposed_variable<CCTK_REAL>(cctkGH, n_timelevels_for_coeffs_in_file, n_timelevels, timelevel, t_fd_order, CCTK_DELTA_TIME, column_gyy[0], sphere_number, *SHR::db[3], SHR::comm_world);
   SHR::C[4] = new SHR::spherical_decomposed_variable<CCTK_REAL>(cctkGH, n_timelevels_for_coeffs_in_file, n_timelevels, timelevel, t_fd_order, CCTK_DELTA_TIME, column_gyz[0], sphere_number, *SHR::db[4], SHR::comm_world);
   SHR::C[5] = new SHR::spherical_decomposed_variable<CCTK_REAL>(cctkGH, n_timelevels_for_coeffs_in_file, n_timelevels, timelevel, t_fd_order, CCTK_DELTA_TIME, column_gzz[0], sphere_number, *SHR::db[5], SHR::comm_world);
   SHR::C[6] = new SHR::spherical_decomposed_variable<CCTK_REAL>(cctkGH, n_timelevels_for_coeffs_in_file, n_timelevels, timelevel, t_fd_order, CCTK_DELTA_TIME, column_shiftx[0], sphere_number, *SHR::db[6], SHR::comm_world);
   SHR::C[7] = new SHR::spherical_decomposed_variable<CCTK_REAL>(cctkGH, n_timelevels_for_coeffs_in_file, n_timelevels, timelevel, t_fd_order, CCTK_DELTA_TIME, column_shifty[0], sphere_number, *SHR::db[7], SHR::comm_world);
   SHR::C[8] = new SHR::spherical_decomposed_variable<CCTK_REAL>(cctkGH, n_timelevels_for_coeffs_in_file, n_timelevels, timelevel, t_fd_order, CCTK_DELTA_TIME, column_shiftz[0], sphere_number, *SHR::db[8], SHR::comm_world);
   SHR::C[9] = new SHR::spherical_decomposed_variable<CCTK_REAL>(cctkGH, n_timelevels_for_coeffs_in_file, n_timelevels, timelevel, t_fd_order, CCTK_DELTA_TIME, column_lapse[0], sphere_number, *SHR::db[9], SHR::comm_world);
   
   SHR::Cr[0] = new SHR::spherical_decomposed_variable<CCTK_REAL>(cctkGH, n_timelevels_for_coeffs_in_file, n_timelevels, timelevel, t_fd_order, CCTK_DELTA_TIME, column_gxx[1], sphere_number, *SHR::db[10], SHR::comm_world);
   SHR::Cr[1] = new SHR::spherical_decomposed_variable<CCTK_REAL>(cctkGH, n_timelevels_for_coeffs_in_file, n_timelevels, timelevel, t_fd_order, CCTK_DELTA_TIME, column_gxy[1], sphere_number, *SHR::db[11], SHR::comm_world);
   SHR::Cr[2] = new SHR::spherical_decomposed_variable<CCTK_REAL>(cctkGH, n_timelevels_for_coeffs_in_file, n_timelevels, timelevel, t_fd_order, CCTK_DELTA_TIME, column_gxz[1], sphere_number, *SHR::db[12], SHR::comm_world);
   SHR::Cr[3] = new SHR::spherical_decomposed_variable<CCTK_REAL>(cctkGH, n_timelevels_for_coeffs_in_file, n_timelevels, timelevel, t_fd_order, CCTK_DELTA_TIME, column_gyy[1], sphere_number, *SHR::db[13], SHR::comm_world);
   SHR::Cr[4] = new SHR::spherical_decomposed_variable<CCTK_REAL>(cctkGH, n_timelevels_for_coeffs_in_file, n_timelevels, timelevel, t_fd_order, CCTK_DELTA_TIME, column_gyz[1], sphere_number, *SHR::db[14], SHR::comm_world);
   SHR::Cr[5] = new SHR::spherical_decomposed_variable<CCTK_REAL>(cctkGH, n_timelevels_for_coeffs_in_file, n_timelevels, timelevel, t_fd_order, CCTK_DELTA_TIME, column_gzz[1], sphere_number, *SHR::db[15], SHR::comm_world);
   SHR::Cr[6] = new SHR::spherical_decomposed_variable<CCTK_REAL>(cctkGH, n_timelevels_for_coeffs_in_file, n_timelevels, timelevel, t_fd_order, CCTK_DELTA_TIME, column_shiftx[1], sphere_number, *SHR::db[16], SHR::comm_world);
   SHR::Cr[7] = new SHR::spherical_decomposed_variable<CCTK_REAL>(cctkGH, n_timelevels_for_coeffs_in_file, n_timelevels, timelevel, t_fd_order, CCTK_DELTA_TIME, column_shifty[1], sphere_number, *SHR::db[17], SHR::comm_world);
   SHR::Cr[8] = new SHR::spherical_decomposed_variable<CCTK_REAL>(cctkGH, n_timelevels_for_coeffs_in_file, n_timelevels, timelevel, t_fd_order, CCTK_DELTA_TIME, column_shiftz[1], sphere_number, *SHR::db[18], SHR::comm_world);
   SHR::Cr[9] = new SHR::spherical_decomposed_variable<CCTK_REAL>(cctkGH, n_timelevels_for_coeffs_in_file, n_timelevels, timelevel, t_fd_order, CCTK_DELTA_TIME, column_lapse[1], sphere_number, *SHR::db[19], SHR::comm_world);
   
   if (time_derivative_in_file) {
      SHR::Ct[0] = new SHR::spherical_decomposed_variable<CCTK_REAL>(cctkGH, n_timelevels_for_coeffs_in_file, n_timelevels, timelevel, t_fd_order, CCTK_DELTA_TIME, column_gxx[2], sphere_number, *SHR::db[20], SHR::comm_world);
      SHR::Ct[1] = new SHR::spherical_decomposed_variable<CCTK_REAL>(cctkGH, n_timelevels_for_coeffs_in_file, n_timelevels, timelevel, t_fd_order, CCTK_DELTA_TIME, column_gxy[2], sphere_number, *SHR::db[21], SHR::comm_world);
      SHR::Ct[2] = new SHR::spherical_decomposed_variable<CCTK_REAL>(cctkGH, n_timelevels_for_coeffs_in_file, n_timelevels, timelevel, t_fd_order, CCTK_DELTA_TIME, column_gxz[2], sphere_number, *SHR::db[22], SHR::comm_world);
      SHR::Ct[3] = new SHR::spherical_decomposed_variable<CCTK_REAL>(cctkGH, n_timelevels_for_coeffs_in_file, n_timelevels, timelevel, t_fd_order, CCTK_DELTA_TIME, column_gyy[2], sphere_number, *SHR::db[23], SHR::comm_world);
      SHR::Ct[4] = new SHR::spherical_decomposed_variable<CCTK_REAL>(cctkGH, n_timelevels_for_coeffs_in_file, n_timelevels, timelevel, t_fd_order, CCTK_DELTA_TIME, column_gyz[2], sphere_number, *SHR::db[24], SHR::comm_world);
      SHR::Ct[5] = new SHR::spherical_decomposed_variable<CCTK_REAL>(cctkGH, n_timelevels_for_coeffs_in_file, n_timelevels, timelevel, t_fd_order, CCTK_DELTA_TIME, column_gzz[2], sphere_number, *SHR::db[25], SHR::comm_world);
      SHR::Ct[6] = new SHR::spherical_decomposed_variable<CCTK_REAL>(cctkGH, n_timelevels_for_coeffs_in_file, n_timelevels, timelevel, t_fd_order, CCTK_DELTA_TIME, column_shiftx[2], sphere_number, *SHR::db[26], SHR::comm_world);
      SHR::Ct[7] = new SHR::spherical_decomposed_variable<CCTK_REAL>(cctkGH, n_timelevels_for_coeffs_in_file, n_timelevels, timelevel, t_fd_order, CCTK_DELTA_TIME, column_shifty[2], sphere_number, *SHR::db[27], SHR::comm_world);
      SHR::Ct[8] = new SHR::spherical_decomposed_variable<CCTK_REAL>(cctkGH, n_timelevels_for_coeffs_in_file, n_timelevels, timelevel, t_fd_order, CCTK_DELTA_TIME, column_shiftz[2], sphere_number, *SHR::db[28], SHR::comm_world);
      SHR::Ct[9] = new SHR::spherical_decomposed_variable<CCTK_REAL>(cctkGH, n_timelevels_for_coeffs_in_file, n_timelevels, timelevel, t_fd_order, CCTK_DELTA_TIME, column_lapse[2], sphere_number, *SHR::db[29], SHR::comm_world);
   }
   
   
   if (initial_time > 0)
   {
      // set starting iteration to current iteration
      // which in case of initial_time!=0 is non-zero!
      SHR::starting_iteration = SHR::C[0]->get_iteration(initial_time+SHR::initial_time_in_file);
      
      ostringstream s;
      s << "Neglecting first t = " << initial_time << " from worldtube file." << endl;
      s << "Starting evolution from time " << initial_time+SHR::initial_time_in_file << " corresponding to timestep n = " << SHR::starting_iteration << " (counting from 0).";
      CCTK_INFO(s.str().c_str());
      if (SHR::starting_iteration < 0) 
         CCTK_WARN(0, "The initial time is not contained in the file!");
   }
   
   // fill all timelevels (only if we do not have time-interpolation)
   // after this step we will have read in timestep 0,1 from input file.
   // this means we will always read two timestep ahead (timelevel==0) because we need it for time-derivatives!
   if (!time_derivative_in_file && !time_interpolate)
   {
      CCTK_INFO("Loading decomposed variables for first timestep...");
      
      if (time_fd_order == 2)
      {
	 for (int i=0; i < 3*NUM_METRIC_COMPONENTS; ++i) {
	    if (SHR::db[i]->n_timesteps() < n_timelevels-1)
	       CCTK_WARN(0, "Not enough timesteps in Cauchy harmonic coefficient file(s) to initialize first characteristic timestep!");
	 }
         for (int i=0; i < NUM_METRIC_COMPONENTS; ++i) {
            SHR::C[i]->set_modes_from_file(1, cctk_iteration+SHR::starting_iteration);   // (timelevel, iteration), timelevel=0 is current timelevel, timelevl=1 is one timelevel in the past etc...
	    SHR::Cr[i]->set_modes_from_file(1, cctk_iteration+SHR::starting_iteration);   // (timelevel, iteration), timelevel=0 is current timelevel, timelevl=1 is one timelevel in the past etc...
	 }
         
         for (int i=0; i < NUM_METRIC_COMPONENTS; ++i) {
            SHR::C[i]->set_modes_from_file(0, cctk_iteration+SHR::starting_iteration+1);
	    SHR::Cr[i]->set_modes_from_file(0, cctk_iteration+SHR::starting_iteration+1);
	 }
      }
      
      if (time_fd_order == 4)
      {
	 for (int i=0; i < 3*NUM_METRIC_COMPONENTS; ++i) {
	    if (SHR::db[i]->n_timesteps() < n_timelevels-3)
	       CCTK_WARN(0, "Not enough timesteps in Cauchy harmonic coefficient file to initialize first characteristic timestep!");
	 }
	 
         for (int i=0; i < NUM_METRIC_COMPONENTS; ++i) {
            SHR::C[i]->set_modes_from_file(3, cctk_iteration+SHR::starting_iteration);
	    SHR::Cr[i]->set_modes_from_file(3, cctk_iteration+SHR::starting_iteration);
	 }
	 
         for (int i=0; i < NUM_METRIC_COMPONENTS; ++i) {
            SHR::C[i]->set_modes_from_file(2, cctk_iteration+SHR::starting_iteration+1);
	    SHR::Cr[i]->set_modes_from_file(2, cctk_iteration+SHR::starting_iteration+1);
	 }
      
         for (int i=0; i < NUM_METRIC_COMPONENTS; ++i) {
            SHR::C[i]->set_modes_from_file(1, cctk_iteration+SHR::starting_iteration+2);
	    SHR::Cr[i]->set_modes_from_file(1, cctk_iteration+SHR::starting_iteration+2);
	 }
	 
         for (int i=0; i < NUM_METRIC_COMPONENTS; ++i) {
            SHR::C[i]->set_modes_from_file(0, cctk_iteration+SHR::starting_iteration+3);
	    SHR::Cr[i]->set_modes_from_file(0, cctk_iteration+SHR::starting_iteration+3);
	 }
      }
   }
   
   // set iteration number in file to current iteration
   // (which in case of recovery and/or initial_time!=0 is non-zero)
   SHR::iteration_in_file = cctk_iteration;
   
}






extern "C" void SphericalHarmonicReconGeneric_Shutdown(CCTK_ARGUMENTS)
{
   DECLARE_CCTK_ARGUMENTS;
   DECLARE_CCTK_PARAMETERS;

   // delete all variables
   for (int i=0; i < NUM_METRIC_COMPONENTS; ++i)
   {
      if (SHR::C[i]) delete SHR::C[i];
      if (SHR::Cr[i]) delete SHR::Cr[i];
      if (time_derivative_in_file && SHR::Ct[i]) delete SHR::Ct[i];
   }

   // there are multiple pointers pointing to the same class,
   // make sure that each class is deleted only once

   for(size_t i=0;i<SHR::db.size();++i) {
     if(SHR::db[i]==0) continue;
     for(size_t j=i+1;j<SHR::db.size();++j)
       if(SHR::db[j]==SHR::db[i]) {
         SHR::db[j]=0;
       }
     delete SHR::db[i]; SHR::db[i]=0;
   }
}


