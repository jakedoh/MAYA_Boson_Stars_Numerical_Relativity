#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "ShiftTracker.hh"

// ShiftTracker_Header():
//   Creates the output files, truncates them if they already exists
//   and outputs the header into them.  It creates one output file
//   per tracker.  They are named ShiftTracker[i].asc.
extern "C" void ShiftTracker_Header(CCTK_ARGUMENTS);

// ShiftTracker_Output():
//   Appends data to the already created output files.  Output data are
//   itteration number, time, and x, y, z coordinates for each tracker.
extern "C" void ShiftTracker_Output(CCTK_ARGUMENTS);

static string output_filename(int i)
{
  // This gets the cactus output directory so that we don't spill data
  // everywhere.
  CCTK_STRING *out_dir  = (CCTK_STRING *) CCTK_ParameterGet("out_dir", 
      "IOUtil", NULL);
  stringstream file_name;

  // This is the name of the output file
  file_name << "ShiftTracker";

  // -1 here means that you will not be appending a number
  // to the end of the filename.  This implies that 
  // num_trackers < output_threshold.  I think there might 
  // be a very clean way to implement this function such that
  // this if/else output statements can be combined.
  if (i != -1) {
    file_name << i;
  }

  file_name << ".asc" << ends;
                                            
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

extern "C" void ShiftTracker_Header(CCTK_ARGUMENTS) 
{ 
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (CCTK_MyProc(cctkGH) ==0 // Output only on first proc
      && cctk_iteration%output_every == 0) { // Only at n*output_every iteration
    if( num_trackers < output_threshold ) {

      for( int i=0; i < num_trackers; i++ ) { 

        // C++ IO is so easy.

        // Open and truncate the file.  Then output the header.
        ofstream outfile;
        outfile.open(output_filename(i).c_str());
        outfile << "# ShiftTracker" << i << ".asc:" << endl
          << "# itt" << '\t'
          << "time" << '\t'
          << "x" << '\t' << "y" << '\t' << "z" << '\t'
          << "vx" << '\t' << "vy" << '\t' << "vz" << '\t'
          << "ax" << '\t' << "ay" << '\t' << "az" << endl
          << "#==========================================="
          << "============================" << endl;

        // Always make sure to close the file when we are done.
        outfile.close();
      }
    } else {
      // This puts all of the stuff in one big file.
      // I have not tested this output to make sure there are no
      // stupid bugs.
      ofstream outfile;
      outfile.open(output_filename(-1).c_str());
      outfile << "# ShiftTracker.asc:" << endl
        << "# itt" << '\t'
        << "time" << '\t';
      for( int i=0; i < num_trackers; i++ ) {
        outfile << "x" << i << '\t'
          << "y" << i << '\t'
          << "z" << i << '\t'
          << "vx" << i << '\t'
          << "vy" << i << '\t'
          << "vz" << i << '\t'
          << "ax" << i << '\t'
          << "ay" << i << '\t'
          << "az" << i << '\t';
      }
      outfile << endl;
      outfile.close();
    }
  }
}
extern "C" void ShiftTracker_Output(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_MyProc(cctkGH) ==0
      && cctk_iteration%output_every == 0 ) {
    if( num_trackers < output_threshold ) {

      for( int i=0; i < num_trackers; i++ ) {

        // Setup dynamic filenames
        stringstream out_name;
        out_name << "ShiftTracker" << i << ".asc" << ends;

        // Open each file and append data to it.
        ofstream outfile;
        outfile.open(output_filename(i).c_str(), 
            ofstream::out | ofstream::app);
        outfile << setprecision(io_precision);
        outfile << cctk_iteration << '\t'
          << cctk_time << '\t'
          << st_x[i] << '\t'
          << st_y[i] << '\t'
          << st_z[i] << '\t'
          << st_vx[i] << '\t'
          << st_vy[i] << '\t'
          << st_vz[i] << '\t'
          << st_ax[i] << '\t'
          << st_ay[i] << '\t'
          << st_az[i] << endl;
        // Always close files.
        outfile.close();
      }
    } else {
      ofstream outfile;
      outfile.open(output_filename(-1).c_str(),
          ofstream::out | ofstream::app);

      outfile << setprecision(io_precision);
      outfile << cctk_iteration << '\t' << cctk_time << '\t';
      for( int i=0; i < num_trackers; i++ ) {
        outfile << st_x[i] << '\t'
          << st_y[i] << '\t'
          << st_z[i] << '\t'
          << st_vx[i] << '\t'
          << st_vy[i] << '\t'
          << st_vz[i] << '\t'
          << st_ax[i] << '\t'
          << st_ay[i] << '\t'
          << st_az[i] << '\t';
      }
      outfile << endl;
      outfile.close();

    }
  }
}
