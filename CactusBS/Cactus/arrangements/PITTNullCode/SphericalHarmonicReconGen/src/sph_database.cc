#include "sph_database.hh"
#include <cassert>

namespace SHR {

using namespace std;

// check return code of HDF5 call and print a warning in case of an error
#define HDF5_ERROR(fn_call) HDF5_CALL(__LINE__, __FILE__, CCTK_THORNSTRING,   \
                                      fn_call, #fn_call)
static long long int HDF5_CALL(const int line, const char* file,
                     const char* thornstring, long long int error_code,
                     const char* fn_call)
{
  // guard against type size changes
  if (sizeof(hid_t) > sizeof(error_code) || sizeof(hsize_t) > sizeof(error_code)) {
    CCTK_VError (__LINE__, __FILE__, CCTK_THORNSTRING,
                 "Incorrect size of hid_t used, sizeof(hid_t) = %d, sizeof(hsize) = %d, sizeof(error_code) = %d",
                 (int)sizeof(hid_t), (int)sizeof(hsize_t), (int)sizeof(error_code));
  }

  if (error_code < 0)
  {
    CCTK_VError (line, file, thornstring,
                "HDF5 call '%s' returned error code %d",
                fn_call, int(error_code));
  }
  return error_code;
}

SPH_db_ASCII::SPH_db_ASCII(const string& fname_,
                         const bool verbose_,
                         const int cached_timesteps,
			 const int col_time_,
			 const int col_iteration_,
			 const int col_radius_) 
	       : SPH_database(cached_timesteps, verbose_), _fname(fname_),
		 _col_time(col_time_), _col_iteration(col_iteration_), _col_radius(col_radius_)
{ 
   file.open(fname_.c_str()); 
   scan(); 
}
            
SPH_db_ASCII::~SPH_db_ASCII() 
{ 
   file.close(); 
}
            

void SPH_db_ASCII::scan()
{
   scan_Carpet_ascii_file();
   
   // check continuity and delta_t of timesteps
   if (_times.size() > 1)
   {
      stringstream str;
      str << "SPH_db_ASCII " << _fname << ": t_0 = " <<  _times[0] << ",    t_final = " << _times.back();
      CCTK_INFO(str.str().c_str());
      _delta_t = 0;
      for (size_t i=1; i < _times.size(); ++i)
      {
         if (fabs(_times[i] - _times[i-1]) < 1e-8 || _times[i] < _times[i-1]) {
            ostringstream str;
            str << "Time-step not strictly monotonic: time[" << i-1 << "] = " << _times[i-1] << ", time[" << i << "] = " << _times[i];
            CCTK_WARN(0, str.str().c_str());
         }
         if (verbose && fabs(_times[i] - _times[i-1] - _delta_t) >= 1e-8)
         {
            stringstream str;
            _delta_t = _times[i]-_times[i-1];
            str << "SPH_db_ASCII " << _fname << ": t = " <<  _times[i-1] << ",    delta_t = " << _delta_t;
            CCTK_INFO(str.str().c_str());
         }
      }
      // set _delta_t to correspond to initial delta_t
      _delta_t = _times[1]-_times[0];
   }
   
}
            


            // read all modes on all extraction spheres for a given timestep and variable number into array
void SPH_db_ASCII::read(const size_t timestep, const size_t varno, 
		  vector<vector<vector<CCTK_COMPLEX> > >& coeff) const
{
   // check cache if data is there and return it
   if (is_in_cache(timestep, varno, coeff)) return;

   // fill cache -------------------------------------
   
   cached_first_timestep = timestep;
   
   // set up cache
   if (setup_cache) {
      for (size_t t=0; t < cache.size(); ++t) {
         cache[t].resize(_n_variables);
         for (int v=0; v < _n_variables; ++v) {
            cache[t][v].resize(_n_spheres);
            for (int s=0; s < _n_spheres; ++s) {
               cache[t][v][s].resize(_lmax+1);
               for (int l=0; l <= _lmax; ++l) {
                  cache[t][v][s][l].resize(2*l+1);
               }
            }
         }
      }
      // cache needs to be setup only at first time
      setup_cache = false;
   }
   
   file.clear();
   
   // move-pointer to byte-offset for requested timestep
   file.seekg(byte_offset[timestep], ios::beg);
   
   char *line_buffer = new char[10000];
   
   istringstream instream;  // input string-stream
   string s;
   
   double val;
   int sphere = 0, l = 0, m = 0;
   int t = 0;
   
   size_t byte_block_end;
   if (timestep+cache.size() >= byte_offset.size())
      byte_block_end = byte_offset[byte_offset.size()-1];
   else
      byte_block_end = byte_offset[timestep+cache.size()];
   
   while (file.good() && (size_t)file.tellg() < byte_block_end)
   {
      if (file.tellg() >= byte_offset[timestep+t+1])
      {
         ++t;
         l=0;
	 m=0;
	 sphere=0;
      }
      // reset input-stringstream (clear from previous errors etc.)
      instream.clear();
      s.clear();
      
      // read line from file (max 10000 chars)
      file.getline(line_buffer, 10000);
      
      // set input -stringstream to line-buffer string
      s = (const char*) line_buffer;
      instream.str(s);
      
      // only read line if it is not a comment
      int this_line_num_columns = 0;
      if (! (instream.str().c_str()[0] == '#'))
      {
	 
	 instream >> ws;  // get rid of white spaces at beginning of stream
	 
	 // try to read values until we hit end of stream (which is end of line)
	 while (!instream.eof())
	 {
	    
	    // if value was successfully extracted, we count this as a column
	    // and extract the data
	    if (instream >> val)
	    {
	       this_line_num_columns++;
	       
	       // based on the current column, we assign the data to the cache
	       cache[t][this_line_num_columns-1][sphere][l][m] = (CCTK_REAL) val;
	    }
	    else
	    {
	       CCTK_WARN(0, "Error in reading spherical coefficient file: not a number!");
	    }
	    instream >> ws; // get rid of whitespaces
	 }
      }
            
      // if this line had no data continue with the next one, then we can assume that we 
      // the next data block is the next extraction sphere.
      if (this_line_num_columns == 0)
      {
	 l=0;
	 m=0;
	 sphere++;
	 continue;
      }
      else
      {
	 m++;
	 
	 if (m >= 2*l+1)
	 {
	    l++;
	    m=0;
	 }
      }
   }
   
   delete [] line_buffer;
   
   
   // finally set coeff array according to cache
   if (!is_in_cache(timestep, varno, coeff)) 
      CCTK_WARN(0, "Variable not in cache even though I have put it there!");
}

            
            
void SPH_db_ASCII::scan_Carpet_ascii_file()
{
   double val;
   
   byte_offset = vector<long>(0);
   
   record_t this_line;  // the extracted data-info of the current line
   record_t prev_line;  // the extracted data-info of the previous line
   
   int datasets_per_line = -1; // the number of the actual data in one line
   int first_line = 1;
   
   if (!file || !file.good()) 
   {
      // error
      CCTK_WARN(0, "Error in opening file: file==NULL or !file.good().");
      return;  
   }
   
   char *line_buffer = new char[10000];
   
   istringstream instream;  // input string-stream
   string s;
   
   int n_blank_lines = 0;
   int n_lines = 0;
   _n_spheres = 1;  // we assume we have at least one sphere in that file
   _radii.resize(_n_spheres);
   
   while (file.good())
   {
      // reset input-stringstream (clear from previous errors etc.)
      instream.clear();
      s.clear();
      
      // store byte offset of this line
      this_line.byte_offset = file.tellg();
      
      // read line from file (max 10000 chars)
      file.getline(line_buffer, 10000);
      
      // set input -stringstream to line-buffer string
      s = (const char*) line_buffer;
      instream.str(s);
      
      // only read line if it is not a comment
      int this_line_num_columns = 0;
      if (! (instream.str().c_str()[0] == '#'))
      {
	 instream >> ws;  // get rid of white spaces at beginning of stream
	 
	 // try to read values until we hit end of stream (which is end of line)
	 while (!instream.eof())
	 {
	    // if value was successfully extracted, we count this as a column
	    // and extract the data
	    if (instream >> val)
	    {
	       this_line_num_columns++;
	       
	       // based on the current column, we assign the data
	       // according to the Carpet column style
	       if (this_line_num_columns == _col_iteration)   // iteration
		  this_line.it = (int) val;
	       if (this_line_num_columns == _col_time)  // time
		  this_line.time = (double) val;
	       if (this_line_num_columns == _col_radius)  // radius
		  this_line.radius = (int) val;
	    }
	    else
	    {
	       CCTK_WARN(0, "Error in reading spherical coefficient file: not a number!");
	    }
	    instream >> ws; // get rid of whitespaces
	 }
      }
      
      // if this line had no data continue with the next one
      if (this_line_num_columns == 0)
      {
	 n_lines = 0;
	 n_blank_lines++;
	 continue;
      }
      
      // Single empty lines separate different spheres
      if (n_blank_lines == 1 && _iterations.size() == 1)
      {
         _radii[_n_spheres-1] = prev_line.radius;
	 _n_spheres++;
	 _radii.resize(_n_spheres);
      }
      
      n_blank_lines = 0;
      n_lines++;
      
      // count number of lines with data until seperated by blank line
      // All contiguous data lines are assumed to count the lm-modes: 00, 1-1, 10, 11, 2-2, 2-1, 20, 21, 22, ...
      // In case there is no blank line between current and next timestep, we assume that all lines belong to one sphere.
      if (_iterations.size() == 1)
      {
         if (sqrt(n_lines)-1 > _lmax) _lmax = sqrt(n_lines)-1;
      }
      
      if (first_line && this_line_num_columns != 0)  // if we are in the first line, set num_columns etc according to non-empty first line!
      {
	 num_columns = this_line_num_columns;
	 
	 if (num_columns < columns_assigned)
	 {
	    CCTK_WARN(0, "The spherical coefficient file contains less columns than I expected!");
	 }
	 // number of datasets corresponds to number of columns
	 datasets_per_line = num_columns;
	 
	 new_dataset(this_line);
      }
      
      if (!first_line && num_columns != this_line_num_columns)
      {
	 CCTK_WARN(0, "Error in reading spherical coefficient file: Number of columns of the current line have to be equal to other lines!");
      }
      
      // compare the data of current line with data of previous line
      // and decide whether we have a new iteration etc...
      if (!first_line)
	 evaluate_data_records(this_line, prev_line);

      prev_line = this_line;
      
      if (first_line)
	 first_line = 0;
   }
   
   // number of columns
   _n_variables = datasets_per_line;
   
   
   delete [] line_buffer;
   
   file.clear();
   
   {
   stringstream str;
   str << "SPH_db_ASCII " << _fname << ": n_spheres   = " << _n_spheres;
   CCTK_INFO(str.str().c_str());
   }
   if (_col_radius >= 0)
      for (int i=0; i < _n_spheres; ++i) {
         stringstream str;
         str << "SPH_db_ASCII " << _fname << ":         radius[" << i << "]   = " << _radii[i];
         CCTK_INFO(str.str().c_str());
      }
   {
   stringstream str;
   str << "SPH_db_ASCII " << _fname << ": n_columns   = " << _n_variables;
   CCTK_INFO(str.str().c_str());
   }
   {
   stringstream str;
   str << "SPH_db_ASCII " << _fname << ": n_variables = " << _n_variables;
   CCTK_INFO(str.str().c_str());
   }
   {
   stringstream str;
   str << "SPH_db_ASCII " << _fname << ": lmax        = " << _lmax;
   CCTK_INFO(str.str().c_str());
   }
   {
   stringstream str;
   str << "SPH_db_ASCII " << _fname << ": timesteps   = " << _iterations.size();
   CCTK_INFO(str.str().c_str());
   }
   
   // we append to byte_offset the end of the file
   byte_offset.push_back(file.tellg());
}
            







//-------------------------------------------------------------------------------------------
//
//
// Dat format goes here...........
//
//
//------------------------------------------------------------------------------------------




SPH_db_DAT::SPH_db_DAT(const string& fname_,
                       const bool verbose_,
                       const int cached_timesteps,
                       const int lmax_,
		       const int col_time_,
		       const int col_iteration_,
		       const int col_radius_,
		       const int col_lmax_,
		       const int col_n_variables_,
		       const int col_data_,
		       const bool inverse_m_) 
	       : SPH_database(cached_timesteps, verbose_), _fname(fname_),
		 _col_time(col_time_), _col_iteration(col_iteration_), _col_radius(col_radius_), _col_lmax(col_lmax_),
		 _col_n_variables(col_n_variables_), _col_data(col_data_), columns_assigned(0), inverse_m(inverse_m_)
{ 
   _lmax = lmax_;
   _n_spheres = 1;

   columns_assigned += _col_time >= 0 ? 1 : 0;
   columns_assigned += _col_radius >= 0 ? 1 : 0;
   columns_assigned += _col_iteration >= 0 ? 1 : 0;
   columns_assigned += _col_lmax >= 0 ? 1 : 0;
   columns_assigned += _col_n_variables >= 0 ? 1 : 0;

   if (columns_assigned > _col_data)
      CCTK_WARN(0, "Number of data column is too small.");

   // at least (0,0) mode of one variable (real and imag part) must be there
   columns_assigned += 2;
   
   stringstream str;
   str << "SPH_db_DAT " << _fname << ": Expecting at least " << columns_assigned << " columns in boundary data file.";
   CCTK_INFO(str.str().c_str());

   file.open(fname_.c_str()); 
   scan(); 
}
            
SPH_db_DAT::~SPH_db_DAT() 
{ 
   file.close(); 
}
            

void SPH_db_DAT::scan()
{
   scan_dat_file();
   
   // check continuity and delta_t of timesteps
   if (_times.size() > 1)
   {
      stringstream str;
      str << "SPH_db_DAT " << _fname << ": t_0 = " <<  _times[0] << ",    t_final = " << _times.back();
      CCTK_INFO(str.str().c_str());
      _delta_t = 0;
      for (size_t i=1; i < _times.size(); ++i)
      {
         if (fabs(_times[i] - _times[i-1]) < 1e-8 || _times[i] < _times[i-1]) {
            ostringstream str;
            str << "Time-step not strictly monotonic: time[" << i-1 << "] = " << _times[i-1] << ", time[" << i << "] = " << _times[i];
            CCTK_WARN(0, str.str().c_str());
         }
         if (verbose && fabs(_times[i] - _times[i-1] - _delta_t) >= 1e-8)
         {
            _delta_t = _times[i]-_times[i-1];
            ostringstream str;
            str << "SPH_db_DAT " << _fname << ": t = " <<  _times[i-1] << ",    delta_t = " << _delta_t;
            CCTK_INFO(str.str().c_str());
         }
      }
      // set _delta_t to correspond to initial delta_t
      _delta_t = _times[1]-_times[0];
   }
}
            
            // read all modes on all extraction spheres for a given timestep and variable number into array
void SPH_db_DAT::read(const size_t timestep, const size_t varno, 
		  vector<vector<vector<CCTK_COMPLEX> > >& coeff) const
{

   // check cache if data is there and return it
   if (is_in_cache(timestep, varno*2+1, coeff)) return;

   // fill cache -------------------------------------
   
   cached_first_timestep = timestep;
   
   // set up cache
   if (setup_cache) {
      for (size_t t=0; t < cache.size(); ++t) {
         cache[t].resize(2*_n_variables);
         for (int v=0; v < 2*_n_variables; ++v) {
            cache[t][v].resize(_n_spheres);
            for (int s=0; s < _n_spheres; ++s) {
               cache[t][v][s].resize(_lmax+1);
               for (int l=0; l <= _lmax; ++l) {
                  cache[t][v][s][l].resize(2*l+1);
               }
            }
         }
      }
      // cache needs to be setup only at first time
      setup_cache = false;
   }

   //cout << "bla" << endl;

   file.clear();
   
   // move-pointer to byte-offset for requested timestep
   file.seekg(byte_offset[timestep], ios::beg);
   
   char *line_buffer = new char[1000000];
   
   istringstream instream;  // input string-stream
   string s;
   
   double val;
   int sphere = 0, l = 0, m = 0;
   int t = 0;
   int v = 0;
   
   size_t byte_block_end;
   if (timestep+cache.size() >= byte_offset.size())
      byte_block_end = byte_offset[byte_offset.size()-1];
   else
      byte_block_end = byte_offset[timestep+cache.size()];
   
   int column = _col_data;
   
   while (file.good() && (size_t)file.tellg() < byte_block_end)
   {
      if (file.tellg() >= byte_offset[timestep+t+1])
      {
         ++t;
         sphere = 0;
         l=0; m=0;
         v=0;
         column = _col_data;
      }
      // reset input-stringstream (clear from previous errors etc.)
      instream.clear();
      s.clear();
      
      // read line from file (max 1000000 chars)
      file.getline(line_buffer, 1000000);
      
      // set input -stringstream to line-buffer string
      s = (const char*) line_buffer;
      instream.str(s);
      
      // only read line if it is not a comment
      int this_line_num_columns = 0;
      if (! (instream.str().c_str()[0] == '#'))
      {
	 
	 instream >> ws;  // get rid of white spaces at beginning of stream
	 
	 // try to read values until we hit end of stream (which is end of line)
	 while (!instream.eof())
	 {
	    
	    // if value was successfully extracted, we count this as a column
	    // and extract the data
	    if (instream >> val)
	    {
	       this_line_num_columns++;
	       
	       // based on the current column, we assign the data
	       // based on the selected varno (column)
	       
	       if (this_line_num_columns == column)  // real part
		  cache[t][v][sphere][l][inverse_m ? 2*l-m : m] = (CCTK_REAL) val;
	       
	       if (this_line_num_columns == column+1)  { // imaginary-part
		  cache[t][v+1][sphere][l][inverse_m ? 2*l-m : m] = (CCTK_REAL) val;
		  
	          column += 2; // next lm-mode
	          m++;
	          if (m >= 2*l+1)
	          {
	             l++;
	             m=0;
	          }
	          if (l > _lmax) {
	             // read next variable
	             v+=2;
	             l=0;
	             m=0;
	             sphere=0;
	             if (v >= _n_variables*2) {
	                // go to next time step since all variables have been read!
                        file.seekg(byte_offset[timestep+t+1], ios::beg);
	                break; // stop reading further!
	             }
	          }
	       }
	    }
	    else
	    {
	       CCTK_WARN(0, "Error in reading spherical coefficient file: not a number!");
	    }
	    instream >> ws; // get rid of whitespaces
	 }
      }
      
   }
   
   delete [] line_buffer;
      
   // finally set coeff array according to cache
   if (!is_in_cache(timestep, varno*2+1, coeff)) 
      CCTK_WARN(0, "Variable not in cache even though I have put it there!");
}
            
            
void SPH_db_DAT::scan_dat_file()
{
   double val;
   
   byte_offset = vector<long>(0);
   
   record_t this_line;  // the extracted data-info of the current line
   record_t prev_line;  // the extracted data-info of the previous line
   
   int datasets_per_line; // the number of the actual data in one line
   int first_line = 1;
   
   if (!file || !file.good()) 
   {
      // error
      CCTK_WARN(0, "Error in opening file: file==NULL or !file.good().");
      return;  
   }
   
   char *line_buffer = new char[1000000];
   
   istringstream instream;  // input string-stream
   string s;
   
   _n_spheres = 1;  // we assume we have at least one sphere in that file
   while (file.good())
   {
      // reset input-stringstream (clear from previous errors etc.)
      instream.clear();
      s.clear();
      
      // store byte offset of this line
      this_line.byte_offset = file.tellg();
      
      // read line from file (max 1000000 chars)
      file.getline(line_buffer, 1000000);
      
      // set input -stringstream to line-buffer string
      s = (const char*) line_buffer;
      instream.str(s);
      
      // only read line if it is not a comment
      int this_line_num_columns = 0;
      if (! (instream.str().c_str()[0] == '#'))
      {
	 instream >> ws;  // get rid of white spaces at beginning of stream
	 
	 // try to read values until we hit end of stream (which is end of line)
	 while (!instream.eof())
	 {
	    // if value was successfully extracted, we count this as a column
	    // and extract the data
	    if (instream >> val)
	    {
	       this_line_num_columns++;
	       
	       // based on the current column, we assign the data
	       // according to the Carpet column style
	       if (this_line_num_columns == _col_iteration)   // iteration
		  this_line.it = (int) val;
	       if (this_line_num_columns == _col_time)  // time
		  this_line.time = (double) val;
	       if (this_line_num_columns == _col_radius)  // radius
		  this_line.radius = (int) val;
	       if (this_line_num_columns == _col_lmax)  // lmax
		  this_line.lmax = (int) val;
	       if (this_line_num_columns == _col_n_variables)  // n_variables
		  this_line.n_variables = (int) val;
	    }
	    else
	    {
	       CCTK_WARN(0, "Error in reading spherical coefficient file: not a number!");
	    }
	    instream >> ws; // get rid of whitespaces
	 }
      }
      
      // if this line had no data continue with the next one
      if (this_line_num_columns == 0)
	 continue;
      
      if (first_line && this_line_num_columns != 0)  // if we are in the first line, set num_columns etc according to non-empty first line!
      {
	 num_columns = this_line_num_columns;
	 
	 if (_col_lmax >= 0 && this_line.lmax != _lmax) {
	    cout << _col_lmax << ", " << this_line.lmax << ", " << _lmax << endl;
	    CCTK_WARN(0, "I found a different lmax in file than what has been specified by user! Stopping.");
	 }
	 if (num_columns < columns_assigned)
	 {
	    CCTK_WARN(0, "The spherical coefficient file contains less columns than I expected!");
	 }
	 // number of datasets corresponds to number of columns
	 datasets_per_line = num_columns;
	 
	 new_dataset(this_line);
      }
      
      
      // compare the data of current line with data of previous line
      // and decide whether we have a new iteration etc...
      if (!first_line)
	 evaluate_data_records(this_line, prev_line);

      prev_line = this_line;
      
      if (first_line)
	 first_line = 0;
   }
   
   // check that there are as many variables as the user has specified!
   _n_variables = (num_columns-_col_data-1) / ((_lmax+1)*(_lmax+1)*2) + 1;
   
   
   delete [] line_buffer;
   
   file.clear();
   
   {
   stringstream str;
   str << "SPH_db_DAT " << _fname << ": n_spheres   = " << _n_spheres;
   CCTK_INFO(str.str().c_str());
   }
   {
   stringstream str;
   str << "SPH_db_DAT " << _fname << ": n_columns   = " << num_columns;
   CCTK_INFO(str.str().c_str());
   }
   {
   stringstream str;
   str << "SPH_db_DAT " << _fname << ": n_variables = " << _n_variables;
   CCTK_INFO(str.str().c_str());
   }
   {
   stringstream str;
   str << "SPH_db_DAT " << _fname << ": lmax        = " << _lmax;
   CCTK_INFO(str.str().c_str());
   }
   {
   stringstream str;
   str << "SPH_db_DAT " << _fname << ": timesteps   = " << _iterations.size();
   CCTK_INFO(str.str().c_str());
   }
   
   // we append to byte_offset the end of the file
   byte_offset.push_back(file.tellg());
}

            













//-------------------------------------------------------------------------------------------
//
//
// SpEC-H5 format goes here...........
//
//
//------------------------------------------------------------------------------------------




const string SPH_db_SpEC_H5::time_attrib_name = "Time";
const string SPH_db_SpEC_H5::step_basename = "Step";

            // table of varibale names: this maps a column number to a group.
            // The first entry per varibale is the main group, the second is a possible subgroup, e.g.
            //                /g/Step00001/xx
            //    for SpEC
const string SPH_db_SpEC_H5::varname_table[30][2] =
                                                 { { "Lapse", "scalar"},
                                                   { "Shift", "x" }, { "Shift", "y" }, { "Shift", "z" },
                                                   { "g", "xx" }, { "g", "xy" }, { "g", "xz" },
                                                   { "g", "yy" }, { "g", "yz" },
                                                   { "g", "zz" },
                                                   
                                                   { "DrLapse", "scalar"},
                                                   { "DrShift", "x" }, { "DrShift", "y" }, { "DrShift", "z" },
                                                   { "Drg", "xx" }, { "Drg", "xy" }, { "Drg", "xz" },
                                                   { "Drg", "yy" }, { "Drg", "yz" },
                                                   { "Drg", "zz" },
                                                   
                                                   { "DtLapse", "scalar"},
                                                   { "DtShift", "x" }, { "DtShift", "y" }, { "DtShift", "z" },
                                                   { "Dtg", "xx" }, { "Dtg", "xy" }, { "Dtg", "xz" },
                                                   { "Dtg", "yy" }, { "Dtg", "yz" },
                                                   { "Dtg", "zz" } };


SPH_db_SpEC_H5::SPH_db_SpEC_H5(const string& fname_,
                               const bool verbose_,
                               const int cached_timesteps) 
	       : SPH_database(cached_timesteps, verbose_), _fname(fname_)
{ 
   // open file
   file = H5Fopen(_fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   if (file < 0) {
     CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING, 
                 "Could not open file \"%s\".", _fname.c_str());
   }
   
   scan(); 
}
            

SPH_db_SpEC_H5::~SPH_db_SpEC_H5() 
{ 
   // close file
   HDF5_ERROR(H5Fclose(file));
}
            


int SPH_db_SpEC_H5::remove_non_monotonic_steps(const int k)
{
   // go backwards from k and remove until we get monotonicity
   const double time = steps[k].time;
   const int size = steps.size();
   int c = 1;
   for (int i=k-1; i >= 0; --i) {
      if (time > steps[i].time) {
         break;
      }
      ++c;
   }
   
   steps.erase(steps.begin()+k-c, steps.begin()+k);
   _times.erase(_times.begin()+k-c, _times.begin()+k);
   _iterations.erase(_iterations.begin()+k-c, _iterations.begin()+k);
   
   // return the number of steps that were removed
   return size-steps.size();
}


void SPH_db_SpEC_H5::scan()
{

   // scan datasets
   scan_HDF5();
   
   
   // check continuity and delta_t of timesteps
   if (_times.size() > 1)
   {
      stringstream str;
      str << "SPH_db_SpEC_H5 " << _fname << ": t_0 = " <<  _times[0] << ",    t_final = " << _times.back();
      CCTK_INFO(str.str().c_str());
      _delta_t = 0;
      for (size_t i=1; i < _times.size(); ++i)
      {
         if (fabs(_times[i] - _times[i-1]) < 1e-8 || _times[i] < _times[i-1]) {
            ostringstream str;
            str << "Removing non-monotonic timesteps: time[" << i-1 << "] = " << _times[i-1] << ", time[" << i << "] = " << _times[i];
            CCTK_WARN(1, str.str().c_str());
            i -= remove_non_monotonic_steps(i);
            continue;
         }
         if (verbose && fabs(_times[i] - _times[i-1] - _delta_t) >= 1e-8)
         {
            _delta_t = _times[i]-_times[i-1];
            stringstream str;
            str << "SPH_db_SpEC_H5 " << _fname << ": t = " <<  _times[i-1] << ",    delta_t = " << _delta_t;
            CCTK_INFO(str.str().c_str());
         }
      }
      // set _delta_t to correspond to initial delta_t
      _delta_t = _times[1]-_times[0];
   }
   
}




void SPH_db_SpEC_H5::read(const size_t timestep, const size_t varno, 
		  vector<vector<vector<CCTK_COMPLEX> > >& coeff) const
{
   if (timestep >= steps.size()) {
      CCTK_WARN(0, "Requested timestep not in worldtube data file. Stopping!");
   }

   hid_t dataset = HDF5_ERROR(H5Dopen(file, (varname_table[varno][0]+"/"+
                                             steps[timestep].name+"/"+
                                             varname_table[varno][1]).c_str()));
   vector<double> buffer((_lmax+1)*(_lmax+1)*2);
   HDF5_ERROR(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &buffer.front()));
   HDF5_ERROR(H5Dclose(dataset));
   
   int c = 0;
   for (int l=0; l <= _lmax; ++l) {
      for (int m=l; m >= -l; --m, ++c) {
         coeff[0][l][m+l] = CCTK_Cmplx(buffer[2*c], buffer[2*c + 1]);
      }
   }
}



herr_t SPH_db_SpEC_H5::H5iter(hid_t loc_id, const char* name, const H5L_info_t* info, void* operator_data)
{
   SPH_db_SpEC_H5* SpEC_H5 = (SPH_db_SpEC_H5*) operator_data;
   
   // Datasets are ignored! We only iterate over groups
   switch (info->type) {
      case H5O_TYPE_GROUP: {
         const string gname = name;
         // we have found a timstep group if the group contains the "step_basename" substring
         if (gname.find(step_basename) != string::npos) {
            // get time attribute
            double time;
            hid_t group = HDF5_ERROR(H5Gopen(loc_id, name));
            hid_t attribute = HDF5_ERROR(H5Aopen(group, time_attrib_name.c_str(), H5P_DEFAULT));
            HDF5_ERROR(H5Aread(attribute, H5T_NATIVE_DOUBLE, &time));
            HDF5_ERROR(H5Aclose(attribute));
            HDF5_ERROR(H5Gclose(group));
            
            SpEC_H5->steps.push_back(timestep_t(name, time));
         }
         break;
      }
      case H5O_TYPE_DATASET:
         // ignore datasets
         return 0;
         break;
      default:
         break;
   }

   return 0;
}



void SPH_db_SpEC_H5::scan_HDF5()
{
   // pick one variable and browse through entire file
   // to find out how many (and what) cycles and times there are
   // (assuming that all variables are present at the same timesteps)
   HDF5_ERROR(H5Literate_by_name (file, varname_table[0][0].c_str(),
                                  H5_INDEX_NAME, H5_ITER_NATIVE, NULL,
                                  H5iter, this, H5P_DEFAULT));

   _iterations.resize(steps.size());
   _times.resize(steps.size());
   for (size_t i=0; i < steps.size(); ++i) {
      _iterations[i] = i;
      _times[i] = steps[i].time;
   }

   // get lmax from first variable and timestep
   hid_t dataset = HDF5_ERROR(H5Dopen(file, (varname_table[0][0]+"/"+
                                             steps[0].name+"/"+
                                             varname_table[0][1]).c_str()));
   hid_t dataspace = HDF5_ERROR(H5Dget_space(dataset));
   hsize_t size = HDF5_ERROR(H5Sget_simple_extent_npoints(dataspace));
   HDF5_ERROR(H5Sclose(dataspace));
   HDF5_ERROR(H5Dclose(dataset));

   assert (size % 2 == 0);

   _lmax = sqrt(size/2)-1;
   _n_variables = 30;
   _n_spheres = 1;
   
   // cross-check that all expected variables are present
   
   
   // display some info
   
   {
   stringstream str;
   str << "SPH_db_SpEC_H5 " << _fname << ": n_spheres   = " << _n_spheres;
   CCTK_INFO(str.str().c_str());
   }
   {
   stringstream str;
   str << "SPH_db_SpEC_H5 " << _fname << ": n_variables = " << _n_variables;
   CCTK_INFO(str.str().c_str());
   }
   {
   stringstream str;
   str << "SPH_db_SpEC_H5 " << _fname << ": lmax        = " << _lmax;
   CCTK_INFO(str.str().c_str());
   }
   {
   stringstream str;
   str << "SPH_db_SpEC_H5 " << _fname << ": timesteps   = " << _iterations.size();
   CCTK_INFO(str.str().c_str());
   }
}




//-------------------------------------------------------------------------------------------
//
//
// SpEC-H5-v2 format goes here...........
//
//
//------------------------------------------------------------------------------------------




// table of varibale names: this maps a column number to a dataset name.
const string SPH_db_SpEC_H5_v2::varname_table[30] =
                                                 { "Lapse.dat",
                                                   "Shiftx.dat", "Shifty.dat", "Shiftz.dat",
                                                   "gxx.dat", "gxy.dat", "gxz.dat",
                                                   "gyy.dat", "gyz.dat",
                                                   "gzz.dat",
                                                   
                                                   "DrLapse.dat",
                                                   "DrShiftx.dat", "DrShifty.dat", "DrShiftz.dat",
                                                   "Drgxx.dat", "Drgxy.dat", "Drgxz.dat",
                                                   "Drgyy.dat", "Drgyz.dat",
                                                   "Drgzz.dat",
                                                   
                                                   "DtLapse.dat",
                                                   "DtShiftx.dat", "DtShifty.dat", "DtShiftz.dat",
                                                   "Dtgxx.dat", "Dtgxy.dat", "Dtgxz.dat",
                                                   "Dtgyy.dat", "Dtgyz.dat",
                                                   "Dtgzz.dat" };



SPH_db_SpEC_H5_v2::SPH_db_SpEC_H5_v2(const string& fname_,
                               const bool verbose_,
                               const int cached_timesteps) 
	       : SPH_database(cached_timesteps, verbose), _fname(fname_)
{ 
   // open file
   file = H5Fopen(_fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   if (file < 0) {
     CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING, 
                 "Could not open file \"%s\".", _fname.c_str());
   }
   
   scan(); 
}
            

SPH_db_SpEC_H5_v2::~SPH_db_SpEC_H5_v2() 
{ 
   // close file
   HDF5_ERROR(H5Fclose(file));
}


int SPH_db_SpEC_H5_v2::remove_non_monotonic_steps(const int k)
{
   // go backwards from k and remove until we get monotonicity
   const double time = steps[k];
   const int size = steps.size();
   int c = 1;
   for (int i=k-1; i >= 0; --i) {
      if (time > steps[i]) {
         break;
      }
      ++c;
   }
   
   steps.erase(steps.begin()+k-c, steps.begin()+k);
   _times.erase(_times.begin()+k-c, _times.begin()+k);
   _iterations.erase(_iterations.begin()+k-c, _iterations.begin()+k);
   
   // return the number of steps that were removed
   return size-steps.size();
}



void SPH_db_SpEC_H5_v2::scan()
{

   // scan datasets
   scan_HDF5();
   
   
   // check continuity and delta_t of timesteps
   if (_times.size() > 1)
   {
      stringstream str;
      str << "SPH_db_SpEC_H5 " << _fname << ": t_0 = " <<  _times[0] << ",    t_final = " << _times.back();
      CCTK_INFO(str.str().c_str());
      _delta_t = 0;
      for (size_t i=1; i < _times.size(); ++i)
      {
         if (fabs(_times[i] - _times[i-1]) < 1e-8 || _times[i] < _times[i-1]) {
            ostringstream str;
            str << "Removing non-monotonic timesteps: time[" << i-1 << "] = " << _times[i-1] << ", time[" << i << "] = " << _times[i];
            CCTK_WARN(1, str.str().c_str());
            i -= remove_non_monotonic_steps(i);
            continue;
         }
         if (verbose && fabs(_times[i] - _times[i-1] - _delta_t) >= 1e-8)
         {
            _delta_t = _times[i]-_times[i-1];
            stringstream str;
            str << "SPH_db_SpEC_H5_v2 " << _fname << ": t = " <<  _times[i-1] << ",    delta_t = " << _delta_t;
            CCTK_INFO(str.str().c_str());
         }
      }
      // set _delta_t to correspond to initial delta_t
      _delta_t = _times[1]-_times[0];
   }
   
}


void SPH_db_SpEC_H5_v2::read(const size_t timestep, const size_t varno, 
		  vector<vector<vector<CCTK_COMPLEX> > >& coeff) const
{
   if (timestep >= steps.size()) {
      CCTK_WARN(0, "Requested timestep not in worldtube data file. Stopping!");
   }

   // check cache if data is there and return it
   if (is_in_cache(timestep, 2*varno+1, coeff)) return;

   // fill cache -------------------------------------
   
   cached_first_timestep = timestep;
   
   // set up cache
   if (setup_cache) {
      for (size_t t=0; t < cache.size(); ++t) {
         cache[t].resize(2*_n_variables);
         for (int v=0; v < 2*_n_variables; ++v) {
            cache[t][v].resize(_n_spheres);
            for (int s=0; s < _n_spheres; ++s) {
               cache[t][v][s].resize(_lmax+1);
               for (int l=0; l <= _lmax; ++l) {
                  cache[t][v][s][l].resize(2*l+1);
               }
            }
         }
      }
      // cache needs to be setup only at first time
      setup_cache = false;
   }

   
   const int nsteps = timestep+cache.size() >= steps.size() ? steps.size()-timestep : cache.size();

   for (int v=0; v < _n_variables; ++v) {

      hid_t dataset = HDF5_ERROR(H5Dopen(file, varname_table[v].c_str()));
      const int nmodes = (_lmax+1)*(_lmax+1)*2;
      vector<double> buffer(nmodes*nsteps);
   
      hsize_t dimsm[1] = { buffer.size() };
      hid_t memspace = HDF5_ERROR(H5Screate_simple (1, dimsm, NULL));
   
      hid_t dataspace = HDF5_ERROR(H5Dget_space(dataset));
      hsize_t dims[2];
      HDF5_ERROR(H5Sget_simple_extent_dims (dataspace, dims, NULL));
   
      hsize_t offset[2] = { timestep, 1 };
      hsize_t count[2]  = { (hsize_t)nsteps, (hsize_t)nmodes };
      HDF5_ERROR(H5Sselect_hyperslab (dataspace, H5S_SELECT_SET, offset, NULL, 
                                      count, NULL));
   
      HDF5_ERROR(H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT,
                 &buffer.front()));
   
      HDF5_ERROR(H5Sclose(memspace));
      HDF5_ERROR(H5Sclose(dataspace));
      HDF5_ERROR(H5Dclose(dataset));
   
      for (int t=0; t < nsteps; ++t) {
   
         int c = 0;
         for (int l=0; l <= _lmax; ++l) {
            for (int m=l; m >= -l; --m, ++c) {
               cache[t][2*v  ][0][l][m+l] = buffer[t*nmodes + 2*c]; 
               cache[t][2*v+1][0][l][m+l] = buffer[t*nmodes + 2*c + 1];
            }
         }
   
      }
   
   }
   
   // finally set coeff array according to cache
   if (!is_in_cache(timestep, 2*varno+1, coeff)) 
      CCTK_WARN(0, "Variable not in cache even though I have put it there!");
}



void SPH_db_SpEC_H5_v2::scan_HDF5()
{  
   // get lmax from first variable and timestep
   hid_t dataset = HDF5_ERROR(H5Dopen(file, "gxx.dat"));
   hid_t dataspace = HDF5_ERROR(H5Dget_space(dataset));
   hsize_t dims[2];
   HDF5_ERROR(H5Sget_simple_extent_dims (dataspace, dims, NULL));
   
   // first dim is number of timesteps
   steps.resize(dims[0]);
   // second dim is the number of lmodes
   _lmax = sqrt(dims[1]/2)-1;
   
   // select hyperslab in dataset
   hsize_t offset[2] = { 0, 0 };
   hsize_t count[2]  = { dims[0], 1 };
   HDF5_ERROR(H5Sselect_hyperslab (dataspace, H5S_SELECT_SET, offset, NULL, 
                                   count, NULL));
   
   // define dataspace of memory
   hsize_t dimsm[1] = { steps.size() };
   hid_t memspace = H5Screate_simple (1, dimsm, NULL);
   
   // read "times"
   HDF5_ERROR(H5Dread (dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                       H5P_DEFAULT, &steps.front()));
   
   HDF5_ERROR(H5Sclose(memspace));
   HDF5_ERROR(H5Sclose(dataspace));
   HDF5_ERROR(H5Dclose(dataset));

   _n_variables = 30;
   _n_spheres = 1;
   
   _iterations.resize(steps.size());
   _times.resize(steps.size());
   for (size_t i=0; i < steps.size(); ++i) {
      _iterations[i] = i;
      _times[i] = steps[i];
   }
   
   // cross-check that all expected variables are present
   
   
   // display some info
   {
   stringstream str;
   str << "SPH_db_SpEC_H5_v2 " << _fname << ": n_spheres   = " << _n_spheres;
   CCTK_INFO(str.str().c_str());
   }
   {
   stringstream str;
   str << "SPH_db_SpEC_H5_v2 " << _fname << ": n_variables = " << _n_variables;
   CCTK_INFO(str.str().c_str());
   }
   {
   stringstream str;
   str << "SPH_db_SpEC_H5_v2 " << _fname << ": lmax        = " << _lmax;
   CCTK_INFO(str.str().c_str());
   }
   {
   stringstream str;
   str << "SPH_db_SpEC_H5_v2 " << _fname << ": timesteps   = " << _iterations.size();
   CCTK_INFO(str.str().c_str());
   }

}





}
