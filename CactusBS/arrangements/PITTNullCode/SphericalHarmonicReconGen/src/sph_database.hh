#include <stdio.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"


#include <cstring>
#include <string>
#include <vector>
#include <deque>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <typeinfo>

#include "mpi.h"

#define H5_USE_16_API 1
#include "hdf5.h"

#ifndef _SPH_database_
#define _SPH_database_



namespace SHR {

using namespace std;



/**
   This abstract base-class represents the input-file for easy access. 
*/
class SPH_database
{
   public :
            SPH_database() 
               : _times(vector<CCTK_REAL>(0)), _iterations(vector<CCTK_INT>(0)), _radii(0), 
                 cache(1), cached_first_timestep(0), setup_cache(true),
                 _delta_t(0), _lmax(0), _n_spheres(0), 
                 _n_variables(0), verbose(false) { };
            
            SPH_database(const size_t cached_timesteps, const bool verbose_)
               : _times(vector<CCTK_REAL>(0)),_iterations(vector<CCTK_INT>(0)), _radii(0), 
                 cache(cached_timesteps), cached_first_timestep(0), setup_cache(true),
                 _delta_t(0), _lmax(0), _n_spheres(0),
                 _n_variables(0), verbose(verbose_)  { assert(cached_timesteps >= 1); };
            
	    virtual ~SPH_database() { };
            
            int lmax() const { return _lmax; }
	    int n_variables() const { return _n_variables; }
	    int n_spheres() const { return _n_spheres; }
	    int n_timesteps() const { return _times.size(); }
	    CCTK_REAL time(size_t n) const { assert(n < _times.size()); return _times[n]; } 
	    CCTK_INT  iteration(size_t n) const { assert(n < _iterations.size()); return _iterations[n]; } 
	    CCTK_REAL delta_t() const { return _delta_t; }
            CCTK_REAL radius(const size_t sphere_number) { if (_radii.size() <= sphere_number) return -1; return _radii[sphere_number]; }
            
            // scan file for timesteps, variables, extraction spheres, l-modes
            virtual void scan() = 0;
            
            // read all modes on all extraction spheres for a given timestep and variable number (i.e. column number) into array
            virtual void read(const size_t timestep, const size_t varno, 
                              vector<vector<vector<CCTK_COMPLEX> > >& coeff) const = 0;
            
            // checks if requested timestep and variable is in cache.
            // If yes, it will be put into coeff array
            bool is_in_cache(const size_t timestep, const size_t varno, 
                              vector<vector<vector<CCTK_COMPLEX> > >& coeff) const
            {
               if (setup_cache) return false;
               if (timestep >= cached_first_timestep && timestep-cached_first_timestep < cache.size())
               {
                  for (int sphere=0; sphere < _n_spheres; ++sphere) {
                     for (int l=0; l <= _lmax; ++l) {
                        for (int m=-l; m <= l; ++m) {
                            CCTK_COMPLEX c = CCTK_Cmplx(cache[timestep-cached_first_timestep][varno-1][sphere][l][m+l], 
                                                        cache[timestep-cached_first_timestep][varno][sphere][l][m+l]);
                            coeff[sphere][l][m+l] = c;
                        }
                     }
                  }
                  return true;
               }
               return false;
            }
            
   protected :
            vector<CCTK_REAL> _times;      // all stored times...
            vector<CCTK_INT> _iterations;  // ...and all corresponding iterations
            vector<CCTK_REAL> _radii;       // radius for each sphere number
            
            // cached timesteps and variables for all spheres and lm modes,
            // we always cache at least one timestep so that the the various variables do not need to be read in again.
            mutable vector<vector<vector<vector<vector<CCTK_REAL> > > > > cache;
            // start timestep of cache
            mutable size_t cached_first_timestep;
            // at first read operation, cache must be setup
            mutable bool setup_cache;
            
            // the timestep size found in file
            CCTK_REAL _delta_t;
            
            // number of l-modes
            int _lmax;
            // number of extraction spheres
            int _n_spheres;
            // number of variables
            int _n_variables;
            
            bool verbose;
};


/**
   An ASCII database
*/
class SPH_db_ASCII : public SPH_database
{
   public :    
            SPH_db_ASCII(const string& fname_,
                         const bool verbose,
                         const int cached_timesteps,
			 const int col_time_,
			 const int col_iteration_,
			 const int col_radius_);
            
	    virtual ~SPH_db_ASCII();
            
            // scans the input file
            virtual void scan();
            
            // read all modes on all extraction spheres for a given timestep and variable number into array
            virtual void read(const size_t timestep, const size_t varno, 
                              vector<vector<vector<CCTK_COMPLEX> > >& coeff) const;
            
   private :
            mutable ifstream file;             // ascii input-file stream
            vector<long> byte_offset;  // the byte-offset for each timestep
            
            string _fname;
            
	    int _col_time;
	    int _col_iteration;
	    int _col_radius;
	    
            // a structure holding the data of one line
            struct record_t
            {
               unsigned long long byte_offset;  // the byte position of this line in file
               int                line_length;  // byte length of line
               
               double time;
               int    it;
	       double radius;
	       
	       record_t() : byte_offset(0), line_length(0), time(-100000), it(0), radius(0) { }
            };
            
            
            int num_columns;
            static const int columns_assigned = 4;     // the coefficient file must have at least 4 columns: iteration, time, data (real and imag)
            
            void scan_Carpet_ascii_file();
            
            
            void evaluate_data_records(const record_t& tl, const record_t& pl)
            {
               // tl = this line
               // pl = previous line
            
               // on the following conditions,
               // we register a new dataset
               
               if (tl.it != pl.it || fabs(tl.time - pl.time) > 1e-8)
                  new_dataset(tl);

            }


            void new_dataset(const record_t& tl)
            {
               _times.push_back(tl.time);
               _iterations.push_back(tl.it);
               byte_offset.push_back(tl.byte_offset);
            }
   
};



/**
   An ASCII "dat" database:
   There is one row per time. Each column carries a given (l,m) mode for n-th variable.
   Currently, we support only ONE radius per row.
*/
class SPH_db_DAT : public SPH_database
{
   public :    
            SPH_db_DAT(const string& fname_,
                       const bool verbose,
                       const int cached_timesteps,
                       const int lmax_,
		       const int col_time_,
		       const int col_iteration_,
		       const int col_radius_,
		       const int col_lmax_,
		       const int col_n_variables_,
		       const int col_data_,
		       const bool inverse_m_);
            
	    virtual ~SPH_db_DAT();
            
            // scans the input file
            virtual void scan();
            
            // read all modes on all extraction spheres for a given timestep and variable number into array
            virtual void read(const size_t timestep, const size_t varno, 
                              vector<vector<vector<CCTK_COMPLEX> > >& coeff) const;
            
   private :
            mutable ifstream file;             // ascii input-file stream
            vector<long> byte_offset;  // the byte-offset for each timestep
            
            string _fname;
            
	    int _col_time;
	    int _col_iteration;
	    int _col_radius;
	    int _col_lmax;
	    int _col_n_variables;
	    int _col_data;
	    
            // a structure holding the data of one line
            struct record_t
            {
               unsigned long long byte_offset;  // the byte position of this line in file
               int                line_length;  // byte length of line
               
               double time;
               int    it;
	       double radius;
	       int    lmax;
	       int    n_variables;
	       
	       record_t() : byte_offset(0), line_length(0), time(-100000), it(0), radius(0),
	                    lmax(-1), n_variables(-1) { }
            };
            
            
            int num_columns;
            int columns_assigned;
            
            bool inverse_m;
            
            void scan_dat_file();
            
            
            void evaluate_data_records(const record_t& tl, const record_t& pl)
            {
               // tl = this line
               // pl = previous line
            
               // on the following conditions,
               // we register a new dataset
            
               if (tl.it != pl.it || fabs(tl.time - pl.time) > 1e-8)
                  new_dataset(tl);

            }


            void new_dataset(const record_t& tl)
            {
               _times.push_back(tl.time);
               _iterations.push_back(tl.it);
               byte_offset.push_back(tl.byte_offset);
            }
   
};




/**
   A SpEC-HDF5 database
*/
class SPH_db_SpEC_H5 : public SPH_database
{
   public :    
            SPH_db_SpEC_H5(const string& fname_,
                      const bool verbose,
                      const int cached_timesteps);
            
	    virtual ~SPH_db_SpEC_H5();
            
            // scans the input file
            virtual void scan();
            
            // read all modes on all extraction spheres for a given timestep and variable number into array
            virtual void read(const size_t timestep, const size_t varno, 
                              vector<vector<vector<CCTK_COMPLEX> > >& coeff) const;
            
   private :
            
            int remove_non_monotonic_steps(const int k);
            void scan_HDF5();
            static herr_t H5iter(hid_t loc_id, const char* name, const H5L_info_t* info, void* operator_data);
            
            
            string _fname;
            
            // the hdf5 file handle
            hid_t file;

            static const string time_attrib_name;
            static const string step_basename;
            
            // table of varibale names: this maps a column number to a group.
            // The first entry per varibale is the main group, the second is a possible subgroup, e.g.
            //                /g/Step00001/xx
            //    for SpEC
            static const string varname_table[30][2];
            
            struct timestep_t
            {
               string name;
               double time;
               timestep_t(const string& name_, const double time_)
                  : name(name_), time(time_)
               {}
               bool operator<(const timestep_t& rhs) const { return time < rhs.time; }
            };
            
            /// the cylces and names of the steps contained in the file
            vector<timestep_t> steps;
   
};


/**
   A SpEC-HDF5-v2 database
*/
class SPH_db_SpEC_H5_v2 : public SPH_database
{
   public :    
            SPH_db_SpEC_H5_v2(const string& fname_,
                      const bool verbose,
                      const int cached_timesteps);
            
	    virtual ~SPH_db_SpEC_H5_v2();
            
            // scans the input file
            virtual void scan();
            
            // read all modes on all extraction spheres for a given timestep and variable number into array
            virtual void read(const size_t timestep, const size_t varno, 
                              vector<vector<vector<CCTK_COMPLEX> > >& coeff) const;
            
   private :
            int remove_non_monotonic_steps(const int k);
            void scan_HDF5();
            
            
            // table of varibale names: this maps a column number to a group.
            // The first entry per varibale is the main group, the second is a possible subgroup, e.g.
            //                /g/Step00001/xx
            //    for SpEC
            static const string varname_table[30];
            
            string _fname;
            
            // the hdf5 file handle
            hid_t file;

            
            
            vector<CCTK_REAL> steps;
            
};





/**
   This class represents a set of spherical coefficients for one variable and over all extraction spheres.
   It can read-in the coefficients from a database and reconstruct the angular dependence by
   recomposing with the basis-functions, i.e. the sYlm.
   In addition it can return radial and time derivatives.
*/
template <typename T> // type of recomposed variable (CCTK_REAL or CCTK_COMPLEX).
class spherical_decomposed_variable
{
   public :
      
            spherical_decomposed_variable(cGH* const cctkGH_, 
                                          const int n_timelevels_for_coeff_in_file_, 
                                          const int n_timelevels_,
                                          const int timelevel_,
					  const int time_fd_order_,
					  const double dt_,
                                          const int varno_,
                                          const int sphere_no_,
                                          const SPH_database& db_,
                                          const MPI_Comm comm_world_) 
              : _lmax(db_.lmax()), 
                _n_timelevels(n_timelevels_), 
                _n_timelevels_for_coeff_in_file(n_timelevels_for_coeff_in_file_),
                _varno(varno_), _sphere_no(sphere_no_), _timelevel(timelevel_), _time_fd_order(time_fd_order_),
                _db(db_),
                _dt(dt_),
                cctkGH(cctkGH_), 
                _comm_world(comm_world_)
            {
	       // we will do 2nd-order centered time-derivative
               if (_time_fd_order == 2)
               {
                  assert(_n_timelevels == 3);
                  assert(_timelevel == 1);
               }
               if (_time_fd_order == 4)
               {
                  assert(_n_timelevels == 5);
                  assert(_timelevel == 2);
               }
	       
	       assert(_sphere_no < _db.n_spheres());
	       assert(_varno < _db.n_variables());
	       
               // set up coefficients array
               coeff = deque<vector<vector<vector<CCTK_COMPLEX> > > >(_n_timelevels);
               coeff_in_file = deque<vector<vector<vector<CCTK_COMPLEX> > > >(_n_timelevels_for_coeff_in_file);
               coeff_in_file_times = deque<CCTK_REAL>(_n_timelevels_for_coeff_in_file);
               
               for (int i=0; i < _n_timelevels; ++i)
               {
                  coeff[i] = vector<vector<vector<CCTK_COMPLEX> > >(_db.n_spheres());
                  
                  for (int j=0; j < _db.n_spheres(); ++j)
                  {
                     coeff[i][j] = vector<vector<CCTK_COMPLEX> >(_db.lmax()+1);
                     
                     for (int l=0; l <= _db.lmax(); ++l)
                     {
                        coeff[i][j][l] = vector<CCTK_COMPLEX>(2*l+1);
                     }
                  }
               }
               
               for (int i=0; i < _n_timelevels_for_coeff_in_file; ++i)
               {
                  coeff_in_file[i] = vector<vector<vector<CCTK_COMPLEX> > >(_db.n_spheres());
                  
                  for (int j=0; j < _db.n_spheres(); ++j)
                  {
                     coeff_in_file[i][j] = vector<vector<CCTK_COMPLEX> >(_db.lmax()+1);
                     
                     for (int l=0; l <= _db.lmax(); ++l)
                     {
                        coeff_in_file[i][j][l] = vector<CCTK_COMPLEX>(2*l+1);
                     }
                  }
               }
               
               CCTK_REAL* tmp = (CCTK_REAL *) (CCTK_VarDataPtrB(cctkGH, 0, CCTK_VarIndex("NullGrid::null_delta"), NULL));
               
               null_delta[0] = tmp[0];
               null_delta[1] = tmp[1];
               
               assert(null_delta != NULL);
            }
            
            virtual ~spherical_decomposed_variable() { }
            
            int varno() const { return _varno; }
            
            int lmax() const { return _lmax; }
            int n_timelevels() const { return _n_timelevels; }
            
            
            // read all modes on all extraction spheres for the given variable-number "var_no" from SPHERICALHARMONICS output file
            // timelevel: for which timelevel in allocated array to store the coeffs
            // iteration: from which iteration to read in input file
            void set_modes_from_file(int timelevel, int timestep)
            {
               _db.read(timestep, _varno, coeff_in_file[timelevel]);
               coeff_in_file_times[timelevel] = _db.time(timestep);

               // remove non-finite values
               for (int s=0; s < _db.n_spheres(); ++s)         
                  for (int l=0; l <= _db.lmax(); ++l)
                     for (int m=-l; m <= l; ++m)
                         if (!finite(coeff_in_file[timelevel][s][l][m+l].real()) || !finite(coeff_in_file[timelevel][s][l][m+l].imag()))
                         {
                            coeff_in_file[timelevel][s][l][m+l] = CCTK_Cmplx(0.0, 0.0);
                         }
            }
            
            // set a mode by hand
            void set_mode(const int l, const int m, const CCTK_COMPLEX val)
            {
               assert(l <= _db.lmax());
               assert(_timelevel < _n_timelevels_for_coeff_in_file);
               assert(_sphere_no < _db.n_spheres());
               assert(m >= -l && m <=l);
               
               coeff_in_file[_timelevel][_sphere_no][l][m+l] = val;
            }
            
            // broadcast coefficients to all processors
            void sync(const int timelevel)
            {
	       int na = _lmax*_lmax + 2*_lmax + 1;
	       
	       vector<double> re(na);
	       vector<double> im(na);
	       
	       for (int l=0; l <= _lmax; ++l)
	       {
		  for (int m=-l; m <= l; ++m)
		  {
		     const int n = l*l + l + m;
		     re[n] = coeff[timelevel][_sphere_no][l][m+l].real();
		     im[n] = coeff[timelevel][_sphere_no][l][m+l].imag();
		  }
	       }
		  
               const int root = 0; // proc 0 sends the data
	       MPI_Bcast(&re.front(), na, MPI_DOUBLE, root, _comm_world);
	       MPI_Bcast(&im.front(), na, MPI_DOUBLE, root, _comm_world);
	       
	       // copy back
	       for (int l=0; l <= _lmax; ++l)
	       {
		  for (int m=-l; m <= l; ++m)
		  {
		     const int n = l*l + l + m;
		     coeff[timelevel][_sphere_no][l][m+l] = CCTK_Cmplx(re[n], im[n]);
		  }
	       }
            }
            
            // get mode
            CCTK_COMPLEX get_mode(const int l, const int m) const
            {
               return coeff[_timelevel][_sphere_no][l][m+l];
            }
            
            // get time of mode in file on timelevel "tl"
            CCTK_REAL get_time_of_coeff_in_file(const int tl) const
            {
               assert(tl < _n_timelevels_for_coeff_in_file);
               return coeff_in_file_times[tl];
            }
            
            
            // get mode
            CCTK_COMPLEX get_mode_dt(const int l, const int m) const
            {
               assert(_time_fd_order == 2 || _time_fd_order ==4);
               
	       CCTK_REAL re = -1e6;
	       CCTK_REAL im = -1e6;
	       if (_time_fd_order == 2)
	       {
	          assert(_n_timelevels == 3);
	          assert(_timelevel == 1);
		  re = (coeff[_timelevel-1][_sphere_no][l][m+l] - coeff[_timelevel+1][_sphere_no][l][m+l]).real() /(2*_dt);
		  im = (coeff[_timelevel-1][_sphere_no][l][m+l] - coeff[_timelevel+1][_sphere_no][l][m+l]).imag() /(2*_dt);
	       }
	       if (_time_fd_order == 4)
	       {
	          assert(_n_timelevels == 5);
	          assert(_timelevel == 2);
		  re = 1.0/(12.0*_dt) * (- coeff[_timelevel-2][_sphere_no][l][m+l] 
					  + coeff[_timelevel+2][_sphere_no][l][m+l]
					  + 8.0*(  coeff[_timelevel-1][_sphere_no][l][m+l] 
						 - coeff[_timelevel+1][_sphere_no][l][m+l] )).real();
		  im = 1.0/(12.0*_dt) * (- coeff[_timelevel-2][_sphere_no][l][m+l] 
					  + coeff[_timelevel+2][_sphere_no][l][m+l]
					  + 8.0*(  coeff[_timelevel-1][_sphere_no][l][m+l] 
						 - coeff[_timelevel+1][_sphere_no][l][m+l] )).imag();
	       }

               return CCTK_Cmplx(re, im);;
            }
            
            /// returns the iteration that is just before time "time"
            int get_iteration(const CCTK_REAL time)
            {
               for (int j=0; j < _db.n_timesteps(); ++j)
                  if (_db.time(j) > time && fabs(_db.time(j) - time) > 1e-8)
                     return j-1;
               
               // error, we cycled to last iteration
               return -1;
            }
            
            // this routine will cycle all timelevels one timestep to the past
            // The last timelevel is therefore abandoned and the first timelevel
            // needs to be filled!
            void cycle_timelevels()
            {
               if (_n_timelevels == 1) return;
               
               coeff_in_file.push_front(coeff_in_file.back());
               coeff_in_file.pop_back();
            }
            
            
            void copy_timelevels()
            {
               //cout << target_dt << endl;
               //cout << current_iteration << endl;
               //cout << iteration_in_file << endl;
               
               assert(_n_timelevels == _n_timelevels_for_coeff_in_file);
               
               {
                  for (int i=0; i < _n_timelevels; ++i)
                  {
                     for (int j=0; j < _db.n_spheres(); ++j)
                     {
                        for (int l=0; l <= _db.lmax(); ++l)
                        {
                           for (int m=-l; m <= l; ++m)
                           {
                              coeff[i][j][l][m+l] = coeff_in_file[i][j][l][m+l];
                           }
                        }
                     }
                  }
               }
            }
            
            // use all timelevels stored in coeffs_in_file to interpolate and
            // store outcome in timelevel "tl" of "coeff"
            // The vector L stores precomputed Langrange interpolation coefficients
            // for a previously determined target time
            void interpolate_to_(const int tl, const vector<CCTK_REAL>& L)
            {
               //cout << target_dt << endl;
               //cout << current_iteration << endl;
               //cout << iteration_in_file << endl;
               
               assert(tl < _n_timelevels);
               assert(_n_timelevels_for_coeff_in_file == 6); // 6 timelevels for 5th-order interpolation
               
               // interpolate
               {
                  for (int j=0; j < _db.n_spheres(); ++j)
                  {
                     for (int l=0; l <= _db.lmax(); ++l)
                     {
                        for (int m=-l; m <= l; ++m)
                        {
                           const CCTK_COMPLEX f[6] = { coeff_in_file[5][j][l][m+l],
                                                       coeff_in_file[4][j][l][m+l],
                                                       coeff_in_file[3][j][l][m+l],
                                                       coeff_in_file[2][j][l][m+l],
                                                       coeff_in_file[1][j][l][m+l],
                                                       coeff_in_file[0][j][l][m+l] };
                           
                           coeff[tl][j][l][m+l] = CCTK_Cmplx(0, 0);
                           for (int p=0; p <= 5; ++p)
                              coeff[tl][j][l][m+l] += CCTK_Cmplx(L[p]*f[p].real(), L[p]*f[p].imag());
                        }
                     }
                  }
               }
            }
            
            
            
            
   private :
            // (interpolated) harmonic coefficients for all l and m's and all allocated timesteps and extraction spheres
            deque<  // per timelevel
            vector<  // per sphere
            vector<  // per l
            vector<  // per m
                   CCTK_COMPLEX> > > > coeff;
            
            // times of interpolated coeffs above
            deque<CCTK_REAL> coeff_in_file_times;
            
            // harmonic coefficients for all l and m's and all allocated timesteps and extraction spheres as read from file
            deque<  // per timelevel
            vector<  // per sphere
            vector<  // per l
            vector<  // per m
                   CCTK_COMPLEX> > > > coeff_in_file;
            // times of coeffs_in_file above
            deque<CCTK_REAL> times;
            
            // the number of modes
            int _lmax;
            // the number of allocated timelevels 
            // (does not necessarily correspond to the number of timesteps in the input file!)
            // This parameter relates to the "coeffs" array
            int _n_timelevels;
            // This parameter relates to the "coeff_in_file" array
            int _n_timelevels_for_coeff_in_file;
            
            
            // variable number
            int _varno;
            
            // the sphere number
            int _sphere_no;
            
            // the timelevel we operate on (important if we need to take time derivatives)
            int _timelevel;
            
	    // order of accuracy of time derivative
	    int _time_fd_order;
	    
            // the underlying database
            const SPH_database& _db;
            
	    // the delta time
	    double _dt;
	    
            cGH* const cctkGH;
            
	    MPI_Comm _comm_world;
	    
            // pointers to Cactus variable "null_delta"
            CCTK_REAL null_delta[2];
};



}


#endif
