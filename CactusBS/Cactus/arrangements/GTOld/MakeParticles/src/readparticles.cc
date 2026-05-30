#include <cmath>
#include <cstdio>
#include <cerrno>
#include <cassert>
#include <cstring>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

using namespace std;

#include "carpet.hh"
extern "C" {
// keep after carpet.hh, otherwise I get strange errors
#define KRANC_C
#include "GenericFD.h"
}

#define SQR(x) ((x)*(x))
#define DIM(x) (sizeof(x)/sizeof((x)[0]))

#ifdef ET_HYDROBASE
#define velx (&vel[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely (&vel[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz (&vel[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#endif

/* information type encodings */
#include "makeparticles_types.h"

static vector<CCTK_REAL> particle_data;
extern "C"
void MakeParticles_ReadParticles(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /* parse parameters into numbers */
  enum info_type_t info_type = info_invalid;   // parameter encoding
  for (int t = 0 ; t < DIM(info) ; t++) {
    if ( CCTK_EQUALS(information, info[t].name) ) {
      info_type = info[t].type;
      break;
    }
  }
  assert(info_type != info_invalid);
  assert(info_type == info_fallback);

  // read in data either from all files that match the file_X pattern or from just the file given
  const bool is_multifile = strstr(read_particles_from_file, ".file_0.") != NULL;
  int filenum = 0;
  do {
    const char* fn;
    ostringstream ofn;
    if(is_multifile) {
      int base_len = strstr(read_particles_from_file, ".file_0.") - read_particles_from_file;
      ofn << string(read_particles_from_file).substr(0,base_len) << ".file_" << filenum << "." << string(read_particles_from_file).substr(base_len+strlen(".file_0."));
      fn = ofn.str().c_str();
    } else {
      fn = read_particles_from_file;
    }

    FILE* fh = fopen(fn, "r");
    if(fh != NULL) {
      size_t file_sz, particles_sz;
      fseek(fh, 0, SEEK_END);
      file_sz = ftell(fh);
      fseek(fh, 0, SEEK_SET);
      if(file_sz % (sizeof(particle_data[0])*info[info_type].nvals) != 0) {
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING, "Invalid file '%s': size not multiple of dataset size (%d)", fn, (int)(info[info_type].nvals*sizeof(CCTK_REAL)));
        return; /* NOTREACHED */
      }
      file_sz /= sizeof(particle_data[0]); // makes expressions below a bit shorter, is all
      particles_sz = particle_data.size();
      particle_data.resize(particles_sz+file_sz);
      if(fread(&particle_data[particles_sz], sizeof(particle_data[0]), file_sz, fh) != file_sz) {
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING, "Could not read %d elements from file '%s': %s", (int)file_sz, fn, strerror(errno));
        return; /* NOTREACHED */
      }
      fclose(fh);
    } else {
      if(filenum == 0) {
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING, "Could not open file '%s' for reading: %s", fn, strerror(errno));
        return; /* NOTREACHED */
      }
      break; // quit do..while loop
    }

    filenum += 1;
  } while(is_multifile); // keep trying until we 'break'
}

extern "C" void MakeParticles_FreeParticleData(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  particle_data.clear();
}

extern "C"
void MakeParticles_SetHydroVariables(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL xmin = x[0] - CCTK_DELTA_SPACE(0)/2.;
  const CCTK_REAL ymin = y[0] - CCTK_DELTA_SPACE(1)/2.;
  const CCTK_REAL zmin = z[0] - CCTK_DELTA_SPACE(2)/2.;

  const CCTK_REAL xmax = x[CCTK_GFINDEX3D(cctkGH, cctk_lsh[0]-1, cctk_lsh[1]-1, cctk_lsh[2]-1)] + CCTK_DELTA_SPACE(0)/2.; 
  const CCTK_REAL ymax = y[CCTK_GFINDEX3D(cctkGH, cctk_lsh[0]-1, cctk_lsh[1]-1, cctk_lsh[2]-1)] + CCTK_DELTA_SPACE(1)/2.; 
  const CCTK_REAL zmax = z[CCTK_GFINDEX3D(cctkGH, cctk_lsh[0]-1, cctk_lsh[1]-1, cctk_lsh[2]-1)] + CCTK_DELTA_SPACE(2)/2.; 

  const CCTK_REAL dxfine = cctk_delta_space[0] / Carpet::spacereffacts[Carpet::maxreflevels-1][0];
  const CCTK_REAL dyfine = cctk_delta_space[0] / Carpet::spacereffacts[Carpet::maxreflevels-1][1];
  const CCTK_REAL dzfine = cctk_delta_space[0] / Carpet::spacereffacts[Carpet::maxreflevels-1][2];

  /* parse parameters into numbers */
  enum info_type_t info_type;               // parameter encoding
  info_type = info_invalid;
  for (int t = 0 ; t < DIM(info) ; t++) {
    if ( CCTK_EQUALS(information, info[t].name) ) {
      info_type = info[t].type;
      break;
    }
  }
  assert(info_type != info_invalid);
  
  assert(info_type == info_fallback);
  enum {i_x ,i_y ,i_z ,i_vx ,i_vy ,i_vz ,i_gxx ,i_gxy ,i_gxz ,i_gyy ,i_gyz ,i_gzz ,i_betax ,i_betay ,i_betaz ,i_alpha ,i_mass ,i_energy_per_unit_mass ,i_rho ,i_eps};

  // first clear to vacuum
  for(size_t idx = 0 ; idx < cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2] ; idx++) {
    rho[idx] = eps[idx] = press[idx] = velx[idx] = vely[idx] = velz[idx] = 0.0;
  }

  // put in values from particles
  int dp=info[info_type].nvals;
#pragma omp parallel for schedule(dynamic)
  for(int p = 0 ; p < int(particle_data.size()) ; p+=dp) {
    CCTK_REAL *part = &particle_data[p];
    CCTK_REAL dx, dy, dz;
    if(part[i_x] > xmin && part[i_x] < xmax &&
       part[i_y] > ymin && part[i_y] < ymax &&
       part[i_z] > zmin && part[i_z] < zmax) {
      // check that the point lies on one of our grid points
      dx = fmod(part[i_x] - CCTK_ORIGIN_SPACE(0), CCTK_DELTA_SPACE(0));
      if(dx > CCTK_DELTA_SPACE(0)/2)
        dx -= CCTK_DELTA_SPACE(0);
      dy = fmod(part[i_y] - CCTK_ORIGIN_SPACE(1), CCTK_DELTA_SPACE(1));
      if(dy > CCTK_DELTA_SPACE(1)/2)
        dy -= CCTK_DELTA_SPACE(1);
      dz = fmod(part[i_z] - CCTK_ORIGIN_SPACE(2), CCTK_DELTA_SPACE(2));
      if(dz > CCTK_DELTA_SPACE(2)/2)
        dz -= CCTK_DELTA_SPACE(2);
      if(fabs(dx) < dxfine/2 && fabs(dy) < dyfine/2 && fabs(dz) < dzfine/2) {
        const CCTK_INT i = CCTK_INT(round((part[i_x] - CCTK_ORIGIN_SPACE(0)) / CCTK_DELTA_SPACE(0))) - cctk_lbnd[0];
        const CCTK_INT j = CCTK_INT(round((part[i_y] - CCTK_ORIGIN_SPACE(1)) / CCTK_DELTA_SPACE(1))) - cctk_lbnd[1];
        const CCTK_INT k = CCTK_INT(round((part[i_z] - CCTK_ORIGIN_SPACE(2)) / CCTK_DELTA_SPACE(2))) - cctk_lbnd[2];
        // final check to make sure we are really really within the grid (since xmin < local_origin)
        if(i >=0 && i < cctk_lsh[0] && j >= 0 && j < cctk_lsh[1] && k >= 0 && k < cctk_lsh[2]) {
          const CCTK_INT idx = CCTK_GFINDEX3D(cctkGH, i,j,k);
          // actual payload code is here
          rho[idx] = part[i_rho];
          eps[idx] = part[i_eps];
          // NB: press is missing
          velx[idx] = part[i_vx];
          vely[idx] = part[i_vy];
          velz[idx] = part[i_vz];
        }
      }
    }
  }
}

extern "C"
void MakeParticles_SelectHydroVariablesBC(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // apply (symmetry) boundary condition to fill in symmetry points
  const char *groups[] = {"HydroBase::rho", "HydroBase::eps", "HydroBase::vel"}; // what about w_lorentz?
  CCTK_INT ierr = 0;
  for(int g = 0 ; g < DIM(groups) ; g++) {
    ierr += Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1, groups[g], "None");
  }
}

extern "C"
void MakeParticles_SetMetric(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL xmin = x[0] - CCTK_DELTA_SPACE(0)/2.;
  const CCTK_REAL ymin = y[0] - CCTK_DELTA_SPACE(1)/2.;
  const CCTK_REAL zmin = z[0] - CCTK_DELTA_SPACE(2)/2.;

  const CCTK_REAL xmax = x[CCTK_GFINDEX3D(cctkGH, cctk_lsh[0]-1, cctk_lsh[1]-1, cctk_lsh[2]-1)] + CCTK_DELTA_SPACE(0)/2.; 
  const CCTK_REAL ymax = y[CCTK_GFINDEX3D(cctkGH, cctk_lsh[0]-1, cctk_lsh[1]-1, cctk_lsh[2]-1)] + CCTK_DELTA_SPACE(1)/2.; 
  const CCTK_REAL zmax = z[CCTK_GFINDEX3D(cctkGH, cctk_lsh[0]-1, cctk_lsh[1]-1, cctk_lsh[2]-1)] + CCTK_DELTA_SPACE(2)/2.; 

  const CCTK_REAL dxfine = cctk_delta_space[0] / Carpet::spacereffacts[Carpet::maxreflevels-1][0];
  const CCTK_REAL dyfine = cctk_delta_space[0] / Carpet::spacereffacts[Carpet::maxreflevels-1][1];
  const CCTK_REAL dzfine = cctk_delta_space[0] / Carpet::spacereffacts[Carpet::maxreflevels-1][2];

  /* parse parameters into numbers */
  enum info_type_t info_type;               // parameter encoding
  info_type = info_invalid;
  for (int t = 0 ; t < DIM(info) ; t++) {
    if ( CCTK_EQUALS(information, info[t].name) ) {
      info_type = info[t].type;
      break;
    }
  }
  assert(info_type != info_invalid);
  
  assert(info_type == info_fallback);
  enum {i_x ,i_y ,i_z ,i_vx ,i_vy ,i_vz ,i_gxx ,i_gxy ,i_gxz ,i_gyy ,i_gyz ,i_gzz ,i_betax ,i_betay ,i_betaz ,i_alpha ,i_mass ,i_energy_per_unit_mass ,i_rho ,i_eps};

  int dp = info[info_type].nvals;
#pragma omp parallel for schedule(dynamic)
  for(int p = 0 ; p < int(particle_data.size()) ; p+=dp) {
    CCTK_REAL *part = &particle_data[p];
    CCTK_REAL dx, dy, dz;
    if(part[i_x] > xmin && part[i_x] < xmax &&
       part[i_y] > ymin && part[i_y] < ymax &&
       part[i_z] > zmin && part[i_z] < zmax) {
      // check that the point lies on one of our grid points
      dx = fmod(part[i_x] - CCTK_ORIGIN_SPACE(0), CCTK_DELTA_SPACE(0));
      if(dx > CCTK_DELTA_SPACE(0)/2)
        dx -= CCTK_DELTA_SPACE(0);
      dy = fmod(part[i_y] - CCTK_ORIGIN_SPACE(1), CCTK_DELTA_SPACE(1));
      if(dy > CCTK_DELTA_SPACE(1)/2)
        dy -= CCTK_DELTA_SPACE(1);
      dz = fmod(part[i_z] - CCTK_ORIGIN_SPACE(2), CCTK_DELTA_SPACE(2));
      if(dz > CCTK_DELTA_SPACE(2)/2)
        dz -= CCTK_DELTA_SPACE(2);
      if(fabs(dx) < dxfine/2 && fabs(dy) < dyfine/2 && fabs(dz) < dzfine/2) {
        const CCTK_INT i = CCTK_INT(round((part[i_x] - CCTK_ORIGIN_SPACE(0)) / CCTK_DELTA_SPACE(0))) - cctk_lbnd[0];
        const CCTK_INT j = CCTK_INT(round((part[i_y] - CCTK_ORIGIN_SPACE(1)) / CCTK_DELTA_SPACE(1))) - cctk_lbnd[1];
        const CCTK_INT k = CCTK_INT(round((part[i_z] - CCTK_ORIGIN_SPACE(2)) / CCTK_DELTA_SPACE(2))) - cctk_lbnd[2];
        // final check to make sure we are really really within the grid (since xmin < local_origin)
        if(i >=0 && i < cctk_lsh[0] && j >= 0 && j < cctk_lsh[1] && k >= 0 && k < cctk_lsh[2]) {
          const CCTK_INT idx = CCTK_GFINDEX3D(cctkGH, i,j,k);
          // actual payload code is here
          gxx[idx] = part[i_gxx];
          gxy[idx] = part[i_gxy];
          gxz[idx] = part[i_gxz];
          gyy[idx] = part[i_gyy];
          gyz[idx] = part[i_gyz];
          gzz[idx] = part[i_gzz];
        }
      }
    }
  }
}

extern "C"
void MakeParticles_SelectMetricBC(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // apply (symmetry) boundary condition to fill in symmetry points
  const char *groups[] = {"ADMBase::metric"};
  CCTK_INT ierr = 0;
  for(int g = 0 ; g < DIM(groups) ; g++) {
    ierr += Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1, groups[g], "None");
  }
}

extern "C"
void MakeParticles_SetLapse(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL xmin = x[0] - CCTK_DELTA_SPACE(0)/2.;
  const CCTK_REAL ymin = y[0] - CCTK_DELTA_SPACE(1)/2.;
  const CCTK_REAL zmin = z[0] - CCTK_DELTA_SPACE(2)/2.;

  const CCTK_REAL xmax = x[CCTK_GFINDEX3D(cctkGH, cctk_lsh[0]-1, cctk_lsh[1]-1, cctk_lsh[2]-1)] + CCTK_DELTA_SPACE(0)/2.; 
  const CCTK_REAL ymax = y[CCTK_GFINDEX3D(cctkGH, cctk_lsh[0]-1, cctk_lsh[1]-1, cctk_lsh[2]-1)] + CCTK_DELTA_SPACE(1)/2.; 
  const CCTK_REAL zmax = z[CCTK_GFINDEX3D(cctkGH, cctk_lsh[0]-1, cctk_lsh[1]-1, cctk_lsh[2]-1)] + CCTK_DELTA_SPACE(2)/2.; 

  const CCTK_REAL dxfine = cctk_delta_space[0] / Carpet::spacereffacts[Carpet::maxreflevels-1][0];
  const CCTK_REAL dyfine = cctk_delta_space[0] / Carpet::spacereffacts[Carpet::maxreflevels-1][1];
  const CCTK_REAL dzfine = cctk_delta_space[0] / Carpet::spacereffacts[Carpet::maxreflevels-1][2];

  /* parse parameters into numbers */
  enum info_type_t info_type;               // parameter encoding
  info_type = info_invalid;
  for (int t = 0 ; t < DIM(info) ; t++) {
    if ( CCTK_EQUALS(information, info[t].name) ) {
      info_type = info[t].type;
      break;
    }
  }
  assert(info_type != info_invalid);
  
  assert(info_type == info_fallback);
  enum {i_x ,i_y ,i_z ,i_vx ,i_vy ,i_vz ,i_gxx ,i_gxy ,i_gxz ,i_gyy ,i_gyz ,i_gzz ,i_betax ,i_betay ,i_betaz ,i_alpha ,i_mass ,i_energy_per_unit_mass ,i_rho ,i_eps};


  int dp=info[info_type].nvals;
#pragma omp parallel for schedule(dynamic)
  for(int p = 0 ; p < int(particle_data.size()) ; p+=dp) {
    CCTK_REAL *part = &particle_data[p];
    CCTK_REAL dx, dy, dz;
    if(part[i_x] > xmin && part[i_x] < xmax &&
       part[i_y] > ymin && part[i_y] < ymax &&
       part[i_z] > zmin && part[i_z] < zmax) {
      // check that the point lies on one of our grid points
      dx = fmod(part[i_x] - CCTK_ORIGIN_SPACE(0), CCTK_DELTA_SPACE(0));
      if(dx > CCTK_DELTA_SPACE(0)/2)
        dx -= CCTK_DELTA_SPACE(0);
      dy = fmod(part[i_y] - CCTK_ORIGIN_SPACE(1), CCTK_DELTA_SPACE(1));
      if(dy > CCTK_DELTA_SPACE(1)/2)
        dy -= CCTK_DELTA_SPACE(1);
      dz = fmod(part[i_z] - CCTK_ORIGIN_SPACE(2), CCTK_DELTA_SPACE(2));
      if(dz > CCTK_DELTA_SPACE(2)/2)
        dz -= CCTK_DELTA_SPACE(2);
      if(fabs(dx) < dxfine/2 && fabs(dy) < dyfine/2 && fabs(dz) < dzfine/2) {
        const CCTK_INT i = CCTK_INT(round((part[i_x] - CCTK_ORIGIN_SPACE(0)) / CCTK_DELTA_SPACE(0))) - cctk_lbnd[0];
        const CCTK_INT j = CCTK_INT(round((part[i_y] - CCTK_ORIGIN_SPACE(1)) / CCTK_DELTA_SPACE(1))) - cctk_lbnd[1];
        const CCTK_INT k = CCTK_INT(round((part[i_z] - CCTK_ORIGIN_SPACE(2)) / CCTK_DELTA_SPACE(2))) - cctk_lbnd[2];
        // final check to make sure we are really really within the grid (since xmin < local_origin)
        if(i >=0 && i < cctk_lsh[0] && j >= 0 && j < cctk_lsh[1] && k >= 0 && k < cctk_lsh[2]) {
          const CCTK_INT idx = CCTK_GFINDEX3D(cctkGH, i,j,k);
          // actual payload code is here
          alp[idx] = part[i_alpha];
        }
      }
    }
  }
}

extern "C"
void MakeParticles_SelectLapseBC(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // apply (symmetry) boundary condition to fill in symmetry points
  const char *groups[] = {"ADMBase::lapse"};
  CCTK_INT ierr = 0;
  for(int g = 0 ; g < DIM(groups) ; g++) {
    ierr += Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1, groups[g], "None");
  }
}

extern "C"
void MakeParticles_SetShift(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL xmin = x[0] - CCTK_DELTA_SPACE(0)/2.;
  const CCTK_REAL ymin = y[0] - CCTK_DELTA_SPACE(1)/2.;
  const CCTK_REAL zmin = z[0] - CCTK_DELTA_SPACE(2)/2.;

  const CCTK_REAL xmax = x[CCTK_GFINDEX3D(cctkGH, cctk_lsh[0]-1, cctk_lsh[1]-1, cctk_lsh[2]-1)] + CCTK_DELTA_SPACE(0)/2.; 
  const CCTK_REAL ymax = y[CCTK_GFINDEX3D(cctkGH, cctk_lsh[0]-1, cctk_lsh[1]-1, cctk_lsh[2]-1)] + CCTK_DELTA_SPACE(1)/2.; 
  const CCTK_REAL zmax = z[CCTK_GFINDEX3D(cctkGH, cctk_lsh[0]-1, cctk_lsh[1]-1, cctk_lsh[2]-1)] + CCTK_DELTA_SPACE(2)/2.; 

  const CCTK_REAL dxfine = cctk_delta_space[0] / Carpet::spacereffacts[Carpet::maxreflevels-1][0];
  const CCTK_REAL dyfine = cctk_delta_space[0] / Carpet::spacereffacts[Carpet::maxreflevels-1][1];
  const CCTK_REAL dzfine = cctk_delta_space[0] / Carpet::spacereffacts[Carpet::maxreflevels-1][2];

  /* parse parameters into numbers */
  enum info_type_t info_type;               // parameter encoding
  info_type = info_invalid;
  for (int t = 0 ; t < DIM(info) ; t++) {
    if ( CCTK_EQUALS(information, info[t].name) ) {
      info_type = info[t].type;
      break;
    }
  }
  assert(info_type != info_invalid);
  
  assert(info_type == info_fallback);
  enum {i_x ,i_y ,i_z ,i_vx ,i_vy ,i_vz ,i_gxx ,i_gxy ,i_gxz ,i_gyy ,i_gyz ,i_gzz ,i_betax ,i_betay ,i_betaz ,i_alpha ,i_mass ,i_energy_per_unit_mass ,i_rho ,i_eps};

  int dp=info[info_type].nvals;
#pragma omp parallel for schedule(dynamic)
  for(int p = 0 ; p < int(particle_data.size()) ; p+=dp) {
    CCTK_REAL *part = &particle_data[p];
    CCTK_REAL dx, dy, dz;
    if(part[i_x] > xmin && part[i_x] < xmax &&
       part[i_y] > ymin && part[i_y] < ymax &&
       part[i_z] > zmin && part[i_z] < zmax) {
      // check that the point lies on one of our grid points
      dx = fmod(part[i_x] - CCTK_ORIGIN_SPACE(0), CCTK_DELTA_SPACE(0));
      if(dx > CCTK_DELTA_SPACE(0)/2)
        dx -= CCTK_DELTA_SPACE(0);
      dy = fmod(part[i_y] - CCTK_ORIGIN_SPACE(1), CCTK_DELTA_SPACE(1));
      if(dy > CCTK_DELTA_SPACE(1)/2)
        dy -= CCTK_DELTA_SPACE(1);
      dz = fmod(part[i_z] - CCTK_ORIGIN_SPACE(2), CCTK_DELTA_SPACE(2));
      if(dz > CCTK_DELTA_SPACE(2)/2)
        dz -= CCTK_DELTA_SPACE(2);
      if(fabs(dx) < dxfine/2 && fabs(dy) < dyfine/2 && fabs(dz) < dzfine/2) {
        const CCTK_INT i = CCTK_INT(round((part[i_x] - CCTK_ORIGIN_SPACE(0)) / CCTK_DELTA_SPACE(0))) - cctk_lbnd[0];
        const CCTK_INT j = CCTK_INT(round((part[i_y] - CCTK_ORIGIN_SPACE(1)) / CCTK_DELTA_SPACE(1))) - cctk_lbnd[1];
        const CCTK_INT k = CCTK_INT(round((part[i_z] - CCTK_ORIGIN_SPACE(2)) / CCTK_DELTA_SPACE(2))) - cctk_lbnd[2];
        // final check to make sure we are really really within the grid (since xmin < local_origin)
        if(i >=0 && i < cctk_lsh[0] && j >= 0 && j < cctk_lsh[1] && k >= 0 && k < cctk_lsh[2]) {
          const CCTK_INT idx = CCTK_GFINDEX3D(cctkGH, i,j,k);
          // actual payload code is here
          betax[idx] = part[i_betax];
          betay[idx] = part[i_betay];
          betaz[idx] = part[i_betaz];
        }
      }
    }
  }
}

extern "C"
void MakeParticles_SelectShiftBC(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // apply (symmetry) boundary condition to fill in symmetry points
  const char *groups[] = {"ADMBase::shift"};
  CCTK_INT ierr = 0;
  for(int g = 0 ; g < DIM(groups) ; g++) {
    ierr += Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1, groups[g], "None");
  }
}

extern "C"
void MakeParticles_Shrink(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // we need some temp storage to avoid a dependence on the order in which we
  // process the grid.
  try {
    vector<bool> remove_cell(size_t(cctk_lsh[2]*cctk_lsh[1]*cctk_lsh[0]), false);

#pragma omp parallel for schedule(dynamic)
    for(CCTK_INT k = shrink_by_cells ; k < cctk_lsh[2]-shrink_by_cells ; k++) {
      for(CCTK_INT j = shrink_by_cells ; j < cctk_lsh[1]-shrink_by_cells ; j++) {
        for(CCTK_INT i = shrink_by_cells ; i < cctk_lsh[0]-shrink_by_cells ; i++) {
          CCTK_INT idx = CCTK_GFINDEX3D(cctkGH, i,j,k);
          if(rho_min < rho[idx] && rho[idx] < shrink_threshold) {
            bool atmo_seen = false;
            for(CCTK_INT kk = -shrink_by_cells ; kk <= shrink_by_cells && !atmo_seen ; kk++) {
              for(CCTK_INT jj = -shrink_by_cells ; jj <= shrink_by_cells && !atmo_seen ; jj++) {
                for(CCTK_INT ii = -shrink_by_cells ; ii <= shrink_by_cells && !atmo_seen ; ii++) {
                  CCTK_INT idxidx = CCTK_GFINDEX3D(cctkGH, i+ii,j+jj,k+kk);
                  if(rho[idxidx] < rho_min) {
                    remove_cell[idx] = true;
                    atmo_seen = true;
                  }
                }
              }
            }
          }
        }
      }
    }

#pragma omp parallel for
    for(CCTK_INT k = shrink_by_cells ; k < cctk_lsh[2]-shrink_by_cells ; k++) {
      for(CCTK_INT j = shrink_by_cells ; j < cctk_lsh[1]-shrink_by_cells ; j++) {
        for(CCTK_INT i = shrink_by_cells ; i < cctk_lsh[0]-shrink_by_cells ; i++) {
          CCTK_INT idx = CCTK_GFINDEX3D(cctkGH, i,j,k);
          if(remove_cell[idx]) {
            rho[idx] = rho_min;
          }
        }
      }
    }
  }
  catch (bad_alloc) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING, "Failed to allocate memory in %s.", __func__);
  }
}
