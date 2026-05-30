
#include <stdio.h> 
#include <string.h> 
#include <assert.h> 

#include "cctk.h" 
#include "cctk_Arguments.h" 
#include "cctk_Parameters.h" 
#include "cctk_Timers.h"

#define DEBUG

static int timer_handle = -1;
static CCTK_REAL last_timeofday;
static CCTK_REAL current_timeofday;
static CCTK_REAL last_coord_time;
static CCTK_REAL current_coord_time;
static CCTK_REAL period; /* in minutes */
static int last_times_valid = 0;
static int stats_valid = 0;

/* Initialise the timer  */
void RunStats_Init(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  timer_handle = CCTK_TimerCreate("Time this run");
  assert(timer_handle >= 0);
  CCTK_TimerStartI(timer_handle);
  last_times_valid = 0;
  stats_valid = 0;
  *cputime_used = 0;
  *speed = 0;
  *maxrss = 0;
  *majflt = 0;
  *arena = 0;
  *ordblks = 0;
  *hblks = 0;
  *hblkhd = 0;
  *uordblks = 0;
  *fordblks = 0;
  *keepcost = 0;

  //  printf("Initialising timer\n");
}

/* Output the statistics  */
void RunStats_Output(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const int buf_size = 1024;
  static int first_time = 1;
  const char *mode = first_time ? "w" : "a";
  char output_name[buf_size];
  CCTK_STRING *out_dir  = (CCTK_STRING *) CCTK_ParameterGet("out_dir", "IOUtil", NULL);
  FILE *fp = NULL;

  if (CCTK_MyProc(cctkGH) != 0) 
    return;

  if (!stats_valid)
    return;

  if (cctk_iteration % calc_every != 0)
    return;

  if (*out_dir == NULL)
  {
    CCTK_WARN(1,"Parameter IOUtil::out_dir not found");
    return;
  }

  if (strlen(*out_dir)+1 > buf_size)
  {
    CCTK_WARN(1,"Parameter IOUtil::out_dir is too long");
    return;
  }

  sprintf(output_name, "%s/runstats.asc", *out_dir);

  fp = fopen(output_name, mode);
  if (fp == NULL)
  {
    CCTK_WARN(1,"Cannot open output file for run statistics");
    return;
  }

  if (first_time)
    fprintf(fp, "# iteration, coord_time, wall_time, speed (hours^-1), period (minutes), cputime (cpu hours)\n");

  first_time = 0;

  fprintf(fp, "%d %f %f %f %f %f\n",
	  cctk_iteration, cctk_time, current_timeofday, *speed, period, *cputime_used);

  fclose(fp);
}

/* Calculate statistics  */
void RunStats_Calculate(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (cctk_iteration % calc_every != 0)
    return;

  assert(timer_handle >= 0);
  CCTK_TimerStopI(timer_handle);

  cTimerData *info = CCTK_TimerCreateData();
  if (info == NULL)
  {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
		"Could not create timer information structure");
    stats_valid = 0;
    last_times_valid = 0;
    return;
  }

  int ierr = CCTK_TimerI(timer_handle,info);
  int nprocs = CCTK_nProcs(cctkGH);


  if (ierr > 0)
  {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
		"Could not access timer information");
    stats_valid = 0;
    last_times_valid = 0;
    return;
  }

  const cTimerVal   *clock = CCTK_GetClockValue("gettimeofday",info);

  if (clock)
  {
    current_timeofday = CCTK_TimerClockSeconds( clock );
  }
  else
  {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
		"Could not access clock value from timer information");
    stats_valid = 0;
    last_times_valid = 0;
    return;
  }
  //  printf("current_timeofday == %f\n", current_timeofday);

  CCTK_TimerDestroyData(info); 

  current_coord_time = cctk_time;
  //  printf("current_coord_time == %f\n", current_coord_time);

  if (last_times_valid)
  {
    CCTK_REAL delta_timeofday = current_timeofday - last_timeofday;
    CCTK_REAL delta_coord_time = current_coord_time - last_coord_time;

    //    printf("delta_timeofday == %f\n", delta_timeofday);
    //    printf("delta_coord_time == %f\n", delta_coord_time);

    *speed = 3600 * delta_coord_time / delta_timeofday;
    //    printf("speed == %f\n", *speed);
    period = delta_timeofday / 60 / delta_coord_time;
    //    printf("period == %f\n", period);
    *cputime_used += nprocs * delta_timeofday / 3600.0;
    stats_valid = 1;
  }


  last_coord_time = current_coord_time;
  last_timeofday = current_timeofday;
  last_times_valid = 1;

  CCTK_TimerStartI(timer_handle);
}
