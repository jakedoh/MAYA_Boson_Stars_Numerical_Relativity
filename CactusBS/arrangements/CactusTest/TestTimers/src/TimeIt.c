#include <stdio.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusTest_TestTimers_TimeIt_c)

#define PRINTERROR                  \
if (ierr < 0)                       \
{                                   \
  printf("Error code is %d\n",ierr);\
}                 

/* Handles for the timers used */
static int index1;
static int index2;

int TimeIt_Startup(void);
void TimeIt(CCTK_ARGUMENTS);
void TimeIt_ShutDown (CCTK_ARGUMENTS);
static void PrintTimerData (const char *timername, const cTimerData *info);


int TimeIt_Startup(void)
{
  /* First timer times each iteration */
  printf("Creating named timer <timemystuff>\n");
  index1 = CCTK_TimerCreate("timemystuff");
  printf("Timer index is %d\n",index1);

  /* Second timer accumulates time */
  printf("Creating unamed timer\n");
  index2 = CCTK_TimerCreateI();
  printf("Timer index is %d\n",index2);
  return 0;
}
 
void TimeIt_ShutDown (CCTK_ARGUMENTS)
{
  int ierr;

  printf("Destroying named timer <timemystuff>\n");
  ierr = CCTK_TimerDestroy("timemystuff");
  PRINTERROR
  printf("Destroying unamed timer\n");
  ierr = CCTK_TimerDestroyI(index2);
  PRINTERROR
}
 

void TimeIt(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  int i;
  int ierr;
  cTimerData *info;
  DECLARE_CCTK_PARAMETERS

  printf("Resetting named timer\n");
  ierr = CCTK_TimerReset("timemystuff");
  PRINTERROR

  printf("Starting named timer\n");
  ierr = CCTK_TimerStart("timemystuff");
  PRINTERROR

  printf("Starting unamed timer\n");
  ierr = CCTK_TimerStartI(index2);
  PRINTERROR

  for (i=0; i<1000000; i++)
  {
    ierr = CCTK_Equals(TestParameter,"testit");
    ierr = CCTK_Equals(TestParameter,"testit");
    ierr = CCTK_Equals(TestParameter,"testit");
    ierr = CCTK_Equals(TestParameter,"testit");
    ierr = CCTK_Equals(TestParameter,"testit");
  } 

  printf("Stopping named timer\n");
  ierr = CCTK_TimerStop("timemystuff");
  PRINTERROR

  printf("Stopping unamed timer\n");
  ierr = CCTK_TimerStopI(index2);
  PRINTERROR

  /* Now print out the data */
  printf("Creating data structure for all clocks\n");
  info = CCTK_TimerCreateData();
 
  printf("Now fill it with measurements for named timer\n");
  ierr = CCTK_Timer("timemystuff",info);
  PRINTERROR

  PrintTimerData ("named", info);

  printf("\nDestroying data structure\n");
  ierr = CCTK_TimerDestroyData(info);
  PRINTERROR
  
  printf("Creating data structure for all clocks\n");
  info = CCTK_TimerCreateData();
 
  printf("\nNow fill it with measurements for unamed timer\n");
  ierr = CCTK_TimerI(index2,info);
  PRINTERROR

  PrintTimerData ("unnamed", info);
  
  printf("\nDestroying data structure\n");
  ierr = CCTK_TimerDestroyData(info);
  PRINTERROR
 
  printf("Creating data structure for all timers\n");
  info = CCTK_TimerCreateData ();
  for (i = 0; i < CCTK_NumTimers (); i++)
  {
    if (CCTK_TimerI (i, info) == 0)
    {
      PrintTimerData (CCTK_TimerName (i), info);
    }
  }

  printf("\nDestroying data structure\n");
  ierr = CCTK_TimerDestroyData(info);
  PRINTERROR
 
}


static void PrintTimerData (const char *timername, const cTimerData *info)
{
  int i;


  printf("\nResults for timer '%s'\n", timername);
  printf("-----------------------\n");
  for (i = 0; i < info->n_vals; i++)
  {
    switch (info->vals[i].type)
      {
      case val_int:
        printf("%s: %d %s\n", 
               info->vals[i].heading,info->vals[i].val.i, 
               info->vals[i].units);
        break;
        
      case val_long:
        printf("%s: %d %s\n", 
               info->vals[i].heading,(int) info->vals[i].val.l, 
               info->vals[i].units);
        break;
        
      case val_double:
        printf("%s: %.3f %s\n", 
               info->vals[i].heading,info->vals[i].val.d, 
               info->vals[i].units);
        break;
        
      default:
        CCTK_WARN(1, "Unknown data type for timer info");
        break;
      }
  }
}
