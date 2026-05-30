 /*@@
   @file      TestSchedule.c
   @date      Wed Aug 18 12:57:41 2004
   @author    Rion Dooley
   @desc 
   Test functions to see if everything is scheduled correctly.
   @enddesc 
   @version $Header$
 @@*/
#include <stdio.h>
#include <stdlib.h>
#include "cctk.h"

#include "cctk.h"

#include "cctk_Arguments.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusTest_TestSchedule_TestSchedule_c);

/********************************************************************
 *********************  Macro Definitions  **************************
 ********************************************************************/

/********************************************************************
 *********************  Local Data Types  ***************************
 ********************************************************************/

/********************************************************************
 *********************  Aliased Routine Prototypes  *****************
 ********************************************************************/

/********************************************************************
 *********************  Scheduled Routine Prototypes  ***************
 ********************************************************************/

int TestSchedule_CCTK_STARTUP(void);
int TestSchedule_CCTK_RECOVER_PARAMETERS(void);
void TestSchedule_CCTK_WRAGH(CCTK_ARGUMENTS);
void TestSchedule_CCTK_PARAMCHECK(CCTK_ARGUMENTS);
void TestSchedule_CCTK_BASEGRID(CCTK_ARGUMENTS);
void TestSchedule_CCTK_INITIAL(CCTK_ARGUMENTS);
void TestSchedule_CCTK_POSTINITIAL(CCTK_ARGUMENTS);
void TestSchedule_CCTK_RECOVER_VARIABLES(CCTK_ARGUMENTS);
void TestSchedule_CCTK_POST_RECOVER_VARIABLES(CCTK_ARGUMENTS);
void TestSchedule_CCTK_CPINITIAL(CCTK_ARGUMENTS);
void TestSchedule_CCTK_CHECKPOINT(CCTK_ARGUMENTS);
void TestSchedule_CCTK_PRESTEP(CCTK_ARGUMENTS);
void TestSchedule_CCTK_EVOL(CCTK_ARGUMENTS);
void TestSchedule_CCTK_POSTSTEP(CCTK_ARGUMENTS);
void TestSchedule_CCTK_POSTRESTRICT(CCTK_ARGUMENTS);
void TestSchedule_CCTK_POSTREGRID(CCTK_ARGUMENTS);
void TestSchedule_CCTK_ANALYSIS(CCTK_ARGUMENTS);
void TestSchedule_CCTK_TERMINATE(CCTK_ARGUMENTS);
int TestSchedule_CCTK_SHUTDOWN(void);

/********************************************************************
 *********************  Fortran Wrapper Prototypes  *****************
 ********************************************************************/

/********************************************************************
 *********************  Local Routine Prototypes  *******************
 ********************************************************************/

/********************************************************************
 *********************  Local Data  *********************************
 ********************************************************************/

/********************************************************************
 *********************  Aliased Routines  ***************************
 ********************************************************************/

/********************************************************************
 *********************  Scheduled Routines  *************************
 ********************************************************************/

 /*@@
   @routine    TestSchedule_CCTK_STARTUP
   @date       Wed Aug 18 12:57:41 2004
   @author     Rion Dooley
   @desc 
   Prints a message this function has has been invoked.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
 
 @@*/
int TestSchedule_CCTK_STARTUP(void)
{
  CCTK_INFO("executed in CCTK_STARTUP");

  return 0;

}

 /*@@
   @routine    TestSchedule_CCTK_STARTUP
   @date       Wed Aug 18 12:57:41 2004
   @author     Rion Dooley
   @desc 
   Prints a message this function has has been invoked.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
 
 @@*/
int TestSchedule_CCTK_RECOVER_PARAMETERS(void)
{
  CCTK_INFO("executed in CCTK_RECOVER_PARAMETERS");

  return 0;
}

 /*@@
   @routine    TestSchedule_CCTK_WRAGH
   @date       Wed Aug 18 12:57:41 2004
   @author     Rion Dooley
   @desc 
   Prints a message this function has has been invoked.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
 
 @@*/
void TestSchedule_CCTK_WRAGH(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;

  CCTK_INFO("executed in CCTK_WRAGH");
}

 /*@@
   @routine    TestSchedule_CCTK_PARAMCHECK
   @date       Wed Aug 18 12:57:41 2004
   @author     Rion Dooley
   @desc 
   Prints a message this function has has been invoked.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
 
 @@*/
void TestSchedule_CCTK_PARAMCHECK(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;

  CCTK_INFO("executed in CCTK_PARAMCHECK");
}

 /*@@
   @routine    TestSchedule_CCTK_BASEGRID
   @date       Wed Aug 18 12:57:41 2004
   @author     Rion Dooley
   @desc 
   Prints a message this function has has been invoked.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
 
 @@*/
void TestSchedule_CCTK_BASEGRID(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;

  CCTK_INFO("executed in CCTK_BASEGRID");
}

 /*@@
   @routine    TestSchedule_CCTK_INITIAL
   @date       Wed Aug 18 12:57:41 2004
   @author     Rion Dooley
   @desc 
   Prints a message this function has has been invoked.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
 
 @@*/
void TestSchedule_CCTK_INITIAL(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;

  CCTK_INFO("executed in CCTK_INITIAL");
}

 /*@@
   @routine    TestSchedule_CCTK_POSTINITIAL
   @date       Wed Aug 18 12:57:41 2004
   @author     Rion Dooley
   @desc 
   Prints a message this function has has been invoked.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
 
 @@*/
void TestSchedule_CCTK_POSTINITIAL(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;

  CCTK_INFO("executed in CCTK_POSTINITIAL");
}

 /*@@
   @routine    TestSchedule_CCTK_RECOVER_VARIABLES
   @date       Wed Aug 18 12:57:41 2004
   @author     Rion Dooley
   @desc 
   Prints a message this function has has been invoked.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
 
 @@*/
void TestSchedule_CCTK_RECOVER_VARIABLES(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;

  CCTK_INFO("executed in CCTK_RECOVER_VARIABLES");
}

 /*@@
   @routine    TestSchedule_CCTK_POST_RECOVER_VARIABLES
   @date       Wed Aug 18 12:57:41 2004
   @author     Rion Dooley
   @desc 
   Prints a message this function has has been invoked.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
 
 @@*/
void TestSchedule_CCTK_POST_RECOVER_VARIABLES(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;

  CCTK_INFO("executed in CCTK_POST_RECOVER_VARIABLES");
}

 /*@@
   @routine    TestSchedule_CCTK_CPINITIAL
   @date       Wed Aug 18 12:57:41 2004
   @author     Rion Dooley
   @desc 
   Prints a message this function has has been invoked.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
 
 @@*/
void TestSchedule_CCTK_CPINITIAL(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;

  CCTK_INFO("executed in CCTK_CPINITIAL");
}

 /*@@
   @routine    TestSchedule_CCTK_CHECKPOINT
   @date       Wed Aug 18 12:57:41 2004
   @author     Rion Dooley
   @desc 
   Prints a message this function has has been invoked.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
 
 @@*/
void TestSchedule_CCTK_CHECKPOINT(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;

  CCTK_INFO("executed in CCTK_CHECKPOINT");
}

 /*@@
   @routine    TestSchedule_CCTK_PRESTEP
   @date       Wed Aug 18 12:57:41 2004
   @author     Rion Dooley
   @desc 
   Prints a message this function has has been invoked.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
 
 @@*/
void TestSchedule_CCTK_PRESTEP(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;

  CCTK_INFO("executed in CCTK_PRESTEP");
}

 /*@@
   @routine    TestSchedule_CCTK_EVOL
   @date       Wed Aug 18 12:57:41 2004
   @author     Rion Dooley
   @desc 
   Prints a message this function has has been invoked.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
 
 @@*/
void TestSchedule_CCTK_EVOL(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;

  CCTK_INFO("executed in CCTK_EVOL");
}

 /*@@
   @routine    TestSchedule_CCTK_POSTSTEP
   @date       Wed Aug 18 12:57:41 2004
   @author     Rion Dooley
   @desc 
   Prints a message this function has has been invoked.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
 
 @@*/
void TestSchedule_CCTK_POSTSTEP(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;

  CCTK_INFO("executed in CCTK_POSTSTEP");
}

 /*@@
   @routine    TestSchedule_CCTK_POSTRESTRICT
   @date       Wed Aug 18 12:57:41 2004
   @author     Rion Dooley
   @desc 
   Prints a message this function has has been invoked.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
 
 @@*/
void TestSchedule_CCTK_POSTRESTRICT(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;

  CCTK_INFO("executed in CCTK_POSTRESTRICT");
}

 /*@@
   @routine    TestSchedule_CCTK_POSTREGRID
   @date       Wed Aug 18 12:57:41 2004
   @author     Rion Dooley
   @desc 
   Prints a message this function has has been invoked.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
 
 @@*/
void TestSchedule_CCTK_POSTREGRID(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;

  CCTK_INFO("executed in CCTK_POSTREGRID");
}

 /*@@
   @routine    TestSchedule_CCTK_ANALYSIS
   @date       Wed Aug 18 12:57:41 2004
   @author     Rion Dooley
   @desc 
   Prints a message this function has has been invoked.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
 
 @@*/
void TestSchedule_CCTK_ANALYSIS(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;

  CCTK_INFO("executed in CCTK_ANALYSIS");
}

 /*@@
   @routine    TestSchedule_CCTK_TERMINATE
   @date       Wed Aug 18 12:57:41 2004
   @author     Rion Dooley
   @desc 
   Prints a message this function has has been invoked.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
 
 @@*/
void TestSchedule_CCTK_TERMINATE(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;

  CCTK_INFO("executed in CCTK_TERMINATE");
}

 /*@@
   @routine    TestSchedule_CCTK_SHUTDOWN
   @date       Wed Aug 18 12:57:41 2004
   @author     Rion Dooley
   @desc 
   Prints a message this function has has been invoked.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
 
 @@*/
int TestSchedule_CCTK_SHUTDOWN(void)
{
  CCTK_INFO("executed in CCTK_SHUTDOWN");
  return 0;
}

/********************************************************************
 *********************  Other External Routines  ********************
 ********************************************************************/

/********************************************************************
 *********************  Local Routines  *****************************
 ********************************************************************/

