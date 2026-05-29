 /*@@
   @file      CRoutines.c
   @date      20th June 2001
   @author    Gabrielle Allen
   @desc
              Test stuff with strings
   @enddesc
   @version   $Id$
 @@*/

#include <stdio.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_FortranString.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusTest_TestStrings_CRoutines_c)


/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

void PrintOneString(int num1, 
                    const int num2,
                    const char *first);
void PrintTwoString(int num1,
                    const int num2,
                    const char *first, 
                    const char *second);
void PrintThreeString(int num1,
                      const int num2,
                      const char *first, 
                      const char *second, 
                      const char *third);


/********************************************************************
 *********************     Fortran Wrappers    **********************
 ********************************************************************/

void CCTK_FCALL CCTK_FNAME (PrintOneString)
     (int *num1, const int *num2, ONE_FORTSTRING_ARG);

void CCTK_FCALL CCTK_FNAME (PrintTwoString)
     (int *num1, const int *num2, TWO_FORTSTRING_ARG);

void CCTK_FCALL CCTK_FNAME (PrintThreeString)
     (int *num1, const int *num2, THREE_FORTSTRING_ARG);


/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/


/*@@
   @routine    PrintOneString
   @date       Wed 20th June 2001
   @author     Gabrielle Allen
   @desc
   Just print a string
   @enddesc
   @calls

   @var     num1
   @vdesc   passed in number
   @vtype   int
   @vio     in
   @endvar
   @var     num2
   @vdesc   passed in number
   @vtype   const int *
   @vio     in
   @endvar
   @var     first
   @vdesc   string to be printed
   @vtype   const char *
   @vio     in
   @endvar

   @returntype void
   @returndesc
   @endreturndesc
@@*/

void PrintOneString(int num1,
                    const int num2,
                    const char *first)
{
  printf("Numbers: %d %d\n",num1, num2);
  printf("PrintOneString: <%s>\n",first);
}

void CCTK_FCALL CCTK_FNAME (PrintOneString)
     (int *num1, const int *num2,
      ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(name);
  PrintOneString(*num1, *num2, name);
  free(name);
}


/*@@
   @routine    PrintTwoString
   @date       Wed 20th June 2001
   @author     Gabrielle Allen
   @desc
   Just print two string
   @enddesc
   @calls

   @var     num1
   @vdesc   passed in number
   @vtype   int
   @vio     in
   @endvar
   @var     num2
   @vdesc   passed in number
   @vtype   const int *
   @vio     in
   @endvar
   @var     first
   @vdesc   string to be printed
   @vtype   const char *
   @vio     in
   @endvar
   @var     second
   @vdesc   string to be printed
   @vtype   const char *
   @vio     in
   @endvar

   @returntype void
   @returndesc
   @endreturndesc
@@*/

void PrintTwoString(int num1, 
                    const int num2, 
                    const char *first, 
                    const char *second)
{
  printf("Numbers: %d %d\n",num1,num2);
  printf("PrintTwoString: <%s> <%s>\n",first,second);
}

void CCTK_FCALL CCTK_FNAME (PrintTwoString)
                           (int *num1, const int *num2,
                            TWO_FORTSTRING_ARG)
{
  TWO_FORTSTRING_CREATE(name1,name2);
  PrintTwoString(*num1,*num2,name1,name2);
  free(name1);
  free(name2);
}


/*@@
   @routine    PrintThreeString
   @date       Wed 20th June 2001
   @author     Gabrielle Allen
   @desc
   Just print three strings
   @enddesc
   @calls

   @var     num1
   @vdesc   passed in number
   @vtype   int
   @vio     in
   @endvar
   @var     num2
   @vdesc   passed in number
   @vtype   const int *
   @vio     in
   @endvar
   @var     first
   @vdesc   string to be printed
   @vtype   const char *
   @vio     in
   @endvar
   @var     second
   @vdesc   string to be printed
   @vtype   const char *
   @vio     in
   @endvar
   @var     third
   @vdesc   string to be printed
   @vtype   const char *
   @vio     in
   @endvar

   @returntype void
   @returndesc
   @endreturndesc
@@*/

void PrintThreeString(int num1, 
                      const int num2,
                      const char *first, 
                      const char *second, 
                      const char *third)
{
  printf("Numbers: %d %d\n",num1, num2);
  printf("PrintThreeString: <%s> <%s> <%s>\n",first,second,third);
}

void CCTK_FCALL CCTK_FNAME (PrintThreeString)
                           (int *num1, const int *num2, THREE_FORTSTRING_ARG)
{
  THREE_FORTSTRING_CREATE(name1,name2,name3);
  PrintThreeString(*num1,*num2,name1,name2,name3);
  free(name1);
  free(name2);
  free(name3);
}


