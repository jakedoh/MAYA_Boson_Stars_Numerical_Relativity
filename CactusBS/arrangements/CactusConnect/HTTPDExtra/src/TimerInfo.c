 /*@@
   @file      TimerInfo.c
   @date      Wed Dec 16 2004
   @author    Andre Werthmann
   @desc
              Pages about CCTK timers
   @enddesc
   @version   $Header$
 @@*/

#include "cctk.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "util_String.h"

#include "httpextra_HostNames.h"
#include "http_Content.h"

/* SW Temporary, while testing the SString module*/
#include "CactusConnect/HTTPD/src/SString.h"
#include "CactusConnect/HTTPD/src/SString_Namespace.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusConnect_HTTPDExtra_TimerInfo_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

#define DECBUFSIZE 64
#define NUMTIMERTYPES 19
/* static field with all know schedule bin names */
static char* TI_ttypes[NUMTIMERTYPES] =
{
        "CCTK_STARTUP",
        "CCTK_WRAGH",
        "CCTK_PARAMCHECK",
        "CCTK_BASEGRID",
        "CCTK_INITIAL",
        "CCTK_POSTREGRID",
        "CCTK_POSTRESTRICTINITIAL",
        "CCTK_POSTINITIAL",
        "CCTK_POSTSTEP",
        "CCTK_RECOVER_VARIABLES",
        "CCTK_POST_RECOVER_VARIABLES",
        "CCTK_CPINITIAL",
        "CCTK_ANALYSIS",
        "CCTK_PRESTEP",
        "CCTK_EVOL",
        "CCTK_POSTRESTRICT",
        "CCTK_CHECKPOINT",
        "CCTK_TERMINATE",
        "CCTK_SHUTDOWN"
};
/* actual number of not resolved schedule bin names */
static int TI_nettypes=0;
/* actual number of allocated entries in the field */
static int TI_maxnettypes=0;
/* dynamic field that holds the names of not resolved schedule bin names */
static char** TI_ettypes=NULL;
/* dynamic field that holds the values of clocks for timers for nicer html output later */
static double** TI_oldtimes=NULL;
static int TI_oldntimers=0;
static int TI_oldnclocks=0;

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

static int TimerInfoPage(const cGH *cctkGH, httpRequest *request, void *data);

/* local list implementation, comments below */
static void TI_initExtraTimerTypeField(void);
static int TI_isElement(char* selem);
static int TI_insertElement(char* selem);
static void TI_deleteExtraTimerTypeField(void);

/* helper functions, comments below */
static int TI_isTimerType(const char* tName, char* tType);
static char* TI_thornName(const char* tn);
static char* TI_funcName(const char* tn);
static char* TI_schedName(const char* tn);

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    HTTPDExtra_RegisterTimerInfoPages
   @date       Thu Dec 09 18:12:43 2004
   @author     Andre Werthmann
   @desc
               Httpd utils registration routine.
   @enddesc
@@*/
int HTTPDExtra_RegisterTimerInfoPages(void);
int HTTPDExtra_RegisterTimerInfoPages(void)
{
  /* Register the group info page. */
  HTTP_RegisterPage("/TimerInfo", TimerInfoPage, NULL);

  HTTP_ContentLink("/TimerInfo/index.html", "Timer Information",
                   "CCTK Timer information",
                   HTTP_QUICKLINK);
  return 0;
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/



/******************************************************************************
 *************************** Groups Page **************************************
 ******************************************************************************/

 /*@@
   @routine    TimerInfoPage
   @date       Thu Dec 09 18:12:43 2004
   @author     Andre Werthmann
   @desc
               Displays the cctk timer description page.
   @enddesc
@@*/
static int TimerInfoPage(const cGH *cctkGH, httpRequest *request, void *data)
{
  int retval = 0;
  int ntimers = 0, nt = 0;
  int nclocks = 0, nc = 0;
  int i, j;
  char *sbuff;
  String *message = String_New();
  cTimerData *tdata;

  HTTP_SendOKHeader( request );

  HTTP_SetDoctype( message );
  HTTP_SendString(request, message);

  /* Start the page */
  HTTP_Send(request, "<html>\n<head>\n");
  HTTP_Send(request, "<title>Cactus Simulation CCTK Timer Information</title>\n");

  HTTP_SetHeadInfo( message);
  HTTP_SendString(request, message );

  HTTP_Send(request, "<style type=\"text/css\">\n");
  HTTP_Send(request, "  th, td { padding-left: 1em; padding-right: 1em; }\n");
  HTTP_Send(request, "  td.name { text-align: left; }\n");
  HTTP_Send(request, "  td.number { text-align: center; }\n");
  HTTP_Send(request, "</style>\n");

  HTTP_Send(request, "</head>\n<body>\n");

  /* HTTP_Write out the header part. */
  HTTP_SetContentHeaderString(cctkGH, 0, message, NULL);
  retval = HTTP_SendString(request, message);

  HTTP_Send(request, "<h1>CCTK Timer Information</h1>\n");

  HTTP_Send(request, "<center>(changed values since last refresh are in bold characters)</center>\n<br>");

  /* init timer data */
  tdata=CCTK_TimerCreateData();
  ntimers = CCTK_NumTimers();
  nclocks = tdata->n_vals;

  /** init the timer clocks buffer
   * for every clock and every timer store an old value to
   * enhance html output of changed values
   * format is: TI_oldtimes[ClockNumber][TimerNumber]
   * memory is allocated dynamically if there are new timers and/or clocks
   * during computation
   */
  if (nclocks != TI_oldnclocks)
  {
    TI_oldtimes=realloc(TI_oldtimes, nclocks*sizeof(double*));
    for (i=TI_oldnclocks; i<nclocks; i++)
    {
      TI_oldtimes[i]=NULL;
    }
    TI_oldnclocks=nclocks;
  }
  if (ntimers != TI_oldntimers)
  {
    for (i=0; i<nclocks; i++)
    {
      TI_oldtimes[i]=realloc(TI_oldtimes[i], ntimers*sizeof(double));
      for (j=TI_oldntimers; j<ntimers; j++)
      {
        TI_oldtimes[i][j]=0.0;
      }
    }
    TI_oldntimers=ntimers;
  }

  /* start of first table */
  HTTP_Send(request, "<h2>Timers which are associated with schedule bins</h2>\n");

  /* set table header */
  SetToCString(message, "<div class=\"centered\">\n<table frame=\"box\" rules=\"all\">\n<tr>");

  /* timer type - table header */
  ConcatCString(message, "<th>Schedule Bin</th>");

  /* thorn name - table header */
  ConcatCString(message, "<th>Thorn Name</th>");

  /* description name - table header */
  ConcatCString(message, "<th>Description</th>");

  /* clock names - table header */
  CCTK_TimerI(0, tdata);
  for (i=0; i<nclocks; i++)
  {
    /** assume every timer has the same number of clocks
     *  take timer nr. 0 to build the table header
     */
    ConcatCString(message, "<th>");
    ConcatCString(message, tdata->vals[i].heading);
    ConcatCString(message, " (");
    ConcatCString(message, tdata->vals[i].units);
    ConcatCString(message, ")</th>");
  }
  ConcatCString(message, "</tr>");
  HTTP_SendString(request, message);

  /* parse schedule bins */
  for (i=0; i<NUMTIMERTYPES; i++)
  {
    for (nt=0; nt<ntimers; nt++)
    {
      if (TI_isTimerType(CCTK_TimerName(nt), TI_ttypes[i]))
      {
        /* fill out timer info for timer number nr into tdata */
        CCTK_TimerI(nt, tdata);
        /* start table row */
        SetToCString(message, "<tr>\n");

        /* timer type, first column */
        ConcatCString(message, "<td class=\"type\">\n");
        sbuff=TI_schedName(TI_ttypes[i]);
        ConcatCString(message, sbuff);
        free(sbuff);
        ConcatCString(message, "</td>\n");

        /* thorn name, second column */
        ConcatCString(message, "<td class=\"name\">");
        sbuff=TI_thornName(CCTK_TimerName(nt));
        ConcatCString(message, sbuff);
        free(sbuff);
        ConcatCString(message, "</td>\n");

        /* description */
        ConcatCString(message, "<td class=\"description\">");
        sbuff=TI_funcName(CCTK_TimerName(nt));
        ConcatCString(message, sbuff);
        free(sbuff);
        ConcatCString(message, "</td>\n");

        /* timer clocks, fourth... (fourth+number clocks) column */
        for (nc=0; nc<nclocks; nc++)
        {
          ConcatCString(message, "<td class=\"");
          ConcatCString(message, tdata->vals[nc].heading);
          ConcatCString(message, "\">");
          switch (tdata->vals[nc].type)
          {
            case val_int:
              ConcatDecimal(message, tdata->vals[nc].val.i);
              break;
            case val_long:
              ConcatDecimal(message, tdata->vals[nc].val.l);
              break;
            case val_double:
              /** if the new value is different from the old one
               * (changed) print it in bold characters
               */
              if ( TI_oldtimes[nc][nt] != tdata->vals[nc].val.d )
              {
                ConcatCString(message, "<span style=\"font-weight:bold\">");
                ConcatDouble(message, tdata->vals[nc].val.d);
                ConcatCString(message, "</span>");
              }
              else
              {
                ConcatDouble(message, tdata->vals[nc].val.d);
              }
              /* old value = new value */
              TI_oldtimes[nc][nt]=tdata->vals[nc].val.d;
              break;
            default:
              ConcatCString(message, "unknown value type");
              break;
          }
          ConcatCString(message, "</td>\n");
        }
        /* end table row */
        ConcatCString(message, "</tr>\n");
        HTTP_SendString(request, message);
      }
    }
  }
  retval = HTTP_Send(request, "</table>\n</div>\n<br>\n");

  HTTP_Send(request, "<h2>Timers which are associated with schedule groups</h2>\n");

  /* init the (char*) timer list of schedule bins which couldn't be associated */
  TI_initExtraTimerTypeField();

  /* start of second table */
  /* set table header */
  SetToCString(message, "<div class=\"centered\">\n<table frame=\"box\" rules=\"all\">\n<tr>");

  /* timer type - table header */
  ConcatCString(message, "<th>Schedule Group</th>");

  /* thorn name - table header */
  ConcatCString(message, "<th>Thorn Name</th>");

  /* description name - table header */
  ConcatCString(message, "<th>Description</th>");

  /* clock names - table header */
  CCTK_TimerI(0, tdata);
  for (i=0; i<nclocks; i++)
  {
    /** assume every timer has the same number of clocks
     *  take timer nr. 0 to build the table header
     */
    ConcatCString(message, "<th>");
    ConcatCString(message, tdata->vals[i].heading);
    ConcatCString(message, " (");
    ConcatCString(message, tdata->vals[i].units);
    ConcatCString(message, ")</th>");
  }
  ConcatCString(message, "</tr>");
  HTTP_SendString(request, message);

  /* second table content */
  for (i=0; i<TI_nettypes; i++)
  {
    for (nt=0; nt<ntimers; nt++)
    {
      if (TI_isTimerType(CCTK_TimerName(nt), TI_ettypes[i]))
      {
        /* fill out timer info for timer number nr into tdata */
        CCTK_TimerI(nt, tdata);
        /* start table row */
        SetToCString(message, "<tr>\n");

        /* timer type, first column */
        ConcatCString(message, "<td class=\"type\">\n");
        sbuff=TI_schedName(TI_ettypes[i]);
        ConcatCString(message, sbuff);
        free(sbuff);
        ConcatCString(message, "</td>\n");

        /* thorn name, second column */
        ConcatCString(message, "<td class=\"name\">");
        sbuff=TI_thornName(CCTK_TimerName(nt));
        ConcatCString(message, sbuff);
        free(sbuff);
        ConcatCString(message, "</td>\n");

        /* description */
        ConcatCString(message, "<td class=\"description\">");
        sbuff=TI_funcName(CCTK_TimerName(nt));
        ConcatCString(message, sbuff);
        free(sbuff);
        ConcatCString(message, "</td>\n");

        /* timer clocks, fourth... (fourth+number clocks) column */
        for (nc=0; nc<nclocks; nc++)
        {
          ConcatCString(message, "<td class=\"");
          ConcatCString(message, tdata->vals[nc].heading);
          ConcatCString(message, "\">");

          switch (tdata->vals[nc].type)
          {
            case val_int:
              ConcatDecimal(message, tdata->vals[nc].val.i);
              break;
            case val_long:
              ConcatDecimal(message, tdata->vals[nc].val.l);
              break;
            case val_double:
              /** if the new value is different from the old one
               * (changed) print it in bold characters
                */
              if ( TI_oldtimes[nc][nt] != tdata->vals[nc].val.d )
              {
                ConcatCString(message, "<span style=\"font-weight:bold\">");
                ConcatDouble(message, tdata->vals[nc].val.d);
                ConcatCString(message, "</span>");
              }
              else
              {
                ConcatDouble(message, tdata->vals[nc].val.d);
              }
              /* old value = new value */
              TI_oldtimes[nc][nt]=tdata->vals[nc].val.d;

              break;
            default:
              ConcatCString(message, "unknown value type");
              break;
          }
          ConcatCString(message, "</td>\n");
        }
        /* end table row */
        ConcatCString(message, "</tr>\n");
        HTTP_SendString(request, message);
      }
    }
  }

  HTTP_Send(request, "</table>\n</div>\n");

  /* free the allocated memory for the extra timer names - list */
  TI_deleteExtraTimerTypeField();

  /* free timer structure memory */
  CCTK_TimerDestroyData(tdata);

  HTTP_SetContentFooterString(cctkGH, 0, message);
  retval = HTTP_SendString(request, message);

  String_Delete( message );

  return retval;
}


/* local list implementation - helper functions */

/** init the local extra timer list, has to be called first
 *  creates a list of all available extra timer types
 */
/*@@
  @routine    TI_initExtraTimerTypeField
  @date       Thu Dec 09 18:12:43 2004
  @author     Andre Werthmann
  @desc
              creates a list of all available extra timer types
  @enddesc
  @calls
              TI_isElement, TI_insertElement
  @calledby
              TimerInfoPage
@@*/
static void TI_initExtraTimerTypeField(void)
{
  int i;
  int nt=CCTK_NumTimers();
  char* temp;
  int diff = strlen("in ");

  TI_nettypes=0;

  for (i=0; i<nt; i++)
  {
    temp=strstr(CCTK_TimerName(i), "in ");
    if (temp)
    {
      if (!strstr(temp,"CCTK_"))
      {
        if (!TI_isElement(temp+diff))
        {
          /* do not insert a timer with no schedule bin name */
          if (strlen(temp+diff)>=2)
          {
            if (TI_insertElement(temp+diff))
            {
              /* should never happen, only if system is out of memory... */
              CCTK_WARN (0, "out of memory !");
            }
          }
        }
      }
    }
  }
}

/* checks if a element is member of the list */
/*@@
  @routine    TI_isElement
  @date       Thu Dec 09 18:12:43 2004
  @author     Andre Werthmann
  @desc
              checks if a element is member of the list
  @enddesc
  @calledby
              TI_initExtraTimerTypeField
 @@*/
static int TI_isElement(char* selem)
{
  int i;

  for (i=0; i<TI_nettypes; i++)
  {
    if (!strcmp(selem, TI_ettypes[i]))
    {
      return 1;
    }
  }
  return 0;
}

/* inserts a new element */
/*@@
  @routine    TI_insertElement
  @date       Thu Dec 09 18:12:43 2004
  @author     Andre Werthmann
  @desc
              inserts a new element
  @enddesc
  @calledby
              TI_initExtraTimerTypeField
@@*/
static int TI_insertElement(char* selem)
{
  /* if the next insert would overflow the array, increase the array size by 30 elements */
  if ( TI_nettypes >= TI_maxnettypes )
  {
    TI_maxnettypes+=30;
    TI_ettypes=realloc(TI_ettypes, TI_maxnettypes*sizeof(char*));
  }
  /* out of memory ? */
  if (TI_ettypes==NULL) return 1;

  TI_ettypes[TI_nettypes]=(char*) strdup(selem);
  TI_nettypes++;
  return 0;
}

/* tests whether a timer is the timer type */
/*@@
  @routine    TI_isTimerType
  @date       Thu Dec 09 18:12:43 2004
  @author     Andre Werthmann
  @desc
              tests whether a timer is the timer type
  @enddesc
  @calledby
              TimerInfoPage
@@*/
static int TI_isTimerType(const char* tName, char* tType)
{
  char* temp=strstr(tName, tType);
  if (temp)
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

/* deletes all temporary allocated memory */
/*@@
  @routine    TI_deleteExtraTimerTypeField
  @date       Thu Dec 09 18:12:43 2004
  @author     Andre Werthmann
  @desc
              deletes all temporary allocated memory
  @enddesc
  @calledby
             TimerInfoPage
@@*/
static void TI_deleteExtraTimerTypeField(void)
{
  int i;

  for (i=0; i<TI_nettypes; i++)
  {
    free(TI_ettypes[i]);
  }
}

/* returns the thornname of a given timername */
/*@@
  @routine    TI_thornName
  @date       Thu Dec 09 18:12:43 2004
  @author     Andre Werthmann
  @desc
              returns the thornname of a given timername
  @enddesc
  @calledby
              TimerInfoPage
@@*/
static char* TI_thornName(const char* tn)
{
  int i;
  int sl=strlen(tn)+1;
  char* temp=(char*) malloc(sl);

  /* cut '[xxxx]' and break before the ':' */
  for (i=7; i<sl; i++)
  {
    if (tn[i]==':') break;
    temp[i-7]=tn[i];
  }
  temp[i-7]='\0';

  return temp;
}

/* returns the functionname of a given timername */
/*@@
  @routine    TI_funcName
  @date       Thu Dec 09 18:12:43 2004
  @author     Andre Werthmann
  @desc
              returns the functionname of a given timername
  @enddesc
  @calledby
              TimerInfoPage
@@*/
static char* TI_funcName(const char* tn)
{
  int i;
  int sl=strlen(tn)+1;
  char* ps=(char*) malloc(sl);
  int diff=0;
  char* temp=(char*) malloc(sl);

  /* break after the ':' and return pos after ': ' */
  for (i=0; i<sl; i++)
  {
    if (tn[i]==':') break;
  }
  i+=2;
  if (i>sl-1)
  {
    free (ps);
    free (temp);
    return strdup("malformed timer name");
  }

  strncpy(temp, tn+i, sl-i);
  temp[sl-i]='\0';

  /* cut schedule bin, sometimes doesn't have a leading 'CCTK_' */
  /*if (strstr(temp, " in CCTK_"))*/
  if (strstr(temp, " in "))
  {
    /*diff = strstr(temp, " in CCTK_")-temp;*/
    diff = strstr(temp, " in ")-temp;
    strncpy(ps, temp, diff);
    ps[diff]='\0';
  }
  else
  {
    free (ps);
    free (temp);
    return strdup("no schedule bin");
  }

  free(temp);

  return ps;
}

/* returns schedule bin name without leading 'CCTK_' */
/*@@
  @routine    TI_schedName
  @date       Thu Dec 09 18:12:43 2004
  @author     Andre Werthmann
  @desc
              returns schedule bin name without leading 'CCTK_'
  @enddesc
  @calledby
              TimerInfoPage
@@*/
static char* TI_schedName(const char* tn)
{
  int sl=strlen(tn)+1;
  int diff=strlen("CCTK_");
  char* temp=(char*) malloc(sl);

  /* if the timername does not have a leading 'CCTK_', just copy the name */
  if (strstr(tn, "CCTK_"))
  {
    /* cut 'CCTK_' */
    strncpy(temp, tn+diff, sl-diff);
    temp[sl-diff]='\0';
  }
  else
  {
    strncpy(temp, tn, sl);
  }

  return temp;
}
