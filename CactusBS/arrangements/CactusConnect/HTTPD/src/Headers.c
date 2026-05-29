 /*@@
   @file      Headers.c
   @date      Wed Sep 17 23:47:43 2000
   @author    Gabrielle Allen
   @desc 
   Functions to return standard headers and footers for HTML pages
   @enddesc
   @version $Header$
 @@*/

#include "cctk.h"

#include "util_String.h"

#include "SString_Namespace.h"
#define EMPTYSTRING { '\0' }

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusConnect_HTTPD_Headers_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

extern struct httpLink *ContentLinks;

struct httpLink
{
  struct httpLink *next;
  char *URL;
  char *name;
  char *description;
  int flags;
};

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/


/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

static const char *cactus_styles =
"<style type=\"text/css\">\n"
"       body { color: black; background-color: white; }\n"
"       h1, h2 { text-align: center; }\n"
"       dfn { font-style: italic; }\n"
"       td { vertical-align: top; }\n"
"       a:link { color: #1B831D; }\n"
"       a:visited { color: #768000; }\n"
"       a:active { color: green; }\n"
"       div.centered { text-align: center; } \n"
"       div.centered table { margin: auto; } \n"
"       span.hilite { color: red; } \n"
"       td.menu { color: black; background-color: #E5FFA2; \n"
"               text-align: left; vertical-align: top; width: 20ex; \n"
"               font-size: small; }\n"
"       td.menu h2 { font-weight: normal; font-size: medium;  \n"
"                    text-align: left; margin-top: 0; } \n"
"       td.menu h3 { font-weight: bold; font-size: small;  \n"
"                    margin-top: 1.2em; margin-bottom: 0; } \n"
"       td.menu span.simulation_name { font-style: italic; } \n"
"       td.menu kbd { font-family: monospace; font-style: normal; } \n"
"       .footer td { font-size: small; vertical-align: top; } \n"
"       .footer td.by { text-align: right; } \n"
"       .footer img { border: 0; } \n"
"       div.banner { text-align: center; } \n"
"       div.banner table { margin: auto; } \n"
"       div.banner img { border: 0; } \n"
"</style>\n";

static const char *cactus_footer =
"\n</td></tr>\n"
"\n<tr><td colspan=\"2\">\n"
"\n<table class=\"footer\" width=\"100%\" border=\"0\">\n"
"<tr><td colspan=\"2\"><hr /></td></tr>\n"
"<tr><td>\n"
"<a href=\"http://www.cactuscode.org\">"
"<img src=\"/Images/www.gif\" alt=\"www.CactusCode.org\" /></a>\n"
"</td><td class=\"by\">"
"Cactus Web Interface by \n"
"<a href=\"mailto:cactusmaint@cactuscode.org\">The Cactus Team</a><br />\n"
"<a href=\"/About.html\"><em>About this Server</em></a>\n"
"</td></tr></table>\n"
"</td></tr></table>\n"
"</body></html>\n";

static const char * cactus_doctype = 
"<?xml version=\"1.0\" encoding=\"ISO-8859-15\"?>\n"
"<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Strict//EN\"\n"
"        \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd\">\n";

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

void HTTP_SetHeadInfo( String *header)
{
  SetToCString( header, cactus_styles );
  ConcatCString( header, "<meta name=\"content-type\" "
                         "content=\"text/html; charset=iso-8859-15\" />\n");
}
void HTTP_SetDoctype( String *header)
{
  SetToCString( header, cactus_doctype );
}
 /*@@
   @routine    HTTP_ContentHeader
   @date       Sat Sep 16 15:22:59 2000
   @author     Gabrielle Allen
   @desc 
   Returns header for HTML pages
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
extern void CCTK_GetRunTitleString( String *s );

int HTTP_SetContentHeaderString(const cGH *GH, int choice, 
                          String *header, const String *menu)
{
  String *title = String_New();
  String *quicklinks = String_New();
  char parfile[200] = EMPTYSTRING;
  char currentdate[50] = EMPTYSTRING;
  char currenttime[50] = EMPTYSTRING;
  struct httpLink *link;
  char *file;
  char *dir;

  Truncate( header,0);

  if (choice == 0)
  {
    ConcatCString( header, "<table cellpadding=\"10\" width=\"100%\">\n"
                           "<tr>\n"
                           "<td class=\"menu\" valign=\"top\">\n"
                           "<h2><a href=\"/\">Master Run Page</a></h2>\n");
  }

  if(ContentLinks)
  {
    SetToCString( quicklinks,"<h3>Options:</h3>\n");
    for(link = ContentLinks; link; link=link->next)
    {
      ConcatCString(quicklinks, "<a href=\"");
      ConcatCString(quicklinks, link->URL);
      ConcatCString(quicklinks, "\">");
      ConcatCString(quicklinks, link->name);
      ConcatCString(quicklinks, "</a><br />\n");
    }
  }
 
  if (choice == 0)
  {
    /* Find strings needed for nonmain-page headers */
    CCTK_GetRunTitleString(title); 
    CCTK_ParameterFilename(sizeof(parfile),parfile);
    Util_SplitFilename(&dir,&file,parfile);

    Util_CurrentDate(sizeof(currentdate),currentdate);
    Util_CurrentTime(sizeof(currenttime),currenttime);

    /* Build the header */
    ConcatCString( header, "\n"
                           "<h3>Environment:</h3>\n"
                           "Time: ");
    ConcatCString( header, currenttime );
    ConcatCString( header, "<br />\n"
                           "Date: " );
    ConcatCString( header, currentdate );
    ConcatCString( header, "<br />\n" );
    ConcatCString( header, "<h3>Simulation:</h3>\n"
                           "<span class=\"simulation_name\">");
    Concat( header, title );
    ConcatCString( header, "</span><br />\n"
                           "<kbd>");
    ConcatCString( header, file );
    ConcatCString( header, "</kbd><br />\n"
                           "Iteration: ");
    ConcatDecimal( header, GH ? GH->cctk_iteration : 0 );
    ConcatCString( header, "<br />\nPhysical time: ");
    ConcatFormattedDouble( header, GH ? GH->cctk_time : 0, 4, 2, SFMT_DEFAULT );
    ConcatCString( header, "<br />\n");
    if (ContentLinks)
    {
      Concat( header, quicklinks);
    }

    if (menu)
    {
      Concat( header, menu);
    }

    ConcatCString( header, 
            "</td>\n<td>\n\n");
  }

  String_Delete( quicklinks );
  String_Delete( title );
  return Length(header);

}

 /*@@
   @routine    HTTP_ContentFooter
   @date       Sat Sep 16 15:22:59 2000
   @author     Tom Goodale
   @desc 
   Returns footer for HTML pages
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
int HTTP_SetContentFooterString(const cGH *GH, int choice, String *footer)
{
  /* avoid compiler warnings about unused parameters */
  GH = GH;
  choice = choice;

  Truncate( footer, 0 );
  SetToCString(footer,cactus_footer);
  return Length(footer);
}

int HTTP_ContentHeader(const cGH *GH, int choice, int len,
                          char *header, const char *menu)
{
  String *menuString = NULL;
  String *headerString = String_New();
  int retval = 0;
  
  if( menu )
    menuString = String_Make( menu );

  retval = HTTP_SetContentHeaderString( GH, choice, headerString, menuString );

  /* SW: original didn't use len, although it was clearly meant
   * to indicate the length of the buffer, and all calls I saw to it
   * passed a buffer length */
  if( len > 0 )
  {
    SetToBuffer( headerString, header, len );
  }
  String_Delete( menuString );
  String_Delete( headerString );
  return retval;
}

int HTTP_ContentFooter(const cGH *GH, int choice, int len, char *header)
{
  String *headerString = String_New();
  int retval = 0;
  
  retval = HTTP_SetContentFooterString( GH, choice, headerString );

  /* note: original didn't use len, although it was clearly meant
   * to indicate the length of the buffer, and all calls I saw to it
   * passed a buffer length */
  if( len > 0 )
  {
    SetToBuffer( headerString, header, len );
  }
  String_Delete( headerString );
  return retval;
}
