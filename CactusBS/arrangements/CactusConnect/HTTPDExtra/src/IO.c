 /*@@
   @file      IO.c
   @date      Sun Sep 17 16:19:17 2000
   @author    Tom Goodale
   @desc
              Stuff for displaying IO callbacks on the web page.
   @enddesc
   @version   $Id$
 @@*/

#include <stdlib.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include <stdlib.h>
#include <fcntl.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include "CactusBase/IOUtil/src/ioutil_AdvertisedFiles.h"

#include "http_Content.h"

#include "HTTPD_FileList.h"

/* SW Temporary, while testing the SString module */
#include "CactusConnect/HTTPD/src/SString_Namespace.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusConnect_HTTPDExtra_IO_c);

/********************************************************************
 ********************    Macro Definitions   ************************
 ********************************************************************/
#ifndef O_BINARY
#define O_BINARY 0
#endif

/********************************************************************
 *********************     Internal Typedefs   **********************
 ********************************************************************/

/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
static int IOFileListener(const cGH *GH, const char *filename,
                          const ioAdvertisedFileDesc *description);
static int AdvertisedFilePage(const cGH *GH, httpRequest *request, void *data);
static int ViewportFilePage(const cGH *GH, httpRequest *request, void *data);
static int SendFilePage(const cGH *GH, httpRequest *request, void *data);


/********************************************************************
 ********************    External Routines   ************************
 ********************************************************************/
int HTTP_util_RegisterIOPages(void);


/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/
static FileList *filelist   = NULL;

 /*@@
   @routine    HTTP_util_RegisterIOPages
   @date       Wed Sep 14 11:29:43 2000
   @author     Gabrielle Allen
   @desc
               Registers HTTPDExtra's Output and ViewPort page.
   @enddesc
   @calls      IOUtil_RegisterAdvertisedFileListener
               HTTP_RegisterPage
               HTTP_ContentLink

   @returntype int
   @returndesc
               0 for success
   @endreturndesc
@@*/
int HTTP_util_RegisterIOPages(void)
{
  ioAdvertisedFileListenerCallbacks listener;

  listener.advertise = IOFileListener;

  IOUtil_RegisterAdvertisedFileListener ( NULL, CCTK_THORNSTRING,
                                         &listener);

  HTTP_RegisterPage( "/Output/index.html", AdvertisedFilePage, NULL);
  HTTP_RegisterPage( "/Output/viewport.html", ViewportFilePage, NULL);

  HTTP_ContentLink( "/Output/index.html", "Files",
                   "Downloadable files",
                   HTTP_QUICKLINK);
  HTTP_ContentLink( "/Output/viewport.html", "Viewport",
                   "Viewport for certain output files",
                   HTTP_QUICKLINK);

  HTTP_RegisterPage( "/Output", SendFilePage, NULL);

  return 0;
}


/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
 /*@@
   @routine    HTTP_util_ExistingFileEntry
   @date       Wed Oct 07 12:49:02 2009
   @author     Frank Loeffler
   @desc
               Returns the pointer to an entry with the given
               filename set or NULL if nothing was found.
   @enddesc
@@*/
static httpFileItem * HTTP_util_ExistingFileEntry(FileList *filelist,
                                                  const char *filename)
{
  if ( !filelist )
    return NULL;

  /* Look for an entry with this filename */
  const size_t n = NumberOfFiles( filelist );
  size_t i;
  httpFileItem *item;
  for (i = 0; i < n; i++)
  {
    item = HTTPD_FileList_Item(filelist, i);
    /* When entry found, return its pointer */
    if (!StringCompareCString(item->filename, filename))
    {
      return item;
    }
  }
  return NULL;
}

/*@@
   @routine    IOFileListener
   @date       Sun Sep 17 17:56:22 2000
   @author     Tom Goodale
   @desc
               Listener for advertised files.
   @enddesc
@@*/
static int IOFileListener(const cGH *GH, const char *filename,
                          const ioAdvertisedFileDesc *description)
{
  size_t position = 0;
  (void) (GH + 0); /* avoid compiler warning about unused parameter */

  if( !filelist )
    filelist = HTTPD_FileList_New();
  httpFileItem *entry = HTTP_util_ExistingFileEntry(filelist, filename);
  /* If there is already an entry, we cannot change it directly (const).
     We _could_ remove it and insert a new entry, but that would probably
     be too expensive. So we do nothing and assume nothing changed
     anyway in the entry */
  if (entry)
  {
    return 0;
  }
  entry = (httpFileItem *)malloc(sizeof( httpFileItem));
  if (entry)
  {
    entry->filename    = String_Make( filename);
    entry->thorn       = String_Make( description->thorn);
    entry->varname     = String_Make( description->varname);
    entry->mimetype    = String_Make( description->mimetype);
    entry->slice       = String_Make( description->slice);
    entry->description = String_Make( description->description);
    String *linknameString = String_Make( filename);
    entry->linkname    = linknameString;
    /* Need to mangle the filename to get a decent linkname */
    while( FindCharFrom( linknameString, '/', &position ) )
      SetNthChar( linknameString, position, '@');
    AppendFile( filelist, entry );
  }
  return 0;
}

static void
IndentAndConcatCString( String *message, const char *c, int depth )
{
  int i;
  for ( i = 0; i < depth; i++ )
    ConcatCString( message, "\t" );
  ConcatCString( message, c );
}
/*
 * SW
 * Build a HTML representation of the contents of filelist.
 * Denis wanted variable name on the topmost level.
 * This is about third iteration. Here, use tables to get nesting effect.
 */
static size_t
buildList( const FileList *list, String * message, size_t itemNo, int depth )
{
  const size_t n = NumberOfFiles( filelist );
  size_t i = 0;
  httpFileItem *lastItem = FileItem( filelist, itemNo );

  for ( i = itemNo; i < n; i++ )
  {
    httpFileItem *file = FileItem( filelist, i );
    if( depth == 0 )
    {
      IndentAndConcatCString( message, "<tr>", depth );
      ConcatCString( message, "<td colspan=\"5\" class=\"variable\">" );
      Concat( message, file->varname );
      ConcatCString( message, "&nbsp;" );
      IndentAndConcatCString( message, "</td></tr>\n", depth );
      i = buildList( list, message, i, 1 );
    }
    if( depth == 1 )
    {
      if( Equals( file->varname, lastItem->varname ) )
      {
        IndentAndConcatCString( message, "<tr>", depth );
        ConcatCString( message, "<td class=\"spacer\">&nbsp;</td>"
                        "<td class=\"thorn\" colspan=\"4\">" );
        Concat( message, file->thorn );
        IndentAndConcatCString( message, "</td></tr>\n", depth );
        i = buildList( list, message, i, 2 );
      }
      else
        break;
    }
    if( depth == 2 )
    {
      if( Equals( file->thorn, lastItem->thorn )
      &&  Equals( file->varname, lastItem->varname ) )
      {
        IndentAndConcatCString( message, "<tr>", depth );
        ConcatCString( message, "<td class=\"spacer\">&nbsp;</td>"
                      "<td class=\"spacer\">&nbsp;</td>\n"
                      "<td class=\"slice\">" );
        Concat( message, file->slice );

        ConcatCString( message, "\n" );
        IndentAndConcatCString( message, "</td><td class=\"link\">\n", depth );
        IndentAndConcatCString( message, "<a href=\"/Output/", depth );
        Concat( message, file->linkname );
        ConcatCString( message, "\">" );
        Concat( message, file->filename );
        ConcatCString( message, "</a>\n" );
        IndentAndConcatCString( message, "</td><td class=\"desc\">\n", depth );
        IndentAndConcatCString( message, "", depth );
        Concat( message, file->description );
        ConcatCString( message, "\n" );
        IndentAndConcatCString( message, "</td></tr>\n", depth );
      }
      else
        break;
    }
  }
  ConcatCString( message, "\n" );
  return i - 1;
}

 /*@@
   @routine    AdvertisedFilePage
   @date       Sun Sep 17 18:26:01 2000
   @author     Tom Goodale
   @desc
               Page to deal with advertised files.
   @enddesc
@@*/
static int AdvertisedFilePage(const cGH *GH, httpRequest *request, void *data)
{
  int retval = 0;
  String *message = String_New();

  data = data; /* avoid compiler warning about unused parameter */

  HTTP_SendOKHeader( request );

  HTTP_SetDoctype( message );
  HTTP_SendString( request, message);

  HTTP_Send( request, "<html><head>\n" );
  HTTP_Send( request, "<title>Cactus Downloadable Files</title>\n");

  HTTP_SetHeadInfo( message);
  HTTP_Send( request, "<style type=\"text/css\">"
                  "\t.files td { text-align: left; font-family: sans-serif;\n"
                  "\t              padding-left: 1ex; padding-right: 1ex; }\n"
                  "\t.files td.variable { background-color: #ff9;\n"
                  "\t              font-weight: bold; }\n"
                  "\t.files td.thorn { background-color: #dff;\n"
                  "\t              color: black; }\n"
                  "\t.files td.slice { background-color: #fde;\n"
                  "\t              font-family: monospace; }\n"
                  "\t.files td.link { \n"
                  "                      }\n"
                  "\t.files td.desc { background-color: #dff;\n"
                  "\t              font-style: italic; }\n"
                  "\t#filescaption { border-bottom: thin black solid;\n"
                  "\t               margin-bottom: 1em; }\n"
                  "\t</style>\n");
  HTTP_SendString( request, message );
 
  HTTP_Send( request, "</head>\n<body>\n");

  HTTP_SetContentHeaderString( GH, 0, message, NULL);
  HTTP_SendString( request, message);

  HTTP_Send( request, "<h1>Downloadable Files</h1>\n");

  HTTP_Send( request,
         "<p>From this page you can download various output files \n"
         "from the simulation. Depending on the software available on your\n"
         " local machine, you can change the browser properties to launch\n"
         " files directly to visualization clients.</p> \n "
         "<p>Some of these output files can be directly viewed on \n"
         "a browser on the simulation \n"
         "<a href=\"/Output/viewport.html\">viewport</a></p>\n"
         "<p>Many IO methods have <dfn>steerable</dfn> parameters which \n"
         "allow you to e.g. add fields and customise behaviour.\n"
         "Depending on your authorisation, you can access the \n"
         "<a href=\"/Parameters/index.html\">parameter steering page</a>"
         "</p>\n");
  HTTP_Send( request,
         "<div class=\"centered\">\n"
         "<table class=\"files\" id=\"filescaption\" "
         "frame=\"below\" cellpadding=\"0\" cellspacing=\"0\">\n"
         "<tr><td class=\"variable\">\n"
         "Variable</td><td class=\"thorn\">Thorn</td>\n"
         "<td class=\"slice\"> Slice</td><td class=\"file\"> File</td>\n"
         "<td class=\"desc\"> Description</td></tr>\n</table>\n"
         "\n"
         );

  Truncate( message, 0 );
  SortFilesAccordingTo( filelist, Variable_Thorn_Slice );

  HTTP_Send( request, "<table cellpadding=\"0\" cellspacing=\"0\""
                  " class=\"files\">\n");
  buildList( filelist, message, 0, 0 );

  HTTP_SendString( request, message);
  HTTP_Send( request, "\n</table>");
  HTTP_Send( request, "\n</div>");

  HTTP_SetContentFooterString( GH, 0, message);
  retval = HTTP_SendString( request, message);

  String_Delete( message );
  return retval;
}

 /*@@
   @routine    SendFilePage
   @date       Sun Sep 17 18:43:40 2000
   @author     Tom Goodale
   @desc
               Sends an advertised file.
   @enddesc
@@*/

static int SendFilePage( const cGH *GH, httpRequest *request, void *data)
{
  int found = 0;

  const size_t n = NumberOfFiles( filelist );
  size_t i;
  /* avoid compiler warning about unused parameters */
  (void) (GH + 0);
  data = data;

  for( i = 0; i < n; i++ )
  {
    httpFileItem *file = FileItem( filelist, i );

    if( !CompareCString( file->linkname, HTTP_Residual( request ) ))
    {
      int filedes = open( GetBuffer( file->filename), O_RDONLY | O_BINARY);
      if(filedes >= 0)
      {
        HTTP_Send( request, "HTTP/1.0 200 OK\r\n");

        HTTP_Send( request, "Content-Type: " );
        HTTP_SendString( request, file->mimetype);
        HTTP_Send( request, "\r\n\r\n" );

        HTTP_ContentSendFromFile( request, filedes);

        close( filedes);
        found = 1;
      }
      else
      {
        HTTP_Send( request,"HTTP/1.0 500 Server Internal Error\r\n");

        HTTP_Send( request, "Content-Type: text/html\r\n\r\n");
        HTTP_Send( request, "<html>\n<head>\n");
        HTTP_Send( request, "<title>Error 500: Internal Error</title>\n");
        HTTP_Send( request, "</head>\n<body>\n");
        HTTP_Send( request, "<div class=\"centered\"><p>Unable to open " );
        HTTP_SendString( request, file->filename );
        HTTP_Send( request, "</p></div>\n");
        HTTP_Send( request, "</body>\n</html>\n" );
      }
      break;
    }
  }

  if(!found)
  {
    HTTP_Send( request,"HTTP/1.0 404 Not Found\r\n");

    HTTP_Send( request,"Content-Type: text/html\r\n\r\n");

    HTTP_Send( request, "<html>\n<head>\n" );
    HTTP_Send( request, "<title>Error 404: Not Found</title>\n" );
    HTTP_Send( request, "</head>\n<body>\n" );
    HTTP_Send( request, "<div class=\"centered\"><p>" );
    HTTP_Send( request, HTTP_URI( request ) );
    HTTP_Send( request, " does not exist</p></div>\n" );
    HTTP_Send( request, "</body>\n</html>\n" );
  }
  return 0;
}

 /*@@
   @routine    AdvertisedFilePage
   @date       Sun Sep 17 18:26:01 2000
   @author     Tom Goodale
   @desc
               Page to deal with advertised files.
   @enddesc
@@*/

static int ViewportFilePage(const cGH *GH, httpRequest *request, void *data)
{
  DECLARE_CCTK_PARAMETERS

  int retval = 0;
  int foundone = 0;
  String *message = String_New();
  const size_t n = NumberOfFiles( filelist );
  size_t i;

  (void) (GH + 0);
  data = data; /* avoid compiler warning about unused parameters */

  HTTP_SendOKRefreshHeader( request, viewport_refresh_seconds );
  
  HTTP_SetDoctype( message );
  HTTP_SendString( request, message);
  
  HTTP_Send( request, "<html>\n<head>\n");
  HTTP_Send( request, "<title>Cactus Downloadable Files</title>\n");

  HTTP_SetHeadInfo( message);

  HTTP_Send( request, "<style type=\"text/css\">"
                      ".files td { text-align: center; vertical-align: middle;"
                      " font-size: smaller; font-family: sans-serif; }\n"
                      ".files .spacer { width: 10ex;\n"
                      "               }\n"
                      ".files .variable { background-color: #ff9;\n"
                      "              font-weight: bold; }\n"
                      ".files .thorn { background-color: #dff;\n"
                      "              color: black; }\n"
                      ".files .slice { background-color: #fde;\n"
                      "              font-family: monospace; }\n"
                      ".files .link { }\n"
                      ".files .desc { font-style: italic;\n"
                      "               background-color: #dff; }\n"
                      ".files img { border: 0;\n"
                      "               }\n"
                      "</style>\n");

  HTTP_SendString( request, message );

  HTTP_Send( request, "</head>\n<body>\n");

  HTTP_SetContentHeaderString( GH, 0, message, NULL);
  HTTP_SendString( request, message);

  HTTP_Send( request, "<h1>Viewport</h1>\n");

  HTTP_Send( request, 
         "<p>This page displays certain types of the output files \n"
         "from the <a href=\"/Output/index.html\">download</a> page \n"
         "as images (currently only JPEGs [mime type image/jpeg]).</p>\n"
         "<p>Many IO methods have <dfn>steerable</dfn> parameters which \n"
         "allow you to e.g. add fields and customise behaviour.\n"
         "Depending on your authorisation, you can access the \n"
         "<a href=\"/Parameters/index.html\">parameter steering page</a></p>\n"
         "<div class=\"centered\">\n");

  SortFilesAccordingTo( filelist, Variable_Thorn_Slice );

  for ( i = 0; i < n; i++ )
  {
    httpFileItem *file = FileItem( filelist, i );
    if ( !CompareCString( file->mimetype, "image/jpeg"))
    {
      if (!foundone)
      {
        HTTP_Send( request, 
               "<table class=\"files\" cellspacing=\"5\" cellpadding=\"5\" "
               "rules=\"groups\" >\n"
               "<thead>\n<tr><th>Variable<br />Slice</th>"
               "<th>Description</th><th>Image</th></tr>\n</thead>\n");
        foundone = 1;
      }
      HTTP_Send( request, "<tr>\n" );
      HTTP_Send( request, "<td><span class=\"variable\">" );
      HTTP_SendString( request, file->varname );
      HTTP_Send( request, "</span><br />\n" );
      HTTP_Send( request, "<span class=\"slice\">" );
      HTTP_SendString( request, file->slice );
      HTTP_Send( request, "</span>" );
      HTTP_Send( request, "</td>\n" );
      HTTP_Send( request, "<td><span class=\"desc\">" );
      HTTP_SendString( request, file->description );
      HTTP_Send( request, "</span></td>\n" );
      HTTP_Send( request, "<td class=\"linkk\"><a href=\"/Output/" );
      HTTP_SendString( request, file->linkname );
      HTTP_Send( request, "\">" );
      HTTP_Send( request, "<img width=\"100\" height=\"100\" src=\"" );
      HTTP_SendString( request, file->linkname );
      HTTP_Send( request, "\" alt=\"" );
      HTTP_SendString( request, file->linkname );
      HTTP_Send( request, "\"/>" );
      HTTP_Send( request, "</a></td>\n" );
      HTTP_Send( request, "</tr>\n" );
    }
  }

  if( foundone)
    HTTP_Send( request, "</table>\n");
  else
  {
    HTTP_Send( request, "<strong>\n<p>No viewable images registered!</p>\n"
           "</strong>\n");
  }

  HTTP_Send( request, "</div>\n");

  HTTP_SetContentFooterString( GH, 0, message);
  retval = HTTP_SendString( request, message);

  String_Delete( message );

  return retval;
}

