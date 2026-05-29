 /*@@
   @file      Groups.c
   @date      Wed Sep 14 23:47:43 2000
   @author    Gabrielle Allen
   @desc
              Pages about groups.
   @enddesc
   @version   $Id$
 @@*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "cctk.h"
#include "util_String.h"
#include "http_Content.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusConnect_HTTPDExtra_Groups_c)


/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/
int HTTPUTILS_RegisterPages (void);


/********************************************************************
 *********************     Internal Routines   **********************
 ********************************************************************/
static int MessagesPage(const cGH *cctkGH, httpRequest *request, void *data);


 /*@@
   @routine    HTTPUTILS_RegisterPages
   @date       Wed Sep 14 11:29:43 2000
   @author     Gabrielle Allen
   @desc
               Httpd utils registration routine.
   @enddesc
   @calls      HTTP_RegisterPage
               HTTP_ContentLink

   @returntype int
   @returndesc
               0 for success
   @endreturndesc
@@*/
int HTTPUTILS_RegisterPages (void)
{
  /* Register the message board page. */
  HTTP_RegisterPage ("/Messages", MessagesPage, NULL);

  HTTP_ContentLink ("/Messages/index.html", "Message Board",
                    "Collaborative simulation notepad",
                    HTTP_QUICKLINK);

  return (0);
}

/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/

 /*@@
   @routine    MessagesPage
   @date       Sat Sep 16
   @author     Gabrielle Allen
   @desc
               Message board for simulation.
   @enddesc
@@*/
static int MessagesPage (const cGH *GH, httpRequest *request, void *data)
{
  int  retval = 0;
  char currtime[64] = EMPTYSTRING, currdate[64] = EMPTYSTRING;
  String  *message = String_New();
  const char *temp;
  static char *message_board = NULL;


  /* avoid compiler warning about unused parameter */
  data = data;

  if (HTTP_NumArguments( request ) > 0)
  {
    const char *name = HTTP_ArgumentValue (request, "name");
    const char *memo = HTTP_ArgumentValue (request, "memo");

    if (name && *name && memo && *memo)
    {
      size_t message_board_len = 0;
      /* concatenate new message, labeled with current date/time, to the
         message board */
      Util_CurrentTime (sizeof (currtime), currtime);
      Util_CurrentDate (sizeof (currdate), currdate);

      if (message_board)
      {
        message_board_len = strlen (message_board);
        message_board = (char *) realloc (message_board,
                                          message_board_len +
                                          strlen (currtime) + strlen (currdate)+
                                          strlen (name) + strlen (memo) + 60);
      }
      else
      {
        message_board_len = 0;
        message_board = (char *) malloc (strlen (currtime) + strlen (currdate) +
                                         strlen (name) + strlen (memo) + 60);
      }
      if (message_board)
      {
        sprintf (message_board + message_board_len,
                 "<p><strong>%s</strong> %s %s<br />\n\n<em>%s</em></p>",
                 name, currtime, currdate, memo);
      }
    }

    /* Now redirect the browser to the normal message board page */
    if (HTTP_MajorVersion( request ) < 1 ||
        (HTTP_MajorVersion( request ) == 1 && HTTP_MinorVersion( request ) < 1))
    {
      /* Older browsers don't understand 303 */
      temp = "HTTP/1.0 302 Found\r\n"
                     "Location: /Messages/index.html\r\n\r\n";
    }
    else
    {
      temp = "HTTP/1.0 303 See Other\r\n"
                     "Location: /Messages/index.html\r\n\r\n";
    }

    HTTP_Send (request, temp);

    return (0);
  }

  /* Status message */
  HTTP_SendOKHeader( request );

  HTTP_SetDoctype( message );
  HTTP_SendString(request, message);

  /* Start the page */
  HTTP_Send (request, "<html><head><title>CACTUS Messages</title>\n" );

  HTTP_SetHeadInfo( message);
  HTTP_SendString(request, message );

  HTTP_Send(request, "<style type=\"text/css\">\n");
  HTTP_Send(request, "  td.nomsg { background-color: #E9F4D3; }\n");
  HTTP_Send(request, "</style>\n");
  HTTP_Send(request, "</head>\n<body>\n");

  /* Write out the header part */
  HTTP_SetContentHeaderString(GH, 0, message, NULL);
  HTTP_SendString(request, message);

  HTTP_Send (request, 
         "<h1>Message Board</h1>\n"
         "<p>This page can be used to post messages during a \n"
         "simulation. At the moment the messages will disappear \n"
         "when the simulation finishes, but soon there will be an \n"
         "option to save them to a file.</p>\n"
         "<div class=\"centered\">"
         "<form action=\"/Messages/\">\n"
         "<table>\n"
         "<tr><td>Name:</td>"
         "<td><input type=\"text\" size=\"40\" maxlength=\"100\" name=\"name\" "
         "value=\"\" /></td></tr>\n"
         "<tr><td valign=\"top\">Message:</td><td>"
         "<textarea name=\"memo\" rows=\"10\" cols=\"40\"></textarea>\n"
         "</td></tr></table>\n"
         "<p>\n<input type=\"submit\" value=\"Submit Message\" />\n</p>\n"
         "<table width=\"80%\"><tr><td>\n"
         "<h2>Messages:</h2></td></tr>\n"
         "<tr><td><table width=\"100%\" cellpadding=\"5\" cellspacing=\"5\">"
         "<tr><td class=\"nomsg\">" );

  temp = message_board ? message_board :
                 "No messages yet ... use the form above to add one";
  HTTP_Send (request, temp );

  HTTP_Send (request, "</td></tr></table>\n" );
  HTTP_Send (request, "</td></tr></table>\n</form>\n</div>\n" );

  /* Write out the footer part. */
  HTTP_SetContentFooterString(GH, 0, message);
  retval = HTTP_SendString(request, message);

  String_Delete( message );
  return retval;
}
