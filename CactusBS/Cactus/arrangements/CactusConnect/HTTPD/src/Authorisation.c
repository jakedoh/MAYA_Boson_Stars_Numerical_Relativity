 /*@@
   @file      Authorisation.c
   @date      Fri Sep 15 12:34:59 2000
   @author    Tom Goodale
   @desc 
   Authorisation stuff for webserver
   @enddesc
   @version $Header
 @@*/

#include "cctk.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_CRYPT_H
#include <crypt.h>
#endif

#include "httpd_Map.h"
#include "util_String.h"

#include "httpRequest.h"
#include "Auth.h"
#include "SString_Namespace.h"

#include "base64.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusConnect_HTTPD_Authorisation_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

struct httpUserData
{
  char *password;
  char *encryption_scheme;
};

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

static int AddUser(uMap database, 
                   const char *name,
                   const char *password,
                   const char *encryption_scheme);

static int VerifyPassword(const char *database, 
                          const char *user, 
                          const char *password);

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

static uMap AuthDatabase = NULL;

#define INITIAL_SIZE 32
#define DECODED_SIZE 100

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    HTTP_AuthAddUser
   @date       Fri Sep 15 12:52:09 2000
   @author     Tom Goodale
   @desc 
   Adds a user to a http authentication database.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
int HTTP_AuthAddUser(const char *database, 
                     const char *name,
                     const char *password,
                     const char *encryption_scheme)
{
  int retcode = -1;
  uMap this_database = NULL;

  /* Create the master database if necessary */
  if(!AuthDatabase)
  {
    AuthDatabase = Httpd_MapCreate();
  }

  if(AuthDatabase)
  {
    /* Does this database exist ? */
    this_database = (uMap)Httpd_MapData(AuthDatabase, strlen(database), database);

    if(!this_database)
    {
      this_database = Httpd_MapCreate();

      if(this_database)
      {
        Httpd_MapStore(AuthDatabase, strlen(database), database, (void *)this_database);
      }
      else
      {
        retcode = -2;
      }
    }
  }

  /* Now add the user to the database */
  if(this_database)
  {
    retcode = AddUser(this_database, name, password, encryption_scheme);
  }

  return retcode;
}


 /*@@
   @routine    HTTP_AuthenticateBasic
   @date       Fri Sep 15 13:12:43 2000
   @author     Tom Goodale
   @desc 
   Authenticates an HTTP request against 
   a particular database.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

   @returntype int
   @returndesc
   The authorisation status.
   +1 means that there was no Authorization header.
   0  succesful authentication
   -1 failed authentication
   @endreturndesc

@@*/
int HTTP_AuthenticateBasic(httpRequest *request,
                           const char *database,
                           char *user,
                           int length)
{
  int retval = -1;

  char *auth_string = NULL;
  char *token;
  
  int decoded_size = 0;
  char decoded[DECODED_SIZE+1] = {'\0'};
  
  char *password = NULL;

  int authorised = 0;

  const char *value = HTTP_HeaderValue(request, "Authorization");

  /* Null terminate the user string */ 
  if(user && length > 0)
  {
    *user = 0;
  }

  /* Ok, there's an authentication string here. */
  if(value)
  {
    auth_string = Util_Strdup(value);
    
    token = strtok(auth_string, " ");

    if(token)
    {
      if(CCTK_Equals(token, "Basic"))
      {
        token = strtok(NULL, " ,\t");
        decoded_size = HTTP_b64_pton(token, (unsigned char *) decoded,
                                     DECODED_SIZE);
        
        /* Null terminate string */
        decoded[decoded_size] = 0;

        password = strchr(decoded, ':');

        if(password)
        {
          *password = 0;
          password++;

          authorised = VerifyPassword(database, decoded, password);
          if(user && (int) strlen(user) < decoded_size)
          {
            sprintf(user,"%s", decoded);
          }
        }
      }
    }

    if(auth_string)
    {
      free(auth_string);
    }

    if(authorised)
    {
      retval = 0;
    }
  }
  else
  {
    /* There's no authentication string here */

    retval = 1;
    
  }

  return retval;
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/

 /*@@
   @routine    AddUser
   @date       Fri Sep 15 12:52:37 2000
   @author     Tom Goodale
   @desc 
   Adds a user to a particular database.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
static int AddUser(uMap database, 
                   const char *name,
                   const char *password,
                   const char *encryption_scheme)
{
  int retcode = -1;

  /* Does this user already exist ? */
  struct httpUserData * this_user = (struct httpUserData *)Httpd_MapData(
                                              database, strlen(name), name);

  if(!this_user)
  {
    /* New user */

    this_user = (struct httpUserData *)malloc(sizeof(struct httpUserData));

    if(this_user)
    {
      this_user->password          = Util_Strdup(password);
      this_user->encryption_scheme = Util_Strdup(encryption_scheme);

      retcode = Httpd_MapStore(database, strlen(name), name, (void *)this_user);
    }
  }
  else
  {
    /* Replace user's current data */
    free(this_user->password);
    free(this_user->encryption_scheme);

    this_user->password          = Util_Strdup(password);
    this_user->encryption_scheme = Util_Strdup(encryption_scheme);

    retcode = 0;
  }

  return retcode;
}

 /*@@
   @routine    VerifyPassword
   @date       Fri Sep 15 13:28:50 2000
   @author     Tom Goodale
   @desc 
   Verifies a user and password against a database.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
static int VerifyPassword(const char *database, 
                          const char *user, 
                          const char *password)
{
  int retcode = 0;

  if(AuthDatabase)
  {
    /* Does this database exist ? */
    uMap this_database = (uMap)Httpd_MapData(AuthDatabase,
                                              strlen(database), database);

    if(this_database)
    {
      struct httpUserData *data = (struct httpUserData *) Httpd_MapData(
                                         this_database, strlen(user), user);

      if(data)
      {
        /* Ok, now verify the password. */
        if(CCTK_Equals(data->encryption_scheme, "none"))
        {
          if(!strcmp(data->password, password))
          {
            retcode = 1;
          }
        }
        else if(CCTK_Equals(data->encryption_scheme, "crypt"))
        {
#ifdef HAVE_CRYPT
          if(!strcmp(data->password, crypt(password, data->password)))
          {
            retcode = 1;
          }
#else
	  CCTK_WARN( 1, "Sorry, crypt(3) not supported in this configuration." );
#endif
        }
        else
        {
          CCTK_VWarn(1, __LINE__,__FILE__,CCTK_THORNSTRING,
                     "Unknown encryption algorithm: %s",
                     data->encryption_scheme  );
        } 
      }
    }
  }

  return retcode;
}
