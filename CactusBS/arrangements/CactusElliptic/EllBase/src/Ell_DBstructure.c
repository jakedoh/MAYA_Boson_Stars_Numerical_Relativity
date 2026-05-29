/*@@
   @file      Ell_DBstructure.c
   @date      
   @author    Gerd Lanfermann
   @desc
   Database for elliptic parameters
   @enddesc
   @version   $Id$
 @@*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_WarnLevel.h"
#include "cctk_FortranString.h"

#include "StoreNamedData.h"

#include "Ell_DBstructure.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusElliptic_EllBase_Ell_DBstructure_c)

/********************************************************************
 ********************    External Routines   ************************
 ********************************************************************/
void CCTK_FCALL CCTK_FNAME(Ell_IsKey)
     (int *ierr, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Ell_DeleteKey)
     (int *ierr, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Ell_SetRealKey)
     (int *ierr, CCTK_REAL *value, ONE_FORTSTRING_ARG); 
void CCTK_FCALL CCTK_FNAME(Ell_SetIntKey)
     (int *ierr, CCTK_INT *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Ell_SetStrKey)
     (int *ierr, TWO_FORTSTRINGS_ARGS);
void CCTK_FCALL CCTK_FNAME(Ell_GetRealKey)
     (int *ierr, CCTK_REAL *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Ell_GetIntKey)
     (int *ierr, CCTK_INT *value, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(Ell_GetStrKey)
     (int *nchar, char **cstring,ONE_FORTSTRING_ARG);

/********************************************************************
 ********************    Internal Typedefs   ************************
 ********************************************************************/

struct t_ellthingy
{
  int type;
  int been_set;
  union 
  {
    CCTK_REAL r;
    CCTK_INT  i;
    char *s;
  } vals;
};

/********************************************************************
 ********************    Static Variables   *************************
 ********************************************************************/

static pNamedData *EllInfoDB;

/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/

int Ell_CreateKey(int vartype, const char *keychain) 
{

  DECLARE_CCTK_PARAMETERS

  struct t_ellthingy*  new;
  int retval;

  if ((struct t_ellthingy*)GetNamedData(EllInfoDB, keychain)) 
  {
    retval = ELLCREATE_TWICE;
  } 
  else 
  {

    new = (struct t_ellthingy*)malloc(sizeof(struct t_ellthingy));
    
    new->type     = vartype;
    new->been_set = 0;
    new->vals.r   = 0.0;
    new->vals.i   = 0;
    new->vals.s   = NULL;
    
    if (StoreNamedData(&EllInfoDB, keychain, new)!=0) 
    {
      if (!CCTK_Equals(elliptic_verbose,"no")) 
      {
        char *msg;
        const char *name=CCTK_VarTypeName(vartype);
        msg = (char *)malloc( 
             (200 + strlen(keychain) + strlen(name) )*sizeof(char) );
        sprintf(msg,"Failed to create %s (%s)",keychain,name);
        CCTK_INFO(msg);
        free(msg);
      }
      retval = ELLCREATE_FAILED;
    }
    else
    {
      if (!CCTK_Equals(elliptic_verbose,"no")) 
      {
        char *msg;
        const char *name=CCTK_VarTypeName(vartype);
        msg = (char *)malloc( 
             (200 + strlen(keychain) + strlen(name) )*sizeof(char) );
        sprintf(msg,"Created %s (%s)",keychain,name);
        CCTK_INFO(msg);
        free(msg);
      }
      retval = 0;
    }

  }
    
  return(retval);

  

} 

int Ell_IsKey(const char *keychain) 
{
  int retval=ELL_ISNOKEY;
  if ((struct t_ellthingy*)GetNamedData(EllInfoDB, keychain)) 
  {
    retval = 0;
  } 
  else
  {
    retval = ELL_ISNOKEY;
  }
  return(retval);
}

void CCTK_FCALL CCTK_FNAME(Ell_IsKey)(int *ierr, ONE_FORTSTRING_ARG) 
{
  ONE_FORTSTRING_CREATE(key)
  *ierr = Ell_IsKey(key);
  free(key);
}

int Ell_UnsetKey(const char *keychain) 
{

  struct t_ellthingy*  getme;
  int retval;

  getme = (struct t_ellthingy*)GetNamedData(EllInfoDB, keychain);
  
  if (!(getme)) 
  {
    CCTK_VWarn(4,__LINE__,__FILE__,CCTK_THORNSTRING,
               "Ell_GetRealKey: Cannot get structure with key %s",keychain);
    retval = ELLGET_NOKEY;
  } 
  else 
  {
    getme->been_set = ELL_NO;
    retval = 0;
  }
  return(retval);
}

int Ell_DeleteKey(const char *keychain) 
{
  int retval;
  /* avoid compiler warnings */
  keychain = keychain;
  CCTK_INFO("Ell_DeleteKey: Routine not implemented yet!");
  retval = 1;
  return(retval);
}

void CCTK_FCALL CCTK_FNAME(Ell_DeleteKey)(int *ierr, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(key);
  *ierr = Ell_DeleteKey(key);
  free(key);
}
   

int Ell_SetRealKey(CCTK_REAL value, const char *keychain) 
{
  int retval;
  struct t_ellthingy *setme;

  setme = GetNamedData(EllInfoDB, keychain);
  if (!(setme)) 
  { 
    retval = ELLSET_FAILED;
  }  
  else if (setme->type!=CCTK_VARIABLE_REAL) 
  {
    retval = ELLSET_BADTYPE;
    CCTK_VWarn(4,__LINE__,__FILE__,CCTK_THORNSTRING,
               "Ell_SetRealKey: Key %s not of type CCTK_REAL (type %d)",
               keychain,setme->type);
  }
  else 
  {
    setme->type     = CCTK_VARIABLE_REAL;
    setme->vals.r   = value;
    setme->been_set = ELL_YES;
    retval=0;
  }

  return(retval);
}

void CCTK_FCALL CCTK_FNAME(Ell_SetRealKey)
     (int *ierr, CCTK_REAL *value, ONE_FORTSTRING_ARG) 
{
  ONE_FORTSTRING_CREATE(key)
  *ierr = Ell_SetRealKey(*value, key);
  free(key);
} 

int Ell_SetIntKey(CCTK_INT value, const char *keychain) 
{
  int retval;
  struct t_ellthingy *setme;

  setme = GetNamedData(EllInfoDB, keychain);
  if (!(setme)) 
  { 
    retval = ELLSET_FAILED;
  }
  else if (setme->type!=CCTK_VARIABLE_INT) 
  {
    retval = ELLSET_BADTYPE;
    printf("Ell_SetIntKey: The key you try to set is not of type CCTK_INT: >%s< (type %d)\n",
           keychain,setme->type);
  }
  else {
    setme->type     = CCTK_VARIABLE_INT;
    setme->vals.i   = value;
    setme->been_set = ELL_YES;
    retval=0;
  }
  return(retval);
}

void CCTK_FCALL CCTK_FNAME(Ell_SetIntKey)
     (int *ierr, CCTK_INT *value, ONE_FORTSTRING_ARG) 
{
  ONE_FORTSTRING_CREATE(key)
  *ierr = Ell_SetIntKey(*value, key);
  free(key);
} 

int Ell_SetStrKey(char *value, const char *keychain) 
{
  int retval;
  struct t_ellthingy *setme;

  setme = GetNamedData(EllInfoDB, keychain);
  if (!(setme)) 
  { 
    retval = ELLSET_FAILED;
  } 
  else if (setme->type!=CCTK_VARIABLE_STRING) 
  {
    retval = ELLSET_BADTYPE;
    CCTK_VWarn(4,__LINE__,__FILE__,CCTK_THORNSTRING,
               "Ell_SetStrKey: Key %s not of type STRING (type %d)",
               keychain,setme->type);
  }
  else 
  {
    setme->type     = CCTK_VARIABLE_STRING;
    setme->vals.s   = strdup(value);
    setme->been_set = ELL_YES;
    retval = 0;
  }

  return(retval);

}

void CCTK_FCALL CCTK_FNAME(Ell_SetStrKey)
     (int *ierr, TWO_FORTSTRINGS_ARGS) 
{
  TWO_FORTSTRINGS_CREATE(value,key)
  *ierr = Ell_SetStrKey(value, key);
  free(value);
  free(key);
} 


int Ell_GetRealKey(CCTK_REAL *value, const char *keychain) 
{

  struct t_ellthingy *getme=NULL;
  int retval;

  getme = (struct t_ellthingy*)GetNamedData(EllInfoDB, keychain);

  if (!(getme)) 
  {
    CCTK_VWarn(4,__LINE__,__FILE__,CCTK_THORNSTRING,
               "Ell_GetRealKey: Cannot get structure for key %s",keychain);
    retval = ELLGET_NOKEY;
  }
  else if (getme->type!=CCTK_VARIABLE_REAL) 
  {
    CCTK_VWarn(4,__LINE__,__FILE__,CCTK_THORNSTRING,
               "Ell_GetRealKey: Not getting a CCTK_REAL value off key %s "
               "type %d (need %d)\n",
               keychain,getme->type,CCTK_VARIABLE_REAL);
    retval = ELLGET_BADTYPE;
  }
  else if (getme->been_set==ELL_NO) 
  {
    CCTK_VWarn(4,__LINE__,__FILE__,CCTK_THORNSTRING,
               "Ell_GetRealKey: Key %s has not been set to any value",
               keychain);
    retval = ELLGET_NOTSET;
  }
  else 
  {
    *value=getme->vals.r;
    retval=0;
  }
  return(retval);
}

void CCTK_FCALL CCTK_FNAME(Ell_GetRealKey)
     (int *ierr, CCTK_REAL *value, ONE_FORTSTRING_ARG) 
{
  ONE_FORTSTRING_CREATE(key)
  *ierr = Ell_GetRealKey(value, key); 
  free(key);
}

int Ell_GetIntKey(CCTK_INT *value,const char *keychain) 
{

  struct t_ellthingy *getme=NULL;
  int retval;

  getme = (struct t_ellthingy*)GetNamedData(EllInfoDB, keychain);

  if (!(getme)) 
  {
    printf("Ell_GetIntKey: Cannot get structure with key >%s< \n",keychain);
    printf("Ell_GetIntKey: Create first!\n");
    retval = ELLGET_NOKEY;
  }
  else if (getme->type!=CCTK_VARIABLE_INT) 
  {
    CCTK_VWarn(4,__LINE__,__FILE__,CCTK_THORNSTRING,
               "Ell_GetIntKey: Not getting a CCTK_INT value from key %s "
               "(type %d)",
               keychain,getme->type);
    retval = ELLGET_BADTYPE;
  } 
  else if (getme->been_set==ELL_NO) 
  {
    CCTK_VWarn(4,__LINE__,__FILE__,CCTK_THORNSTRING,
               "Ell_GetIntKey: Key %s has not been set to any value (type %d)",
               keychain,getme->type);
    retval = ELLGET_NOTSET;
  }  
  else 
  {
    *value=getme->vals.i;
    retval=0;
  }  

  return(retval);
}

void CCTK_FCALL CCTK_FNAME(Ell_GetIntKey)
     (int *ierr, CCTK_INT *value, ONE_FORTSTRING_ARG) 
{
  ONE_FORTSTRING_CREATE(key)
  *ierr = Ell_GetIntKey(value, key);
  free(key);
}

int Ell_GetStrKey(char **value, const char *keychain) 
{

  struct t_ellthingy *getme=NULL;
  int retval;

  *value = NULL;

  getme = (struct t_ellthingy*)GetNamedData(EllInfoDB, keychain);

  if (!(getme)) 
  {
    CCTK_VWarn(4,__LINE__,__FILE__,CCTK_THORNSTRING,
               "Ell_GetStrKey: Cannot get structure with key %s",keychain);
    retval = ELLGET_NOKEY;
  }
  else if (getme->type!=CCTK_VARIABLE_STRING) 
  {
    CCTK_VWarn(4,__LINE__,__FILE__,CCTK_THORNSTRING,
               "Ell_GetStrKey: Not getting a CCTK_STRING from this key %s "
               "type: %d", keychain,getme->type);
    retval = ELLGET_BADTYPE;
  }
  else if (getme->been_set==ELL_NO) 
  {
    CCTK_VWarn(4,__LINE__,__FILE__,CCTK_THORNSTRING,
               "Ell_GetStrKey: Key %s has not been set", keychain);
    retval = ELLGET_NOTSET;
  }  
  else 
  {
    *value=strdup(getme->vals.s);
    retval=0;
  }
  return(retval);
}

void CCTK_FCALL CCTK_FNAME(Ell_GetStrKey)
     (int *nchar, char **cstring,ONE_FORTSTRING_ARG)
{   
   int i;
  ONE_FORTSTRING_CREATE(fstring)
  ONE_FORTSTRING_PTR(fptr)

  if (strlen(*cstring) > cctk_strlen1) 
  {
    char *message;
    message = (char *)malloc( (200+strlen(*cstring))*sizeof(char) );
    sprintf(message,"Cannot output %s to char* of length %zu\n",
            *cstring,(size_t)cctk_strlen1);
    CCTK_Warn (1,__LINE__,__FILE__,"Cactus",message);
    free(message);
    *nchar = -1;
  }

  for (i=0;i<(int)strlen(*cstring);i++) 
  {
    fptr[i] = (*cstring)[i];
  }

  for (i=(int)strlen(*cstring);i<(int)cctk_strlen1;i++)
  {
    fptr[i] = ' ';
  }

  fptr[strlen(*cstring)] = '\0';

  *nchar = strlen(*cstring);

  free(fstring);
}
