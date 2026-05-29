#ifndef _ELL_DBSTRUCTURE_H_
#define _ELL_DBSTRUCTURE_H_



#define ELL_NO  0
#define ELL_YES 1
#define ELL_ISNOKEY -1

#define ELLCREATE_FAILED  -1
#define ELLCREATE_TWICE    1
 
#define ELLGET_NOKEY      -1
#define ELLGET_BADTYPE    -2
#define ELLGET_NOTSET     -3

#define ELLSET_FAILED     -1
#define ELLSET_BADTYPE    -2


#ifdef CCODE

#ifdef __cplusplus 
extern "C" {
#endif

int Ell_CreateKey(int vartype, const char *keychain);
int Ell_IsKey(const char *keychain);
int Ell_UnsetKey(const char *keychain);
int Ell_DeleteKey(const char *keychain);

int Ell_SetRealKey(CCTK_REAL value, const char *keychain);
int Ell_SetIntKey(CCTK_INT value, const char *keychain);
int Ell_SetStrKey(char *value, const char *keychain);

int Ell_GetRealKey(CCTK_REAL *value, const char *keychain);
int Ell_GetIntKey(CCTK_INT *value,const char *keychain);
int Ell_GetStrKey(char **value, const char *keychain);


#ifdef __cplusplus 
}
#endif

#endif /* CCODE */



#endif /* _ELL_DBSTRUCTURE_H_ */
