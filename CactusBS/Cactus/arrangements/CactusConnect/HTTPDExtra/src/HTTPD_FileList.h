#ifndef _HTTPDFILELIST_H_
#define _HTTPDFILELIST_H_ 1

/* SW temporary explicit paths while testing SString module */
#include "CactusConnect/HTTPD/src/SString.h"
#include "CactusConnect/HTTPD/src/SString_Namespace.h"

typedef struct httpFileItem_tag
{
  const String *thorn;
  const String *varname;
  const String *mimetype;
  const String *slice;
  const String *description;
  const String *filename;
  const String *linkname;
} httpFileItem;

typedef enum {FLFALSE, FLTRUE} FLBOOL;

/* a wrapper module around PtrList for type safety */

typedef struct PtrList_tag FileList;

#ifdef __cplusplus
extern "C"
{
#endif

size_t          HTTPD_FileList_NumberOfItems( const FileList * );
FileList *      HTTPD_FileList_New( void );
void            HTTPD_FileList_Delete( FileList * );

void            HTTPD_FileList_FreeItemsInListAndEmpty( FileList * );
void            HTTPD_FileList_Append( FileList *, httpFileItem * item );
httpFileItem *  HTTPD_FileList_Item( const FileList *, size_t index );

typedef int     (*HTTPD_FileListSortComparison)( const httpFileItem *,
                                                const httpFileItem * );
int             HTTPD_FileListCompare_Var_Thorn_Slice( const httpFileItem *,
                                                const httpFileItem * );

void            HTTPD_FileList_SortAccordingTo( FileList *, HTTPD_FileListSortComparison comparison );

#ifdef __cplusplus
}
#endif


#define NumberOfFiles( a ) \
                HTTPD_FileList_NumberOfItems( a )
#define FileItem( a, b ) \
                HTTPD_FileList_Item( a, b )
#define AppendFile( a, b ) \
                HTTPD_FileList_Append( a, b )
#define SortFilesAccordingTo( a, b ) \
                HTTPD_FileList_SortAccordingTo( a, b )
#define FreeFileItemsInListAndEmpty( a ) \
                HTTPD_FileList_FreeItemsInListAndEmpty( a )
#define Variable_Thorn_Slice \
                HTTPD_FileListCompare_Var_Thorn_Slice
#endif
