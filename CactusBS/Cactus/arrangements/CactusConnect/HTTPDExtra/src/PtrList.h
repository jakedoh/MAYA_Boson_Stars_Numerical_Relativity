#ifndef __PTRLIST_HH__
#define __PTRLIST_HH__

#include <stddef.h>

typedef enum { PLFALSE, PLTRUE } PLBOOL;

typedef struct PtrList_tag PtrList;

#ifdef __cplusplus
extern "C"
{
#endif
                                /* Creation and deletion */
PtrList *       List_New( void );
PtrList *       List_NewWithPageSize( size_t pagesize );
PtrList *       List_MakeCopy( const PtrList * other );
void            List_Delete( PtrList * );
                                /* Counts the items */
size_t          List_NumberOfItems( const PtrList * );
                                /* Item accessors */
void *          List_Item( const PtrList *, size_t index );
void            List_SetItem( PtrList *, size_t index, void * ptr );
                                /* List manipulation */
void            List_Append( PtrList *, void * item );
void            List_Insert( PtrList *, size_t index, void * item );
void *          List_RemoveItem( PtrList *, size_t index );
void            List_SwapItems( PtrList *, size_t a_index, size_t b_index );
                                /* Actions on pointer of particular value */
void            List_Remove( PtrList *, void * item );
PLBOOL          List_GetIndexOf( const PtrList *, const void * item,
                                size_t * index );
                                /* Remove all items from list */
void            List_Empty( PtrList * list );
                                /* Copy another list */
void            List_CopyList( PtrList * list, const PtrList * other );
                                /* Special freeing utility */
void            List_FreeItemsInListAndEmpty( PtrList * );
                                /* Sort and Search */
typedef int     (*ListSortComparison)( const void *, const void * );

void            List_SortAccordingTo( PtrList *, ListSortComparison );

typedef PLBOOL  (*ListCondition)( const void * );

void *          List_FirstItemSuchThat( const PtrList *, ListCondition  );
PLBOOL          List_FindFirstIndexSuchThat( const PtrList *,
                                ListCondition condition, size_t * index );
#ifdef __cplusplus
}
#endif

#endif
