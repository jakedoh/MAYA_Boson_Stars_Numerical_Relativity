#ifndef __PTRLIST_NAMESPACE_HH__
#define __PTRLIST_NAMESPACE_HH__

#include "PtrList.h"

/* A poor man's namespace for the PtrList module */

#define NumberOfItems( a ) \
                List_NumberOfItems( a )
#define Remove( a, b ) \
                List_Remove( a, b )
#define Item( a, b ) \
                List_Item( a, b )
#define Insert( a, b, c ) \
                List_Insert( a, b, c )
#define SetItem( a, b, c ) \
                List_SetItem( a, b, c )
#define Append( a, b ) \
                List_Append( a, b )
#define RemoveItem( a, b ) \
                List_RemoveItem( a, b )
#define SwapItems( a, b, c ) \
                List_SwapItems( a, b, c )
#define Remove( a, b ) \
                List_Remove( a, b )
#define GetIndexOf( a, b, c ) \
                List_GetIndexOf( a, b, c )
#define DeleteItemsAndEmpty( a ) \
                List_DeleteItemsAndEmpty( a )
#define SortAccordingTo( a, b ) \
                List_SortAccordingTo( a, b )
#define FirstItemSuchThat( a, b ) \
                List_FirstItemSuchThat( a, b )
#define FindFirstIndexSuchThat( a, b ) \
                List_FindFirstIndexSuchThat( a, b )

#endif
