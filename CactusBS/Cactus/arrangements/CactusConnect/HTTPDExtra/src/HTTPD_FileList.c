#include "HTTPD_FileList.h"
#include "PtrList.h"

/* wrapper functions to make PtrList type-safe for FileList items */

size_t
HTTPD_FileList_NumberOfItems( const FileList * list )
{
        return List_NumberOfItems( list );
}

FileList *
HTTPD_FileList_New()
{
        return (FileList *)List_New();
}

void
HTTPD_FileList_Delete( FileList * list )
{
        List_Delete( list );
}

void
HTTPD_FileList_FreeItemsInListAndEmpty( FileList * list )
{
        List_FreeItemsInListAndEmpty( list );
}

void
HTTPD_FileList_Append( FileList * list, httpFileItem * item )
{
        List_Append( list, item );
}

httpFileItem *
HTTPD_FileList_Item( const FileList * list, size_t index )
{
        return List_Item( list, index );
}

void
HTTPD_FileList_SortAccordingTo( FileList * list,
                HTTPD_FileListSortComparison comparison )
{
        List_SortAccordingTo( list, (ListSortComparison)comparison );
}

int
HTTPD_FileListCompare_Var_Thorn_Slice(
                const httpFileItem * a, const httpFileItem * b )
{
        int order = Compare( a->varname, b->varname );
        if( order == 0 )
        {
                order = Compare( a->thorn, b->thorn );
                if( order == 0 )
                        return Compare( a->slice, b->slice );
                else
                        return order;
        }
        else
                return order;
}
