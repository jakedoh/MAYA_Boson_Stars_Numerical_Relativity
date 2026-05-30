#include "PtrList.h"
#include <stdlib.h>
#include <memory.h>

#ifdef macintosh
#include <string.h>
#endif

#define ITEMSONFULLPAGE 62
#define TAKEDEFAULTSIZE 0
#define MIN(a,b) (a<b?a:b)

typedef struct ListPage_tag ListPage;

typedef struct PtrList_tag
{
        ListPage *firstPage;
        size_t  itemsPerPage;
} PtrList_PLACEHOLDER;

/* ===========================================================================
 * ____________________________ ListPage _____________________________________
 *  For efficiency, the pointers in the list are arranged into blocks,
 *  called ListPage's.
 * ======================================================================== */

typedef struct ListPage_tag
{
        size_t          numberOfItems;
        size_t          maxItems;
        ListPage *      next;
        void *          *items;
} ListPage_PLACEHOLDER;

static ListPage *
ListPage_New( size_t maxItems )
{
        ListPage *new = (ListPage *)calloc( 1, sizeof( ListPage ) );
        if( new )
        {
                new->items = malloc( maxItems * sizeof( void * ) );
                new->maxItems = maxItems;
        }
        return new;
}

#define ListPage_IsFull( page ) \
        ( ((const ListPage *)page)->numberOfItems >= page->maxItems )


static ListPage *
ListPage_Dup( const ListPage * other );

static void *
ListPage_Append( ListPage * this, void * item );

static void *
ListPage_RemoveItem( ListPage * this, size_t index );

static void *
ListPage_InsertItem( ListPage * this, size_t index, void * item );

                /* ListPage_Dup does not copy the next pointer */
ListPage *
ListPage_Dup( const ListPage * other )
{
        ListPage * this = ListPage_New( other->maxItems );
        if( other->numberOfItems > 0 )
                memcpy( this->items, other->items,
                                other->numberOfItems * sizeof( void * ) );

        this->maxItems = other->maxItems;
        this->numberOfItems = other->numberOfItems;

        return this;
}

        /* Returns the input pointer if it doesn't fit on page 
         * (for similarity with InsertItem) */
void *
ListPage_Append( ListPage *this, void * item )
{
        if( this->numberOfItems < this->maxItems )      /* If page not full */
        {
                *( this->items + this->numberOfItems ) = item; /* item to end*/
                this->numberOfItems++;  /* Update number of items on page */
                return NULL;
        }
        return item;
}

        /* Returns a pointer to the item removed, or NULL if no item found */
void *
ListPage_RemoveItem( ListPage *this, size_t index )
{
        if( index < this->numberOfItems )
        {
                const size_t    lastItemToShift = this->numberOfItems - 1;
                size_t          i;
                void            *item = *( this->items + index );

                this->numberOfItems--;

                for( i = index; i < lastItemToShift; i++ )
                        *( this->items + i ) = *( this->items + i + 1 );

                return item;
        }
        return NULL;
}

        /* Returns a pointer to any item bumped off the page due to
         * the page being full. (could be the input item) */
void *
ListPage_InsertItem( ListPage *this, size_t index, void * item )
{
        void * remainder = NULL;

        if( index < this->numberOfItems )
        {
                const size_t    firstItemToShift = this->numberOfItems - 1;
                size_t          i;

                if( this->numberOfItems >= this->maxItems )
                        remainder = *( this->items + firstItemToShift );

                this->numberOfItems = MIN( this->numberOfItems + 1,
                                this->maxItems );

                for( i = firstItemToShift; i > index; i-- )
                        *( this->items + i ) = *( this->items + i - 1 );

                *( this->items + index ) = item;
        }
        else
                remainder = ListPage_Append( this, item );

        return remainder;
}

/* ===========================================================================
   ______________________________ PtrList functions _________________________

   ======================================================================== */

static void
List_HandleAddressingError(void)
{
        /* To taste */
}

static ListPage * List_Xerox( const PtrList *other);
static void List_FindItemAddress( const PtrList *this, size_t itemNo,
                ListPage ** thePage, size_t *thePageIndex );

        /*_____________________________________________________________________
         **Constructor & Destructor**__________________________________________
         *_____________________________________________________________________
         */

PtrList *
List_New( void )
{
        return List_NewWithPageSize( TAKEDEFAULTSIZE  );
}

PtrList *
List_NewWithPageSize( size_t itemsOnFullPage  )
{
        PtrList *new = (PtrList *)calloc( 1, sizeof( PtrList ) );
        if( new )
        {
                if( itemsOnFullPage == TAKEDEFAULTSIZE )
                        new->itemsPerPage = ITEMSONFULLPAGE;
                else
                        new->itemsPerPage = itemsOnFullPage;
        }
        return new;
}

PtrList *
List_MakeCopy( const PtrList *other )
{
        PtrList *this = (PtrList *)calloc( 1, sizeof( PtrList ) );
        this->itemsPerPage = other->itemsPerPage;
        this->firstPage = List_Xerox( other );
        return this;
}

void
List_Delete( PtrList *list )
{
        List_Empty( list );
        free( list );
}

        /*_____________________________________________________________________
         **Empty**____________________________________________________________
         * Runs through the list of ListPages deleting the (previous) page.
         *_____________________________________________________________________
         */
void
List_Empty( PtrList * this )
{
        if( this->firstPage != NULL )
        {
                ListPage        *prevPage = this->firstPage,
                                *startPage = this->firstPage->next;
                ListPage        *page;

                this->firstPage = NULL;         /* Prefer to do this first,
                                                   so list is consistent */
                for( page = startPage; page != NULL; page = page->next )
                {
                        free( prevPage );
                        prevPage = page;
                }
                free( prevPage );
        }
}
        /*_____________________________________________________________________
        **NumberOfItems**______________________________________________________
        * Runs through the list of ListPages totalling the items on each page.
        * -> This number could be cached.
        *______________________________________________________________________
        */
size_t
List_NumberOfItems( const PtrList *this )
{
        size_t          number = 0;
        ListPage        *page;

        for( page = this->firstPage; page != NULL; page = page->next )
                number += page->numberOfItems;

        return number;
}
        /*_____________________________________________________________________
        **Append**_____________________________________________________________
        *  This and Insert the only PtrList method that allocates ListPages.
        *  If there are no pages with space on them, it allocates an new
        *  page and links the page the the previous last page.
        *  If there is space on the page, it adds the item to the end of
        *  the page and updates the number of items on the page.
        *______________________________________________________________________
        */
void
List_Append( PtrList *this, void * newItem )
{
        ListPage        *page;

        if( this->firstPage == NULL ) /* Add first page if necessary*/
                this->firstPage = ListPage_New( this->itemsPerPage );

        page = this->firstPage;

        while( page != NULL )                   /* Cycle through pages */
        {                                       /* Try to put on page*/
                if( ListPage_Append( page, newItem ) == NULL )
                        return;          /* Go home with a smile */
                else if( page->next == NULL )   /* If next page doesn't exist*/
                        page->next = ListPage_New( this->itemsPerPage );
                                                /* make new page and */
                                                /* link to this*/

                page = page->next;              /* Try the next page*/
        }
        List_HandleAddressingError();
}

        /*_____________________________________________________________________
         * Internal routine for List_Insert.
         * Assumes check has already been done for item address and that
         * it exists on the given page.
         *
         * It is inserted on the page.  The question is, what to do if the
         * page was almost empty, and an item had to be bumped off the page
         * to accomodate it.
         *
         * Could shift all elemtents in the list up by one.  Might be
         * expensive
         *
         * Strategey here will be to add a new page after the present one
         * if the next page is either full or non-existant.  That way, only
         * elements on the next page are affected, and if another item is
         * added to the present page (as often happens) very little more work
         * need be done.
        *______________________________________________________________________
        */
static void
List_InsertExistingItem( PtrList *list, ListPage *page, size_t pageIndex,
                void *item )
{
        void * remnant = ListPage_InsertItem( page, pageIndex, item );
        if( remnant ) /* An item got bumped off when inserting into page */
        {
                ListPage *next = page->next;
                if( next == NULL || ListPage_IsFull( next ) )
                {
                        ListPage *new = ListPage_New( list->itemsPerPage );
                        ListPage_Append( new, item ); /* Item on new */
                        new->next = next;   /* link new page */
                        page->next = new;
                }
                else
                        ListPage_InsertItem( next, 0, remnant );
        }
}
        /*_____________________________________________________________________
        **List_Insert**
        * Inserts the item at the specified index, shifting indices of
        * existing list items as necessary.  If index equals the number
        * of items in the list less one, the item is appended to the list.
        *______________________________________________________________________
        */

void
List_Insert( PtrList * list, size_t index, void * item )
{
        if( index >= List_NumberOfItems( list ) )
                List_Append( list, item );
        else
        {
                ListPage        *page = NULL;
                size_t          pageIndex;

                List_FindItemAddress( list, index, &page, &pageIndex );

                if( page != NULL )
                        List_InsertExistingItem( list, page, pageIndex, item );
                else
                        List_HandleAddressingError();
                /* probably should handle an addressing error */
        }
}

        /*_____________________________________________________________________
        **Item**_______________________________________________________________
        * Returns the item of given index; returns NULL if no such index.
        *______________________________________________________________________
        */

void *
List_Item( const PtrList *this, size_t requestedIndex )
{
        ListPage        *itsPage = NULL;
        size_t          itsPageIndex;

        List_FindItemAddress( this, requestedIndex, &itsPage, &itsPageIndex );

        if( itsPage != NULL )
                return *( itsPage->items + itsPageIndex );

        return NULL;
}
        /*_____________________________________________________________________
        **SetItem**____________________________________________________________
        * Sets the item of given index.  Calls user-supplied addressing
        * handler if no such item exists.
        *______________________________________________________________________
        */
void
List_SetItem( PtrList *this, size_t requestedIndex, void *value )
{
        ListPage        *itsPage = NULL;
        size_t          itsPageIndex;

        List_FindItemAddress( this, requestedIndex, &itsPage, &itsPageIndex );

        if( itsPage != NULL )
                *( itsPage->items + itsPageIndex ) = value;
        else
                List_HandleAddressingError();
}

        /*_____________________________________________________________________
        **CopyList**_________________________________________________________
        * Empties this list and puts items from other list in it.
        *______________________________________________________________________
        */
void
List_CopyList( PtrList * this, const PtrList * other )
{
        List_Empty( this );
        this->itemsPerPage = other->itemsPerPage;
        this->firstPage = List_Xerox(  other );
}
        /*_____________________________________________________________________
        **RemoveItem**_________________________________________________________
        * Removes the item of given index from the list by shifting all the
        * following items back by one and reducing itsNumberOfItems by one.
        * -> It would kind of be nice if this also deleted empty pages
        *______________________________________________________________________
        */
void *
List_RemoveItem( PtrList *this, const size_t requestedIndex )
{
        ListPage        *page = NULL;
        size_t          pageIndex = 0;

        List_FindItemAddress( this, requestedIndex, &page, &pageIndex );

        if( page != NULL )
                return ListPage_RemoveItem( page, pageIndex );
        else
                List_HandleAddressingError();

        return NULL;
}
        /*_____________________________________________________________________
        **Remove**_____________________________________________________________
        * Removes the first item in the list with matching pointer.
        * Note:  No check is done for duplicate pointers in the list.
        *______________________________________________________________________
        */
void
List_Remove( PtrList *this, void * item )
{
        size_t  index = 0;

        if( List_GetIndexOf( this, item, &index ) )
                List_RemoveItem( this, index );
}
        /*_____________________________________________________________________
        **GetIndexOf**_________________________________________________________
        * Inverse of Item(), if pointers in list are unique.
        * Note:  No check is done for duplicate pointers in the list.
        * Returns TRUE and sets the indexi argument to the first item in the
        * list with matching pointer, if such an item exists; otherwise
        * returns FALSE.
        *______________________________________________________________________
        */
PLBOOL
List_GetIndexOf( const PtrList *this, const void * item, size_t * index )
{
        size_t  numItemsOnPreviousPages = 0;
        ListPage        *page = NULL;

        for( page = this->firstPage; page != NULL; page = page->next )
        {
                const size_t    numItemsOnPage = page->numberOfItems;
                size_t  i;

                for( i = 0; i < numItemsOnPage; i++ )
                        if( *( page->items + i ) == item )
                        {
                                *index = numItemsOnPreviousPages + i;
                                return PLTRUE;
                        }
                numItemsOnPreviousPages += numItemsOnPage;
        }
        return PLFALSE;
}
#if 0
static void
List_BubbleSort( PtrList * list, ListSortComparison comparison )
{
        const size_t    n = List_NumberOfItems( list );
        size_t          i, j;

        for( i = 0; i < n; i++ )
        {
                void    *item_i = List_Item( list, i );

                for( j = i + 1; j < n; j++ )
                {
                        void    *item_j = List_Item( list, j );

                        if( comparison( item_i, item_j ) > 0 )
                        {
                                List_SwapItems( list, i, j );
                                item_i = item_j;
                        }
                }
        }
}
#endif
void
List_SwapItems( PtrList * v, size_t a, size_t b )
{
        ListPage        *aPage = NULL, *bPage = NULL;
        size_t          aPageInd, bPageInd;

        List_FindItemAddress( v, a, &aPage, &aPageInd );
        List_FindItemAddress( v, b, &bPage, &bPageInd );

        if( aPage != NULL && bPage != NULL )
        {
                void *temp = *( aPage->items + aPageInd );
                *( aPage->items + aPageInd ) = *( bPage->items + bPageInd );
                *( bPage->items + bPageInd ) = temp;
        }
}
        /* Adapted from Kernighan & Ritchie */
static void
List_Qsort( PtrList * v, size_t left, size_t right,
                          ListSortComparison comparison )
{
        size_t i, last;

        if( left >= right )     /* do nothing if array contains */
                return;         /* fewer than two elements */
        List_SwapItems( v, left, ( left + right ) / 2); /* move partition */
        last = left;                            /* elem to v[0] */
        for( i = left + 1; i <= right; i++ )
        {/* SW this could be improved a lot */
                void *  item_i = List_Item( v, i );
                void *  item_left = List_Item( v, left );
                if( comparison( item_i, item_left ) < 0 )
                        List_SwapItems( v, ++last, i );
        }
        List_SwapItems( v, left, last );        /* restore partition elem */
/* SW NEED TO THINK ABOUT THIS */
        if( last > 0 )  /* K&R use int, so this isn't an issue for them */
                List_Qsort( v, left, last - 1, comparison );
        List_Qsort( v, last + 1, right, comparison );
}

        /*_____________________________________________________________________
        **SortAccordingTo**__________________________________________________
        * Based on the comparison of any two items provided by the comparison 
        * function argument, sorts the list.
        * Note here it uses a bubble sort, which is slow but easy to 
        * understand.  Want something faster?  Write your own!
        *______________________________________________________________________
        */
void
List_SortAccordingTo( PtrList * list, ListSortComparison comparison )
{
        size_t n = List_NumberOfItems( list );
        if( n > 1 )
                List_Qsort( list, 0, n - 1, comparison );
}

        /*_____________________________________________________________________
        **FirstItemSuchThat**__________________________________________________
        * Takes a ListCondition function, returns the pointer value that 
        * satisfies the condition, otherwise returns NULL,
        *______________________________________________________________________
        */
void *
List_FirstItemSuchThat( const PtrList *list, ListCondition condition )
{
        const size_t    n = List_NumberOfItems( list );
        size_t          i;

        for( i = 0; i < n; i++ )
                if( condition( List_Item( list, i ) ) )
                        return List_Item( list, i );

        return NULL;
}

        /*_____________________________________________________________________
        **FindFirstIndexSuchThat**_____________________________________________
        * Takes a ListCondition function, sets the index pointer to the index
        * of the first pointer in the list that satisfies condition.
        * Otherwise returns PLFALSE.
        *______________________________________________________________________
        */
PLBOOL
List_FindFirstIndexSuchThat( const PtrList *list, ListCondition condition,
                                        size_t * index )
{
        const size_t    n = List_NumberOfItems( list );
        size_t          i;

        for( i = 0; i < n; i++ )
                if( condition( List_Item( list, i ) ) )
                {
                        *index = i;
                        return PLTRUE;
                }

        return PLFALSE;
}

        /*_____________________________________________________________________
        **FreeItemsInListAndEmpty**____________________________________________
        * Lists don't usually delete the things they refer to.
        * This function takes the list and uses it to free the listed items.
        *______________________________________________________________________
        */
void
List_FreeItemsInListAndEmpty( PtrList * this )
{
        size_t  i = List_NumberOfItems( this );
        void    *item = NULL;

        while( i > 0 )  /* Remove items in reverse order for efficeincy */
        {
                i--;
                item = List_Item( this, i );
                List_RemoveItem( this, i ); /* Remove item before deleting */
                                                /* (for consistency) */
                free( item );
        }
}
        /*_____________________________________________________________________
        **FindItemAddress**____________________________________________________
        * Private utility that gets the page and index relative to the page
        * given the item index.  Sets *thePage to NULL if no such index.
        *______________________________________________________________________
        */
void
List_FindItemAddress( const PtrList *this, size_t itemNo,
                ListPage ** thePage, size_t *thePageIndex ) 
{
        ListPage *page = NULL;
        size_t  numItemsIncludingThisPage = 0,
                numItemsOnPreviousPages = 0;

        *thePage = NULL;
        *thePageIndex = 0;

        for( page = this->firstPage; page != NULL; page = page->next )
        {
                numItemsIncludingThisPage += page->numberOfItems;

                if( numItemsIncludingThisPage > itemNo )
                {
                        *thePage = page;
                        *thePageIndex = itemNo - numItemsOnPreviousPages;
                        return;
                }
                numItemsOnPreviousPages = numItemsIncludingThisPage;
        }
}

        /*_____________________________________________________________________
        **Xerox**______________________________________________________________
        * Private utility that makes a duplicate of the pages in the list.
        *______________________________________________________________________
        */
ListPage *
List_Xerox( const PtrList *list)
{
        ListPage *copy = NULL, *newPage = NULL, *page = list->firstPage;

        if( list->firstPage != NULL )
                copy = newPage = ListPage_Dup( list->firstPage );

        if( newPage != NULL )
                do
                {
                        if( page->next != NULL )
                        {
                                page = page->next;
                                newPage->next = ListPage_Dup( page );
                        }
                        newPage = newPage->next;
                }
                while( newPage != NULL );
        return copy;
}
