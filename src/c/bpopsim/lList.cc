#include "lList.h"

// ============================================================================

// queue -- this function puts an item at the end of the list
void lList::queue( void *pData )
{
	Node *pNewNode;

	// create the new node
	pNewNode = new Node( pData );

	// place the node in the list
	if( pHead == NULL )
	{
		// list is empty
		pHead = pNewNode;
		pTail = pNewNode;
	} // end if -- empty list
	else
	{
		// list has at least one item in it
		pTail->SetNext(pNewNode );
		pTail = pNewNode;
	} // end else
} // end queue

// ----------------------------------------------------------------------------

// push - this function puts an item at the front of the list
void lList::push( void *pData )
{
	Node *pNewNode;
	
	// create the new node
	pNewNode = new Node( pData );

	// place the node on the list
	pNewNode->SetNext( pHead );
	pHead = pNewNode;

	// if the list was empty then this is both the head and tail 
	if( pTail == NULL )
		pTail = pNewNode;

} // end push

// ----------------------------------------------------------------------------

// pop - this function removes an item from the front of the list
void *lList::pop()
{
	void *pData;

	// if the list is empty return NULL
	if( pHead == NULL )
		return( NULL );

	// get the head node's data
	pData = pHead->Data();

	// if there is only one item in the list 
	if( pHead == pTail )
	{
		// remove old node
		delete pHead;

		// set head and tail pointers
		pHead = NULL;
		pTail = NULL;
	} // end if -- one item in list
	else
	{
		// remove old node
		delete pHead;
	
		// set head pointer
		pHead = pHead->Next();
	} // end else -- were at least two items in list
	
	return( pData );
} // end pop

// ----------------------------------------------------------------------------

// destroy - deletes entire list
void lList::destroy()
{
	void *pData;
	Node *pNode, *pTempNode;

	pNode = pHead;
	while( pNode != NULL )
	{
		pTempNode = pNode->Next();
		pData = pNode->Data();
		delete pData;
		delete pNode;
		
		pNode = pTempNode;
	} // end while

	pHead = NULL;
	pTail = NULL;
	pCurrent = NULL;
} // end destroy

// ----------------------------------------------------------------------------

// next - sets the current pointer to the next item and returns a pointer to
//        its data object; if current pointer is NULL then returns the first
//        item in the list returns NULL if no next item
//
void *lList::next()
{
	// if the list is empty return NULL
	if( pHead == NULL )
		return( NULL );

	// change the current pointer
	if( pCurrent == NULL )
		pCurrent = pHead;
	else
		pCurrent = pCurrent->Next();

	if( pCurrent == NULL )
		return( NULL );

	return( pCurrent->Data() );
} // end next

// ----------------------------------------------------------------------------

