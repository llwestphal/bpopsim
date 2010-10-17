#include <stdio.h>

class Node {

public:
	Node( void *_pData )
	{
		pData = _pData;
		pNext = NULL;
	} // end constructor

	Node *Next() { return( pNext ); }
	void SetNext( Node *pNode ) { pNext = pNode; }
	void *Data() { return( pData ); }

private:
	void *pData;
	Node *pNext;
}; 
class lList
{
public:
	lList()
	{
		pHead = NULL;
		pTail = NULL;
		pCurrent = NULL;
	}
	~lList() { this->destroy(); }

	void push( void *pData );
	void queue( void *pData );
	void *pop();
	void destroy();

	void top() { pCurrent = NULL; }
	void *first()
	{
		if( pHead == NULL )
			return( NULL );

		return( pHead->Data() );
	} // end first

	void *last()
	{
		if( pTail == NULL )
			return( NULL );

		return( pTail->Data() ); 
	} // end last

	void *next();

private:
	Node *pHead;
	Node *pTail;
	Node *pCurrent;
};

