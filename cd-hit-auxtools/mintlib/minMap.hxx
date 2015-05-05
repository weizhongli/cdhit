//=================================================================
// This file is a part of the Minimum Template Library.
// By Limin Fu (phoolimin@gmail.com, lmfu@ucsd.edu)
// Note: Adapted from Dao Virtual Machine source codes (daovm.net).
// Released under the MIT license.
//================================================================= 

#ifndef __MIN_MAP_HXX__
#define __MIN_MAP_HXX__

#include<stdlib.h>
#include<string.h>
#include<stdint.h>
#include"minBase.hxx"

BEGIN_NS_MIN

#define RB_RED    0
#define RB_BLACK  1

typedef enum{ KEY_EQ=0, KEY_ML/*Max Less*/, KEY_MG/*Min Great*/ } KeySearchType;

struct DNode
{
	unsigned int color :  1;
	unsigned int hash  : 31;
	
	DNode  *parent;
	DNode  *left;
	DNode  *right;
	DNode  *next;
	
	DNode(){ parent = left = right = NULL; }
	virtual ~DNode(){}
	
	DNode* First()const;
	DNode* Next()const;
	
	virtual void HashKey( int T ){}
	virtual void CopyValue( DNode *other ){}
};

struct DTree
{
	DNode **table;
	DNode  *root;
	DNode  *first;
	DNode  *last;
	int     size;
	int     tsize;
	bool    hashing;
	
	DTree( bool hashing = false );
	virtual ~DTree();
	
	void Clear();
	int Size()const{ return size; }
	
protected:
	
	virtual int CompareNode( DNode *left, DNode *right )=0;
	virtual void SwapNode( DNode *node, DNode *extreme ){}
	
	DNode* First();
	DNode* Next( const DNode *node );
	const DNode* First()const;
	const DNode* Next( const DNode *node )const;

	DNode* Insert( DNode *node );
	DNode* SimpleInsert( DNode *node );
	
	void DeleteTree( DNode *node );
	void RotateLeft( DNode *child );
	void RotateRight( DNode *parent );
	void InsertNode( DNode *node );
	void EraseChild( DNode *node );
	void EraseNode( DNode *node );
	void InsertTree( DNode *node );
	void ResetTable();
};

unsigned int MakeHash( unsigned int k );
unsigned int MakeHash( int k );
unsigned int MakeHash( void *p );
unsigned int MakeHash( uint64_t k );
int Compare( int left, int right );
int Compare( unsigned int left, unsigned int right );
int Compare( void *left, void *right );
int Compare( uint64_t left, uint64_t right );


template<typename KeyType, typename ValueType>
struct DKeyValue : public DNode
{
	KeyType   key;
	ValueType value;
	
	DKeyValue( const KeyType & k, const ValueType & v ){
		key = k;
		value = v;
		parent = left = right = next = NULL;
	}
	void HashKey( int T ){ hash = MakeHash( key ) % T; }
	void CopyValue( DNode *other ){ value = ((DKeyValue*)other)->value; }
};


template<typename KeyType, typename ValueType>
struct DMap : public DTree
{
	typedef Min::DKeyValue<KeyType,ValueType> DKeyValue;
	
	DMap(){}
	DMap( const DMap & other ) : DTree( other.hashing ){ this->operator=( other ); }
	
	DMap& operator=(const DMap & other){
		const DKeyValue *node;
		Clear();
		hashing = other.hashing;
		for(node=other.First(); node; node=other.Next(node)) Insert( node->key, node->value );
		return *this;
	}
	
	void Erase( const KeyType & key ){
		EraseNode( FindNode( key, KEY_EQ ) );
	}
	DKeyValue* Insert( const KeyType & key, const ValueType & value ){
		DKeyValue *node = new DKeyValue( key, value );
		return (DKeyValue*) this->DTree::Insert( (DNode*) node );
	}
	DKeyValue* Find( const KeyType & key ){ return FindNode( key, KEY_EQ ); }
	DKeyValue* FindML( const KeyType & key ){ return FindNode( key, KEY_ML ); }
	DKeyValue* FindMG( const KeyType & key ){ return FindNode( key, KEY_MG ); }
	
	DKeyValue* First(){ return (DKeyValue*) DTree::First(); }
	DKeyValue* Next( const DKeyValue *node ){ return (DKeyValue*) DTree::Next( node ); }
	const DKeyValue* First()const{ return (DKeyValue*) DTree::First(); }
	const DKeyValue* Next( const DKeyValue *node )const{ return (DKeyValue*) DTree::Next( node ); }
	
	ValueType& operator[]( KeyType key ){
		DKeyValue *node = Find( key );
		if( node ) return node->value;
		node = Insert( key, ValueType() );
		return node->value;
	}
	
protected:
	
	DMap( bool hashing ) : DTree( hashing ){}
	
	virtual int CompareNode( DNode *left0, DNode *right0 ){
		DKeyValue *left = (DKeyValue*) left0;
		DKeyValue *right = (DKeyValue*) right0;
		return Compare( left->key, right->key );
	}
	void SwapNode( DNode *node0, DNode *extreme0 )
	{
		DKeyValue *node = (DKeyValue*) node0;
		DKeyValue *extreme = (DKeyValue*) extreme0;
		KeyType key = extreme->key;
		ValueType value = extreme->value;
		int hash = extreme->hash;
		extreme->hash = node->hash;
		extreme->key = node->key;
		extreme->value = node->value;
		node->hash = hash;
		node->key = key;
		node->value = value;
	}
	DKeyValue* FindChild( DKeyValue *root, KeyType key, KeySearchType type )
	{
		DKeyValue *p = root;
		DKeyValue *m = 0;
		int compare;
		
		if( root == NULL ) return NULL;
		
		for(;;){
			compare = Compare( key, p->key );
			if( compare == 0 ) return p;
			
			if( compare < 0 ){
				if( type == KEY_MG ) m = p;
				if( p->left ) p = (DKeyValue*) p->left; else break;
			}else{
				if( type == KEY_ML ) m = p;
				if( p->right ) p = (DKeyValue*) p->right; else break;
			}
		}
		return m;
	}
	DKeyValue* FindNode( KeyType key, KeySearchType type )
	{
		DKeyValue *root = (DKeyValue*) this->root;
		int id;
		if( hashing ){
			id = MakeHash( key ) % this->tsize;
			root = (DKeyValue*) this->table[id];
			if( root == NULL ) return NULL;
		}
		return FindChild( root, key, type );
	}
};

template<typename KeyType,typename ValueType>
struct DHash : public DMap<KeyType,ValueType>
{
	DHash() : DMap<KeyType,ValueType>( true ){}
	DHash( const DHash & other ) : DMap<KeyType,ValueType>( other ){}
};

#define Map  DMap
#define Hash DHash
#define Node DKeyValue

END_NS_MIN

#endif
