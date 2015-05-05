//=================================================================
// This file is a part of the Minimum Template Library.
// By Limin Fu (phoolimin@gmail.com, lmfu@ucsd.edu)
// Note: Adapted from Dao Virtual Machine source codes (daovm.net).
// Released under the MIT license.
//================================================================= 

#include"minMap.hxx"

BEGIN_NS_MIN

DNode* DNode::First()const
{
	DNode *p = (DNode*) this;
	if( p ) while( p->left ) p = p->left;
	return p;
}
DNode* DNode::Next()const
{
	DNode* p = (DNode*) this->right;
	if( this->right ){
		while( p->left ) p = p->left;
	}else if( this->parent ){
		if( this == this->parent->left ){
			p = this->parent;
		}else{
			p = this->parent;
			while( p->parent && p==p->parent->right ) p = p->parent;
			p = p->parent;
		}
	}
	return p;
}

DTree::DTree( bool hs )
{
	root = NULL;
	table = NULL;
	size = 0;
	tsize = 0;
	hashing = hs;
	if( hashing ){
		tsize = 4;
		table = (DNode**) calloc( tsize, sizeof(DNode*) );
	}
}
DTree::~DTree()
{
	Clear();
	if( table ) free( table );
}
void DTree::Clear()
{
	int i;
	if( hashing ){
		for(i=0; i<tsize; i++) DeleteTree( table[i] );
		if( table ) free( table );
		tsize = 4;
		table = (DNode**) calloc( tsize, sizeof(DNode*) );
	}else{
		DeleteTree( root );
	}
	root = NULL;
	size = 0;
}
DNode* DTree::Insert( DNode *node )
{
	DNode *p;
	int id = 0;
	if( hashing ){
		node->HashKey( this->tsize );
		id = node->hash;
		this->root = this->table[id];
		if( this->root ==NULL ){
			this->size += 1;
			this->table[id] = node;
			node->color = RB_BLACK;
			return node;
		}
	}
	p = SimpleInsert( node );
	if( p == node ){ /* key not exist: */
		InsertNode( node );
		if( hashing ){
			this->table[id] = this->root;
			if( this->size >= this->tsize ) ResetTable();
		}
	}else{
		p->CopyValue( node );
		delete node;
	}
	return p;
}
DNode* DTree::First()
{
	const DTree *tree = this;
	return (DNode*) tree->First();
}
DNode* DTree::Next( const DNode *node )
{
	const DTree *tree = this;
	return (DNode*) tree->Next( node );
}
const DNode* DTree::First()const
{
	const DNode *node = NULL;
	int i = 0;
	if( hashing ){
		while( i < this->tsize && this->table[i] == NULL ) i += 1;
		if( i < this->tsize ) node = this->table[i]->First();
	}
	if( node == NULL && this->root ) node = this->root->First();
	return node;
}
const DNode* DTree::Next( const DNode *node )const
{
	const DNode *next;
	if( node == NULL ) return NULL;
	next = node->Next();
	if( next == NULL && hashing ){
		int i = node->hash + 1;
		while( i < this->tsize && this->table[i] == NULL ) i += 1;
		if( i < this->tsize ) next = this->table[i]->First();
	}
	return next;
}

#define HASH_SEED  0xda0

void DTree::DeleteTree( DNode *node )
{
	if( node == NULL ) return;
	DeleteTree( node->left );
	DeleteTree( node->right );
	delete node;
}
void DTree::RotateLeft( DNode *child )
{
	DNode *grandpa = child->parent;
	DNode *parent = child->right;
	
	if( grandpa ){
		if( child == grandpa->right )
			grandpa->right = parent;
		else
			grandpa->left = parent;
	}else{
		this->root = parent;
	}
	parent->parent = grandpa;
	
	child->right = parent->left;
	if( child->right ) child->right->parent = child;
	
	parent->left = child;
	child->parent = parent;
}
void DTree::RotateRight( DNode *parent )
{
	DNode *grandpa = parent->parent;
	DNode *child = parent->left;
	
	if( grandpa ){
		if( parent == grandpa->right )
			grandpa->right = child;
		else
			grandpa->left = child;
	}else{
		this->root = child;
	}
	child->parent = grandpa;
	
	parent->left = child->right;
	if( parent->left ) parent->left->parent = parent;
	
	child->right = parent;
	parent->parent = child;
}
void DTree::InsertNode( DNode *node )
{
	DNode *grandpa, *parent, *uncle;
	
	node->color = RB_RED;
	this->size ++;
	
	while( node->parent != NULL ){
		parent = node->parent;
		grandpa = parent->parent;
		if( parent->color == RB_RED ){ /* insert_case2() */
			/* grandpa can't be NULL, since parent is RED and can't be root. */
			uncle = ( parent == grandpa->left ? grandpa->right : grandpa->left );
			if( uncle != NULL && uncle->color == RB_RED ){ /* insert_case3(): */
				parent->color = RB_BLACK;
				uncle->color  = RB_BLACK;
				grandpa->color = RB_RED;
				node = grandpa;
				continue; /* insert_case1() */
			}else{
				if( node == parent->right && parent == grandpa->left ){
					RotateLeft( parent );
					node = node->left;
				}else if( node == parent->left && parent == grandpa->right ){
					/* rotate right around parent: */
					RotateRight( parent );
					node = node->right;
				}
				/* node changed, update parent and grandpa. */
				parent = node->parent;
				grandpa = parent->parent;
				/* insert_case5() */
				
				parent->color = RB_BLACK;
				grandpa->color = RB_RED;
				if( node == parent->left && parent == grandpa->left )
					RotateRight( grandpa );
				else
					RotateLeft( grandpa );
			}
		}
		break;
	}
	/* insert_case1() as in Wikipedia term: Red-black tree. */
	if( node->parent == NULL){
		this->root = node;
		node->color = RB_BLACK;
	}
}
void DTree::EraseChild( DNode *node )
{
	DNode *extreme = node;
	DNode *child = 0;
	
	if( node == NULL ) return;
	this->size --;
	
	/* deletion by coping */
	if( node->left ){
		extreme = node->left;
		while( extreme->right ) extreme = extreme->right;
		child = extreme->left;
	}else if( node->right ){
		extreme = node->right;
		while( extreme->left ) extreme = extreme->left;
		child = extreme->right;
	}
	SwapNode( node, extreme );
	
	if( child ){
		/* replace node */
		child->parent = extreme->parent;
		if( extreme->parent ){
			if( extreme == extreme->parent->left )
				extreme->parent->left = child;
			else
				extreme->parent->right = child;
		}
		if( extreme->color == RB_BLACK ){
			if( child->color == RB_RED )
				child->color = RB_BLACK;
			else{
				node = child;
				while( node->parent ){ /* delete_case1() */
					
					DNode *parent = node->parent;
					DNode *sibling = ( node == parent->left ? parent->right : parent->left );
					if( sibling && sibling->color == RB_RED ){ /* delete_case2() */
						parent->color = RB_RED;
						sibling->color = RB_BLACK;
						if( node == parent->left )
							RotateLeft( parent );
						else
							RotateRight( parent );
					}
					/* node relationship changed, update parent and sibling: */
					parent = node->parent;
					sibling = ( node == parent->left ? parent->right : parent->left );
					if( ! sibling ) break;
					/* delete_case3() */
					if( parent->color == RB_BLACK && sibling->color == RB_BLACK
						&& ( ! sibling->left || sibling->left->color == RB_BLACK )
						&& ( ! sibling->right|| sibling->right->color == RB_BLACK)){
						sibling->color = RB_RED;
						node = node->parent;
						continue; /* delete_case1() */
					}else{
						/* delete_case4() */
						if( parent->color == RB_RED && sibling->color == RB_BLACK
							&& ( ! sibling->left || sibling->left->color == RB_BLACK )
							&& ( ! sibling->right|| sibling->right->color == RB_BLACK)){
							sibling->color = RB_RED;
							parent->color = RB_BLACK;
						}else{
							/* delete_case5() */
							if( node == parent->left && sibling->color == RB_BLACK
								&&( sibling->left && sibling->left->color == RB_RED )
								&&( !sibling->right|| sibling->right->color == RB_BLACK)){
								sibling->color = RB_RED;
								sibling->left->color = RB_BLACK;
								RotateRight( sibling );
							}else if( node == parent->right && sibling->color == RB_BLACK
								&&( sibling->right && sibling->right->color == RB_RED )
								&&( !sibling->left || sibling->left->color == RB_BLACK)){
								sibling->color = RB_RED;
								sibling->right->color = RB_BLACK;
								RotateLeft( sibling );
							}
							/* node relationship changed, update parent and sibling: */
							parent = node->parent;
							sibling = ( node==parent->left ? parent->right:parent->left );
							/* delete_case6() */
							sibling->color = parent->color;
							parent->color = RB_BLACK;
							if( node == parent->left ){
								sibling->right->color = RB_BLACK;
								RotateLeft( parent );
							}else{
								sibling->left->color = RB_BLACK;
								RotateRight( parent );
							}
						} /* end if */
					} /* end if */
				} /* end while */
			}
		}
	}else if( extreme->parent ){
		if( extreme == extreme->parent->left )
			extreme->parent->left = NULL;
		else
			extreme->parent->right = NULL;
	}else{
		this->root = NULL;
	}
	delete extreme;
}
void DTree::EraseNode( DNode *node )
{
	if( node == NULL ) return;
	if( this->hashing ){
		unsigned int hash = node->hash;
		this->root = this->table[ node->hash ];
		if( this->root == NULL ) return;
		EraseChild( node );
		this->table[ hash ] = this->root;
		if( this->size < 0.25*this->tsize ) ResetTable();
	}else{
		EraseChild( node );
	}
}
void DTree::InsertTree( DNode *node )
{
	DNode *left = node->left;
	DNode *right = node->right;
	node->HashKey( this->tsize );
	node->parent = node->left = node->right = NULL;
	this->root = this->table[ node->hash ];
	if( this->root == NULL ){
		node->color = RB_BLACK;
		this->table[ node->hash ] = node;
		this->size += 1;
	}else{
		SimpleInsert( node );
		InsertNode( node );
		this->table[ node->hash ] = this->root;
	}
	if( left ) InsertTree( left );
	if( right ) InsertTree( right );
}
void DTree::ResetTable()
{
	DNode **nodes = this->table;
	int i, tsize = this->tsize;
	
	if( hashing ==0 ) return;
	this->tsize = 2 * this->size + 1;
	this->table = (DNode**)calloc( this->tsize, sizeof(DNode*) );
	this->size = 0;
	for(i=0; i<tsize; i++) if( nodes[i] ) InsertTree( nodes[i] );
	if( nodes ) free( nodes );
}
DNode* DTree::SimpleInsert( DNode *node )
{
	DNode *p = this->root;
	int compare;
	node->color = RB_RED;
	if( this->root == NULL ) return node;
	for(;;){
		compare = CompareNode( node, p );
		if( compare == 0 ){
			return p;
		}else if( compare < 0 ){
			if( p->left ){
				p = p->left;
			}else{
				p->left = node;
				node->parent = p;
				break;
			}
		}else{
			if( p->right ){
				p = p->right;
			}else{
				p->right = node;
				node->parent = p;
				break;
			}
		}
	}
	return node;
}

static unsigned int MurmurHash2_Int32( unsigned int k )
{
	/* 'm' and 'r' are mixing constants generated offline.
	They're not really 'magic', they just happen to work well. */
	const unsigned int m = 0x5bd1e995;
	const int r = 24;
	
	/* Initialize the hash to a 'random' value */
	unsigned int h = HASH_SEED ^ 4;
	
	/* Mix 4 bytes into the hash */
	k *= m; 
	k ^= k >> r; 
	k *= m; 
	
	h *= m; 
	h ^= k;
	
	/* Do a few final mixes of the hash to ensure the last few
	       bytes are well-incorporated. */
	h ^= h >> 13;
	h *= m;
	h ^= h >> 15;
	
	return h;
}
unsigned int MurmurHash2( const void *key, int len, unsigned int seed );

unsigned int MakeHash( unsigned int k )
{
	return MurmurHash2_Int32( k );
}
unsigned int MakeHash( int k )
{
	return MurmurHash2_Int32( (unsigned int) k );
}
unsigned int MakeHash( void *p )
{
	return MurmurHash2( & p, sizeof(void*), HASH_SEED );
}
unsigned int MakeHash( uint64_t k )
{
	return MurmurHash2( & k, sizeof(uint64_t), HASH_SEED );
}
int Compare( int left, int right )
{
	return left == right ? 0 : (left < right ? -1 : 1);
}
int Compare( unsigned int left, unsigned int right )
{
	return left == right ? 0 : (left < right ? -1 : 1);
}
int Compare( void *left, void *right )
{
	if( left == right ) return 0;
	if( left < right ) return -1;
	return 1;
}
int Compare( uint64_t left, uint64_t right )
{
	return left == right ? 0 : (left < right ? -1 : 1);
}


/*
PUBLIC DOMAIN CODES
http://sites.google.com/site/murmurhash/
http://www.burtleburtle.net/bob/hash/doobs.html
*/

/* -----------------------------------------------------------------------------
   MurmurHash2, by Austin Appleby

   Note - This code makes a few assumptions about how your machine behaves -
   1. We can read a 4-byte value from any address without crashing
   2. sizeof(int) == 4

   And it has a few limitations -
   1. It will not work incrementally.
   2. It will not produce the same results on little-endian and big-endian
   machines.
 */
unsigned int MurmurHash2( const void *key, int len, unsigned int seed )
{
	/* 'm' and 'r' are mixing constants generated offline.
	They're not really 'magic', they just happen to work well. */
	const unsigned int m = 0x5bd1e995;
	const int r = 24;
	
	/* Initialize the hash to a 'random' value */
	unsigned int h = seed ^ len;
	
	/* Mix 4 bytes at a time into the hash */
	const unsigned char * data = (const unsigned char *)key;
	
	while(len >= 4) {
		unsigned int k = *(unsigned int *)data;
		
		k *= m; 
		k ^= k >> r; 
		k *= m; 
		
		h *= m; 
		h ^= k;
		
		data += 4;
		len -= 4;
	}
	
	/* Handle the last few bytes of the input array */
	switch(len)
	{
		case 3: h ^= data[2] << 16;
		case 2: h ^= data[1] << 8;
		case 1: h ^= data[0];
		h *= m;
	};
	
	/* Do a few final mixes of the hash to ensure the last few
	bytes are well-incorporated. */
	h ^= h >> 13;
	h *= m;
	h ^= h >> 15;
	
	return h;
}

END_NS_MIN
