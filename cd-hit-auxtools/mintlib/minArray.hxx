//=================================================================
// This file is a part of the Minimum Template Library.
// By Limin Fu (phoolimin@gmail.com, lmfu@ucsd.edu)
// Released under the MIT license.
//================================================================= 

#ifndef __MIN_ARRAY_HXX__
#define __MIN_ARRAY_HXX__

#include<string.h>
#include<stdint.h>
#include<stdlib.h>
#include<assert.h>
#include<new>
#include"minBase.hxx"

BEGIN_NS_MIN

template<typename item_t>
class Array
{
	item_t *items;
	item_t *buf;
	size_t  size;
	size_t bufsize;

	public:
	Array( size_t n=0, item_t it=item_t() ){
		items = buf = NULL;
		size = bufsize = 0;
		Resize( n, it );
	}
	Array( const Array & other ){
		items = buf = NULL;
		size = bufsize = 0;
		this->operator=( other );
	}
	~Array(){
		Clear();
	}

	size_t Size()const{ return size; }
	item_t* Data(){ return items; }

	item_t& operator[]( size_t i ){ return items[i]; }
	item_t  operator[]( size_t i )const{ return items[i]; }

	Array& operator=( const Array & other ){
		size_t i;
		Clear();
		size = other.size;
		bufsize = other.bufsize;
		items = buf = (item_t*) calloc( bufsize, sizeof(item_t) );
		for(i=0; i<size; i++) new (items+i) item_t( other.items[i] );
		return *this;
	}

	void Swap( Array & other ){
		item_t *it = items, *bu = buf;
		size_t s = size, bs = bufsize;
		items = other.items;
		buf = other.buf;
		size = other.size;
		bufsize = other.bufsize;
		other.items = it;
		other.buf = bu;
		other.size = s;
		other.bufsize = bs;
	}
	void Clear(){
		size_t i;
		for(i=0; i<size; i++) items[i].~item_t();
		if( buf ) free( buf );
		items = buf = NULL;
		size = bufsize = 0;
	}
	void ResetSize( size_t n = 0 ){ if( n > bufsize ) Resize( n ); else size = n; }
	void ResetBuffer(){ // remove buffer from front:
		if( items == buf ) return;
		memmove( buf, items, size*sizeof(item_t) );
		items = buf;
	}
	void Reserve( size_t n, bool extra=false ){
		size_t front = (intptr_t) (items - buf);
		if( bufsize >= (n + front) ) return; // enough space at back
		if( bufsize >= n ){ // front > 0
			ResetBuffer();
			return;
		}
		void *old = buf;
		if( extra ) n += n/5;
		buf = (item_t*) malloc( n*sizeof(item_t) );
		memmove( buf, items, size*sizeof(item_t) );
		items = buf;
		bufsize = n;
		free( old );
	}
	void Resize( size_t n, item_t it=item_t() ){
		size_t i;
		ResetBuffer();
		if( bufsize != n ){
			items = buf = (item_t*) realloc( buf, n*sizeof(item_t) );
			bufsize = n;
		}
		for(i=size; i<n; i++) new (items+i) item_t( it );
		size = bufsize = n;
	}
	item_t& Front(){
		assert( size );
		return items[0];
	}
	item_t& Back(){
		assert( size );
		return items[size-1];
	}
	void PushFront( const item_t & it ){
		size_t front = (intptr_t) (items - buf);
		if( front ){
			items -= 1;
		}else{
			front = bufsize/5 + 5;
			bufsize += front;
			buf = (item_t*) malloc( bufsize*sizeof(item_t) );
			memmove( buf + front, items, size*sizeof(item_t) );
			free( items );
			items = buf + front - 1;
		}
		new (items) item_t( it );
		size += 1;
	}
	void PushBack( const item_t & it ){
		Reserve( size + 1, true );
		new (items+size) item_t( it );
		size += 1;
	}
	void Append( const item_t & it ){ PushBack( it ); }
	void Erase( size_t start, size_t n=-1 ){
		if( size == 0 ) return;
		if( start < 0 or start >= size ) return;
		if( n < 0 ) n = size - start;
		n += start;
		if( n > size ) n = size;

		size_t i;
		for(i=start; i<n; i++){
			items[i].~item_t();
			memset( items+i, 0, sizeof(item_t) );
		}
		for(i=0; i<size-n; i++){
			memcpy( items+start+i, items+n+i, sizeof(item_t) );
		}
		size -= n - start;
	}
	void PopBack(){ if( size ) Erase( size-1, 1 ); }

	void BubbleSort(){
		size_t i, sorted=size;
		while( sorted ){
			size_t swap = 0;
			for(i=1; i<sorted; i++){
				if( items[i-1] > items[i] ){
					item_t tmp = items[i-1];
					items[i-1] = items[i];
					items[i] = tmp;
					swap = i;
				}
			}
			sorted = swap;
		}
	}
	void QuickSort( size_t partial=0, size_t sub=0 ){
		if( size <= 1 ) return;
		if( partial ==0 ) partial = size;
		if( sub ==0 ) sub = size;
		PartialQuickSort( items, 0, sub-1, partial );
	}
	/* Quick Sort.
	 * Adam Drozdek: Data Structures and Algorithms in C++, 2nd Edition.
	 */
	void PartialQuickSort( item_t *data, size_t first, size_t last, size_t partial )
	{
		size_t lower=first+1, upper=last;
		item_t pivot;
		item_t val;
		if( first >= last ) return;
		val = data[first];
		data[first] = data[ (first+last)/2 ];
		data[ (first+last)/2 ] = val;
		pivot = data[ first ];

		while( lower <= upper ){
			while( lower <= last && data[lower] < pivot ) lower ++;
			while( pivot < data[upper] ) upper --;
			if( lower < upper ){
				val = data[lower];
				data[lower] = data[upper];
				data[upper] = val;
				upper --;
			}
			lower ++;
		}
		val = data[first];
		data[first] = data[upper];
		data[upper] = val;
		if( upper && first < upper-1 ) PartialQuickSort( data, first, upper-1, partial );
		if( upper >= partial ) return;
		if( upper+1 < last ) PartialQuickSort( data, upper+1, last, partial );
	}
	item_t QuickXth( item_t *data, size_t first, size_t last, size_t X )
	{
		size_t lower=first+1, upper=last;
		item_t pivot;
		item_t val;
		if( first >= last ) return data[last];
		val = data[first];
		data[first] = data[ (first+last)/2 ];
		data[ (first+last)/2 ] = val;
		pivot = data[ first ];

		while( lower <= upper ){
			while( lower <= last && data[lower] < pivot ) lower ++;
			while( pivot < data[upper] ) upper --;
			if( lower < upper ){
				val = data[lower];
				data[lower] = data[upper];
				data[upper] = val;
				upper --;
			}
			lower ++;
		}
		val = data[first];
		data[first] = data[upper];
		data[upper] = val;
		if( first < upper-1 and upper > X ) return QuickXth( data, first, upper-1, X );
		if( upper+1 < last and upper < X ) return QuickXth( data, upper+1, last, X );
		return val;
	}
	item_t QuickMedian(){
		Array<item_t> tmp( *this );
		return QuickXth( tmp.items, 0, size-1, size/2 );
	}
	item_t QuickXth( size_t xth ){
		Array<item_t> tmp( *this );
		return QuickXth( tmp.items, 0, size-1, xth );
	}
};

END_NS_MIN

#endif
