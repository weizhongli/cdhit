//=================================================================
// This file is a part of cdhit-dup.
// By Limin Fu (phoolimin@gmail.com, lmfu@ucsd.edu)
//================================================================= 

#include <stdint.h>
#include <stdio.h>
#include "minArray.hxx"
#include "minString.hxx"

#define BEGIN_NS_BIO namespace Bio {
#define END_NS_BIO }

BEGIN_NS_BIO

class Sequence
{
	uint32_t     id;
	Min::String  des;
	Min::String  seq;
	Min::String  qs;

	public:
	Sequence(){ id = 0; }
#if 0
	Sequence( const Sequence & other ){
		id = other.id;
		des = other.des;
		seq = other.seq;
		qs = other.qs;
	}
#endif

	void Copy( const Sequence & other ){
		des = other.des;
		seq = other.seq;
		qs = other.qs;
	}

	int Length()const{ return seq.Size(); }

	Min::String& Description(){ return des; }
	Min::String& SequenceData(){ return seq; }
	Min::String& QualityScore(){ return qs; }
	Min::String GetDescription( int deslen = 0 );

	void Reset(){ des.Reset(); seq.Reset(); qs.Reset(); }
	void ToReverseComplimentary();

	void Print( FILE *fout = stdout, int width = 0 );
};

class SequenceCache
{
	Min::Array<Sequence*>  sequences;

	FILE  *inputFile;

	bool  getdes;
	
	Sequence  buffer;
	
	Min::String  *current;

	public:
	SequenceCache( const Min::String & file ){
		inputFile = fopen( file.Data(), "r" );
		current = NULL;
		getdes = true;
	}
	~SequenceCache(){
		if( inputFile ) fclose( inputFile );
		Clear();
	}
	void Clear(){
		for(int i=0,n=sequences.Size(); i<n; i++) delete sequences[i];
		sequences.Clear();
	}

	int Update( int max = 10000, bool idonly = false );

	size_t Count()const{ return sequences.Size(); }
	Sequence* operator[]( int i )const{ return sequences[i]; }
};

class SequenceList
{
	Min::Array<Sequence*>  sequences;
	
	int  max_length;
	int  min_length;

	public:
	SequenceList(){ max_length = min_length = 0; }
	~SequenceList(){
		Clear();
	}

	void Clear(){
		for(int i=0,n=sequences.Size(); i<n; i++) delete sequences[i];
		sequences.Clear();
	}

	bool ReadFastAQ( const Min::String & file, bool idonly=false, FILE *fout_log=stdout );
	void SortByLength();

	size_t Count()const{ return sequences.Size(); }
	int MaxLength(){
		int i, N = sequences.Size();
		if( N == 0 ) return 0;
		for(i=1, max_length = sequences[0]->Length(); i<N; i++){
			if( sequences[i]->Length() > max_length ) max_length = sequences[i]->Length();
		}
		return max_length;
	}
	int MinLength(){
		int i, N = sequences.Size();
		if( N == 0 ) return 0;
		for(i=1, min_length = sequences[0]->Length(); i<N; i++){
			if( sequences[i]->Length() < min_length ) min_length = sequences[i]->Length();
		}
		return min_length;
	}
	int RemoveEmptySequences();

	void Shuffle();

	Sequence* operator[]( int i )const{ return sequences[i]; }
};


END_NS_BIO
