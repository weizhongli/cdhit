//=================================================================
// This file is a part of cdhit-dup.
// By Limin Fu (phoolimin@gmail.com, lmfu@ucsd.edu)
//================================================================= 

#include <ctype.h>
#include "bioSequence.hxx"

using namespace Min;
using namespace Bio;

Min::String Bio::Sequence::GetDescription( int deslen )
{
	String des = Description();
	int i = 0;
	if( des[i] == '>' || des[i] == '@' || des[i] == '+' ) i += 1;
	if( des[i] == ' ' || des[i] == '\t' ) i += 1;
	if( deslen == 0 || deslen > des.Size() ) deslen = des.Size();
	while( i < deslen and ! isspace( des[i] ) ) i += 1;
	des.Resize( i );
	return des;
}
const char *dna_reverse_complimentary = "TBGDEFCHIJKLMNOPQRSAAVWXYZ";
void Bio::Sequence::ToReverseComplimentary()
{
	int i, j;
	int n = seq.Size();
	int m = (n/2) + (n%2);
	for(i=0, j=n-1; i<m; i++, j--){
		char L = seq[i];
		char R = seq[j];
		seq[i] = dna_reverse_complimentary[ R - 'A' ];
		seq[j] = dna_reverse_complimentary[ L - 'A' ];
	}
	if( qs.Size() != n ) return;
	for(i=0, j=n-1; i<m; i++, j--){
		char L = qs[i];
		char R = qs[j];
		qs[i] = R;  qs[j] = L;
	}
}

void Bio::Sequence::Print( FILE *fout, int width )
{
	if( qs.Size() ){
		assert( seq.Size() == qs.Size() );
		if( des.Size() ){
			fprintf( fout, "@%s\n", des.Data() );
		}else{
			fprintf( fout, "@%i\n", id );
		}
		fprintf( fout, "%s\n", seq.Data() );
		if( des.Size() ){
			fprintf( fout, "+%s\n", des.Data() );
		}else{
			fprintf( fout, "+%i\n", id );
		}
		fprintf( fout, "%s\n", qs.Data() );
	}else{
		if( des.Size() ){
			fprintf( fout, ">%s\n", des.Data() );
		}else{
			fprintf( fout, ">%i\n", id );
		}
		fprintf( fout, "%s\n", seq.Data() );
	}
}

#define BUFSIZE 1024

int Bio::SequenceCache::Update( int max, bool idonly )
{
	char buf[ BUFSIZE+1 ];
	char *res = NULL;
	FILE *fin = inputFile;
	String & ds = buffer.Description();
	String & ss = buffer.SequenceData();
	String & qs = buffer.QualityScore();

	Clear();

	if( fin == NULL ) return 0;
	while( feof( fin ) == 0 ){
		if( (res=fgets( buf, BUFSIZE, fin )) == NULL ) break;
		if( getdes && (buf[0] == '>' || buf[0] == '@' || buf[0] == '+') ){
			getdes = 0;
			if( buf[0] == '+' ){
				current = & qs;
			}else{
				if( ss.Size() ){
					ss.ToUpper();
					sequences.Append( new Sequence( buffer ) );
					buffer.Reset();
					if( sequences.Size() % 100000 == 0 )
						printf( "Read %9i sequences ...\n", (int)sequences.Size() );
				}

				current = & ss;

				ds = buf + 1;
				while( ds[ds.Size()-1] != '\n' ){
					if( (res=fgets( buf, BUFSIZE, fin )) == NULL ) break;
					ds += buf;
				}
				ds.Chop();
				if( idonly ){
					int d = 0;
					while( d < ds.Size() && isspace( ds[d] ) ==0 ) d += 1;
					ds.Resize( d );
				}
				//printf( "%s\n", ds->Data() );
				if( max > 0 && (int)sequences.Size() >= max ) return sequences.Size();
			}
		}else{
			assert( current != NULL );
			*current += buf;
			while( current->Size() && isspace( (*current)[current->Size()-1] ) ) current->Chop();
			getdes = current == & ss || (current == & qs && qs.Size() == ss.Size());
		}
	}
	if( ss.Size() ){
		sequences.Append( new Sequence( buffer ) );
		buffer.Reset();
	}
	return sequences.Size();
}

bool Bio::SequenceList::ReadFastAQ( const Min::String & file, bool idonly, FILE *fout_log )
{
	char buf[ BUFSIZE+1 ];
	FILE *fin = fopen( file.Data(), "r" );
	bool wasEmpty = sequences.Size() == 0;
	char *res = NULL;
	int max = 0, min = 0x7fffffff;
	int getdes = 1;
	Sequence one;
	String & ds = one.Description();
	String & ss = one.SequenceData();
	String & qs = one.QualityScore();
	String *data=NULL;

	if( fin == NULL ) return 0;
	while( feof( fin ) == 0 ){
		if( (res=fgets( buf, BUFSIZE, fin )) == NULL ) break;
		if( getdes && (buf[0] == '>' || buf[0] == '@' || buf[0] == '+') ){
			getdes = 0;
			if( buf[0] == '+' ){
				data = & qs;
			}else{
				if( ss.Size() ){
					ss.ToUpper();
					if( ss.Size() > max ) max = ss.Size();
					if( ss.Size() < min ) min = ss.Size();
					sequences.Append( new Sequence( one ) );
					one.Reset();
					if( sequences.Size() % 100000 == 0 )
						fprintf( fout_log, "Read %9i sequences ...\n", (int)sequences.Size() );
				}

				data = & ss;

				ds = buf + 1;
				while( ds[ds.Size()-1] != '\n' ){
					if( (res=fgets( buf, BUFSIZE, fin )) == NULL ) break;
					ds += buf;
				}
				ds.Chop();
				if( idonly ){
					int d = 0;
					while( d < ds.Size() && isspace( ds[d] ) ==0 ) d += 1;
					ds.Resize( d );
				}
				//printf( "%s\n", ds->Data() );
			}
		}else{
			*data += buf;
			while( data->Size() && isspace( (*data)[data->Size()-1] ) ) data->Chop();
			getdes = data == & ss || (data == & qs && qs.Size() == ss.Size());
		}
	}
	fclose( fin );
	if( ss.Size() ){
		if( ss.Size() > max ) max = ss.Size();
		if( ss.Size() < min ) min = ss.Size();
		sequences.Append( new Sequence( one ) );
	}
	if( wasEmpty ){
		max_length = max;
		min_length = min;
	}else{
		if( max > max_length ) max_length = max;
		if( min < min_length ) min_length = min;
	}
	return true;
}
void Bio::SequenceList::SortByLength()
{
	if( max_length == min_length ) return;

	int i;
	int N = sequences.Size();
	int M = max_length - min_length + 1;
	Min::Array<int> count( M, 0 ); // count for each size = max_len - i
	Min::Array<int> accum( M, 0 ); // count for all size > max_len - i
	Min::Array<int> offset( M, 0 ); // offset from accum[i] when filling sorting
	Min::Array<Sequence*> sorting( N );

	for (i=0; i<N; i++) count[ max_length - sequences[i]->Length() ] ++;
	for (i=1; i<M; i++) accum[i] = accum[i-1] + count[i-1];
	for (i=0; i<N; i++){
		int len = max_length - sequences[i]->Length();
		int id = accum[len] + offset[len];
		//sequences[i].index = id;
		sorting[id] = sequences[i];
		offset[len] ++;
	}
	sequences.Swap( sorting );
}
int Bio::SequenceList::RemoveEmptySequences()
{
	int i, M=0, N = sequences.Size();
	for(i=0; i<N; i++){
		if( sequences[i]->Length() == 0 ){
			delete sequences[i];
			continue;
		}
		sequences[M++] = sequences[i];
	}
	sequences.Resize( M );
}
void Bio::SequenceList::Shuffle()
{
	int i, N = sequences.Size();
	if( N <= 1 ) return;
	for(i=N-1; i>=0; i--){
		int rnd = (int)(i * (rand() / (1.0 + RAND_MAX)));
		Sequence *seq = sequences[i];
		sequences[i] = sequences[rnd];
		sequences[rnd] = seq;
	}
}
