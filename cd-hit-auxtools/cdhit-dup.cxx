//=================================================================
// This file is a part of cdhit-dup.
// By Limin Fu (phoolimin@gmail.com, lmfu@ucsd.edu)
//================================================================= 

#include <math.h>
#include <ctype.h>
#include "minString.hxx"
#include "minArray.hxx"
#include "minMap.hxx"
#include "bioSequence.hxx"

using namespace Min;
using namespace Bio;

class SequenceCluster : public Array<Sequence*>
{
	int  id;
	int  abundance;
	int  chiHead;
	int  chiTail;

	public:
	SequenceCluster( Sequence *rep = NULL ){
		id = 0;
		abundance = 0;
		chiHead = chiTail = 0;
		if( rep ) Append( rep );
	}

	String GetDescription( Sequence *seq, int deslen=0 );

	int GetID()const{ return id; }
	void SetID( int i ){ id = i; }

	int GetAbundance()const{ return abundance; }
	void SetAbundance( int ab ){ abundance = ab; }

	int GetChimericParent1()const{ return chiHead; }
	int GetChimericParent2()const{ return chiTail; }
	void SetChimericParent( int head, int tail ){ chiHead = head; chiTail = tail; }

	void Write( FILE *fout = stdout, int id = 0, int deslen=0, const char *des=NULL );

	bool operator<( const SequenceCluster & other ){
		assert( Size() && other.Size() );
		return (*this)[0]->Length() < other[0]->Length();
	}
};
String SequenceCluster::GetDescription( Sequence *seq, int deslen )
{
	String des = seq->Description();
	int i = 0;
	if( des[i] == '>' || des[i] == '@' || des[i] == '+' ) i += 1;
	if( des[i] == ' ' || des[i] == '\t' ) i += 1;
	if( deslen == 0 || deslen > des.Size() ) deslen = des.Size();
	while( i < deslen and ! isspace( des[i] ) ) i += 1;
	des.Resize( i );
	return des;
}
void SequenceCluster::Write( FILE *fout, int id, int deslen, const char *cdes )
{
	//Array<Sequence*> & seqs = *this;
	SequenceCluster & seqs = *this;
	String des = GetDescription( seqs[0], deslen );
	int i = 0, n = Size();

	fprintf( fout, ">Cluster %i%s\n", id, cdes ? cdes : "" );
	fprintf( fout, "%i\t%int, >%s... *\n", 0, seqs[0]->Length(), des.Data() );
	for(i=1; i<n; i++){
		String des = GetDescription( seqs[i], deslen );
		int len = seqs[i]->Length();
		int mm = CountMismatch( seqs[0]->SequenceData(), seqs[i]->SequenceData() );
		fprintf( fout, "%i\t%int, >%s... at 1:%i:1:%i/+/%.2f%%\n",
				i, len, des.Data(), len, len, 100*(len-mm)/(float)len );
	}
}

void WriteClusters( Array<SequenceCluster> & clusters, const String & name = "temp.txt", int deslen = 0 )
{
	String cfile = name + ".clstr";
	String cfile2 = name + "2.clstr";
	FILE *fout1 = fopen( name.Data(), "w" );
	FILE *fout2 = fopen( cfile.Data(), "w" );
	FILE *fout3 = fopen( cfile2.Data(), "w" );
	char cdes[200];
	int i, n = clusters.Size();
	int k1 = 0, k2 = 0;
	for(i=0; i<n; i++){
		SequenceCluster & cluster = clusters[i];
		int head = cluster.GetChimericParent1();
		int tail = cluster.GetChimericParent2();
		if( cluster.Size() == 0 ) continue;
		if( head == tail ){
			cluster[0]->Print( fout1 );
			cluster.SetID( k1 );
			cluster.Write( fout2, k1++, deslen );
		}else{
			head = clusters[head].GetID();
			tail = clusters[tail].GetID();
			sprintf( cdes, " chimeric_parent1=%i,chimeric_parent2=%i", head, tail );
			cluster.Write( fout3, k2++, deslen, cdes );
		}
	}
	fclose( fout1 );
	fclose( fout2 );
	fclose( fout3 );
}

//liwz
//skip fout1 since R1 being modified due to padding
void WriteClusters_clstronly( Array<SequenceCluster> & clusters, const String & name = "temp.txt", int deslen = 0 )
{
	String cfile = name + ".clstr";
	String cfile2 = name + "2.clstr";
//	FILE *fout1 = fopen( name.Data(), "w" );
	FILE *fout2 = fopen( cfile.Data(), "w" );
	FILE *fout3 = fopen( cfile2.Data(), "w" );
	char cdes[200];
	int i, n = clusters.Size();
	int k1 = 0, k2 = 0;
	for(i=0; i<n; i++){
		SequenceCluster & cluster = clusters[i];
		int head = cluster.GetChimericParent1();
		int tail = cluster.GetChimericParent2();
		if( cluster.Size() == 0 ) continue;
		if( head == tail ){
			//cluster[0]->Print( fout1 );
			cluster.SetID( k1 );
			cluster.Write( fout2, k1++, deslen );
		}else{
			head = clusters[head].GetID();
			tail = clusters[tail].GetID();
			sprintf( cdes, " chimeric_parent1=%i,chimeric_parent2=%i", head, tail );
			cluster.Write( fout3, k2++, deslen, cdes );
		}
	}
//	fclose( fout1 );
	fclose( fout2 );
	fclose( fout3 );
}

void WriteClusters_seqonly( Array<SequenceCluster> & clusters, const String & name = "temp.txt", int deslen = 0 )
{
	FILE *fout1 = fopen( name.Data(), "w" );
	int i, n = clusters.Size();
	for(i=0; i<n; i++){
		SequenceCluster & cluster = clusters[i];
		int head = cluster.GetChimericParent1();
		int tail = cluster.GetChimericParent2();
		if( cluster.Size() == 0 ) continue;
		if( head == tail ){
			cluster[0]->Print( fout1 );
		}
	}
	fclose( fout1 );
}


void SortByAbundance( Array<SequenceCluster> & clusters )
{
	int i, max = 0, min = 0;
	int N = clusters.Size();
	if( N <= 1 ) return;
	max = min = clusters[0].GetAbundance();
	for(i=1; i<N; i++){
		int ab = clusters[i].GetAbundance();
		if( ab > max ) max = ab;
		if( ab < min ) min = ab;
	}

	int M = max - min + 1;
	Min::Array<int> count( M, 0 ); // count for each size = max_len - i
	Min::Array<int> accum( M, 0 ); // count for all size > max_len - i
	Min::Array<int> offset( M, 0 ); // offset from accum[i] when filling sorting
	Array<SequenceCluster> sorted( N );

	for (i=0; i<N; i++) count[ max - clusters[i].GetAbundance() ] ++;
	for (i=1; i<M; i++) accum[i] = accum[i-1] + count[i-1];
	for (i=0; i<N; i++){
		int len = max - clusters[i].GetAbundance();
		int id = accum[len] + offset[len];
		//clusters[i].index = id;
		sorted[id] = clusters[i];
		offset[len] ++;
	}
	clusters.Swap( sorted );
}


int HashingDepth( int len, int min )
{
	assert( len >= min );
	return (int)sqrt( (len - min) / 10);
}
int HashingLength( int dep, int min )
{
	return min + 10 * dep * dep;
}
struct ChimericSource
{
	int  index;
	int  head;
	int  tail;

	ChimericSource( int i=0, int h=0, int t=0 ){ index = i; head = h; tail = t; }
};
struct HashHit
{
	int  index;
	int  offset;

	HashHit( int i = 0, int o = 0 ){ index = i; offset = o; }
};
#define MAX_OFFSETS 5
struct HashHit3
{
	int  index;
	int  size;
	int  offsets[MAX_OFFSETS];
	int  counts[MAX_OFFSETS]; // number of word hits with the corresponding offsets;

	HashHit3( int i = 0, int o = 0 ){
		memset( offsets, 0, MAX_OFFSETS*sizeof(int) );
		memset( counts, 0, MAX_OFFSETS*sizeof(int) );
		index = i;
		size = 1;
		offsets[0] = o;
		counts[0] = 1;
	}
	void Update( int o ){
		int imin = 0, min = counts[0];
		for(int i=0; i<size; i++){
			if( o == offsets[i] ){
				counts[i] += 1;
				return;
			}else if( counts[i] < min ){
				min = counts[i];
				imin = i;
			}
		}
		if( size == MAX_OFFSETS ){
			offsets[imin] = o;
			counts[imin] = 1;
		}else{
			offsets[size] = o;
			counts[size] = 1;
			size += 1;
		}
	}
	int BestOffset(){
		int imax = 0, max = counts[0];
		for(int i=0; i<size; i++){
			if( counts[i] > max ){
				max = counts[i];
				imax = i;
			}
		}
		//printf( "BestOffset: %i %i %i %i\n", max, imax, size, offsets[imax] );
		return offsets[imax];
	}
};
typedef Hash<unsigned int,Array<int> >  HashHitTable;

void UpdateHashTables( Array<Array<HashHitTable> > & middles, String & seq, int id, int shared, int min, int primer )
{
	unsigned int hash = MakeHash( seq, shared );
	int j, k;
	for(j=primer; (j+shared)<=seq.Size(); j+=shared){
		int dep = j ? 0 : HashingDepth( seq.Size()-j, min );
		for(k=0; k<=dep; k++){
			hash = MakeHash( seq, j, HashingLength( k, min ) );
			middles[j][k][hash].Append( id );
		}
	}
}

void ClusterDuplicate( SequenceList & seqlist, Array<SequenceCluster> & clusters, bool mlen, int errors = 0, float errors2 = 0.0 )
{
	int primer = 0;
	int maxmm = errors;
	unsigned int hash;
	int i, j, k, K, K2, M, N = seqlist.Count();
	int shared, min, count = 0;
	int max = seqlist.MaxLength();
	Hash<unsigned int,Array<int> > wholes;
	Array<Array<Hash<unsigned int,Array<int> > > > middles;
	Node<unsigned int,Array<int> > *node;

	String & first = seqlist[0]->SequenceData();
	primer = first.Size();
	for(i=1; i<N; i++){
		primer = ComparePrefix( first, 0, seqlist[i]->SequenceData(), 0, primer, 0 );
		if( primer < 15 ){
			primer = 0;
			break;
		}
	}
	if( errors == 0 ) maxmm = max * errors2;
	printf( "primer = %i\n", primer );
	shared = (seqlist.MinLength() - primer) / (maxmm + 1);
	if( shared < 30 ) shared = 30;
	min = shared;

	middles.Resize( max );
	for(i=0; (i+shared)<=max; i++) middles[i].Resize( HashingDepth( max - i, min ) + 1 );
	for(i=0; i<N; i++){
		String & seq = seqlist[i]->SequenceData();
		int qlen = seq.Size();
		bool clustered = false;

		maxmm = errors;
		if( errors == 0 ) maxmm = qlen * errors2;
		if( i && (i%10000 == 0) )
			printf( "Clustered %9i sequences with %9i clusters ...\n", i, (int)clusters.Size() );

		hash = MakeHash( seq );
		node = wholes.Find( hash );
		if( node != NULL ){
			Array<int> & hits = node->value;
			for(j=0, M=hits.Size(); j<M; j++){
				SequenceCluster & clust = clusters[hits[j]];
				String & rep = clust[0]->SequenceData();
				if( qlen == rep.Size() && Compare( seq, rep ) == 0 ){
					clust.Append( seqlist[i] );
					clustered = true;
					break;
				}
			}
			if( clustered ) continue;
			for(j=0, M=hits.Size(); j<M; j++){
				SequenceCluster & clust = clusters[hits[j]];
				String & rep = clust[0]->SequenceData();
				if( mlen && rep.Size() != qlen ) continue;
				K = ComparePrefix( seq, rep, maxmm );
				if( K == qlen ){
					clust.Append( seqlist[i] );
					clustered = true;
					break;
				}
			}
			if( clustered ) continue;
		}

		int dep = HashingDepth( qlen - primer, min );
		int hashlen = HashingLength( dep, min );

		hash = MakeHash( seq, primer, hashlen );
		node = middles[primer][ dep ].Find( hash );
		if( node ){
			Array<int> & hits = node->value;
			for(j=0, M=hits.Size(); j<M; j++){
				int hit = hits[j];

				SequenceCluster & clust = clusters[hit];
				String & rep = clust[0]->SequenceData();

				if( mlen && rep.Size() != qlen ) continue;
				K = ComparePrefix( seq, rep, maxmm );
				if( K == qlen ){
					clustered = true;
					clust.Append( seqlist[i] );
					break;
				}
			}
		}

		//seqlists[i]->Print();

		int start = primer;
		dep = 0;//HashingDepth( shared, min );
		while( maxmm == 0 && clustered == false && (start+2*shared) <= qlen && start <= primer + (maxmm+1) * shared ){
			hash = MakeHash( seq, start, shared );
			node = middles[start][ dep ].Find( hash );
			start += shared;
			if( node == NULL ) continue;

			Array<int> & hits = node->value;
			for(j=0, M=hits.Size(); j<M; j++){
				int hit = hits[j];
				SequenceCluster & clust = clusters[hit];
				String & rep = clust[0]->SequenceData();

				if( mlen && rep.Size() != qlen ) continue;
				K = ComparePrefix( seq, rep, maxmm );
				if( K == qlen ){
					clustered = true;
					clust.Append( seqlist[i] );
					break;
				}
			}
		}
		if( clustered == false ){
			hash = MakeHash( seq );
			UpdateHashTables( middles, seq, clusters.Size(), shared, min, primer );
			wholes[ hash ].Append( clusters.Size() );
			clusters.Append( SequenceCluster( seqlist[i] ) );
		}
	}
}

void PrintChimeric( FILE *fout, Sequence *SQ, Sequence *SA, Sequence *SB, int IX, int GAP = 0 )
{
	String A( SA->SequenceData() ), B( SB->SequenceData() ), Q( SQ->SequenceData() );

	fprintf( fout, "\nQuery    = %s\n", SQ->Description().Data() );
	fprintf( fout, "Parent A = %s\n", SA->Description().Data() );
	fprintf( fout, "Parent B = %s\n", SB->Description().Data() );
	fprintf( fout, "Crossover = %i\n", IX );
	//if( offset != 0 ) fprintf( fout, "Offset = %i\n", offset );
	//fprintf( fout, "LQA = %i\n", LQA );
	//fprintf( fout, "LQB = %i\n", LQB );
	//fprintf( fout, "RQA = %i\n", RQA );
	//fprintf( fout, "RQB = %i\n", RQB );

	if( GAP > 0 ){
		while( GAP-- > 0 ) B.Insert( '-', IX );
	}else if( GAP < 0 ){
		while( GAP++ < 0 ){
			A.Insert( '-', IX );
			Q.Insert( '-', IX );
		}
	}
	A.Insert( "|X|", IX );
	Q.Insert( "|X|", IX );
	B.Insert( "|X|", IX );
	int I, J, W = 80;
	int LA = A.Size();
	int LB = B.Size();
	int LQ = Q.Size();
	int L = LA;
	if( L < LB ) L = LB;
	if( L < LQ ) L = LQ;
	for(I=0; I<L; I+=W){
		int M = (I + W) < L ? (I + W) : L;
		fprintf( fout, "A:  " );
		for(J=I; J<M; J++) fprintf( fout, "%c", J < LA ? A[J] : '-' );
		fprintf( fout, "\n" );
		fprintf( fout, " :  " );
		for(J=I; J<M; J++){
			bool match = J < LQ && J < LA && Q[J] == A[J] && Q[J] != '-';
			fprintf( fout, "%c", match ? (J==(IX+1) ? 'X' : '|') : ' ' );
		}
		fprintf( fout, "\n" );
		fprintf( fout, "Q:  " );
		for(J=I; J<M; J++) fprintf( fout, "%c", J < LQ ? Q[J] : '-' );
		fprintf( fout, "\n" );
		fprintf( fout, " :  " );
		for(J=I; J<M; J++){
			bool match = J < LQ && J < LB && Q[J] == B[J] && Q[J] != '-';
			fprintf( fout, "%c", match ? (J==(IX+1) ? 'X' : '|') : ' ' );
		}
		fprintf( fout, "\n" );
		fprintf( fout, "B:  " );
		for(J=I; J<M; J++) fprintf( fout, "%c", J < LB ? B[J] : '-' );
		fprintf( fout, "\n\n" );
	}
}
void PrintChimeric2( FILE *fout, Sequence *SQ, Sequence *SA, Sequence *SB, int IX, int offset )
{
	String A( SA->SequenceData() ), B( SB->SequenceData() ), Q( SQ->SequenceData() );

	fprintf( fout, "\nQuery    = %s\n", SQ->Description().Data() );
	fprintf( fout, "Parent A = %s\n", SA->Description().Data() );
	fprintf( fout, "Parent B = %s\n", SB->Description().Data() );
	fprintf( fout, "Crossover = %i\n", IX );
	//if( offset != 0 ) fprintf( fout, "Offset = %i\n", offset );
	//fprintf( fout, "LQA = %i\n", LQA );
	//fprintf( fout, "LQB = %i\n", LQB );
	//fprintf( fout, "RQA = %i\n", RQA );
	//fprintf( fout, "RQB = %i\n", RQB );

	if( offset > 0 ){
		while( offset-- > 0 ){
			A.Insert( '-', 0 );
			Q.Insert( '-', 0 );
		}
	}else if( offset < 0 ){
		while( offset++ < 0 ) B.Insert( '-', 0 );
	}
	A.Insert( "|X|", IX );
	Q.Insert( "|X|", IX );
	B.Insert( "|X|", IX );
	int I, J, W = 80;
	int LA = A.Size();
	int LB = B.Size();
	int LQ = Q.Size();
	int L = LA;
	if( L < LB ) L = LB;
	if( L < LQ ) L = LQ;
	for(I=0; I<L; I+=W){
		int M = (I + W) < L ? (I + W) : L;
		fprintf( fout, "A:  " );
		for(J=I; J<M; J++) fprintf( fout, "%c", J < LA ? A[J] : '-' );
		fprintf( fout, "\n" );
		fprintf( fout, " :  " );
		for(J=I; J<M; J++){
			bool match = J < LQ && J < LA && Q[J] == A[J] && Q[J] != '-';
			fprintf( fout, "%c", match ? (J==(IX+1) ? 'X' : '|') : ' ' );
		}
		fprintf( fout, "\n" );
		fprintf( fout, "Q:  " );
		for(J=I; J<M; J++) fprintf( fout, "%c", J < LQ ? Q[J] : '-' );
		fprintf( fout, "\n" );
		fprintf( fout, " :  " );
		for(J=I; J<M; J++){
			bool match = J < LQ && J < LB && Q[J] == B[J] && Q[J] != '-';
			fprintf( fout, "%c", match ? (J==(IX+1) ? 'X' : '|') : ' ' );
		}
		fprintf( fout, "\n" );
		fprintf( fout, "B:  " );
		for(J=I; J<M; J++) fprintf( fout, "%c", J < LB ? B[J] : '-' );
		fprintf( fout, "\n\n" );
	}
}

struct ParentInfo
{
	int index;
	int start;
	int len;
	int offset;
	ParentInfo( int index0=0, int start0=0, int len0=0, int offset0=0 ){
		index = index0;
		start = start0;
		len = len0;
		offset = offset0;
	}
	bool operator<( const ParentInfo & other )const{
		return (start - offset) < (other.start - other.offset);
	}
};
struct OffsetCount
{
	int offset;
	int count;

	OffsetCount( int offset = 0, int count = 0 ){
		this->offset = offset;
		this->count = count;
	}
	bool operator<( const OffsetCount & other )const{
		if( count != other.count ) return count > other.count;
		return offset < other.offset;
	}
};
void DetectChimeric( Array<SequenceCluster> & clusters, Array<ChimericSource> & chistat, int max, int shared, float percent, float abratio = 1, int minabu = 2 )
{
	int primer = 20;
	unsigned int hash;
	int i, j, k, K, K2, M, N = clusters.Size();
	int min, count = 0, maxmm = 2, maxmm2 = maxmm + maxmm;
	Hash<unsigned int,Array<HashHit> > hashTable;
	Node<unsigned int,Array<HashHit> > *node;
	Hash<int,int> chimap;
	Array<int> hitMapping;
	Array<unsigned int> hashes;
	Array<HashHit3> hitList;
	Array<ParentInfo> parentAList;
	Array<ParentInfo> parentBList;
	Array<OffsetCount> offsetCounts;
	FILE *debug = NULL;

	if( N <= 2 ) return;

	if( shared < 20 ) shared = 20;
	min = shared;

	//debug = fopen( "debug.txt", "w+" );

	hitList.Resize( N );
	hitMapping.Resize( N );

	String & first = clusters[0][0]->SequenceData();
	primer = first.Size();
	for(i=1; i<N; i++){
		primer = ComparePrefix( first, 0, clusters[i][0]->SequenceData(), 0, primer, 0 );
		if( primer < 15 ){
			primer = 0;
			break;
		}
	}
	printf( "primer = %i\n", primer );

	printf( "Searching for chimeric clusters ...\n" );
	for(i=0; i<N; i++){
		String & seq = clusters[i][0]->SequenceData();
		int qlen = seq.Size();
		int qn = clusters[i].GetAbundance();
		int maxmm = (int)(0.01*qlen*percent);
		int maxmm2 = (int)(0.015*qlen*percent);

		if( i && i % 1000 ==0 ) printf( "Checked %9i clusters, detected %9i chimeric clusters\n", i, (int)chistat.Size() );
		if( qn < minabu ) break;

		hashes.ResetSize( primer );
		for(j=primer; j<=qlen-shared; j++) hashes.Append( MakeHash( seq, j, shared ) );

		if( qlen < (2*shared + primer) ){
			for(j=primer,k=hashes.Size(); j<k; j+=1) hashTable[hashes[j]].Append( HashHit( i, j ) );
			continue;
		}
		//clusters[i][0]->Print();

		HashHit3 *hit = & hitList[0];
		for(j=0,M=hitList.Size(); j<M; j++, hit++) hitMapping[hit->index] = 0;
		hitList.ResetSize( 0 );

		bool detected = false;
		bool duplicate = false;
		int start = primer;
		int dep = 0;

		int proto = -1;
		int protomm = max;
		while( (start + shared) <= qlen ){
			hash = hashes[start];
			node = hashTable.Find( hash );
			start += 1;
			if( node == NULL ) continue;

			Array<HashHit> & hits = node->value;
			for(j=0, M=hits.Size(); j<M; j++){
				int hit = hits[j].index;
				int offset = hits[j].offset - start + 1; // current start is "start - 1";
				if( hitMapping[hit] == 0 ){
					hitList.Append( HashHit3( hit, offset ) );
					hitMapping[hit] = hitList.Size();
				}else{
					hitList[ hitMapping[hit] - 1 ].Update( offset );
				}
			}
		}
		parentAList.ResetSize( 0 );
		parentBList.ResetSize( 0 );
		Array<HashHit3> & hits = hitList;
		for(j=0, M=hits.Size(); j<M; j++){
			int hit = hits[j].index;
			SequenceCluster & clust = clusters[hit];
			String & rep = clust[0]->SequenceData();

			if( clusters[hit].GetAbundance() < qn*abratio ) continue;

			offsetCounts.ResetSize( 0 );
			for(k=0; k<hits[j].size; k++)
				offsetCounts.Append( OffsetCount( hits[j].offsets[k], hits[j].counts[k] ) );

			offsetCounts.QuickSort( 2 );
			//int best = hits[j].BestOffset();

			int tail = 0, head = ComparePrefix( seq, rep, maxmm );
			if( head >= (primer + shared) ) parentAList.Append( ParentInfo( hit, 0, head, 0 ) );

			int best = offsetCounts[0].offset;
			for(k=0; k<2; k++){ // check the best two offsets:
				best = offsetCounts[k].offset;
				if( qlen + best < rep.Size() ){
					tail = CompareSuffix( seq, qlen-1, rep, qlen+best-1, 0, maxmm );
					if( tail >= shared ){
						parentBList.Append( ParentInfo( hit, qlen - tail + best, tail, best ) );
						break;
					}
				}else{
					tail = CompareSuffix( seq, rep.Size()-best-1, rep, rep.Size()-1, 0, maxmm );
					if( tail >= shared ){
						parentBList.Append( ParentInfo( hit, rep.Size() - tail, tail, best ) );
						break;
					}
				}
				tail = 0;
			}

			if( (head + tail) >= qlen ){
				int X = (head + qlen - tail) / 2;
				int mm = CountMismatch2( seq, rep, 0, X, maxmm2 );
				mm += CountMismatch3( seq, X, rep, X+best, qlen, maxmm2 );
				//printf( "\noffset = %5i;  mm = %5i\n", best, mm );
				//PrintChimeric2( stdout, clusters[i][0], clusters[hit][0], clusters[hit][0], X, best );
				if( mm <= maxmm2 /*&& rep.Size() >= qlen*/ ){
					if( mm < protomm ){
						protomm = mm;
						proto = hit;
					}
				}
			}else{
				int mm = CountMismatch2( seq, rep, 0, qlen, maxmm2 );
				if( mm <= maxmm2 /*&& rep.Size() >= qlen*/ ){
					if( mm < protomm ){
						protomm = mm;
						proto = hit;
					}
				}
			}
		}
		if( protomm <= maxmm ){
			if( chimap.Find( proto ) ){
				ChimericSource cs = chistat[chimap[proto]];
				chimap[i] = chistat.Size();
				chistat.Append( ChimericSource( i, cs.head, cs.tail ) );
			}
			continue;
		}
		//printf( "debug: %9i %9i %9i\n", (int)hitList.Size(), (int)parentAList.Size(), (int)parentBList.Size() );
		parentBList.QuickSort();
		for(j=0, M=parentAList.Size(); j<M; j++){
			ParentInfo & infoA = parentAList[j];
			String & repA = clusters[infoA.index][0]->SequenceData();
			for(k=0; k<parentBList.Size(); k++){
				ParentInfo & infoB = parentBList[k];
				String & repB = clusters[infoB.index][0]->SequenceData();
				if( infoA.index == infoB.index ) continue;
				if( (infoB.start - infoB.offset) > infoA.len ) break;
				int XA = (infoA.len + infoB.start - infoB.offset) / 2;
				int XB = XA + infoB.offset;
				int LQA = CountMismatch( seq, repA, 0, XA );
				int RQA = CountMismatch( seq, repA, XA, qlen );
				int LQB = CountMismatch3( seq, 0, repB, 0, XA, qlen );
				int RQB = CountMismatch3( seq, XA, repB, XB, qlen - XA, qlen );
				if( proto < 0 ) protomm = maxmm2;
				if( LQB > (protomm + maxmm) && RQA > (protomm + maxmm) && (LQA + RQB) <= protomm )
				{
					//PrintChimeric2( stdout, clusters[i][0], clusters[infoA.index][0], clusters[infoB.index][0], XA, infoB.offset );
					detected = true;
					chimap[i] = chistat.Size();
					chistat.Append( ChimericSource( i, infoA.index, infoB.index ) );
					count += 1;
					break;
				}
			}
			if( detected ) break;
		}
		

		//if( (i+1) % 1000 ==0 ) printf( "Checked %9i clusters, detected %9i chimeric clusters\n", i+1, chistat.Size() );
		//if( detected == false && duplicate == false ){
		if( detected == false ){
			for(j=primer; j<hashes.Size(); j+=1) hashTable[hashes[j]].Append( HashHit( i, j ) );
		}
	}
	if( debug ) fclose( debug );
}


const char *help =
"CD-HIT-DUP\n\n"
"Usage:\n"
"cd-hit-dup -i input.fa -o output [other options] (for single reads FASTQ)\n"
"cd-hit-dup -i input.fq -o output [other options] (for single reads FASTA)\n"
"cd-hit-dup -i R1.fq -i2 R2.fq -o output -o2 output-R2 [other options] (for PE reads FASTQ)\n"
"cd-hit-dup -i R1.fa -i2 R2.fa -o output -o2 output-R2 [other options] (for PE reads FASTA)\n\n"
"Options:\n"
"    -i        Input file (FASTQ or FASTA);\n"
"    -i2       Second input file (FASTQ or FASTA);\n"
"    -o        Output file;\n"
"    -o2       Output file for R2;\n"
"    -d        Description length (default 0, truncate at the first whitespace character)\n"
"    -u        Length of prefix to be used in the analysis (default 0, for full/maximum length);\n"
"    -m        Match length (true/false, default true);\n"
"    -e        Maximum number of mismatches allowd;\n"
"    -f        Filter out chimeric clusters (true/false, default false);\n"
"    -s        Minimum length of common sequence shared between a chimeric read\n"
"              and each of its parents (default 30, minimum 20);\n"
"    -a        Abundance cutoff (default 1 without chimeric filtering, 2 with chimeric filtering);\n"
"    -b        Abundance ratio between a parent read and chimeric read (default 1);\n"
"    -p        Dissimilarity control for chimeric filtering (default 1);\n"
;

int main( int argc, char *argv[] )
{
	if( argc < 5 ){
		printf( "%s\n", help );
		return 1;
	}
	String input, input2, output, output2;
	bool matchLength = true;
	bool nochimeric = false;
	int abundance = -1;
	int deslen = 0;
	int shared = 30;
	int uselen = 0;
	int errors = 0;
	float errors2 = 0;
	float abratio = 1;
	float percent = 1;
	int i, k, m, N;

	for(i=1; i<argc; i+=2){
		if( i+1 == argc ){
			printf( "Incomplete argument %s\n", argv[i] );
			printf( "\n%s\n", help );
			return 1;
		}
		if( strcmp( argv[i], "-i" ) == 0 ) input = argv[i+1];
		else if( strcmp( argv[i], "-i2" ) == 0 ) input2 = argv[i+1];
		else if( strcmp( argv[i], "-o" ) == 0 ) output = argv[i+1];
		else if( strcmp( argv[i], "-o2" ) == 0 ) output2 = argv[i+1];
		else if( strcmp( argv[i], "-a" ) == 0 ) abundance = strtol( argv[i+1], NULL, 10 );
		else if( strcmp( argv[i], "-d" ) == 0 ) deslen = strtol( argv[i+1], NULL, 10 );
		else if( strcmp( argv[i], "-u" ) == 0 ) uselen = strtol( argv[i+1], NULL, 10 );
		else if( strcmp( argv[i], "-b" ) == 0 ) abratio = strtod( argv[i+1], NULL );
		else if( strcmp( argv[i], "-p" ) == 0 ) percent = strtod( argv[i+1], NULL );
		else if( strcmp( argv[i], "-m" ) == 0 ) matchLength = strcmp( argv[i+1], "true" ) ==0;
		else if( strcmp( argv[i], "-f" ) == 0 ) nochimeric = strcmp( argv[i+1], "true" ) ==0;
		else if( strcmp( argv[i], "-e" ) == 0 ){
			errors = strtol( argv[i+1], NULL, 10 );
			errors2 = strtod( argv[i+1], NULL );
		}else if( strcmp( argv[i], "-s" ) == 0 ){
			shared = strtol( argv[i+1], NULL, 10 );
			nochimeric = true;
			printf( "Chimeric cluster filtering is automatically enabled by \"-s\" parameter!\n" );
		}else{
			printf( "Unknown argument %s\n", argv[i] );
			printf( "\n%s\n", help );
			return 1;
		}
	}
	if( nochimeric && input2.Size() ){
		printf( "ERROR: Chimeric filtering can only be done with single (-i) input file!\n" );
		return 1;
	}
        if (input2.Size() && ( ! output2.Size())) {
                output2 = output + ".2";
        }

	SequenceList seqlist;
	SequenceList seqlist2;
        int seqlist_modified = 0;
        int seqlist2_modified = 0;
	if( not seqlist.ReadFastAQ( input ) ){
		printf( "File openning failed: %s\n", input.Data() );
		return 1;
	}
	if( input2.Size() and not seqlist2.ReadFastAQ( input2 ) ){
		printf( "File openning failed: %s\n", input2.Data() );
		return 1;
	}

	printf( "From input: %s\n", input.Data() );
	printf( "Total number of sequences: %i\n", seqlist.Count() );
	printf( "Longest: %i\n", seqlist.MaxLength() );
	printf( "Shortest: %i\n", seqlist.MinLength() );

	if( input2.Size() ){
		if( seqlist2.Count() != seqlist.Count() ){
			printf( "Error: the pair end files contain different number of reads!\n" );
			return 1;
		}
		printf( "\nFrom input: %s\n", input2.Data() );
		printf( "Total number of sequences: %i\n", seqlist2.Count() );
		printf( "Longest: %i\n", seqlist2.MaxLength() );
		printf( "Shortest: %i\n", seqlist2.MinLength() );

		int max1 = seqlist.MaxLength();
		int max2 = seqlist2.MaxLength();
		int max = max1 > max2 ? max1 : max2;

		if( uselen <= 0 ) uselen = max;
		if( uselen > max ) uselen = max;

		N = seqlist.Count();
		for(k=0; k<N; k++){
			Sequence *first = seqlist[k];
			Sequence *second = seqlist2[k];

			int qs = first->QualityScore().Size();
			String & s1 = first->SequenceData();
			String & s2 = second->SequenceData();
			String & qs1 = first->QualityScore();
			String & qs2 = second->QualityScore();
			int n1 = s1.Size();
			int n2 = s2.Size();
			int i, j, m = 1;

			if( n1 < uselen ){
				s1.Resize( uselen, 'N' );
				qs1.Resize( uselen, 0 );
			}
			if( n2 < uselen ){
				s2.Resize( uselen, 'N' );
				qs2.Resize( uselen, 0 );
                                seqlist2_modified = 1;
			}

			seqlist[k]->SequenceData().Chop( n1 - uselen );
			seqlist[k]->SequenceData() += s2.Data();
			seqlist[k]->SequenceData().Chop( n2 - uselen );
			if( qs ){
				seqlist[k]->QualityScore().Chop( n1 - uselen );
				seqlist[k]->QualityScore() += qs2.Data();
				seqlist[k]->QualityScore().Chop( n2 - uselen );
			}
		}
                seqlist_modified = 1;
	}else if( uselen > 0 ){
		N = seqlist.Count();
		for(k=0; k<N; k++){
			Sequence *first = seqlist[k];
			int qs = first->QualityScore().Size();
			String & s1 = first->SequenceData();
			String & qs1 = first->QualityScore();
			int n1 = s1.Size();
			int i, j, m = 1;

			if( n1 < uselen ) continue;

                	seqlist_modified = 1;

			seqlist[k]->SequenceData().Chop( n1 - uselen );
			if( qs ) seqlist[k]->QualityScore().Chop( n1 - uselen );
		}
	}
	if( seqlist.Count() == 0 ){
		printf( "Abort: no sequence available for the analysis!\n" );
		return 1;
	}

	if( seqlist.MaxLength() != seqlist.MinLength() ){
		seqlist.SortByLength();
		printf( "Sorted by length ...\n" );
	}

	Array<SequenceCluster> clusters;
	Array<ChimericSource> chistat;
	printf( "Start clustering duplicated sequences ...\n" );
	ClusterDuplicate( seqlist, clusters, matchLength, errors, errors2 );
	printf( "Number of reads: %i\n", seqlist.Count() );
	printf( "Number of clusters found: %i\n", clusters.Size() );

	for(i=0, m=clusters.Size(); i<m; i++){
		clusters[i].SetAbundance( clusters[i].Size() );
	}
	if( nochimeric ){
		if( abundance < 0 ) abundance = 2;
		if( seqlist.Count() == clusters.Size() ){
			for(i=0, m=clusters.Size(); i<m; i++){
				Sequence *rep = clusters[i][0];
				String & des = rep->Description();
				int ab = 1;
				int pos = des.Find( "_abundance_" );
				if( pos ) ab = strtol( des.Data() + pos + 10, NULL, 10 );
				clusters[i].SetAbundance( ab );
			}
		}

		SortByAbundance( clusters );
		DetectChimeric( clusters, chistat, seqlist.MaxLength(), shared, percent, abratio, abundance );
		for(i=0, m=0; i<chistat.Size(); i++){
			clusters[chistat[i].index].SetChimericParent( chistat[i].head, chistat[i].tail );
		}
		printf( "Number of chimeric clusters found: %i\n", chistat.Size() );
	}
	if( abundance < 0 ) abundance = 1;
	int above = 0;
	int below = 0;
	for(i=0, m=clusters.Size(); i<m; i++){
		if( clusters[i].GetChimericParent1() == clusters[i].GetChimericParent2() ){
			bool ab = clusters[i].GetAbundance() < abundance;
			above += ab == false;
			below += ab == true;
		}
	}
	printf( "Number of clusters with abundance above the cutoff (=%i): %i\n", abundance, above );
	printf( "Number of clusters with abundance below the cutoff (=%i): %i\n", abundance, below );
	printf( "Writing clusters to files ...\n" );
	// Write .clstr file and merged output, the later will be overwritten if R1 is modified
	// need to write .clstr file now so that the .clstr has correct info, e.g. identity %

	// by liwz - avoid modified R1 being write to file. padded R1 may crash in Print(fout1)
	// WriteClusters( clusters, output, deslen );
        seqlist_modified ? WriteClusters_clstronly( clusters, output, deslen ) : WriteClusters( clusters, output, deslen );
        // liwz:
        // 1. for PE reads, previous code connect R1 and R2 and output R1 and R2 together
        //    now we need to output R1 and R2 in two different files
        // 2. with -u option, the original reads are cut, previous code output cut reads
        //    now we need to output full-length reads
        // since seqlist2 is already read in, 
        // we output R2 first
        if( input2.Size() ){
		//re-read seqlist2 if modified due to padding
		if (seqlist2_modified) {
                	seqlist2.Clear();
			if( input2.Size() and not seqlist2.ReadFastAQ( input2 ) ){
				printf( "File openning failed: %s\n", input2.Data() );
				return 1;
			}
                }
		
		// copy seqlist2 => seqlist
		N = seqlist.Count();
		for(k=0; k<N; k++){
			int qs = seqlist[k]->QualityScore().Size();
			seqlist[k]->SequenceData().Reset();
                        seqlist[k]->SequenceData() += seqlist2[k]->SequenceData();
			seqlist[k]->Description().Reset();
                        seqlist[k]->Description() += seqlist2[k]->Description();
			if( qs ){
				seqlist[k]->QualityScore().Reset();
				seqlist[k]->QualityScore() += seqlist2[k]->QualityScore();
			}
			 seqlist_modified = 1;
		}

           // write R2
	   printf( "Write R2 reads\n" );
           WriteClusters_seqonly( clusters, output2, deslen );
        } 
	if (seqlist_modified){
		// since seqlist was nodified, re-read from file, this time use seqlist2 as storage
		seqlist2.Clear();
		if( not seqlist2.ReadFastAQ( input ) ){
			printf( "File openning failed: %s\n", input.Data() );
			return 1;
		}
		// copy seqlist2 => seqlist
		N = seqlist.Count();
		for(k=0; k<N; k++){
			int qs = seqlist[k]->QualityScore().Size();
			seqlist[k]->SequenceData().Reset();
                        seqlist[k]->SequenceData() += seqlist2[k]->SequenceData();
			seqlist[k]->Description().Reset();
                        seqlist[k]->Description() += seqlist2[k]->Description();
			if( qs ){
				seqlist[k]->QualityScore().Reset();
				seqlist[k]->QualityScore() += seqlist2[k]->QualityScore();
			}
		}
		// overwrite R1 if R1 was modified
	   	printf( "Write R1 reads\n" );
        	WriteClusters_seqonly( clusters, output, deslen );
	}

	printf( "Done!\n" );
	return 0;
}
