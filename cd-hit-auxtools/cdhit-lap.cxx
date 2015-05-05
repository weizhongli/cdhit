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

char base_mapping[128] = {0};
char rev_comp_mapping[128] = {0};
uint64_t base_powers[33];

FILE *fout_log = stdout;
FILE *fout_rep = NULL;
FILE *fout_clstr = NULL;

uint64_t EncodeWord( const char *word, unsigned char W )
{
	unsigned char i;
	uint64_t *power = base_powers;
	uint64_t code = 0;
	for(i=0; i<W; i++, word++, power++) code += base_mapping[ (unsigned char)*word ] * *power;
	return code;
}

struct ClusterMember
{
	Sequence *seq;
	short offset;  // offset on query
	short offset2; // offset on representative
	bool revcomp;

	ClusterMember( Sequence *s = NULL, int o = 0, int o2 = 0, bool rc = false ){
		seq = s;
		offset = o;
		offset2 = o2;
		revcomp = rc;
	}
};

class SequenceCluster : public Array<ClusterMember>
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
		if( rep ) Append( ClusterMember( rep ) );
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
		return (*this)[0].seq->Length() < other[0].seq->Length();
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
	String des = GetDescription( seqs[0].seq, deslen );
	String & rep = seqs[0].seq->SequenceData();
	int len2 = seqs[0].seq->Length();
	int i = 0, n = Size();

	fprintf( fout, ">Cluster %i%s\n", id, cdes ? cdes : "" );
	fprintf( fout, "%i\t%int, >%s... *\n", 0, seqs[0].seq->Length(), des.Data() );
	for(i=1; i<n; i++){
		String & seq = seqs[i].seq->SequenceData();
		String des = GetDescription( seqs[i].seq, deslen );
		int len = seqs[i].seq->Length();
		int startQ = seqs[i].offset + 1;
		int startR = seqs[i].offset2 + 1;
		int endQ = len;
		int endR = len2;
		int mm = 0;
		if( (endQ - startQ) < (endR - startR) ){
			endR = startR + endQ - startQ;
		}else{
			endQ = startQ + endR - startR;
		}
		if( seqs[i].revcomp ){
			startQ = len - startQ + 1;
			endQ = len - endQ + 1;
		}
		fprintf( fout, "%i\t%int, >%s... at %i:%i:%i:%i/%c/%.2f%%\n", 
				i, len, des.Data(), startQ, endQ, startR, endR, 
				seqs[i].revcomp ? '-' : '+', 100*(len-mm)/(float)len );
	}
}
void WriteClusters( Array<SequenceCluster> & clusters, const String & name = "temp.txt", int deslen = 0 )
{
	char cdes[200];
	int i, n = clusters.Size();
	int k1 = 0, k2 = 0;
	for(i=0; i<n; i++){
		SequenceCluster & cluster = clusters[i];
		if( cluster.Size() == 0 ) continue;
		cluster[0].seq->Print( fout_rep );
		cluster.SetID( k1 );
		cluster.Write( fout_clstr, k1++, deslen );
	}
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
#define USE_WORD
void ClusterOverlap( SequenceList & seqlist, Array<SequenceCluster> & clusters, int min, float minper )
{
	String reverse;
	Array<SequenceCluster*> clusters2;
#ifdef USE_WORD
	Hash<uint64_t,Array<SequenceCluster*> > headHashes;
	Hash<uint64_t,Array<SequenceCluster*> > tailHashes;
	Node<uint64_t,Array<SequenceCluster*> > *node;
#else
	Hash<unsigned int,Array<SequenceCluster*> > headHashes;
	Hash<unsigned int,Array<SequenceCluster*> > tailHashes;
	Node<unsigned int,Array<SequenceCluster*> > *node;
#endif
	int i, i2, j, k, m, N = seqlist.Count();
	int WORD = min;
#ifndef USE_WORD
	if( WORD > 30 ) WORD = 30;
#endif
	for(i=0; i<N; i++){
		Sequence *seq = seqlist[i];
		int len = seq->Length();
		int min2 = len * minper;
		bool clustered = false;

		if( min2 < min ) min2 = min;

		String *ss = & seq->SequenceData();
		reverse.ResetSize( len );
		for(j=0; j<len; j++) reverse[j] = rev_comp_mapping[ ss->Data()[len - j - 1] ];

		for(j=(len-WORD-1); j>=min2-WORD; j--){
			for(i2=0; i2<2; i2++){
				String *ss = & seq->SequenceData();
				if( i2 ) ss = & reverse;
				//if( i2 ) break;
#ifdef USE_WORD
				uint64_t hash = EncodeWord( ss->Data() + j, WORD );
#else
				unsigned int hash = MakeHash( *ss, j, WORD );
#endif
				node = tailHashes.Find( hash );
				if( node == NULL ) continue;
				Array<SequenceCluster*> & clusts = node->value;
				for(int j2=0, m=clusts.Size(); j2<m; j2++){
					SequenceCluster & clust = *clusts[j2];
					String & rep = clust[0].seq->SequenceData();
					if( rep.Size() < (j + WORD) ) continue;
					const char *r = rep.Data() + (rep.Size() - j - WORD);
					const char *q = ss->Data();
					int M = j + WORD;
					for(k=0; k<M; k++, q++, r++) if( *q != *r ) break;
					if( k == M ){
						clust.Append( ClusterMember( seq, 0, rep.Size() - j - WORD, i2 == 1 ) );
						clustered = true;
						break;
					}
				}
				if( clustered ) break;
			}
			if( clustered ) break;
		}
		if( clustered == false ){
			for(j=0; j<=len-min2; j++){
				for(i2=0; i2<2; i2++){
					String *ss = & seq->SequenceData();
					if( i2 ) ss = & reverse;
					//if( i2 ) break;
#ifdef USE_WORD
					uint64_t hash = EncodeWord( ss->Data() + j, WORD );
#else
					unsigned int hash = MakeHash( *ss, j, WORD );
#endif
					node = headHashes.Find( hash );
					if( node == NULL ) continue;
					Array<SequenceCluster*> & clusts = node->value;
					for(int j2=0, m=clusts.Size(); j2<m; j2++){
						SequenceCluster & clust = *clusts[j2];
						String & rep = clust[0].seq->SequenceData();
						if( rep.Size() < (len - j) ) continue;
						const char *r = rep.Data();
						const char *q = ss->Data() + j;
						int M = len - j;
						for(k=0; k<M; k++, q++, r++) if( *q != *r ) break;
						if( k == M ){
							clust.Append( ClusterMember( seq, j, 0, i2 == 1 ) );
							clustered = true;
							break;
						}
					}
					if( clustered ) break;
				}
				if( clustered ) break;
			}
		}
		if( not clustered ){
#ifdef USE_WORD
			uint64_t head = EncodeWord( seq->SequenceData().Data(), WORD );
			uint64_t tail = EncodeWord( seq->SequenceData().Data() + (len - WORD), WORD );
#else
			unsigned int head = MakeHash( seq->SequenceData(), WORD );
			unsigned int tail = MakeHash( seq->SequenceData(), (len - WORD), WORD );
#endif
			clusters2.Append( new SequenceCluster( seq ) );
			headHashes[head].Append( clusters2.Back() );
			tailHashes[tail].Append( clusters2.Back() );
		}
		if( (i+1)%100000 == 0 )
			fprintf( fout_log, "Clustered %9i sequences with %9i clusters ...\n", i+1, clusters2.Size() );
	}
	for(i=0, m=clusters2.Size(); i<m; i++){
		clusters.Append( SequenceCluster() );
		clusters.Back().Swap( *clusters2[i] );
	}
}

const char *help =
"Options:\n"
"    -i        Input file;\n"
"    -o        Output file;\n"
"    -m        Minimum length of overlapping part (default 20);\n"
"    -p        Minimum percentage of overlapping part (default 0, any percentage);\n"
"    -d        Description length (default 0, truncate at the first whitespace character)\n"
"    -s        Random number seed for shuffling (default 0, no shuffling; shuffled before sorting by length);\n"
"    -stdout   Standard output type (default \"log\", other options \"rep\", \"clstr\");\n"
;

int main( int argc, char *argv[] )
{
	if( argc < 5 ){
		printf( "%s\n", help );
		return 1;
	}
	String input, output, stdout_type = "log";
	unsigned int seed = 0;
	int minlen = 20;
	int deslen = 0;
	float minper = 0;
	int i, m;

	base_mapping[(unsigned)'a'] = 0;
	base_mapping[(unsigned)'A'] = 0;
	base_mapping[(unsigned)'c'] = 1;
	base_mapping[(unsigned)'C'] = 1;
	base_mapping[(unsigned)'g'] = 2;
	base_mapping[(unsigned)'G'] = 2;
	base_mapping[(unsigned)'t'] = 3;
	base_mapping[(unsigned)'T'] = 3;

	memset( rev_comp_mapping, 'A', 128 );
	rev_comp_mapping[(unsigned)'a'] = 't';
	rev_comp_mapping[(unsigned)'A'] = 'T';
	rev_comp_mapping[(unsigned)'c'] = 'g';
	rev_comp_mapping[(unsigned)'C'] = 'G';
	rev_comp_mapping[(unsigned)'g'] = 'c';
	rev_comp_mapping[(unsigned)'G'] = 'C';
	rev_comp_mapping[(unsigned)'t'] = 'a';
	rev_comp_mapping[(unsigned)'T'] = 'A';
	base_powers[0] = 1;
	for(i=1; i<=32; i++) base_powers[i] = base_powers[i-1] * 4;

	for(i=1; i<argc; i+=2){
		if( i+1 == argc ){
			printf( "Incomplete argument %s\n", argv[i] );
			printf( "\n%s\n", help );
			return 1;
		}
		if( strcmp( argv[i], "-i" ) == 0 ) input = argv[i+1];
		else if( strcmp( argv[i], "-o" ) == 0 ) output = argv[i+1];
		else if( strcmp( argv[i], "-stdout" ) == 0 ) stdout_type = argv[i+1];
		else if( strcmp( argv[i], "-m" ) == 0 ) minlen = strtol( argv[i+1], NULL, 10 );
		else if( strcmp( argv[i], "-p" ) == 0 ) minper = strtod( argv[i+1], NULL );
		else if( strcmp( argv[i], "-d" ) == 0 ) deslen = strtol( argv[i+1], NULL, 10 );
		else if( strcmp( argv[i], "-s" ) == 0 ) seed = strtol( argv[i+1], NULL, 10 );
		else{
			printf( "Unknown argument %s\n", argv[i] );
			printf( "\n%s\n", help );
			return 1;
		}
	}
	if( stdout_type == "log" ){
		String cfile = output + ".clstr";
		fout_rep = fopen( output.Data(), "w" );
		fout_clstr = fopen( cfile.Data(), "w" );
	}else if( stdout_type == "rep" ){
		String cfile = output + ".log";
		String cfile2 = output + ".clstr";
		fout_log = fopen( cfile.Data(), "w" );
		fout_clstr = fopen( cfile2.Data(), "w" );
		fout_rep = stdout;
	}else if( stdout_type == "clstr" ){
		String cfile = output + ".log";
		fout_log = fopen( cfile.Data(), "w" );
		fout_rep = fopen( output.Data(), "w" );
		fout_clstr = stdout;
	}else{
		printf( "Unknown standard output type %s\n", stdout_type.Data() );
		return 1;
	}

	SequenceList seqlist;
	if( not seqlist.ReadFastAQ( input, false, fout_log ) ){
		printf( "File openning failed: %s\n", input.Data() );
		return 1;
	}

	fprintf( fout_log, "Total number of sequences: %i\n", seqlist.Count() );
	fprintf( fout_log, "Longest: %i\n", seqlist.MaxLength() );
	fprintf( fout_log, "Shortest: %i\n", seqlist.MinLength() );

	if( seed ){
		srand( seed );
		seqlist.Shuffle();
	}
	if( seqlist.MaxLength() != seqlist.MinLength() ){
		seqlist.SortByLength();
		fprintf( fout_log, "Sorted by length ...\n" );
	}

	Array<SequenceCluster> clusters;
	fprintf( fout_log, "Start clustering duplicated sequences ...\n" );
	ClusterOverlap( seqlist, clusters, minlen, minper );
	fprintf( fout_log, "Number of clusters found: %i\n", clusters.Size() );

	WriteClusters( clusters, output, deslen );
	fprintf( fout_log, "Done!\n" );
	fclose( fout_log );
	fclose( fout_rep );
	fclose( fout_clstr );
	return 0;
}
