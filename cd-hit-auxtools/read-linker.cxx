#include "bioSequence.hxx"

using namespace Min;
using namespace Bio;

class SequenceLinker
{
	int  minOverlap;
	int  tolerance;
	int  maxlen;

	int  overlap;
	int  error;
	int  count;

	Array<Sequence> sequences;

	public:
	SequenceLinker( int minOverlap = 10, int tolerance = 1 ){
		this->minOverlap = minOverlap;
		this->tolerance = tolerance;
		maxlen = 0;
		error = 0;
		count = 0;
		overlap = 0;
	}

	void SetMaxLen( int max ){ maxlen = max; }

	int Overlap()const{ return overlap; }
	int Error()const{ return error; }
	int Count()const{ return count; }

	int Check( const String & first, const String & second );
	int Link( Sequence *first, Sequence *second );
	int Link2( Sequence *first, Sequence *second );
	void Reset(){ count = 0; sequences.Clear(); }
	void Write( FILE *fout ){
		for(int i=0,n=sequences.Size(); i<n; i++) {
                   //printf("-->%d\t%s\n", i, sequences[i].Description().Data());
                   sequences[i].Print( fout );
                }
	}
};


int SequenceLinker::Check( const String & first, const String & second )
{
	int n1 = first.Size();
	int n2 = second.Size();
	int m = n1 < n2 ? n1 : n2;
	int i, j, k;
	for(k=m; k>=minOverlap; k--){
		error = 0;
		for(i=n1-k, j=0; j<k; i++, j++){
			error += (first[i] != second[j]) || (first[i] == 'N' && second[j] == 'N');
			if( error > tolerance ) break;
		}
		overlap = k;
		if( error <= tolerance ) return k;
	}
	overlap = 0;
	return 0;
}
int SequenceLinker::Link( Sequence *first, Sequence *second )
{
	int qs = first->QualityScore().Size();
	const String & s1 = first->SequenceData();
	const String & s2 = second->SequenceData();
	const String & qs1 = first->QualityScore();
	const String & qs2 = second->QualityScore();
	int n1 = s1.Size();
	int n2 = s2.Size();
	int i, j, m = 1, O = Check( s1, s2 );
	int error2 = 0;
	if( O == 0 ) return 0;
	String locs;
	for(i=n1-O, j=0; j<O; i++, j++){
		error2 += s1[i] != s2[j] && s1[i] != 'N' && s2[j] != 'N';
		if( s1[i] != s2[j] ){
			if( locs.Size() ) locs += ',';
			locs += String::FromNumber( i+1 );
		}
	}
	m = 1 << error2;
	while( (count + m) > (int)sequences.Size() ) sequences.Append( Sequence() );
	for(i=0; i<m; i++){
		Sequence & seq = sequences[count+i];
		char tag[50];
		sprintf( tag, ".contig.%i length=%i overlap=%i mismatch_no=%i", (i+1), n1+n2-O, O, error );
		seq.Reset();
		seq.Description() = first->GetDescription();
		seq.Description() += tag;
		if( locs.Size() ) seq.Description() += " mismatch_pos=" + locs;
		seq.SequenceData() += s1;
		seq.SequenceData() += s2.Data() + O;
		if( qs ){
			seq.QualityScore() += qs1;
			seq.QualityScore() += qs2.Data() + O;
		}
	}
        
	int step = m;
	for(i=n1-O, j=0; j<O; i++, j++){
		if( s1[i] == s2[j] ) continue;
		if( s1[i] == 'N' ){
			for(int k=0; k<m; k++){
				Sequence & seq = sequences[count+k];
				seq.SequenceData()[i] = s2[j];
				if( qs ) seq.QualityScore()[i] = qs2[j];
			}
		}else if( s2[j] == 'N' ){
			for(int k=0; k<m; k++){
				Sequence & seq = sequences[count+k];
				seq.SequenceData()[i] = s1[i];
				if( qs ) seq.QualityScore()[i] = qs1[i];
			}
		}else{
			step = step >> 1;
			for(int k=0; k<m; k++){
				Sequence & seq = sequences[count+k];
				int odd = (k/step) % 2;
				seq.SequenceData()[i] = odd ? s2[j] : s1[i];
				if( qs ) seq.QualityScore()[i] = odd ? qs2[j] : qs1[i];
			}
		}
	}
	count += m;
	return m;
}
int SequenceLinker::Link2( Sequence *first, Sequence *second )
{
	int qs = first->QualityScore().Size();
	const String & s1 = first->SequenceData();
	const String & s2 = second->SequenceData();
	const String & qs1 = first->QualityScore();
	const String & qs2 = second->QualityScore();
	int n1 = s1.Size();
	int n2 = s2.Size();
	int i, j, m = 1;

	if( n1 < minOverlap || n2 < minOverlap ) return 0;
	if( n2 < maxlen ) return 0;

	Sequence joint = Sequence();
	joint.Description() = first->GetDescription();
	joint.SequenceData() += s1;
	joint.SequenceData().Chop( n1 - minOverlap );
	joint.SequenceData() += s2.Data() + (maxlen - minOverlap);
	joint.SequenceData().Chop( n2 - maxlen );
	if( qs ){
		joint.QualityScore() = first->QualityScore();
		joint.QualityScore().Chop( n1 - minOverlap );
		joint.QualityScore() += qs2.Data() + (maxlen - minOverlap);
		joint.QualityScore().Chop( n2 - maxlen );
	}
	sequences.Append( joint );
	return 1;
}

const char *help =
"Options:\n"
"    -1 file       Input file, first end;\n"
"    -2 file       Input file, second end;\n"
"    -o file       Output file;\n"
"    -l number     Minimum overlapping length (default 10);\n"
//"                  Or the selection length, if \"number\" is negative,\n"
//"                  to simply join the two ends;\n"
"    -e number     Maximum number of errors (mismatches, default 1);\n"
//"    -m number     Cutoff length when \"-l\" is negative;\n"
;

int main( int argc, char *argv[] )
{
	String first;
	String second;
	String output;
	int min = 10;
	int error = 1;
	int maxlen = 0;
	int i;

	if( argc < 7 ){
		printf( "%s\n", help );
		return 1;
	}
	for(i=1; i<argc; i+=2){
		if( i+1 == argc ){
			printf( "Incomplete argument %s\n", argv[i] );
			printf( "\n%s\n", help );
			return 1;
		}
		if( strcmp( argv[i], "-1" ) == 0 ) first = argv[i+1];
		else if( strcmp( argv[i], "-2" ) == 0 ) second = argv[i+1];
		else if( strcmp( argv[i], "-o" ) == 0 ) output = argv[i+1];
		else if( strcmp( argv[i], "-l" ) == 0 ) min = strtol( argv[i+1], NULL, 10 );
		else if( strcmp( argv[i], "-e" ) == 0 ) error = strtol( argv[i+1], NULL, 10 );
		else if( strcmp( argv[i], "-m" ) == 0 ) maxlen = strtol( argv[i+1], NULL, 10 );
		else{
			printf( "Unknown argument %s\n", argv[i] );
			printf( "\n%s\n", help );
			return 1;
		}
	}

	SequenceLinker linker( abs( min ), error );
	SequenceCache cache1( first );
	SequenceCache cache2( second );
	FILE *fout = fopen( output.Data(), "w" );

	int omin = 0x7fffffff;
	int omax = 0;
	int count = 0;
	int countInputPairs = 0;
	int countUsedPairs = 0;
	int countAllContig = 0;
	Array<int> countPairs( error+1, 0 );
	Array<int> countContigs( error+1, 0 );

	linker.SetMaxLen( maxlen );
	while(1){
                //cache1.Clear();
                //cache2.Clear();
		int n1 = cache1.Update();
		int n2 = cache2.Update();
		int n = n1 < n2 ? n1 : n2;
		if( n1 != n2 ){
			printf( "Warning: the pair end files contain different number of reads!\n" );
		}
		if( n == 0 ) break;
                //printf("==>%d\t%d\n", n1, n2);
		if( min < 0 && maxlen == 0 ){
			maxlen = cache2[0]->Length();
			linker.SetMaxLen( maxlen );
		}
		//Array<Sequence> sequences( n );
		for(int i=0; i<n; i++){
			cache2[i]->ToReverseComplimentary();
			int c = 0;
			if( min > 0 ){
				c = linker.Link( cache1[i], cache2[i] );
				int o = linker.Overlap();
				if( o > omax ) omax = o;
				if( o && o < omin ) omin = o;
			}else{
				c = linker.Link2( cache1[i], cache2[i] );
			}
			countUsedPairs += c != 0;
			countPairs[ linker.Error() ] += c != 0;
			countContigs[ linker.Error() ] += c;
		}
		countInputPairs += n;
		countAllContig += linker.Count();;
		linker.Write( fout );
		linker.Reset();
		count += n;
		printf( "handled: %9i\n", count );
	}
	printf( "Total input pairs of read: %i\n", countInputPairs );
	printf( "Total pairs of read used: %i\n", countUsedPairs );
	printf( "Total contigs: %i\n", countAllContig );
	if( min > 0 ){
		for(i=0; i<=error; i++)
			printf( "%2i mismatch: %9i pairs of reads => %9i contigs\n", i, countPairs[i], countContigs[i] );

		printf( "Overlap range %i-%i\n", omin, omax );
	}
	fclose( fout );
	return 0;
}
