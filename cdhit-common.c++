// =============================================================================
// CD-HI/CD-HIT
//
// Cluster Database at High Identity
//
// CD-HIT clusters protein sequences at high identity threshold.
// This program can remove the high sequence redundance efficiently.
//
// program written by
//                    Weizhong Li
//                    UCSD, San Diego Supercomputer Center
//                    La Jolla, CA, 92093
//                    Email liwz@sdsc.edu
//
//                 at
//                    Adam Godzik's lab
//                    The Burnham Institute
//                    La Jolla, CA, 92037
//                    Email adam@burnham-inst.org
//
// modified by:
//                    Limin Fu
//                    Center for Research in Biological Systems (CRBS), UCSD
//                    La Jolla, CA, 92093
//                    Email: l2fu@ucsd.edu, fu@daovm.net
// =============================================================================

#include "cdhit-common.h"
#include<valarray>
#include<stdint.h>
#include<assert.h>
#include<limits.h>

#ifndef NO_OPENMP

#include<omp.h>

#define WITH_OPENMP " (+OpenMP)"

#else

#define WITH_OPENMP ""

#define omp_set_num_threads(T) (T = T)
#define omp_get_thread_num() 0

#endif

//class function definition
const char aa[] = {"ARNDCQEGHILKMFPSTWYVBZX"};
//{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,2,6,20};
int aa2idx[] = {0, 2, 4, 3, 6, 13,7, 8, 9,20,11,10,12, 2,20,14,
                5, 1,15,16,20,19,17,20,18, 6};
    // idx for  A  B  C  D  E  F  G  H  I  J  K  L  M  N  O  P
    //          Q  R  S  T  U  V  W  X  Y  Z
    // so  aa2idx[ X - 'A'] => idx_of_X, eg aa2idx['A' - 'A'] => 0,
    // and aa2idx['M'-'A'] => 12

int BLOSUM62[] = {
  4,                                                                  // A
 -1, 5,                                                               // R
 -2, 0, 6,                                                            // N
 -2,-2, 1, 6,                                                         // D
  0,-3,-3,-3, 9,                                                      // C
 -1, 1, 0, 0,-3, 5,                                                   // Q
 -1, 0, 0, 2,-4, 2, 5,                                                // E
  0,-2, 0,-1,-3,-2,-2, 6,                                             // G
 -2, 0, 1,-1,-3, 0, 0,-2, 8,                                          // H
 -1,-3,-3,-3,-1,-3,-3,-4,-3, 4,                                       // I
 -1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,                                    // L
 -1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,                                 // K
 -1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5,                              // M
 -2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,                           // F
 -1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,                        // P
  1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4,                     // S
  0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,                  // T
 -3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11,               // W
 -2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,            // Y
  0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,         // V
 -2,-1, 3, 4,-3, 0, 1,-1, 0,-3,-4, 0,-3,-3,-2, 0,-1,-4,-3,-3, 4,      // B
 -1, 0, 0, 1,-3, 3, 4,-2, 0,-3,-3, 1,-1,-3,-1, 0,-1,-3,-2,-2, 1, 4,   // Z
  0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-1,-1,-1 // X
//A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X
//0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19  2  6 20
};


int na2idx[] = {0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 3, 3, 4, 4, 4, 4, 4};
    // idx for  A  B  C  D  E  F  G  H  I  J  K  L  M  N  O  P
    //          Q  R  S  T  U  V  W  X  Y  Z
    // so aa2idx[ X - 'A'] => idx_of_X, eg aa2idx['A' - 'A'] => 0,
    // and aa2idx['M'-'A'] => 4
int BLOSUM62_na[] = {
  2,                  // A
 -2, 2,               // C
 -2,-2, 2,            // G
 -2,-2,-2, 2,         // T
 -2,-2,-2, 1, 2,      // U
 -2,-2,-2,-2,-2, 1,   // N
  0, 0, 0, 0, 0, 0, 1 // X
//A  C  G  T  U  N  X
//0  1  2  3  3  4  5
};

void setaa_to_na();

struct TempFile
{
	FILE *file;
	char buf[512];

	TempFile( const char *dir = NULL ){
		int len = dir ? strlen( dir ) : 0;
		assert( len < 400 );
		buf[0] = 0;
		if( len ){
			strcat( buf, dir );
			if( buf[len-1] != '/' && buf[len-1] != '\\' ){
				buf[len] = '/';
				len += 1;
			}
		}
		strcat( buf, "cdhit.temp." );
		len += 11;
		sprintf( buf + len, "%p", this );
		file = fopen( buf, "w+" );
	}
	~TempFile(){
		if( file ){
			fclose( file );
			remove( buf );
		}
	}
};

struct TempFiles
{
	NVector<TempFile*> files;

	~TempFiles(){ Clear(); }

	void Clear(){
		int i;
		#pragma omp critical
		{
			for(i=0; i<files.size; i++) if( files[i] ) delete files[i];
			files.Clear();
		}
	}
};

const char *temp_dir = "";
TempFiles temp_files;

FILE* OpenTempFile( const char *dir = NULL )
{
	TempFile *file = new TempFile( dir );
	#pragma omp critical
	{
		temp_files.files.Append( file );
	}
	return file->file;
}
static void CleanUpTempFiles()
{
	if( temp_files.files.Size() ) printf( "Clean up temporary files ...\n" );
	temp_files.Clear();
}

int NAA1 ;
int NAA2 ;
int NAA3 ;
int NAA4 ;
int NAA5 ;
int NAA6 ;
int NAA7 ;
int NAA8 ;
int NAA9 ;
int NAA10;
int NAA11;
int NAA12;
int NAAN_array[13] = { 1 };

void InitNAA( int max )
{
	NAA1  = NAAN_array[1]  = max;
	NAA2  = NAAN_array[2]  = NAA1 * NAA1;
	NAA3  = NAAN_array[3]  = NAA1 * NAA2;
	NAA4  = NAAN_array[4]  = NAA2 * NAA2;
	NAA5  = NAAN_array[5]  = NAA2 * NAA3;
	NAA6  = NAAN_array[6]  = NAA3 * NAA3;
	NAA7  = NAAN_array[7]  = NAA3 * NAA4;
	NAA8  = NAAN_array[8]  = NAA4 * NAA4;
	NAA9  = NAAN_array[9]  = NAA4 * NAA5;
	NAA10 = NAAN_array[10] = NAA5 * NAA5;
	NAA11 = NAAN_array[11] = NAA5 * NAA6;
	NAA12 = NAAN_array[12] = NAA6 * NAA6;
}

extern Options options;
ScoreMatrix mat;
Vector<int> Comp_AAN_idx;

void make_comp_iseq(int len, char *iseq_comp, char *iseq)
{
	int i, c[6] = {3,2,1,0,4,5};
	for (i=0; i<len; i++) iseq_comp[i] = c[ (int)iseq[len-i-1] ];
} // make_comp_iseq

bool Options::SetOptionCommon( const char *flag, const char *value )
{
	int intval = atoi( value );
	if      (strcmp(flag, "-i" ) == 0) input = value;
	else if (strcmp(flag, "-j" ) == 0) input_pe = value;
	else if (strcmp(flag, "-o" ) == 0) output = value;
	else if (strcmp(flag, "-op") == 0) output_pe = value;
	else if (strcmp(flag, "-M" ) == 0) max_memory  = atoll(value) * 1000000;
	else if (strcmp(flag, "-l" ) == 0) min_length  = intval;
	else if (strcmp(flag, "-c" ) == 0) cluster_thd  = atof(value), useIdentity = true;
	else if (strcmp(flag, "-D" ) == 0) distance_thd  = atof(value), useDistance = true;
	else if (strcmp(flag, "-b" ) == 0) band_width  = intval;
	else if (strcmp(flag, "-n" ) == 0) NAA       = intval;
	else if (strcmp(flag, "-d" ) == 0) des_len   = intval;
	else if (strcmp(flag, "-s" ) == 0) diff_cutoff  = atof(value);
	else if (strcmp(flag, "-S" ) == 0) diff_cutoff_aa  = intval;
	else if (strcmp(flag, "-B" ) == 0) store_disk  = intval;
	else if (strcmp(flag, "-P" ) == 0) PE_mode  = intval;
	else if (strcmp(flag, "-cx") == 0) trim_len = intval;
	else if (strcmp(flag, "-cy") == 0) trim_len_R2 = intval;
	else if (strcmp(flag, "-ap") == 0) align_pos = intval;
	else if (strcmp(flag, "-sc") == 0) sort_output = intval;
	else if (strcmp(flag, "-sf") == 0) sort_outputf = intval;
	else if (strcmp(flag, "-p" ) == 0) print  = intval;
	else if (strcmp(flag, "-g" ) == 0) cluster_best  = intval;
	else if (strcmp(flag, "-G" ) == 0) global_identity  = intval;
	else if (strcmp(flag, "-aL") == 0) long_coverage = atof(value);
	else if (strcmp(flag, "-AL") == 0) long_control = intval;
	else if (strcmp(flag, "-aS") == 0) short_coverage = atof(value);
	else if (strcmp(flag, "-AS") == 0) short_control = intval;
	else if (strcmp(flag, "-A" ) == 0) min_control  = intval;
	else if (strcmp(flag, "-uL") == 0) long_unmatch_per = atof(value);
	else if (strcmp(flag, "-uS") == 0) short_unmatch_per = atof(value);
	else if (strcmp(flag, "-U") == 0) unmatch_len = intval;
	else if (strcmp(flag, "-tmp" ) == 0) temp_dir  = value;
	else if (strcmp(flag, "-bak" ) == 0) backupFile  = intval;
	else if (strcmp(flag, "-T" ) == 0){
#ifndef NO_OPENMP
		int cpu = omp_get_num_procs();
		threads  = intval;
		if( threads > cpu ){
			threads = cpu;
			printf( "Warning: total number of CPUs in the system is %i\n", cpu );
		}else if( threads < 0 ){
			threads += cpu;
			if( threads < 0 ) threads = 0;
		}
		if( threads == 0 ){
			threads = cpu;
			printf( "total number of CPUs in the system is %i\n", cpu );
		}
		if( threads != intval ) printf( "Actual number of CPUs to be used: %i\n\n", threads );
#else
		printf( "Option -T is ignored: multi-threading with OpenMP is NOT enabled!\n" );
#endif
	}else return false;
	return true;
}
bool Options::SetOption( const char *flag, const char *value )
{
	if( is454 ){
		if( strcmp(flag, "-s") == 0 ) return false;
		else if( strcmp(flag, "-S") == 0 ) return false;
		else if( strcmp(flag, "-G") == 0 ) return false;
		else if( strcmp(flag, "-A") == 0 ) return false;
		else if( strcmp(flag, "-r") == 0 ) return false;
		else if( strcmp(flag, "-D") == 0 ){ max_indel = atoi(value); return true; }
	}
	if( SetOptionCommon( flag, value ) ) return true;
	if (strcmp(flag, "-t" ) == 0) tolerance = atoi(value);
	else if (strcmp(flag, "-F" ) == 0) frag_size = atoi(value);
	else if (has2D && SetOption2D( flag, value ) ) return true;
	else if (isEST && SetOptionEST( flag, value ) ) return true;
	else return false;
	return true;
}
bool Options::SetOption2D( const char *flag, const char *value )
{
	if( SetOptionCommon( flag, value ) ) return true;
	if (strcmp(flag, "-i2" ) == 0) input2 = value;
        else if (strcmp(flag, "-j2" ) == 0) input2_pe = value;
	else if (strcmp(flag, "-s2") == 0) diff_cutoff2 = atof(value);
	else if (strcmp(flag, "-S2") == 0) diff_cutoff_aa2 = atoi(value);
	else return false;
	return true;
}
bool Options::SetOptionEST( const char *flag, const char *value )
{
	NAA_top_limit = 12;
	if( SetOptionCommon( flag, value ) ) return true;
	if (strcmp(flag, "-r" ) == 0) option_r  = atoi(value); 
	else if (strcmp(flag, "-gap") == 0) mat.gap = MAX_SEQ * atoi(value);
	else if (strcmp(flag, "-gap-ext") == 0) mat.ext_gap = MAX_SEQ * atoi(value);
	else if (strcmp(flag, "-match") == 0) mat.set_match( atoi(value) );
	else if (strcmp(flag, "-mismatch") == 0) mat.set_mismatch( atoi(value) );
	else if (strcmp(flag, "-mask") == 0){
		string letters = value;
		int i, n = letters.size();
		for(i=0; i<n; i++){
			char ch = toupper( letters[i] );
			if( ch < 'A' || ch > 'Z' ) continue;
			na2idx[ ch - 'A' ] = 5;
		}
		setaa_to_na();
	}else return false;
	return true;
}
bool Options::SetOptions( int argc, char *argv[], bool twod, bool est )
{
	int i, n;
	char date[100];
	strcpy( date, __DATE__ );
	n = strlen( date );
	for(i=1; i<n; i++) if( date[i-1] == ' ' and date[i] == ' ' ) date[i] = '0';
	printf( "================================================================\n" );
	printf( "Program: CD-HIT, V" CDHIT_VERSION WITH_OPENMP ", %s, " __TIME__ "\n", date );
	printf( "Command:" );
	n = 9;
	for(i=0; i<argc; i++){
		n += strlen( argv[i] ) + 1;
		if( n >= 64 ){
			printf( "\n         %s", argv[i] );
			n = strlen( argv[i] ) + 9;
		}else{
			printf( " %s", argv[i] );
		}
	}
	printf( "\n\n" );
	time_t tm = time(NULL);
	printf( "Started: %s", ctime(&tm) );
	printf( "================================================================\n" );
	printf( "                            Output                              \n" );
	printf( "----------------------------------------------------------------\n" );
	has2D = twod;
	isEST = est;
	for (i=1; i+1<argc; i+=2) if ( SetOption( argv[i], argv[i+1] ) == 0) return false;
	if( i < argc ) return false;

	atexit( CleanUpTempFiles );
	return true;
}
void Options::Validate()
{
	if( useIdentity and useDistance ) bomb_error( "can not use both identity cutoff and distance cutoff" );
	if( useDistance ){
		if ((distance_thd > 1.0) || (distance_thd < 0.0)) bomb_error("invalid distance threshold");
	}else if( isEST ){
		if ((cluster_thd > 1.0) || (cluster_thd < 0.8)) bomb_error("invalid clstr threshold, should >=0.8");
	}else{
		if ((cluster_thd > 1.0) || (cluster_thd < 0.4)) bomb_error("invalid clstr");
	}

        if (input.size()  == 0) bomb_error("no input file");
        if (output.size() == 0) bomb_error("no output file");
        if (PE_mode) {
          if (input_pe.size()  == 0) bomb_error("no input file for R2 sequences in PE mode");
          if (output_pe.size() == 0) bomb_error("no output file for R2 sequences in PE mode");
        }
        if (isEST && (align_pos==1)) option_r = 0; 

	if (band_width < 1 ) bomb_error("invalid band width");
	if (NAA < 2 || NAA > NAA_top_limit) bomb_error("invalid word length");
	if (des_len < 0 ) bomb_error("too short description, not enough to identify sequences");
	if (not isEST && (tolerance < 0 || tolerance > 5) ) bomb_error("invalid tolerance");
	if ((diff_cutoff<0) || (diff_cutoff>1)) bomb_error("invalid value for -s");
	if (diff_cutoff_aa<0) bomb_error("invalid value for -S");
	if( has2D ){
		if ((diff_cutoff2<0) || (diff_cutoff2>1)) bomb_error("invalid value for -s2");
		if (diff_cutoff_aa2<0) bomb_error("invalid value for -S2");
          if (PE_mode) {
            if (input2_pe.size()  == 0) bomb_error("no input file for R2 sequences for 2nd db in PE mode");
          }
	}
	if (global_identity == 0) print = 1;
	if (short_coverage < long_coverage) short_coverage = long_coverage;
	if (short_control > long_control) short_control = long_control;
	if ((global_identity == 0) && (short_coverage == 0.0) && (min_control == 0))
		bomb_error("You are using local identity, but no -aS -aL -A option");
	if (frag_size < 0) bomb_error("invalid fragment size");

#if 0
	if( useDistance ){
		/* when required_aan becomes zero */
		if( distance_thd * NAA >= 1 )
			bomb_warning( "word length is too long for the distance cutoff" );
	}else{
		/* when required_aan becomes zero */
		if( cluster_thd <= 1.0 - 1.0 / NAA )
			bomb_warning( "word length is too long for the identity cutoff" );
	}
#endif

	const char *message = "Your word length is %i, using %i may be faster!\n";
	if ( not isEST && tolerance ) {
		int i, clstr_idx = (int) (cluster_thd * 100) - naa_stat_start_percent;
		int tcutoff = naa_stat[tolerance-1][clstr_idx][5-NAA];

		if (tcutoff < 5 )
			bomb_error("Too low cluster threshold for the word length.\n"
					"Increase the threshold or the tolerance, or decrease the word length.");
		for ( i=5; i>NAA; i--) {
			if ( naa_stat[tolerance-1][clstr_idx][5-i] > 10 ) {
				printf( message, NAA, i );
				break;
			}
		}
	} else if( isEST ) {
		if      ( cluster_thd > 0.9  && NAA < 8 ) printf( message, NAA, 8 );
		else if ( cluster_thd > 0.87 && NAA < 5 ) printf( message, NAA, 5 );
		else if ( cluster_thd > 0.80 && NAA < 4 ) printf( message, NAA, 4 );
		else if ( cluster_thd > 0.75 && NAA < 3 ) printf( message, NAA, 3 );
	} else {
		if      ( cluster_thd > 0.85 && NAA < 5 ) printf( message, NAA, 5 );
		else if ( cluster_thd > 0.80 && NAA < 4 ) printf( message, NAA, 4 );
		else if ( cluster_thd > 0.75 && NAA < 3 ) printf( message, NAA, 3 );
	}

	if ( (min_length + 1) < NAA ) bomb_error("Too short -l, redefine it");
}
void Options::Print()
{
	printf( "isEST = %i\n", isEST );
	printf( "has2D = %i\n", has2D );
	printf( "NAA = %i\n", NAA );
	printf( "NAA_top_limit = %i\n", NAA_top_limit );
	printf( "min_length = %i\n", min_length );
	printf( "cluster_best = %i\n", cluster_best );
	printf( "global_identity = %i\n", global_identity );
	printf( "cluster_thd = %g\n", cluster_thd );
	printf( "diff_cutoff = %g\n", diff_cutoff );
	printf( "diff_cutoff_aa = %i\n", diff_cutoff_aa );
	printf( "tolerance = %i\n", tolerance );
	printf( "long_coverage = %g\n", long_coverage );
	printf( "long_control = %i\n", long_control );
	printf( "short_coverage = %g\n", short_coverage );
	printf( "short_control = %i\n", short_control );
	printf( "frag_size = %i\n", frag_size );
	printf( "option_r = %i\n", option_r );
	printf( "print = %i\n", print );
}

void bomb_error(const char *message)
{
	fprintf( stderr, "\nFatal Error:\n%s\nProgram halted !!\n\n", message );
	temp_files.Clear();
	exit (1);
} // END void bomb_error

void bomb_error(const char *message, const char *message2)
{
	fprintf( stderr, "\nFatal Error:\n%s %s\nProgram halted !!\n\n", message, message2 );
	temp_files.Clear();
	exit (1);
} // END void bomb_error


void bomb_warning(const char *message)
{
	fprintf( stderr, "\nWarning:\n%s\nNot fatal, but may affect results !!\n\n", message );
} // END void bomb_warning


void bomb_warning(const char *message, const char *message2)
{
	fprintf( stderr, "\nWarning:\n%s %s\nNot fatal, but may affect results !!\n\n", message, message2 );
} // END void bomb_warning

void format_seq(char *seq)
{
	int i, j;
	char c1;
	int len = strlen(seq);

	for (i=0,j=0; i<len; i++) {
		c1 = toupper(seq[i]);
		if ( isalpha(c1) ) seq[j++] = c1;
	}
	seq[j] = 0;
} // END void format_seq

void strrev(char *p)
{
  char *q = p;
  while(q && *q) ++q;
  for(--q; p < q; ++p, --q)
    *p = *p ^ *q,
    *q = *p ^ *q,
    *p = *p ^ *q;
}

////For smiple len1 <= len2, len2 is for existing representative
////walk along all diag path of two sequences,
////find the diags with most aap
////return top n diags
////added on 2006 11 13
////band 0                      XXXXXXXXXXXXXXXXXX               seq2, rep seq
////                            XXXXXXXXXXXXXXX                  seq1
////band 1                      XXXXXXXXXXXXXXXXXX               seq2, rep seq
////                             XXXXXXXXXXXXXXX                 seq1
////extreme right (+)           XXXXXXXXXXXXXXXXXX               seq2, rep seq
////    band = len2-1                            XXXXXXXXXXXXXXX seq1
////band-1                      XXXXXXXXXXXXXXXXXX               seq2, rep seq
////                           XXXXXXXXXXXXXXX                   seq1
////extreme left (-)            XXXXXXXXXXXXXXXXXX               seq2, rep seq
////              XXXXXXXXXXXXXXX   band = -(len1-1)             seq1
////index of diag_score = band+len1-1;
int diag_test_aapn(int NAA1, char iseq2[], int len1, int len2, WorkingBuffer & buffer,
		int &best_sum, int band_width, int &band_left, int &band_center, int &band_right, int required_aa1)
{
	int i, i1, j, k;
	int *pp;
	int nall = len1+len2-1;
	Vector<int> & taap = buffer.taap;
	Vector<INTs> & aap_begin = buffer.aap_begin;
	Vector<INTs> & aap_list = buffer.aap_list;
	Vector<int> & diag_score = buffer.diag_score;
	Vector<int> & diag_score2 = buffer.diag_score2;

	if (nall > MAX_DIAG) bomb_error("in diag_test_aapn, MAX_DIAG reached");
	for (pp=&diag_score[0], i=nall; i; i--, pp++) *pp=0;
	for (pp=&diag_score2[0], i=nall; i; i--, pp++) *pp=0;

	int c22, cpx;
	INTs *bip;
	int len11 = len1-1;
	int len22 = len2-1;
	i1 = len11;
	for (i=0; i<len22; i++,i1++) {
		c22 = iseq2[i]*NAA1+ iseq2[i+1];
		cpx = 1 + (iseq2[i] != iseq2[i+1]);
		if ( (j=taap[c22]) == 0) continue;
		int m = aap_begin[c22];
		for(int k=0; k<j; k++){
			diag_score[ i1 - aap_list[m+k] ] ++;
			diag_score2[ i1 - aap_list[m+k] ] += cpx;
		}
	}

	//find the best band range
	//  int band_b = required_aa1;
	int band_b = required_aa1-1 >= 0 ? required_aa1-1:0;  // on dec 21 2001
	int band_e = nall - band_b;

	int band_m = ( band_b+band_width-1 < band_e ) ? band_b+band_width-1 : band_e;
	int best_score=0;
	int best_score2=0;
	int max_diag = 0;
	int max_diag2 = 0;
	int imax_diag = 0;
	for (i=band_b; i<=band_m; i++){
		best_score += diag_score[i];
		best_score2 += diag_score2[i];
		if( diag_score2[i] > max_diag2 ){
			max_diag2 = diag_score2[i];
			max_diag = diag_score[i];
			imax_diag = i;
		}
	}
	int from=band_b;
	int end =band_m;
	int score = best_score;
	int score2 = best_score2;
	for (k=from, j=band_m+1; j<band_e; j++, k++) {
		score -= diag_score[k]; 
		score += diag_score[j]; 
		score2 -= diag_score2[k]; 
		score2 += diag_score2[j]; 
		if ( score2 > best_score2 ) {
			from = k + 1;
			end  = j;
			best_score = score;
			best_score2 = score2;
			if( diag_score2[j] > max_diag2 ){
				max_diag2 = diag_score2[j];
				max_diag = diag_score[j];
				imax_diag = j;
			}
		}
	}
	int mlen = imax_diag;
	if( imax_diag > len1 ) mlen = nall - imax_diag;
	int emax = int((1.0 - options.cluster_thd) * mlen) + 1;
	for (j=from; j<imax_diag; j++) { // if aap pairs fail to open gap
		if ( (imax_diag - j) > emax || diag_score[j] < 1 ) {
			best_score -= diag_score[j]; from++;
		} else break;
	}
	for (j=end; j>imax_diag; j--) { // if aap pairs fail to open gap
		if ( (j - imax_diag) > emax || diag_score[j] < 1 ) {
			best_score -= diag_score[j]; end--;
		} else break;
	}

	//  delete [] diag_score;
	band_left = from - len1 + 1; 
	band_right= end - len1 + 1;
	band_center = imax_diag - len1 + 1;
	best_sum = best_score;
	return OK_FUNC;
}
// END diag_test_aapn
 

int diag_test_aapn_est(int NAA1, char iseq2[], int len1, int len2, WorkingBuffer & buffer, 
        int &best_sum, int band_width, int &band_left, int &band_center, int &band_right, int required_aa1)
{
	int i, i1, j, k;
	int *pp, *pp2;
	int nall = len1+len2-1;
	int NAA2 = NAA1 * NAA1;
	int NAA3 = NAA2 * NAA1;
	Vector<int> & taap = buffer.taap;
	Vector<INTs> & aap_begin = buffer.aap_begin;
	Vector<INTs> & aap_list = buffer.aap_list;
	Vector<int> & diag_score = buffer.diag_score;
	Vector<int> & diag_score2 = buffer.diag_score2;

	if (nall > MAX_DIAG) bomb_error("in diag_test_aapn_est, MAX_DIAG reached");
	pp = & diag_score[0];
	pp2 = & diag_score2[0];
	for (i=nall; i; i--, pp++, pp2++) *pp = *pp2 =0;

	INTs *bip;
	int c22, cpx;
	int len22 = len2-3;
	i1 = len1-1;
	for (i=0; i<len22; i++,i1++,iseq2++) {
		unsigned char c0 = iseq2[0];
		unsigned char c1 = iseq2[1];
		unsigned char c2 = iseq2[2];
		unsigned char c3 = iseq2[3];
		if ((c0>=4) || (c1>=4) || (c2>=4) || (c3>=4)) continue; //skip N

		c22 = c0*NAA3+ c1*NAA2 + c2*NAA1 + c3;
		if ( (j=taap[c22]) == 0) continue;
		cpx = 1 + (c0 != c1) + (c1 != c2) + (c2 != c3);
		bip = & aap_list[ aap_begin[c22] ];     //    bi = aap_begin[c22];
		for (; j; j--, bip++) { 
			diag_score[i1 - *bip]++;
			diag_score2[i1 - *bip] += cpx;
		}
	}
#if 0
	int mmax = 0;
	int immax = 0;
	for(i=0; i<=nall; i++){
		if( i%len2 ==0 or i == nall ) printf( "\n" );
		printf( "%3i ", diag_score[i] );
		if( diag_score[i] > mmax ){
			mmax = diag_score[i];
			immax = i;
		}
	}
#endif
	
	//find the best band range
	//  int band_b = required_aa1;
	int band_b = required_aa1-1 >= 0 ? required_aa1-1:0;  // on dec 21 2001
	int band_e = nall - band_b;

	if( options.is454 ){
		band_b = len1 - band_width;
		band_e = len1 + band_width;
		if( band_b < 0 ) band_b = 0;
		if( band_e > nall ) band_e = nall;
	}

	int band_m = ( band_b+band_width-1 < band_e ) ? band_b+band_width-1 : band_e;
	int best_score=0;
	int best_score2=0;
	int max_diag = 0;
	int max_diag2 = 0;
	int imax_diag = 0;
	for (i=band_b; i<=band_m; i++){
		best_score += diag_score[i];
		best_score2 += diag_score2[i];
		if( diag_score2[i] > max_diag2 ){
			max_diag2 = diag_score2[i];
			max_diag = diag_score[i];
			imax_diag = i;
		}
	}
	int from=band_b;
	int end =band_m;
	int score = best_score;  
	int score2 = best_score2;  

	for (k=from, j=band_m+1; j<band_e; j++, k++) {
		score -= diag_score[k]; 
		score += diag_score[j]; 
		score2 -= diag_score2[k]; 
		score2 += diag_score2[j]; 
		if ( score2 > best_score2 ) {
			from = k + 1;
			end  = j;
			best_score = score;
			best_score2 = score2;
			if( diag_score2[j] > max_diag2 ){
				max_diag2 = diag_score2[j];
				max_diag = diag_score[j];
				imax_diag = j;
			}
		}
	}
#if 0
	printf( "%i %i\n", required_aa1, from );
	printf( "max=%3i  imax=%3i; band:  %3i  %3i  %i\n", max_diag, imax_diag, band_b, band_e, band_m );
	printf( "best: %i\n", best_score );
	printf( "from: %i, end: %i,  best: %i\n", from, end, best_score );
#endif
	int mlen = imax_diag;
	if( imax_diag > len1 ) mlen = nall - imax_diag;
	int emax = int((1.0 - options.cluster_thd) * mlen) + 1;
	for (j=from; j<imax_diag; j++) { // if aap pairs fail to open gap
		if ( (imax_diag - j) > emax || diag_score[j] < 1 ) {
			best_score -= diag_score[j]; from++;
		} else break;
	}
	for (j=end; j>imax_diag; j--) { // if aap pairs fail to open gap
		if ( (j - imax_diag) > emax || diag_score[j] < 1 ) {
			best_score -= diag_score[j]; end--;
		} else break;
	}

	band_left = from-len1+1; 
	band_right= end-len1+1;
	band_center = imax_diag - len1 + 1;
	best_sum = best_score;
	if( options.is454 ){
		if( band_left > 0 ) best_sum = 0;
		if( band_right < 0 ) best_sum = 0;
	}
#if 0
	printf( "%3i:  best: %i,  %i  %i  %i\n", required_aa1, best_score, band_left, band_right, band_width );
	printf( "max=%3i  imax=%3i; band:  %3i  %3i  %i\n", mmax, immax, band_b, band_e, band_m );
#endif
	return OK_FUNC;
}
// END diag_test_aapn_est


/*
local alignment of two sequence within a diag band
for band 0 means direction (0,0) -> (1,1)
         1 means direction (0,1) -> (1,2)
        -1 means direction (1,0) -> (2,1)
added on 2006 11 13
band 0                      XXXXXXXXXXXXXXXXXX               seq2, rep seq
                            XXXXXXXXXXXXXXX                  seq1
band 1                      XXXXXXXXXXXXXXXXXX               seq2, rep seq
                             XXXXXXXXXXXXXXX                 seq1
extreme right (+)           XXXXXXXXXXXXXXXXXX               seq2, rep seq
    band = len2-1                            XXXXXXXXXXXXXXX seq1
band-1                      XXXXXXXXXXXXXXXXXX               seq2, rep seq
                           XXXXXXXXXXXXXXX                   seq1
extreme left (-)            XXXXXXXXXXXXXXXXXX               seq2, rep seq
              XXXXXXXXXXXXXXX   band = -(len1-1)             seq1
iseq len are integer sequence and its length,
mat is matrix, return ALN_PAIR class

       band:  -101   seq2 len2 = 17
                \\\1234567890123456
              0  \xxxxxxxxxxxxxxxxx
              1   xxxxxxxxxxxxxxxxx\ most right band = len2-1
              2   xxxxxxxxxxxxxxxxx
    seq1      3   xxxxxxxxxxxxxxxxx
    len1 = 11 4   xxxxxxxxxxxxxxxxx
              5   xxxxxxxxxxxxxxxxx
              6   xxxxxxxxxxxxxxxxx
              7   xxxxxxxxxxxxxxxxx
              8   xxxxxxxxxxxxxxxxx
              9   xxxxxxxxxxxxxxxxx
              0   xxxxxxxxxxxxxxxxx
                  \
                   most left band = -(len1-1)

*/

int local_band_align( char iseq1[], char iseq2[], int len1, int len2, ScoreMatrix &mat, 
		int &best_score, int &iden_no, int &alnln, float &dist, int *alninfo,
		int band_left, int band_center, int band_right, WorkingBuffer & buffer)
{
	int i, j, k, j1;
	int jj, kk;
	int iden_no1;
	int64_t best_score1;
	iden_no = 0;

	if ( (band_right >= len2 ) ||
			(band_left  <= -len1) ||
			(band_left  > band_right) ) return FAILED_FUNC;

	// allocate mem for score_mat[len1][len2] etc
	int band_width = band_right - band_left + 1;
	int band_width1 = band_width + 1;

    // score_mat, back_mat [i][j]: i index of seqi (0 to len(seqi)-1), j index of band (0 to band_width-1)
	MatrixInt64 & score_mat = buffer.score_mat;
	MatrixInt   & back_mat = buffer.back_mat;

	//printf( "%i  %i\n", band_right, band_left );

	if( score_mat.size() <= len1 ){
		VectorInt   row( band_width1, 0 );
		VectorInt64 row2( band_width1, 0 );
		while( score_mat.size() <= len1 ){
			score_mat.Append( row2 );
			back_mat.Append( row );
		}
	}
	for(i=0; i<=len1; i++){
		if( score_mat[i].Size() < band_width1 ) score_mat[i].Resize( band_width1 );
		if( back_mat[i].Size() < band_width1 ) back_mat[i].Resize( band_width1 );
	}

	best_score = 0;
	/* seq1 is query, seq2 is rep
                  seq2    len2 = 17       seq2    len2 = 17    seq2    len2 = 17
                  01234567890123456       01234567890123456    01234567890123456
       0          xxxxxxxxxxxxxxxxx \\\\\\XXXxxxxxxxxxxxxxx    xXXXXXXXxxxxxxxxx
       1     \\\\\Xxxxxxxxxxxxxxxxx  \\\\\Xxx\xxxxxxxxxxxxx    xx\xxxxx\xxxxxxxx
       2      \\\\X\xxxxxxxxxxxxxxx   \\\\Xxxx\xxxxxxxxxxxx    xxx\xxxxx\xxxxxxx
  seq1 3       \\\Xx\xxxxxxxxxxxxxx    \\\Xxxxx\xxxxxxxxxxx    xxxx\xxxxx\xxxxxx
  len1 4        \\Xxx\xxxxxxxxxxxxx     \\Xxxxxx\xxxxxxxxxx    xxxxx\xxxxx\xxxxx
  = 11 5         \Xxxx\xxxxxxxxxxxx      \Xxxxxxx\xxxxxxxxx    xxxxxx\xxxxx\xxxx
       6          Xxxxx\xxxxxxxxxxx       Xxxxxxxx\xxxxxxxx    xxxxxxx\xxxxx\xxx
       7          x\xxxx\xxxxxxxxxx       x\xxxxxxx\xxxxxxx    xxxxxxxx\xxxxx\xx
       8          xx\xxxx\xxxxxxxxx       xx\xxxxxxx\xxxxxx    xxxxxxxxx\xxxxx\x
       9          xxx\xxxx\xxxxxxxx       xxx\xxxxxxx\xxxxx    xxxxxxxxxx\xxxxx\
       0          xxxx\xxxx\xxxxxxx       xxxx\xxxxxxx\xxxx    xxxxxxxxxxx\xxxxx
                  band_left < 0 (-6)      band_left < 0 (-6)   band_left >=0 (1)
                  band_right < 0 (-1)     band_right >=0 (2)   band_right >=0(7)
                  band_width 6            band_width 9         band_width 7
       init score_mat, and iden_mat (place with upper 'X')
     */

	if (band_left < 0) {  //set score to left border of the matrix within band
		int tband = (band_right < 0) ? band_right : 0;
		//for (k=band_left; k<tband; k++)
		for (k=band_left; k<=tband; k++) { // fixed on 2006 11 14
			i = -k;
			j1 = k-band_left;
			// penalty for leading gap opening = penalty for gap extension
            // each of the left side query hunging residues give ext_gap (-1)
			score_mat[i][j1] =  mat.ext_gap * i;
			back_mat[i][j1] = DP_BACK_TOP;
		}
		back_mat[-tband][tband-band_left] = DP_BACK_NONE;
	}

	if (band_right >=0) { //set score to top border of the matrix within band
		int tband = (band_left > 0) ? band_left : 0;
		for (j=tband; j<=band_right; j++) {
			j1 = j-band_left;
			score_mat[0][j1] = mat.ext_gap * j;
			back_mat[0][j1] = DP_BACK_LEFT;
		}
		back_mat[0][tband-band_left] = DP_BACK_NONE;
	}

	int gap_open[2] = { mat.gap, mat.ext_gap };
	int max_diag = band_center - band_left;
	int extra_score[4] = { 4, 3, 2, 1 };
	for (i=1; i<=len1; i++) {
		int J0 = 1 - band_left - i;
		int J1 = len2 - band_left - i;
		if( J0 < 0 ) J0 = 0;
		if( J1 >= band_width ) J1 = band_width;
		for (j1=J0; j1<=J1; j1++){
			j = j1+i+band_left;

			int ci = iseq1[i-1];
			int cj = iseq2[j-1];
			int sij = mat.matrix[ci][cj];
			//int iden_ij = (ci == cj);
			int s1, k0, back;

			/* extra score according to the distance to the best diagonal */
			int extra = extra_score[ abs(j1 - max_diag) & 3 ]; // max distance 3
			sij += extra * (sij>0);

			back = DP_BACK_LEFT_TOP;
			best_score1 = score_mat[i-1][j1] + sij;
			int gap0 = gap_open[ (i == len1) | (j == len2) ];
			int gap = 0;
			int64_t score;

			if( j1 > 0 ){
				gap = gap0;
				if( back_mat[i][j1-1] == DP_BACK_LEFT ) gap = mat.ext_gap;
				if( (score = score_mat[i][j1-1] + gap) > best_score1 ){
					back = DP_BACK_LEFT;
					best_score1 = score;
				}
			}
			if(j1+1<band_width){
				gap = gap0;
				if( back_mat[i-1][j1+1] == DP_BACK_TOP ) gap = mat.ext_gap;
				if( (score = score_mat[i-1][j1+1] + gap) > best_score1 ){
					back = DP_BACK_TOP;
					best_score1 = score;
				}
			}
			score_mat[i][j1] = best_score1;
			back_mat[i][j1]  = back;
			//printf( "%2i(%2i) ", best_score1, iden_no1 );

		}
		//printf( "\n" );
	}
	i = j = 0;
	if( len2 - band_left < len1 ){
		i = len2 - band_left;
		j = len2;
	}else if( len1 + band_right < len2 ){
		i = len1;
		j = len1 + band_right;
	}else{
		i = len1;
		j = len2;
	}
	j1 = j - i - band_left;
	best_score = score_mat[i][j1];
	best_score1 = score_mat[i][j1];

#if 1
	const char *letters = "acgtnx";
	const char *letters2 = "ACGTNX";
#else
	const char *letters = "arndcqeghilkmfpstwyvbzx";
	const char *letters2 = "ARNDCQEGHILKMFPSTWYVBZX";
#endif
	int back = back_mat[i][j1];
	int last = back;
	int count = 0, count2 = 0, count3 = 0;
	int match, begin1, begin2, end1, end2;
	int gbegin1=0, gbegin2=0, gend1=0, gend2=0;
	int64_t score, smin = best_score1, smax = best_score1 - 1;
	int posmin, posmax, pos = 0;
	int bl, dlen = 0, dcount = 0;
	posmin = posmax = 0;
	begin1 = begin2 = end1 = end2 = 0;

#ifdef PRINT
#define PRINT
	printf( "%i %i\n", best_score, score_mat[i][j1] );
	printf( "%i %i %i\n", band_left, band_center, band_right );
	printf( "%i %i %i %i\n", i, j, j1, len2 );
#endif
#ifdef MAKEALIGN
#define MAKEALIGN
	char AA[ MAX_SEQ ], BB[ MAX_SEQ ];
	int NN = 0;
	int IA, IB;
	for(IA=len1;IA>i;IA--){
		AA[NN] = letters[ iseq1[IA-1] ];
		BB[NN++] = '-';
	}
	for(IB=len2;IB>j;IB--){
		AA[NN] = '-';
		BB[NN++] = letters[ iseq2[IB-1] ];
	}
#endif

	int masked = 0;
	int indels = 0;
	int max_indels = 0;
	while( back != DP_BACK_NONE ){
		switch( back ){
		case DP_BACK_TOP  :
#ifdef PRINT
			printf( "%5i: %c %c %9i\n", pos, letters[ iseq1[i-1] ], '|', score_mat[i][j1] );
#endif
#ifdef MAKEALIGN
			AA[NN] = letters[ iseq1[i-1] ];
			BB[NN++] = '-';
#endif
			bl = (last != back) & (j != 1) & (j != len2);
			dlen += bl;
			dcount += bl;
			score = score_mat[i][j1];
			if( score < smin ){
				count2 = 0;
				smin = score;
				posmin = pos - 1;
				begin1 = i;
				begin2 = j;
			}
			i -= 1;
			j1 += 1;
			break;
		case DP_BACK_LEFT :
#ifdef PRINT
			printf( "%5i: %c %c %9i\n", pos, '|', letters[ iseq2[j-1] ], score_mat[i][j1] );
#endif
#ifdef MAKEALIGN
			AA[NN] = '-';
			BB[NN++] = letters[ iseq2[j-1] ];
#endif
			bl = (last != back) & (i != 1) & (i != len1);
			dlen += bl;
			dcount += bl;
			score = score_mat[i][j1];
			if( score < smin ){
				count2 = 0;
				smin = score;
				posmin = pos - 1;
				begin1 = i;
				begin2 = j;
			}
			j1 -= 1;
			j -= 1;
			break;
		case DP_BACK_LEFT_TOP :
#ifdef PRINT
			if( iseq1[i-1] == iseq2[j-1] ){
				printf( "%5i: %c %c %9i\n", pos, letters2[ iseq1[i-1] ], letters2[ iseq2[j-1] ], score_mat[i][j1] );
			}else{
				printf( "%5i: %c %c %9i\n", pos, letters[ iseq1[i-1] ], letters[ iseq2[j-1] ], score_mat[i][j1] );
			}
#endif
#ifdef MAKEALIGN
			if( iseq1[i-1] == iseq2[j-1] ){
				AA[NN] = letters2[ iseq1[i-1] ];
				BB[NN++] = letters2[ iseq2[j-1] ];
			}else{
				AA[NN] = letters[ iseq1[i-1] ];
				BB[NN++] = letters[ iseq2[j-1] ];
			}
#endif
			if( alninfo && options.global_identity ){
				if( i == 1 || j == 1 ){
					gbegin1 = i-1;
					gbegin2 = j-1;
				}else if( i == len1 || j == len2 ){
					gend1 = i-1;
					gend2 = j-1;
				}
			}
			score = score_mat[i][j1];
			i -= 1;
			j -= 1;
			match = iseq1[i] == iseq2[j];
			if( score > smax ){
				count = 0;
				smax = score;
				posmax = pos;
				end1 = i;
				end2 = j;
			}
			if( options.isEST && (iseq1[i] > 4 || iseq2[j] > 4) ){
				masked += 1;
			}else{
				dlen += 1;
				dcount += ! match;
				count += match;
				count2 += match;
				count3 += match;
			}
			if( score < smin ){
				int mm = match == 0;
				count2 = 0;
				smin = score;
				posmin = pos - mm;
				begin1 = i + mm;
				begin2 = j + mm;
			}
			break;
		default : printf( "%i\n", back ); break;
		}
		if( options.is454 ){
			if( back == DP_BACK_LEFT_TOP ){
				if( indels > max_indels ) max_indels = indels;
				indels = 0;
			}else{
				if( last == DP_BACK_LEFT_TOP ){
					indels = 1;
				}else if( indels ){
					indels += 1;
				}
			}
		}
		pos += 1;
		last = back;
		back = back_mat[i][j1];
	}
	if( options.is454 and max_indels > options.max_indel ) return FAILED_FUNC;
	iden_no = options.global_identity ? count3 : count - count2;
	alnln = posmin - posmax + 1 - masked;
	dist = dcount/(float)dlen;
	//dist = - 0.75 * log( 1.0 - dist * 4.0 / 3.0 );
	int umtail1 = len1 - 1 - end1;
	int umtail2 = len2 - 1 - end2;
	int umhead = begin1 < begin2 ? begin1 : begin2;
	int umtail = umtail1 < umtail2 ? umtail1 : umtail2;
	int umlen = umhead + umtail;
	if( umlen > options.unmatch_len ) return FAILED_FUNC;
	if( umlen > len1 * options.short_unmatch_per ) return FAILED_FUNC;
	if( umlen > len2 * options.long_unmatch_per ) return FAILED_FUNC;
	if( alninfo ){
		alninfo[0] = begin1;
		alninfo[1] = end1;
		alninfo[2] = begin2;
		alninfo[3] = end2;
		alninfo[4] = masked;
		if( options.global_identity ){
			alninfo[0] = gbegin1;
			alninfo[1] = gend1;
			alninfo[2] = gbegin2;
			alninfo[3] = gend2;
		}
	}
#ifdef PRINT
	printf( "%6i %6i:  %4i %4i %4i %4i\n", alnln, iden_no, begin1, end1, begin2, end2 );
	printf( "%6i %6i:  %4i %4i\n", posmin, posmax, posmin - posmax, count - count2 );
	printf( "smin = %9i, smax = %9i\n", smin, smax );
	printf( "dlen = %5i, dcount = %5i, dist = %.3f\n", dlen, dcount, dcount/(float)dlen );
#endif
#ifdef MAKEALIGN
	float identity = iden_no / (float)( options.global_identity ? (len1 - masked) : alnln);
	if( identity < options.cluster_thd ) return OK_FUNC;
	while(i--){
		AA[NN] = letters[ iseq1[i-1] ];
		BB[NN++] = '-';
	}
	while(j--){
		AA[NN] = '-';
		BB[NN++] = letters[ iseq2[j-1] ];
	}
	AA[NN] = '\0';
	BB[NN] = '\0';
	for(i=0; i<NN/2; i++){
		char aa = AA[i], bb = BB[i];
		AA[i] = AA[NN-i-1];
		BB[i] = BB[NN-i-1];
		AA[NN-i-1] = aa;
		BB[NN-i-1] = bb;
	}
	static int fcount = 0; 
	fcount += 1;
	FILE *fout = fopen( "alignments.txt", "a" );
	if( fout == NULL ){
		if( fcount <= 1 ) printf( "alignment files open failed\n" );
		return OK_FUNC;
	}
	fprintf( fout, "\n\n######################################################\n" );
	fprintf( fout, "# length X = %i\n", len2 );
	fprintf( fout, "# length Y = %i\n", len1 );
	fprintf( fout, "# best align X: %i-%i\n", begin2+1, end2+1 );
	fprintf( fout, "# best align Y: %i-%i\n", begin1+1, end1+1 );
	if( alninfo ){
		fprintf( fout, "# align X: %i-%i\n", alninfo[2]+1, alninfo[3]+1 );
		fprintf( fout, "# align Y: %i-%i\n", alninfo[0]+1, alninfo[1]+1 );
	}
	fprintf( fout, "# alignment length: %i\n", alnln );
	fprintf( fout, "# identity count: %i\n", iden_no );
	fprintf( fout, "# identity: %g\n", identity );
	fprintf( fout, "# distance: %g\n", dist );
	if( options.is454 ) fprintf( fout, "# max indel: %i\n", max_indels );
#if 0
	fprintf( fout, "%i %s\n", seq1->index, AA );
	fprintf( fout, "%i %s\n", seq2->index, BB );
#else
	bool printaa = true;
	IB = IA = 0;
	fprintf( fout, "\n\nX " );
	while( IA < NN ){
		if( printaa ){
			fprintf( fout, "%c", BB[IB] );
			IB += 1;
			if( IB % 75 ==0 or IB == NN ) printaa = false, fprintf( fout, "\nY " );
		}else{
			fprintf( fout, "%c", AA[IA] );
			IA += 1;
			if( IA % 75 ==0 ) printaa = true, fprintf( fout, "\n\nX " );
		}
	}
#endif
	fclose( fout );
#endif

	return OK_FUNC;
} // END int local_band_align



void setaa_to_na()
{
	int i;
	for (i=0; i<26; i++) aa2idx[i]   = na2idx[i];
} // END void setaa_to_na


/////////////////
ScoreMatrix::ScoreMatrix()
{
	init();
}

void ScoreMatrix::init() 
{
	set_gap( -11, -1 );
	set_matrix( BLOSUM62 );
}

void ScoreMatrix::set_gap(int gap1, int ext_gap1)
{
	int i;
	gap = MAX_SEQ * gap1;
	ext_gap = MAX_SEQ * ext_gap1;
}

void ScoreMatrix::set_matrix(int *mat1)
{
	int i, j, k;
	k = 0;
	for ( i=0; i<MAX_AA; i++)
		for ( j=0; j<=i; j++)
			matrix[j][i] = matrix[i][j] = MAX_SEQ * mat1[ k++ ];
}

void ScoreMatrix::set_to_na()
{
	set_gap( -6, -1 );
	set_matrix( BLOSUM62_na );
}
// Only for est
void ScoreMatrix::set_match( int score )
{
	int i;
	for ( i=0; i<5; i++) matrix[i][i] = MAX_SEQ * score;
	//matrix[3][4] = matrix[4][3] = MAX_SEQ * score;
}
// Only for est
void ScoreMatrix::set_mismatch( int score )
{
	int i, j;
	for ( i=0; i<MAX_AA; i++)
		for ( j=0; j<i; j++)
			matrix[j][i] = matrix[i][j] = MAX_SEQ * score;
	matrix[3][4] = matrix[4][3] = MAX_SEQ;
}

WordTable::WordTable( int naa, int naan )
{
	NAA      = 0;
	NAAN     = 0;
	is_aa    = 1;
	size = 0;
	frag_count = 0;
	Init( naa, naan );
}

void WordTable::SetDNA()
{
	is_aa = 0;
}

void WordTable::Init(int naa, int naan)
{
	NAA  = naa;
	NAAN = naan;
	indexCounts.resize( NAAN );
}

void WordTable::Clear()
{
	int i;
#if 0
	int n1 = 0, n2 = 0, n3 = 0, ns = 0;
	for(i=0; i<NAAN; i++){
		NVector<IndexCount> & ics = indexCounts[i];
		for(int j=0; j<ics.size; j++){
			IndexCount ic = ics[j];
			n1 += ic.count == 1;
			n2 += ic.count == 2;
			n3 += ic.count == 3;
			ns += ic.count >= 4;
		}
	}
	printf( "%9i %9i %9i %9i\n", n1, n2, n3, ns );
#endif
	size = 0;
	frag_count = 0;
	sequences.clear();
	for (i=0; i<NAAN; i++) indexCounts[i].size = 0;//Clear();
}

int WordTable::AddWordCounts( NVector<IndexCount> & counts, Sequence *seq, bool skipN)
{
	int aan_no = counts.Size();
	int i, j, k, idx = sequences.size();
	for (i=0; i<aan_no; i++) {
		IndexCount ic = counts[i];
		if ( (k=ic.count) ) {
			j = ic.index;
			if ( skipN && j<0) continue; // for those has 'N'
			NVector<IndexCount> & row = indexCounts[j];
			ic.index = idx;
			row.Append( ic );
			size += 1;
		}
	}
	sequences.Append( seq );
	return OK_FUNC;
}
int WordTable::AddWordCountsFrag( NVector<IndexCount> & counts, int frag, int frag_size, int repfrag )
{
	return 0;
}
int WordTable::AddWordCounts(int aan_no, Vector<int> & word_encodes, Vector<INTs> & word_encodes_no, int idx, bool skipN)
{
	int i, j, k;
	//printf( "seq %6i: ", idx );
	for (i=0; i<aan_no; i++) {
		if ( (k=word_encodes_no[i]) ) {
			j = word_encodes[i];
			if( skipN && j<0) continue; // for those has 'N'
			NVector<IndexCount> & row = indexCounts[j];
			row.Append( IndexCount( idx, k ) );
			size += 1;
			//if( k >1 ) printf( " %3i", k );
		}
	}
	//printf( "\n" );
	return OK_FUNC;
}

int WordTable::AddWordCountsFrag( int aan_no, Vector<int> & word_encodes,
    Vector<INTs> & word_encodes_no, int frag, int frag_size )
{
	int i, j, k, i1, k1, fra;

	for (i=0; i<frag; i++) {
		k = (i+1)*frag_size < aan_no ? (i+1)*frag_size-1: aan_no-1;
		//quick_sort(&word_encodes[0], i*frag_size, k);
		std::sort( word_encodes.begin() + i*frag_size, word_encodes.begin() + k + 1 );
	}
	for(j=aan_no-1; j; j--) {
		if (word_encodes[j] == word_encodes[j-1]) {
			word_encodes_no[j-1] += word_encodes_no[j];
			word_encodes_no[j]=0;
		}
	}
	// END check_word_encodes

	for (i=0; i<aan_no; i+=frag_size) {
		k = frag_size < (aan_no-i) ? frag_size : (aan_no -i);
		fra = i / frag_size;
		//AddWordCounts(k, word_encodes+i, word_encodes_no+i, NR90f_no+fra);
		for (i1=i; i1<i+k; i1++) {
			if ( (k1=word_encodes_no[i1]) ) {
				j = word_encodes[i1];
				NVector<IndexCount> & row = indexCounts[j];
				row.Append( IndexCount( frag_count + fra, k1 ) );
				size += 1;
			}
		}
	}
	frag_count += frag;

	return 0;
}


void WordTable::PrintAll()
{
	int  i, j, k;
	int cols = 0;
	long long total_words = 0;
	k = 0;
	for (i=0; i<NAAN; i++) {
		int size = indexCounts[i].Size();
		if ( size == 0 ) continue;
		cols++;
		cout << k << "\t" << i << "\tsize:" << size << "\t";
		for (j=0; j<size; j++) {
			cout << indexCounts[i][j].index << "," << indexCounts[i][j].count << " ";
			total_words += indexCounts[i][j].count;
		}
		cout << endl;
		k++;
	}

	cout << "total cols: " << cols << " total words: " << total_words << endl;
}


/* Quick Sort.
 * Adam Drozdek: Data Structures and Algorithms in C++, 2nd Edition.
 */
void PartialQuickSort( IndexCount *data, int first, int last, int partial )
{
	int lower=first+1, upper=last;
	IndexCount pivot;
	IndexCount val;
	if( first >= last ) return;
	val = data[first];
	data[first] = data[ (first+last)/2 ];
	data[ (first+last)/2 ] = val;
	pivot = data[ first ];

	while( lower <= upper ){
		while( lower <= last && data[lower].count < pivot.count ) lower ++;
		while( pivot.count < data[upper].count ) upper --;
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
	if( first < upper-1 ) PartialQuickSort( data, first, upper-1, partial );
	if( upper >= partial ) return;
	if( upper+1 < last ) PartialQuickSort( data, upper+1, last, partial );
}
int WordTable::CountWords(int aan_no, Vector<int> & word_encodes, Vector<INTs> & word_encodes_no,
    NVector<IndexCount> &lookCounts, NVector<uint32_t> & indexMapping, 
	bool est, int min)
{
	int S = frag_count ? frag_count : sequences.size();
	int  j, k, j0, j1, k1, m;
	int ix1, ix2, ix3, ix4;
	IndexCount tmp;

	IndexCount *ic = lookCounts.items;
	for(j=0; j<lookCounts.size; j++, ic++) indexMapping[ ic->index ] = 0;
	lookCounts.size = 0;

	int *we = & word_encodes[0];
	j0 = 0;
	if( est ) while( *we <0 ) j0++, we++; // if met short word has 'N'
	INTs *wen = & word_encodes_no[j0];
	//printf( "\nquery : " );
	for (; j0<aan_no; j0++, we++, wen++) {
		j  = *we;
		j1 = *wen;
		//if( j1 >1 ) printf( " %3i", j1 );
		if( j1==0 ) continue;
		NVector<IndexCount> & one = indexCounts[j];
		k1 = one.Size();
		IndexCount *ic = one.items;

		int rest = aan_no - j0 + 1;
		for (k=0; k<k1; k++, ic++){
			int c = ic->count < j1 ? ic->count : j1;
			uint32_t *idm = indexMapping.items + ic->index;
			if( *idm ==0 ){
				if( rest < min ) continue;
				IndexCount *ic2 = lookCounts.items + lookCounts.size;
				lookCounts.size += 1;
				*idm = lookCounts.size;
				ic2->index = ic->index;
				ic2->count = c;
			}else{
				lookCounts[ *idm - 1 ].count += c;
			}
		}
	}
	//printf( "%6i %6i\n", S, lookCounts.size );
	lookCounts[ lookCounts.size ].count = 0;
	//printf( "\n\n" );
	return OK_FUNC;
}

Sequence::Sequence()
{
	memset( this, 0, sizeof( Sequence ) );
	distance = 2.0;
}
Sequence::Sequence( const Sequence & other )
{
	int i;
	//printf( "new: %p  %p\n", this, & other );
	memcpy( this, & other, sizeof( Sequence ) );
	distance = 2.0;
	if( other.data ){
		size = bufsize = other.size;
                size_R2 = 0;
		data = new char[size+1];
		//printf( "data: %p  %p\n", data, other.data );
		data[size] = 0;
		memcpy( data, other.data, size );
		//for (i=0; i<size; i++) data[i] = other.data[i];
	}
	if( other.identifier ){
		int len = strlen( other.identifier );
		identifier = new char[len+1];
		memcpy( identifier, other.identifier, len );
		identifier[len] = 0;
	}
}

// back to back merge for PE
// R1 -> XXXXXXABC ------------------- NMLYYYYYY <--R2
// >R1           >R2
// XXXXXXABC     YYYYYYLMN =====> Merge into
// >R12
// NMLYYYYYYXXXXXXABC
Sequence::Sequence( const Sequence & other, const Sequence & other2, int mode )
{
	int i;
        if (mode != 1) bomb_error("unknown mode");

	//printf( "new: %p  %p\n", this, & other );
	memcpy( this, & other, sizeof( Sequence ) );
	distance = 2.0;

	if( other.data && other2.data ){
		size = bufsize = (other.size + other2.size);
                size_R2 = other2.size;
		data = new char[size+1];
		//printf( "data: %p  %p\n", data, other.data );
		data[size] = 0;     
                data[size_R2] = 0;  
                memcpy( data, other2.data, size_R2); // copy R2 first
                strrev( data );                      // reverse R2 on data
		memcpy( data+size_R2, other.data, size-size_R2 ); // copy R1 to end of R2
		//for (i=0; i<size; i++) data[i] = other.data[i];
		des_begin2 = other2.des_begin;
                tot_length2= other2.tot_length;
	}
        else if ( other.data || other2.data ) {
                bomb_error("Not both PE sequences have data");
        }

	if( other.identifier ){ // only use R1
		int len = strlen( other.identifier );
		identifier = new char[len+1];
		memcpy( identifier, other.identifier, len );
		identifier[len] = 0;
	}
}


Sequence::~Sequence()
{
	//printf( "delete: %p\n", this );
	if( data ) delete[] data;
	if( identifier ) delete[] identifier;
}

void Sequence::Clear()
{
	if( data ) delete[] data;
	/* do not set size to zero here, it is need for writing output */
	bufsize = 0;
	data = NULL;
}

void Sequence::operator=( const char *s )
{
	size = 0; // avoid copying;
	Resize( strlen( s ) );
	strcpy( data, s );
}
void Sequence::operator+=( const char *s )
{
	int i, m = size, n = strlen( s );
	Reserve( m + n );
	memcpy( data+m, s, n );
}
void Sequence::Resize( int n )
{
	int i, m = size < n ? size : n;
	size = n;
	if( size != bufsize ){
		char *old = data;
		bufsize = size;
		data = new char[ bufsize + 1 ];
		if ( data == NULL ) bomb_error( "Memory" );
		if ( old ){
			memcpy( data, old, m );
			delete []old;
		}
		if( size ) data[size] = 0;
	}
}
void Sequence::Reserve( int n )
{
	int i, m = size < n ? size : n;
	size = n;
	if( size > bufsize ){
		char *old = data;
		bufsize = size + size/5 + 1;
		data = new char[ bufsize + 1 ];
		if ( data == NULL ) bomb_error( "Memory" );
		if ( old ){
			memcpy( data, old, m );
			delete []old;
		}
	}
	if( size ) data[size] = 0;
}
void Sequence::trim(int trim_len) {
    if (trim_len >= size) return;
    size = trim_len;
    if (size) data[size]=0;
}
void Sequence::ConvertBases()
{
	int i;
	for(i=0; i<size; i++) data[i] = aa2idx[data[i] - 'A'];
}

void Sequence::Swap( Sequence & other )
{
	Sequence tmp;
	memcpy( & tmp, this, sizeof( Sequence ) );
	memcpy( this, & other, sizeof( Sequence ) );
	memcpy( & other, & tmp, sizeof( Sequence ) );
	memset( & tmp, 0, sizeof( Sequence ) );
}
int Sequence::Format()
{
	int i, j=0, m = 0;
	while( size && isspace( data[size-1] ) ) size --;
	if( size && data[size-1] == '*' ) size --;
	if( size ) data[size] = 0;
	for (i=0; i<size; i++){
		char ch = data[i];
		m += ! (isalpha( ch ) | isspace( ch ));
	}
	if( m ) return m;
	for (i=0; i<size; i++){
		char ch = data[i];
		if ( isalpha( ch ) ) data[j++] = toupper( ch );
	}
	data[j] = 0;
	size = j;
	return 0;
}

void Sequence::SwapIn()
{
	if( data ) return;
	if( swap == NULL ) bomb_error( "Can not swap in sequence" );
	Resize( size );
	fseek( swap, offset, SEEK_SET );
	if( fread( data, 1, size, swap ) ==0 ) bomb_error( "Can not swap in sequence" );
	data[size] = 0;
}
void Sequence::SwapOut()
{
	if( swap && data ){
		delete[] data;
		bufsize = 0;
		data = NULL;
	}
}
void Sequence::PrintInfo( int id, FILE *fout, const Options & options, char *buf )
{
	const char *tag = options.isEST ? "nt" : "aa";
	bool print = options.print != 0;
	bool strand = options.isEST;
	fprintf( fout, "%i\t%i%s, >%s...", id, size, tag, identifier+1 );
	if( identity ){
		int *c = coverage;
		fprintf( fout, " at " );
		if (print) fprintf( fout, "%i:%i:%i:%i/", c[0], c[1], c[2], c[3] );
		if (strand) fprintf( fout, "%c/", (state & IS_MINUS_STRAND) ? '-' : '+' );
		fprintf( fout, "%.2f%%", identity*100 );
		if( options.useDistance ) fprintf( fout, "/%.2f%%", distance*100 );
		fprintf( fout, "\n" );
	}else{
		fprintf( fout, " *\n" );
	}
}

// by liwz
// disable swap option
// change des_begin, des_length, des_length2, dat_length => des_begin, tot_length
// where des_begin is the FILE pointer of sequence record start
//       tot_length is the total bytes of sequence record 
void SequenceDB::Read( const char *file, const Options & options )
{
    Sequence one;
    Sequence des;
    FILE *fin = fopen( file, "rb" );
    char *buffer = NULL;
    char *res = NULL;
    int option_l = options.min_length;
    if( fin == NULL ) bomb_error( "Failed to open the database file" );
    Clear();
    buffer = new char[ MAX_LINE_SIZE+1 ];

    while (not feof( fin ) || one.size) { /* do not break when the last sequence is not handled */
        buffer[0] = '>';
        if ( (res=fgets( buffer, MAX_LINE_SIZE, fin )) == NULL && one.size == 0) break;
        if( buffer[0] == '+' ){
            int len = strlen( buffer );
            int len2 = len;
            while( len2 && buffer[len2-1] != '\n' ){
                if ( (res=fgets( buffer, MAX_LINE_SIZE, fin )) == NULL ) break;
                len2 = strlen( buffer );
                len += len2;
            }
            one.tot_length += len;

            // read next line quality score
            if ( (res=fgets( buffer, MAX_LINE_SIZE, fin )) == NULL ) bomb_error("can not read quality score after");
            len = strlen( buffer );
            len2 = len;
            while( len2 && buffer[len2-1] != '\n' ){
                if ( (res=fgets( buffer, MAX_LINE_SIZE, fin )) == NULL ) break;
                len2 = strlen( buffer );
                len += len2;
            }
            one.tot_length += len;
        }else if (buffer[0] == '>' || buffer[0] == '@' || (res==NULL && one.size)) {
            if ( one.size ) { // write previous record
                if( one.identifier == NULL || one.Format() ){
                    printf( "Warning: from file \"%s\",\n", file );
                    printf( "Discarding invalid sequence or sequence without identifier and description!\n\n" );
                    if( one.identifier ) printf( "%s\n", one.identifier );
                    printf( "%s\n", one.data );
                    one.size = 0;
                }
                one.index = sequences.size();
                if( one.size > option_l ) {
                    if (options.trim_len    > 0) one.trim(options.trim_len);
                    sequences.Append( new Sequence( one ) ); 
                }
            }
            one.size = 0;
            one.tot_length = 0;

            int len = strlen( buffer );
            int len2 = len;
            des.size = 0;
            des += buffer;
            while( len2 && buffer[len2-1] != '\n' ){
                if ( (res=fgets( buffer, MAX_LINE_SIZE, fin )) == NULL ) break;
                des += buffer;
                len2 = strlen( buffer );
                len += len2;
            }
            size_t offset = ftell( fin );
            one.des_begin = offset - len;
            one.tot_length += len;              // count first line

            int i = 0;
            if( des.data[i] == '>' || des.data[i] == '@' || des.data[i] == '+' ) i += 1;
            if( des.data[i] == ' ' or des.data[i] == '\t' ) i += 1;
            if( options.des_len and options.des_len < des.size ) des.size = options.des_len;
            while( i < des.size and ! isspace( des.data[i] ) ) i += 1;
            des.data[i] = 0;
            one.identifier = des.data;
        } else {
            one.tot_length += strlen(buffer);  one += buffer;
        }
    }
#if 0
    int i, n = 0;
    for(i=0; i<sequences.size(); i++) n += sequences[i].bufsize + 4;
    cout<<n<<"\t"<<sequences.capacity() * sizeof(Sequence)<<endl;
    int i;
    scanf( "%i", & i );
#endif
    one.identifier = NULL;
    delete[] buffer;
    fclose( fin );
}

// PE reads liwz, disable swap option
void SequenceDB::Read( const char *file, const char *file2, const Options & options )
{
    Sequence one, two;
    Sequence des;
    FILE *fin = fopen( file, "rb" );
    FILE *fin2= fopen( file2,"rb" );
    char *buffer = NULL;
    char *buffer2= NULL;
    char *res = NULL;
    char *res2= NULL;
    int option_l = options.min_length;
    if( fin == NULL ) bomb_error( "Failed to open the database file" );
    if( fin2== NULL ) bomb_error( "Failed to open the database file" );
    Clear();
    buffer = new char[ MAX_LINE_SIZE+1 ];
    buffer2= new char[ MAX_LINE_SIZE+1 ];

    while (((not feof( fin )) && (not feof( fin2)) ) || (one.size && two.size)) { /* do not break when the last sequence is not handled */
        buffer[0] = '>'; res =fgets( buffer,  MAX_LINE_SIZE, fin  );
        buffer2[0]= '>'; res2=fgets( buffer2, MAX_LINE_SIZE, fin2 );

        if ( (res      == NULL) && (res2     != NULL)) bomb_error( "Paired input files have different number sequences" );
        if ( (res      != NULL) && (res2     == NULL)) bomb_error( "Paired input files have different number sequences" );
        if ( (one.size == 0   ) && (two.size >     0)) bomb_error( "Paired input files have different number sequences" );
        if ( (one.size >  0   ) && (two.size ==    0)) bomb_error( "Paired input files have different number sequences" );
        if ( (res      == NULL) && (one.size ==    0)) break;

        if( buffer[0] == '+' ){ // fastq 3rd line
            // file 1
            int len = strlen( buffer ); 
            int len2 = len;
            while( len2 && buffer[len2-1] != '\n' ){ // read until the end of the line
                if ( (res=fgets( buffer, MAX_LINE_SIZE, fin )) == NULL ) break;
                len2 = strlen( buffer );
                len += len2;
            }
            one.tot_length += len;

            // read next line quality score
            if ( (res=fgets( buffer, MAX_LINE_SIZE, fin )) == NULL ) bomb_error("can not read quality score after");
            len = strlen( buffer );
            len2 = len;
            while( len2 && buffer[len2-1] != '\n' ){
                if ( (res=fgets( buffer, MAX_LINE_SIZE, fin )) == NULL ) break;
                len2 = strlen( buffer );
                len += len2;
            }
            one.tot_length += len;

            // file 2
            len = strlen( buffer2 );
            len2 = len;
            while( len2 && buffer2[len2-1] != '\n' ){ // read until the end of the line
                if ( (res2=fgets( buffer2, MAX_LINE_SIZE, fin2 )) == NULL ) break;
                len2 = strlen( buffer2 );
                len += len2;
            }
            two.tot_length += len;

            // read next line quality score
            if ( (res2=fgets( buffer2, MAX_LINE_SIZE, fin2 )) == NULL ) bomb_error("can not read quality score after");
            len = strlen( buffer2 );
            len2 = len;
            while( len2 && buffer2[len2-1] != '\n' ){
                if ( (res2=fgets( buffer2, MAX_LINE_SIZE, fin2 )) == NULL ) break;
                len2 = strlen( buffer2 );
                len += len2;
            }
            two.tot_length += len;

        }else if (buffer[0] == '>' || buffer[0] == '@' || (res==NULL && one.size)) {
            if ( one.size && two.size ) { // write previous record
                if( one.identifier == NULL || one.Format() ){
                    printf( "Warning: from file \"%s\",\n", file );
                    printf( "Discarding invalid sequence or sequence without identifier and description!\n\n" );
                    if( one.identifier ) printf( "%s\n", one.identifier );
                    printf( "%s\n", one.data );
                    one.size=0; two.size=0;
                }
                if( two.identifier == NULL || two.Format() ){
                    printf( "Warning: from file \"%s\",\n", file2 );
                    printf( "Discarding invalid sequence or sequence without identifier and description!\n\n" );
                    if( two.identifier ) printf( "%s\n", two.identifier );
                    printf( "%s\n", two.data );
                    one.size=0; two.size = 0;
                }
                one.index = sequences.size();
                if( (one.size + two.size)> option_l ) {
                    if (options.trim_len    > 0) one.trim(options.trim_len);
                    if (options.trim_len_R2 > 0) two.trim(options.trim_len_R2);
                    sequences.Append( new Sequence( one, two, 1 ) ); 
                }
            }
            // R1
            one.size = 0;
            one.tot_length = 0;

            int len = strlen( buffer );
            int len2 = len;
            des.size = 0;
            des += buffer;
            while( len2 && buffer[len2-1] != '\n' ){
                if ( (res=fgets( buffer, MAX_LINE_SIZE, fin )) == NULL ) break;
                des += buffer;
                len2 = strlen( buffer );
                len += len2;
            }
            size_t offset = ftell( fin );    
            one.des_begin = offset - len; // offset of ">" or "@" 
            one.tot_length += len;              // count first line

            int i = 0;
            if( des.data[i] == '>' || des.data[i] == '@' || des.data[i] == '+' ) i += 1;
            if( des.data[i] == ' ' or des.data[i] == '\t' ) i += 1;
            if( options.des_len and options.des_len < des.size ) des.size = options.des_len;
            while( i < des.size and ! isspace( des.data[i] ) ) i += 1;
            des.data[i] = 0;                   // find first non-space letter
            one.identifier = des.data;

            // R2
            two.size = 0;
            two.tot_length = 0;

            len = strlen( buffer2 );
            len2 = len;
            while( len2 && buffer2[len2-1] != '\n' ){
                if ( (res=fgets( buffer2, MAX_LINE_SIZE, fin2 )) == NULL ) break;
                len2 = strlen( buffer2 );
                len += len2;
            }
            offset = ftell( fin2 );
            two.des_begin = offset - len;
            two.tot_length += len;              // count first line
            two.identifier = des.data;
        } else {
            one.tot_length += strlen(buffer);  one += buffer;
            two.tot_length+= strlen(buffer2); two+= buffer2;
        }
    }
#if 0
    int i, n = 0;
    for(i=0; i<sequences.size(); i++) n += sequences[i].bufsize + 4;
    cout<<n<<"\t"<<sequences.capacity() * sizeof(Sequence)<<endl;
    int i;
    scanf( "%i", & i );
#endif
    one.identifier = NULL;
    two.identifier = NULL;
    delete[] buffer;
    fclose( fin );
    delete[] buffer2;
    fclose( fin2 );
}

#if 0
void SequenceDB::Sort( int first, int last )
{
	int lower=first+1, upper=last;
	Sequence *pivot;
	Sequence *val;
	if( first >= last ) return;
	val = sequences[first];
	sequences[first] = sequences[ (first+last)/2 ];
	sequences[ (first+last)/2 ] = val;
	pivot = sequences[ first ];

	while( lower <= upper ){
		while( lower <= last && sequences[lower]->stats < pivot->stats ) lower ++;
		while( pivot->stats < sequences[upper]->stats ) upper --;
		if( lower < upper ){
			val = sequences[lower];
			sequences[lower] = sequences[upper];
			sequences[upper] = val;
			upper --;
		}
		lower ++;
	}
	val = sequences[first];
	sequences[first] = sequences[upper];
	sequences[upper] = val;
	if( first < upper-1 ) Sort( first, upper-1 );
	if( upper+1 < last ) Sort( upper+1, last );
}
#endif
void SequenceDB::SortDivide( Options & options, bool sort )
{
	int i, j, k, len;
	int N = sequences.size();
	total_letter=0;
	total_desc=0;
	max_len = 0;
	min_len = (size_t)-1;
	for (i=0; i<N; i++) {
		Sequence *seq = sequences[i];
		len = seq->size;
		total_letter += len;
		if (len > max_len) max_len = len;
		if (len < min_len) min_len = len;
		if (seq->swap == NULL) seq->ConvertBases();
		if( seq->identifier ) total_desc += strlen( seq->identifier );
	}
	options.max_entries = max_len * MAX_TABLE_SEQ;
	if (max_len >= 65536 and sizeof(INTs) <=2) 
		bomb_warning("Some seqs longer than 65536, you may define LONG_SEQ");

	if (max_len > MAX_SEQ ) 
		bomb_warning("Some seqs are too long, please rebuild the program with make parameter "
				"MAX_SEQ=new-maximum-length (e.g. make MAX_SEQ=10000000)");

	cout << "longest and shortest : " << max_len << " and " << min_len << endl;
	cout << "Total letters: " << total_letter << endl;
	// END change all the NR_seq to iseq

	len_n50 = (max_len + min_len) / 2; // will be properly set, if sort is true;
	if( sort ){
		// **************************** Form NR_idx[], Sort them from Long to short
		long long sum = 0;
		int M = max_len - min_len + 1;
		Vector<int> count( M, 0 ); // count for each size = max_len - i
		Vector<int> accum( M, 0 ); // count for all size > max_len - i
		Vector<int> offset( M, 0 ); // offset from accum[i] when filling sorting
		Vector<Sequence*> sorting( N ); // TODO: use a smaller class if this consumes to much memory!

		for (i=0; i<N; i++) count[ max_len - sequences[i]->size ] ++;
		for (i=1; i<M; i++) accum[i] = accum[i-1] + count[i-1];
		for (i=0; i<M; i++){
			sum += (max_len - i) * count[i];
			if( sum >= (total_letter>>1) ){
				len_n50 = max_len - i;
				break;
			}
		}
		for (i=0; i<N; i++){
			int len = max_len - sequences[i]->size;
			int id = accum[len] + offset[len];
			//sequences[i].index = id;
			sorting[id] = sequences[i];
			offset[len] ++;
		}
		options.max_entries = 0;
		for (i=0; i<N; i++){
			sequences[i] = sorting[i];
			if( i < MAX_TABLE_SEQ ) options.max_entries += sequences[i]->size;
		}
#if 0
		if( options.isEST ){
			int start = 0;
			for (i=0; i<M; i++){
				Sort( start, accum[i] );
				start = accum[i];
			}
		}
#endif
		cout << "Sequences have been sorted" << endl;
		// END sort them from long to short
	}
}// END sort_seqs_divide_segs

void SequenceDB::DivideSave( const char *db, const char *newdb, int n, const Options & options )
{
	if( n == 0 or sequences.size() ==0 ) return;

	size_t max_seg = total_letter / n + sequences[0]->size;
	if( max_seg >= MAX_BIN_SWAP ) max_seg = (size_t) MAX_BIN_SWAP;

	FILE *fin = fopen( db, "rb" );
	char *buf = new char[MAX_LINE_SIZE+1];
	char outfile[512];
	size_t seg_size = 0;
	int i, j, count, rest, seg = 0;
	sprintf( outfile, "%s-%i", newdb, 0 );
	FILE *fout = fopen( outfile, "w+" );
	n = sequences.size();
	for (i=0; i<n; i++){
		Sequence *seq = sequences[i];
		fseek( fin, seq->des_begin, SEEK_SET );

		seg_size += seq->size;
		if( seg_size >= max_seg ){
			seg += 1;
			sprintf( outfile, "%s-%i", newdb, seg );
			fclose( fout );
			fout = fopen( outfile, "w+" );
			seg_size = seq->size;
		}

		count = seq->tot_length / MAX_LINE_SIZE;
		rest  = seq->tot_length % MAX_LINE_SIZE;
		//printf( "count = %6i,  rest = %6i\n", count, rest );
		for (j=0; j<count; j++){
			if( fread( buf, 1, MAX_LINE_SIZE, fin ) ==0 ) bomb_error( "Can not swap in sequence" );
			fwrite( buf, 1, MAX_LINE_SIZE, fout );
		}
		if( rest ){
			if( fread( buf, 1, rest, fin ) ==0 ) bomb_error( "Can not swap in sequence" );
			fwrite( buf, 1, rest, fout );
		}
	}
	fclose( fin );
	fclose( fout );
	delete []buf;
}
void SequenceDB::WriteClusters( const char *db, const char *newdb, const Options & options )
{
	FILE *fin = fopen( db, "rb" );
	FILE *fout = fopen( newdb, "w+" );
	int i, j, n = rep_seqs.size();
	int count, rest;
	char *buf = new char[MAX_LINE_SIZE+1];
	vector<uint64_t> sorting( n );
	if( fin == NULL || fout == NULL ) bomb_error( "file opening failed" );
	for (i=0; i<n; i++) sorting[i] = ((uint64_t)sequences[ rep_seqs[i] ]->index << 32) | rep_seqs[i];
	std::sort( sorting.begin(), sorting.end() );
	for (i=0; i<n; i++){
		Sequence *seq = sequences[ sorting[i] & 0xffffffff ];
		fseek( fin, seq->des_begin, SEEK_SET );

		count = seq->tot_length / MAX_LINE_SIZE;
		rest  = seq->tot_length % MAX_LINE_SIZE;
		//printf( "count = %6i,  rest = %6i\n", count, rest );
		for (j=0; j<count; j++){
			if( fread( buf, 1, MAX_LINE_SIZE, fin ) ==0 ) bomb_error( "Can not swap in sequence" );
			fwrite( buf, 1, MAX_LINE_SIZE, fout );
		}
		if( rest ){
			if( fread( buf, 1, rest, fin ) ==0 ) bomb_error( "Can not swap in sequence" );
			fwrite( buf, 1, rest, fout );
		}
	}
	fclose( fin );
	fclose( fout );
	delete []buf;
}
// liwz PE output
void SequenceDB::WriteClusters( const char *db, const char *db_pe, const char *newdb, const char *newdb_pe, const Options & options )
{
	FILE *fin = fopen( db, "rb" );
	FILE *fout = fopen( newdb, "w+" );
	FILE *fin_pe = fopen( db_pe, "rb" );
	FILE *fout_pe = fopen( newdb_pe, "w+" );
	int i, j, n = rep_seqs.size();
	int count, rest;
	char *buf = new char[MAX_LINE_SIZE+1];
	vector<uint64_t> sorting( n );
	if( fin == NULL || fout == NULL ) bomb_error( "file opening failed" );
	if( fin_pe == NULL || fout_pe == NULL ) bomb_error( "file opening failed" );
	for (i=0; i<n; i++) sorting[i] = ((uint64_t)sequences[ rep_seqs[i] ]->index << 32) | rep_seqs[i];
	std::sort( sorting.begin(), sorting.end() );

        //sort fasta / fastq
        int *clstr_size;
        int *clstr_idx1;
        if (options.sort_outputf) {
            clstr_size = new int[n];
            clstr_idx1 = new int[n];
            for (i=0; i<n; i++) { 
                clstr_size[i] = 0;
                clstr_idx1[i]  = i;
            }

            int N = sequences.size();
            for (i=0; i<N; i++) { 
                int id = sequences[i]->cluster_id;
                if (id < 0) continue;
                if (id >=n) continue;
                clstr_size[id]++;
            }
            quick_sort_idxr(clstr_size, clstr_idx1, 0, n-1);
        }

	for (i=0; i<n; i++){
		Sequence *seq = sequences[ sorting[i] & 0xffffffff ];
                if (options.sort_outputf) seq = sequences[  rep_seqs[ clstr_idx1[i] ] ];
                //R1
		fseek( fin, seq->des_begin, SEEK_SET );

		count = seq->tot_length / MAX_LINE_SIZE;
		rest  = seq->tot_length % MAX_LINE_SIZE;
		//printf( "count = %6i,  rest = %6i\n", count, rest );
		for (j=0; j<count; j++){
			if( fread( buf, 1, MAX_LINE_SIZE, fin ) ==0 ) bomb_error( "Can not swap in sequence" );
			fwrite( buf, 1, MAX_LINE_SIZE, fout );
		}
		if( rest ){
			if( fread( buf, 1, rest, fin ) ==0 ) bomb_error( "Can not swap in sequence" );
			fwrite( buf, 1, rest, fout );
		}

                //R2
		fseek( fin_pe, seq->des_begin2, SEEK_SET );

		count = seq->tot_length2 / MAX_LINE_SIZE;
		rest  = seq->tot_length2 % MAX_LINE_SIZE;
		//printf( "count = %6i,  rest = %6i\n", count, rest );
		for (j=0; j<count; j++){
			if( fread( buf, 1, MAX_LINE_SIZE, fin_pe ) ==0 ) bomb_error( "Can not swap in sequence" );
			fwrite( buf, 1, MAX_LINE_SIZE, fout_pe );
		}
		if( rest ){
			if( fread( buf, 1, rest, fin_pe ) ==0 ) bomb_error( "Can not swap in sequence" );
			fwrite( buf, 1, rest, fout_pe );
		}

	}
	fclose( fin );
	fclose( fout );
	fclose( fin_pe );
	fclose( fout_pe );
	delete []buf;
}

void SequenceDB::WriteExtra1D( const Options & options )
{
	string db_clstr = options.output + ".clstr";
	string db_clstr_bak = options.output + ".bak.clstr";
	int i, i0, k, N = sequences.size();
	vector<long long> sorting( N );
	for (i=0; i<N; i++) sorting[i] = ((long long)sequences[i]->index << 32) | i;
	std::sort( sorting.begin(), sorting.end() );

	FILE *fout;
	char *buf = new char[ MAX_DES + 1 ];

	if( options.backupFile ){
		fout = fopen( db_clstr_bak.c_str(), "w+" );
		for (i=0; i<N; i++) {
			Sequence *seq = sequences[ sorting[i] & 0xffffffff ];
			seq->PrintInfo( seq->cluster_id, fout, options, buf );
		}
		fclose( fout );
	}

	cout << "writing clustering information" << endl;
	int M = rep_seqs.size();
	Vector<Vector<int> > clusters( M );
	for (i=0; i<N; i++){
		int k = sorting[i] & 0xffffffff;
		int id = sequences[k]->cluster_id;
		clusters[id].Append( k );
	}

	fout = fopen( db_clstr.c_str(), "w+" );

        if (options.sort_output) {
            int *clstr_size = new int[M];
            int *clstr_idx1 = new int[M];

            for (i=0; i<M; i++) { 
                clstr_size[i] = (int)clusters[i].size();
                clstr_idx1[i]  = i;
            }
            quick_sort_idxr(clstr_size, clstr_idx1, 0, M-1);

  	    for (i=0; i<M; i++) {
                i0 = clstr_idx1[i];
		fprintf( fout, ">Cluster %i\n", i );
		for (k=0; k<(int)clusters[i0].size(); k++)
			sequences[ clusters[i0][k] ]->PrintInfo( k, fout, options, buf );
	    }   
        }
        else {
  	    for (i=0; i<M; i++) {
		fprintf( fout, ">Cluster %i\n", i );
		for (k=0; k<(int)clusters[i].size(); k++)
			sequences[ clusters[i][k] ]->PrintInfo( k, fout, options, buf );
	    }   

        }

	delete []buf;
}
void SequenceDB::WriteExtra2D( SequenceDB & other, const Options & options )
{
	string db_clstr = options.output + ".clstr";
	string db_clstr_bak = options.output + ".bak.clstr";
	int i, k, N = other.sequences.size();
	int N2 = sequences.size();
	vector<long long> sorting( N );
	for (i=0; i<N; i++) sorting[i] = ((long long)other.sequences[i]->index << 32) | i;
	std::sort( sorting.begin(), sorting.end() );

	FILE *fout;
	char *buf = new char[ MAX_DES + 1 ];
	if( options.backupFile ){
		fout = fopen( db_clstr_bak.c_str(), "w+" );
		for (i=0; i<N; i++) {
			Sequence *seq = other.sequences[ sorting[i] & 0xffffffff ];
			seq->PrintInfo( seq->cluster_id, fout, options, buf );
		}
		for (i=0; i<N2; i++) {
			Sequence *seq = sequences[i];
			if( seq->state & IS_REDUNDANT ) seq->PrintInfo( seq->cluster_id, fout, options, buf );
		}
		fclose( fout );
	}

	cout << "writing clustering information" << endl;
	Vector<Vector<int> > clusters( N );
	for (i=0; i<N2; i++){
		int id = sequences[i]->cluster_id;
		if( sequences[i]->state & IS_REDUNDANT ) clusters[id].Append( i );
	}

	fout = fopen( db_clstr.c_str(), "w+" );
	for (i=0; i<N; i++) {
		Sequence *seq = other.sequences[ i ];
		fprintf( fout, ">Cluster %i\n", i );
		seq->PrintInfo( 0, fout, options, buf );
		for (k=0; k<(int)clusters[i].size(); k++)
			sequences[ clusters[i][k] ]->PrintInfo( k+1, fout, options, buf );
	}
	delete []buf;
}
void WorkingParam::ControlShortCoverage( int len, const Options & options )
{
	len_eff = len;
	aln_cover_flag = 0;
	if ((options.short_coverage > 0.0) || (options.min_control>0) ) { // has alignment coverage control
		aln_cover_flag = 1;
		min_aln_lenS = (int) (double(len) * options.short_coverage);
		if ( len-options.short_control > min_aln_lenS) min_aln_lenS = len-options.short_control;
		if ( options.min_control > min_aln_lenS) min_aln_lenS = options.min_control;
	}
	if (options.global_identity == 0) len_eff = min_aln_lenS; //global_identity==0
}
void WorkingParam::ControlLongCoverage( int len2, const Options & options )
{
	if (aln_cover_flag) {
		min_aln_lenL = (int) (double(len2) * options.long_coverage);
		if ( len2-options.long_control > min_aln_lenL) min_aln_lenL = len2-options.long_control;
		if ( options.min_control > min_aln_lenL) min_aln_lenL = options.min_control;
	}
}


// when alignment coverage such as -aL is specified
// if a existing rep is too long, it won't be qulified 
int upper_bound_length_rep(int len, double opt_s, int opt_S, double opt_aL, int opt_AL )
{
	int len_upper_bound = 99999999;
	double r1 = (opt_s > opt_aL) ? opt_s : opt_aL;
	int    a2 = (opt_S < opt_AL) ? opt_S : opt_AL;
	if (r1 > 0.0) len_upper_bound = (int) ( ((float) len)  / r1);
	if ((len+a2) < len_upper_bound)  len_upper_bound = len+a2;

	return len_upper_bound;
} // END upper_bound_length_rep
int upper_bound_length_rep(int len, const Options & options )
{
	double opt_s = options.diff_cutoff;
	int    opt_S = options.diff_cutoff_aa;
	double opt_aL = options.long_coverage;
	int    opt_AL = options.long_control;
	return upper_bound_length_rep( len, opt_s, opt_S, opt_aL, opt_AL );
}


void cal_aax_cutoff(double &aa1_cutoff, double &aa2_cutoff, double &aan_cutoff,
		double cluster_thd, int tolerance, int naa_stat_start_percent,
		int naa_stat[5][61][4], int NAA)
{
	aa1_cutoff = cluster_thd;
	aa2_cutoff = 1 - (1-cluster_thd)*2;
	aan_cutoff = 1 - (1-cluster_thd)*NAA;
	if (tolerance==0) return; 

	int clstr_idx = (int) (cluster_thd * 100) - naa_stat_start_percent;
	if (clstr_idx <0) clstr_idx = 0;
	double d2  = ((double) (naa_stat[tolerance-1][clstr_idx][3]     )) / 100;
	double dn  = ((double) (naa_stat[tolerance-1][clstr_idx][5-NAA] )) / 100;
	aa2_cutoff = d2 > aa2_cutoff ? d2 : aa2_cutoff;
	aan_cutoff = dn > aan_cutoff ? dn : aan_cutoff;
	return;
} // END cal_aax_cutoff


void update_aax_cutoff(double &aa1_cutoff, double &aa2_cutoff, double &aan_cutoff,
		int tolerance, int naa_stat_start_percent,
		int naa_stat[5][61][4], int NAA, double cluster_thd)
{
	if (cluster_thd > 1.0) cluster_thd = 1.00;

	double aa1_t, aa2_t, aan_t;
	cal_aax_cutoff(aa1_t, aa2_t, aan_t, cluster_thd, tolerance, naa_stat_start_percent,
			naa_stat, NAA);
	if (aa1_t > aa1_cutoff) aa1_cutoff = aa1_t;
	if (aa2_t > aa2_cutoff) aa2_cutoff = aa2_t;
	if (aan_t > aan_cutoff) aan_cutoff = aan_t;
	return;  
} // END update_aax_cutoff

void WorkingParam::ComputeRequiredBases( int NAA, int ss, const Options & option )
{
	// d: distance, fraction of errors;
	// e: number of errors;
	// g: length of the maximum gap;
	// m: word length;
	// n: sequence length;
	// alignment length = n - g + 1;
	// d = e / (n - g + 1);
	// e >= 1, so that, g <= n + 1 - 1/d
	// word count = (n - g - m + 1) - (e - 1)*m;
	//            = (n - g - m + 1) - (d*(n - g + 1) - 1)*m
	//            = (n - g + 1) - d*m*(n - g + 1)
	//            = (n - g + 1)*(1 - d*m)
	// minimum word count is reached when g == n + 1 - 1/d
	// so, minimum word count = 1/d - m.
	// if g == band_width: word count = (n - band + 1)*(1 - d*m);
	if( options.useDistance ){
		int band = options.band_width + 1;
		int invd = int( 1.0 / (options.distance_thd + 1E-9) );
		int k = len_eff < invd ? len_eff : invd;
		int ks = len_eff - ss + 1;
		int kn = len_eff - NAA + 1;
		int ks2 = invd - ss;
		int kn2= invd - NAA;
		int ks3 = int((len_eff - band + 1.0)*(1.0 - options.distance_thd * ss));
		int kn3 = int((len_eff - band + 1.0)*(1.0 - options.distance_thd * NAA));
		//if( ks3 > ks2 ) ks2 = ks3;
		//if( kn3 > kn2 ) kn2 = kn3;
		required_aa1 = required_aas = (ks2 < ks ? ks2 : ks);
		required_aan = kn2 < kn ? kn2 : kn;
		if( required_aa1 <=0 ) required_aa1 = required_aas = 1;
		if( required_aan <=0 ) required_aan = 1;
		//required_aa1 = required_aas = required_aan = 0;
		return;
	}
	// (N-K)-K*(1-C)*N = C*K*N-(K-1)*N-K = (C*K-K+1)*N-K
	required_aa1 = (len_eff - ss) - int(ss * ceil( (1.0 - aa1_cutoff) * len_eff ));
	if( required_aa1 < 0 ) required_aa1 = 0;
	required_aas = required_aa1;
	required_aan = (len_eff - NAA) - int(NAA * ceil( (1.0 - aa1_cutoff) * len_eff ));
	//printf( "%i %i\n", required_aa1, required_aan );
	if( required_aan < 0 ) required_aan = 0;

	int aa1_old = int (aa1_cutoff* (double) len_eff) - ss + 1;
	int aas_old = int (aas_cutoff* (double) len_eff);
	int aan_old = int (aan_cutoff* (double) len_eff);

	double thd = option.cluster_thd;
	//double rest = (len_eff - ss) / double(len_eff * ss);
	double rest = (len_eff - NAA) / double(len_eff * NAA);
	double thd0 = 1.0 - rest;
	double fnew = 0;
	double fold = 1;
	if( thd > thd0 ){
		fnew = (thd - thd0) / rest;
		fold = 1.0 - fnew;
	}
	//printf( "%g %g %g\n", thd, thd0, fnew );

	required_aa1 = (int)(fnew*required_aa1 + fold*aa1_old);
	required_aas = (int)(fnew*required_aas + fold*aas_old);
	required_aan = (int)(fnew*required_aan + fold*aan_old);
}
int WorkingBuffer::EncodeWords( Sequence *seq, int NAA, bool est )
{
	char *seqi = seq->data;
	int len = seq->size;
	// check_word_encodes
	int aan_no = len - NAA + 1;
	int i, j, i0, i1;
	int skip = 0;
	unsigned char k, k1;
	for (j=0; j<aan_no; j++) {
		char *word = seqi + j;
		int encode = 0;
		for (k=0, k1=NAA-1; k<NAA; k++, k1--) encode += word[k] * NAAN_array[k1];
		word_encodes[j] = word_encodes_backup[j] = encode;
	}

	if( est ){
		for (j=0; j<len; j++){
			if ( seqi[j] >= 4 ) {                      // here N is 4
				i0 = (j-NAA+1 > 0) ? j-NAA+1 : 0;
				i1 = j < aan_no ? j : aan_no - 1;
				for (i=i0; i<=i1; i++) word_encodes[i]=-1;
			}
		}
		for (j=0; j<aan_no; j++) skip += (word_encodes[j] == -1);
	}

	std::sort( word_encodes.begin(), word_encodes.begin() + aan_no );
	for(j=0; j<aan_no; j++) word_encodes_no[j]=1;
	for(j=aan_no-1; j; j--) {
		if (word_encodes[j] == word_encodes[j-1]) {
			word_encodes_no[j-1] += word_encodes_no[j];
			word_encodes_no[j]=0;
		}
	}
	return skip;
	// END check_word_encodes
}
void WorkingBuffer::ComputeAAP( const char *seqi, int size )
{
	int len1 = size - 1;
	int sk, j1, mm, c22;
	for (sk=0; sk<NAA2; sk++) taap[sk] = 0;
	for (j1=0; j1<len1; j1++) {
		c22= seqi[j1]*NAA1 + seqi[j1+1];
		taap[c22]++;
	}
	for (sk=0,mm=0; sk<NAA2; sk++) {
		aap_begin[sk] = mm; mm+=taap[sk]; taap[sk] = 0;
	}
	for (j1=0; j1<len1; j1++) {
		c22= seqi[j1]*NAA1 + seqi[j1+1];
		aap_list[aap_begin[c22]+taap[c22]++] =j1;
	}
}
void WorkingBuffer::ComputeAAP2( const char *seqi, int size )
{
	int len1 = size - 3;
	int sk, j1, mm, c22;
	for (sk=0; sk<NAA4; sk++) taap[sk] = 0;
	for (j1=0; j1<len1; j1++) {
		if ((seqi[j1]>=4) || (seqi[j1+1]>=4) || (seqi[j1+2]>=4) || (seqi[j1+3]>=4)) continue; //skip N
		c22 = seqi[j1]*NAA3 + seqi[j1+1]*NAA2 + seqi[j1+2]*NAA1 + seqi[j1+3];
		taap[c22]++;
	}
	for (sk=0,mm=0; sk<NAA4; sk++) {
		aap_begin[sk] = mm;  mm += taap[sk];  taap[sk] = 0;
	}
	for (j1=0; j1<len1; j1++) {
		if ((seqi[j1]>=4) || (seqi[j1+1]>=4) || (seqi[j1+2]>=4) || (seqi[j1+3]>=4)) continue; //skip N
		c22 = seqi[j1]*NAA3 + seqi[j1+1]*NAA2 + seqi[j1+2]*NAA1 + seqi[j1+3];
		aap_list[aap_begin[c22]+taap[c22]++] =j1;
	}
}

void SequenceDB::ClusterOne( Sequence *seq, int id, WordTable & table,
		WorkingParam & param, WorkingBuffer & buffer, const Options & options )
{
	if (seq->state & IS_REDUNDANT) return;
	int frag_size = options.frag_size;
	int NAA = options.NAA;
	int len = seq->size;
	int len_bound = upper_bound_length_rep(len, options);
	param.len_upper_bound = len_bound;
	int flag = CheckOne( seq, table, param, buffer, options );

	if( flag == 0 ){
		if ((seq->identity>0) && (options.cluster_best)) {
			// because of the -g option, this seq is similar to seqs in old SEGs
			seq->state |= IS_REDUNDANT ;
			seq->Clear();
		} else {                  // else add to NR90 db
			int aan_no = len - NAA + 1;
			int size = rep_seqs.size();
			rep_seqs.Append( id );
			seq->cluster_id = size;
			seq->identity = 0;
			seq->state |= IS_REP;
			if (frag_size){ /* not used for EST */
				int frg1 = (len - NAA ) / frag_size + 1;
				table.AddWordCountsFrag( aan_no, buffer.word_encodes_backup, 
						buffer.word_encodes_no, frg1, frag_size );
			}else{
				table.AddWordCounts(aan_no, buffer.word_encodes, buffer.word_encodes_no, table.sequences.size(), options.isEST);
			}
			table.sequences.Append( seq );
			if( frag_size ){
				while( table.sequences.size() < table.frag_count )
					table.sequences.Append( seq );
			}
		}
	}
	if ( (id+1) % 1000 == 0 ) {
		int size = rep_seqs.size();
		printf( "." );
		fflush( stdout );
		if ( (id+1) % 10000 == 0 ) printf( "\r..........%9i  finished  %9i  clusters\n", id+1, size );
	}
}
#include<assert.h>
size_t SequenceDB::MinimalMemory( int frag_no, int bsize, int T, const Options & options, size_t extra )
{
	int N = sequences.size();
	int F = frag_no < MAX_TABLE_SEQ ? frag_no : MAX_TABLE_SEQ;
	size_t mem_need = 0;
	size_t mem, mega = 1000000;
	int table = T > 1 ? 2 : 1;

	printf( "\nApproximated minimal memory consumption:\n" );
	mem = N*sizeof(Sequence) + total_desc + N + extra;
	if( options.store_disk == false ) mem += total_letter + N;
	printf( "%-16s: %zuM\n", "Sequence", mem/mega );
	mem_need += mem;

	mem = bsize;
	printf( "%-16s: %i X %zuM = %zuM\n", "Buffer", T, mem/mega, T*mem/mega );
	mem_need += T*mem;

	mem = F*(sizeof(Sequence*) + sizeof(IndexCount)) + NAAN*sizeof(NVector<IndexCount>);
	printf( "%-16s: %i X %zuM = %zuM\n", "Table", table, mem/mega, table*mem/mega );
	mem_need += table*mem;

	mem = sequences.capacity()*sizeof(Sequence*) + N*sizeof(int);
	mem += Comp_AAN_idx.size()*sizeof(int);
	printf( "%-16s: %zuM\n", "Miscellaneous", mem/mega );
	mem_need += mem;

	printf( "%-16s: %zuM\n\n", "Total", mem_need/mega );

	if(options.max_memory and options.max_memory < mem_need + 50*table ){
		char msg[200];
		sprintf( msg, "not enough memory, please set -M option greater than %zu\n", 
				50*table + mem_need/mega );
		bomb_error(msg);
	}
	return mem_need;
}
size_t MemoryLimit( size_t mem_need, const Options & options )
{
	size_t mem_limit = (options.max_memory - mem_need) / sizeof(IndexCount);

	//printf( "Table limit with the given memory limit:\n" );
	if( options.max_memory == 0 ){
		mem_limit = options.max_entries;
		if( mem_limit > MAX_TABLE_SIZE ) mem_limit = MAX_TABLE_SIZE;
	}
	//printf( "Max number of representatives: %zu\n", mem_limit );
	//printf( "Max number of word counting entries: %zu\n\n", mem_limit );
	return mem_limit;
}
void Options::ComputeTableLimits( int min_len, int max_len, int typical_len, size_t mem_need )
{
//liwz Fri Jan 15 15:44:47 PST 2016
//T=1 scale=1
//T=2 scale=0.6035
//T=4 scale=0.375
//T=8 scale=0.2392
//T=16 scale=0.1562
//T=32 scale=0.104
//T=64 scale=0.0703

	double scale = 0.5/threads + 0.5/sqrt(threads);
	max_sequences = (size_t)(scale * MAX_TABLE_SEQ);
	max_entries = (size_t)(scale * (500*max_len + 500000*typical_len + 50000000));
	if( max_memory ){
		double frac = max_sequences / (double) max_entries;
		max_entries = (options.max_memory - mem_need) / sizeof(IndexCount);
		max_sequences = (size_t)(max_entries * frac);
		if( max_sequences < MAX_TABLE_SEQ / 100 ) max_sequences = MAX_TABLE_SEQ / 100;
		if( max_sequences > MAX_TABLE_SEQ ) max_sequences = MAX_TABLE_SEQ;
	}
	printf( "Table limit with the given memory limit:\n" );
	printf( "Max number of representatives: %zu\n", max_sequences );
	printf( "Max number of word counting entries: %zu\n\n", max_entries );
}
void SequenceDB::DoClustering( int T, const Options & options )
{
	int i, j, k;
	int NAA = options.NAA;
	double aa1_cutoff = options.cluster_thd;
	double aas_cutoff = 1 - (1-options.cluster_thd)*4;
	double aan_cutoff = 1 - (1-options.cluster_thd)*options.NAA;
	int seq_no = sequences.size();
	int frag_no = seq_no;
	int frag_size = options.frag_size;
	int len, len_bound;
	int flag;
	valarray<size_t>  letters(T);

	//printf( "%li\n", options.mem_limit );

	if (frag_size){ 
		frag_no = 0;
		for (i=0; i<seq_no; i++) frag_no += (sequences[i]->size - NAA) / frag_size + 1;
	}

	if( not options.isEST )
		cal_aax_cutoff(aa1_cutoff, aas_cutoff, aan_cutoff, options.cluster_thd,
				options.tolerance, naa_stat_start_percent, naa_stat, NAA);

	Vector<WorkingParam> params(T);
	Vector<WorkingBuffer> buffers(T);
	for(i=0; i<T; i++){
		params[i].Set( aa1_cutoff, aas_cutoff, aan_cutoff );
		buffers[i].Set( frag_no, max_len, options );
	}

	// word_table as self comparing table and table buffer:
	WordTable word_table( options.NAA, NAAN );

	WordTable last_table( options.NAA, NAAN );

	int N = sequences.size();
	int K = N - 100 * T;
	size_t mem_need = MinimalMemory( frag_no, buffers[0].total_bytes, T, options );
	size_t mem_limit = MemoryLimit( mem_need, options );
	size_t mem, mega = 1000000;
	size_t tabsize = 0;
	int remaining = 0;

	Options opts( options );
	opts.ComputeTableLimits( min_len, max_len, len_n50, mem_need );

	omp_set_num_threads(T);
	for(i=0; i<N; ){
		int start = i;
		int m = i;
		size_t sum = remaining;
		float redundancy = (rep_seqs.size() + 1.0) / (i + 1.0);
		size_t max_items = opts.max_entries;
		size_t max_seqs = opts.max_sequences;
		size_t items = 0;
		if( i == 0 && max_seqs > 1000 ){ // first SCB with small size
			max_items /= 8;
			max_seqs /= 8;
		}
		while( m < N && (sum*redundancy) < max_seqs && items < max_items ){
			Sequence *seq = sequences[m];
			if( ! (seq->state & IS_REDUNDANT) ){
				if ( options.store_disk ) seq->SwapIn();
				//items += seq->size;
				items += (size_t)(seq->size * redundancy);
				sum += 1;
			}
			m ++;
		}
		if( (m > i + 1E4) && (m > i + (N - i) / (2+T)) ) m = i + (N - i) / (2+T);
		if( m == i || m >= N ){
			m = N;
			if( m > i + 1E3 ) m = i + (N - i) / (2+T);
		}
		//printf( "m = %i  %i,  %i\n", i, m, m-i );
		printf( "\r# comparing sequences from  %9i  to  %9i\n", i, m );
		if( last_table.size ){
			int print = (m-i)/20 + 1;
			#pragma omp parallel for schedule( dynamic, 1 )
			for(int j=i; j<m; j++){
				Sequence *seq = sequences[j];
				if (seq->state & IS_REDUNDANT) continue;
				int tid = omp_get_thread_num();
				CheckOne( seq, last_table, params[tid], buffers[tid], options );
				if ( options.store_disk && (seq->state & IS_REDUNDANT) ) seq->SwapOut();
				if( j%print==0 ){
					printf( "." ); fflush( stdout );
				}
			}
			int may_stop = 0;
			int self_stop = 0;
			float p0 = 0;
			int min = last_table.sequences[ last_table.sequences.size()-1 ]->size;
			int m0 = m;
			bool stop = false;
			#pragma omp parallel for schedule( dynamic, 1 )
			for(int j=m-1; j<N; j++){
				#pragma omp flush (stop)
				if( ! stop ){
					if( j+1 == N ) may_stop = 1;
					if( j == (m0-1) ){ // use m0 to avoid other iterations satisfying the condition:
						int tid = omp_get_thread_num();
						for(int ks=i; ks<m; ks++){
							Sequence *seq = sequences[ks];
							i = ks + 1;
							if (seq->state & IS_REDUNDANT) continue;
							ClusterOne( seq, ks, word_table, params[tid], buffers[tid], options );
							if ( options.store_disk && (seq->state & IS_REDUNDANT) ) seq->SwapOut();
							if( may_stop and word_table.sequences.size() >= 100 ) break;
							if( word_table.size >= max_items ) break;
							int tmax = max_seqs - (frag_size ? seq->size / frag_size + 1 : 0);
							if( word_table.sequences.size() >= tmax ) break;
						}
						self_stop = 1;
					}else{
						Sequence *seq = sequences[j];
						if (seq->state & IS_REDUNDANT) continue;
						if ( options.store_disk ){
							#pragma omp critical
							seq->SwapIn();
						}
						int tid = omp_get_thread_num();
						CheckOne( seq, last_table, params[tid], buffers[tid], options );
						if ( options.store_disk && (seq->state & IS_REDUNDANT) ) seq->SwapOut();
						if( min > params[tid].len_upper_bound ){
							may_stop = 1;
							stop = true;
							#pragma omp flush (stop)
						}
						if( self_stop && tid ==1 ){
							float p = (100.0*j)/N;
							if( p > p0+1E-1 ){ // print only if the percentage changed
								printf( "\r%4.1f%%", p );
								fflush( stdout );
								p0 = p;
							}
						}

					}
				}
			}
		}
		if( i == start || m == N ){
			//printf( "comparing the first or last or very small group ...\n" ); fflush( stdout );
			for(k=i; k<m; ){
				int kk, mm = k, sum = 0;
				while( mm < m && sum < 1E5 ){
					if( ! (sequences[mm]->state & IS_REDUNDANT) ) sum += sequences[mm]->size;
					mm += 1;
				}
				if( mm < k + 1000 ) mm = k + 1000;
				if( mm > m ) mm = m;
				#pragma omp parallel for schedule( dynamic, 1 )
				for(kk=k; kk<mm; kk++){
					Sequence *seq = sequences[kk];
					if (seq->state & IS_REDUNDANT) continue;
					int tid = omp_get_thread_num();
					CheckOne( seq, word_table, params[tid], buffers[tid], options );
					if ( options.store_disk && (seq->state & IS_REDUNDANT) ) seq->SwapOut();
				}
				bool bk = false;
				for(int ks=k; ks<mm; ks++){
					Sequence *seq = sequences[ks];
					i = k = ks + 1;
					if (seq->state & IS_REDUNDANT) continue;
					ClusterOne( seq, ks, word_table, params[0], buffers[0], options );
					bk = true;
					if ( options.store_disk && (seq->state & IS_REDUNDANT) ) seq->SwapOut();
					if( word_table.size >= max_items ) break;
					int tmax = max_seqs - (frag_size ? seq->size / frag_size + 1 : 0);
					if( word_table.sequences.size() >= tmax ) break;
					bk = false;
				}
				if( bk ) break;
			}
		}else if( i < m ){
			remaining = remaining/2 + (m - i);
			printf( "\r---------- %6i remaining sequences to the next cycle\n", m-i );
		}
		printf( "---------- new table with %8i representatives\n", word_table.sequences.size() );
		if( (last_table.size + word_table.size) > tabsize )
			tabsize = last_table.size + word_table.size;
		last_table.Clear();
		last_table.sequences.swap( word_table.sequences );
		last_table.indexCounts.swap( word_table.indexCounts );
		last_table.size = word_table.size;
		word_table.size = 0;
	}
	printf( "\n%9i  finished  %9i  clusters\n", sequences.size(), rep_seqs.size() );
	mem = (mem_need + tabsize*sizeof(IndexCount))/mega;
	printf( "\nApprixmated maximum memory consumption: %zuM\n", mem );
	last_table.Clear();
	word_table.Clear();
}

int SequenceDB::CheckOne( Sequence *seq, WordTable & table, WorkingParam & param, WorkingBuffer & buf, const Options & options )
{
	int len = seq->size;
	param.len_upper_bound = upper_bound_length_rep(len, options);
	if( options.isEST ) return CheckOneEST( seq, table, param, buf, options );
	return CheckOneAA( seq, table, param, buf, options );
}
int SequenceDB::CheckOneAA( Sequence *seq, WordTable & table, WorkingParam & param, WorkingBuffer & buf, const Options & options )
{
	NVector<IndexCount> & lookCounts = buf.lookCounts;
	NVector<uint32_t> & indexMapping = buf.indexMapping;
	Vector<INTs> & word_encodes_no = buf.word_encodes_no;
	Vector<INTs> & aap_list = buf.aap_list;
	Vector<INTs> & aap_begin = buf.aap_begin;
	Vector<int>  & word_encodes = buf.word_encodes;
	Vector<int>  & taap = buf.taap;
	double aa1_cutoff = param.aa1_cutoff;
	double aa2_cutoff = param.aas_cutoff;
	double aan_cutoff = param.aan_cutoff;

	char *seqi = seq->data;
	int j, k, j1, len = seq->size;
	int flag = 0;
	int frag_size = options.frag_size;
	int & aln_cover_flag = param.aln_cover_flag;
	int & required_aa1 = param.required_aa1;
	int & required_aa2 = param.required_aas;
	int & required_aan = param.required_aan;
	int & min_aln_lenS = param.min_aln_lenS;
	int & min_aln_lenL = param.min_aln_lenL;

	int NAA = options.NAA;
	int S = table.sequences.size();
	int len_eff = len;

	if( S ){
		int min = table.sequences[S-1]->size;
		if( min < len ){
			if( len * options.diff_cutoff2 > min ) min = (int)(len * options.diff_cutoff2);
			if( (len - options.diff_cutoff_aa2) > min ) min = len - options.diff_cutoff_aa2;
			len_eff = min;
		}
	}

    //liwz 2016 01, seq is too short for the shortest (longer) seq in word_table to satisfy -aL option
    //longer seqeunce * -aL -band_width
    if ( S ) {
		int min = table.sequences[S-1]->size;
		int min_red = min * options.long_coverage - options.band_width;
		if (len < min_red) return 0; // return flag=0
	} 

	param.ControlShortCoverage( len_eff, options );
	param.ComputeRequiredBases( options.NAA, 2, options );

	buf.EncodeWords( seq, options.NAA, false );

	// if minimal alignment length > len, return
	// I can not return earlier, because I need to calc the word_encodes etc
	if (options.min_control>len) return 0; // return flag=0

	// lookup_aan
	int aan_no = len - options.NAA + 1;
	int M = frag_size ? table.frag_count : S;
	table.CountWords(aan_no, word_encodes, word_encodes_no, lookCounts, indexMapping, false, required_aan);

	// contained_in_old_lib()
	int len_upper_bound = param.len_upper_bound;
	int len_lower_bound = param.len_lower_bound;
	int band_left, band_right, best_score, band_width1, best_sum, len2, alnln, len_eff1;
	int tiden_no, band_center;
	float tiden_pc, distance=0;
	int talign_info[5];
	int best1, sum;
	INTs *lookptr;
	char *seqj;
	int frg2 = frag_size ? (len - NAA + options.band_width ) / frag_size + 1 + 1 : 0;
	int lens;
	int has_aa2 = 0;

	IndexCount *ic = lookCounts.items;
	ic = lookCounts.items;
	for(; ic->count; ic++){
		if( ! frag_size ){
			indexMapping[ ic->index ] = 0;
			if ( ic->count < required_aan ) continue;
		}

		Sequence *rep = table.sequences[ ic->index ];
		len2 = rep->size;
		if (len2 > len_upper_bound ) continue;
		if (options.has2D && len2 < len_lower_bound ) continue;
		if( frag_size ){
			uint32_t *ims = & indexMapping[ ic->index ];
			int count = ic->count;
			k = (len2 - NAA) / frag_size + 1;
			sum = 0;
			for (j1=0; j1<frg2; j1++){
				uint32_t im = ims[j1];
				if( im ) sum += lookCounts[im-1].count;
			}
			count = sum;
			for (j1=frg2; j1<k; j1++) {
				uint32_t im1 = ims[j1];
				uint32_t im2 = ims[j1-frg2];
				if( im1 ) sum += lookCounts[im1-1].count;
				if( im2 ) sum -= lookCounts[im2-1].count;
				if (sum > count) count = sum;
			}
			if ( count < required_aan ) continue;
		}

		param.ControlLongCoverage( len2, options );

		if ( has_aa2 == 0 )  { // calculate AAP array
			buf.ComputeAAP( seqi, seq->size );
			has_aa2 = 1;
		}
		seqj = rep->data; //NR_seq[NR90_idx[j]];

		band_width1 = (options.band_width < len+len2-2 ) ? options.band_width : len+len2-2;
		diag_test_aapn(NAA1, seqj, len, len2, buf, best_sum,
				band_width1, band_left, band_center, band_right, required_aa1);
		if ( best_sum < required_aa2 ) continue;

		int rc = FAILED_FUNC;
		if (options.print || aln_cover_flag) //return overlap region
			rc = local_band_align(seqi, seqj, len, len2, mat,
					best_score, tiden_no, alnln, distance, talign_info,
					band_left, band_center, band_right, buf);
		else
			rc = local_band_align(seqi, seqj, len, len2, mat,
					best_score, tiden_no, alnln, distance, talign_info, 
					band_left, band_center, band_right, buf);
		if ( rc == FAILED_FUNC ) continue;
		if ( tiden_no < required_aa1 ) continue;
		lens = len;
		if( options.has2D && len > len2 ) lens = len2;
		len_eff1 = (options.global_identity == 0) ? alnln : (lens - talign_info[4]);
		tiden_pc = tiden_no / (float) len_eff1;
		if( options.useDistance ){
			if (distance > options.distance_thd ) continue;
			if (distance >= seq->distance) continue; // existing distance
		}else{
			if (tiden_pc < options.cluster_thd) continue;
			if (tiden_pc <= seq->identity) continue; // existing iden_no
		}
		if (aln_cover_flag) {
			if ( talign_info[3]-talign_info[2]+1 < min_aln_lenL) continue;
			if ( talign_info[1]-talign_info[0]+1 < min_aln_lenS) continue;
		}
		if( options.has2D ) seq->state |= IS_REDUNDANT ;
		flag = 1; seq->identity = tiden_pc; seq->cluster_id = rep->cluster_id;
		seq->distance = distance;
		seq->coverage[0] = talign_info[0] +1;
		seq->coverage[1] = talign_info[1] +1;
		seq->coverage[2] = talign_info[2] +1;
		seq->coverage[3] = talign_info[3] +1;
		if (not options.cluster_best) break;
		update_aax_cutoff(aa1_cutoff, aa2_cutoff, aan_cutoff,
				options.tolerance, naa_stat_start_percent, naa_stat, NAA, tiden_pc);
		param.ComputeRequiredBases( options.NAA, 2, options );
	}
	if( frag_size ) ic = lookCounts.items;
	while( ic->count ){
		indexMapping[ ic->index ] = 0;
		ic += 1;
	}
	lookCounts.size = 0;
	if (flag == 1) { // if similar to old one delete it
		if (! options.cluster_best) {
			seq->Clear();
			seq->state |= IS_REDUNDANT ;
		}
	}
	return flag;

}
int SequenceDB::CheckOneEST( Sequence *seq, WordTable & table, WorkingParam & param, WorkingBuffer & buf, const Options & options )
{
	NVector<IndexCount> & lookCounts = buf.lookCounts;
	NVector<uint32_t> & indexMapping = buf.indexMapping;
	Vector<INTs> & word_encodes_no = buf.word_encodes_no;
	Vector<INTs> & aap_list = buf.aap_list;
	Vector<INTs> & aap_begin = buf.aap_begin;
	Vector<int>  & word_encodes = buf.word_encodes;
	Vector<int>  & taap = buf.taap;
	Vector<int> & aan_list_comp = buf.aan_list_comp;
	char *seqi_comp = & buf.seqi_comp[0];

	int & aln_cover_flag = param.aln_cover_flag;
	int & required_aa1 = param.required_aa1;
	int & required_aas = param.required_aas;
	int & required_aan = param.required_aan;
	int & min_aln_lenS = param.min_aln_lenS;
	int & min_aln_lenL = param.min_aln_lenL;

	char *seqi = seq->data;
	int j, len = seq->size;
	int flag = 0;
	int S = table.sequences.size();
	int len_eff = len;
	if( S ){
		int min = table.sequences[S-1]->size;
		if( min < len ){
			if( len * options.diff_cutoff2 > min ) min = (int)(len * options.diff_cutoff2);
			if( (len - options.diff_cutoff_aa2) > min ) min = len - options.diff_cutoff_aa2;
			len_eff = min;
		}
	}


        //liwz 2016 01, seq is too short for the shortest (longer) seq in word_table to satisfy -aL option
        //longer seqeunce * -aL -band_width
        if ( S ) {
		int min = table.sequences[S-1]->size;
		int min_red = min * options.long_coverage - options.band_width;
		if (len < min_red) return 0; // return flag=0
	} 


	param.ControlShortCoverage( len_eff, options );
	param.ComputeRequiredBases( options.NAA, 4, options );
	int skip = buf.EncodeWords( seq, options.NAA, true );
	required_aan -= skip;
	required_aas -= skip;
	required_aa1 -= skip;
	if( required_aan <= 0 ) required_aan = 1;
	if( required_aas <= 0 ) required_aas = 1;
	if( required_aa1 <= 0 ) required_aa1 = 1;

	// if minimal alignment length > len, return
	// I can not return earlier, because I need to calc the word_encodes etc
	if (options.min_control>len) return 0; // return flag=0

	int aan_no = len - options.NAA + 1;

	// contained_in_old_lib()
	int len_upper_bound = param.len_upper_bound;
	int len_lower_bound = param.len_lower_bound;
	int band_left, band_right, best_score, band_width1, best_sum, len2, alnln, len_eff1;
	int tiden_no, band_center;
	float tiden_pc, distance=0;
	int talign_info[5];
	int j0, comp, lens;
	char *seqj;

	for(comp=0; comp<2; comp++){
		if( comp ){
			for (j0=0; j0<aan_no; j0++) {
				j = word_encodes[j0];
				if ( j<0 ) aan_list_comp[j0] = j;
				else       aan_list_comp[j0] = Comp_AAN_idx[j];
			}
			make_comp_iseq(len, seqi_comp, seqi);
			seqi = seqi_comp;
		}
		int has_aas = 0;

		if( comp ){
			table.CountWords(aan_no, aan_list_comp, word_encodes_no, lookCounts, indexMapping, true, required_aan );
		}else{
			table.CountWords(aan_no, word_encodes, word_encodes_no, lookCounts, indexMapping, true, required_aan ); 
		}

		IndexCount *ic = lookCounts.items;
		ic = lookCounts.items;
		for(; ic->count; ic++){
			indexMapping[ ic->index ] = 0;
			if ( ic->count < required_aan ) continue;
			Sequence *rep = table.sequences[ic->index];

			len2 = rep->size;
			if (len2 > len_upper_bound ) continue;
			if (options.has2D && len2 < len_lower_bound ) continue;

			seqj = rep->data;

			param.ControlLongCoverage( len2, options );

			if ( has_aas == 0 )  { // calculate AAP array
				buf.ComputeAAP2( seqi, seq->size );
				has_aas = 1;
			}

			band_width1 = (options.band_width < len+len2-2 ) ? options.band_width : len+len2-2;
			diag_test_aapn_est(NAA1, seqj, len, len2, buf, best_sum,
					band_width1, band_left, band_center, band_right, required_aa1);
			if ( best_sum < required_aas ) continue;
			//if( comp and flag and (not options.cluster_best) and j > rep->cluster_id ) goto Break;

			int rc = FAILED_FUNC;
			if (options.print || aln_cover_flag){ //return overlap region
				rc = local_band_align(seqi, seqj, len, len2, mat,
						best_score, tiden_no, alnln, distance, talign_info,
						band_left, band_center, band_right, buf);
				if( comp ){
					talign_info[0] = len - talign_info[0] - 1;
					talign_info[1] = len - talign_info[1] - 1;
				}
			}else{
				//printf( "%5i %5i %5i %5i\n", band_width1, band_right-band_left, band_left, band_right );
				rc = local_band_align(seqi, seqj, len, len2, mat,
						best_score, tiden_no, alnln, distance, talign_info,
						band_left, band_center, band_right, buf);
			}
			if ( rc == FAILED_FUNC ) continue;
			//printf( "%i  %i  %i\n", best_score, tiden_no, required_aa1 );
			if ( tiden_no < required_aa1 ) continue;
			if ( options.is454 ){
				if (talign_info[2] != talign_info[0]) continue; // same start
				if (talign_info[0] > 1) continue; // one mismatch allowed at beginning
				if ((len-talign_info[1]) > 2) continue; // one mismatch allowed at end
			}

			lens = len;
			if( options.has2D && len > len2 ) lens = len2;
			len_eff1 = (options.global_identity == 0) ? alnln : (lens - talign_info[4]);
			tiden_pc = tiden_no / (float)len_eff1;
			//printf( "%i %f\n", tiden_no, tiden_pc );
			if( options.useDistance ){
				if (distance > options.distance_thd ) continue;
				if (options.cluster_best and distance >= seq->distance) continue; // existing distance
			}else{
				if (tiden_pc < options.cluster_thd) continue;
				if (options.cluster_best and tiden_pc < seq->identity) continue; // existing iden_no
			}
			if (aln_cover_flag) {
				if ( talign_info[3]-talign_info[2]+1 < min_aln_lenL) continue;
				if( comp ){
					if ( talign_info[0]-talign_info[1]+1 < min_aln_lenS) continue;
				}else{
					if ( talign_info[1]-talign_info[0]+1 < min_aln_lenS) continue;
				}
			}
			if( options.cluster_best and fabs(tiden_pc - seq->identity) < 1E-9 and rep->cluster_id >= seq->cluster_id ) continue;
			if( (not options.cluster_best) and flag !=0 and rep->cluster_id >= seq->cluster_id ) continue;
			flag = comp ? -1 : 1;
			seq->identity = tiden_pc;
			seq->distance = distance;
			seq->cluster_id = rep->cluster_id;
			seq->coverage[0] = talign_info[0] +1;
			seq->coverage[1] = talign_info[1] +1;
			seq->coverage[2] = talign_info[2] +1;
			seq->coverage[3] = talign_info[3] +1;
			if (not options.cluster_best) break;
		}
		while( ic->count ){
			indexMapping[ ic->index ] = 0;
			ic += 1;
		}
		lookCounts.size = 0;
		if (not options.option_r ) break;
	}
	if ((flag == 1) || (flag == -1)) { // if similar to old one delete it
		if (! options.cluster_best) {
			seq->Clear();
			seq->state |= IS_REDUNDANT ;
		}
		if( flag == -1 )
			seq->state |= IS_MINUS_STRAND; 
		else
			seq->state &= ~IS_MINUS_STRAND; 
	}
	return flag;
}
void SequenceDB::ComputeDistance( const Options & options )
{
	int i, j, N = sequences.size();
	int best_score, best_sum;
	int band_width1, band_left, band_center, band_right, required_aa1;
	int tiden_no, alnln;
	int talign_info[5];
	float distance;
	WorkingBuffer buf( N, max_len, options );

	Vector<NVector<float> > dists( N, NVector<float>(N) );

	Sequence comseq( *sequences[0] );

	for(i=0; i<N; i++){
		Sequence *seq = sequences[i];
		char *seqi = seq->data;
		int len = seq->size;
		buf.EncodeWords( seq, options.NAA, false );
		buf.ComputeAAP2( seqi, seq->size );
		dists[i][i] = 0.0;
		if((i+1)%1000 ==0) printf( "%9i\n", (i+1) );
		for(j=0; j<i; j++){
			Sequence *rep = sequences[j];
			char *seqj = rep->data;
			int len2 = rep->size;
			band_width1 = (options.band_width < len+len2-2 ) ? options.band_width : len+len2-2;
			diag_test_aapn_est(NAA1, seqj, len, len2, buf, best_sum,
					band_width1, band_left, band_center, band_right, 0);
			local_band_align(seqi, seqj, len, len2, mat,
					best_score, tiden_no, alnln, distance, talign_info,
					band_left, band_center, band_right, buf);
			dists[seq->index][rep->index] = dists[rep->index][seq->index] = distance;
		}
		if (not options.option_r ) break;
		comseq.index = seq->index;
		comseq.size = len;
		for(j=0; j<len; j++) comseq.data[i] = seq->data[len-i-1];
		seqi = comseq.data;
		buf.EncodeWords( &comseq, options.NAA, false );
		buf.ComputeAAP2( seqi, seq->size );
		for(j=0; j<i; j++){
			Sequence *rep = sequences[j];
			char *seqj = rep->data;
			int len2 = rep->size;
			band_width1 = (options.band_width < len+len2-2 ) ? options.band_width : len+len2-2;
			diag_test_aapn_est(NAA1, seqj, len, len2, buf, best_sum,
					band_width1, band_left, band_center, band_right, 0);
			local_band_align(seqi, seqj, len, len2, mat,
					best_score, tiden_no, alnln, distance, talign_info,
					band_left, band_center, band_right, buf);
			if( distance < dists[seq->index][rep->index] )
				dists[seq->index][rep->index] = dists[rep->index][seq->index] = distance;
		}
	}
	std::string output = options.output + ".dist";
	FILE *fout = fopen( output.c_str(), "w+" );
	fprintf( fout, "1" );
	for(i=1; i<N; i++) fprintf( fout, "\t%i", i+1 );
	fprintf( fout, "\n" );
	for(i=0; i<N; i++){
		fprintf( fout, "%g", dists[i][0] );
		for(j=1; j<N; j++) fprintf( fout, "\t%g", dists[i][j] );
		fprintf( fout, "\n" );
	}
	fclose( fout );
}
void SequenceDB::DoClustering( const Options & options )
{
	int i;
	int NAA = options.NAA;
	double aa1_cutoff = options.cluster_thd;
	double aas_cutoff = 1 - (1-options.cluster_thd)*4;
	double aan_cutoff = 1 - (1-options.cluster_thd)*options.NAA;
	int seq_no = sequences.size();
	int frag_no = seq_no;
	int frag_size = options.frag_size;
	int len, len_bound;
	int flag;

#if 0
	ComputeDistance( options );
	return;
#endif

	if( options.threads > 1 ){
		DoClustering( options.threads, options );
		temp_files.Clear();
		return;
	}

	if (frag_size){ 
		frag_no = 0;
		for (i=0; i<seq_no; i++) frag_no += (sequences[i]->size - NAA) / frag_size + 1;
	}

	if( not options.isEST )
		cal_aax_cutoff(aa1_cutoff, aas_cutoff, aan_cutoff, options.cluster_thd,
				options.tolerance, naa_stat_start_percent, naa_stat, NAA);

	WorkingParam param( aa1_cutoff, aas_cutoff, aan_cutoff );
	WorkingBuffer buffer( frag_no, max_len, options );

	WordTable word_table( options.NAA, NAAN );

	size_t mem_need = MinimalMemory( frag_no, buffer.total_bytes, 1, options );
	size_t mem_limit = MemoryLimit( mem_need, options );
	size_t mem, mega = 1000000;
	int N = sequences.size();

	size_t total_letters = total_letter;
	size_t tabsize = 0;

	Options opts( options );
	opts.ComputeTableLimits( min_len, max_len, len_n50, mem_need );

	for(i=0; i<N; ){
		float redundancy = (rep_seqs.size() + 1.0) / (i + 1.0);
		int m = i;
		size_t sum = 0;
		size_t max_items = opts.max_entries;
		size_t max_seqs = opts.max_sequences;
		size_t items = 0;

// find a block from i to m, so that this block can fit into a word table
//     ...
//  i  ++++++++++++++++++++++++++
//     ++++++++++++++++++++
//     ++++++++++++++++
//  m  +++++++++++++
//     ...
		while( m < N && (sum*redundancy) < max_seqs && items < max_items ){
			Sequence *seq = sequences[m];
			if( ! (seq->state & IS_REDUNDANT) ){
				if ( options.store_disk ) seq->SwapIn();
				items += (size_t)(seq->size * redundancy);
				sum += 1;
			}
			m ++;
		}
		if( m > N ) m = N;
		printf( "\rcomparing sequences from  %9i  to  %9i\n", i, m );
		fflush( stdout );
		for(int ks=i; ks<m; ks++){ // clustering this block
			Sequence *seq = sequences[ks];
			i = ks + 1;
			if (seq->state & IS_REDUNDANT) continue;
			ClusterOne( seq, ks, word_table, param, buffer, options );
			total_letters -= seq->size;
			if( options.store_disk && (seq->state & IS_REDUNDANT) ) seq->SwapOut();
			if( word_table.size >= max_items ) break;
			int tmax = max_seqs - (frag_size ? seq->size / frag_size + 1 : 0);
			if( word_table.sequences.size() >= tmax ) break;
		} // finishing word table from this block
		m = i;
		if( word_table.size == 0 ) continue;
		float p0 = 0;
		for(int j=m; j<N; j++){ // use this word table to screen rest sequences m->N
			Sequence *seq = sequences[j];
			if (seq->state & IS_REDUNDANT) continue;
			if ( options.store_disk ) seq->SwapIn();
			CheckOne( seq, word_table, param, buffer, options );
			total_letters -= seq->size;
			if ( options.store_disk && (seq->state & IS_REDUNDANT) ) seq->SwapOut();
			int len_bound = param.len_upper_bound;
			if( word_table.sequences[ word_table.sequences.size()-1 ]->size > len_bound ){
				break;
			}
			float p = (100.0*j)/N;
			if( p > p0+1E-1 ){ // print only if the percentage changed
				printf( "\r%4.1f%%", p );
				fflush( stdout );
				p0 = p;
			}
		}
		if( word_table.size > tabsize ) tabsize = word_table.size;
		//if( i && i < m ) printf( "\r---------- %6i remaining sequences to the next cycle\n", m-i );
		word_table.Clear();
	}
	printf( "\n%9i  finished  %9i  clusters\n", sequences.size(), rep_seqs.size() );
	mem = (mem_need + tabsize*sizeof(IndexCount))/mega;
	printf( "\nApprixmated maximum memory consumption: %liM\n", mem );
	temp_files.Clear();
	word_table.Clear();

#if 0
	int zeros = 0;
	for(i=0; i<word_table.indexCounts.size(); i++) zeros += word_table.indexCounts[i].Size() ==0;
	printf( "%9i  empty entries out of  %9i\n", zeros, word_table.indexCounts.size() );
#endif
}

void SequenceDB::ClusterTo( SequenceDB & other, const Options & options )
{
	int i, is, js, flag;
	int len, len_tmp, len_lower_bound, len_upper_bound;
	int NR2_red_no = 0;
	int aan_no = 0;
	char *seqi;
	int NAA = options.NAA;
	double aa1_cutoff = options.cluster_thd;
	double aas_cutoff = 1 - (1-options.cluster_thd)*4;
	double aan_cutoff = 1 - (1-options.cluster_thd)*options.NAA;
	Vector<int>  word_encodes( MAX_SEQ );
	Vector<INTs> word_encodes_no( MAX_SEQ );

	if( not options.isEST ){
		cal_aax_cutoff(aa1_cutoff, aas_cutoff, aan_cutoff, options.cluster_thd,
				options.tolerance, naa_stat_start_percent, naa_stat, NAA);
	}


	int N = other.sequences.size();
	int M = sequences.size();
	int T = options.threads;

	valarray<size_t>  counts(T);
	Vector<WorkingParam> params(T);
	Vector<WorkingBuffer> buffers(T);
	WorkingParam & param = params[0];
	WorkingBuffer & buffer = buffers[0];
	for(i=0; i<T; i++){
		params[i].Set( aa1_cutoff, aas_cutoff, aan_cutoff );
		buffers[i].Set( N, max_len, options );
	}
	if( T >1 ) omp_set_num_threads(T);

	size_t mem_need = MinimalMemory( N, buffer.total_bytes, T, options, other.total_letter+other.total_desc );
	size_t mem_limit = MemoryLimit( mem_need, options );

	Options opts( options );
	opts.ComputeTableLimits( min_len, max_len, len_n50, mem_need );

	WordTable word_table( options.NAA, NAAN );

	size_t max_items = opts.max_entries;
	size_t max_seqs = opts.max_sequences;
	for(i=0; i<N; ){
		size_t items = 0;
		size_t sum = 0;
		int m = i;
		while( m < N && sum < max_seqs && items < max_items ){
			Sequence *seq = other.sequences[m];
			if( ! (seq->state & IS_REDUNDANT) ){
				if ( options.store_disk ) seq->SwapIn();
				items += seq->size;
				sum += 1;
			}
			m ++;
		}
		if( m > N ) m = N;
		//printf( "m = %i  %i,  %i\n", i, m, m-i );
		for(int ks=i; ks<m; ks++){
			Sequence *seq = other.sequences[ks];
			len = seq->size;
			seqi = seq->data;
			calc_ann_list(len, seqi, NAA, aan_no, word_encodes, word_encodes_no, options.isEST);
			word_table.AddWordCounts(aan_no, word_encodes, word_encodes_no, ks-i, options.isEST);
			word_table.sequences.Append( seq );
			seq->cluster_id = ks;
			seq->state |= IS_REP;
			if ( (ks+1) % 1000 == 0 ) {
				printf( "." );
				fflush( stdout );
				if ( (ks+1) % 10000 == 0 ) printf( "%9i  finished\n", ks+1 );
			}  
		}
		float p0 = 0;
		if( T > 1 ){
			int JM = M;
			counts = 0;
			#pragma omp parallel for schedule( dynamic, 1 )
			for(int j=0; j<JM; j++){
				Sequence *seq = sequences[j];
				if( seq->state & IS_REDUNDANT ) continue;
				int len = seq->size;
				char *seqi = seq->data;
				int len_upper_bound = upper_bound_length_rep(len,options);
				int len_lower_bound = len - options.diff_cutoff_aa2;

				int len_tmp = (int) ( ((double)len) * options.diff_cutoff2);
				if (len_tmp < len_lower_bound) len_lower_bound = len_tmp;

				int tid = omp_get_thread_num();
				params[tid].len_upper_bound = len_upper_bound;
				params[tid].len_lower_bound = len_lower_bound;

				if( word_table.sequences[ word_table.sequences.size()-1 ]->size > len_upper_bound ){
					JM = 0;
					continue;
				}

				int flag = other.CheckOne( seq, word_table, params[tid], buffers[tid], options );
				if ((flag == 1) || (flag == -1)) { // if similar to old one delete it
					if (! options.cluster_best) {
						seq->Clear();
						seq->state |= IS_REDUNDANT ;
						counts[tid] ++;
					}
					if( flag == -1 ) seq->state |= IS_MINUS_STRAND; // for EST only
				}
				float p = (100.0*j)/M;
				if( p > p0+1E-1 ){ // print only if the percentage changed
					printf( "\r%4.1f%%", p );
					fflush( stdout );
					p0 = p;
				}
			}
			for(int j=0; j<T; j++) NR2_red_no += counts[j];
		}else{
			for(int j=0; j<M; j++){
				Sequence *seq = sequences[j];
				if( seq->state & IS_REDUNDANT ) continue;
				len = seq->size;
				seqi = seq->data;
				len_upper_bound = upper_bound_length_rep(len,options);
				len_lower_bound = len - options.diff_cutoff_aa2;

				len_tmp = (int) ( ((double)len) * options.diff_cutoff2);
				if (len_tmp < len_lower_bound) len_lower_bound = len_tmp;
				param.len_upper_bound = len_upper_bound;
				param.len_lower_bound = len_lower_bound;

				if( word_table.sequences[ word_table.sequences.size()-1 ]->size > len_upper_bound ){
					break;
				}

				flag = other.CheckOne( seq, word_table, param, buffer, options );
				if ((flag == 1) || (flag == -1)) { // if similar to old one delete it
					if (! options.cluster_best) {
						seq->Clear();
						seq->state |= IS_REDUNDANT ;
						NR2_red_no ++;
					}
					if( flag == -1 ) seq->state |= IS_MINUS_STRAND; // for EST only
				}
				float p = (100.0*j)/M;
				if( p > p0+1E-1 ){ // print only if the percentage changed
					printf( "\r%4.1f%%", p );
					fflush( stdout );
					p0 = p;
				}
			}
		}
		printf( "\r..........%9i  compared  %9i  clusters\n", i, NR2_red_no );
		word_table.Clear();
		word_table.size = 0;
		i = m;
	}

	if (options.cluster_best) {//delete redundant sequences in options.cluster_best mode
		for (i=0; i<(int)sequences.size(); i++){
			Sequence *seq = sequences[i];
			if (seq->identity > 0 ){
				seq->state |= IS_REDUNDANT;
				NR2_red_no ++;
			}
		}
	}
	for (i=0; i<(int)sequences.size(); i++){
		Sequence *seq = sequences[i];
		if( seq->identity <0 ) seq->identity *= -1;
		if( not(seq->state & IS_REDUNDANT) ) rep_seqs.Append( i );
	}

	cout << endl;
	cout << sequences.size() << " compared\t" << NR2_red_no << " clustered" << endl;
	temp_files.Clear();
}

int calc_ann_list(int len, char *seqi, int NAA, int& aan_no, Vector<int> & aan_list, Vector<INTs> & aan_list_no, bool est) 
{
	int i, j, k, i0, i1, k1;

	// check_aan_list 
	aan_no = len - NAA + 1;
	for (j=0; j<aan_no; j++) {
		aan_list[j] = 0;
		for (k=0, k1=NAA-1; k<NAA; k++, k1--) aan_list[j] += seqi[j+k] * NAAN_array[k1];
	}
	if( est ){
		// for the short word containing 'N', mask it to '-1'
		for (j=0; j<len; j++){
			if ( seqi[j] >= 4 ) {                      // here N is 4
				i0 = (j-NAA+1 > 0) ? j-NAA+1 : 0;
				i1 = j < aan_no ? j : aan_no - 1;
				for (i=i0; i<=i1; i++) aan_list[i]=-1;
			}
		}
	}

	std::sort(aan_list.begin(), aan_list.begin() + aan_no);
	for(j=0; j<aan_no; j++) aan_list_no[j]=1;
	for(j=aan_no-1; j; j--) {
		if (aan_list[j] == aan_list[j-1]) {
			aan_list_no[j-1] += aan_list_no[j];
			aan_list_no[j]=0;
		}
	}
	return OK_FUNC;
} // END calc_ann_list

void make_comp_short_word_index(int NAA, int *NAAN_array, Vector<int> &Comp_AAN_idx) {
	int i, j, k, icomp, k1;
	int c[4] = {3,2,1,0};
	unsigned char short_word[32]; //short_word[12] is enough

	int NAA1 = NAAN_array[1];
	int NAAN = NAAN_array[NAA];

	for (i=0; i<NAAN; i++) {
		// decompose i back to short_word
		for (k=i, j=0; j<NAA; j++) {
			short_word[j] = (unsigned char) (k % NAA1);
			k  = k / NAA1;
		}

		// calc_comp_aan_list
		icomp=0;
		for (k=0, k1=NAA-1; k<NAA; k++, k1--) icomp += c[short_word[k1]] * NAAN_array[k];

		Comp_AAN_idx[i] = icomp;
	}
} // make_comp_short_word_index

//quick_sort_idx calling (a, idx, 0, no-1)
//sort a with another array idx
//so that idx rearranged
int quick_sort_idx (int *a, int *idx, int lo0, int hi0 ) {
  int lo = lo0;
  int hi = hi0;
  int mid;
  int tmp;

  if ( hi0 > lo0) {
    mid = a[ ( lo0 + hi0 ) / 2 ];

    while( lo <= hi ) {
      while( ( lo < hi0 ) && ( a[lo] < mid ) ) lo++;
      while( ( hi > lo0 ) && ( a[hi] > mid ) ) hi--;
      if( lo <= hi ) {
        tmp=a[lo];   a[lo]=a[hi];     a[hi]=tmp;
        tmp=idx[lo]; idx[lo]=idx[hi]; idx[hi]=tmp;
        lo++; hi--;
      }
    } // while

    if( lo0 < hi ) quick_sort_idx(a, idx, lo0, hi );
    if( lo < hi0 ) quick_sort_idx(a, idx, lo, hi0 );
  } // if ( hi0 > lo0)
  return 0;
} // quick_sort_idx


//decreasing can not use reverse of quick_sort_idx due to tie
//quick_sort_idxr calling (a, idx, 0, no-1)
//sort a with another array idx
//so that idx rearranged
int quick_sort_idxr (int *a, int *idx, int lo0, int hi0 ) {
  int lo = lo0;
  int hi = hi0;
  int mid;
  int tmp;

  if ( hi0 > lo0) {
    mid = a[ ( lo0 + hi0 ) / 2 ];

    while( lo <= hi ) {
      while( ( lo < hi0 ) && ( a[lo] > mid ) ) lo++;
      while( ( hi > lo0 ) && ( a[hi] < mid ) ) hi--;
      if( lo <= hi ) {
        tmp=a[lo];   a[lo]=a[hi];     a[hi]=tmp;
        tmp=idx[lo]; idx[lo]=idx[hi]; idx[hi]=tmp;
        lo++; hi--;
      }
    } // while

    if( lo0 < hi ) quick_sort_idxr(a, idx, lo0, hi );
    if( lo < hi0 ) quick_sort_idxr(a, idx, lo, hi0 );
  } // if ( hi0 > lo0)
  return 0;
} // quick_sort_idxr

/////////////////////////// END ALL ////////////////////////

int naa_stat_start_percent = 40;
int naa_stat[5][61][4] = {

	// cover 0.99
	{
		// N=5   N=4   N=3   N=2
		{  0,    0,    0,    7,  },  // 40%
		{  0,    0,    0,    8,  },  // 41%
		{  0,    0,    0,    9,  },  // 42%
		{  0,    0,    0,    9,  },  // 43%
		{  0,    0,    1,   10,  },  // 44%
		{  0,    0,    1,   11,  },  // 45%
		{  0,    0,    1,   12,  },  // 46%
		{  0,    0,    2,   13,  },  // 47%
		{  0,    0,    2,   14,  },  // 48%
		{  0,    0,    4,   16,  },  // 49%
		{  0,    0,    4,   16,  },  // 50%
		{  0,    0,    5,   17,  },  // 51%
		{  0,    0,    5,   18,  },  // 52%
		{  0,    0,    7,   20,  },  // 53%
		{  0,    1,    7,   21,  },  // 54%
		{  0,    1,    7,   21,  },  // 55%
		{  0,    2,    8,   23,  },  // 56%
		{  0,    2,    8,   25,  },  // 57%
		{  0,    2,   10,   25,  },  // 58%
		{  0,    3,   10,   26,  },  // 59%
		{  0,    4,   13,   28,  },  // 60%
		{  0,    5,   13,   30,  },  // 61%
		{  0,    5,   14,   30,  },  // 62%
		{  1,    6,   15,   33,  },  // 63%
		{  2,    7,   17,   34,  },  // 64%
		{  2,    7,   17,   35,  },  // 65%
		{  2,    9,   20,   37,  },  // 66%
		{  4,   10,   20,   37,  },  // 67%
		{  4,   11,   22,   40,  },  // 68%
		{  5,   12,   24,   41,  },  // 69%
		{  5,   12,   25,   42,  },  // 70%
		{  6,   16,   27,   43,  },  // 71%
		{  8,   16,   27,   45,  },  // 72%
		{  9,   17,   29,   47,  },  // 73%
		{ 10,   18,   31,   47,  },  // 74%
		{ 10,   20,   32,   50,  },  // 75%
		{ 12,   20,   32,   51,  },  // 76%
		{ 14,   22,   36,   54,  },  // 77%
		{ 15,   24,   37,   55,  },  // 78%
		{ 17,   26,   41,   58,  },  // 79%
		{ 18,   29,   41,   59,  },  // 80%
		{ 20,   30,   45,   60,  },  // 81%
		{ 24,   35,   48,   62,  },  // 82%
		{ 26,   36,   48,   64,  },  // 83%
		{ 27,   38,   51,   65,  },  // 84%
		{ 31,   43,   54,   68,  },  // 85%
		{ 35,   43,   55,   70,  },  // 86%
		{ 36,   48,   60,   71,  },  // 87%
		{ 36,   50,   61,   73,  },  // 88%
		{ 40,   50,   61,   75,  },  // 89%
		{ 45,   54,   65,   75,  },  // 90%
		{ 52,   60,   70,   79,  },  // 91%
		{ 53,   62,   71,   81,  },  // 92%
		{ 57,   66,   75,   84,  },  // 93%
		{ 57,   66,   76,   85,  },  // 94%
		{ 64,   71,   78,   85,  },  // 95%
		{ 70,   75,   82,   89,  },  // 96%
		{ 77,   81,   86,   92,  },  // 97%
		{ 82,   86,   90,   94,  },  // 98%
		{ 83,   87,   91,   95,  },  // 99%
		{ 91,   93,   95,   97,  },  // 100%
	},
	// cover 0.95
	{
		// N=5   N=4   N=3   N=2
		{  0,    0,    1,    9,  },  // 40%
		{  0,    0,    2,   10,  },  // 41%
		{  0,    0,    2,   11,  },  // 42%
		{  0,    0,    3,   12,  },  // 43%
		{  0,    0,    3,   12,  },  // 44%
		{  0,    0,    4,   14,  },  // 45%
		{  0,    0,    4,   14,  },  // 46%
		{  0,    1,    5,   16,  },  // 47%
		{  0,    1,    6,   17,  },  // 48%
		{  0,    2,    7,   19,  },  // 49%
		{  0,    2,    8,   19,  },  // 50%
		{  0,    2,    8,   20,  },  // 51%
		{  0,    2,    9,   21,  },  // 52%
		{  0,    4,   10,   23,  },  // 53%
		{  1,    4,   11,   24,  },  // 54%
		{  1,    4,   11,   24,  },  // 55%
		{  1,    5,   13,   26,  },  // 56%
		{  2,    5,   13,   27,  },  // 57%
		{  2,    6,   15,   29,  },  // 58%
		{  2,    7,   15,   30,  },  // 59%
		{  3,    8,   16,   31,  },  // 60%
		{  4,    8,   18,   32,  },  // 61%
		{  4,    9,   18,   33,  },  // 62%
		{  5,   11,   20,   36,  },  // 63%
		{  6,   12,   22,   37,  },  // 64%
		{  6,   12,   22,   38,  },  // 65%
		{  8,   14,   24,   40,  },  // 66%
		{  8,   15,   25,   41,  },  // 67%
		{ 10,   16,   27,   42,  },  // 68%
		{ 10,   18,   28,   45,  },  // 69%
		{ 11,   18,   29,   45,  },  // 70%
		{ 14,   21,   31,   47,  },  // 71%
		{ 14,   22,   32,   48,  },  // 72%
		{ 14,   22,   33,   50,  },  // 73%
		{ 17,   24,   36,   52,  },  // 74%
		{ 17,   25,   36,   52,  },  // 75%
		{ 18,   27,   39,   54,  },  // 76%
		{ 20,   29,   41,   56,  },  // 77%
		{ 21,   31,   42,   58,  },  // 78%
		{ 21,   31,   46,   60,  },  // 79%
		{ 27,   35,   46,   60,  },  // 80%
		{ 28,   37,   50,   63,  },  // 81%
		{ 31,   38,   50,   64,  },  // 82%
		{ 34,   43,   53,   66,  },  // 83%
		{ 36,   45,   54,   67,  },  // 84%
		{ 41,   50,   60,   70,  },  // 85%
		{ 43,   51,   60,   71,  },  // 86%
		{ 45,   54,   63,   74,  },  // 87%
		{ 48,   55,   64,   75,  },  // 88%
		{ 54,   60,   68,   78,  },  // 89%
		{ 55,   62,   71,   80,  },  // 90%
		{ 56,   63,   71,   80,  },  // 91%
		{ 64,   70,   76,   84,  },  // 92%
		{ 69,   74,   80,   86,  },  // 93%
		{ 73,   78,   83,   88,  },  // 94%
		{ 74,   78,   84,   89,  },  // 95%
		{ 80,   84,   87,   91,  },  // 96%
		{ 83,   86,   90,   93,  },  // 97%
		{ 86,   89,   92,   95,  },  // 98%
		{ 91,   93,   95,   97,  },  // 99%
		{ 92,   93,   95,   97,  },  // 100%
	},
	// cover 0.9
	{
		// N=5   N=4   N=3   N=2
		{  0,    0,    2,   11,  },  // 40%
		{  0,    0,    3,   12,  },  // 41%
		{  0,    0,    3,   12,  },  // 42%
		{  0,    1,    4,   13,  },  // 43%
		{  0,    1,    5,   14,  },  // 44%
		{  0,    1,    5,   15,  },  // 45%
		{  0,    1,    6,   16,  },  // 46%
		{  0,    2,    7,   18,  },  // 47%
		{  0,    2,    7,   18,  },  // 48%
		{  0,    3,    9,   20,  },  // 49%
		{  1,    4,    9,   20,  },  // 50%
		{  1,    4,   10,   21,  },  // 51%
		{  1,    4,   11,   23,  },  // 52%
		{  2,    5,   12,   24,  },  // 53%
		{  2,    5,   12,   25,  },  // 54%
		{  2,    6,   13,   26,  },  // 55%
		{  3,    7,   14,   28,  },  // 56%
		{  3,    7,   15,   28,  },  // 57%
		{  4,    8,   16,   30,  },  // 58%
		{  5,    9,   17,   31,  },  // 59%
		{  5,   10,   18,   32,  },  // 60%
		{  6,   11,   20,   35,  },  // 61%
		{  6,   11,   20,   35,  },  // 62%
		{  7,   13,   22,   38,  },  // 63%
		{  8,   14,   23,   39,  },  // 64%
		{  8,   15,   24,   39,  },  // 65%
		{ 10,   16,   26,   42,  },  // 66%
		{ 10,   17,   27,   42,  },  // 67%
		{ 12,   19,   29,   44,  },  // 68%
		{ 13,   20,   30,   46,  },  // 69%
		{ 13,   21,   31,   47,  },  // 70%
		{ 16,   23,   33,   48,  },  // 71%
		{ 18,   25,   34,   50,  },  // 72%
		{ 18,   26,   36,   51,  },  // 73%
		{ 19,   28,   38,   53,  },  // 74%
		{ 20,   29,   38,   53,  },  // 75%
		{ 23,   30,   41,   56,  },  // 76%
		{ 24,   33,   43,   57,  },  // 77%
		{ 26,   34,   45,   59,  },  // 78%
		{ 28,   37,   48,   61,  },  // 79%
		{ 30,   37,   48,   62,  },  // 80%
		{ 33,   42,   52,   64,  },  // 81%
		{ 35,   43,   53,   65,  },  // 82%
		{ 38,   47,   56,   68,  },  // 83%
		{ 40,   47,   56,   68,  },  // 84%
		{ 44,   53,   61,   71,  },  // 85%
		{ 45,   53,   62,   73,  },  // 86%
		{ 50,   58,   66,   75,  },  // 87%
		{ 51,   58,   66,   76,  },  // 88%
		{ 57,   63,   71,   79,  },  // 89%
		{ 60,   66,   72,   81,  },  // 90%
		{ 62,   68,   75,   83,  },  // 91%
		{ 70,   74,   80,   85,  },  // 92%
		{ 74,   78,   82,   88,  },  // 93%
		{ 85,   87,   90,   92,  },  // 94%
		{ 86,   88,   90,   92,  },  // 95%
		{ 87,   89,   91,   93,  },  // 96%
		{ 87,   89,   92,   94,  },  // 97%
		{ 89,   91,   93,   96,  },  // 98%
		{ 93,   94,   96,   97,  },  // 99%
		{ 94,   95,   97,   98,  },  // 100%
	},
	// cover 0.8
	{
		// N=5   N=4   N=3   N=2
		{  0,    1,    4,   13,  },  // 40%
		{  0,    1,    5,   13,  },  // 41%
		{  0,    1,    5,   14,  },  // 42%
		{  0,    2,    6,   15,  },  // 43%
		{  0,    2,    6,   16,  },  // 44%
		{  0,    2,    7,   17,  },  // 45%
		{  1,    3,    8,   18,  },  // 46%
		{  1,    4,    9,   20,  },  // 47%
		{  1,    4,    9,   20,  },  // 48%
		{  2,    5,   11,   22,  },  // 49%
		{  2,    5,   11,   22,  },  // 50%
		{  2,    6,   12,   24,  },  // 51%
		{  3,    6,   13,   25,  },  // 52%
		{  3,    7,   14,   26,  },  // 53%
		{  4,    8,   14,   27,  },  // 54%
		{  4,    8,   15,   28,  },  // 55%
		{  5,    9,   17,   30,  },  // 56%
		{  5,    9,   17,   30,  },  // 57%
		{  6,   11,   19,   32,  },  // 58%
		{  7,   12,   20,   34,  },  // 59%
		{  8,   12,   20,   34,  },  // 60%
		{  9,   14,   22,   37,  },  // 61%
		{  9,   14,   23,   37,  },  // 62%
		{ 10,   16,   25,   39,  },  // 63%
		{ 11,   17,   26,   41,  },  // 64%
		{ 12,   18,   27,   41,  },  // 65%
		{ 13,   20,   28,   43,  },  // 66%
		{ 14,   21,   30,   45,  },  // 67%
		{ 15,   22,   31,   46,  },  // 68%
		{ 17,   24,   33,   48,  },  // 69%
		{ 17,   24,   34,   48,  },  // 70%
		{ 19,   26,   36,   50,  },  // 71%
		{ 20,   27,   37,   51,  },  // 72%
		{ 21,   29,   39,   53,  },  // 73%
		{ 23,   31,   41,   55,  },  // 74%
		{ 23,   31,   41,   55,  },  // 75%
		{ 26,   34,   44,   58,  },  // 76%
		{ 28,   36,   46,   59,  },  // 77%
		{ 29,   37,   47,   60,  },  // 78%
		{ 34,   41,   50,   62,  },  // 79%
		{ 34,   42,   51,   63,  },  // 80%
		{ 38,   45,   55,   66,  },  // 81%
		{ 39,   46,   55,   67,  },  // 82%
		{ 44,   51,   60,   70,  },  // 83%
		{ 44,   51,   60,   70,  },  // 84%
		{ 49,   56,   64,   73,  },  // 85%
		{ 50,   57,   64,   74,  },  // 86%
		{ 57,   63,   69,   77,  },  // 87%
		{ 58,   64,   70,   78,  },  // 88%
		{ 68,   71,   76,   82,  },  // 89%
		{ 68,   72,   77,   83,  },  // 90%
		{ 75,   79,   81,   85,  },  // 91%
		{ 86,   87,   89,   90,  },  // 92%
		{ 88,   89,   90,   92,  },  // 93%
		{ 90,   91,   92,   93,  },  // 94%
		{ 91,   92,   93,   94,  },  // 95%
		{ 92,   94,   94,   95,  },  // 96%
		{ 93,   94,   95,   96,  },  // 97%
		{ 94,   95,   95,   96,  },  // 98%
		{ 94,   95,   96,   98,  },  // 99%
		{ 95,   96,   97,   98,  },  // 100%
	},
	// cover 0.6
	{
		// N=5   N=4   N=3   N=2
		{  1,    2,    6,   15,  },  // 40%
		{  1,    3,    7,   16,  },  // 41%
		{  1,    3,    8,   17,  },  // 42%
		{  2,    4,    9,   18,  },  // 43%
		{  2,    4,    9,   19,  },  // 44%
		{  2,    5,   10,   20,  },  // 45%
		{  3,    5,   10,   21,  },  // 46%
		{  3,    6,   12,   22,  },  // 47%
		{  3,    6,   12,   23,  },  // 48%
		{  4,    8,   14,   25,  },  // 49%
		{  4,    8,   14,   25,  },  // 50%
		{  5,    8,   15,   26,  },  // 51%
		{  5,    9,   16,   27,  },  // 52%
		{  6,   10,   17,   29,  },  // 53%
		{  6,   11,   18,   30,  },  // 54%
		{  7,   11,   18,   31,  },  // 55%
		{  8,   12,   20,   32,  },  // 56%
		{  8,   13,   20,   33,  },  // 57%
		{ 10,   14,   22,   35,  },  // 58%
		{ 10,   15,   23,   37,  },  // 59%
		{ 11,   16,   24,   37,  },  // 60%
		{ 12,   18,   26,   39,  },  // 61%
		{ 13,   18,   26,   40,  },  // 62%
		{ 14,   20,   28,   42,  },  // 63%
		{ 16,   22,   30,   43,  },  // 64%
		{ 16,   22,   31,   44,  },  // 65%
		{ 17,   23,   32,   45,  },  // 66%
		{ 18,   25,   33,   47,  },  // 67%
		{ 19,   26,   35,   48,  },  // 68%
		{ 21,   27,   36,   50,  },  // 69%
		{ 22,   29,   37,   51,  },  // 70%
		{ 24,   30,   39,   52,  },  // 71%
		{ 25,   32,   41,   53,  },  // 72%
		{ 26,   33,   42,   55,  },  // 73%
		{ 29,   35,   44,   57,  },  // 74%
		{ 29,   36,   45,   57,  },  // 75%
		{ 32,   39,   48,   60,  },  // 76%
		{ 34,   41,   50,   61,  },  // 77%
		{ 36,   43,   51,   62,  },  // 78%
		{ 40,   46,   54,   65,  },  // 79%
		{ 40,   46,   54,   65,  },  // 80%
		{ 46,   52,   59,   68,  },  // 81%
		{ 46,   52,   60,   69,  },  // 82%
		{ 53,   59,   65,   73,  },  // 83%
		{ 54,   60,   66,   73,  },  // 84%
		{ 63,   67,   73,   78,  },  // 85%
		{ 68,   71,   75,   79,  },  // 86%
		{ 78,   80,   82,   85,  },  // 87%
		{ 79,   81,   83,   85,  },  // 88%
		{ 83,   85,   86,   87,  },  // 89%
		{ 85,   86,   87,   89,  },  // 90%
		{ 86,   88,   89,   90,  },  // 91%
		{ 88,   89,   90,   91,  },  // 92%
		{ 90,   90,   91,   92,  },  // 93%
		{ 91,   92,   92,   93,  },  // 94%
		{ 92,   93,   94,   94,  },  // 95%
		{ 94,   94,   95,   95,  },  // 96%
		{ 95,   95,   96,   96,  },  // 97%
		{ 95,   96,   97,   97,  },  // 98%
		{ 96,   96,   97,   98,  },  // 99%
		{ 97,   98,   98,   99,  },  // 100%
	},
};

