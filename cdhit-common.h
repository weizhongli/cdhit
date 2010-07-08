// =============================================================================
// CD-HI/CD-HIT
//
// Cluster Database at High Identity Threshold
//
// CD-HIT clusters protein sequence database at high sequence identity threshold.
// This program can remove the high sequence redundance efficiently.
//
// program written by
//                                      Weizhong Li
//                                      UCSD, San Diego Supercomputer Center
//                                      La Jolla, CA, 92093
//                                      Email liwz@sdsc.edu
//
//                 at
//                                      Adam Godzik's lab
//                                      The Burnham Institute
//                                      La Jolla, CA, 92037
//                                      Email adam@burnham-inst.org
//
// modified by:
//                                      Limin Fu
//                                      Calit2, UCSD
//                                      La Jolla, CA, 92093
//                                      Email: l2fu@ucsd.edu, fu@daovm.net
// =============================================================================

#include<iostream>
#include<fstream>
#include<iomanip>
#include<cstdlib>
#include<stdio.h>
#include<string.h>
#include<ctype.h>
#include<sys/times.h>

#include<valarray>
#include<vector>
#include<map>

#define MAX_AA 23
#define MAX_NA 6
#define MAX_UAA 21
#define MAX_SEQ 655360
#define MAX_DIAG (MAX_SEQ<<1)              // MAX_DIAG be twice of MAX_SEQ
#define MAX_GAP MAX_SEQ                    // MAX_GAP <= MAX_SEQ
#define MAX_DES 300000
#define MAX_LINE_SIZE 300000
#define MAX_FILE_NAME 1280
#define MAX_SEG 50
#define MAX_BIN_SWAP 2E9
#define CLOCK_TICKS 100
#define FAILED_FUNC 1
#define OK_FUNC 0

#define IS_REP 1
#define IS_REDUNDANT 2
#define IS_PROCESSED 16
#define IS_MINUS_STRAND 32

#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))

typedef unsigned int UINT4;
typedef unsigned short UINT2;

#define LONG_SEQ

//if the longset sequence is longer than 65535, I use INT4
#ifdef LONG_SEQ
#define INTs UINT4
#else
#define INTs UINT2
#endif

using namespace std;

// the parent containter must guarantee continuous memory allocation.
// std::valarray could be used instead of std::vector.
template<class TYPE>
class Vector : public vector<TYPE>
{
  public:
  Vector() : vector<TYPE>(){}
  Vector( size_t size ) : vector<TYPE>( size ){}
  Vector( size_t size, const TYPE & deft ) : vector<TYPE>( size, deft ){}

  void Append( const TYPE & item ){
    size_t n = this->size();
    if( n + 1 >= this->capacity() ) this->reserve( n + n/5 + 1 );
    push_back( item );
  }
};

// for primitive types only
template<class TYPE>
class NVector
{
  public:
	TYPE   *items;
	int     size;
	int     capacity;

  NVector(){ size = capacity = 0; items = NULL; }
  NVector( int n, const TYPE & v=TYPE() ){ 
    size = capacity = 0; items = NULL; 
    Resize( n, v );
  }
  NVector( const NVector & other ){
    size = capacity = 0; items = NULL; 
    if( other.items ){
      Resize( other.size );
      memcpy( items, other.items, other.size * sizeof(TYPE) );
    }
  }

  ~NVector(){ if( items ) free( items ); }

  int  Size()const{ return size; }
  void Clear(){
    if( items ) free( items );
    size = capacity = 0; items = NULL; 
  }

  void Resize( int n, const TYPE & value=TYPE() ){
    if( n == size && capacity > 0 ) return;
    int i;
    // When resize() is called, probably this is the intended size,
    // and will not be changed frequently.
    if( n != capacity ){
      capacity = n;
      items = (TYPE*)realloc( items, (capacity+1)*sizeof(TYPE) );
    }
    for(i=size; i<n; i++ ) items[i] = value;
    size = n;
  }
  void Append( const TYPE & item ){
    if( size + 1 >= capacity ){
      capacity = size + size/5 + 1;
      items = (TYPE*)realloc( items, (capacity+1)*sizeof(TYPE) );
    }
    items[size] = item;
    size ++;
  }

	TYPE& operator[]( const int i ){
    //if( i <0 or i >= size ) printf( "out of range\n" );
    return items[i];
  }
	TYPE& operator[]( const int i )const{
    //if( i <0 or i >= size ) printf( "out of range\n" );
    return items[i];
  }
};
typedef NVector<int>      VectorInt;
typedef Vector<VectorInt> MatrixInt;

////////// Class definition //////////
class ScoreMatrix { //Matrix
  private:

  public:
    int matrix[MAX_AA][MAX_AA];
    int gap, ext_gap;
    VectorInt gap_array;

    ScoreMatrix();
    void init();
    void set_gap(int gap1, int ext_gap1);
    void set_matrix(int *mat1);
    void set_to_na();
}; // END class ScoreMatrix


typedef NVector<INTs>      VectorIntX;
typedef Vector<VectorIntX> MatrixIntX;


extern int NAA1 ;
extern int NAA2 ;
extern int NAA3 ;
extern int NAA4 ;
extern int NAA5 ;
extern int NAA6 ;
extern int NAA7 ;
extern int NAA8 ;
extern int NAA9 ;
extern int NAA10;
extern int NAA11;
extern int NAA12;
extern int NAAN_array[13];

void InitNAA( int max );


extern int naa_stat_start_percent;
extern int naa_stat[5][61][4];

struct IndexCount
{
  int index;
  int count;

  IndexCount( int i=0, int c=0 ){ index = i, count = c; }
};

struct Sequence;

class WordTable
{
  private:
  public:
    Vector<NVector<IndexCount> > indexCounts; // hold index and word counts of seqs
    Vector<Sequence*>            sequences;
    int     NAA;                // length of word
    int     NAAN;               // rows of table
    char    is_aa;              // aa is for prot
    size_t  size;
    int     frag_count;

  public:
    WordTable( int naa=0, int naan=0 );
    void Init(int, int);
    void Clear();
    void SetDNA();
    int  AddWordCounts( NVector<IndexCount> & counts, Sequence *seq, bool skipN=false);
    int  AddWordCountsFrag( NVector<IndexCount> & counts, int frag, int frag_size, int repfrag );
    int  CountWords( NVector<IndexCount> & counts, Vector<INTs> & look_and_count, bool est=false);

    int  AddWordCounts(int aan_no, Vector<int> & word_encodes, 
        Vector<INTs> & word_encodes_no, int idx, bool skipN=false);
    int AddWordCountsFrag( int aan_no, Vector<int> & word_encodes, 
        Vector<INTs> & word_encodes_no, int frag, int frag_size );
    int CountWords(int aan_no, Vector<int> & aan_list,
        Vector<INTs> & aan_list_no, Vector<INTs> & look_and_count, bool est=false);
    void PrintAll();
}; // END class INDEX_TBL
struct Options
{
  int     NAA;
  int     NAAN;
  int     NAA_top_limit;

	size_t  max_memory; // -M: 400,000,000 in bytes
	int     min_length; // -l: 10 bases
	bool    cluster_best;  // -g: 0, the first; 1, the best
	bool    global_identity; // -G:
	bool    store_disk; // -B:
	int     band_width; // -b: 20
  int     cluster_thd100;
  double  cluster_thd; // -c
	double  diff_cutoff; // -s: 0.0
	double  diff_cutoff2; // -s2: 1.0
	int     diff_cutoff_aa; // -S: 999999
	int     diff_cutoff_aa2; // -S2: 0
	int     tolerance; // -t: 2
	double  long_coverage; // -aL:
	int     long_control; // -AL:
	double  short_coverage; // -aS:
	int     short_control; // -AS:
	int     min_control; // -A:
  int     print;
  int     des_len;
  int     frag_size;
  int     option_r;
  int     threads;

  size_t  mem_limit;
  size_t  total_letters;

  bool    has2D;
  bool    isEST;

  string  input;
  string  input2;
  string  output;

	Options(){
    has2D = false;
    isEST = false;
    NAA = 5;
    NAA_top_limit = 5;
    cluster_thd = 0.9;
		max_memory = 800000000;
		min_length = 10;
		cluster_best = false;
		global_identity = true;
		store_disk = false;
		band_width = 20;
		diff_cutoff = 0.0;
		diff_cutoff2 = 1.0;
		diff_cutoff_aa = 999999;
		diff_cutoff_aa2 = 0;
		tolerance = 2;
		long_coverage = 0.0;
		long_control = 99999999;
		short_coverage = 0.0;
		short_control = 99999999;
		min_control = 0;
    print = 0;
    option_r  = 0;
    frag_size = 0;
    des_len = 20;
    threads = 1;
    mem_limit = 100000000;
    total_letters = 0;
	};

  bool SetOptionCommon( const char *flag, const char *value );
  bool SetOption( const char *flag, const char *value );
  bool SetOption2D( const char *flag, const char *value );
  bool SetOptionEST( const char *flag, const char *value );
  bool SetOptions( int argc, char *argv[], bool twodata=false, bool est=false );

  void Validate();

  void Print();
};

void bomb_error(const char *message);


struct Sequence
{
  // real sequence, if it is not stored swap file:
  char *data;
  // length of the sequence:
  int   size;
  int   bufsize;

  // if swap != NULL, the sequence is stored in file.
  // swap is opened as temporary file, which will be deleted automatically
  // after the program is finished:
  FILE *swap;
  // stream offset of the sequence:
  int   offset;

  // stream offset of the description string in the database:
  size_t   des_begin;
  // length of the description:
  int   des_length;
  // length of data in fasta file, including line wrapping:
  int   dat_length;

  // index of the sequence in the original database:
  int   index;
  short state;
  short xletter; // number of un-encoded letters
  int   cluster_id;
  int   fragment; // the first fragment, for representatives only
  int   identity;
  int   coverage[4];

  Sequence();
  Sequence( const Sequence & other );
  ~Sequence();

  void Clear();

  void operator=( const char *s );
  void operator+=( const char *s );

  void Resize( int n );
  void Reserve( int n );

  void Swap( Sequence & other );
  void Format();

  void SwapIn();
  void SwapOut();
  void PrintInfo( int id, FILE *fin, FILE *fout, const Options & options );
};

struct WorkingParam
{
  double aa1_cutoff;
  double aas_cutoff; /* or aa2 */
  double aan_cutoff;
  int    len_upper_bound;
  int    len_lower_bound;

  WorkingParam( double a1=0, double a2=0, double an=0 ){
    Set( a1, a2, an );
  }
  void Set( double a1=0, double a2=0, double an=0 ){
    aa1_cutoff = a1;
    aas_cutoff = a2;
    aan_cutoff = an;
    len_upper_bound = 0;
    len_lower_bound = 0;
  }

  int len_eff;
  int aln_cover_flag;
  int min_aln_lenS;
  int min_aln_lenL;
  int required_aa1;
  int required_aas; /* or aa2 */
  int required_aan;

  void ControlShortCoverage( int len, const Options & option );
  void ControlLongCoverage( int len, const Options & option );
  void ComputeRequiredBases( int NAA, int ss=2 );
};
struct WorkingBuffer
{
  Vector<int>  taap;
  Vector<int>  word_encodes;
  Vector<int>  word_encodes_backup;
  Vector<INTs> word_encodes_no;
  Vector<INTs> aap_list;
  Vector<INTs> aap_begin;
  Vector<INTs> look_and_count;
  Vector<IndexCount>  indexCounts;
  Vector<int>  diag_score;
  Vector<int> aan_list_comp;
  char seqi_comp[MAX_SEQ];

  WorkingBuffer( int frag=0, bool est=false ){
    Set( frag, est );
  }
  void Set( int frag, bool est=false ){
    int m = MAX_UAA*MAX_UAA;
    if( est ) m = m * m;
    taap.resize( m );
    aap_list.resize( MAX_SEQ );
    aap_begin.resize( m );
    indexCounts.resize( MAX_SEQ );
    word_encodes.resize( MAX_SEQ );
    word_encodes_no.resize( MAX_SEQ );
    word_encodes_backup.resize( MAX_SEQ );
    look_and_count.resize( frag );
    diag_score.resize( MAX_DIAG );
    aan_list_comp.resize( MAX_SEQ );
  }

  int EncodeWords( Sequence *seq, int NA, bool est = false );
  void ComputeAAP( const char *seqi, int size );
  void ComputeAAP2( const char *seqi, int size );
};
extern Vector<int>  Comp_AAN_idx;
extern ScoreMatrix  mat;


class SequenceDB
{
public:

  int NAAN;
  Vector<Sequence*>  sequences;
  Vector<int>        rep_seqs;

  void Clear(){
    for(int i=0; i<sequences.size(); i++) delete sequences[i];
    sequences.clear(); rep_seqs.clear();
  }

  SequenceDB(){ }
  ~SequenceDB(){ Clear(); }

  void Read( const char *file, const Options & options );
  void WriteClusters( const char *db, const char *newdb, const Options & options );
  void WriteExtra1D( const Options & options );
  void WriteExtra2D( SequenceDB & other, const Options & options );
  void DivideSave( const char *db, const char *newdb, int n, const Options & options );

  void SwapIn( int seg, bool reponly=false );
  void SwapOut( int seg );

  void SortDivide( Options & options, bool sort=true );
  void MakeWordTable( const Options & optioins );

  void ClusterOne( Sequence *seq, int id, WordTable & table,
      WorkingParam & param, WorkingBuffer & buf, const Options & options );

  //void SelfComparing( int start, int end, WordTable & table, 
  //    WorkingParam & param, WorkingBuffer & buf, const Options & options );

  void DoClustering( const Options & options );
  void DoClustering( int T, const Options & options );
  void ClusterTo( SequenceDB & other, const Options & optioins );
  int  CheckOne( Sequence *seq, WordTable & tab, WorkingParam & par, WorkingBuffer & buf, const Options & opt );
  int  CheckOneEST( Sequence *seq, WordTable & tab, WorkingParam & par, WorkingBuffer & buf, const Options & opt );
  int  CheckOneAA( Sequence *seq, WordTable & tab, WorkingParam & par, WorkingBuffer & buf, const Options & opt );
};


int print_usage (char *arg);
void bomb_error(const char *message);
void bomb_error(const char *message, const char *message2);
void bomb_warning(const char *message);
void bomb_warning(const char *message, const char *message2);
void format_seq(char *seq);
int diag_test_aapn(int NAA1, char iseq2[], int len1, int len2, 
        WorkingBuffer & buffer, int &best_sum,
        int band_width, int &band_left, int &band_right, int required_aa1);
int diag_test_aapn_est(int NAA1, char iseq2[], int len1, int len2, 
        WorkingBuffer & buffer, int &best_sum,
        int band_width, int &band_left, int &band_right, int required_aa1);
int local_band_align(char iseq1[], char iseq2[], int len1, int len2,
                     ScoreMatrix &mat, int &best_score, int &iden_no,
                     int band_left, int band_right);
int local_band_align2(char iseq1[], char iseq2[], int len1, int len2,
                     ScoreMatrix &mat, int &best_score, int &iden_no,
                     int band_left, int band_right,
                     int &from1, int &end1, int &from2, int &end2, int &alnln);
int setiseq(char *seq, int len);

int print_usage_2d (char *arg);
int print_usage_est (char *arg);
int print_usage_div (char *arg);
int print_usage_est_2d (char *arg);

int upper_bound_length_rep(int len, const Options & options );
void cal_aax_cutoff (double &aa1_cutoff, double &aa2_cutoff, double &aan_cutoff,
                     double NR_clstr, int tolerance, int naa_stat_start_percent,
                     int naa_stat[5][61][4], int NAA);
void update_aax_cutoff(double &aa1_cutoff, double &aa2_cutoff, double &aan_cutoff,
                     int tolerance, int naa_stat_start_percent,
                     int naa_stat[5][61][4], int NAA, int iden);
void show_cpu_time(tms &CPU_begin, tms &CPU_end);

int calc_ann_list(int len, char *seqi, int NAA, int& aan_no, Vector<int> & aan_list, Vector<INTs> & aan_list_no, bool est=false);

