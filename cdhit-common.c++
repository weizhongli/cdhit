// =============================================================================
// CD-HI/CD-HIT
//
// Cluster Database at High Identity
//
// CD-HIT clusters protein sequences at high identity threshold.
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

#include "cdhit-common.h"
#include<valarray>
#include<stdint.h>
#include<assert.h>

#ifndef NO_OPENMP

#include<omp.h>

#else

#define omp_set_num_threads(T) (T = T)
#define omp_get_thread_num() 0

#endif

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
    sprintf( buf + len, "%x", this );
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

ScoreMatrix mat;
Vector<int> Comp_AAN_idx;

void make_comp_iseq(int len, char *iseq_comp, char *iseq)
{
  int i, c[5] = {3,2,1,0,4};
  for (i=0; i<len; i++) iseq_comp[i] = c[ (int)iseq[len-i-1] ];
} // make_comp_iseq

bool Options::SetOptionCommon( const char *flag, const char *value )
{
  int intval = atoi( value );
  if      (strcmp(flag, "-i" ) == 0) input = value;
  else if (strcmp(flag, "-o" ) == 0) output = value;
  else if (strcmp(flag, "-M" ) == 0) max_memory  = atoll(value) * 1000000;
  else if (strcmp(flag, "-l" ) == 0) min_length  = intval;
  else if (strcmp(flag, "-c" ) == 0) cluster_thd  = atof(value);
  else if (strcmp(flag, "-b" ) == 0) band_width  = intval;
  else if (strcmp(flag, "-n" ) == 0) NAA       = intval;
  else if (strcmp(flag, "-d" ) == 0) des_len   = intval;
  else if (strcmp(flag, "-s" ) == 0) diff_cutoff  = atof(value);
  else if (strcmp(flag, "-S" ) == 0) diff_cutoff_aa  = intval;
  else if (strcmp(flag, "-B" ) == 0) store_disk  = intval;
  else if (strcmp(flag, "-p" ) == 0) print  = intval;
  else if (strcmp(flag, "-g" ) == 0) cluster_best  = intval;
  else if (strcmp(flag, "-G" ) == 0) global_identity  = intval;
  else if (strcmp(flag, "-aL") == 0) long_coverage = atof(value);
  else if (strcmp(flag, "-AL") == 0) long_control = intval;
  else if (strcmp(flag, "-aS") == 0) short_coverage = atof(value);
  else if (strcmp(flag, "-AS") == 0) short_control = intval;
  else if (strcmp(flag, "-A" ) == 0) min_control  = intval;
  else if (strcmp(flag, "-tmp" ) == 0) temp_dir  = value;
  else if (strcmp(flag, "-T" ) == 0){
#ifndef NO_OPENMP
    int cpu = omp_get_num_procs();
    threads  = intval;
    if( threads > cpu ){
      threads = cpu;
      printf( "Warning: total number of CPUs in the system is %i\n", cpu );
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
  else return false;
  return true;
}
bool Options::SetOptions( int argc, char *argv[], bool twod, bool est )
{
  int i;
  has2D = twod;
  isEST = est;
  for (i=1; i+1<argc; i+=2) if ( SetOption( argv[i], argv[i+1] ) == 0) return false;
  if( i < argc ) return false;
  cluster_thd100 = (int)(cluster_thd * 100);

  atexit( CleanUpTempFiles );
  return true;
}
void Options::Validate()
{
  if( isEST ){
    if ((cluster_thd > 1.0) || (cluster_thd < 0.8)) bomb_error("invalid clstr threshold, should >=0.8");
  }else{
    if ((cluster_thd > 1.0) || (cluster_thd < 0.4)) bomb_error("invalid clstr");
  }
  if (band_width < 1 ) bomb_error("invalid band width");
  if (NAA < 2 || NAA > NAA_top_limit) bomb_error("invalid word length");
  if (des_len < 0 ) bomb_error("too short description, not enough to identify sequences");
  if (not isEST && (tolerance < 0 || tolerance > 5) ) bomb_error("invalid tolerance");
  if ((diff_cutoff<0) || (diff_cutoff>1)) bomb_error("invalid value for -s");
  if (diff_cutoff_aa<0) bomb_error("invalid value for -S");
  if( has2D ){
    if ((diff_cutoff2<0) || (diff_cutoff2>1)) bomb_error("invalid value for -s2");
    if (diff_cutoff_aa2<0) bomb_error("invalid value for -S2");
  }
  if (global_identity == 0) print = 1;
  if (short_coverage < long_coverage) short_coverage = long_coverage;
  if (short_control > long_control) short_control = long_control;
  if ((global_identity == 0) && (short_coverage == 0.0) && (min_control == 0))
    bomb_error("You are using local identity, but no -aS -aL -A option");
  if (frag_size < 0) bomb_error("invalid fragment size");

  const char *message = "Your word length is %i, using %i may be faster!\n";
  if ( not isEST && tolerance ) {
    int i, clstr_idx = (int) (cluster_thd * 100) - naa_stat_start_percent;
    int tcutoff = naa_stat[tolerance-1][clstr_idx][5-NAA];

    if (tcutoff < 5 )
      bomb_error("Too short word length, increase it or the tolerance");
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

  if ( min_length <= NAA ) bomb_error("Too short -l, redefine it");
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
  printf( "cluster_thd100 = %i\n", cluster_thd100 );
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

void format_seq(char *seq) {
  int i, j;
  char c1;
  int len = strlen(seq);

  for (i=0,j=0; i<len; i++) {
    c1 = toupper(seq[i]);
    if ( isalpha(c1) ) seq[j++] = c1;
  }
  seq[j] = 0;
} // END void format_seq


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
        int &best_sum, int band_width, int &band_left, int &band_right, int required_aa1)
{
  int i, i1, j, k;
  int *pp;
  int nall = len1+len2-1;
  Vector<int> & taap = buffer.taap;
  Vector<INTs> & aap_begin = buffer.aap_begin;
  Vector<INTs> & aap_list = buffer.aap_list;
  Vector<int> & diag_score = buffer.diag_score;

  if (nall > MAX_DIAG) bomb_error("in diag_test_aapn, MAX_DIAG reached");
  for (pp=&diag_score[0], i=nall; i; i--, pp++) *pp=0;

  int c22;
  INTs *bip;
  int len11 = len1-1;
  int len22 = len2-1;
  i1 = len11;
  for (i=0; i<len22; i++,i1++) {
//    c22 = iseq2[i]*NAA1 + iseq2[i+1];
    c22 = iseq2[i]*NAA1+ iseq2[i+1];
    if ( (j=taap[c22]) == 0) continue;
#if 0
    int bi, bj;
    bip = aap_list+ aap_begin[c22];     //    bi = aap_begin[c22];
    for (; j; j--, bip++) {  //  for (j=0; j<taap[c22]; j++,bi++) {
      diag_score[i1 - *bip]++;
    }
#endif
    int m = aap_begin[c22];
    for(int k=0; k<j; k++) diag_score[ i1 - aap_list[m+k] ] ++;
  }

  //find the best band range
//  int band_b = required_aa1;
  int band_b = required_aa1-1 >= 0 ? required_aa1-1:0;  // on dec 21 2001
  int band_e = nall - required_aa1;
  int band_m = ( band_b+band_width-1 < band_e ) ? band_b+band_width-1 : band_e;
  int best_score=0;
  for (i=band_b; i<=band_m; i++) best_score += diag_score[i];
  int from=band_b;
  int end =band_m;
  int score = best_score;  
  for (k=from, j=band_m+1; j<band_e; j++) {
    score -= diag_score[k++]; 
    score += diag_score[j]; 
    if ( score > best_score ) {
      from = k;
      end  = j;
      best_score = score;
    }
  }
  for (j=from; j<=end; j++) { // if aap pairs fail to open gap
    if ( diag_score[j] < 5 ) { best_score -= diag_score[j]; from++;}
    else break;
  }
  for (j=end; j>=from; j--) { // if aap pairs fail to open gap
    if ( diag_score[j] < 5 ) { best_score -= diag_score[j]; end--;}
    else break;
  }

//  delete [] diag_score;
  band_left = from-len1+1; 
  band_right= end-len1+1;
  best_sum = best_score;
  return OK_FUNC;
}
// END diag_test_aapn
 

int diag_test_aapn_est(int NAA1, char iseq2[], int len1, int len2, WorkingBuffer & buffer, 
        int &best_sum, int band_width, int &band_left, int &band_right, int required_aa1)
{
  int i, i1, j, k;
  int *pp;
  int nall = len1+len2-1;
  int NAA2 = NAA1 * NAA1;
  int NAA3 = NAA2 * NAA1;
  Vector<int> & taap = buffer.taap;
  Vector<INTs> & aap_begin = buffer.aap_begin;
  Vector<INTs> & aap_list = buffer.aap_list;
  Vector<int> & diag_score = buffer.diag_score;

  if (nall > MAX_DIAG) bomb_error("in diag_test_aapn_est, MAX_DIAG reached");
  for (pp=&diag_score[0], i=nall; i; i--, pp++) *pp=0;

  int c22;
  INTs *bip;
  int len22 = len2-3;
  i1 = len1-1;
  for (i=0; i<len22; i++,i1++) {
    if ((iseq2[i]==4) || (iseq2[i+1]==4) || (iseq2[i+2]==4) || (iseq2[i+3]==4)) continue; //skip N

    c22 = iseq2[i]*NAA3+ iseq2[i+1]*NAA2 + iseq2[i+2]*NAA1 + iseq2[i+3];
    if ( (j=taap[c22]) == 0) continue;
    bip = & aap_list[ aap_begin[c22] ];     //    bi = aap_begin[c22];
    for (; j; j--, bip++) {  //  for (j=0; j<taap[c22]; j++,bi++) {
      diag_score[i1 - *bip]++;
    }
  }
#if 0
  int mmax = 0;
  int immax = 0;
  for(i=0; i<=nall; i++){
    printf( "%3i\t", diag_score[i] );
    if( i%10 ==0 or i == nall ) printf( "\n" );
    if( diag_score[i] > mmax ){
      mmax = diag_score[i];
      immax = i;
    }
  }
#endif

  //find the best band range
//  int band_b = required_aa1;
  int band_b = required_aa1-1 >= 0 ? required_aa1-1:0;  // on dec 21 2001
  int band_e = nall - required_aa1;
  int band_m = ( band_b+band_width-1 < band_e ) ? band_b+band_width-1 : band_e;
  int best_score=0;
  for (i=band_b; i<=band_m; i++) best_score += diag_score[i];
  int from=band_b;
  int end =band_m;
  int score = best_score;  
  
#if 0
  printf( "%i\n", required_aa1 );
  printf( "max=%3i  imax=%3i; band:  %3i  %3i  %i\n", mmax, immax, band_b, band_e, band_m );
  printf( "best: %i\n", best_score );
#endif
  for (k=from, j=band_m+1; j<band_e; j++) {
    score -= diag_score[k++]; 
    score += diag_score[j]; 
    if ( score > best_score ) {
      from = k;
      end  = j;
      best_score = score;
    }
  }
  //printf( "best: %i\n", best_score );
  for (j=from; j<=end; j++) { // if aap pairs fail to open gap
    if ( diag_score[j] < 5 ) { best_score -= diag_score[j]; from++;}
    else break;
  }
  for (j=end; j>=from; j--) { // if aap pairs fail to open gap
    if ( diag_score[j] < 5 ) { best_score -= diag_score[j]; end--;}
    else break;
  }
  //printf( "best: %i\n", best_score );

//  delete [] diag_score;
  band_left = from-len1+1; 
  band_right= end-len1+1;
  best_sum = best_score;
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

int local_band_align(char iseq1[], char iseq2[], int len1, int len2,
                     ScoreMatrix &mat, int &best_score, int &iden_no,
                     int band_left, int band_right) {
  int i, j, k, j1;
  int jj, kk;
  int best_score1, iden_no1;
  iden_no = 0;

  if ( (band_right >= len2 ) ||
       (band_left  <= -len1) ||
       (band_left  > band_right) ) return FAILED_FUNC;

  // allocate mem for score_mat[len1][len2] etc
  int band_width = band_right - band_left + 1;
  VectorInt row( band_width, 0 );
  MatrixInt score_mat( len1, row );
  MatrixInt iden_mat( len1, row );

  VectorInt & gap_array  = mat.gap_array;
  best_score = 0;
  /*
              seq2 len2 = 17            seq2 len2 = 17      seq2 len2 = 17
              01234567890123456       01234567890123456    01234567890123456
        0     xxxxxxxxxxxxxxxxx \\\\\\XXXxxxxxxxxxxxxxx    xXXXXXXXxxxxxxxxx
        1\\\\\Xxxxxxxxxxxxxxxxx  \\\\\Xxx\xxxxxxxxxxxxx    xx\xxxxx\xxxxxxxx
        2 \\\\X\xxxxxxxxxxxxxxx   \\\\Xxxx\xxxxxxxxxxxx    xxx\xxxxx\xxxxxxx
   seq1 3  \\\Xx\xxxxxxxxxxxxxx    \\\Xxxxx\xxxxxxxxxxx    xxxx\xxxxx\xxxxxx
   len1 4   \\Xxx\xxxxxxxxxxxxx     \\Xxxxxx\xxxxxxxxxx    xxxxx\xxxxx\xxxxx
   = 11 5    \Xxxx\xxxxxxxxxxxx      \Xxxxxxx\xxxxxxxxx    xxxxxx\xxxxx\xxxx
        6     Xxxxx\xxxxxxxxxxx       Xxxxxxxx\xxxxxxxx    xxxxxxx\xxxxx\xxx
        7     x\xxxx\xxxxxxxxxx       x\xxxxxxx\xxxxxxx    xxxxxxxx\xxxxx\xx
        8     xx\xxxx\xxxxxxxxx       xx\xxxxxxx\xxxxxx    xxxxxxxxx\xxxxx\x
        9     xxx\xxxx\xxxxxxxx       xxx\xxxxxxx\xxxxx    xxxxxxxxxx\xxxxx\
        0     xxxx\xxxx\xxxxxxx       xxxx\xxxxxxx\xxxx    xxxxxxxxxxx\xxxxx
                  band_left < 0           band_left < 0        band_left >=0
                  band_right < 0          band_right >=0       band_right >=0
     init score_mat, and iden_mat (place with upper 'X')
   */

  if (band_left < 0) {  //set score to left border of the matrix within band
    int tband = (band_right < 0) ? band_right : 0;
    //for (k=band_left; k<tband; k++) {
    for (k=band_left; k<=tband; k++) { // fixed on 2006 11 14
      i = -k;
      j1 = k-band_left;
      if ( ( score_mat[i][j1] = mat.matrix[iseq1[i]][iseq2[0]] ) > best_score) 
        best_score = score_mat[i][j1];
      iden_mat[i][j1] = (iseq1[i] == iseq2[0]) ? 1 : 0;
    }
  }

  if (band_right >=0) { //set score to top border of the matrix within band
    int tband = (band_left > 0) ? band_left : 0;
    for (i=0,j=tband; j<=band_right; j++) {
      j1 = j-band_left;
      if ( ( score_mat[i][j1] = mat.matrix[iseq1[i]][iseq2[j]] ) > best_score)
        best_score = score_mat[i][j1];
      iden_mat[i][j1] = (iseq1[i] == iseq2[j]) ? 1 : 0;
    }
  }

  for (i=1; i<len1; i++) {
    for (j1=0; j1<band_width; j1++) {
      j = j1+i+band_left;
      if ( j<1 ) continue;
      if ( j>=len2) continue;

      int sij = mat.matrix[iseq1[i]][iseq2[j]];
      int iden_ij = (iseq1[i] == iseq2[j] ) ? 1 : 0;
      int s1, k0;

      // from (i-1,j-1)
      if ( (best_score1 = score_mat[i-1][j1] )> 0 ) {
        iden_no1 = iden_mat[i-1][j1];
      }
      else {
        best_score1 = 0;
        iden_no1 = 0;
      }

      // from last row
      VectorInt & mat_row = score_mat[i-1];
      k0 = (-band_left+1-i > 0) ? -band_left+1-i : 0;
      for (k=j1-1, kk=0; k>=k0; k--, kk++) {
        if ( (s1 = mat_row[k] + gap_array[kk] ) > best_score1 ){
           best_score1 = s1;
           iden_no1 = iden_mat[i-1][k];
        }
      }

      k0 = (j-band_right-1 > 0) ? j-band_right-1 : 0;
      for(k=i-2, jj=j1+1,kk=0; k>=k0; k--,kk++,jj++) {
        if ( (s1 = score_mat[k][jj] + gap_array[kk] ) > best_score1 ){
           best_score1 = s1;
           iden_no1 = iden_mat[k][jj];
        }
      }

      best_score1 += sij;
      iden_no1    += iden_ij;
      score_mat[i][j1] = best_score1;
      iden_mat[i][j1]  = iden_no1;

      if ( best_score1 > best_score ) {
        best_score = best_score1;
        iden_no = iden_no1;
      }
    } // END for (j=1; j<len2; j++)
  } // END for (i=1; i<len1; i++)

  return OK_FUNC;
} // END int local_band_align


////local alignment of two sequence within a diag band
////for band 0 means direction (0,0) -> (1,1)
////         1 means direction (0,1) -> (1,2)
////        -1 means direction (1,0) -> (2,1)
////iseq len are integer sequence and its length,
////mat is matrix, return ALN_PAIR class
////copied from local_band_align, but also return alignment position
int local_band_align2(char iseq1[], char iseq2[], int len1, int len2,
                     ScoreMatrix &mat, int &best_score, int &iden_no,
                     int band_left, int band_right, 
                     int &from1, int &end1, int &from2, int &end2, int &alnln) {
  int i, j, k, j1;
  int jj, kk;
  int best_score1, iden_no1;
  int best_from1, best_from2, best_alnln;
  iden_no = 0; from1=0; from2=0;

  if ( (band_right >= len2 ) ||
       (band_left  <= -len1) ||
       (band_left  > band_right) ) return FAILED_FUNC;

  // allocate mem for score_mat[len1][len2] etc
  int band_width = band_right - band_left + 1;
  VectorInt row( band_width, 0 );
  MatrixInt score_mat( len1, row );
  MatrixInt iden_mat( len1, row );
  MatrixInt from1_mat( len1, row );
  MatrixInt from2_mat( len1, row );
  MatrixInt alnln_mat( len1, row );

  VectorInt & gap_array  = mat.gap_array;
  best_score = 0;

  if (band_left < 0) {  //set score to left border of the matrix within band
    int tband = (band_right < 0) ? band_right : 0;
    //for (k=band_left; k<tband; k++) {
    for (k=band_left; k<=tband; k++) { // fixed on 2006 11 14
      i = -k;
      j1 = k-band_left;
      if ( ( score_mat[i][j1] = mat.matrix[iseq1[i]][iseq2[0]] ) > best_score) {
        best_score = score_mat[i][j1];
        from1 = i; from2 = 0; end1 = i; end2 = 0; alnln = 1;
      }
      iden_mat[i][j1] = (iseq1[i] == iseq2[0]) ? 1 : 0;
      from1_mat[i][j1] = i;
      from2_mat[i][j1] = 0;
      alnln_mat[i][j1] = 1;
    }
  }

  if (band_right >=0) { //set score to top border of the matrix within band
    int tband = (band_left > 0) ? band_left : 0;
    for (i=0,j=tband; j<=band_right; j++) {
      j1 = j-band_left;
      if ( ( score_mat[i][j1] = mat.matrix[iseq1[i]][iseq2[j]] ) > best_score) {
        best_score = score_mat[i][j1];
        from1 = i; from2 = j; end1 = i; end2 = j; alnln = 0;
      }
      iden_mat[i][j1] = (iseq1[i] == iseq2[j]) ? 1 : 0;
      from1_mat[i][j1] = i;
      from2_mat[i][j1] = j;
      alnln_mat[i][j1] = 1;
    }
  }

  for (i=1; i<len1; i++) {
    for (j1=0; j1<band_width; j1++) {
      j = j1+i+band_left;
      if ( j<1 ) continue;
      if ( j>=len2) continue;

      int sij = mat.matrix[iseq1[i]][iseq2[j]];
      int iden_ij = (iseq1[i] == iseq2[j] ) ? 1 : 0;
      int s1, k0;

      // from (i-1,j-1)
      if ( (best_score1 = score_mat[i-1][j1] )> 0 ) {
        iden_no1 = iden_mat[i-1][j1];
        best_from1 = from1_mat[i-1][j1];
        best_from2 = from2_mat[i-1][j1];
        best_alnln = alnln_mat[i-1][j1] + 1;
      }
      else {
        best_score1 = 0;
        iden_no1 = 0;
        best_from1 = i;
        best_from2 = j;
        best_alnln = 1;
      }

      // from last row
      VectorInt & mat_row = score_mat[i-1];
      k0 = (-band_left+1-i > 0) ? -band_left+1-i : 0;
      for (k=j1-1, kk=0; k>=k0; k--, kk++) {
        if ( (s1 = mat_row[k] + gap_array[kk] ) > best_score1 ){
           best_score1 = s1;
           iden_no1 = iden_mat[i-1][k];
           best_from1 = from1_mat[i-1][k];
           best_from2 = from2_mat[i-1][k];
           best_alnln = alnln_mat[i-1][k]+kk+2;
        }
      }

      k0 = (j-band_right-1 > 0) ? j-band_right-1 : 0;
      for(k=i-2, jj=j1+1,kk=0; k>=k0; k--,kk++,jj++) {
        if ( (s1 = score_mat[k][jj] + gap_array[kk] ) > best_score1 ){
           best_score1 = s1;
           iden_no1 = iden_mat[k][jj];
           best_from1 = from1_mat[k][jj];
           best_from2 = from2_mat[k][jj];
           best_alnln = alnln_mat[k][jj]+kk+2;
        }
      }

      best_score1 += sij;
      iden_no1    += iden_ij;
      score_mat[i][j1] = best_score1;
      iden_mat[i][j1]  = iden_no1;
      from1_mat[i][j1] = best_from1;
      from2_mat[i][j1] = best_from2;
      alnln_mat[i][j1] = best_alnln;
      if ( best_score1 > best_score ) {
        best_score = best_score1;
        iden_no = iden_no1;
        end1 = i; end2 = j;
        from1 = best_from1; from2 = best_from2; alnln = best_alnln;
      }
    } // END for (j=1; j<len2; j++)
  } // END for (i=1; i<len1; i++)

  return OK_FUNC;
} // END int local_band_align2


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
  1,               // A
 -2, 1,            // C
 -2,-2, 1,         // G
 -2,-2,-2, 1,      // T
 -2,-2,-2, 1, 1,   // U
 -2,-2,-2,-2,-2, 1 // N
//A  C  G  T  U  N
//0  1  2  3  3  4
};

void setaa_to_na() {
  int i;
  for (i=0; i<26; i++) aa2idx[i]   = na2idx[i];
} // END void setaa_to_na


int setiseq(char *seq, int len)
{
  int m = 0;
  for (int i=0; i<len; i++) {
    seq[i] = aa2idx[seq[i] - 'A'];
    //m += (seq[i] == NAA1);
  }
  //printf( "NAA1 = %i %i\n", NAA1, m );
  return m;
} // END void SEQ::seq2iseq()


/////////////////
ScoreMatrix::ScoreMatrix() : gap_array( MAX_GAP )
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
  gap = gap1; ext_gap = ext_gap1;
  for (i=0; i<MAX_GAP; i++)  gap_array[i] = gap + i * ext_gap;
}

void ScoreMatrix::set_matrix(int *mat1)
{
  int i, j, k;
  k = 0;
  for ( i=0; i<MAX_AA; i++)
    for ( j=0; j<=i; j++)
      matrix[j][i] = matrix[i][j] = mat1[ k++ ];
}

void ScoreMatrix::set_to_na()
{
	set_gap( -6, -1 );
	set_matrix( BLOSUM62_na );
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
  size = 0;
  frag_count = 0;
  sequences.clear();
  for (i=0; i<NAAN; i++) indexCounts[i].Clear();
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
      ic.index = idx;
      indexCounts[j].Append( ic );
      size ++;
    }
  }
  sequences.Append( seq );
  return OK_FUNC;
}
int WordTable::AddWordCountsFrag( NVector<IndexCount> & counts, int frag, int frag_size, int repfrag )
{
}
int WordTable::AddWordCounts(int aan_no, Vector<int> & word_encodes, Vector<INTs> & word_encodes_no, int idx, bool skipN)
{
  int i, j, k;
  for (i=0; i<aan_no; i++) {
    if ( (k=word_encodes_no[i]) ) {
      j = word_encodes[i];
      if ( skipN && j<0) continue; // for those has 'N'
      indexCounts[j].Append( IndexCount( idx, k ) );
      size ++;
    }
  }
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
        indexCounts[j].Append( IndexCount( frag_count + fra, k1 ) );
        size ++;
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

int WordTable::CountWords( NVector<IndexCount> & counts, Vector<INTs> & look_and_count, bool est)
{
  int aan_no = counts.Size();
  int  j, k, j0, j1, k1;

  for (j0=0; j0<aan_no; j0++) {
    IndexCount & ec = counts[j0];
    if ( (j1=ec.count) ) {
      j = ec.index;
      if (est && j<0) continue; // if met short word has 'N'
      NVector<IndexCount> & one = indexCounts[j];
      k1 = one.Size();
      for (k=0; k<k1; k++){
        IndexCount & ic = one[k];
        look_and_count[ ic.index ] += ( j1 < ic.count ) ? j1 : ic.count ;
      }
    }
  }
                                                                                
  return OK_FUNC;
}
int WordTable::CountWords(int aan_no, Vector<int> & word_encodes,
                           Vector<INTs> & word_encodes_no, Vector<INTs> &look_and_count, bool est)
{
  int  j, k, j0, j1, k1;

  //j0 = 0;
  //if( est ) while( word_encodes[j0] <0 ) ++ j0; // if met short word has 'N'
  for (j0=0; j0<aan_no; j0++) {
    if ( (j1=word_encodes_no[j0]) ) {
      j = word_encodes[j0];
      if (est && j<0) continue; // if met short word has 'N'
      NVector<IndexCount> & one = indexCounts[j];
      k1 = one.Size();
      for (k=0; k<k1; k++){
        IndexCount & ic = one[k];
        look_and_count[ ic.index ] += ( j1 < ic.count ) ? j1 : ic.count ;
      }
    }
  }
                                                                                
  return OK_FUNC;
}

Sequence::Sequence()
{
  memset( this, 0, sizeof( Sequence ) );
}
Sequence::Sequence( const Sequence & other )
{
  int i;
  //printf( "new: %p  %p\n", this, & other );
  memcpy( this, & other, sizeof( Sequence ) );
  if( other.data ){
    size = bufsize = other.size;
    data = new char[size+1];
    //printf( "data: %p  %p\n", data, other.data );
    data[size] = 0;
    memcpy( data, other.data, size );
    //for (i=0; i<size; i++) data[i] = other.data[i];
  }
}
Sequence::~Sequence()
{
  //printf( "delete: %p\n", this );
  if( data ) delete[] data;
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
    if( size ) data[size] = 0;
  }
}

void Sequence::Swap( Sequence & other )
{
  Sequence tmp;
  memcpy( & tmp, this, sizeof( Sequence ) );
  memcpy( this, & other, sizeof( Sequence ) );
  memcpy( & other, & tmp, sizeof( Sequence ) );
  memset( & tmp, 0, sizeof( Sequence ) );
}
void Sequence::Format()
{
  int i, j=0;
  for (i=0; i<size; i++){
    char ch = data[i];
    if ( isalpha( ch ) ) data[j++] = toupper( ch );
  }
  data[j] = 0;
  size = j;
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
void Sequence::PrintInfo( int id, FILE *fin, FILE *fout, const Options & options )
{
  const char *tag = options.isEST ? "nt" : "aa";
  char *buf = new char[ MAX_DES + 1 ];
  bool print = options.print != 0;
  bool strand = options.isEST;
  int i, len = options.des_len ? options.des_len : MAX_DES;
  fseek( fin, des_begin, SEEK_SET );
  if( (len=fread( buf, 1, len, fin )) ==0 ) bomb_error( "Can not swap in sequence" );
  i = 0;
  if( buf[i] == '>' ) i += 1;
  if( buf[i] == ' ' or buf[i] == '\t' ) i += 1;
  while( i < len and ! isspace( buf[i] ) ) i += 1;
  buf[i] = 0;
#if 0
  if( options.des_len ==0 ){
    len = 0;
    while( len < MAX_DES && ! isspace( buf[len] ) ) len ++;
  }
  buf[ len ] = 0;
#endif
  fprintf( fout, "%i\t%i%s, %s...", id, size, tag, buf );
  if( identity ){
    int *c = coverage;
    fprintf( fout, " at " );
    if (print) fprintf( fout, "%i:%i:%i:%i/", c[0], c[1], c[2], c[3] );
    if (strand) fprintf( fout, "%c/", (state & IS_MINUS_STRAND) ? '-' : '+' );
    fprintf( fout, "%i%%\n", identity );
  }else{
    fprintf( fout, " *\n" );
  }
  delete []buf;
}

void SequenceDB::Read( const char *file, const Options & options )
{
  Sequence one;
  Sequence dummy;
  FILE *swap = NULL;
  FILE *fin = fopen( file, "r" );
  char *buffer = NULL;
  char *res = NULL;
  size_t swap_size = 0;
  int option_l = options.min_length;
  if( fin == NULL ) bomb_error( "Failed to open the database file" );
  if( options.store_disk ) swap = OpenTempFile( temp_dir );
  Clear();
  dummy.swap = swap;
  buffer = new char[ MAX_LINE_SIZE+1 ];

  while (not feof( fin ) || one.size) { /* do not break when the last sequence is not handled */
    if ( (res=fgets( buffer, MAX_LINE_SIZE, fin )) == NULL && one.size == 0) break;
    if (buffer[0] == '>' || (res==NULL && one.size)) {
      if ( one.size ) { // write previous record
        one.dat_length = dummy.dat_length = one.size;
        one.Format();
        one.index = dummy.index = sequences.size();
        if ( one.size > option_l ) {
          if ( swap ) {
            swap_size += one.size;
            // so that size of file < MAX_BIN_SWAP about 2GB
            if ( swap_size >= MAX_BIN_SWAP) {
              dummy.swap = swap = OpenTempFile( temp_dir );
              swap_size = one.size;
            }
            dummy.size = one.size;
            dummy.offset = ftell( swap );
            sequences.Append( new Sequence( dummy ) ); 
            setiseq( one.data, one.size );
            fwrite( one.data, 1, one.size, swap );
          }else{
            //printf( "==================\n" );
            sequences.Append( new Sequence( one ) ); 
            //printf( "------------------\n" );
            //if( sequences.size() > 10 ) break;
          }
        }
        one.size = 0;
      }
      int len = strlen( buffer );
      int len2 = len;
      while( len2 && buffer[len2-1] != '\n' ){
        if ( (res=fgets( buffer, MAX_LINE_SIZE, fin )) == NULL ) break;
        len2 = strlen( buffer );
        len += len2;
      }
      size_t offset = ftell( fin );
      one.des_begin = dummy.des_begin = offset - len;
      one.des_length = dummy.des_length = len;
    } else {
      one += buffer;
    }
  }
#if 0
  int i, n = 0;
  for(i=0; i<sequences.size(); i++) n += sequences[i].bufsize + 4;
  cout<<n<<"\t"<<sequences.capacity() * sizeof(Sequence)<<endl;
  int i;
  scanf( "%i", & i );
#endif
  delete[] buffer;
  fclose( fin );
}

void SequenceDB::SortDivide( Options & options, bool sort )
{
  int i, j, k, len;
  long long total_letter=0;
  int max_len = 0, min_len = 99999999;
  int N = sequences.size();
  for (i=0; i<N; i++) {
    Sequence *seq = sequences[i];
    len = seq->size;
    total_letter += len;
    if (len > max_len) max_len = len;
    if (len < min_len) min_len = len;
    if (seq->swap == NULL) setiseq( seq->data, len );
  }
  options.total_letters = total_letter;
  if (max_len >= 65536) 
    bomb_warning("Some seqs longer than 65536, you may define LONG_SEQ");
  cout << "longest and shortest : " << max_len << " and " << min_len << endl;
  cout << "Total letters: " << total_letter << endl;
  // END change all the NR_seq to iseq

  if( sort ){
    // **************************** Form NR_idx[], Sort them from Long to short
    int M = max_len - min_len + 1;
    Vector<int> count( M, 0 ); // count for each size = max_len - i
    Vector<int> accum( M, 0 ); // count for all size > max_len - i
    Vector<int> offset( M, 0 ); // offset from accum[i] when filling sorting
    Vector<Sequence*> sorting( N ); // TODO: use a smaller class if this consumes to much memory!

    for (i=0; i<N; i++) count[ max_len - sequences[i]->size ] ++;
    for (i=1; i<M; i++) accum[i] = accum[i-1] + count[i-1];
    for (i=0; i<N; i++){
      int len = max_len - sequences[i]->size;
      int id = accum[len] + offset[len];
      //sequences[i].index = id;
      sorting[id] = sequences[i];
      offset[len] ++;
    }
    for (i=0; i<N; i++) sequences[i] = sorting[i];
    cout << "Sequences have been sorted" << endl;
    // END sort them from long to short
  }
}// END sort_seqs_divide_segs

void SequenceDB::DivideSave( const char *db, const char *newdb, int n, const Options & options )
{
  if( n == 0 or sequences.size() ==0 ) return;

  size_t max_seg = options.total_letters / n + sequences[0]->size;
  if( max_seg >= MAX_BIN_SWAP ) max_seg = (size_t) MAX_BIN_SWAP;

  FILE *fin = fopen( db, "r" );
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

    count = (seq->des_length + seq->dat_length) / MAX_LINE_SIZE;
    rest = (seq->des_length + seq->dat_length) % MAX_LINE_SIZE;
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
  FILE *fin = fopen( db, "r" );
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

    count = (seq->des_length + seq->dat_length) / MAX_LINE_SIZE;
    rest = (seq->des_length + seq->dat_length) % MAX_LINE_SIZE;
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
void SequenceDB::WriteExtra1D( const Options & options )
{
  string db_clstr = options.output + ".clstr";
  string db_clstr_bak = options.output + ".bak.clstr";
  int i, k, N = sequences.size();
  vector<long long> sorting( N );
  for (i=0; i<N; i++) sorting[i] = ((long long)sequences[i]->index << 32) | i;
  std::sort( sorting.begin(), sorting.end() );

  FILE *fin = fopen( options.input.c_str(), "r" );
  FILE *fout = fopen( db_clstr_bak.c_str(), "w+" );
  for (i=0; i<N; i++) {
    Sequence *seq = sequences[ sorting[i] & 0xffffffff ];
    seq->PrintInfo( seq->cluster_id, fin, fout, options );
  }
  fclose( fout );

  cout << "writing clustering information" << endl;
  int M = rep_seqs.size();
  Vector<Vector<int> > clusters( M );
  for (i=0; i<N; i++){
    int k = sorting[i] & 0xffffffff;
    int id = sequences[k]->cluster_id;
    clusters[id].Append( k );
  }

  fout = fopen( db_clstr.c_str(), "w+" );
  for (i=0; i<M; i++) {
    fprintf( fout, ">Cluster %i\n", i );
    for (k=0; k<(int)clusters[i].size(); k++)
      sequences[ clusters[i][k] ]->PrintInfo( k, fin, fout, options );
  }
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

  FILE *fin = fopen( options.input.c_str(), "r" );
  FILE *fin2 = fopen( options.input2.c_str(), "r" );
  FILE *fout = fopen( db_clstr_bak.c_str(), "w+" );
  for (i=0; i<N; i++) {
    Sequence *seq = other.sequences[ sorting[i] & 0xffffffff ];
    seq->PrintInfo( seq->cluster_id, fin, fout, options );
  }
  for (i=0; i<N2; i++) {
    Sequence *seq = sequences[i];
    if( seq->state & IS_REDUNDANT ) seq->PrintInfo( seq->cluster_id, fin2, fout, options );
  }
  fclose( fout );

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
    seq->PrintInfo( 0, fin, fout, options );
    for (k=0; k<(int)clusters[i].size(); k++)
      sequences[ clusters[i][k] ]->PrintInfo( k+1, fin2, fout, options );
  }
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

void show_cpu_time(tms &CPU_begin, tms &CPU_end) {
  int  ClockTicksPerSecond, total_seconds;
//  ClockTicksPerSecond = (int)sysconf(_SC_CLK_TCK);
//  ClockTicksPerSecond = (int)(100);
  ClockTicksPerSecond = CLOCK_TICKS;

  total_seconds = (CPU_end.tms_utime - CPU_begin.tms_utime) 
                  / ClockTicksPerSecond;

  cout << "Total CPU time " << total_seconds << endl;
} // END  show_current_cpu_time

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
                     int naa_stat[5][61][4], int NAA) {
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
                     int naa_stat[5][61][4], int NAA, int iden) {
  double cluster_thd;
  cluster_thd = ((double)(iden)) / 100.0;
  if (cluster_thd > 1.0) cluster_thd = 1.00;

  double aa1_t, aa2_t, aan_t;
  cal_aax_cutoff(aa1_t, aa2_t, aan_t, cluster_thd, tolerance, naa_stat_start_percent,
                 naa_stat, NAA);
  if (aa1_t > aa1_cutoff) aa1_cutoff = aa1_t;
  if (aa2_t > aa2_cutoff) aa2_cutoff = aa2_t;
  if (aan_t > aan_cutoff) aan_cutoff = aan_t;
  return;  
} // END update_aax_cutoff

void WorkingParam::ComputeRequiredBases( int NAA, int ss )
{
  required_aa1 = int (aa1_cutoff* (double) len_eff);
  required_aas = (aa1_cutoff > 0.95) ?
    len_eff-ss  +1-(len_eff-required_aa1)*ss   :
    int (aas_cutoff* (double) len_eff);
  required_aan = (aa1_cutoff > 0.95) ?
    len_eff-NAA+1-(len_eff-required_aa1)*NAA :
    int (aan_cutoff* (double) len_eff);
}
int WorkingBuffer::EncodeWords( Sequence *seq, int NAA, bool est )
{
  char *seqi = seq->data;
  int len = seq->size;
  // check_word_encodes
  int aan_no = len - NAA + 1;
  int i, j, k, i0, i1, k1;
  int skip = 0;
  for (j=0; j<aan_no; j++) {
    word_encodes[j] = 0;
    for (k=0, k1=NAA-1; k<NAA; k++, k1--) word_encodes[j] += seqi[j+k] * NAAN_array[k1];
    word_encodes_backup[j] = word_encodes[j];
  }

  if( est ){
    for (j=0; j<len; j++){
      if ( seqi[j] == 4 ) {                      // here N is 4
        i0 = (j-NAA+1 > 0)      ? j-NAA+1 : 0;
        for (i=i0; i<=j; i++) word_encodes[i]=-1;
      }
    }
    for (j=0; j<len; j++) skip += (word_encodes[j] == -1);
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
    if ((seqi[j1]==4) || (seqi[j1+1]==4) || (seqi[j1+2]==4) || (seqi[j1+3]==4)) continue; //skip N
    c22 = seqi[j1]*NAA3 + seqi[j1+1]*NAA2 + seqi[j1+2]*NAA1 + seqi[j1+3];
    taap[c22]++;
  }
  for (sk=0,mm=0; sk<NAA4; sk++) {
    aap_begin[sk] = mm; mm+=taap[sk]; taap[sk] = 0;
  }
  for (j1=0; j1<len1; j1++) {
    if ((seqi[j1]==4) || (seqi[j1+1]==4) || (seqi[j1+2]==4) || (seqi[j1+3]==4)) continue; //skip N
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
        seq->fragment = table.frag_count;
        table.AddWordCountsFrag( aan_no, buffer.word_encodes_backup, 
            buffer.word_encodes_no, frg1, frag_size );
        // the first fragment for the representative
      }else{
        table.AddWordCounts(aan_no, buffer.word_encodes, buffer.word_encodes_no, table.sequences.size(), options.isEST);
      }
      table.sequences.Append( seq );
    }
  }
  if ( (id+1) % 100 == 0 ) {
    int size = rep_seqs.size();
    printf( "." );
    fflush( stdout );
    if ( (id+1) % 1000 == 0 ) printf( "\r..........%9i  finished  %9i  clusters\n", id+1, size );
  }
}
#include<assert.h>
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
  valarray<int>  indices(T);

  //printf( "%li\n", options.mem_limit );

  if (frag_size){ 
    frag_no = 0;
    for (i=0; i<seq_no; i++) frag_no += (sequences[i]->size - NAA) / frag_size + 1;
  }

  if( not options.isEST )
    cal_aax_cutoff(aa1_cutoff, aas_cutoff, aan_cutoff, options.cluster_thd,
        options.tolerance, naa_stat_start_percent, naa_stat, NAA);

  WorkingParam param( aa1_cutoff, aas_cutoff, aan_cutoff );
  WorkingBuffer buffer( frag_no, options.isEST );

  Vector<WorkingParam> params(T);
  Vector<WorkingBuffer> buffers(T);
  for(i=0; i<T; i++){
    params[i].Set( aa1_cutoff, aas_cutoff, aan_cutoff );
    buffers[i].Set( frag_no, options.isEST );
  }
  omp_set_num_threads(T);

  // word_table as self comparing table and table buffer:
  WordTable word_table( options.NAA, NAAN );

  WordTable last_table( options.NAA, NAAN );

  size_t mem_limit = options.max_memory / sizeof(IndexCount);
  size_t mem_limit2 = mem_limit / 50;
  int N = sequences.size();
  int K = N - 100 * T;

  size_t total_letters = options.total_letters;
  letters = 0;

  for(i=0; i<N; ){
    for(j=0; j<T; j++) total_letters -= letters[j];
    letters = 0;

    int start = i;
    int m = i;
    size_t sum = 0;
    size_t lim = total_letters / T; // preferred Self-Comparing Block(SCB) size
    assert( total_letters );
    if( lim > mem_limit ) lim = mem_limit; // SCB size has upper limit
    if( i ==0 ) lim /= 8; // first SCB with small size
    if( lim < mem_limit2 ) lim = (lim + mem_limit2) / 2; // SCB size has lower limit
    while( m < N && sum < lim ){
      Sequence *seq = sequences[m];
      if( ! (seq->state & IS_REDUNDANT) ){
        if ( options.store_disk ) seq->SwapIn();
        sum += seq->size;
      }
      m ++;
    }
    if( m > N ) m = N;
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
      int JN = N;
      int p0 = 0;
      int min = last_table.sequences[ last_table.sequences.size()-1 ]->size;
      int m0 = m;
      indices = 0;
      #pragma omp parallel for schedule( dynamic, 1 )
      for(int j=m-1; j<JN; j++){
        if( j+1 == N ){
          //printf( "stoping\n" );
          //#pragma omp atomic
          may_stop = 1;
        }
        int tid = omp_get_thread_num();
        if( j == (m0-1) ){ // use m0 to avoid other iterations satisfying the condition:
          indices[tid] = N; // just in case if this tid is not zero
          for(int ks=i; ks<m; ks++){
            if( ks == (m-1) ){
              int done = N;
              for(int t=0; t<T; t++) if( indices[t] < done ) done = indices[t];
              //printf( "K = %9i,  done = %9i\n", K, done );
              if( done < K && done > m){
                int Q = m + (N-done)/T + 1000;
                if( Q > done ) Q = done;
                //printf( "update m: %9i\n", m );
                while( m < Q && sum < mem_limit ){
                  Sequence *seq = sequences[m];
                  if( ! (seq->state & IS_REDUNDANT) ){
                    if ( options.store_disk ){
                      #pragma omp critical
                      seq->SwapIn();
                    }
                  }
                  sum += seq->size;
                  m ++;
                }
                if( m > done ) m = done;
                printf( "\r---------- extending the current cycle to sequence  %9i (others=%2i%)\n", m, (int)(100*done/N) );
              }
            }
            Sequence *seq = sequences[ks];
            i = ks + 1;
            if (seq->state & IS_REDUNDANT) continue;
            ClusterOne( seq, ks, word_table, params[tid], buffers[tid], options );
            total_letters -= seq->size;
            if ( options.store_disk && (seq->state & IS_REDUNDANT) ) seq->SwapOut();
            if( may_stop ) break;
            //if( ks % 100 ==0 ) printf( "done = %i\n", done );
            //if( may_stop && (ks+10*T) < m ) break;
          }
          self_stop = 1;
        }else{
          Sequence *seq = sequences[j];
          if (seq->state & IS_REDUNDANT) continue;
          if ( options.store_disk ){
            #pragma omp critical
            seq->SwapIn();
          }
          CheckOne( seq, last_table, params[tid], buffers[tid], options );
          letters[tid] -= (seq->state & IS_REDUNDANT) ? seq->size : 0;
          if ( options.store_disk && (seq->state & IS_REDUNDANT) ) seq->SwapOut();
          if( min > params[tid].len_upper_bound ){
            may_stop = 1;
            #pragma omp critical
            JN = j;
            continue;
          }
          if( self_stop && tid ==1 ){
            int p = (int)(100*j/N);
            if( p > p0 ){ // print only if the percentage changed
              printf( "\r%2i%%", p );
              fflush( stdout );
              p0 = p;
            }
          }

        }
        indices[tid] = j;
      }
    }
    if( i == start || m == N || word_table.size < mem_limit2 ){
      for(k=i; k<m; ){
        int kk, mm = k, sum = 0;
        while( mm < m && sum < 1E5 ){
          if( ! (sequences[mm]->state & IS_REDUNDANT) ) sum += sequences[mm]->size;
          mm += T;
        }
        if( mm < k + 100*T ) mm = k + 100*T;
        if( mm > m ) mm = m;
        #pragma omp parallel for schedule( static, 1 )
        for(kk=k; kk<mm; kk++){
          Sequence *seq = sequences[kk];
          if (seq->state & IS_REDUNDANT) continue;
          int tid = omp_get_thread_num();
          CheckOne( seq, word_table, params[tid], buffers[tid], options );
          letters[tid] -= (seq->state & IS_REDUNDANT) ? seq->size : 0;
          if ( options.store_disk && (seq->state & IS_REDUNDANT) ) seq->SwapOut();
        }
        for(int ks=k; ks<mm; ks++){
          Sequence *seq = sequences[ks];
          if (seq->state & IS_REDUNDANT) continue;
          ClusterOne( seq, ks, word_table, param, buffer, options );
          total_letters -= seq->size;
          if ( options.store_disk && (seq->state & IS_REDUNDANT) ) seq->SwapOut();
        }
        k = mm;
      }
      i = m;
    }else if( i < m ){
      printf( "\r---------- %6i remaining sequences to the next cycle\n", m-i );
    }
    last_table.Clear();
    last_table.sequences.swap( word_table.sequences );
    last_table.indexCounts.swap( word_table.indexCounts );
    last_table.size = word_table.size;
    word_table.size = 0;
  }
  printf( "\n%9li  finished  %9li  clusters\n", sequences.size(), rep_seqs.size() );
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
  Vector<INTs> & look_and_count = buf.look_and_count;
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

  param.ControlShortCoverage( len, options );
  param.ComputeRequiredBases( options.NAA, 2 );

  buf.EncodeWords( seq, options.NAA, false );

  // if minimal alignment length > len, return
  // I can not return earlier, because I need to calc the word_encodes etc
  if (options.min_control>len) return 0; // return flag=0

  // lookup_aan
  if( frag_size )
    for (j=0; j<table.frag_count; j++) look_and_count[j]=0;
  else
    for (j=0; j<S; j++) look_and_count[j]=0;

  int aan_no = len - options.NAA + 1;
  table.CountWords(aan_no, word_encodes, word_encodes_no, look_and_count);

  // contained_in_old_lib()
  int len_upper_bound = param.len_upper_bound;
  int len_lower_bound = param.len_lower_bound;
  int band_left, band_right, best_score, band_width1, best_sum, len2, alnln, len_eff1;
  int tiden_no;
  int talign_info[5];
  int best1, sum;
  INTs *lookptr;
  char *seqj;
  int frg2 = frag_size ? (len - NAA + options.band_width ) / frag_size + 1 + 1 : 0;
  int lens;
  int has_aa2 = 0;
  for (j=0; j<S; j++) {
    if ( frag_size ==0 && look_and_count[j] < required_aan ) continue;

    Sequence *rep = table.sequences[j];
    len2 = rep->size;
    if (len2 > len_upper_bound ) continue;
    if (options.has2D && len2 < len_lower_bound ) continue;
    if ( frag_size ){
      k = (len2 - NAA) / frag_size + 1;
      lookptr = & look_and_count[ rep->fragment ];

      if ( frg2 >= k ) {
        best1=0;
        for (j1=0; j1<k; j1++) best1 += lookptr[j1];
      } else {
        sum = 0;
        for (j1=0; j1<frg2; j1++) sum += lookptr[j1];
        best1 = sum;
        for (j1=frg2; j1<k; j1++) {
          sum += lookptr[j1] - lookptr[j1-frg2];
          if (sum > best1) best1 = sum;
        }
      }

      if ( best1 < required_aan ) continue;
    }

    param.ControlLongCoverage( len2, options );

    if ( has_aa2 == 0 )  { // calculate AAP array
      buf.ComputeAAP( seqi, seq->size );
      has_aa2 = 1;
    }
    seqj = rep->data; //NR_seq[NR90_idx[j]];

    band_width1 = (options.band_width < len+len2-2 ) ? options.band_width : len+len2-2;
    diag_test_aapn(NAA1, seqj, len, len2, buf, best_sum,
        band_width1, band_left, band_right, required_aa1-1);
    if ( best_sum < required_aa2 ) continue;

    if (options.print || aln_cover_flag) //return overlap region
      local_band_align2(seqi, seqj, len, len2, mat,
                        best_score, tiden_no, band_left, band_right,
                        talign_info[1],talign_info[2],
                        talign_info[3],talign_info[4],alnln);
    else
      local_band_align(seqi, seqj, len, len2, mat,
                             best_score, tiden_no, band_left, band_right);
    if ( tiden_no < required_aa1 ) continue;
    lens = len;
    if( options.has2D && len > len2 ) lens = len2;
    len_eff1 = (options.global_identity == 0) ? alnln : lens;
    tiden_no = (tiden_no * 100) / len_eff1;
    if (tiden_no < options.cluster_thd100) continue;
    if (tiden_no <= seq->identity) continue; // existing iden_no
    if (aln_cover_flag) {
      if ( talign_info[4]-talign_info[3]+1 < min_aln_lenL) continue;
      if ( talign_info[2]-talign_info[1]+1 < min_aln_lenS) continue;
    }
    if( options.has2D ) seq->state |= IS_REDUNDANT ;
    flag = 1; seq->identity = tiden_no; seq->cluster_id = rep->cluster_id;
    seq->coverage[0] = talign_info[1] +1;
    seq->coverage[1] = talign_info[2] +1;
    seq->coverage[2] = talign_info[3] +1;
    seq->coverage[3] = talign_info[4] +1;
    if (not options.cluster_best) break;
    update_aax_cutoff(aa1_cutoff, aa2_cutoff, aan_cutoff,
        options.tolerance, naa_stat_start_percent, naa_stat, NAA, tiden_no);
    param.ComputeRequiredBases( options.NAA, 2 );
  }
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
  Vector<INTs> & look_and_count = buf.look_and_count;
  Vector<INTs> & word_encodes_no = buf.word_encodes_no;
  Vector<INTs> & aap_list = buf.aap_list;
  Vector<INTs> & aap_begin = buf.aap_begin;
  Vector<int>  & word_encodes = buf.word_encodes;
  Vector<int>  & taap = buf.taap;
  Vector<int> & aan_list_comp = buf.aan_list_comp;
  char *seqi_comp = buf.seqi_comp;

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

  param.ControlShortCoverage( len, options );
  param.ComputeRequiredBases( options.NAA, 4 );
  int skip = buf.EncodeWords( seq, options.NAA, true );
  required_aan -= skip;

  // if minimal alignment length > len, return
  // I can not return earlier, because I need to calc the word_encodes etc
  if (options.min_control>len) return 0; // return flag=0

  int aan_no = len - options.NAA + 1;

  // contained_in_old_lib()
  int len_upper_bound = param.len_upper_bound;
  int len_lower_bound = param.len_lower_bound;
  int band_left, band_right, best_score, band_width1, best_sum, len2, alnln, len_eff1;
  int tiden_no;
  int talign_info[5];
  int j0, comp;
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
    for (j=0; j<S; j++) look_and_count[j]=0;
    if( comp ){
      table.CountWords(aan_no, aan_list_comp, word_encodes_no, look_and_count, true);
    }else{
      table.CountWords(aan_no, word_encodes, word_encodes_no, look_and_count, true);
    }
    for (j=0; j<S; j++){
      //printf( "1: %4i %4i %3i\n", look_and_count[j], required_aan, seq->xletter );
      if ( look_and_count[j] < required_aan ) continue;

      Sequence *rep = table.sequences[j];
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
          band_width1, band_left, band_right, required_aa1-3);
      //printf( "a:  %5i  %5i\n", best_sum, required_aas );
      if ( best_sum < required_aas ) continue;

      if (options.print || aln_cover_flag){ //return overlap region
        local_band_align2(seqi, seqj, len, len2, mat,
            best_score, tiden_no, band_left, band_right,
            talign_info[1],talign_info[2],
            talign_info[3],talign_info[4],alnln);
        if( comp ){
          talign_info[1] = len - talign_info[1] - 1;
          talign_info[2] = len - talign_info[2] - 1;
        }
      }else{
        //printf( "b:  %5i  %5i\n", band_left, band_right );
        local_band_align(seqi, seqj, len, len2, mat,
            best_score, tiden_no, band_left, band_right);
      }
      if ( tiden_no < required_aa1 ) continue;
      len_eff1 = (options.global_identity == 0) ? alnln : len;
      tiden_no = (tiden_no * 100) / len_eff1;
      if (tiden_no < options.cluster_thd100) continue;
      if (options.cluster_best and tiden_no < seq->identity) continue; // existing iden_no
      if (aln_cover_flag) {
        if ( talign_info[4]-talign_info[3]+1 < min_aln_lenL) continue;
        if( comp ){
          if ( talign_info[1]-talign_info[2]+1 < min_aln_lenS) continue;
        }else{
          if ( talign_info[2]-talign_info[1]+1 < min_aln_lenS) continue;
        }
      }
      if( options.cluster_best and tiden_no == seq->identity and rep->cluster_id >= seq->cluster_id ) continue;
      if( (not options.cluster_best) and flag !=0 and rep->cluster_id >= seq->cluster_id ) continue;
      flag = comp ? -1 : 1;
      seq->identity = tiden_no;
      seq->cluster_id = rep->cluster_id;
      seq->coverage[0] = talign_info[1] +1;
      seq->coverage[1] = talign_info[2] +1;
      seq->coverage[2] = talign_info[3] +1;
      seq->coverage[3] = talign_info[4] +1;
      if (not options.cluster_best) break;
    }
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
  WorkingBuffer buffer( frag_no, options.isEST );

  WordTable word_table( options.NAA, NAAN );

  size_t mem_limit = options.max_memory / sizeof(IndexCount);
  int N = sequences.size();

  size_t total_letters = options.total_letters;

  for(i=0; i<N; ){
    size_t sum = 0;
    int m = i;
    while( m < N && sum < mem_limit ){
      Sequence *seq = sequences[m];
      if( ! (seq->state & IS_REDUNDANT) ){
        if ( options.store_disk ) seq->SwapIn();
        sum += seq->size;
      }
      m ++;
    }
    if( m > N ) m = N;
    printf( "\rcomparing sequences from  %9i  to  %9i\n", i, m );
    for(int ks=i; ks<m; ks++){
      Sequence *seq = sequences[ks];
      if (seq->state & IS_REDUNDANT) continue;
      ClusterOne( seq, ks, word_table, param, buffer, options );
      total_letters -= seq->size;
      if ( options.store_disk && (seq->state & IS_REDUNDANT) ) seq->SwapOut();
    }
    i = m;
    if( word_table.size == 0 ) continue;
    int p0 = 0;
    for(int j=m; j<N; j++){
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
      int p = (int)((100*j)/N);
      if( p > p0 ){ // print only if the percentage changed
        printf( "\r%2i%%", p );
        fflush( stdout );
        p0 = p;
      }
    }
    //if( i && i < m ) printf( "\r---------- %6i remaining sequences to the next cycle\n", m-i );
    word_table.Clear();
  }
  printf( "\n%9li  finished  %9li  clusters\n", sequences.size(), rep_seqs.size() );
  temp_files.Clear();

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

  WorkingParam param( aa1_cutoff, aas_cutoff, aan_cutoff );
  WorkingBuffer buffer( other.sequences.size(), options.isEST );

  size_t mem_limit = options.max_memory / sizeof(IndexCount);
  int N = other.sequences.size();
  int M = sequences.size();

  int T = options.threads;
  valarray<size_t>  counts(T);
  Vector<WorkingParam> params(T);
  Vector<WorkingBuffer> buffers(T);
  for(i=0; i<T; i++){
    params[i].Set( aa1_cutoff, aas_cutoff, aan_cutoff );
    buffers[i].Set( N, options.isEST );
  }
  if( T >1 ) omp_set_num_threads(T);


  size_t total_letters = 0;
  for(i=0; i<N; i++) total_letters += other.sequences[i]->size;

  WordTable word_table( options.NAA, NAAN );

  for(i=0; i<N; ){
    size_t sum = 0;
    int m = i;
    while( m < N && sum < mem_limit ){
      Sequence *seq = other.sequences[m];
      if( ! (seq->state & IS_REDUNDANT) ){
        if ( options.store_disk ) seq->SwapIn();
        sum += seq->size;
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
      if ( (ks+1) % 100 == 0 ) {
        printf( "." );
        fflush( stdout );
        if ( (ks+1) % 1000 == 0 ) printf( "%9i  finished\n", ks+1 );
      }  
    }
    int p0 = 0;
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
        int p = (int)((100*j)/N);
        if( p > p0 ){ // print only if the percentage changed
          printf( "\r%2i%%", p );
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
        int p = (int)((100*j)/N);
        if( p > p0 ){ // print only if the percentage changed
          printf( "\r%2i%%", p );
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
      if ( seqi[j] == 4 ) {                      // here N is 4
        i0 = (j-NAA+1 > 0)      ? j-NAA+1 : 0;
        for (i=i0; i<=j; i++) aan_list[i]=-1;
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

