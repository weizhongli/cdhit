// =============================================================================
// CD-HIT
// http://bioinformatics.burnham-inst.org/cd-hi
// 
// program written by 
//                                      Weizhong Li
//                                      UCSD, San Diego Supercomputer Center
//                                      La Jolla, CA, 92093
//                                      Email liwz@sdsc.edu
//                 at
//                                      Adam Godzik's lab
//                                      The Burnham Institute
//                                      La Jolla, CA, 92037
//                                      Email adam@burnham-inst.org
// =============================================================================

#include "cd-hi.h"
#include "cd-hi-init.h"

int db_read_and_div(ifstream &in1, int SEG_no, int option_l,
                    char *NR_seg, char db_swap[MAX_SEG][MAX_FILE_NAME]);
int sort_seqs_divide_segs_div(int NR_no, int *NR_len, int *NR_idx, char *NR_seg,
                           int SEG_no, char db_swap[MAX_SEG][MAX_FILE_NAME],
                           char db_out[]);
int db_read_in_lenf(ifstream &in1,  int & NR_no, int *NR_len);
int db_read_and_divf(ifstream &in1, int SEG_no,
                    char *NR_seg, char db_swap[MAX_SEG][MAX_FILE_NAME]);

////////////////////////////////////  MAIN /////////////////////////////////////
int main(int argc, char **argv) {
  int i, j, k, i1, j1, k1, i0, j0, k0, sg_i, sg_j;
  int si, sj, sk;
  char db_in[MAX_FILE_NAME];
  char db_out[MAX_FILE_NAME];

  times(&CPU_begin);

  // ***********************************    parse command line and open file
  if (argc < 5) print_usage_div(argv[0]);

  DB_no = 30000000;
  for (i=1; i<argc; i++) {
    if      (strcmp(argv[i], "-i"    )==0) strncpy(db_in,  argv[++i], MAX_FILE_NAME-1);
    else if (strcmp(argv[i], "-o"    )==0) strncpy(db_out, argv[++i], MAX_FILE_NAME-1);
    else if (strcmp(argv[i], "-div"  )==0) SEG_no    = atoi(argv[++i]);
    else if (strcmp(argv[i], "-dbmax")==0) DB_no     = atoi(argv[++i]);
    else                                   print_usage_div(argv[0]);
  }

  if ((NR_len      = new int   [DB_no]) == NULL) bomb_error("Memory");
  if ((NR_idx      = new int   [DB_no]) == NULL) bomb_error("Memory");
  if ((NR_seg      = new char  [DB_no]) == NULL) bomb_error("Memory");

  ifstream in1(db_in); if (! in1) bomb_error("Can not open file", db_in); 
  db_read_in_lenf(in1, NR_no, NR_len);

  in1.close(); 
  cout << "total seq: " << NR_no << endl;

  sort_seqs_divide_segs_div(NR_no, NR_len, NR_idx, NR_seg,
                        SEG_no, db_swap, db_out);

  ifstream in1b(db_in); if (! in1b) bomb_error("Can not open file", db_in); 
  db_read_and_divf(in1b, SEG_no, NR_seg, db_swap);
  in1b.close();

  cout << "program completed !" << endl << endl;

  times(&CPU_end);
  show_cpu_time(CPU_begin, CPU_end);
  return 0;
} // END int main

///////////////////////FUNCTION of common tools////////////////////////////

int db_read_and_div(ifstream &in1, int SEG_no,
                    int option_l, char *NR_seg, 
                    char db_swap[MAX_SEG][MAX_FILE_NAME] ) {

  char raw_seq[MAX_SEQ], raw_des[MAX_DES], raw_seq1[MAX_SEQ];
  char buffer1[MAX_LINE_SIZE];
  raw_seq[0] = raw_des[0] = buffer1[0] = 0;
  int read_in = 0;
  int NR_no1 = 0;
  int i, j, k;

  ofstream out1[255];
  for (i=0; i<SEG_no; i++) {
    out1[i].open(db_swap[i]);
    if (! out1[i] ) bomb_error("Can not open file", db_swap[i]);
  }

  while(1) {
    if ( in1.eof()) break;
    in1.getline(buffer1, MAX_LINE_SIZE-2, '\n');
    if ( buffer1[0] == '>' || buffer1[0] == ';') {
      if ( read_in ) { // write last record
         strcpy(raw_seq1, raw_seq);
         format_seq(raw_seq1);
         if ( strlen(raw_seq1) > option_l ) {
           out1[ NR_seg[NR_no1] ] << raw_des << "\n" << raw_seq;
           NR_no1++;
         }
      }
      strncpy(raw_des, buffer1, MAX_DES-2);
      raw_seq[0] = 0;
    }
    else {
      read_in = 1;
      strcat(raw_seq, buffer1); strcat(raw_seq,"\n");
    }
  } // END while(1);

  if (1) { // the last record
    strcpy(raw_seq1, raw_seq);
    format_seq(raw_seq1);
    if ( strlen(raw_seq1) > option_l ) {
      out1[ NR_seg[NR_no1] ] << raw_des << "\n" << raw_seq;
      NR_no1++;
    }
  }

  for (i=0; i<SEG_no; i++) out1[i].close();
  return 0;
} // END db_read_and_write

int db_read_and_divf(ifstream &in1, int SEG_no, char *NR_seg, 
                    char db_swap[MAX_SEG][MAX_FILE_NAME] ) {

  char raw_seq[MAX_SEQ], raw_des[MAX_DES], raw_seq1[MAX_SEQ];
  char buffer1[MAX_LINE_SIZE];
  raw_seq[0] = raw_des[0] = buffer1[0] = 0;
  int read_in = 0;
  int NR_no1 = 0;
  int i, j, k;

  ofstream out1[255];
  for (i=0; i<SEG_no; i++) {
    out1[i].open(db_swap[i]);
    if (! out1[i] ) bomb_error("Can not open file", db_swap[i]);
  }

  NR_no1 = -1;
  while(1) {
    if ( in1.eof()) break;
    in1.getline(buffer1, MAX_LINE_SIZE-2, '\n');
    if ( buffer1[0] == '>' || buffer1[0] == ';') {
      NR_no1++;
    }
    out1[ NR_seg[NR_no1] ] << buffer1 << "\n";
  }

  for (i=0; i<SEG_no; i++) out1[i].close();
  return 0;
} // END db_read_and_writef

int sort_seqs_divide_segs_div(int NR_no, int *NR_len, int *NR_idx, char *NR_seg,
                           int SEG_no, char db_swap[MAX_SEG][MAX_FILE_NAME],
                           char db_out[]) {
  int i, j, k, i1;

  int len, len1, len2, len22;
  long long total_letter=0, letter_per_seg;
  int max_len = 0, min_len = 9999999;
  for (i=0; i<NR_no; i++) {
    len = NR_len[i];
    total_letter += len;
    if (len > max_len) max_len = len;
    if (len < min_len) min_len = len;
  }
  letter_per_seg = (int) (total_letter/SEG_no);
  cout << "longest and shortest : " << max_len << " and " << min_len << endl;
  cout << "Total letters: " << total_letter << endl;

  // **************************** Form NR_idx[], Sort them from Long to short
  int *size_no;
  int *size_begin;
  if ((size_no = new int[max_len-min_len+1]) == NULL ) bomb_error("Memory");
  if ((size_begin = new int[max_len-min_len+1]) == NULL ) bomb_error("Memory");

  for (i=max_len; i>=min_len; i--) {
    size_no[max_len - i] = 0;
    size_begin[max_len - i] = 0;
  }
  for (i=0; i<NR_no; i++)  size_no[max_len - NR_len[i]]++;
  for (i=max_len; i>=min_len; i--)
    for (j=max_len; j>i; j--)
      size_begin[max_len-i] += size_no[max_len-j];
  for (i=max_len; i>=min_len; i--) size_no[max_len - i] = 0;
  for (i=0; i<NR_no; i++) {
    j = max_len-NR_len[i];
    NR_idx[ size_begin[j] + size_no[j]] = i;
    size_no[j]++;
  }
  delete []size_no; delete []size_begin;
  cout << "Sequences have been sorted" << endl;
  // END sort them from long to short

  long long letter_this_seg;

  letter_this_seg = 0;
  k=0; // current seg
  db_swap[k][0] = 0;
  sprintf(db_swap[k], "%s-%d",db_out,k);
  cout << "db " << k << " " << db_swap[k];
  max_len = 0;
  min_len = 9999999;
  for (i1=0; i1<NR_no; i1++) {
    i = NR_idx[i1];
    len = NR_len[i];
    if (len>max_len) max_len = len;
    if (len<min_len) min_len = len;
    letter_this_seg += len;

    NR_seg[i] = k;
    if ( letter_this_seg > letter_per_seg) {
      cout << " length " << max_len << " to " << min_len << endl;
      k++;
      if (k == SEG_no) k--; // k can not be over SEG_no
      db_swap[k][0] = 0;
      sprintf(db_swap[k], "%s-%d",db_out,k);
      cout << "db " << k << " " << db_swap[k];
      letter_this_seg = 0;
      max_len = 0;
      min_len = 9999999;
    }
  }
  cout << " length " << max_len << " to " << min_len << endl;
  return 0;
}// END sort_seqs_divide_segs


/////////////////////////// END ALL ////////////////////////
