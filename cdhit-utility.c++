
#include<stdlib.h>
#include<iostream>
#include"cdhit-common.h"

using namespace std;

// information
char cd_hit_ver[]  = "\t\t====== CD-HIT version " CDHIT_VERSION " (built on " __DATE__ ") ======";
char cd_hit_ref1[] = "\"Clustering of highly homologous sequences to reduce thesize of large protein database\", Weizhong Li, Lukasz Jaroszewski & Adam Godzik. Bioinformatics, (2001) 17:282-283";
char cd_hit_ref2[] = "\"Tolerating some redundancy significantly speeds up clustering of large protein databases\", Weizhong Li, Lukasz Jaroszewski & Adam Godzik. Bioinformatics, (2002) 18:77-82";
char cd_hit_ref3[] = "\"Cd-hit: a fast program for clustering and comparing large sets of protein or nucleotide sequences\", Weizhong Li & Adam Godzik. Bioinformatics, (2006) 22:1658-1659";
char cd_hit_ref4[] = "\"Beifang Niu, Limin Fu, Shulei Sun and Weizhong Li. Artificial and natural duplicates in pyrosequencing reads of metagenomic data. BMC Bioinformatics (2010) 11:187";
//

char contacts[] =
  "    Questions, bugs, contact Limin Fu at l2fu@ucsd.edu, or Weizhong Li at liwz@sdsc.edu\n"
  "    For updated versions and information, please visit: http://cd-hit.org\n\n"
  "    cd-hit web server is also available from http://cd-hit.org\n\n"
  "    If you find cd-hit useful, please kindly cite:\n\n";

char txt_option_i[] = "\tinput filename in fasta format, required\n";
char txt_option_i_2d[] = "\tinput filename for db1 in fasta format, required\n";
char txt_option_i2[] = "\tinput filename for db2 in fasta format, required\n";
char txt_option_o[] = "\toutput filename, required\n";
char txt_option_c[] = 
"\tsequence identity threshold, default 0.9\n \
\tthis is the default cd-hit's \"global sequence identity\" calculated as:\n \
\tnumber of identical amino acids in alignment\n \
\tdivided by the full length of the shorter sequence\n";
char txt_option_G[] = 
"\tuse global sequence identity, default 1\n \
\tif set to 0, then use local sequence identity, calculated as :\n \
\tnumber of identical amino acids in alignment\n \
\tdivided by the length of the alignment\n \
\tNOTE!!! don't use -G 0 unless you use alignment coverage controls\n \
\tsee options -aL, -AL, -aS, -AS\n";
char txt_option_g[] =
"\t1 or 0, default 0\n \
\tby cd-hit's default algorithm, a sequence is clustered to the first \n \
\tcluster that meet the threshold (fast cluster). If set to 1, the program\n \
\twill cluster it into the most similar cluster that meet the threshold\n \
\t(accurate but slow mode)\n \
\tbut either 1 or 0 won't change the representatives of final clusters\n";
char txt_option_b[] = "\tband_width of alignment, default 20\n";
char txt_option_M[] = "\tmemory limit (in MB) for the program, default 800; 0 for unlimitted;\n";
char txt_option_n[] = "\tword_length, default 5, see user's guide for choosing it\n";
char txt_option_n_est[] = "\tword_length, default 10, see user's guide for choosing it\n";
char txt_option_l[] = "\tlength of throw_away_sequences, default 10\n";
char txt_option_t[] = "\ttolerance for redundance, default 2\n";
char txt_option_T[] = "\tnumber of threads, default 1; with 0, all CPUs will be used\n";
char txt_option_d[] =
"\tlength of description in .clstr file, default 20\n \
\tif set to 0, it takes the fasta defline and stops at first space\n";
char txt_option_s[] =
"\tlength difference cutoff, default 0.0\n \
\tif set to 0.9, the shorter sequences need to be\n \
\tat least 90% length of the representative of the cluster\n";
char txt_option_S[] =
"\tlength difference cutoff in amino acid, default 999999\n \
\tif set to 60, the length difference between the shorter sequences\n \
\tand the representative of the cluster can not be bigger than 60\n";
char txt_option_s2[] =
"\tlength difference cutoff for db1, default 1.0\n \
\tby default, seqs in db1 >= seqs in db2 in a same cluster\n \
\tif set to 0.9, seqs in db1 may just >= 90% seqs in db2\n";
char txt_option_S2[] =
"\tlength difference cutoff, default 0\n \
\tby default, seqs in db1 >= seqs in db2 in a same cluster\n \
\tif set to 60, seqs in db2 may 60aa longer than seqs in db1\n";
char txt_option_aL[] = 
"\talignment coverage for the longer sequence, default 0.0\n \
\tif set to 0.9, the alignment must covers 90% of the sequence\n";
char txt_option_AL[] = 
"\talignment coverage control for the longer sequence, default 99999999\n \
\tif set to 60, and the length of the sequence is 400,\n \
\tthen the alignment must be >= 340 (400-60) residues\n";
char txt_option_aS[] = 
"\talignment coverage for the shorter sequence, default 0.0\n \
\tif set to 0.9, the alignment must covers 90% of the sequence\n";
char txt_option_AS[] = 
"\talignment coverage control for the shorter sequence, default 99999999\n \
\tif set to 60, and the length of the sequence is 400,\n \
\tthen the alignment must be >= 340 (400-60) residues\n";
char txt_option_A[] = 
"\tminimal alignment coverage control for the both sequences, default 0\n \
\talignment must cover >= this value for both sequences \n"; 
char txt_option_B[] =
"\t1 or 0, default 0, by default, sequences are stored in RAM\n \
\tif set to 1, sequence are stored on hard drive\n \
\tit is recommended to use -B 1 for huge databases\n";
char txt_option_uL[] = 
"\tmaximum unmatched percentage for the longer sequence, default 1.0\n \
\tif set to 0.1, the unmatched region (excluding leading and tailing gaps)\n \
\tmust not be more than 10% of the sequence\n";
char txt_option_uS[] = 
"\tmaximum unmatched percentage for the shorter sequence, default 1.0\n \
\tif set to 0.1, the unmatched region (excluding leading and tailing gaps)\n \
\tmust not be more than 10% of the sequence\n";
char txt_option_U[] = 
"\tmaximum unmatched length, default 99999999\n \
\tif set to 10, the unmatched region (excluding leading and tailing gaps)\n \
\tmust not be more than 10 bases\n";
char txt_option_p[] =
"\t1 or 0, default 0\n \tif set to 1, print alignment overlap in .clstr file\n";
char txt_option_r[] =
"\t1 or 0, default 1, by default do both +/+ & +/- alignments\n \
\tif set to 0, only +/+ strand alignment\n";

char txt_option_match[] = "\tmatching score, default 2 (1 for T-U and N-N)\n";
char txt_option_match2[] = "\tmatching score, default 2\n";
char txt_option_mismatch[] = "\tmismatching score, default -2\n";
char txt_option_mismatch2[] = "\tmismatching score, default -1\n";
char txt_option_gap[] = "\tgap opening score, default -6\n";
char txt_option_gap2[] = "\tgap opening score, default -3\n";
char txt_option_gap_ext[] = "\tgap extension score, default -1\n";

int print_usage (char *arg) {
  cout << cd_hit_ver << "\n\n" ;
  cout << "Usage: "<< arg << " [Options] \n\nOptions\n\n";
  cout << "    -i" << txt_option_i;
  cout << "    -o" << txt_option_o;
  cout << "    -c" << txt_option_c;
  cout << "    -G" << txt_option_G;
  cout << "    -b" << txt_option_b;
  cout << "    -M" << txt_option_M;
  cout << "    -T" << txt_option_T;
  cout << "    -n" << txt_option_n;
  cout << "    -l" << txt_option_l;
  cout << "    -t" << txt_option_t;
  cout << "    -d" << txt_option_d;
  cout << "    -s" << txt_option_s;
  cout << "    -S" << txt_option_S;
  cout << "    -aL" << txt_option_aL;
  cout << "    -AL" << txt_option_AL;
  cout << "    -aS" << txt_option_aS;
  cout << "    -AS" << txt_option_AS;
  cout << "    -A" << txt_option_A;
  cout << "    -uL" << txt_option_uL;
  cout << "    -uS" << txt_option_uS;
  cout << "    -U" << txt_option_U;
  cout << "    -B" << txt_option_B;
  cout << "    -p" << txt_option_p;
  cout << "    -g" << txt_option_g;
  cout << "    -h print this help\n\n";
  cout << contacts;
  cout << "    " << cd_hit_ref1 << "\n";
  cout << "    " << cd_hit_ref2 << "\n\n\n";
  exit(1);
} // END print_usage



int print_usage_2d (char *arg) {
  cout << cd_hit_ver << "\n\n" ;
  cout << "Usage: "<< arg << " [Options] \n\nOptions\n\n";
  cout << "    -i" << txt_option_i_2d;
  cout << "    -i2"<< txt_option_i2;
  cout << "    -o" << txt_option_o;
  cout << "    -c" << txt_option_c;
  cout << "    -G" << txt_option_G;
  cout << "    -b" << txt_option_b;
  cout << "    -M" << txt_option_M;
  cout << "    -T" << txt_option_T;
  cout << "    -n" << txt_option_n;
  cout << "    -l" << txt_option_l;
  cout << "    -t" << txt_option_t;
  cout << "    -d" << txt_option_d;
  cout << "    -s" << txt_option_s;
  cout << "    -S" << txt_option_S;
  cout << "    -s2" << txt_option_s2;
  cout << "    -S2" << txt_option_S2;
  cout << "    -aL" << txt_option_aL;
  cout << "    -AL" << txt_option_AL;
  cout << "    -aS" << txt_option_aS;
  cout << "    -AS" << txt_option_AS;
  cout << "    -A" << txt_option_A;
  cout << "    -uL" << txt_option_uL;
  cout << "    -uS" << txt_option_uS;
  cout << "    -U" << txt_option_U;
  cout << "    -B" << txt_option_B;
  cout << "    -p" << txt_option_p;
  cout << "    -g" << txt_option_g;
  cout << "    -h print this help\n\n";
  cout << "    Questions, bugs, contact Weizhong Li at liwz@sdsc.edu\n\n";
  cout << "    If you find cd-hit useful, please kindly cite:\n\n";
  cout << "    " << cd_hit_ref1 << "\n";
  cout << "    " << cd_hit_ref3 << "\n\n\n";
  exit(1);
} // END print_usage_2d


int print_usage_est (char *arg) {
  cout << cd_hit_ver << "\n\n" ;
  cout << "Usage: "<< arg << " [Options] \n\nOptions\n\n";
  cout << "    -i" << txt_option_i;
  cout << "    -o" << txt_option_o;
  cout << "    -c" << txt_option_c;
  cout << "    -G" << txt_option_G;
  cout << "    -b" << txt_option_b;
  cout << "    -M" << txt_option_M;
  cout << "    -T" << txt_option_T;
  cout << "    -n" << txt_option_n_est;
  cout << "    -l" << txt_option_l;
  cout << "    -d" << txt_option_d;
  cout << "    -s" << txt_option_s;
  cout << "    -S" << txt_option_S;
  cout << "    -aL" << txt_option_aL;
  cout << "    -AL" << txt_option_AL;
  cout << "    -aS" << txt_option_aS;
  cout << "    -AS" << txt_option_AS;
  cout << "    -A" << txt_option_A;
  cout << "    -uL" << txt_option_uL;
  cout << "    -uS" << txt_option_uS;
  cout << "    -U" << txt_option_U;
  cout << "    -B" << txt_option_B;
  cout << "    -p" << txt_option_p;
  cout << "    -g" << txt_option_g;
  cout << "    -r" << txt_option_r;
  cout << "    -match" << txt_option_match;
  cout << "    -mismatch" << txt_option_mismatch;
  cout << "    -gap" << txt_option_gap;
  cout << "    -gap-ext" << txt_option_gap_ext;
  cout << "    -h print this help\n\n";
  cout << contacts;
  cout << "    " << cd_hit_ref1 << "\n";
  cout << "    " << cd_hit_ref3 << "\n\n\n";
  exit(1);
} // END print_usage_est


int print_usage_est_2d (char *arg) {
  cout << cd_hit_ver << "\n\n" ;
  cout << "Usage: "<< arg << " [Options] \n\nOptions\n\n";
  cout << "    -i" << txt_option_i_2d;
  cout << "    -i2"<< txt_option_i2;
  cout << "    -o" << txt_option_o;
  cout << "    -c" << txt_option_c;
  cout << "    -G" << txt_option_G;
  cout << "    -b" << txt_option_b;
  cout << "    -M" << txt_option_M;
  cout << "    -T" << txt_option_T;
  cout << "    -n" << txt_option_n_est;
  cout << "    -l" << txt_option_l;
  cout << "    -d" << txt_option_d;
  cout << "    -s" << txt_option_s;
  cout << "    -S" << txt_option_S;
  cout << "    -s2" << txt_option_s2;
  cout << "    -S2" << txt_option_S2;
  cout << "    -aL" << txt_option_aL;
  cout << "    -AL" << txt_option_AL;
  cout << "    -aS" << txt_option_aS;
  cout << "    -AS" << txt_option_AS;
  cout << "    -A" << txt_option_A;
  cout << "    -uL" << txt_option_uL;
  cout << "    -uS" << txt_option_uS;
  cout << "    -U" << txt_option_U;
  cout << "    -B" << txt_option_B;
  cout << "    -p" << txt_option_p;
  cout << "    -g" << txt_option_g;
  cout << "    -r" << txt_option_r;
  cout << "    -match" << txt_option_match;
  cout << "    -mismatch" << txt_option_mismatch;
  cout << "    -gap" << txt_option_gap;
  cout << "    -gap-ext" << txt_option_gap_ext;
  cout << "    -h print this help\n\n";
  cout << contacts;
  cout << "    " << cd_hit_ref1 << "\n";
  cout << "    " << cd_hit_ref3 << "\n\n\n";
  exit(1);
} // END print_usage_est_2d


int print_usage_div (char *arg) {
  cout << cd_hit_ver << "\n\n" ;
  cout << "Usage: "<< arg << " [Options] \n\nOptions\n\n";
  cout << "Options " << endl << endl;
  cout << "    -i in_dbname, required" << endl;
  cout << "    -o out_dbname, required" << endl;
  cout << "    -div number of divide, required " << endl;
//  cout << "    -dbmax max size of your db\n\n\n";
  exit(1);
} // END print_usage_div

char mytxt_option_c[] =
"\tsequence identity threshold, default 0.98\n \
\tthis is a \"global sequence identity\" calculated as :\n \
\tnumber of identical amino acids in alignment\n \
\tdivided by the full length of the shorter sequence + gaps\n";
char mytxt_option_b[] = "\tband_width of alignment, default 10\n";
char mytxt_option_n_est[] = "\tword_length, default 10, see user's guide for choosing it\n";
char mytxt_option_D[] = "\tmax size per indel, default 1\n";
int print_usage_454 (char *arg)
{
  cout << cd_hit_ver << "\n\n" ;
  cout << "Usage: "<< arg << " [Options] \n\nOptions\n\n";
  cout << "    -i" << txt_option_i;
  cout << "    -o" << txt_option_o;
  cout << "    -c" << mytxt_option_c;
  cout << "    -b" << mytxt_option_b;
  cout << "    -M" << txt_option_M;
  cout << "    -T" << txt_option_T;
  cout << "    -n" << mytxt_option_n_est;
  cout << "    -aL" << txt_option_aL;
  cout << "    -AL" << txt_option_AL;
  cout << "    -aS" << txt_option_aS;
  cout << "    -AS" << txt_option_AS;
  cout << "    -B" << txt_option_B;
  cout << "    -g" << txt_option_g;
  cout << "    -D" << mytxt_option_D;
  cout << "    -match" << txt_option_match2;
  cout << "    -mismatch" << txt_option_mismatch2;
  cout << "    -gap" << txt_option_gap2;
  cout << "    -gap-ext" << txt_option_gap_ext;
  cout << "    -h print this help\n\n";
  cout << "    Questions, bugs, contact Weizhong Li at liwz@sdsc.edu\n\n";
  cout << "    If you find cd-hit useful, please kindly cite:\n\n";
  cout << "    " << cd_hit_ref1 << "\n";
  cout << "    " << cd_hit_ref3 << "\n";
  cout << "    " << cd_hit_ref4 << "\n\n\n";
  exit(1);
}
