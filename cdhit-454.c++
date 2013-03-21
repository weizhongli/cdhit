// =============================================================================
// CD-HIT
// http://cd-hit.org/
// http://bioinformatics.burnham-inst.org/cd-hi
//
// program written by
//                    Weizhong Li
//                    UCSD, San Diego Supercomputer Center
//                    La Jolla, CA, 92093
//                    Email liwz@sdsc.edu
//                 at
//                    Adam Godzik's lab
//                    The Burnham Institute
//                    La Jolla, CA, 92037
//                    Email adam@burnham-inst.org
//
// Modified by:
//                    Limin Fu
//                    Center for Research in Biological Systems (CRBS), UCSD
//                    La Jolla, CA, 92093
//                    Email: l2fu@ucsd.edu, fu@daovm.net
// =============================================================================

#include "cdhit-common.h"
#include "cdhit-utility.h"

#undef MAX_UAA
#define MAX_UAA 4

Options options;
SequenceDB seq_db;

// matrix below is for highly similar seqs
// with this matrix, alignment can be maintained if there is a mismatch
// at very end of alignment
// |x||||||||
//  ^
//  mismatch will cost -2 with BLOSUM62_na, therefore the first 2 bases won't be in alignment
//  but with BLOSUM62_na2, the first 2 bases are in alignment
int myBLOSUM62_na2[] = {
  2,               // A
 -1, 2,            // C
 -1,-1, 2,         // G
 -1,-1,-1, 2,      // T
 -1,-1,-1, 2, 2,   // U
 -1,-1,-1,-1,-1, 2 // N
//A  C  G  T  U  N
//0  1  2  3  3  4
};

void setaa_to_na();

////////////////////////////////////  MAIN /////////////////////////////////////
int main(int argc, char *argv[])
{
	string db_in;
	string db_out;

	options.isEST = 1;
	options.is454 = 1;
	options.NAA = 10;
	options.NAAN = NAA8;
	options.NAA_top_limit = 12;
	seq_db.NAAN = NAA8;

	options.cluster_thd = 0.98;
	options.band_width = 10;
	options.print = 1;
	options.des_len = 0;
	options.option_r = 0;

	setaa_to_na();
	mat.set_gap(-3,-1); //instead of -6 -1 to maintain maxium length of alignment
	mat.set_matrix(myBLOSUM62_na2);

	float begin_time = current_time();
	float end_time;

	// ***********************************    parse command line and open file
	if (argc < 5) print_usage_454(argv[0]);
	if (options.SetOptions( argc, argv, false, true ) == 0) print_usage_454(argv[0]);
	options.Validate();

	db_in = options.input;
	db_out = options.output;

	InitNAA( MAX_UAA );
	options.NAAN = NAAN_array[options.NAA];
	seq_db.NAAN = NAAN_array[options.NAA];

	//printf( "%i  %i  %i\n", sizeof(NVector<IndexCount>), seq_db.NAAN, sizeof(NVector<IndexCount>) * seq_db.NAAN );

	seq_db.Read( db_in.c_str(), options );
	cout << "total seq: " << seq_db.sequences.size() << endl;
	seq_db.SortDivide( options );

	seq_db.DoClustering( options );

	printf( "writing new database\n" );
	seq_db.WriteClusters( db_in.c_str(), db_out.c_str(), options );

	// write a backup clstr file in case next step crashes
	seq_db.WriteExtra1D( options );
	cout << "program completed !" << endl << endl;
	end_time = current_time();
	printf( "Total CPU time %.2f\n", end_time - begin_time );
	return 0;
} // END int main
