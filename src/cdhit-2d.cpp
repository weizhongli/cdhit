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

// next two control how if seqs in db2 is longer than reps in db1
// by deault, only seqs in db2 that are shorter than rep in db1 
// are clustered to the rep in db1

Options options;
SequenceDB seq_db;
SequenceDB seq_db2;

////////////////////////////////////  MAIN /////////////////////////////////////
int main(int argc, char **argv)
{
	string db_in;
	string db_in2;
	string db_out;

	float begin_time = current_time();
	float end_time;

	// ***********************************    parse command line and open file
	if (argc < 7) print_usage_2d(argv[0]);
	if (options.SetOptions( argc, argv, true ) == 0) print_usage_2d(argv[0]);
	options.Validate();

	db_in = options.input;
	db_in2 = options.input2;
	db_out = options.output;

	InitNAA( MAX_UAA );
	options.NAAN = NAAN_array[options.NAA];
	seq_db.NAAN = seq_db2.NAAN = NAAN_array[options.NAA];

	seq_db.Read( db_in.c_str(), options );
	cout << "total seq in db1: " << seq_db.sequences.size() << endl;

	seq_db2.Read( db_in2.c_str(), options );
	cout << "total seq in db2: " << seq_db2.sequences.size() << endl;

	seq_db.SortDivide( options );
	seq_db2.SortDivide( options, false );
	seq_db2.ClusterTo( seq_db, options );

	cout << "writing non-redundant sequences from db2" << endl;
	seq_db2.WriteClusters( db_in2.c_str(), db_out.c_str(), options );

	// write a backup clstr file in case next step crashes
	seq_db2.WriteExtra2D( seq_db, options );
	cout << "program completed !" << endl << endl;
	end_time = current_time();
	printf( "Total CPU time %.2f\n", end_time - begin_time );
	return 0;
} // END int main

