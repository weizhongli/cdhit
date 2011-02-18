// =============================================================================
// CD-HI-EST
// http://cd-hit.org/
// Cluster Database at High Identity (EST version)
// modified from CD-HI
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
//over-write some defs in cd-hi.h for est version
#undef MAX_UAA
#define MAX_UAA 4

//over-write some defs in cd-hi-init.h for est version

void setaa_to_na();
void make_comp_iseq(int len, char *iseq_comp, char *iseq);
void make_comp_short_word_index(int NAA, int *NAAN_array, Vector<int> & Comp_AAN_idx);

Options options;
SequenceDB seq_db;
SequenceDB seq_db2;
struct tms CPU_current, CPU_begin, CPU_end;

// next two control how if seqs in db2 is longer than reps in db1
// by deault, only seqs in db2 that are shorter than rep in db1
// are clustered to the rep in db1


////////////////////////////////////  MAIN /////////////////////////////////////
int main(int argc, char **argv)
{
	string db_in;
	string db_in2;
	string db_out;

	options.cluster_thd = 0.95;
	options.NAA = 10;
	options.NAAN = NAA8;
	seq_db.NAAN = NAA8;
	options.NAA_top_limit = 12;
	setaa_to_na();
	mat.set_to_na(); //mat.set_gap(-6,-1);

	times(&CPU_begin);

	// ***********************************    parse command line and open file
	if (argc < 7) print_usage_est_2d(argv[0]);
	if (options.SetOptions( argc, argv, true, true ) == 0) print_usage_est_2d(argv[0]);
	options.Validate();

	db_in = options.input;
	db_in2 = options.input2;
	db_out = options.output;

	InitNAA( MAX_UAA );
	options.NAAN = NAAN_array[options.NAA];
	seq_db.NAAN = NAAN_array[options.NAA];
	seq_db2.NAAN = NAAN_array[options.NAA];

	if ( options.option_r ) {
		Comp_AAN_idx.resize( seq_db.NAAN );
		make_comp_short_word_index(options.NAA, NAAN_array, Comp_AAN_idx);
	}

	seq_db.Read( db_in.c_str(), options );
	cout << "total seq in db1: " << seq_db.sequences.size() << endl;

	seq_db2.Read( db_in2.c_str(), options );
	cout << "total seq in db2: " << seq_db2.sequences.size() << endl;

	seq_db.SortDivide( options );
	seq_db2.SortDivide( options, false );
	seq_db2.ClusterTo( seq_db, options );

	cout << "writing non-redundant sequences from db2" << endl;
	seq_db2.WriteClusters( db_in2.c_str(), db_out.c_str(), options );

	seq_db2.WriteExtra2D( seq_db, options );
	cout << "program completed !" << endl << endl;
	times(&CPU_end);
	show_cpu_time(CPU_begin, CPU_end);
	return 0;
} // END int main

