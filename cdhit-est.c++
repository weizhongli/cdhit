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
void make_comp_short_word_index(int NAA, int *NAAN_array, Vector<int> & Comp_AAN_idx);
void make_comp_iseq(int len, char *iseq_comp, char *iseq);


Options options;
SequenceDB seq_db;

////////////////////////////////////  MAIN /////////////////////////////////////
int main(int argc, char **argv) 
{
	string db_in;
	string db_out;
	string db_in_pe;
	string db_out_pe;

	options.cluster_thd = 0.95;
	options.NAA = 10;
	options.NAAN = NAA8;
	seq_db.NAAN = NAA8;
	options.NAA_top_limit = 12;
	setaa_to_na();
	mat.set_to_na(); //mat.set_gap(-6,-1);

	float begin_time = current_time();
	float end_time;

	// ***********************************    parse command line and open file
	if (argc < 5) print_usage_est(argv[0]);
	if (options.SetOptions( argc, argv, false, true ) == 0) print_usage_est(argv[0]);
	options.Validate();

	db_in     = options.input;
	db_in_pe  = options.input_pe;
	db_out    = options.output;
	db_out_pe = options.output_pe;

	InitNAA( MAX_UAA );
	seq_db.NAAN = NAAN_array[options.NAA];

	if ( options.option_r ) {
		Comp_AAN_idx.resize( seq_db.NAAN );
		make_comp_short_word_index(options.NAA, NAAN_array, Comp_AAN_idx);
	}

        if ( options.PE_mode ) {seq_db.Read( db_in.c_str(), db_in_pe.c_str(), options );}
        else                   {seq_db.Read( db_in.c_str(),                   options );} 

	cout << "total seq: " << seq_db.sequences.size() << endl;
	seq_db.SortDivide( options );
	seq_db.DoClustering( options );

	printf( "writing new database\n" );
        if ( options.PE_mode ) { seq_db.WriteClusters( db_in.c_str(), db_in_pe.c_str(), db_out.c_str(), db_out_pe.c_str(), options ); }
        else                   { seq_db.WriteClusters( db_in.c_str(),                   db_out.c_str(),                    options ); }

	// write a backup clstr file in case next step crashes
	seq_db.WriteExtra1D( options );
	cout << "program completed !" << endl << endl;
	end_time = current_time();
	printf( "Total CPU time %.2f\n", end_time - begin_time );
	return 0;
}
