// =============================================================================
// CD-HIT
// http://cd-hit.org/
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
//
// Modified by:
//                                      Limin Fu
//                                      Center for Research in Biological Systems (CRBS), UCSD
//                                      La Jolla, CA, 92093
//                                      Email: l2fu@ucsd.edu, fu@daovm.net
// =============================================================================

#include "cdhit-common.h"

Options options;
SequenceDB seq_db;

////////////////////////////////////  MAIN /////////////////////////////////////
int main(int argc, char *argv[])
{
	string db_in;
	string db_out;
	int i, div = 1;

	float begin_time = current_time();
	float end_time;

	// ***********************************    parse command line and open file
	if (argc < 5) print_usage_div(argv[0]);
	for (i=1; i<argc; i++) {
		if      (strcmp(argv[i], "-i"    )==0) db_in =  argv[++i];
		else if (strcmp(argv[i], "-o"    )==0) db_out = argv[++i];
		else if (strcmp(argv[i], "-div"  )==0) div    = atoi(argv[++i]);
		else                                   print_usage_div(argv[0]);
	}
	if( div <= 1 ){
		printf( "Warning: -div must be greater than 1.\n" );
		printf( "Warning: no database is writen.\n" );
		return 0;
	}
	options.store_disk = 1;

	//printf( "%i  %i  %i\n", sizeof(NVector<IndexCount>), seq_db.NAAN, sizeof(NVector<IndexCount>) * seq_db.NAAN );

	seq_db.Read( db_in.c_str(), options );
	cout << "total seq: " << seq_db.sequences.size() << endl;
	seq_db.SortDivide( options );

	printf( "writing new databases\n" );
	seq_db.DivideSave( db_in.c_str(), db_out.c_str(), div, options );

	end_time = current_time();
	printf( "Total CPU time %.2f\n", end_time - begin_time );

	cout << "program completed !" << endl << endl;
	return 0;
} // END int main

