#!/usr/bin/perl -w

our $script_name = $0;
our $script_dir = $0;
   $script_dir =~ s/[^\/]+$//;
   $script_dir = "./" unless ($script_dir);
require "$script_dir/psi-cd-hit-local.pl";

parse_para_etc(@ARGV);
open_LOG();

our @seqs    = ();
our @dess    = ();
our @idens   = ();
our @lens    = ();
our @covs    = ();
our @passeds = ();
our @NR_clstr_nos = ();
our @in_bg   = ();
our @NR_idx = ();
our $NR_no  = 0;
my ($i, $j, $k, $i0, $j0, $k0, $ll);

read_db();

our $NR_passed = 0;
our $formatdb_no = $NR_no;;

@NR_idx = (0..($NR_no-1));
@NR_idx = sort { $lens[$b] <=> $lens[$a] } @NR_idx unless (-e $restart_in);
start_deamon_master();

our $NR90_no = 0;
our @NR90_seq = ();

$i0 = 0;
if ( -e $restart_in) { $i0 = read_restart(); } ##  restart after crash

blast_formatdb();
for (; $i0<$NR_no; $i0++) {
  $i = $NR_idx[$i0];
  my $fout = ($host_no>0) ? "$bl_dir/$i.out" : "$bl_dir/$i";
  if ( ($host_no>0)     and (not (-e $fout)) and (not $passeds[$i] )) {
    run_batch_blast3($i0) unless ($in_bg[$i]);
  }

  if ( not $passeds[$i] ) { # this is a new representative
    $NR_passed++;
    $NR_clstr_nos[$i]   = $NR90_no;
    $idens[$i]          = "*";
    $covs[$i]           = "100%";
    $passeds[$i]        = 1;
    $NR90_seq[$NR90_no] = [$i];
    fish_other_homolog($i);
    $NR90_no++;
  }

  $seqs[$i] = "" unless ($print_db);

  watch_progress($i0, $NR90_no, $NR_passed, $NR_no, 0);

  if (($i0+1) % $restart_seg == 0 ) { 
    write_restart(); write_db_clstr(); remove_raw_blout($i0);
  }
  if ($formatdb_no - ($NR_no-$NR_passed) >= $reformat_seg) {blast_formatdb(); }
}
## END for ($i=0; $i<$NR_no; $i++)
watch_progress($NR_no-1, $NR90_no, $NR_passed, $NR_no, 1);

if ($print_db) {
  open(DBOUT, "> $db_out")   || die "Can not write $db_out";
  for ($i=0; $i<$NR_no; $i++) {
    next unless ($idens[$i] eq "*");
    my $seq = $seqs[$i];
    $seq =~ s/(.{70})/$1\n/g;
    $seq =~ s/\n$//;
    print DBOUT "$dess[$i]\n$seq\n";
  }
  close(DBOUT);
}

write_restart();
write_db_clstr();
remove_blast_db();
close_LOG();


