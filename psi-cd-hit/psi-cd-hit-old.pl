#!/usr/bin/perl -w
################################################################################
######### PSI-cd-hit written by Weizhong Li at http://cd-hit.org
################################################################################

our $script_name = $0;
our $script_dir = $0;
   $script_dir =~ s/[^\/]+$//;
   $script_dir = "./" unless ($script_dir);
require "$script_dir/psi-cd-hit-local-old.pl";

parse_para_etc(@ARGV);
open_LOG();

our @seqs    = ();
our @dess    = ();
our @idens   = ();
our @lens    = ();
our @passeds = ();
our @NR_clstr_nos = ();
our @in_bg   = ();
our @NR_idx = ();
our $NR_no  = 0;
our $DB_no  = 0;
our $DB_len = 0;
our $DB_len0 = 0;
our $DB_len_reduced = 0;
our $DB_len_reduced2 = 0;  #### for write_restart etc purpose

our $opt_aL_upper_band = 0; #### sequences < this length will not be submitted unless reformatdb
our $opt_al_upper_bandi= 0;
our $opt_aL_lower_band = 0; #### sequences < this length don't need to be searched
my ($i, $j, $k, $i0, $j0, $k0, $ll);

read_db();

our $NR_passed = 0;
our $formatdb_no = $NR_no;;

@NR_idx = (0..($NR_no-1));
@NR_idx = sort { $lens[$b] <=> $lens[$a] } @NR_idx unless (-e $restart_in);

our $NR90_no = 0;
our @NR90_seq = ();

$i0 = 0;
if ( -e $restart_in) { $i0 = read_restart(); } ##  restart after crash
elsif ($skip_long > 0) { #### skip very long seqs
  for (; $i0<$NR_no; $i0++) {
    $i = $NR_idx[$i0];
    last if ($lens[$i] < $skip_long);
    
    $NR_passed++;
    $NR_clstr_nos[$i]   = $NR90_no;
    $idens[$i]          = "*";
    $passeds[$i]        = 1;
    $NR90_seq[$NR90_no] = [$i];
    $NR90_no++;
    $DB_len_reduced += $lens[$i];
  }
}

#### set init opt_aL_uppper/lower_bands
if ( ($opt_aL > 0.3) ) {
  die ("option -aL > 1.0") if ($opt_aL > 1.0);

####################
###################
##################
#################
################
###############  <-upper band 
##############  <- seq below not submit, unless band change
#############
############
###########
##########  <- lower band
######### <- seq below not in format db
########
#######
#####
####
###
##
#
  my $total_jobs = $batch_no_per_node * $num_qsub * $para_no;
  my $space = ($total_jobs > $restart_seg) ? $total_jobs : $restart_seg;
  my $d1 = $i0+$space;
     $d1 = ($NR_no-1) if ($d1 >= $NR_no-1); 
  $opt_aL_upper_band = $lens[ $NR_idx[$d1] ];
  $opt_aL_lower_band = int($opt_aL_upper_band * $opt_aL);
  $opt_aL_upper_bandi= $d1;
  write_LOG("set opt_aL_band $opt_aL_upper_band($opt_aL_upper_bandi) $opt_aL_lower_band");
}


($DB_no, $DB_len) = blast_formatdb();
$DB_len0 = $DB_len;
$DB_len_reduced = 0;
$DB_len_reduced2 = 0;
for (; $i0<$NR_no; $i0++) {
  $i = $NR_idx[$i0];
  run_batch_blast3($i0) unless ($in_bg[$i] or (-e "$bl_dir/$i.out") or $passeds[$i]);

  if ( not $passeds[$i] ) { # this is a new representative
    $NR_passed++;
    $NR_clstr_nos[$i]   = $NR90_no;
    $idens[$i]          = "*";
    $passeds[$i]        = 1;
    $NR90_seq[$NR90_no] = [$i];
    fish_other_homolog($i);
    $NR90_no++;
    $DB_len_reduced += $lens[$i];
    $DB_len_reduced2 += $lens[$i];
  }

  watch_progress($i0, $NR90_no, $NR_passed, $NR_no, 0);

  if ((($i0+1) % $restart_seg == 0) or ($DB_len_reduced2 > $DB_len0/10) ) { 
    write_restart(); write_db_clstr(); remove_raw_blout_bg($i0);
    $DB_len_reduced2 = 0;
  }

  my $opt_aL_format_flag = 0;
  if ( ($opt_aL > 0.3) ) { #### formatdb maybe needed if current length of seq.i0 close to opt_aL_upper_band
    my $total_jobs = $batch_no_per_node * $num_qsub * $para_no;
    if ( ($opt_aL_upper_bandi - $i0) < $total_jobs ) { #### seqs left for possible submission < total_jobs

      my $space = ($total_jobs > $restart_seg) ? $total_jobs : $restart_seg;
      my $d1 = $i0+$space;
         $d1 = ($NR_no-1) if ($d1 >= $NR_no-1); 
      $opt_aL_upper_band = $lens[ $NR_idx[$d1] ];
      $opt_aL_lower_band = int($opt_aL_upper_band * $opt_aL);
      $opt_aL_upper_bandi= $d1;
      $opt_aL_format_flag = 1;
      write_LOG("set opt_aL_band $opt_aL_upper_band($opt_aL_upper_bandi) $opt_aL_lower_band");
    }
  }
  if ((($i0+1) % (int($NR_no/10)) == 0) or ($DB_len_reduced > $DB_len/10) or $opt_aL_format_flag ) {
    ($DB_no, $DB_len) = blast_formatdb();
    $DB_len_reduced = 0;
  }
  #if ($formatdb_no - ($NR_no-$NR_passed) >= $reformat_seg) {blast_formatdb(); }
}
## END for ($i=0; $i<$NR_no; $i++)
watch_progress($NR_no-1, $NR90_no, $NR_passed, $NR_no, 1);

if (1) { ### print NR db
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


