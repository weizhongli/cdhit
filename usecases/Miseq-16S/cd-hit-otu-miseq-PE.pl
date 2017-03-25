#!/usr/bin/perl

use Getopt::Std;
my $script_name = $0;
my $script_dir = $0;
   $script_dir =~ s/[^\/]+$//;
   chop($script_dir);
   $script_dir = "./" unless ($script_dir);

getopts("i:j:o:p:c:s:t:m:e:Z:a:f:d:R:",\%opts);
die usage() unless ($opts{i} and $opts{o});

my $input      = $opts{i};
my $input2     = $opts{j};
my $dir        = $opts{o};
my $abs_cutoff = $opts{a}; $abs_cutoff = 0.00005 unless ($abs_cutoff); #5e-5
my $otu_cutoff = $opts{c}; $otu_cutoff = 0.97    unless ($otu_cutoff);
my $chimera_f  = $opts{m}; $chimera_f  = "true"  unless ($chimera_f);
my $debug_mode = $opts{Z};
my $fast_mode  = $opts{f}; #### use cd-hit-dup for stage 1 and 2 clustering
my $cdhit_opt  = $opts{d}; 
my $restart_n  = $opts{R}; $restart_n  = 0       unless (defined($restart_n));
my $LOGf = "$dir/OTU.log";
my $cd_hit_dup    = "$script_dir/../../cd-hit-auxtools/cd-hit-dup"; die "no $cd_hit_dup" unless (-e $cd_hit_dup);
my $cd_hit_est    = "$script_dir/../../cd-hit-est";                 die "no $cd_hit_est" unless (-e $cd_hit_est);

my ($i, $j, $k, $str, $cmd, $ll);
$cmd = `mkdir -p $dir`;
open(LOG, "> $LOGf") || die "can not write to $LOGf";
my $f2 = "$dir/seq";

################################################################################
#### Stage 0 ----------- clustering at 100% - stage 0
################################################################################
my $clstr = "$f2.dup.clstr";
my $clstr2 = "$f2.dup2.clstr";
if ($restart_n <= 0) {
  nice_run("$cd_hit_dup -i $input -i2 $input2 -o $f2.dup -o2 $f2.dup.2 -u 100 -d 0 -m false -f $chimera_f > $f2.dup2.log");
  nice_run("cat $f2.dup.clstr $f2.dup2.clstr > $f2-stage0.clstr.tmp");
  nice_run("$script_dir/cd-hit/clstr_sort_by.pl < $f2-stage0.clstr.tmp > $f2-stage0.clstr; rm -f $f2-stage0.clstr.tmp");
  nice_run("$script_dir/clstr_sort_rep.pl $f2-stage0.clstr $input > $f2-stage0-rep.fa");
#
# /home/oasis/data/etc/git/cdhit/cd-hit-auxtools/cd-hit-dup -i qc/R1.fa -i2 qc/R2.fa -o otu/seq.dup -o2 otu/seq.dup.2 -u 100 -d 0 -f true > otu/seq.dup.log # no work 
# /home/oasis/data/etc/git/cdhit/cd-hit-auxtools/cd-hit-dup -i qc/R1.fa -i2 qc/R2.fa -o otu/seq.dup -o2 otu/seq.dup.2 -u 100 -d 0 > otu/seq.dup.log 
#
# what if cd-hit-est
# /home/oasis/data/etc/git/cdhit/cd-hit-est -i qc/R1.fa   -j qc/R2.fa     -o otu/seq.nr -op otu/seq.nr.2 -sf 1 -sc 1 -P 1 -r 0 -cx 100 -cy 100 -c 1.0  -n 10 -G 1 -b 1  -T 1 -M 8000  -d 0 -p 1 > otu/seq.nr.log
# /home/oasis/data/etc/git/cdhit/cd-hit-est -i otu/seq.nr -o otu/seq.nr.R1                                                -r 0 -cx 100         -c 1.0  -n 10 -G 1 -b 1  -T 1 -M 8000  -d 0 -p 1 > otu/seq.nr.R1.log
# /home/oasis/data/etc/git/cdhit/cd-hit-est -i otu/seq.nr.2 -o otu/seq.nr.R2                                              -r 0 -cx 100         -c 1.0  -n 10 -G 1 -b 1  -T 1 -M 8000  -d 0 -p 1 > otu/seq.nr.R2.log

# /home/oasis/data/etc/git/cdhit/cd-hit-est -i otu/seq.nr -j otu/seq.nr.2 -o otu/seq.99 -op otu/seq.99.2             -P 1 -r 0 -cx 100 -cy 100 -c 0.99 -n 10 -G 1 -b 1  -T 1 -M 8000  -d 0 -p 1 > otu/seq.99.log
# /home/oasis/data/etc/git/cdhit/cd-hit-est -i otu/seq.99 -j otu/seq.99.2 -o otu/seq.97 -op otu/seq.97.2             -P 1 -r 0 -cx 100 -cy 100 -c 0.97 -n 10 -G 1 -b 5  -T 1 -M 8000  -d 0 -p 1 > otu/seq.97.log
# do not sort 99.clstr, always trust cd-hit-dup ordered sequences
# /home/oasis/data/etc/git/cdhit/clstr_rev.pl otu/seq.nr.clstr otu/seq.99.clstr | /home/oasis/data/etc/git/cdhit/clstr_sort_by.pl > otu/seq.99-full.clstr
# /home/oasis/data/etc/git/cdhit/clstr_rev.pl otu/seq.99-full.clstr otu/seq.97.clstr | /home/oasis/data/etc/git/cdhit/clstr_sort_by.pl > otu/seq.97-full.clstr
#
# combine ref
# /home/oasis/data/etc/git/cdhit/cd-hit-est -i seq.99.wref.R1 -o seq.97.wref.R1only -r 0 -cx 100 -c 0.97 -n 10 -b 5 -T 1 -M 8000 -d 1 -p 1 -G 0 -A 50 -g 1 
#
}
if (not $debug_mode) {
  my $no1       = count_seqs_from_fasta_file($input);
  my $no_clstr  = count_clstrs_from_clstr_file($clstr);
  my $no_clstr2 = count_clstrs_from_clstr_file($clstr2);
  print LOG "Number_contigs\t$no1\n";
  print LOG "Number_unique_contigs\t$no_clstr\n";
  print LOG "Number_unique_chimaric_contigs\t$no_clstr2\n";
}

################################################################################
#### Stage 1 ---------- clustering at 99.25%  #### distance 0.75%
################################################################################
my $seq_n = `grep -c "^>" $input`; $seq_n =~ s/\D//g;
my $cutoff = int($seq_n * $abs_cutoff);
my $c1 = 0.9925; 
if ($restart_n <= 1) {
  if ($fast_mode) {
    nice_run("$script_dir/cd-hit-auxtools/cd-hit-dup -i $f2-stage0-rep.fa -o $f2-stage1 -d 0 -m false -e 3 > $f2-stage1.log");
  }
  else {
    nice_run("$script_dir/cd-hit/cd-hit-est -i $f2-stage0-rep.fa -o $f2-stage1 -c $c1 -n 10 -l 11 -p 1 -d 0 -g 1 -b 3 $cdhit_opt > $f2-stage1.log");
  }
  nice_run("$script_dir/cd-hit/clstr_rev.pl  $f2-stage0.clstr     $f2-stage1.clstr | $script_dir/cd-hit/clstr_sort_by.pl > $f2-stage1-all.clstr");
  nice_run("$script_dir/clstr_sort_rep.pl    $f2-stage1-all.clstr $f2-stage1 > $f2-stage1-rep.fa");
}
if (not $debug_mode) {
  $no_clstr = count_clstrs_from_clstr_file("$f2-stage1.clstr");
  print LOG "Stage1 clustering at $c1\n";
  print LOG "Number_clusters_stage1\t$no_clstr\n";
}

################################################################################
#### Stage 2 ---------- clustering at 98.50%  #### distance 1.50%
################################################################################
   $c1 = 0.985;
if ($restart_n <= 2) {
  if ($fast_mode) {
    nice_run("$script_dir/cd-hit-auxtools/cd-hit-dup -i $f2-stage1-rep.fa -o $f2-stage2 -d 0 -m false -e 6 > $f2-stage2.log");
  }
  else {
    nice_run("$script_dir/cd-hit/cd-hit-est -i $f2-stage1-rep.fa -o $f2-stage2 -c $c1 -n 10 -l 11 -p 1 -d 0 -g 1 -b 3 $cdhit_opt > $f2-stage2.log");
  }
  nice_run("$script_dir/cd-hit/clstr_rev.pl  $f2-stage1-all.clstr $f2-stage2.clstr | $script_dir/cd-hit/clstr_sort_by.pl > $f2-stage2-all.clstr");
  nice_run("$script_dir/clstr_sort_rep.pl    $f2-stage2-all.clstr $f2-stage2 > $f2-stage2-rep.fa");
}
if (not $debug_mode) {
  $no_clstr = count_clstrs_from_clstr_file("$f2-stage2.clstr");
  print LOG "Stage2 clustering at $c1\n";
  print LOG "Number_clusters_stage2\t$no_clstr\n";
}


################################################################################
#### Stage pre-3 ---------- filtering 
################################################################################

if ($restart_n <= 3) {
  nice_run("$script_dir/clstr_select_rep.pl size $cutoff 999999999 < $f2-stage2-all.clstr > $f2-stage2-rep-big.ids");
  nice_run("$script_dir/fetch_fasta_by_ids.pl      $f2-stage2-rep-big.ids  $f2-stage2-rep.fa   > $f2-stage2-rep-big.fa");
  nice_run("$script_dir/fetch_fasta_exclude_ids.pl $f2-stage2-rep-big.ids  $f2-stage2-rep.fa   > $f2-stage2-rep-small.fa");

  if (-s $clstr2) {
    nice_run("$script_dir/clstr_select_rep.pl size 1 999999999 < $clstr2 > $dir/chimaric.ids"); ## save chimaric ids
    nice_run("$script_dir/fetch_fasta_exclude_ids.pl $dir/chimaric.ids $f2-stage2-rep-big.fa > $f2-stage2-rep-big-good.fa"); ## exclude chimaric reads from $t1-pri-rep.fa
    nice_run("rm -f $f2-stage2-rep-big.fa");

    nice_run("$script_dir/fetch_fasta_exclude_ids.pl $dir/chimaric.ids $f2-stage2-rep-small.fa > $f2-stage2-rep-small-good.fa");
    nice_run("rm -f $f2-stage2-rep-small.fa");
  }
  else {
    nice_run("mv $f2-stage2-rep-big.fa   $f2-stage2-rep-big-good.fa");
    nice_run("mv $f2-stage2-rep-small.fa $f2-stage2-rep-small-good.fa");
  }
}

if (not $debug_mode) {
  print LOG "Min_clstr_size\t$cutoff\n";
  my $no_seq = count_seqs_from_fasta_file("$f2-stage2-rep-big-good.fa");
  print LOG "Number_clstrs_above_min_size\t$no_seq\n";
}

################################################################################
#### Stage 3 ---------- clustering at 97%
################################################################################
   $c1 = $otu_cutoff;
if ($restart_n <= 3) {
  nice_run("$script_dir/cd-hit/cd-hit-est -i $f2-stage2-rep-big-good.fa -o $f2-stage3 -c $c1 -n  8 -l 11 -p 1 -d 0 -g 1 -b 5 $cdhit_opt > $f2-stage3.log");
  nice_run("$script_dir/cd-hit/clstr_rev.pl  $f2-stage2-all.clstr $f2-stage3.clstr | $script_dir/cd-hit/clstr_sort_by.pl > $f2-stage3-all.clstr");
  nice_run("$script_dir/clstr_sort_rep.pl    $f2-stage3-all.clstr $f2-stage3 > $f2-stage3-rep.fa");
  nice_run("mv -f $f2-stage3-all.clstr $dir/OTU.clstr");
  nice_run("$script_dir/cd-hit-otu-table-faa.pl -i $dir/OTU.clstr -s $f2-stage3-rep.fa -o $dir/OTU-dist.txt -f $dir/OTU.fa");
}

if (not $debug_mode) {
  $no_clstr = count_clstrs_from_clstr_file("$dir/OTU.clstr");
  $no_seq   = count_seqs_from_clstr_file("$dir/OTU.clstr");
  print LOG "OTU clustering at $c1\n";
  print LOG "Number_OTUs\t$no_clstr\n";
  print LOG "Number_seqs_in_OTUs\t$no_seq\n";
  my ($tu,$ts,$cu,$cs)=times(); my $tt=$tu+$ts+$cu+$cs;
  print LOG "Total_CPU_time\t$tt\n";
}
close(LOG);


sub usage {
<<EOF
Usage:
$script_name -i contig_fasta_file -o output_dir -a abundance_cutoff -c OTU_cutoff -m check_chimera_flag

Options:
    -i input fasta file of contig
    -o output dir
    -c OTU cutoff, default 0.97
    -m whether to perform chimera checking (true/false), default true
    -a abundance cutoff, default 0.00005
       small clusters < this size will be considiered as noise and will be removed
       if total input sequence is 50,000, then clusters < 2 (i.e. singletons) are removed
    -f 1 or 0, default 0
       if set to 1, then use cd-hit-dup instead of cd-hit-est for stage 1 and 2 clustering
       which is very fast
    -R restart flag, if re-run at different abundance cutoff value or something,
       with this parameter, program can skip the first n step and restart at certain step
       values:
       0    default, start from the scratch cd-hit-dup
       1    cd-hit-est at 99.25
       2    cd-hit-est at 98.75
       3    filtering and cd-hit-est at 97%

EOF
}
###### END usage

sub nice_run {
  my $str = shift;
  print STDERR "$str\n";
  my $cmd = `$str` unless ($debug_mode);
  return $cmd;
}
##########

sub count_clstrs_from_clstr_file {
  my $clstr = shift;
  my $n = `grep -c "^>" $clstr`;
     $n =~ s/\s//g;
  return $n;
}

sub count_seqs_from_clstr_file {
  my $clstr = shift;
  my $n = `grep -cv "^>" $clstr`;
     $n =~ s/\s//g;
  return $n;
}

sub count_seqs_from_fasta_file {
  my $faa  = shift;
  my $n = `grep -c "^>" $faa`;
     $n =~ s/\s//g;
  return $n;
}

