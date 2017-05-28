#!/usr/bin/perl

use Getopt::Std;
my $script_name = $0;
my $script_dir = $0;
   $script_dir =~ s/[^\/]+$//;
   chop($script_dir);
   $script_dir = "./" unless ($script_dir);

getopts("k:i:j:o:p:c:s:t:m:e:Z:a:f:d:R:g:",\%opts);
die usage() unless ($opts{k} and $opts{i} and $opts{j} and $opts{a} and $opts{f} and $opts{g} and $opts{o});

my $input0     = $opts{k}; ## nr.clstr
my $input      = $opts{i}; ## R1 only clstr
my $input2     = $opts{j}; ## R2 only clstr
my $clstr_99   = $opts{a}; ## seq.99.clstr         #### can be any 2nd -preclustering e.g. 98.5%
my $seq_99     = $opts{f}; ## seq.99 - fasta file R1
my $seq_992    = $opts{g}; ## seq.99 - fasta file R2
my $output     = $opts{o}; ## seq.99f
my $abs_cutoff = $opts{c}; $abs_cutoff = 0.0001  unless ($abs_cutoff);
my $output_2   = "$output.2";     ## seq.99f.2 -- R2
my $output_cls = "$output.clstr"; ## seq.99f.clstr
my $output_log = "$output.log";   ## seq.99f.log

my ($i, $j, $k, $str, $cmd, $ll);

my $num_total_seq;
my %seq_nr_size; 
my %seqs_of_nr;
open(LOG, "> $output_log") || die "can not open $output_log";
open(TMP, $input0) || die "can not open $input0";
if (1) {
  my $rep;
  while($ll=<TMP>){
    if ($ll =~ /^>/) {
      $rep = "";
    }
    else {
      chop($ll);
      my $id;
      if ($ll =~ /\d+(aa|nt), >(.+)\.\.\./) {
        $id = $2;
        $num_total_seq++;
        if ($ll =~ /\*$/) { $rep=$id; $seq_nr_size{$rep}=0; $seqs_of_nr{$rep} = [];}
        $seq_nr_size{$rep}++ if ($rep);
        push(@{$seqs_of_nr{$rep}}, $id) if ($rep);
      }
    }
  }
}
close(TMP);

my %seq_R1_clstr;
my %seq_R2_clstr;
foreach my $f (($input, $input2)) {
  open(TMP, $f) || die "can not open $f";
  my $rep;

  while($ll=<TMP>){
    if ($ll =~ /^>/) {
      $rep = "";
    }
    else {
      chop($ll);
      my $id;
      if ($ll =~ /\d+(aa|nt), >(.+)\.\.\./) {
        $id = $2;
        if ($ll =~ /\*$/) {
          $rep=$id;
        }
        if ($rep) {
          if ($f eq $input) { $seq_R1_clstr{$id} = $rep;}
          else              { $seq_R2_clstr{$id} = $rep;}
        }
      }
    }
  }
  close(TMP);
}

#### open $clstr_99 first time
open(TMP, $clstr_99) || die "can not open $clstr_99";
%rep_2_otu = ();
$OTU = -1;
while($ll=<TMP>){
  if ($ll =~ /^>/) {
    $OTU++;
  }
  else {
    my $id;
    if ($ll =~ /\d+(aa|nt), >(.+)\.\.\./) {
      $id = $2;
      $rep_2_otu{$id} = $OTU;
    }
  }
}
close(TMP);

my %chimeric_ids = (); 
#### those ids are candidates, if they are recurited by other non-chimeric clusters, 
#### then they are not chimeric anymore
foreach $i (keys %seq_R1_clstr) {
  my $rep1 = $seq_R1_clstr{$i};
  my $rep2 = $seq_R2_clstr{$i};

  next if ($rep1 eq $rep2);
  next unless ($seq_nr_size{$rep1} >= $seq_nr_size{$i}*2);
  next unless ($seq_nr_size{$rep2} >= $seq_nr_size{$i}*2);

  my $OTU1 = $rep_2_otu{$rep1};
  my $OTU2 = $rep_2_otu{$rep2};
  next if ($OTU1 eq $OTU2);
  $chimeric_ids{$i} = 1;
}

#### parse seq.99.clstr
my $cutoff_clstr_size = int($num_total_seq * $abs_cutoff);
   $cutoff_clstr_size = 1 unless ($cutoff_clstr_size >= 1); #### singleton will be removed
#print LOG "cutoff_clstr_size\t$cutoff_clstr_size\n";

open(TMP, $clstr_99) || die "can not open $clstr_99";
open(OUT, "> $output_cls") || die "can not write to $output_cls";
my %good_ids = ();
my @seqs_this_cls = ();
if (1) {
  my $clstr_txt = "";
  my $clstr_size = 0;
  my $rep;

  while($ll=<TMP>){
    if ($ll =~ /^>/) {
      if ($clstr_txt) {
        if (($clstr_size > $cutoff_clstr_size) and (not $chimeric_ids{$rep})) {
          print OUT $clstr_txt;
          $good_ids{$rep} = 1;
        }
        elsif ( $chimeric_ids{$rep} ) {
          foreach $j (@seqs_this_cls) {
            foreach $i ( @{ $seqs_of_nr{$j} } ) {
              print LOG "$i\tChimeric_cluster\t$rep\t$clstr_size\tP1:$seq_R1_clstr{$rep}\tP2:$seq_R2_clstr{$rep}\tOTU1:$rep_2_otu{$seq_R1_clstr{$rep}}\tOTU2:$rep_2_otu{$seq_R2_clstr{$rep}}\n";
            }
          }
        }
        else {
          foreach $j (@seqs_this_cls) {
            foreach $i ( @{ $seqs_of_nr{$j} } ) {
              print LOG "$i\tSmall_cluster\t$rep\t$clstr_size\n";
            }
          }
        }
      }
      $clstr_size = 0;
      $clstr_txt = $ll;
      $rep = "";
      @seqs_this_cls=();
    }
    else {
      $clstr_txt .= $ll;
      chop($ll);
      my $id;
      if ($ll =~ /\d+(aa|nt), >(.+)\.\.\./) {
        $id = $2;
        $clstr_size += $seq_nr_size{$id};
        $rep=$id if ($ll =~ /\*$/);
        push(@seqs_this_cls, $id);
      }
    }
  }
      if ($clstr_txt) {
        if (($clstr_size > $cutoff_clstr_size) and (not $chimeric_ids{$rep})) {
          print OUT $clstr_txt;
          $good_ids{$rep} = 1;
        }
        elsif ( $chimeric_ids{$rep} ) {
          foreach $j (@seqs_this_cls) {
            foreach $i ( @{ $seqs_of_nr{$j} } ) {
              print LOG "$i\tChimeric_cluster\t$rep\t$clstr_size\tP1:$seq_R1_clstr{$rep}\tP2:$seq_R2_clstr{$rep}\tOTU1:$rep_2_otu{$seq_R1_clstr{$rep}}\tOTU2:$rep_2_otu{$seq_R2_clstr{$rep}}\n";
            }
          }
        }
        else {
          foreach $j (@seqs_this_cls) {
            foreach $i ( @{ $seqs_of_nr{$j} } ) {
              print LOG "$i\tSmall_cluster\t$rep\t$clstr_size\n";
            }
          }
        }
      }
}
close(TMP);
close(OUT);

foreach my $f (($seq_99, $seq_992)) {
  my $fout = ($f eq $seq_99) ? $output : $output_2;

  open(TMP, $f) || die "can not open $f";
  open(OUT, ">$fout") || die "can not write to $fout";

  my $flag = 0;
  while($ll = <TMP>) {
    if ($ll =~ /^>/) {
      $gi = substr($ll,1);
      chop($gi);
      $gi =~ s/\s.+$//;
      $flag = ( $good_ids{$gi} ) ? 1 : 0;
    }
    print OUT $ll  if ($flag);
  }

  close(TMP);
  close(OUT);
}


close(LOG);

sub usage {
<<EOF
Usage:
$script_name -k seq.nr.clstr -i seq.nr.R1.clstr -j seq.nr.R2.clstr -c 0.0001 -a seq.99.clstr -f seq.99 -g seq.99.2 -o seq.99f

Options:
    -k input seq.nr.clstr
    -i input seq.nr.R1.clstr
    -j input seq.nr.R2.clstr
    -a input seq.99.clstr
    -f input seq.99
    -g input seq.99.2
    -o output 
    -c abundance cutoff, default $abs_cutoff
       small clusters < this size will be considiered as noise and will be removed
       if total input sequence is 50,000, then clusters < 2 (i.e. singletons) are removed

EOF
}
###### END usage

