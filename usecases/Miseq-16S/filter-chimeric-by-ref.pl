#!/usr/bin/perl

use Getopt::Std;
my $script_name = $0;
my $script_dir = $0;
   $script_dir =~ s/[^\/]+$//;
   chop($script_dir);
   $script_dir = "./" unless ($script_dir);

getopts("k:i:j:o:p:c:s:t:m:e:Z:a:f:d:R:g:",\%opts);
die usage() unless ($opts{i} and $opts{j} and $opts{a} and $opts{f} and $opts{g} and $opts{o});

my $input      = $opts{i}; ## R1 only clstr
my $input2     = $opts{j}; ## R2 only clstr
my $clstr_99   = $opts{a}; ## seq.97f-full.clstr        #### can be any 2nd -preclustering e.g. 98.5%
my $seq_99     = $opts{f}; ## seq.99 - fasta file R1
my $seq_992    = $opts{g}; ## seq.99 - fasta file R2
my $output     = $opts{o}; ## seq.99f
my $abs_cutoff = $opts{c}; $abs_cutoff = 0.01  unless ($abs_cutoff); #### small cluster will be checked for chimeric 
my $output_2   = "$output.2";     ## seq.99f.2 -- R2
my $output_cls = "$output.clstr"; ## seq.99f.clstr
my $output_log = "$output.log";   ## seq.99f.log

my ($i, $j, $k, $str, $cmd, $ll);

my $num_total_seq;
my %seq_nr_size;
my %seqs_of_rep;
open(LOG, "> $output_log") || die "can not open $output_log";
open(TMP, $clstr_99) || die "can not open $clstr_99";
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
        $num_total_seq++ if ($id =~ /^Sample/);
        if ($ll =~ /\*$/) { $rep=$id; $seq_nr_size{$rep}=0; $seqs_of_rep{$rep} = [];}
        $seq_nr_size{$rep}++ if ($rep);
        push(@{$seqs_of_rep{$rep}}, $id) if ($rep);
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
        if ($rep and ($id =~ /^Sample/) ) { 
          if ($f eq $input) { $seq_R1_clstr{$id} = $rep;}
          else              { $seq_R2_clstr{$id} = $rep;}
        }
      }
    }
  }
  close(TMP);
}

my $cutoff_clstr_size = int($num_total_seq * $abs_cutoff); 
   $cutoff_clstr_size = 1 unless ($cutoff_clstr_size >= 1); 
#print LOG "cutoff_clstr_size\t$cutoff_clstr_size\n";

my %chimeric_ids = (); 
#### those ids are candidates, if they are recurited by other non-chimeric clusters, 
#### then they are not chimeric anymore
foreach $i (keys %seq_nr_size) {
  next unless ($i =~ /^Sample/);
  my $rep1 = $seq_R1_clstr{$i};
  my $rep2 = $seq_R2_clstr{$i};
  next unless ($rep1 and $rep2);

  next if ($rep1 eq $rep2);
  next if ($rep1 eq $i);
  next if ($rep2 eq $i);
  next if ($seq_nr_size{$i} > $cutoff_clstr_size);
  if (defined($seq_nr_size{$rep1})) { next unless ($seq_nr_size{$rep1} >= $seq_nr_size{$i}*2); }
  if (defined($seq_nr_size{$rep2})) { next unless ($seq_nr_size{$rep2} >= $seq_nr_size{$i}*2); }

  $chimeric_ids{$i} = 1;
}

#### parse seq.97fwref.clstr
#### do chimeric checking for sample-only clusters 
open(TMP, $clstr_99) || die "can not open $clstr_99";
open(OUT, "> $output_cls") || die "can not write to $output_cls";
my %good_ids = ();
if (1) {
  my $clstr_txt = "";
  my $clstr_size = 0;
  my $rep;
  my $refonly = 1;

  while($ll=<TMP>){
    if ($ll =~ /^>/) {
      if ($clstr_txt) {
        if ( not $refonly ) {
          if (not $chimeric_ids{$rep}) {
            print OUT $clstr_txt;
            $good_ids{$rep} = 1;
          }
          elsif ( $chimeric_ids{$rep} ) {
            foreach $i ( @{ $seqs_of_rep{$rep} }) {
              print LOG "Chimeric_cluster\t$i\t$rep\t$clstr_size\tP1:$seq_R1_clstr{$rep}\tP2:$seq_R2_clstr{$rep}\n";
            }
          }
        }
      }
      $clstr_size = 0;
      $clstr_txt = $ll;
      $rep = "";
      $refonly = 1;
    }
    else {
      $clstr_txt .= $ll;
      chop($ll);
      my $id;
      if ($ll =~ /\d+(aa|nt), >(.+)\.\.\./) {
        $id = $2;
        $clstr_size++;
        $rep=$id if ($ll =~ /\*$/);
        $refonly = 0 if ($id =~ /^Sample/);
      }
    }
  }
      if ($clstr_txt) {
        if ( not $refonly ) {
          if (not $chimeric_ids{$rep}) {
            print OUT $clstr_txt;
            $good_ids{$rep} = 1;
          }
          elsif ( $chimeric_ids{$rep} ) {
            foreach $i ( @{ $seqs_of_rep{$rep} }) {
              print LOG "Chimeric_cluster\t$i\t$rep\t$clstr_size\tP1:$seq_R1_clstr{$rep}\tP2:$seq_R2_clstr{$rep}\n";
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
$script_name -i seq.nr.R1.clstr -j seq.nr.R2.clstr -c 0.0001 -a seq.97f-full.clstr -f seq.99 -g seq.99.2 -o seq.99f

Options:
    -i input seq.nr.R1.clstr
    -j input seq.nr.R2.clstr
    -a input seq.97f-full.clstr
    -f input seq.99
    -g input seq.99.2
    -o output cluster without chimeric cluster, without ref-only cluster
    -c abundance cutoff, default $abs_cutoff
       small clusters < this size will be checked for chimeric and be removed if is chimeric
       if total input sequence is 50,000, then clusters < 2 (i.e. singletons) are checked

EOF
}
###### END usage

