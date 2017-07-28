#!/usr/bin/perl

use Getopt::Std;
my $script_name = $0;
my $script_dir = $0;
   $script_dir =~ s/[^\/]+$//;
   chop($script_dir);
   $script_dir = "./" unless ($script_dir);

getopts("i:s:o:f:c:",\%opts);
die usage() unless ($opts{i} and $opts{s} and $opts{o} and $opts{f});

my $input      = $opts{i}; ## nr.clstr
my $fa         = $opts{s}; ## R1.fa
my $output     = $opts{o}; ## nr-filtered.clstr
my $output_fa  = $opts{f}; ## nr-filtered
my $cutoff     = $opts{c}; $cutoff = 1  unless ($cutoff);

my ($i, $j, $k, $str, $cmd, $ll);

open(TMP, $input) || die "can not open $input";
open(OUT, "> $output") || die "can not write to $output";
$no = 0;
$clstr = "";
$rep = "";
my %good_ids = ();

while($ll=<TMP>){
  if ($ll =~ /^>/) {
    if ($no > $cutoff) {
      print OUT $clstr;
      $good_ids{$rep}=1;
    }
    $clstr = $ll;
    $no = 0;
  }
  else {
    $clstr .= $ll;
    chop($ll);
    if ($ll =~ /\*$/) {
      $rep = "";
      if ($ll =~ /\d+(aa|nt), >(.+)\.\.\./) {
        $rep = $2;
      }
      else {
        die "format error $ll";
      }
    }
    $no++;
  }
}
    if ($no > $cutoff) {
      print OUT $clstr;
      $good_ids{$rep}=1;
    }
close(OUT);
close(TMP);

open(TMP, $fa) || die "can not open $fa";
open(OUT, ">$output_fa") || die "can not write to $output_fa";

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



sub usage {
<<EOF
Usage:
$script_name -i seq.nr.clstr -s R1.fa -o output.clstr -f output.fa -c 1 

Options:
    -i input seq.nr.clstr
    -s R1.fa
    -o output.clstr
    -f output.fa
    -c abundance cutoff, default $cutoff
       small clusters <= this size will be considiered as noise and will be removed

EOF
}
###### END usage

