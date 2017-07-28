#!/usr/bin/perl
#
use Getopt::Std;
getopts("i:s:S:o:f:j:",\%opts);

my $input             = $opts{i}; $input   = "Cluster.clstr" unless $input;
my $output            = $opts{o}; $output  = "Cluster.txt" unless ($output);
my ($i, $j, $k, $str, $cmd, $ll);

my %count = ();
my %count_t = ();
my %count_s = ();
my $Cluster_2_ann = ();
# >4360486|k__Bacteria;.p__Firmicutes;.c__Clostridia;.o__Clostridiales;.f__Lachnospiraceae;.g__Roseburia;.s__faecis
open(TMP, $input) || die "can not open $input";
my $Cluster=0;
while($ll=<TMP>){
  if ($ll =~ /^>/) {
    $Cluster++;
  }
  else {
    chop($ll);
    if ($ll =~ /\d+(aa|nt), >(.+)\.\.\./) {
      my $id = $2;
      if ($id =~ /^Sample\|([^\|]+)\|/) {
        $sample_id = $1;
        $sample_id{$sample_id}=1;
        $count{$Cluster}{$sample_id}++;
        $count_t{$Cluster}++;
        $count_s{$sample_id}++;
      }
      else {
        $Cluster_2_ann{$Cluster} .= "$id\t";
      }
    }
    else {
      die "format error $ll";
    }
  }
}
close(TMP);

my @sample_ids = sort keys %sample_id;

open(OUT1, "> $output") || die "can not write $output";
print OUT1 "Cluster";
foreach $sample_id (@sample_ids){
  print OUT1 "\t$sample_id";
}
#print OUT1 "\tTotal\n";
print OUT1 "\tAnnotation\tSpike\n";

for ($i=1; $i<=$Cluster; $i++){
  $ann = "None";
  if ($Cluster_2_ann{$i}) { $ann = $Cluster_2_ann{$i}; }
  print OUT1 "Cluster$i";
  foreach $sample_id (@sample_ids){
    $k = $count{$i}{$sample_id}? $count{$i}{$sample_id} : 0;
    print OUT1 "\t$k";
  }
  #print OUT1 "\t$count_t{$i}";
  print OUT1 "\t$ann\n";
}
close(OUT1);


