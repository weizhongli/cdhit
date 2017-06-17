#!/usr/bin/perl
#
use Getopt::Std;
getopts("i:s:S:o:f:j:",\%opts);

my $input             = $opts{i}; $input   = "OTU.clstr" unless $input;
my $output            = $opts{o}; $output  = "OTU.txt" unless ($output);
my ($i, $j, $k, $str, $cmd, $ll);

my %count = ();
my %count_t = ();
my %count_s = ();
my $OTU_2_ann = ();

open(TMP, $input) || die "can not open $input";
my $OTU=0;
while($ll=<TMP>){
  if ($ll =~ /^>/) {
    $OTU++;
  }
  else {
    chop($ll);
    if ($ll =~ /\d+(aa|nt), >(.+)\.\.\./) {
      my $id = $2;
      if ($id =~ /^Sample\|([^\|]+)\|/) {
        $sample_id = $1;
        $sample_id{$sample_id}=1;
        $count{$OTU}{$sample_id}++;
        $count_t{$OTU}++;
        $count_s{$sample_id}++;
      }
      else {
        $OTU_2_ann{$OTU} = $id;
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
print OUT1 "OTU";
foreach $sample_id (@sample_ids){
  print OUT1 "\t$sample_id";
}
#print OUT1 "\tTotal\n";
print OUT1 "\tAnnotation\n";

for ($i=1; $i<=$OTU; $i++){
  $ann = "None";
  if ($OTU_2_ann{$i}) { $ann = $OTU_2_ann{$i}; }
  print OUT1 "OTU$i";
  foreach $sample_id (@sample_ids){
    $k = $count{$i}{$sample_id}? $count{$i}{$sample_id} : 0;
    print OUT1 "\t$k";
  }
  #print OUT1 "\t$count_t{$i}";
  print OUT1 "\t$ann\n";
}
close(OUT1);


