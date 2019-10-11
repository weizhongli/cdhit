#!/usr/bin/env perl
#
use Getopt::Std;
getopts("i:s:S:o:f:j:b:h:",\%opts);

die usage() unless($opts{o});

my $input             = $opts{i}; $input   = "OTU.clstr" unless $input;
my $output            = $opts{o};
my $output_pre        = $output;
   $output_pre        =~ s/\.([^\.])+$//;
my $output_meta       = "$output_pre-sample-meta.txt";
my $output_feature    = "$output_pre-feature.txt";
my $output_short      = "$output_pre-short.txt";
my $output_biom       = "$output_pre.biom";
my $biom_exe          = $opts{b};


my ($i, $j, $k, $str, $cmd, $ll);

my %count = ();
my %count_t = ();
my %count_s = ();
my $OTU_2_ann = ();
my $tree_flag = 0;  #### for greengene header format
# >4360486|k__Bacteria;.p__Firmicutes;.c__Clostridia;.o__Clostridiales;.f__Lachnospiraceae;.g__Roseburia;.s__faecis
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
        $id =~ s/^([^\|]+)\|//;
        $id =~ s/;\./;/g;
        $id = "Root;$id";
        $id =~ s/;D_0__/;k__/;
        $id =~ s/;D_1__/;p__/;
        $id =~ s/;D_2__/;c__/;
        $id =~ s/;D_3__/;o__/;
        $id =~ s/;D_4__/;f__/;
        $id =~ s/;D_5__/;g__/;
        $id =~ s/;D_6__/;s__/;
        $OTU_2_ann{$OTU} = $id;
        $tree_flag = 1 if ($id =~ /;k__Bacteria/);
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
open(OUT2, "> $output_short") || die "can not write $output_short";
open(OUT3, "> $output_feature") || die "can not write $output_feature";

print OUT1 "#OTUID";
print OUT2 "#OTUID";
print OUT3 "#OTUID";
foreach $sample_id (@sample_ids){
  print OUT1 "\t$sample_id";
  print OUT2 "\t$sample_id";
}
if ($tree_flag) {
  print OUT1 "\t", join("\t", qw/Kingdom Phylum Class Order Family Genus Species/);
}
#print OUT1 "\tTotal\n";
print OUT1 "\tAnnotation\n";
print OUT2 "\n";
print OUT3 "\ttaxonomy\tconfidence\n";

for ($i=1; $i<=$OTU; $i++){
  $ann = "";
  if ($OTU_2_ann{$i}) { $ann = $OTU_2_ann{$i}; }
  print OUT1 "OTU$i";
  print OUT2 "OTU$i";
  print OUT3 "OTU$i";
  foreach $sample_id (@sample_ids){
    $k = $count{$i}{$sample_id}? $count{$i}{$sample_id} : 0;
    print OUT1 "\t$k";
    print OUT2 "\t$k";
  }
  if ($tree_flag) {
    my ($tax_k, $tax_p, $tax_c, $tax_o, $tax_f, $tax_g, $tax_s);
    if ($ann =~ /k__(\w+)/) {$tax_k = $1} else {$tax_k = "";}
    if ($ann =~ /p__(\w+)/) {$tax_p = $1} else {$tax_p = "";}
    if ($ann =~ /c__(\w+)/) {$tax_c = $1} else {$tax_c = "";}
    if ($ann =~ /o__(\w+)/) {$tax_o = $1} else {$tax_o = "";}
    if ($ann =~ /f__(\w+)/) {$tax_f = $1} else {$tax_f = "";}
    if ($ann =~ /g__(\w+)/) {$tax_g = $1} else {$tax_g = "";}
    if ($ann =~ /s__(\w+)/) {$tax_s = $1} else {$tax_s = "";}
    print OUT1 "\t", join("\t", ($tax_k, $tax_p, $tax_c, $tax_o, $tax_f, $tax_g, $tax_s));
  }
  #print OUT1 "\t$count_t{$i}";
  print OUT1 "\t$ann\n";
  print OUT2 "\n";
  print OUT3 "\t$ann\t1.0\n";
}
close(OUT1);
close(OUT2);
close(OUT3);

open(OUT, ">$output_meta") || die "can not write to $output_meta";
print OUT "#SampleID\tGroup\n";
foreach $sample_id (@sample_ids){
  print OUT "$sample_id\tnogroup\n";
}
close(OUT);

if (-e $biom_exe) {
  $cmd = `$biom_exe convert -i $output_short -o $output_biom --to-hdf5 --observation-metadata-fp $output_feature --sample-metadata-fp $output_meta`;
}

sub usage {
<<EOF
This script converts OTU clusters to tsv format, and if biom is available,
convert .biom file

Usage:
$script_name -i OTU.clstr -o OTU.txt

Options:
    -i input OTU.clstr, by cd-hit-otu-miseq
    -o output OTU table
    -b path to biom executible, if provided, will make .biom file

EOF
}

