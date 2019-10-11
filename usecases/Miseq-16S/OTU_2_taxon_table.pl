#!/usr/bin/env perl
## =========================== NGS tools ==========================================
## NGS tools for metagenomic sequence analysis
## May also be used for other type NGS data analysis
##
##                                      Weizhong Li, UCSD
##                                      liwz@sdsc.edu
## http://weizhongli-lab.org/
## ================================================================================

use Getopt::Std;
getopts("i:o:a:t:r:N:c:P:",\%opts);
die usage() unless ($opts{o} and $opts{i} and $opts{t});


my $otu_table   = $opts{i}; ### e.g. OTU-short.txt
my $taxon_file  = $opts{t}; ### e.g. OTU-feature.txt
my $output      = $opts{o};

my ($i, $j, $k, $ll, $cmd);

my @samples = ();
my @otus = ();
my %otu_mat = ();

my $fh;
if ($otu_table eq "-") { $fh = "STDIN";}
else {
  open(TMP, $otu_table) || die "can not open $otu_table";
  $fh = "TMP";
}

$ll = <$fh>; chop($ll);
my ($t1, @lls) = split(/\t/, $ll);
@samples = @lls;
my $num_samples = $#samples+1;

while($ll=<$fh>){
  next if ($ll =~ /^#/);
  next unless ($ll =~ /^otu/i);
  chop($ll);
  my ($otu, @v) = split(/\t/, $ll);
  push(@otus, $otu);
  for ($i=0; $i<$num_samples; $i++) {
    $otu_mat{$otu}{$samples[$i]} = $v[$i];
  }
}

open(TMP, $taxon_file) || die "can not open $taxon_file";
my %taxon_info = ();
while($ll=<TMP>) {
  chop($ll);
  next if ($ll =~ /^#/);
  my ($otu, $taxon, $c) = split(/\t/, $ll);
  # next unless ($taxon =~ /__/); #### skip unknown OTUs

#OTUID	taxonomy	confidence
#OTU1	Root;k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__	1.0
#OTU2	Root;k__Bacteria;p__TM7;c__TM7-3;o__CW040;f__F16;g__;s__	1.0
#OTU3	Root;k__Bacteria;p__Proteobacteria;c__Epsilonproteobacteria;o__Campylobacterales;f__Helicobacteraceae;g__Flexispira;s__rappini	1.0
#OTU4	Root;k__Bacteria;p__Proteobacteria;c__Deltaproteobacteria;o__Desulfovibrionales;f__Desulfovibrionaceae;g__Desulfovibrio;s__C21_c20	1.0

  my $k = "unclassified";
  my $p = "unclassified";
  my $c = "unclassified";
  my $o = "unclassified";
  my $f = "unclassified";
  my $g = "unclassified";
  my $s = "unclassified";

  $j = $taxon;
  if ($j =~ /^k__([^;]+)/) {$k = $1;}
  if ($j =~ /;k__([^;]+)/) {$k = $1;}
  if ($j =~ /;p__([^;]+)/) {$p = $1;}
  if ($j =~ /;c__([^;]+)/) {$c = $1;}
  if ($j =~ /;o__([^;]+)/) {$o = $1;}
  if ($j =~ /;f__([^;]+)/) {$f = $1;}
  if ($j =~ /;g__([^;]+)/) {$g = $1;}
  if ($j =~ /;s__([^;]+)/) {$s = $1;}

  if ($j =~ /^D_0__([^;]+)/) {$k = $1;}
  if ($j =~ /;D_0__([^;]+)/) {$k = $1;}
  if ($j =~ /;D_1__([^;]+)/) {$p = $1;}
  if ($j =~ /;D_2__([^;]+)/) {$c = $1;}
  if ($j =~ /;D_3__([^;]+)/) {$o = $1;}
  if ($j =~ /;D_4__([^;]+)/) {$f = $1;}
  if ($j =~ /;D_5__([^;]+)/) {$g = $1;}
  if ($j =~ /;D_6__([^;]+)/) {$s = $1;}

  if (($g ne "unclassified") and ($s ne "unclassified")) {
    if ( substr($s, 0, length($g)) ne $g) { #### if species name doesn't contain genus name, add
      $s = "$g $s";
    }
  } 
  $taxon_info{$otu} = [$k,$p,$c,$o,$f,$g,$s];

}
close(TMP);

my @ranks = qw/kingdom phylum class order family genus species/;
my %rank_col = qw/kingdom 0 phylum 1 class 2 order 3 family 4 genus 5 species 6/;
foreach $rank (@ranks) {
  next if ($rank eq "kingdom");

  my $c = $rank_col{$rank};
  my $out = "$output.$rank.txt";
  open(OUT, "> $out") || die "can not write to $out";
  #### print table header
  print OUT "#", join("\t", @ranks[0..$c]);
  print OUT "\t", join("\t", @samples), "\n";

  my %rank_ti_info = ();
  my %rank_mat = ();
  my %ti_sum = ();
  foreach $otu (@otus) {
    my @ann = @{$taxon_info{$otu}};
    my $ti = join("|", @ann[0 .. $c] );

    if (not defined($rank_ti_info{$ti})) {
      $rank_ti_info{$ti} = [ @ann[0 .. $c] ];
    }
    foreach $sample (@samples) {
      $rank_mat{$ti}{$sample} += $otu_mat{$otu}{$sample}; 
      $ti_sum{$ti} += $otu_mat{$otu}{$sample};
    }
  }
  my @tis = keys %rank_mat;
     @tis = sort {$ti_sum{$b} <=> $ti_sum{$a} } @tis;

  foreach $ti (@tis) {
    print OUT join("\t", @{ $rank_ti_info{$ti} } );
    foreach $sample (@samples) {
      print OUT "\t", $rank_mat{$ti}{$sample};
    }
    print OUT "\n";
  }
  close(OUT);
}


sub usage {
<<EOD
Given cd-hit-otu pipeline output, this script generate
taxonomy abundance tables

usage:
  $0  -i otu_table -t taxon_info_file -o output

  options
    -i input OTU table from cd-hit-otu 16S pipeline, default STDIN
       e.g. OTU-short.txt
    -t taxon info, e.g. OTU-features.txt
    -o output file

EOD
}
########## END usage

