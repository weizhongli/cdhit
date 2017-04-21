#!/usr/bin/perl
## =========================== NGS tools ==========================================
## NGS tools for metagenomic sequence analysis
## May also be used for other type NGS data analysis
##
##                                      Weizhong Li, UCSD
##                                      liwz@sdsc.edu
## http://weizhongli-lab.org/
## ================================================================================

use Getopt::Std;
getopts("i:j:o:r:e:p:q:c:d:N:t:u:d:M:T:S:",\%opts);
die usage() unless ($opts{i} and $opts{j} and $opts{o});
my ($i, $j, $k, $cmd);
my ($ll, $lla, $llb, $id, $ida, $idb, $seq, $seqa, $seqb, $qua, $quaa, $quab);
my ($len, $lena, $lenb);

my $file1 = $opts{i};
my $fasta = $opts{j};
my $output  = $opts{o};

my %id_2_ann;
open(TMP, $file1) || die "can not open $file1";
while($ll=<TMP>){
  chop($ll);
  my ($id, $txt) = split(/\s+/, $ll, 2);
  $txt =~ s/ /./g;
  $id_2_ann{$id} = $txt;
}
close(TMP);

my %id_2_seq = ();
my $id = "";
open(TMP, $fasta) || die "can not open $fasta";
while($ll=<TMP>){
  if ($ll =~ /^>(\d+)/) {
    chop($ll);
    $id = $1;
    $ann = $id_2_ann{$id};
    $id = "$id|$ann" if ($ann);
  }
  else {
    $id_2_seq{$id} .= $ll;
  }
}

close(TMP);

my @ids = keys %id_2_seq;
   @ids = sort {length($b) <=> length($a) } @ids;

open(OUT, "> $output") || die "can not write to $output";
foreach $id (@ids) {
  print OUT ">$id\n$id_2_seq{$id}";
}
close(OUT);



sub usage {
<<EOD;
This script formats Greengene FASTA file for CD-HIT-OTU-MiSeq. You should download Greengene sequences
from http://greengenes.secondgenome.com/downloads, or ftp://greengenes.microbio.me/.
download file like greengenes_release/gg_13_5/gg_13_5_otus.tar.gz, unpack the tar file. You may find
gg_13_5_otus/taxonomy/99_otu_taxonomy.txt and gg_13_5_otus/rep_set/99_otus.fasta

Run this script as $0 -i gg_13_5_otus/taxonomy/99_otu_taxonomy.txt -j gg_13_5_otus/rep_set/99_otus.fasta -o gg_13_5_processed.fasta

Options:
======================
        -i path for gg_13_5_otus/taxonomy/99_otu_taxonomy.txt
        -j path for gg_13_5_otus/rep_set/99_otus.fasta
        -o output FASTA file of formatted Greengene reference DB
EOD
}
