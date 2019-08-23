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
  if ($ll =~ /^>(\S+)/) {
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
This script formats Qiime comtatible Silva FASTA file for CD-HIT-OTU-MiSeq. You should download Silva sequences
from https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_xxx_release.zip
unzip this file

Run this script as $0 -i SILVA_132_QIIME_release/taxonomy/16S_only/99/consensus_taxonomy_7_levels.txt \\
  -j SILVA_132_QIIME_release/rep_set/rep_set_16S_only/99/silva_132_99_16S.fna -o silva_132_99_16S_processed.fna

Options:
======================
        -i path for SILVA_132_QIIME_release/taxonomy/16S_only/99/consensus_taxonomy_7_levels.txt
        -j path for SILVA_132_QIIME_release/rep_set/rep_set_16S_only/99/silva_132_99_16S.fna -o silva_132_99_16S_processed.fna
        -o output FASTA file of formatted Silva reference DB
EOD
}
