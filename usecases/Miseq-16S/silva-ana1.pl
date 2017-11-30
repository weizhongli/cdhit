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
die usage() unless ($opts{j} and $opts{o});
my ($i, $j, $k, $cmd);
my ($ll, $lla, $llb, $id, $ida, $idb, $seq, $seqa, $seqb, $qua, $quaa, $quab);
my ($len, $lena, $lenb);

my $fasta = $opts{j};
my $output  = $opts{o};

my %id_2_seq = ();
my $id = "";
my $ann;
open(TMP, $fasta) || die "can not open $fasta";
while($ll=<TMP>){
  if ($ll =~ /^>/) {
    chop($ll);
    ($id, $ann) = split(/\s+/, substr($ll,1), 2);
    $ann =~ s/\s/_/g;
    $id = "$id|$ann";
  }
  else {
    $ll =~ s/U/T/g;
    $ll =~ s/u/T/g;
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
This script formats Silva FASTA file for CD-HIT-OTU-MiSeq. You should download Silva sequences
from https://www.arb-silva.de/fileadmin/silva_databases/release_128/Exports/SILVA_128_SSURef_Nr99_tax_silva.fasta.gz
or similar file, gunzip the file

Run this script as $0 -j SILVA_128_SSURef_Nr99_tax_silva.fasta -o silva_128_SSURef_processed.fasta

Options:
======================
        -j path for SILVA_128_SSURef_Nr99_tax_silva.fasta
        -o output FASTA file of formatted Silva reference DB
EOD
}
