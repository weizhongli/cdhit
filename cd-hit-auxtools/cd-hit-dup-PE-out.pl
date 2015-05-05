#!/usr/bin/perl

my $script_name = $0;
my $script_dir = $0;
   $script_dir =~ s/[^\/]+$//;
   chop($script_dir);
   $script_dir = "./" unless ($script_dir);


use Getopt::Std;
getopts("i:j:o:p:c:",\%opts);
die usage() unless ($opts{i} and $opts{j} and $opts{o} and $opts{p} and $opts{c} );
my ($i, $j, $k, $cmd);
my ($ll, $lla, $llb, $id, $ida, $idb, $seq, $seqa, $seqb, $qua, $quaa, $quab);
my ($len, $lena, $lenb);

my %FHZ=();

my $fastq1        = $opts{i};
my $fastq2        = $opts{j};
my $clstr_file    = $opts{c};
my $outfile1      = $opts{o};
my $outfile2      = $opts{p};
my $format        = input_test($fastq1); #fasta or fastq

my %rep_ids = ();
open(TMP, $clstr_file) || die "can not open $clstr_file";
while($ll = <TMP>){
  if ($ll =~ /^>/) {
    next;
  }
  else {
    chop($ll);
    if ($ll =~ /\*$/) {
      if ($ll =~ /\s(\d+)(aa|nt), >(.+)\.\.\./) {
        my $id = $3;
        $rep_ids{$id} = 1;
      }
    }
  }
}
close(TMP);

open(INPUTa, $fastq1) || die "can not open $fastq1"; open(OUTa, "> $outfile1") || die "can not write to $outfile1";
open(INPUTb, $fastq2) || die "can not open $fastq2"; open(OUTb, "> $outfile2") || die "can not write to $outfile2";

if ($format eq "fasta") {
  my $flag;
  while($lla = <INPUTa> and $llb=<INPUTb>) {
    if ( ($lla =~ /^>/) and ($llb =~ /^>/) ) {
      my $ida = substr($lla,1); chop($ida); $ida =~ s/\s.+$//;
      my $idb = substr($llb,1); chop($idb); $idb =~ s/\s.+$//;
      $flag = (($rep_ids{$ida}) or ($rep_ids{$idb})) ? 1 : 0;
    }
    if ($flag) {
      print OUTa $lla;
      print OUTb $llb;
    }
  }
}
else {
  my $flag;
  while($lla = <INPUTa> and $llb=<INPUTb>) {
    if ( ($lla =~ /^@/) and ($llb =~ /^@/) ) {
      my $ida = substr($lla,1); chop($ida); $ida =~ s/\s.+$//;
      my $idb = substr($llb,1); chop($idb); $idb =~ s/\s.+$//;
      $flag = (($rep_ids{$ida}) or ($rep_ids{$idb})) ? 1 : 0;
      if ($flag) {
        print OUTa $lla;
        print OUTb $llb;
        $lla = <INPUTa>; print OUTa $lla;
        $lla = <INPUTa>; print OUTa $lla;
        $lla = <INPUTa>; print OUTa $lla;
        $llb = <INPUTb>; print OUTb $llb;
        $llb = <INPUTb>; print OUTb $llb;
        $llb = <INPUTb>; print OUTb $llb;
      }
    }
  }
}

close(INPUTa); close(OUTa);
close(INPUTb); close(OUTb);


sub usage {
<<EOD
This script exports the representative PE reads into two seperate files after running
cd-hit-dup
 
         -i fasta or fastq file of PE read 1
         -j fasta or fastq file of PE read 2
         -c .clstr file produced by cd-hit-dup
         -o output file of representative reads, PE read 1
         -p output file of representative reads, PE read 2

EOD
}
######### END usage

sub input_test {
  my $f = shift;
  open(TTT, $f) || die "can not open $f\n";
  my $ll = <TTT>;
  close(TTT);

  my $c = substr($ll,0,1);
  if ($c eq ">") {return "fasta";}
  else           {return "fastq";}
}

