#!/usr/bin/perl
## =========================== NGS tools ==========================================
## NGS tools for metagenomic sequence analysis
## May also be used for other type NGS data analysis
##
##                                      Weizhong Li, UCSD
##                                      liwz@sdsc.edu
## http://weizhongli-lab.org/
## ================================================================================

my $script_name = $0;
my $script_dir = $0;
   $script_dir =~ s/[^\/]+$//;
   chop($script_dir);
   $script_dir = "./" unless ($script_dir);

use Getopt::Std;
getopts("i:j:o:r:e:p:q:c:d:N:t:u:d:M:T:S:",\%opts);
die usage() unless ($opts{i} and $opts{j} and $opts{o} and $opts{d});
my ($i, $j, $k, $cmd);
my ($ll, $lla, $llb, $id, $ida, $idb, $seq, $seqa, $seqb, $qua, $quaa, $quab);
my ($len, $lena, $lenb);

my $fastq         = $opts{i};
my $fastq2        = $opts{j};
my $ref           = $opts{d};
my $output        = $opts{o};
my $trim_R1       = $opts{p}; $trim_R1 = 100 unless ($trim_R1);
my $trim_R2       = $opts{q}; $trim_R2 = 100 unless ($trim_R2);
my $clstr_cutoff  = $opts{c}; #### post clustering
my $full_frag     = $opts{S};
my $prime_len     = 45;
my $output_R1     = "$output-R1";
my $output_R2     = "$output-R2";
my $session       = "OTU-session-$$";
my $output_S      = "$output-single";
my $consensus_db  = "$output-consensus";
my $cd_hit_2d     = "$script_dir/../../cd-hit-est-2d";  die "no $cd_hit_2d"  unless (-e $cd_hit_2d);
my $cd_hit_est    = "$script_dir/../../cd-hit-est";     die "no $cd_hit_est" unless (-e $cd_hit_est);
my $format        = input_test($fastq); #fasta or fastq
my $cdhit_opt_M   = $opts{M}; $cdhit_opt_M = 16000 unless defined($cdhit_opt_M);

if (defined($clstr_cutoff)) {
  die "Clustering cutoff $clstr_cutoff is not reasonable, should be <=1.0 and >= 0.97" unless (($clstr_cutoff <=1.0) and ($clstr_cutoff>=0.97));
}

my %FHZ=();

my %ref_map = ();
foreach my $f (($fastq, $fastq2)) {
  my $R = ( $f eq $fastq ) ? "R1" : "R2";
  open(OUT, "> $consensus_db.$R") || die "can not write to $consensus_db.$R";

  my %con = ();
  my $num_seq = 0;
  open_files_z_safe("TTTa", $f);

  if ($format eq "fastq") {
    while(1) {
      ($lla, $ida, $seqa, $quaa, $lena) = read_next_fastq("TTTa");
      last unless ($lla);
      for ($i=0; $i<$prime_len; $i++) {
        $c=uc(substr($seqa, $i, 1));
        $con{$i}{$c}++;
      }
      $num_seq++;
    }
  }
  else { #### fasta
    my $seqa = "";
    while($ll = <TTTa>) {
      if ($ll =~ /^>/) {
        if ($seqa) {
          for ($i=0; $i<$prime_len; $i++) {
            $c=uc(substr($seqa, $i, 1));
            $con{$i}{$c}++;
          }
          $num_seq++;
        }
        chop($ll);
        $seqa = "";
      }
      else {
        chop($ll);
        $seqa .= $ll;
      }
    }
        if ($seqa) {
          for ($i=0; $i<$prime_len; $i++) {
            $c=uc(substr($seqa, $i, 1));
            $con{$i}{$c}++;
          }
          $num_seq++;
        }
  } #### END fasta

  close(TTTa);

  my @cons = ();   #which letter
  my @cons_v = (); #abundance
  for ($i=0; $i<$prime_len; $i++) {
     my %t = %{ $con{$i} };
     my @k = keys %t;
        @k = sort { $t{$b} <=> $t{$a} } @k;
     push(@cons,   $k[0]);
     push(@cons_v, $t{ $k[0] } / $num_seq);
  }
  ## set minimal consensus to be 30
  for ($i=33; $i<$prime_len; $i++) {
    if ( ($cons_v[$i  ] <0.75) and
         ($cons_v[$i-1] <0.75) and 
         ($cons_v[$i-2] <0.75) ) {
      $i = $i-2; last;
    } 
  }
  my $trim_len_new = $i;

  print OUT ">$R\n";
  for ($i=0; $i<$trim_len_new; $i++) {
    print OUT $cons[$i];
  }
  print OUT "\n";
  close(OUT);

  my $cmd_line = "$cd_hit_2d -i $consensus_db.$R -i2 $ref -d 0 -c 0.8 -n 5 -r 1 -p 1 -b 5 -o $session.$R-vs-ref -G 0 -A 30 -s2 0.01 -M $cdhit_opt_M > $session.$R-vs-ref.log";
  print "running $cmd_line\n";
  $cmd = `$cmd_line`;

  my $parse_template=<<EOD;
>Cluster 0
0 45nt, >R1... *
1 1479nt, >1111882... at 1:42:4:45/+/95.24%
2 1500nt, >1111856... at 1:42:4:45/+/88.10%
3 1426nt, >1111848... at 2:44:3:45/+/90.70%
4 1530nt, >1111847... at 1:42:4:45/+/85.71%
5 1497nt, >1111839... at 1:41:5:45/+/85.37%
6 1492nt, >1111819... at 1:42:4:45/+/88.10%
7 1482nt, >1111782... at 1:42:4:45/+/88.10%
8 1496nt, >1111776... at 1:42:4:45/+/88.10%
9 1500nt, >1111768... at 1:42:4:45/+/85.71%
...
>Cluster 0
0 45nt, >R2... *
1 1428nt, >1111883... at 483:440:2:45/-/84.09%
2 1479nt, >1111882... at 511:468:2:45/-/88.64%
3 1336nt, >1111879... at 435:399:2:38/-/86.49%
4 1402nt, >1111874... at 469:426:2:45/-/84.09%
5 1500nt, >1111856... at 513:470:2:45/-/84.09%
6 1530nt, >1111847... at 532:489:2:45/-/86.36%
7 1497nt, >1111839... at 509:473:2:38/-/86.49%
8 1492nt, >1111819... at 514:471:2:45/-/88.64%
9 1482nt, >1111782... at 502:464:2:40/-/84.62%
10  1496nt, >1111776... at 516:473:2:45/-/84.09%
EOD

  open(TMP, "$session.$R-vs-ref.clstr") || die "can not open $session.$R-vs-ref.clstr";
  while($ll=<TMP>){
    next if ($ll =~ /^>/);
    next if ($ll =~ /^0/);
    chop($ll);
    if ($ll =~ /^\d+\s+\d+(aa|nt), >(.+)\.\.\./) {
      my $id = $2;
      my @lls = split(/\s+/, $ll);
      my @lls = split(/\//, $lls[-1]); ##516:473:2:45/-/84.09%
      my ($query_start, $query_end, $rep_star, $rep_end) = split(/:/, $lls[0]);
      $ref_map{$id}{$R}=[$query_start, $query_end, $rep_star, $rep_end, $lls[1]];
    }
  }
  close(TMP); 
}

my %ref_cut;
foreach $id (keys %ref_map) {
  next unless (defined $ref_map{$id}{"R1"});
  next unless (defined $ref_map{$id}{"R2"});

  my @R1_info = @{$ref_map{$id}{"R1"}};
  my @R2_info = @{$ref_map{$id}{"R2"}};

  next unless ($R1_info[4] eq "+");
  next unless ($R2_info[4] eq "-");

  my $p1 = $R1_info[0] -  ($R1_info[2]-1);   #### 1-based, can be -1 value for V1
  my $p2 = $R2_info[0] +  ($R2_info[2]-1);   #### 1-based, can be longer than len($seq)
  $ref_cut{$id} = [$p1, $p2];
}

open(TMP, $ref) || die "can not open $ref";
open(OUT1, "> $output_R1") || die "can not write to $output_R1";
open(OUT2, "> $output_R2") || die "can not write to $output_R2";
if ($full_frag) {
  open(OUT3, "> $output_S") || die "can not write to $output_S";
}
my $seq;
my $des;
$id = "";

while($ll = <TMP>) {
  if ($ll =~ /^>/) {
    if ($seq) {
      if ($ref_cut{$id}) {
        $seq =~ s/\s//g;
        my ($p1, $p2) = @{$ref_cut{$id}};
        my $len = length($seq);
        my $seq1 = "";
        my $seq2 = "";
        if ($p1>=1) {
          $seq1 = substr($seq, $p1-1, $trim_R1);
        }
        else {
          my $pad = 1 - $p1; #### add NNN at 5'
          $seq1 = "N" x $pad;
          $seq1 .= substr($seq, 0, $trim_R1-$pad);
        }

        if ($p2 <= $len) {
          my $p2a = $p2 - $trim_R2;  #### 0 - based substr idx
          if ($p2a < 0) { #### not long enough
            $seq2 = substr($seq, 0, $p2);
          }
          else {
            $seq2 = substr($seq, $p2a, $trim_R2);
          }
        }
        else {  #### add NNN at 5'
          my $pad = $p2 - $len;
          my $trim_t2_t = $trim_R2 - $pad;
          $seq2 = "N" x $pad;

          my $p2a = $len - $trim_R2_t;  #### 0 - based substr idx
          if ($p2a < 0) { #### not long enough
            $seq2.= $seq;
          }
          else {
            $seq2 .= substr($seq, $p2a, $trim_R2_t);
          }
        }
        $seq2 = reverse_complement($seq2);
        ### now have $seq1 $seq2
        print OUT1 "$des loc=$p1 len=", length($seq1), "\n$seq1\n";
        print OUT2 "$des loc=$p2 len=", length($seq2), "\n$seq2\n";
        if ($full_frag) {
          if ($p1 < 1   ) {$p1 = 1;   }
          if ($p2 > $len) {$p2 = $len;}
          my $eff_len = $p2-$p1+1;
          my $seq1 = substr($seq, $p1-1, $eff_len);
          print OUT3 "$des loc=$p1:$p2 len=$eff_len\n$seq1\n";
        }
      }
    }
    chop($ll);
    $des = $ll;
    $id = substr($ll,1);
    $id =~ s/\s.+$//;
    $seq = "";
  }
  else {
    $seq .= $ll;
  }
}

    if ($seq) {
      if ($ref_cut{$id}) {
        $seq =~ s/\s//g;
        my ($p1, $p2) = @{$ref_cut{$id}};
        my $len = length($seq);
        my $seq1 = "";
        my $seq2 = "";
        if ($p1>=1) {
          $seq1 = substr($seq, $p1-1, $trim_R1);
        }
        else {
          my $pad = 1 - $p1; #### add NNN at 5'
          $seq1 = "N" x $pad;
          $seq1 .= substr($seq, 0, $trim_R1-$pad);
        }

        if ($p2 <= $len) {
          my $p2a = $p2 - $trim_R2;  #### 0 - based substr idx
          if ($p2a < 0) { #### not long enough
            $seq2 = substr($seq, 0, $p2);
          }
          else {
            $seq2 = substr($seq, $p2a, $trim_R2);
          }
        }
        else {  #### add NNN at 5'
          my $pad = $p2 - $len;
          my $trim_t2_t = $trim_R2 - $pad;
          $seq2 = "N" x $pad;

          my $p2a = $len - $trim_R2_t;  #### 0 - based substr idx
          if ($p2a < 0) { #### not long enough
            $seq2.= $seq;
          }
          else {
            $seq2 .= substr($seq, $p2a, $trim_R2_t);
          }
        }
        $seq2 = reverse_complement($seq2);
        ### now have $seq1 $seq2
        print OUT1 "$des loc=$p1 len=", length($seq1), "\n$seq1\n";
        print OUT2 "$des loc=$p2 len=", length($seq2), "\n$seq2\n";
        if ($full_frag) {
          if ($p1 < 1   ) {$p1 = 1;   }
          if ($p2 > $len) {$p2 = $len;}
          my $eff_len = $p2-$p1+1;
          my $seq1 = substr($seq, $p1-1, $eff_len);
          print OUT3 "$des loc=$p1:$p2 len=$eff_len\n$seq1\n";
        }
      }
    }

close(OUT1);
close(OUT2);
if ($full_frag) { close(OUT3); }
close(TMP);

if (defined($clstr_cutoff)) {
  my $output_R1_tmp = "$output_R1.$$";
  my $output_R2_tmp = "$output_R2.$$";

  my $cmd_line = "$cd_hit_est -i $output_R1 -j $output_R2 -d 0 -c $clstr_cutoff -n 10 -p 1 -b 5" .
                 " -o $output_R1_tmp -op $output_R2_tmp -G 1 -g 1 -M $cdhit_opt_M -P 1 -l 11 -sc 1       > $output_R1_tmp.log";
  print "running $cmd_line\n";
  $cmd = `$cmd_line`;

  die "Can not run $cd_hit_est" unless (-e "$output_R1_tmp.clstr");
  $cmd = `mv $output_R1_tmp $output_R1`;
  $cmd = `mv $output_R2_tmp $output_R2`;
  $cmd = `mv $output_R1_tmp.clstr $output.clstr`;

  if ($full_frag) {
    my $output_S_tmp = "$output_S.$$";
    my $cmd_line = "$cd_hit_est -i $output_S -d 0 -c $clstr_cutoff -n 10 -p 1 -b 5" .
                   " -o $output_S_tmp  -G 1 -g 1 -M $cdhit_opt_M -l 11 -sc 1 > $output_S_tmp.log";
    print "running $cmd_line\n";
    $cmd = `$cmd_line`;
    die "Can not run $cd_hit_est" unless (-e "$output_S_tmp.clstr");
    $cmd = `mv $output_S_tmp $output_S`;
    $cmd = `mv $output_S_tmp.clstr $output_S.clstr`;
  }
}

$cmd = `rm -f $session*`;

# need %FHZ
# open one or more files including zipped files
# above open_files_z may have broken pipe problem
# so this safe sub, open each file individually
sub open_files_z_safe {
  my ($fh, @files) = @_;
  my ($i, $j, $k);

  my $no = $#files+1;

  $FHZ{$fh} = {
    'files'      => [@files],
    'no'         => $no,
    'open_idx'   => 0,
  };

  my $f0 = $files[0];
  if    ($f0 =~ /\.gz$/ ) { open($fh, "gunzip -c $f0 |") || die "can not gunzip -c $f0\n"; }
  elsif ($f0 =~ /\.bz2$/) { open($fh, "bzcat     $f0 |") || die "can not bzcat $f0\n"; }
  else                    { open($fh,            $f0   ) || die "can not open $f0\n"; }
  return 0;
}
########## END open_files_z_safe


sub read_FHZ {
  my $fh = shift;
  my $ll;

  $ll = <$fh>;
  if ($ll) { return $ll;} ##### read from existing opened file

  #otherwise, last opened file reaches EOF
  if ($FHZ{$fh}->{open_idx} < $FHZ{$fh}->{no} -1 ) { ### still file not opened yet
    close($fh); #### close last open file

    $FHZ{$fh}->{open_idx}++;
    my $f0 = $FHZ{$fh}->{files}->[ $FHZ{$fh}->{open_idx} ];

    if    ($f0 =~ /\.gz$/ ) { open($fh, "gunzip -c $f0 |") || die "can not gunzip -c $f0\n"; }
    elsif ($f0 =~ /\.bz2$/) { open($fh, "bzcat     $f0 |") || die "can not bzcat $f0\n"; }
    else                    { open($fh,            $f0   ) || die "can not open $f0\n"; }

    $ll = <$fh>;
    return $ll;
  }
  else { #### no more file to open, return undef
    return undef;
  }
}
########### END read_FHZ


########## read_next_fastq
sub read_next_fastq {
  my $fh = shift;
  my ($lla, $seqa, $lla2, $quaa, $ida, $lena);
  $lla  = read_FHZ($fh); return unless ($lla);
  chop($lla); $lla =~ s/\s.+$//;
  $ida = substr($lla,1);
  $seqa  = read_FHZ($fh); $seqa =~ s/\s+$//g; $lena = length($seqa);
  $lla2  = read_FHZ($fh); #read ID
  $quaa  = read_FHZ($fh); $quaa =~ s/\s+$//g;
  return ($lla, $ida, $seqa, $quaa, $lena);
}
########## END read_next_fastq


sub reverse_complement {
    my ($in_seq) = @_;
    my $opposite = reverse $in_seq;
    $opposite =~ tr/ACGT/TGCA/;
    return("$opposite");
}


sub input_test {
  my $f = shift;
  open(TTT, $f) || die "can not open $f\n";
  my $ll = <TTT>;
  close(TTT);

  my $c = substr($ll,0,1);
  if ($c eq ">") {return "fasta";}
  else           {return "fastq";}
}
########## END input_test


sub usage {
<<EOD;
This script takes a paired-end (PE) read files (Fastq or Fasta) for a 16S dataset, e.g. from V3-V4
region, it also takes a Fasta file of full-length 16S reference database, e.g. Greengene.
this script identifies the sequencing region on the reference sequencs and it cuts the forward
and reverse segments and outputs them in PE fasta files. The output PE reference database can be used
to cluster together with 16S datasets

Options:
======================
        -i input fasta or fastq file for R1
        -j input fasta or fastq file for R2
        -d 16S reference sequence file in fasta format
        -o output prefix
        -p lenght of forward sequence in output file
        -q length of reverse sequence in output file
        -S also output full fragment
        -c cutoff for clustering the output PE files to remove redundant reference seqeunces. 
           Suggested cutoffs: 1.00, 0.99, 0.98 and 0.97
           The script will not cluster the output unless user specifies this cutoff.
        -M available memory to use, default 16000, means 16000MB. This option will be passed to cd-hit.
EOD
}
