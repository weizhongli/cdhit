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
getopts("i:j:o:r:e:p:q:c:d:N:t:u:d:",\%opts);
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
my $prime_len     = 30;
my $output_R1     = "$output-R1";
my $output_R2     = "$output-R2";
my $consensus_db  = "$$-consensus";
my $cd_hit_2d     = "$script_dir/../../cd-hit-2d";

my %FHZ=();


open(OUT, "> $consensus_db") || die "can not write to $consensus_db";

foreach my $f (($fastq, $fastq2)) {

  my %con = ();
  my $num_seq = 0;
  open_files_z_safe("TTTa", $f);
  while(1) {
    ($lla, $ida, $seqa, $quaa, $lena) = read_next_fastq("TTTa");
    last unless ($lla);
    for ($i=0; $i<$prime_len; $i++) {
      $c=uc(substr($seqa, $i, 1));
      $con{$i}{$c}++;
    }
    $num_seq++;
  }
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
  ## set minimal consensus to be 15
  for ($i=15; $i<$prime_len; $i++) {
    last if ($cons_v[$i] <0.75);
  }
  my $trim_len_new = $i;

  my $name = ($f eq $fastq) ? "R1" : "R2";
  print OUT ">$name\n";
  for ($i=0; $i<$trim_len_new; $i++) {
    print OUT $cons[$i];
  }
  print OUT "\n";
}

close(OUT);


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



