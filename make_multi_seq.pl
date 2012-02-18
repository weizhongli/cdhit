#!/usr/bin/perl

#note you have to use "-d 0" in the cd-hit run
#note you better to use "-g 1" in the cd-hit run

# this script read the .clstr file, it generate a seperate fasta file
# for each cluster over certain size
# for example, if you run cd-hit -i db -o dbout -c 0.6 -n 4 -d 0 -g 1
# then you will have a dbout.clstr
# now run this script:
# ./make_multi_seq.pl db dbout.clstr multi-seq 20
# 

my $fasta = shift;
my $clstr = shift;
my $out_dir = shift;
my $size_cutoff = shift;

die unless (-e $fasta);
die unless (-e $clstr);
die unless ($out_dir);
$size_cutoff = 1 unless ($size_cutoff);

if (not -e $out_dir) {my $cmd = `mkdir $out_dir`;}


open(TMP, $clstr) || die;
my %id2cid=();
my $cid = "";
my @ids = ();
my $no = 0;
while($ll = <TMP>) {
  if ($ll =~ /^>/) {
    if ($no >= $size_cutoff) {
      foreach $i (@ids) { $id2cid{$i} = $cid;}
    }

    if ($ll =~ /^>Cluster (\d+)/) { $cid = $1; }
    else { die "Wrong format $ll"; }
    @ids = ();
    $no = 0;
  }
  else {
    if ($ll =~ /(aa|nt), >(.+)\.\.\./) { push(@ids, $2); $no++; }
    else { die "Wrong format $ll"; }
  }
}
close(TMP);
    if ($no >= $size_cutoff) {
      foreach $i (@ids) { $id2cid{$i} = $cid;}
    }


open(FASTA, $fasta) || die;
my $outfile_open = 0;
my $flag = 0;
while($ll=<FASTA>) {
  if ($ll=~ /^>(\S+)/) {
    my $id = $1;
    my $cid = $id2cid{$id};
    if (defined($cid)) {
      close(OUT) if ($outfile_open);
      open(OUT, ">> $out_dir/$cid") || die "can not open file to write $out_dir/$cid";
      $outfile_open = 1;
      $flag = 1;
    }
    else {
      $flag = 0;
    }
  }
  if ($flag) {print OUT $ll;}
}
close(FASTA);
      close(OUT) if ($outfile_open);


