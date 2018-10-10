#!/usr/bin/perl

#not like cd-hit-div, this script do not sort input
#or throw away seq

$in  = shift; $in  or die "no input file";
$out = shift; $out or die "no output file";
$div = shift; $div or die "no div number";

@fhs = ();
for ($i=0; $i<$div; $i++) {
  my $tf = "$out-$i";
  $tfh = "FH". $i;
  open($tfh, "> $tf") || die "can not open output files for write";
  push(@fhs, $tfh);
}

open(TMP, $in) || die "can not open input file";
my $seq;
my $des;
my $no = 0;
while($ll = <TMP>) {
  if ($ll =~ /^>/) {
    if ($seq) {
      $i = $no % $div;
      $tfh = $fhs[$i];
      $no++;
      print $tfh $des, $seq;
    }
    $des = $ll;
    $seq = "";
  }
  else {
    $seq .= $ll;
  }
}
    if ($seq) {
      $i = $no % $div;
      $tfh = $fhs[$i];
      $no++;
      print $tfh $des, $seq;
    }

close(TMP);

for ($i=0; $i<$div; $i++) {
  $tfh = $fhs[$i];
  close($tfh);
}
