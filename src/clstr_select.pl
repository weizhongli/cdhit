#!/usr/bin/perl

#my $by = shift;
my $min;
my $max;
#if ($by eq "size") {
  $min = shift;
  $max = shift;
#}

$no = 0;
$clstr = "";

while($ll=<>){
  if ($ll =~ /^>/) {
    if (($no >= $min) and ($no <= $max)) {
      print $clstr;
    }
    $clstr = $ll;
    $no = 0;
  }
  else {
    $clstr .= $ll;
    $no++;
  }
}
    if (($no >= $min) and ($no <= $max)) {
      print $clstr;
    }
