#!/usr/bin/perl

#my $by = shift;
my $min;
my $max;
#if ($by eq "size") {
  $min = shift;
  $max = shift;
#}

$rep = "";
$no = 0;

while($ll=<>){
  if ($ll =~ /^>/) {
    if (($no >= $min) and ($no <= $max)) {
      print "$rep\n";
    }
    $rep = "";
    $no = 0;
  }
  else {
    chop($ll);
    if ($ll =~ /\*$/) {
      $rep = "";
      if ($ll =~ /\d+(aa|nt), >(.+)\.\.\./) {
        $rep = $2;
      }
      else {
        die "format error $ll";
      }
    }
    $no++;
  }
}
    if (($no >= $min) and ($no <= $max)) {
      print "$rep\n";
    }
