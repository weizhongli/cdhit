#!/usr/bin/perl

$rep = "";
$no = 0;
while($ll=<>){
  if ($ll =~ /^>/) {
    if ($no) {
      print "$cid\t$rep\t$no\n";
    }
    $cid = "";
    if    ($ll =~ /^>Cluster (\d+)/) {$cid = $1;}
    elsif ($ll =~ /^>Clstr (\d+)/)   {$cid = $1;}
    else  {die "format error $ll"}
    $rep = "";
    $no = 0;
  }
  else {
    if ($ll =~ /\*$/) {
      $rep = "";
      if ($ll =~ /aa, >(.+)\.\.\./) {
        $rep = $1;
      }
      else {
        die "format error $ll";
      }
    }
    $no++;
  }
}
    if ($no) {
      print "$cid\t$rep\t$no\n";
    }
