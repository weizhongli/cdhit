#!/usr/bin/perl

#keep only top $no proteins in cluster

my $no_cutoff = shift;
$no_cutoff or die "no number\n";

my $clstr= "";
my $no = 0;
my $repout=1;

while($ll=<>){
  if ($ll =~ /^>/) {
    if ($no) {
      print $clstr;
    }
    $clstr = $ll;
    $no = 0;
    $repout=1;
  }
  else {
    if ($no < $no_cutoff-$repout or $ll =~ /\*$/) {
      if($ll=~/\*$/){$repout=0;}
      $clstr .= $ll;
      $no++;
    }
  }
}
    if ($no) {
      print $clstr;
    }
