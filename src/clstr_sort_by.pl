#!/usr/bin/perl

my $sort_by_what = shift;
   $sort_by_what = "no" unless $sort_by_what;

my @clstr = ();
my $no = 0;
my $len = 0;
my $clstr = "";

while($ll = <>) {
  if ($ll =~ /^>/) {
    if ($clstr) {
      push(@clstr, [$len, $no, $clstr]);
    }
    $no = 0;
    $len = 0;
    $clstr = "";
  }
  else {
    $clstr .= $ll;
    $no++;
    chop($ll);
    if ($ll =~ /(\d+)(aa|nt),/) {
      my $this_len = $1;
      if ($this_len > $len) {$len = $this_len;}
    }
  }
}
    if ($clstr) {
      push(@clstr, [$len, $no, $clstr]);
    }

if ($sort_by_what eq "no") {
  @clstr = sort {$b->[1] <=> $a->[1]} @clstr;
}
elsif ($sort_by_what eq "len") {
  @clstr = sort {$b->[0] <=> $a->[0]} @clstr;
}

my $clstr_no = 0;
foreach $c (@clstr) {
  print ">Cluster $clstr_no\n";
  print $c->[-1];
  $clstr_no++;
}

