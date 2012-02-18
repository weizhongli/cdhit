#!/usr/bin/perl

$clstr = shift;
$fr = shift; # for nr80.clstr $fr = 0.8
$fra = 1 - $fr;

print "(\n";
open(TMP, $clstr) || die;
my $ll;
my $no = 0;
my @ids = ();
my $rep = "";
while($ll=<TMP>) {
  if ($ll =~ /^>/) {
    if ($no) {
      if ($no ==1 ) { print "$rep:1.0,\n";}
      else {
        my @mms = ();
        foreach my $tid (@ids) {
          push(@mms, "$tid:$fra");
        }
        my $mm = join(",\n", @mms);
        print "(\n$mm\n):$fr,\n";
      }
    }
    $no = 0;
    @ids = ();
    $rep = "";
  }
  else {
    my $tid= "";
    if ($ll =~ /(aa|nt), >(.+)\.\.\./) { $tid=$2;}
    push(@ids,$tid);
    if ($ll =~ /\*/) {$rep = $tid;}
    $no++;
  }
}
    if ($no) {
      if ($no ==1 ) { print "$rep:1.0\n";}
      else {
        my @mms = ();
        foreach my $tid (@ids) {
          push(@mms, "$tid:$fra");
        }
        my $mm = join(",\n", @mms);
        print "(\n$mm\n):$fr\n";
      }
    }
close(TMP);

print ");\n";
