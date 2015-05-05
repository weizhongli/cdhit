#!/usr/bin/perl

my ($gi_file, $seq_file) = @ARGV;

my %gi1 = ();
my $gi;

open(TMP, $gi_file) || die;
while($ll = <TMP>) {
  chop($ll);
  $ll =~ s/\s.+$//;
  $ll =~ s/^>//;
  $gi1{$ll} = 1;
}
close(TMP);

my $flag = 0;
open(TMP, $seq_file) || die;
while($ll = <TMP>) {
  if ($ll =~ /^>/) {
    $gi = substr($ll,1);
    chop($gi);
    $gi =~ s/\s.+$//;
    $flag = ( $gi1{$gi} ) ? 1 : 0;
  }
  print $ll  if ($flag);
}
close(TMP);

