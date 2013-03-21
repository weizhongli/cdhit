#!/usr/bin/perl

my $no = 0;
my $clstr_no = "";
my @this_cluster = ();

print "id	clstr	clstr_size	length	clstr_rep	clstr_iden	clstr_cov\n";
while($ll = <>) {
  if ($ll =~ /^>/) {
    if ($no>0) {
      process_this_cluster();
    }
    if ($ll =~ /^>Cluster (\d+)/) {
      $clstr_no = $1;
    }
    $no = 0;
    @this_cluster = ();
  }
  else {
    my ($id, $len, $rep, $iden);
    if ($ll =~ /\d+\t(\d+)[a-z]{2}, >(.+)\.\.\. \*/) {
      $len = $1;
      $id  = $2;
      $rep = 1;
      $iden = 100;
      
    }
    elsif ($ll =~ /\d+\t(\d+)[a-z]{2}, >(.+)\.\.\./) {
      $len = $1;
      $id  = $2;
      $rep = 0;
      $ll=~/(\d+%|\d+\.\d+%)$/;
      $iden = $1;
    }
    else {
      print STDERR "***********\n";
    }
    push(@this_cluster, [($id, $len, $rep, $iden)]);
    $no++;
  }
}
    if ($no>0) {
      process_this_cluster();
    }

sub process_this_cluster {
  my ($i, $j, $k);
  my @t = sort { ($b->[2] <=> $a->[2]) or ( $b->[1] <=> $a->[1])  } @this_cluster;
  my $longest = 0;
  foreach $i (@t) {
    $longest = $i->[1] if ($i->[2]);
  }
  foreach $i (@t) {
    my $cov = int ( $i->[1]/$longest * 100);
    print "$i->[0]\t$clstr_no\t$no\t$i->[1]\t$i->[2]\t$i->[3]\t$cov\%\n";
  }
}
