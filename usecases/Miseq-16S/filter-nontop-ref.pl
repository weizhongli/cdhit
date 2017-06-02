#!/usr/bin/perl

use Getopt::Std;
my $script_name = $0;
my $script_dir = $0;
   $script_dir =~ s/[^\/]+$//;
   chop($script_dir);
   $script_dir = "./" unless ($script_dir);

my ($i, $j, $k, $str, $cmd, $ll);

my $clstr = "";
my $best_ref = "";
my $best_score = 0;

my $refonly = 1;
while($ll=<>){
  if ($ll =~ /^>/) {
    if ($clstr) {
      print $clstr;
      print $best_ref if ($best_ref);
    }

    $clstr = $ll;
    $best_ref = "";
    $best_score = 0;
  }
  else {
    if ($ll =~ /\d+(aa|nt), >(.+)\.\.\./) {
      my $id = $2;
      if ($id =~ /^Sample/) {
        $clstr .= $ll;
      }
      elsif ( $ll =~ /\/([\d|\.]+)%$/) {
        my $iden = $1;
        if ($iden > $best_score) {
          $best_score = $iden;
          $best_ref = $ll;
        }
      }
    }
    else {
      print STDERR "format err: $ll";
    }
  }
}

    if ($clstr) {
      print $clstr;
      print $best_ref if ($best_ref);
    }
