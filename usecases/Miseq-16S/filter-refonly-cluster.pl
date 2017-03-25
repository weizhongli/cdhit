#!/usr/bin/perl

use Getopt::Std;
my $script_name = $0;
my $script_dir = $0;
   $script_dir =~ s/[^\/]+$//;
   chop($script_dir);
   $script_dir = "./" unless ($script_dir);

my ($i, $j, $k, $str, $cmd, $ll);

my $num_total_seq;
my %seq_nr_size; 

if (1) {
  my $clstr = "";
  my $refonly = 1;
  while($ll=<>){
    if ($ll =~ /^>/) {
      print $clstr unless ($refonly);
      $clstr = $ll;
      $refonly = 1;
    }
    else {
      $clstr .= $ll;
      my $id;
      if ($ll =~ /\d+(aa|nt), >(.+)\.\.\./) {
        $id = $2;
        $refonly = 0 if ($id =~ /^Sample/);
      }
    }
  }
}

