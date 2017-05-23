#!/usr/bin/perl

use Storable;
use strict;

my $input_file = shift;
my $output_file = shift;
my $sort_by_what = shift;
   $sort_by_what = "no" unless $sort_by_what;

my @clstr = values %{retrieve($input_file)};


if ($sort_by_what eq "no") {

  ### Added by liwz sort by No. sequences instead of No. nodes
  my %rep2size = ();
  my $clstr_no = scalar(@clstr);
  my ($i);


  for ($i=0; $i<$clstr_no; $i++){
    my $node_size = 0;
    foreach my $seq1 (@{$clstr[$i][3]}) {
      if ($$seq1[2]) { # can be futher expanded
        foreach my $seq2(@{$$seq1[3]}) {
          if ($$seq2[2]) { $node_size += scalar(@{$$seq2[3]}); }
          else           { $node_size++; }
        }
      }
      else {
        $node_size++;
      }
    }
    $rep2size{ $clstr[$i][0] } = $node_size;
  }
  ### END

  #@clstr = sort {scalar(@{$b->[3]}) <=> scalar(@{$a->[3]})} @clstr;
  @clstr = sort {$rep2size{$b->[0]} <=> $rep2size{$a->[0]}} @clstr;
}
elsif ($sort_by_what eq "len") {
  @clstr = sort {$b->[1] <=> $a->[1]} @clstr;
}
elsif ($sort_by_what eq "des") {
  @clstr = sort {$a->[0] cmp $b->[0]} @clstr;
}

store \@clstr, $output_file;


