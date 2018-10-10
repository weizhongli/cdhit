#!/usr/bin/perl

use Storable;
use strict;
#my $sort_by_what = shift;
#  $sort_by_what = "no" unless $sort_by_what;

my $clstr_file = shift;
my $store_file = shift;

my %clstr = ();  # an array of hashes for all the cluster
my $rep_len = 0;
my $rep_acc = "";
my @cur_sequences = (); # array of hashes for all sequences in a cluster
my $ll = "";
my @record = ();

open(TMP, $clstr_file) || die;
while($ll = <TMP>) { # read .clstr files
        if ($ll =~ /^>/) { # the begin of a cluster
                if (scalar(@cur_sequences)) {  # not the first cluster, therefore collect the information of last clstr
                        #@cur_sequences = sort {$$b{"seq_len"} <=> $$a{"seq_len"}} @cur_sequences;
                        @cur_sequences = sort {$$b[1] <=> $$a[1]} @cur_sequences;
                        @record = ($rep_acc, $rep_len, 1, [@cur_sequences], "");
                        $clstr{$rep_acc} = [@record];
                }
    @cur_sequences=();
        }
        else { # the sequence line
                chop($ll);
                if ($ll =~ /^(\d+)\s+(\d+)(aa|nt),\s+>(.+)\.\.\./) {
                        @record = ($4, $2, 0, [], "");
                        if ($ll =~ /\*$/) { # representative sequence or not
                                $rep_acc = $record[0];
                                $rep_len = $record[1];
                                $record[4] = "100%";
                        }
#                       elsif ($ll =~ / at (\d.+)$/ ) { 
            elsif ($ll =~ / at (.+\d.+)$/ ) {# because cd-hit-est have strand info 
                                $record[4] = $1;
                        }
                }
                push(@cur_sequences, [@record]);
        }
}
if (scalar(@cur_sequences)) {
        #@cur_sequences = sort {$$b{"seq_len"} <=> $$a{"seq_len"}} @cur_sequences;
        @cur_sequences = sort {$$b[1] <=> $$a[1]} @cur_sequences;
        @record = ($rep_acc, $rep_len, 1, [@cur_sequences], "");
        $clstr{$rep_acc} = [@record];
}
close(TMP);

if (-e $store_file){ # already have a cluster file
        my %old_clstr = %{retrieve($store_file)};
        foreach my $rep_acc (keys %clstr){
                my $seqs = $clstr{$rep_acc}[3]; # $seqs a reference to the sequences;
                my $tmp_size = scalar(@{$seqs});  # how many sequences in a top level cluster, each sequence should be a representative sequence for lower level cluster
                #print "$rep_acc, $tmp_size\n";
                my $i;
                for $i  (0..($tmp_size-1)){
                        my $seq = $$seqs[$i];
                        if ($old_clstr{$$seq[0]}){
                                $clstr{$rep_acc}[3][$i][3] = [@{$old_clstr{$$seq[0]}[3]}];
                                $clstr{$rep_acc}[3][$i][2] = 1;
                        }
                }
        }
}

store \%clstr, $store_file;

#~ my $size = scalar(keys %clstr);
#~ print "$size\n";

#~ my $acc = 'D8F4YGO02FSTQP|range|2:370|frame|2|len|123';

#~ my $temp = $clstr{$acc}[1];
#~ print "$temp\n";

#~ my $temp = scalar(@{$clstr{$acc}[3]});
#~ print "$temp\n";

#~ my $x;
#~ for $x (@{$clstr{$acc}[3]} ){
        #~ my $tmp_1 = scalar(@{$x->[3]});
        #~ print "$x->[2], $x->[4], $x->[0], $x->[1], $tmp_1\n";
#~ }

