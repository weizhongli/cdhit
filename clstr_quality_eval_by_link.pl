#!/usr/bin/perl

## calculate the sensitivity and specificity of clusters
## if the input fasta file has pre-defined classification term
## such as COG, pfam, TAXID etc.
## these terms can be used as benchmark to calcuate
## the sensitivity and specificity of clustering

## sensitivity is the ratio of number of correct pairs in cd-hit to number of pairs in benchmark
## specificity is the ratio of number of correct pairs in cd-hit to number of all cd-hit pairs
## here a pair is a pair of seqs that in same cluster

## the fasta file defline should be ">id||term"
## such as ">id1||COG00001";

## NOTE !!!!!!!!!!
## this script use independant link instead of pair
## for example
## if a cluster has 8 sequences, it has 28 pairs, but only 7 pairs are independant
## these 7 pairs are independ links

my %bench_clusters = ();
my $total_bench_links = 0;
my $total_cdhit_links = 0;
my $correct_links     = 0;

my %clstr_by_ben = ();
my $t_no = 0;
my $readin = 0;
while($ll=<>) {
  if ($ll =~ /^>/) {
    if ($readin) {
      if ($t_no > 1) {
        foreach my $ben_id (keys %clstr_by_ben) {
          $correct_links += $clstr_by_ben{$ben_id}-1;
        }
        $total_cdhit_links += $t_no-1;
      }
    }
    $t_no = 0;
    %clstr_by_ben = ();
  }
  else {
    $readin=1;
    chop($ll);
    if ($ll =~ /(\d+)(aa|nt), >(.+)\.\.\./) {
      $id = $3;
      my ($seq_id, $ben_id) = split(/\|\|/, $id);
      if ($ben_id) {
        if (not defined($bench_clusters{$ben_id})) {
          $bench_clusters{$ben_id} = [];
        }
        push(@{$bench_clusters{$ben_id}},$seq_id);
        $clstr_by_ben{$ben_id}++;
        $t_no++ ;
      }
    }
  }
}
      if ($t_no > 1) {
        foreach my $ben_id (keys %clstr_by_ben) {
          $correct_links += $clstr_by_ben{$ben_id}-1;
        }
        $total_cdhit_links += $t_no-1;
      }

#process links in benchmark
foreach my $family (keys %bench_clusters) {
  my @seqs = @{$bench_clusters{$family}};
  my $N = $#seqs+1;
  $total_bench_links += $N-1;
  #print "$family $N\n";
}

my $sen = $correct_links / $total_bench_links;
my $spe = $correct_links / $total_cdhit_links;

print "Total benchmark links\t$total_bench_links\n";
print "Total cd-hit links\t$total_cdhit_links\n";
print "Total correct links\t$correct_links\n";
print "Sensitivity\t$sen\n";
print "Specificity\t $spe\n";



