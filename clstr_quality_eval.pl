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
## for example: ">id1||COG00001";

my %pair_status = ();
my %bench_clusters = ();
my %cdhit_clusters = ();
my %seqs = ();
my $seq_idx = 0;

my @clstr = ();
my $readin = 0;
my $i = 0;
while($ll=<>) {
  if ($ll =~ /^>/) {

    if ($readin) {
      if (@clstr > 1) { #skip singleton
        $cdhit_clusters{$i} = [@clstr];
        $i++;
      }
    }
    @clstr = ();
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
        push(@{$bench_clusters{$ben_id}},$seq_idx);
      }
      push(@clstr, $seq_idx);
      push(@seqs,$seq_id); $seq_idx++;
    }
  }
}
    if ($readin) {
      if (@clstr > 1) { #skip singleton
        $cdhit_clusters{$i} = [@clstr];
        $i++;
      }
    }

#process pairs in cdhit 
my $total_cdhit_pairs=0;
foreach my $family (keys %cdhit_clusters) {
  my @seqs = @{$cdhit_clusters{$family}};
  @seqs = sort {$a <=> $b} @seqs;
  my $N = $#seqs+1;

  for ($i=0; $i<$N-1; $i++) {
    my $idi = $seqs[$i];
    for ($j=$i+1; $j<$N; $j++){
      my $idj = $seqs[$j];
      $pair_status{"$idi|$idj"} = "C"; #cdhit pairs
      $total_cdhit_pairs++;
    }
  }
}

#process pairs in benchmark
my $total_bench_pairs = 0;
my $correct_pairs=0;
foreach my $family (keys %bench_clusters) {
  my @seqs = @{$bench_clusters{$family}};
  next unless (@seqs > 1); #skip singleton

  @seqs = sort {$a <=> $b} @seqs;
  my $N = $#seqs+1;

  for ($i=0; $i<$N-1; $i++) {
    my $idi = $seqs[$i];
    for ($j=$i+1; $j<$N; $j++){
      my $idj = $seqs[$j];

      if (defined( $pair_status{"$idi|$idj"} )) {
        $correct_pairs++;
        $pair_status{"$idi|$idj"} = "T"; #True pair
      }
      else {
        $pair_status{"$idi|$idj"} = "B"; #bench pairs
      }
      $total_bench_pairs++;
    }
  }
}



my $sen = $correct_pairs / $total_bench_pairs;
my $spe = $correct_pairs / $total_cdhit_pairs;

print "Total benchmark pairs\t$total_bench_pairs\n";
print "Total cd-hit pairs\t$total_cdhit_pairs\n";
print "Total correct pairs\t$correct_pairs\n";
print "Sensitivity\t$sen\n";
print "Specificity\t $spe\n";

print "\n\nPairs in benchmark but not in cd-hit\n";
foreach $pair (keys %pair_status) { 
  next unless ($pair_status{$pair} eq "B");
  my ($idx1, $idx2) = split(/\|/,$pair);
  print "$seqs[$idx1]\t$seqs[$idx2]\n";
}

print "\n\nPairs in cd-hit but not in benchmark\n";
foreach $pair (keys %pair_status) { 
  next unless ($pair_status{$pair} eq "C");
  my ($idx1, $idx2) = split(/\|/,$pair);
  print "$seqs[$idx1]\t$seqs[$idx2]\n";
}

print "\n\nPairs in both cd-hit and benchmark\n";
foreach $pair (keys %pair_status) { 
  next unless ($pair_status{$pair} eq "T");
  my ($idx1, $idx2) = split(/\|/,$pair);
  print "$seqs[$idx1]\t$seqs[$idx2]\n";
}


