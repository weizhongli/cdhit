#!/usr/bin/perl
#

my $rep;
my @non_reps = ();
my @aln_info = ();
while($ll=<>){
  if ($ll =~ /^>/ ) {
    if (@non_reps) {
      for ($i=0; $i<@non_reps; $i++) {
        print "$non_reps[$i]\t$rep\t$aln_info[$i]\n";
      }
    }
    $rep = "";
    @non_reps = ();
    @aln_info = ();
  }
  else {
    chop($ll);
    if ($ll =~ />(.+)\.\.\./ ) {
      my $id = $1;
      if ($ll =~ /\*$/) { $rep = $id}
      else              { 
        push(@non_reps, $id);
        my @lls = split(/\s+/, $ll);
        my ($a, $iden) = split(/\//, $lls[-1]);
        chop($iden); ### removing % sign
        my ($qb, $qe, $sb, $se) = split(/:/, $a);
        my $alnln = $qe-$qb+1;
        my $mis = int($alnln * (100-$iden) / 100);
        my $gap = 0;
        my $exp = 0;
        my $bit = $alnln*3 - $mis*6; #### very rough
        push(@aln_info, "$iden\t$alnln\t$mis\t$gap\t$qb\t$qe\t$sb\t$se\t$exp\t$bit");
      }
    }
#>Cluster 582
#0       6671aa, >gi|302514050|uid|51... *
#1       339aa, >SRS584717|scaffold|... at 2:332:4020:4356/89.32%
#2       182aa, >SRS584717|scaffold|... at 1:182:6490:6671/100.00%
#3       367aa, >SRS584717|scaffold|... at 1:332:4543:4874/90.66%
#4       109aa, >SRS584717|scaffold|... at 1:108:5782:5889/97.22%

  }
}
    if (@non_reps) {
      for ($i=0; $i<@non_reps; $i++) {
        print "$non_reps[$i]\t$rep\t$aln_info[$i]\n";
      }
    }

#query                          subject         %       alnln   mis     gap     q_b     q_e     s_b     s_e     expect  bits
##0                              1               2       3       4       5       6       7       8       9       10      11
#mHE-SRS012902|scaffold|86.16    gnl|CDD|226997  47.62   42      17      2       164     201     210     250     5e-04   37.6
#mHE-SRS012902|scaffold|109.23   gnl|CDD|225183  47.46   236     122     1       1       236     475     708     1e-92    284
#mHE-SRS012902|scaffold|109.23   gnl|CDD|224055  44.35   239     130     2       1       239     332     567     2e-84    259
#mHE-SRS012902|scaffold|109.23   gnl|CDD|227321  39.50   238     140     3       1       238     324     557     9e-69    218

