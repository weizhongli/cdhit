#!/usr/bin/perl
# output single fasta file
# for each cluster output at least $cutoff seqs


$clstr_file = shift;
$fasta_file = shift;
$cutoff = shift;

die "" unless ((-e $clstr_file) and (-e $fasta_file) and ($cutoff>0));
my ($i, $j, $k, $cmd);
my %rep_ids = (); #rep_ids to be skipped


open(TMP, $clstr_file) || die "Can not open file";
my @gis = ();
my $clstr_id = "";
while(my $ll=<TMP>) {
  if ($ll =~ /^>/ ) {
    if (@gis) {
      my $no = $#gis+1;
      for ($i=$cutoff; $i<$no; $i++) { $rep_ids{$gis[$i]}=1; }
    }
    @gis = ();
  }
  else {
    chop($ll);
    if ($ll =~ />(.+)\.\.\./ ) {
      my $id = $1;
      if ($ll =~ /\*$/) { unshift(@gis, $id);}
      else              {    push(@gis, $id);}
    }
  }
}
    if (@gis) {
      my $no = $#gis+1;
      for ($i=$cutoff; $i<$no; $i++) { $rep_ids{$gis[$i]}=1; }
    }
close(TMP);


#########################################################

my $flag = 0;

open(TMP, $fasta_file) || die;
while($ll = <TMP>) {
    chomp($ll);$ll=~s/\r$//;
  if ($ll =~ /^>/) {
    $gi = substr($ll,1);
    $gi =~ s/\s.+$//;
    $flag = ($rep_ids{$gi}) ? 0 : 1;
  }
  if ($flag) {print "$ll\n";}
}
close(TMP);

