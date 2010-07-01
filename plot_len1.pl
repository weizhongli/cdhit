#!/usr/bin/perl

$file90 = shift;
$segs = shift;
@segs = split(/,/, $segs);
$len_segs = shift;
@len_segs = split(/,/,$len_segs);

my @clstr_nos = ();
my @clstr_len = ();
open(TMP, $file90) || die "Can not open file";
$readin = 0;
my $this_no = 0;
my $this_len = 0;
my $max_no = 0;

while(my $ll=<TMP>) {
  if ($ll =~ /^>/ ) {
    if ($readin) { 
      $clstr_nos[$this_no]++; 
      $max_no = $this_no if ($this_no>$max_no);
      if (not defined($clstr_len[$this_no])) {
        $clstr_len[$this_no] = [];
      }
      push(@{$clstr_len[$this_no]}, $this_len)
    }
    $this_no=0;
  }
  else {
    $readin = 1;
    $this_no++;
    chop($ll);
    if ($ll =~ /\*$/) {
      if ($ll =~ /(\d+)(aa|nt), /) { $this_len=$1;}
    }
  }
}
close(TMP);

if ($readin) { 
  $clstr_nos[$this_no]++; 
  $max_no = $this_no if ($this_no>$max_no);
  if (not defined($clstr_len[$this_no])) {
    $clstr_len[$this_no] = [];
  }
  push(@{$clstr_len[$this_no]}, $this_len)
}

print "Size\tNo. seq\tNo. clstr";
my @tlen_nos = ();
for ($j=0; $j<@len_segs; $j++) {
  $len_seg = $len_segs[$j];
  print "\t$len_seg";
  $tlen_nos[$j] = 0;
}
print "\n";
my $tno = 0;
my $tno1 = 0;
for ($i=0; $i<@segs; $i++) {
  $seg = $segs[$i];

  my @lens = ();
  if ($seg =~ /-/) {
    $no = 0;
    $no1 = 0;
    ($b, $e) = split(/-/, $seg);
    $e = $max_no if ($e =~ /up/i);
    for($j=$b; $j<=$e; $j++) {
      $no += $j * $clstr_nos[$j];
      $no1+= $clstr_nos[$j];
      push(@lens, @{$clstr_len[$j]});
    }
    $tno += $no; $tno1 += $no1;
    print "$seg\t$no\t$no1";
  }
  else {
    $tno += $seg * $clstr_nos[$seg];
    $tno1 += $clstr_nos[$seg];
    push(@lens, @{$clstr_len[$seg]});
    print "$seg\t", $seg * $clstr_nos[$seg], "\t$clstr_nos[$seg]";
  }

  for ($j=0; $j<@len_segs; $j++) {
    $len_seg = $len_segs[$j];
    $no1 = 0;
    my ($b, $e);
    if ($len_seg =~ /-/) {
      ($b, $e) = split(/-/, $len_seg);
    }
    else {
       $b = $e = $len_seg;
    }
    foreach $tlen (@lens) {
      $no1++ if (($tlen>=$b) and ($tlen<=$e));
    }
    print "\t$no1";
    $tlen_nos[$j] += $no1;
  }
  print "\n";
}
print "Total\t$tno\t$tno1";
for ($j=0; $j<@len_segs; $j++) {
  print "\t$tlen_nos[$j]";
}
print "\n";
