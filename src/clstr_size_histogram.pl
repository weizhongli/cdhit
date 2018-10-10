#!/usr/bin/perl

if(@ARGV==0){
   print "Usage:\n\tclstr_size_histogram.pl [-bin N] clstr_file\n";
   exit(1);
}

#$file90 = shift;
$step = 100;
if($ARGV[0] eq "-bin"){$step=$ARGV[1];$file90=$ARGV[2];}
else{$file90=$ARGV[0];}

my @clstr_nos = ();

open(TMP, $file90) || die "Can not open file";
$readin = 0;
my $this_no = 0;
while(my $ll=<TMP>) {
  if ($ll =~ /^>/ ) {
    $this_no = int ( ($this_no-1) / $step);
    if ($readin) { $clstr_nos[$this_no]++; }
    $this_no=0;
  }
  else {
    $readin = 1;
    $this_no++;
  }
}
close(TMP);

$this_no = int ( ($this_no-1) / $step);
$clstr_nos[$this_no]++;

print "bin_size\tNo_of_clusters\n";
for ($i=0; $i<@clstr_nos; $i++) {
  if(!$clstr_nos[$i]){$clstr_nos[$i]=0;}
  print $i*$step+1, "-", $i*$step+$step,"\t$clstr_nos[$i]\n";
}
