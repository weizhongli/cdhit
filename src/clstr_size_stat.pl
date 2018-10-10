#!/usr/bin/perl

if(@ARGV==0){
   print "Usage:\n\tclstr_size_stat.pl clstr_file\n";
   exit(1);
}

$file90 = shift;

my @clstr_nos = ();
my $max_size = 0;
open(TMP, $file90) || die "Can not open file";
$readin = 0;
my $this_no = 0;
while(my $ll=<TMP>) {
  if ($ll =~ /^>/ ) {
    if ($readin) { 
      $clstr_nos[$this_no]++; 
      if ($this_no > $max_size) {$max_size=$this_no;} 
    }
    $this_no=0;
  }
  else {
    $readin = 1;
    $this_no++;
  }
}
close(TMP);

$clstr_nos[$this_no]++;
if ($this_no > $max_size) {$max_size=$this_no;} 

print "size	No.clstr	No.seq\n";
for ($i=1; $i<=$max_size; $i++) {
  $noc = $clstr_nos[$i] ? $clstr_nos[$i] : 0;
  next unless $noc;
  $nos = $noc * $i;
  print "$i	$noc	$nos\n";
}

