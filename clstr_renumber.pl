#!/usr/bin/perl
$no = 0;
while($ll=<>){
  if ($ll =~ /^>Cluster (\d+)/) {
    print ">Cluster $no\n"; $no++;
    $cno  = 0;
  }
  else {
    $ll =~ s/^\d+/$cno/;
    print $ll;
    $cno++;
  }

}
