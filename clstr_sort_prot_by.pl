#!/usr/bin/perl

my $sort_by = shift;
   $sort_by = "len" unless ($sort_by);


my @clstr = ();
my $readin = 0;
my $head = "";
while($ll=<>) {
  if ($ll =~ /^>/) {

    if ($readin) {
      print $head;
      if ($sort_by eq "len") {
        @clstr = sort { $b->[1]<=>$a->[1] } @clstr;
      }
      elsif ($sort_by eq "id") {
        @clstr = sort {$a->[2] cmp $b->[2] or $b->[1]<=>$a->[1]} @clstr;
      }
      else {
        @clstr = sort {$b->[1]<=>$a->[1]  or $a->[2] cmp $b->[2]} @clstr;
      }
      for ($i=0; $i<@clstr; $i++) {
        print "$i\t$clstr[$i]->[0]\n";
      }
    }
    @clstr = ();
    $head = $ll;
  }
  else {
    $readin=1;
    chop($ll);
    if ($ll =~ /(\d+)(aa|nt), >(.+)\.\.\./) {
      $len = $1;
      $id = $3;

      $len = 99999999 if ($ll =~ /\*/);
      $ll =~ s/^\d+\t//;
      push(@clstr, [$ll, $len, $id]);

    }
  }
}
    if ($readin) {
      print $head;
      if ($sort_by eq "len") {
        @clstr = sort { $b->[1]<=>$a->[1] } @clstr;
      }
      elsif ($sort_by eq "id") {
        @clstr = sort {$a->[2] cmp $b->[2] or $b->[1]<=>$a->[1]} @clstr;
      }
      else {
        @clstr = sort {$b->[1]<=>$a->[1]  or $a->[2] cmp $b->[2]} @clstr;
      }
      for ($i=0; $i<@clstr; $i++) {
        print "$i\t$clstr[$i]->[0]\n";
      }
    }
