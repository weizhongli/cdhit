#!/usr/bin/perl
# if nr90 from nr100 and
#    nr80 from nr90, so I have nr90.clstr and nr80.clstr
# but, in nr80.clstr, some gi numbers whose from nr100 are there
# use this script, I create a new nr80new.clstr, as it is clustered from nr100

$file90 = shift;
$file80 = shift;

my %gi2clstr = ();
open(TMP, $file90) || die "Can not open file";
$readin = 0;
my $gi = "";
my $clstr = "";
my $this_no = 0;
while(my $ll=<TMP>) {
  if ($ll =~ /^>/ ) {

    if ($readin and $gi and $this_no>1 ) {
      $gi2clstr{$gi} = $clstr;
    }
    $gi="";
    $clstr="";
    $this_no=0;
  }
  else {
    $readin = 1;
    $clstr .= $ll;
    if ($ll =~ /\*/ and $ll =~ />(.+)\.\.\./ ) { $gi = $1; }
    $this_no++;
  }
}
close(TMP);
if ($readin and $gi and $this_no>1 ) {
  $gi2clstr{$gi} = $clstr;
}

my $no = 0;
open(TMP, $file80) || die "Can not open file";
while( $ll = <TMP>) {
  if ($ll =~ /^>/ ) {
    print $ll;
    $no = 0;
  }
  elsif ($ll =~ />(.+)\.\.\./ ) {
    $gi = $1; 
    chop($ll);
    $rep = ( $ll =~ /\*$/) ? 1 : 0;
    $iden = "";
    if ($ll =~ / at (.+)$/) { $iden = $1; }
    else                  { $iden = "100%"; }

    if ( $gi2clstr{$gi} ) {
      $aa = $gi2clstr{$gi};
      @aa = split(/\n/, $aa);

      foreach $a (@aa) {
        $a =~ s/^\d+/$no/;
        if (not $rep) {
          if ($a =~ /\*$/) {
            $a =~ s/\*/at $iden/;
          }
          else {
            $a =~ s/at (.+)$/at $iden,$1/;
          }
        }
        print "$a\n";
        $no++;
      }
    }
    else { 
     $ll =~ s/^\d+/$no/;
     print "$ll\n";
     $no++;
    }
  }
  else {
    print $ll;
  }
}

close(TMP);

