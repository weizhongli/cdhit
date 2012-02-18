#!/usr/bin/perl

# order of clusters don't need to be the same
# but then I have to read everything into memory

my ($master_clstr, @clstr) = @ARGV;
my $clstr_file_no = $#clstr+1;

my %slave_clstr = ();
foreach $file (@clstr) {
  my $rep = "";
  my $no  = "";
  my @members = ();
  open(TC, $file) || die;
  while($ll=<TC>){
    if ($ll =~ /^>/) {
      if ($no) {
        die "format error, no rep before cluster $ll" unless ($rep);
        if (not defined($slave_clstr{$rep})) {
          $slave_clstr{$rep}=[];
        }
        push(@{$slave_clstr{$rep}}, @members);
      }
      $rep = "";
      $no  = "";
      @members = ();
    }
    else {
      my $id = "";
      if ($ll =~ /(aa|nt), >(.+)\.\.\./) {
        $id = $2;
      }
      else {
        die "format error at $ll\n";
      }
      chop($ll);
      if ($ll =~ /\*$/) { $rep = $id; }
      else {
        push(@members, $ll); $no++;
      }
    }
  }
      if ($no) {
        die "format error, no rep before cluster $ll" unless ($rep);
        if (not defined($slave_clstr{$rep})) {
          $slave_clstr{$rep}=[];
        }
        push(@{$slave_clstr{$rep}}, @members);
      }
  close(TC);
}
##########

my $master_rep = "";
my $master_seq = "";
my $rep_no = 0;
open(TMP, $master_clstr) || die "can not open $master_clstr";
while($ll = <TMP>) {
  if ($ll =~ /^>/) {
    if ($master_rep) {
      print $master_seq;
      if (defined( $slave_clstr{$master_rep} )) {
        foreach $i (@{$slave_clstr{$master_rep}}) {
          $i =~ s/^\d+/$rep_no/;
          print $i, "\n";
          $rep_no++;
        }
      }
    }

    $master_rep = "";
    $master_seq = $ll;
    $rep_no     = 0;
  }
  else {
    $master_seq .= $ll;
    $rep_no++;
    chop($ll);
    if ($ll =~ /\*$/) {
      if ($ll =~ /(aa|nt), >(.+)\.\.\./) {
        $master_rep = $2;
      }
      else {
        die "format error $ll";
      }
    }
  }
}
    if ($master_rep) {
      print $master_seq;
      if (defined( $slave_clstr{$master_rep} )) {
        foreach $i (@{$slave_clstr{$master_rep}}) {
          $i =~ s/^\d+/$rep_no/;
          print $i, "\n";
          $rep_no++;
        }
      }
    }
close(TMP);

