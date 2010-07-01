#!/usr/bin/perl

# the order of clusters need to be identical
my ($master_clstr, @clstr) = @ARGV;
my $clstr_file_no = $#clstr+1;

my @fhs = ();
my @div_reps = ();
my @div_seqs = ();
my @div_rep_no = ();
for ($i=0; $i<$clstr_file_no; $i++) {
  $fh = "FH" . $i;
  open($fh, $clstr[$i]) || die "can not open $clstr[$i]";
  $div_reps[$i] = "";
  $div_seqs[$i] = "";
  $div_rep_no[$i] = 0;
}

my $master_rep = "";
my $master_seq = "";
my $rep_no = 0;
open(TMP, $master_clstr) || die "can not open $master_clstr";
while($ll = <TMP>) {
  if ($ll =~ /^>/) {
    if ($master_rep) {
      print $master_seq;
      foreach ($i=0; $i<$clstr_file_no; $i++) {
        $this_no = process_this($i, $master_rep, $rep_no);
        $rep_no += $this_no;
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
      $rep = "";
      if ($ll =~ /(aa|nt), >(.+)\.\.\./) {
        $rep = $2; 
        $master_rep = $rep;
      }
      else {
        die "format error $ll";
      }
    }
  }
}
    if ($master_rep) {
      print $master_seq;
      foreach ($i=0; $i<$clstr_file_no; $i++) {
        $this_no = process_this($i, $master_rep, $rep_no);
        $rep_no += $this_no;
      }
    }
close(TMP);

for ($i=0; $i<$clstr_file_no; $i++) {
  $fh = "FH" . $i;
  close($fh);
}

sub process_this {
  my ($i, $master_rep, $rep_no) = @_;
  my $ll;
  my ($j, $k);
  $fh = "FH" . $i;

  while($ll = <$fh>) {
    if ($ll =~ /^>/) {

      if ($div_reps[$i] eq $master_rep) {

        if ($div_rep_no[$i] > 1) {
          $j = $rep_no;
          my @lls = split(/\n/,$div_seqs[$i]);
          foreach $k (@lls) {
            next if ($k =~ /\*$/);
            $k =~ s/^\d+/$j/;
            print $k, "\n";
            $j++;
          }
        }

        $div_reps[$i] = "";
        $div_seqs[$i] = "";
        my $t1 = $div_rep_no[$i];
        $div_rep_no[$i] = 0;
 
        return ($t1-1);
        #return ($div_rep_no[$i]-1);
      }
      else {
        $div_reps[$i] = "";
        $div_seqs[$i] = "";
        $div_rep_no[$i] = 0;
      }
    }
    else {
      $div_seqs[$i] .= $ll;
      $div_rep_no[$i]++;
      chop($ll);
      if ($ll =~ /\*$/) {
        my $rep = "";
        if ($ll =~ /(aa|nt), >(.+)\.\.\./) {
          $rep = $2; 
          $div_reps[$i] = $rep;
        }
        else {
          die "format error $ll";
        }
      }
    }
  }

      if ($div_reps[$i] eq $master_rep) {
                                                                                
        if ($div_rep_no[$i] > 1) {
          $j = $rep_no;
          my @lls = split(/\n/,$div_seqs[$i]);
          foreach $k (@lls) {
            next if ($k =~ /\*$/);
            $k =~ s/^\d+/$j/;
            print $k, "\n";
            $j++;
          }
        }
                                                                                
        $div_reps[$i] = "";
        $div_seqs[$i] = "";
        my $t1 = $div_rep_no[$i];
        $div_rep_no[$i] = 0;
                                                                                
        return ($t1-1);
        #return ($div_rep_no[$i]-1);
      }



}
