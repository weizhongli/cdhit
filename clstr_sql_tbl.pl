#!/usr/bin/perl

if(@ARGV==0){
  print "Usage:\n\tclstr_sql_tbl.pl clstr_file tbl_file\n";
  exit(1);
}
#usage clstr_sql_tbl.pl clstr_file tbl_file
  # 1 if tbl_file does not exist
  # create a new tbl_file

  # 2 if tbl_file exists
  # add 2 new columns to it

#format prot_id prot_len cid_level1 rep_level1 cid_level2 rep_level2 ... cid_leveln rep_leveln
# if I have db90.clstr db60.clstr then db30.clstr
#
# run clstr_sql_tbl.pl db90.clstr db_sql.txt #create first 4 columns
#     clstr_sql_tbl.pl db60.clstr db_sql.txt #adding next  2 columns
#     clstr_sql_tbl.pl db30.clstr db_sql.txt #adding next  2 columns

my $clstr_file = shift;
my $tbl_file = shift;
my $tbl_file_tmp = shift;
   $tbl_file_tmp = "$tbl_file.$$" unless ($tbl_file_tmp);
my ($i, $j, $k, $ll);



if (not (-e $tbl_file)) {
  open(TMP, $clstr_file) || die;
  open(OUT, "> $tbl_file") || die;
  my $cid = -1;
  while($ll=<TMP>){
    if ($ll =~ /^>/) {
      $cid++;
    }
    else {
      chop($ll);
      if ($ll =~ /\s(\d+)(aa|nt), >(.+)\.\.\./) {
        my $len = $1;
        my $id = $3;
        my $rep = 0;
        if ($ll =~ /\*$/) { $rep = 1; }
        print OUT "$id\t$len\t$cid\t$rep\n";
      }
      else {
        die "format error $ll";
      }
    }
  }
  close(OUT);
  close(TMP);
}
else {
  add_info_to_tbl();
}

sub add_info_to_tbl {
  my ($i, $j, $k, $ll);

  my %id2cid = ();
  my %idisrep = ();
  open(TMP, $clstr_file) || die;
  my $cid = -1;
  while($ll=<TMP>){
    if ($ll =~ /^>/) {
      $cid++;
    }
    else {
      chop($ll);
      if ($ll =~ /\s(\d+)(aa|nt), >(.+)\.\.\./) {
        my $id = $3;
        if ($ll =~ /\*$/) { $idisrep{$id} = 1;}
        $id2cid{$id} = $cid;
      }
      else {
        die "format error $ll";
      }
    }
  }
  close(TMP);
  print STDERR "done reading $clstr_file\n";

  my @last_cid_2_id = ();
  open(TMP0, $tbl_file) || die;
  while($ll=<TMP0>) {
    chop($ll);
    my @lls = split(/\t/, $ll);
    my $id = $lls[0];
    my $last_level_cid = $lls[-2];
    my $last_level_rep = $lls[-1];
    next unless ($last_level_rep == 1);
    $last_cid_2_id[$last_level_cid] = $id;
  }
  close(TMP0);

  open(TMP2, $tbl_file) || die;
  open(OUT, "> $tbl_file_tmp") || die;
  while($ll=<TMP2>) {
    chop($ll);
    my @lls = split(/\t/, $ll);
    my $id = $lls[0];
    my $last_level_cid = $lls[-2];
    my $last_level_rep = $lls[-1];

    my $this_level_cid = "-1";
    my $this_level_rep = 0;

    if (defined($idisrep{$id})) { $this_level_rep = $idisrep{$id} ;}
    if (defined($id2cid{$id})) { $this_level_cid = $id2cid{$id}; }
    else {
      my $rep1 = $last_cid_2_id[$last_level_cid];
      if (defined($id2cid{$rep1})) { $this_level_cid = $id2cid{$rep1}; }
      else { die "at $ll";}
    }

    print OUT "$ll\t$this_level_cid\t$this_level_rep\n";
  }
  close(OUT);
  close(TMP2);

  if(-e $tbl_file_tmp){
    system("cp $tbl_file_tmp $tbl_file");
    system("rm -f $tbl_file_tmp");
  }
}

