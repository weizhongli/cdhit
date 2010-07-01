#!/usr/bin/perl -w

our $pid         = $$;

our $db_in     = "";
our $db_in2    = "";
our $NR_clstr  = 0.3;
our $NR_clstre = -1;
our $db_out    = "";
our $coverage  = 0.0;
our $coverageR = 0.0;
our $use_prof  = 0;
our $g_iden    = 1;
our $des_len   = 30;
our $len_t     = 10;
our $prof_db   = "nr80";
our $prof_para = "-d $prof_db -j 3 -F F -e 0.001 -b 500 -v 500";
our $bl_expect = "0.000001";
our $bl_para   = "-F F -e $bl_expect -b 100000 -v 100000";
our $psi_blast = "blastpgp ";
our $blastp    = "blastall -a 2 -p blastp";
our $formatdb  = "formatdb";
our $reformat_seg = 200000;
our $para      = 0;
our $host_no   = 0;
our @hosts     = ();
our $bl_suffix = "bl";
our $restart_seg = 5000;
our $keep_bl     = 1;
our $print_db    = 1;
our $job         = "";
our $job_file    = "";
our $local_db    = ""; # in case of para, local dbs on each node maybe faste r
our $date        = `date`;
our $restart_in  = "";
our $blast_with  = "ssh"; # deamon or ssh
our $pbs         = 0;
our $batch_no_per_node = 50; # number of blast jobs per node each submission
our $blast_type  = "blastp";
our $pwd = `pwd`; chop($pwd);
our $db2_len_over = 0;
our $div_no = 16;

our $db_clstr;
our $db_log;
our $db_out1;
our $pl_dir;
our $prof_dir;
our $seq_dir;
our $bl_dir;
our $restart_file;
our $tmp_db;
our $remote_perl_script;



sub parse_para_etc {
  my ($arg, $cmd);
  while($arg = shift) {
    if    ($arg eq "-i")          { $db_in     = shift; }
    elsif ($arg eq "-i2")         { $db_in2    = shift; }
    elsif ($arg eq "-o")          { $db_out    = shift; }
    elsif ($arg eq "-c")          { $NR_clstr  = shift; }
    elsif ($arg eq "-ce")         { $NR_clstre = shift; }
    elsif ($arg eq "-p")          { $prof_para = shift; }
    elsif ($arg eq "-L")          { $coverage  = shift; }
    elsif ($arg eq "-M")          { $coverageR = shift; }
    elsif ($arg eq "-s")          { $bl_para   = shift; }
    elsif ($arg eq "-d")          { $des_len   = shift; }
    elsif ($arg eq "-l")          { $len_t     = shift; }
    elsif ($arg eq "-R")          { $use_prof  = shift; }
    elsif ($arg eq "-G")          { $g_iden    = shift; }
    elsif ($arg eq "-k")          { $keep_bl   = shift; }
    elsif ($arg eq "-b")          { $para      = shift; }
    elsif ($arg eq "-be")         { $bl_expect = shift;
                                    $bl_para =~ s/-e \S+/-e $bl_expect/;}
    elsif ($arg eq "-bfdb")       { $prof_db   = shift;
                                    $prof_para =~ s/-d \S+/-d $prof_db/;}
    elsif ($arg eq "-rs")         { $restart_seg = shift; }
    elsif ($arg eq "-rf")         { $reformat_seg= shift; }
    elsif ($arg eq "-lo")         { $db2_len_over= shift; }
    elsif ($arg eq "-pbs")        { $pbs       = shift; $host_no = $pbs;}
    elsif ($arg eq "-div")        { $div_no    = shift; }
    elsif ($arg eq "-restart")    { $restart_in= shift; }
    elsif ($arg eq "-local")      { $local_db  = shift; }
    elsif ($arg eq "-blastn")     { $blast_type= "blastn"; 
                                    $formatdb  = "formatdb -p F";
                                    $blastp    = "megablast -a 2 -D 2"; }
    elsif ($arg eq "-J") { $job       = shift; $job_file = shift; }
    else {
      ($script_name =~ /cd-hit-2d/) ? print_usage_2d() : print_usage();
      exit();
    }
  }

  # speical jobs
  if ($job eq "parse_blout") { job_parse_blout(); exit();}
  if ($job eq "deamon")      { start_deamon(); exit();}

  (-e $db_in) || die "No input";
  ($db_out)   || die "No output";
  die "pbs or ssh, I'm confused" if ($para and $pbs);

  $db_clstr  = "$db_out.clstr";
  $db_log    = "$db_out.log";
  $db_out1   = "$db_out.out";
  $pl_dir    = "$db_out-perl-$pid";
  $prof_dir  = "$db_in-prof";
  $seq_dir   = "$db_in-seq";
  $bl_dir    = "$db_in-$bl_suffix".(($use_prof) ? 1 : "");
  $restart_file   =" $db_out.restart";

  $tmp_db    = ($db_in2) ? "$db_in2.$pid": "$db_in.$pid";
  $remote_perl_script = "$tmp_db-bl.pl";

  $cmd = `mkdir $prof_dir` if ( $use_prof and (not -e $prof_dir));
  $cmd = `mkdir $pl_dir`   if (($blast_with eq "deamon") and (not -e $pl_dir));
  $cmd = `mkdir $bl_dir`   unless (-e $bl_dir);
  $cmd = `mkdir $seq_dir`  unless (-e $seq_dir);

  read_host() if ($para);
  write_remote_perl_script() if ($host_no);
}
########## END parse_para_etc


sub read_host {
  my $ll;
  open(PARA, "$para") || die "can not open $para";
  while($ll= <PARA>){
    chop($ll); $ll =~ s/\s//g;
    next unless ($ll);
    push(@hosts, $ll); $host_no++;
  }
  close(PARA);
}
########## END read_host


sub read_db {
  my $des = "";
  my $seq = "";
  my $ll;

  open(DBIN, $db_in)         || die "Can not open $db_in";
  while($ll=<DBIN>){
    chop($ll);
    if ($ll =~ /^>/) {
      $seq =~ s/\s//g;
      if (length($seq) > $len_t) { add_seq($des, $seq); }
      $des = $ll; $seq = "";
    }
    else { $seq .= $ll; }
  }
  $seq =~ s/\s//g;
  if (length($seq) > $len_t) { add_seq($des, $seq); }
  close(DBIN);

  ($NR_no >=1 ) || die "No sequence readin";

  print OUTT "Total seqs $NR_no in $db_in\n";
}
########## END read_db


sub read_db2 {
  my $des = "";
  my $seq = "";
  my $ll;
  open(DBIN2, $db_in2)       || die "Can not open $db_in2";
  
  while($ll=<DBIN2>){
    chop($ll);
    if ($ll =~ /^>/) {
      $seq =~ s/\s//g;
      if (length($seq) > $len_t) { add_seq_db2($des, $seq); }
      $des = $ll; $seq = "";
    }
    else { $seq .= $ll; }
  }
  $seq =~ s/\s//g;
  if (length($seq) > $len_t) { add_seq_db2($des, $seq); }
  close(DBIN2);

  ($NR2_no >=1 ) || die "No sequence readin";
  print OUTT "Total seqs $NR2_no in $db_in2\n";
}
########## END read_db2


sub add_seq {
  my ($des, $seq) = @_;
  push(@seqs,   $seq);
  if ($print_db) { push(@dess,   $des); }
  else           { 
    if ($des_len > 0) { push(@dess, substr($des, 0, $des_len)); }
    else {              push(@dess, (split(/\s+/,$des))[0]); }
  }
  push(@lens,   length($seq));
  push(@idens,  0);
  push(@covs,   0);
  push(@passeds,0);
  push(@NR_clstr_nos,0); 
  push(@in_bg, 0) if ($host_no>0);
  $NR_no++; 
}
########## END add_seq


sub add_seq_db2 {
  my ($des, $seq) = @_;
  push(@seqs_db2,   $seq);
  if ($print_db) { push(@dess_db2,   $des); }
  else           { 
    if ($des_len > 0) { push(@dess_db2, substr($des, 0, $des_len)); }
    else {              push(@dess_db2, (split(/\s+/,$des))[0]); }
  }
  push(@lens_db2,   length($seq));
  push(@idens_db2,  0);
  push(@covs_db2,   0);
  push(@passeds_db2,0);
  push(@NR2_clstr_nos,0);
  $NR2_no++;
}
########## END add_seq_db2


sub read_db_no_seq {
  my $des = "";
  my $seq = "";
  my $ll;

  open(DBIN, $db_in)         || die "Can not open $db_in";
  while($ll=<DBIN>){
    chop($ll);
    if ($ll =~ /^>/) {
      $seq =~ s/\s//g;
      if (length($seq) > $len_t) { add_seq_ez($des, $seq); }
      $des = $ll; $seq = "";
    }
    else { $seq .= $ll; }
  }
  $seq =~ s/\s//g;
  if (length($seq) > $len_t) { add_seq_ez($des, $seq); }
  close(DBIN);

  ($NR_no >=1 ) || die "No sequence readin";

  print OUTT "Total seqs $NR_no in $db_in\n";
}
########## END read_db_no_seq



sub add_seq_ez {
  my ($des, $seq) = @_;
  if ($print_db) { push(@dess,   $des); }
  else           { 
    if ($des_len > 0) { push(@dess, substr($des, 0, $des_len)); }
    else {              push(@dess, (split(/\s+/,$des))[0]); }
  }
  push(@lens,   length($seq));
  $NR_no++; 
}
########## END add_seq_ez


sub open_LOG {
  open(OUTT, ">> $db_out1") || die "can not open $db_out1";
  print OUTT "Started $date";

  open(LOG, ">> $db_log")     || die "Can not open $db_log";
}
########## END open_LOG


{## use static variables
my $last_NR90_no=0;
my $last_NR_passed=0;
sub watch_progress {
  my ($i0, $NR90_no, $NR_passed, $NR_no, $flag) = @_;
  my $i1 = $i0+1;

  if ( $i1 % 10 == 0 ) {
    print OUTT ".";
    $flag = 1 if ( $i1 % 100 == 0 );
  }

  if ($flag) {
    my $t1 = (int($NR_passed/$NR_no*10000)) / 100;
    my $t90 = $NR90_no - $last_NR90_no;
    my $tno = $NR_passed - $last_NR_passed; 
    my ($tu, $ts, $cu, $cs) = times();
    my $tt = $tu + $ts + $cu + $cs;
    print OUTT 
      "$i1 finished $NR90_no clusters $NR_passed passed $t90/$tno clstr/passed $t1% done $tt cpu\n";
    $last_NR90_no = $NR90_no;
    $last_NR_passed = $NR_passed;
  }
}
}
########## END add_seq_db2


sub close_LOG {
  my $date = `date`; print OUTT "Completed $date\n";
  if ($host_no) {
    my $total_cpu = total_remote_cpu();
    print OUTT "Total CPUs on remote hosts: $total_cpu\n";
  }
  close(OUTT);
  close(LOG);
}
########## END close_LOG


sub total_remote_cpu {
  my ($i, $j, $k, $ll);
  my $tt = 0;
  for ($j=0; $j<$host_no; $j++) {
    open(TCPU, "$seq_dir/host.$j.cpu") || next;
    while($ll = <TCPU>) {
      chop($ll);
      $tt += $ll;
    }
    close(TCPU);
  }
  return $tt;
}
########## END total_remote_cpu


sub job_parse_blout {
  my ($i, $j, $k);
  my @hits = ($blast_type eq "blastn") ?
             process_blout_blastn($job_file) : process_blout($job_file);


  open(BLOUT2, "> $job_file.out") || return;
  foreach $i (@hits) {
    print BLOUT2 join("\t", @{$i}), "\n";
  }
  print BLOUT2 "#\n";
  close(BLOUT2);
}
########## END job_parse_blout


sub start_deamon_master {
  return unless ($host_no and ($blast_with eq "deamon"));

  my ($i, $j, $k, $cmd);

  my %cpu_of_host = ();
  foreach $i (@hosts) {
    my $host1 = $i;
    if (not defined($cpu_of_host{$host1})) { $cpu_of_host{$host1} = 0; }
    else                                   { $cpu_of_host{$host1}++; }
    my $t1 = "$pl_dir,$host1,$cpu_of_host{$host1}";
    print LOG "run deamon on node $t1\n";
    $cmd = `ssh -xqf $host1 'cd $pwd; $script_name -J deamon $t1 >/dev/null 2>&1 &'`;
  }

}
########## END start_deamon_master


sub start_deamon {
  my ($perl_dir, $node, $cpu) = split(/,/, $job_file);
  my ($i, $j, $k);
  my $signal_f = "$perl_dir/$node.$cpu";
  $i=0;
  while(1) {
    die unless (-e $perl_dir);

    opendir(DIR1, $perl_dir) || next;
    my @files = grep { /\.pl$/} readdir(DIR1);
    foreach $i (@files) {
      if ($i =~ /^$node\.$cpu/) {
        my $cmd = `perl $perl_dir/$i`;
           $cmd = `rm -f $perl_dir/$i`;
      }
    }
    closedir(DIR1);
    sleep(1); $i++;
    if ($i % 30 == 0) {
      my $cmd = `touch $signal_f`;
      $i = 0;
    }
  }
} 
########## END start_deamon


sub write_restart {
  my ($i0, $i, $j, $k);
  open(RES, "> $restart_file") || die;

  for ($i0=0; $i0<$NR_no; $i0++) {
    $i = $NR_idx[$i0];
    print RES "$i\t$NR_clstr_nos[$i]\t$idens[$i]\t$covs[$i]\t$passeds[$i]\n";
  }

  close(RES);
}
########## END write_restart


sub write_restart_db2 {
  my ($i0, $i, $j, $k);
  open(RES, "> $restart_file") || die;

  for ($i0=0; $i0<$NR_no; $i0++) {
    $i = $NR_idx[$i0];
    print RES "$i\t$NR_clstr_nos[$i]\t$idens[$i]\t$covs[$i]\t$passeds[$i]\n";
  }

  #db2 part
  print RES ">>>>\n";
  for ($i0=0; $i0<$NR2_no; $i0++) {
    $i = $i0;
    print RES "$i\t$NR2_clstr_nos[$i]\t$idens_db2[$i]\t$covs_db2[$i]\t$passeds_db2[$i]\n";
  }
  close(RES);
}
########## END write_restart_db2


sub read_restart {
  my ($ii, $i0, $i, $j, $k, $ll);
  my @lls;
  open(RESIN, $restart_in) || die;

  $NR_passed = 0;
  $NR90_no   = 0;
  @NR90_seq  = ();
  $ii = -1;
  $i0 = 0;
  while($ll = <RESIN>) {
    chop($ll);
    @lls = split(/\t/,$ll);
    $i = $lls[0];
    $NR_clstr_nos[$i] = $lls[1];
    $idens[$i]       = $lls[2];
    $covs[$i]        = $lls[3];
    $passeds[$i]     = $lls[4];
    $NR_passed++   if ($lls[4]);

    if ($lls[2] eq "*") { #rep
      $NR90_seq[$lls[1]] = [$i];
      $NR90_no++; 
      $ii = $i0 if ($lls[4]);
    }
    else {
      push( @{$NR90_seq[$lls[1]]}, $i) if ($lls[4]);
    }
    $NR_idx[$i0] = $i;
    $i0++; # idx of sorted , see write_restart 
  }
  close(RESIN);

  $ii++; # $ii to be last rep processed
  return $ii;
}
########## END read_restart


sub read_restart_db2 {
  my ($ii, $i0, $i, $j, $k, $ll);
  my @lls;
  open(RESIN, $restart_in) || die;

  $NR_passed = 0;
  $NR90_no   = 0;
  @NR90_seq  = ();
  $ii = -1;
  $i0 = 0;
  while($ll = <RESIN>) {
    last if ($ll =~ /^>>>>/);
    chop($ll);
    @lls = split(/\t/,$ll);
    $i = $lls[0];
    $NR_clstr_nos[$i] = $lls[1];
    $idens[$i]       = $lls[2];
    $covs[$i]        = $lls[3];
    $passeds[$i]     = $lls[4];
    $NR_passed++   if ($lls[4]);

    if ($lls[2] eq "*") { #rep
      $NR90_seq[$lls[1]] = [$i];
      $NR90_no++; 
      $ii = $i0 if ($lls[4]);
    }
    else {
      push( @{$NR90_seq[$lls[1]]}, $i) if ($lls[4]);
    }
    $NR_idx[$i0] = $i;
    $i0++; # idx of sorted , see write_restart 
  }

  #seqs in db2
  while($ll = <RESIN>) {
    last if ($ll =~ /^>>>>/);
    chop($ll);
    @lls = split(/\t/,$ll);
    $i = $lls[0];
    $NR2_clstr_nos[$i] = $lls[1];
    $idens_db2[$i]       = $lls[2];
    $covs_db2[$i]        = $lls[3];
    $passeds_db2[$i]     = $lls[4];
    $NR2_passed++   if ($lls[4]);
    push( @{$NR90_seq[$lls[1]]}, $i) if ($lls[4]);
  }
  close(RESIN);

  $ii++; # $ii to be last rep processed
  return $ii;
}# END read_restart_db2


sub write_db_clstr {
  my ($i, $j, $k);

  open(DBCLS, "> $db_clstr") || die "Can not write $db_clstr";
  for ($i=0; $i<$NR90_no; $i++) {
    print DBCLS ">Cluster $i\n";
    $k = 0;
    foreach $j (@{ $NR90_seq[$i] }) {
      my $des = ($des_len>0) ? substr($dess[$j], 0, $des_len) : $dess[$j];
      print DBCLS "$k\t$lens[$j]"."aa, $des... ";
      if ($idens[$j] eq "*") { print DBCLS "*\n"; }
      else                   { print DBCLS "at $idens[$j]/$covs[$j]\n";}
      $k++;
    }
  }
  close(DBCLS);
}
########## END write_db_clstr


sub write_db_clstr_db2 {
  my ($i, $j, $k);

  open(DBCLS, "> $db_clstr") || die "Can not write $db_clstr";
  for ($i=0; $i<$NR90_no; $i++) {
    print DBCLS ">Cluster $i\n";
    $k = 0;

    foreach $j (@{ $NR90_seq[$i] }) {
      if ( $k == 0 ) { # for reps in db1
        my $des = ($des_len>0) ? substr($dess[$j], 0, $des_len) : $dess[$j];
        print DBCLS "$k\t$lens[$j]"."aa, $des... ";
        print DBCLS "*\n";
        $k++;
      }
      else {# for redundant seq in db2
        my $des = ($des_len>0) ? substr($dess_db2[$j], 0, $des_len) : $dess_db2[$j];
        print DBCLS "$k\t$lens_db2[$j]"."aa, $des... ";
        print DBCLS "at $idens_db2[$j]/$covs_db2[$j]\n";
        $k++;
      }
    }
  }
  close(DBCLS);
}
########## END write_db_clstr_db2


sub write_db_clstr_db3 {
  my ($i0, $i, $j, $k);

  open(DBCLS, "> $db_clstr") || die "Can not write $db_clstr";

  for ($i0=0; $i0<$NR_no; $i0++) {
    $i = $NR_idx[$i0];
    print DBCLS ">Cluster $i0\n";
    my $des = substr($dess[$i],1);
    $k = 0;
    print DBCLS "$k\t$lens[$i]"."aa, >$des... ";
    print DBCLS "*\n";
    $k++;

    foreach $j ( @{$redundant_in_db2{$des}} ) {
      print DBCLS "$k\t$j\n";
      $k++;
    }

  }
  close(DBCLS);
}
########## END write_db_clstr_db3


sub remove_raw_blout {
  my $NR_sofar = shift;
  my ($i0, $i, $j, $k);
  return if ($keep_bl);

  for ($i0=$NR_sofar; $i0>=0; $i0--) {
    $i = $NR_idx[$i0];
    next unless $passeds[$i];
    next unless ($idens[$i] eq "*"); #only reps have blout
    my $fout = "$bl_dir/$i";
    last unless (-e $fout);          #removed from last call
    my $cmd = `rm -f $fout`;
       $cmd = `rm -f $bl_dir/$i.out` if ($host_no>0);
  }
}
########## END remove_raw_blout


sub fish_other_homolog {
  my ($i, $j, $k, $i0, $j0, $k0);
  $id = shift; # real idx, not sorted idx
  my @hits = ();

  if ($host_no>0) { 
    wait_blast_out("$bl_dir/$id.out");
    open(BLPOUT, "$bl_dir/$id.out") || return;
    while($i=<BLPOUT>) {
      last if ($i =~ /^#/);
      chop($i);
      push(@hits, [split(/\t/,$i)]);
    }
    close(BLPOUT);
  }
  else { 
    run_this_blast($id); 
    return unless (-e "$bl_dir/$id");
    @hits = process_blout("$bl_dir/$id");
  }

  foreach $i (@hits) {
    my $id1 = $i->[0];
    next unless ($id1 < $NR_no); 
    next if ( $passeds[$id1] );
    my $coverage1 = $i->[3];

    next unless (($coverage/$lens[$id]) >= $coverageR);

    $idens[$id1]   = ($NR_clstre==-1) ? $i->[1]. "%" : "$i->[4]/$i->[1]%";
    $covs[$id1]    = $i->[2]. "%";
    $passeds[$id1] = 1;
    $NR_clstr_nos[$id1] = $NR90_no;
    push( @{$NR90_seq[$NR90_no]}, $id1);
    $NR_passed++;
  }
}
########## END fish_other_homolog


sub fish_other_homolog_db2 {
  my ($i, $j, $k, $i0, $j0, $k0);
  $id = shift; # real idx, not sorted idx
  my @hits = ();

  if ($host_no>0) { 
    wait_blast_out("$bl_dir/$id.out");
    open(BLPOUT, "$bl_dir/$id.out") || return;
    while($i=<BLPOUT>) {
      last if ($i =~ /^#/);
      chop($i);
      push(@hits, [split(/\t/,$i)]);
    }
    close(BLPOUT);
  }
  else { 
    run_this_blast($id); 
    return unless (-e "$bl_dir/$id");
    @hits = ($blast_type eq "blastn") ? 
            process_blout_blastn("$bl_dir/$id") : process_blout("$bl_dir/$id");
  }

  # if seq in db 2 is longer then $len_to_skip
  # it won't be grouped into db1 cluster
  # therefore be treated as singletons
  my $len_to_skip = $lens[$id] + $db2_len_over;
  foreach $i (@hits) {
    my $id1 = $i->[0];
    next unless ($id1 < $NR2_no);
    next if ( $passeds_db2[$id1] );
    next if ( ($db2_len_over >= 0) and ($lens_db2[$id1]>$len_to_skip));
    my $coverage1 = $i->[3];

    my $len_longer = ($lens[$id] > $lens_db2[$id1] ) ? 
                      $lens[$id] : $lens_db2[$id1];
    next unless (($coverage/$len_longer) >= $coverageR);

    if ($blast_type eq "blastn") {
      $idens_db2[$id1]   = ($NR_clstre==-1) ? "$i->[5]/$i->[1]%" : 
                                      "$i->[4]/$i->[5]/$i->[1]%";
    }
    else {
      $idens_db2[$id1]   = ($NR_clstre==-1) ? "$i->[1]%" : 
                                      "$i->[4]/$i->[1]%";
    }

    $covs_db2[$id1]    = $i->[2]. "%";
    $passeds_db2[$id1] = 1;
    $NR2_clstr_nos[$id1] = $NR90_no;
    push( @{$NR90_seq[$NR90_no]}, $id1);
    $NR2_passed++;
  }
}
########## END fish_other_homolog_db2

sub process_multi_bl {
  my ($i, $j, $k, $i0, $j0, $k0, $ll);
  my $blout = shift;

  open(TMPBL, $blout) || die;
  while($ll = <TMPBL>) {
    if ($ll =~ /^Query=\s+(\S+)/) {
      my $query_id = $1;
      my $query_len = 0;
      while($ll = <TMPBL>) {
        if ($ll =~ /^\s+.(\d+)\s+letters/) {
          $query_len = $1; last;
        }
      }

      process_multi_bl_each_g("TMPBL", $query_id, $query_len);

    }
  }
  close(TMPBL);
}
########## END process_multi_bl


sub process_multi_bl_each {
  my ($fh1, $query_id, $query_len) = @_;
  my ($i, $j, $k, $i0, $j0, $k0);

  my $bl = readblast_fh("", $fh1, $blast_type);

  for ($i=0; $i<$bl->{no}; $i++) {
    my $p = $bl->{sbj}->[$i];
    my $id1 = (split(/\s+/, $p->{sbj} ))[0]; $id1 = substr($id1,1);

    my $qaln  = $p->{qaln};
    my $saln  = $p->{saln};
    my $alnln = length($qaln);
    my $iden  = $p->{iden};
    if ($g_iden) {
      my $iden_no = 0;
      for ($j=0; $j<$alnln; $j++) {
        my $c1 = substr($qaln,$j,1);
        my $c2 = substr($saln,$j,1);
        $iden_no++ if ($c1 eq $c2);
      }
      $iden = int( 100 * $iden_no / $query_len );
    }
    my $coverage1 = $p->{qend} - $p->{qfrom} + 1;
    my $exp1 = $p->{expect};
    my $len = $query_len;

    if (($iden > ($NR_clstr*100) or $exp1<$NR_clstre) and 
        ($coverage1/$len     >= $coverage)      and
        ($coverage1/$p->{ln} >= $coverageR)) {

      my $lla = "$len"."aa, >$query_id... at";
      my $overlap = "$p->{qfrom}:$p->{qend}:$p->{sfrom}:$p->{send}";
      my $iden2 = ($NR_clstre==-1) ? "$iden%" : "$exp1/$iden%";
      my $cov2  = int($coverage1/$len * 100) . "%";
      push(@{ $redundant_in_db2{$id1} }, "$lla $overlap/$iden2/$cov2");
      last; #just top one hit
    }
  }
}
########## END process_multi_bl_each


#copied from above, but merge hsp to longer one
sub process_multi_bl_each_g {
  my ($fh1, $query_id, $query_len) = @_;
  my ($i, $j, $k, $i0, $j0, $k0);

  my $bl = readblast_fh("", $fh1, $blast_type);

  my $last_id = "";
  my ($exp1, $qfrom, $qend, $sfrom, $send, $hit_ln, $iden, $iden_no, $coverage1);
  my $readin = 0;

  for ($i=0; $i<$bl->{no}; $i++) {
    my $p = $bl->{sbj}->[$i];
    my $id1 = (split(/\s+/, $p->{sbj} ))[0]; $id1 = substr($id1,1);

    if ($id1 ne $last_id) {
      if ($readin) { #process last record

        $iden = int(100 * $iden_no / $query_len ) if ($g_iden);
        my $len = $query_len;
        if (($iden > ($NR_clstr*100) or $exp1<$NR_clstre) and 
            ($coverage1/$len     >= $coverage)      and
            ($coverage1/$hit_ln  >= $coverageR)) {

          my $lla = "$len"."aa, >$query_id... at";
          my $overlap = "$qfrom:$qend:$sfrom:$send";
          my $iden2 = ($NR_clstre==-1) ? "$iden%" : "$exp1/$iden%";
          my $cov2  = int($coverage1/$len * 100) . "%";
          push(@{ $redundant_in_db2{$last_id} }, "$lla $overlap/$iden2/$cov2");
          #last; #just top one hit
          return;
        }
      }

      # readin the first HSP of a new hit
      $readin = 1;
      $exp1   = $p->{expect};
      $hit_ln = $p->{ln};
      $qfrom  = $p->{qfrom};
      $qend   = $p->{qend};
      $sfrom  = $p->{sfrom};
      $send   = $p->{send};
      $iden_no= 0;
      $iden   = $p->{iden};
      $coverage1 = $qend - $qfrom + 1;
      
      if ($g_iden) {
        my $qaln  = $p->{qaln};
        my $saln  = $p->{saln};
        my $alnln = length($qaln);
        for ($j=0; $j<$alnln; $j++) {
          my $c1 = substr($qaln,$j,1);
          my $c2 = substr($saln,$j,1);
          $iden_no++ if ($c1 eq $c2);
        }
      }
      $last_id = $id1;
    }
    else { # same id as last HSP
      # readin following HSPs of same id
      $qfrom1  = $p->{qfrom};
      $qend1   = $p->{qend};
      $sfrom1  = $p->{sfrom};
      $send1   = $p->{send};
      if ( ($qfrom1 > $qend ) and ($sfrom1 > $send) ) { #this HSP after last HSP
         $qend = $qend1; $send = $send1;
      }
      elsif ( ($qend1 < $qfrom) and ($send1 < $sfrom) ) { #this HSP before last HSP
         $qfrom = $qfrom1; $sfrom = $sfrom1;
      }
      else {
        next; # only consider non-overlap, non-crossing HSPs
      }
      $coverage1 += $p->{qend} - $p->{qfrom} + 1;
      
      if ($g_iden) {
        my $qaln  = $p->{qaln};
        my $saln  = $p->{saln};
        my $alnln = length($qaln);
        for ($j=0; $j<$alnln; $j++) {
          my $c1 = substr($qaln,$j,1);
          my $c2 = substr($saln,$j,1);
          $iden_no++ if ($c1 eq $c2);
        }
      }
      $last_id = $id1;
    }
  } #for

      #last record
      if ($readin) { #process last record

        $iden = int(100 * $iden_no / $query_len ) if ($g_iden);
        my $len = $query_len;
        if (($iden > ($NR_clstr*100) or $exp1<$NR_clstre) and 
            ($coverage1/$len     >= $coverage)      and
            ($coverage1/$hit_ln  >= $coverageR)) {

          my $lla = "$len"."aa, >$query_id... at";
          my $overlap = "$qfrom:$qend:$sfrom:$send";
          my $iden2 = ($NR_clstre==-1) ? "$iden%" : "$exp1/$iden%";
          my $cov2  = int($coverage1/$len * 100) . "%";
          push(@{ $redundant_in_db2{$last_id} }, "$lla $overlap/$iden2/$cov2");
        }
      }
}
########## END process_multi_bl_each


sub process_blout {
  my ($i, $j, $k, $i0, $j0, $k0);
  my $blout = shift;
  my @blhits = ();

  my $bl = readblast("", $blout, $blast_type);
  my %hits_i = ();
  my %hits_a = ();
  my $hit_no = 0;
  my $len_hit = ();
  my %hits_2_exp = ();
  for ($i=0; $i<$bl->{no}; $i++) {
    my $p = $bl->{sbj}->[$i];
    my $id1 = (split(/\s+/, $p->{sbj} ))[0]; $id1 = substr($id1,1);

    if (not defined( $hits_i{$id1} )) {
      $hits_i{$id1} = [];
      $hits_a{$id1} = [];
      for ($j=0; $j<$p->{ln}; $j++) {
        push(@{$hits_i{$id1}}, 0); 
        push(@{$hits_a{$id1}}, 0);
      }
      $len_hit{$id1} = $p->{ln};
      $hit_no++;
    }

    if (not defined(   $hits_2_exp{$id1})) { $hits_2_exp{$id1} = 99999; }
    if ($p->{expect} < $hits_2_exp{$id1}) { $hits_2_exp{$id1} = $p->{expect}; }

    my $qaln  = $p->{qaln};
    my $saln  = $p->{saln};
    my $alnln = length($qaln);
    my $si    = $p->{sfrom}-1;
    for ($j=0; $j<$alnln; $j++) {
      my $c1 = substr($qaln,$j,1);
      my $c2 = substr($saln,$j,1);
      if (not $hits_a{$id1}->[$si] ) { # if not covered by pervious HSP
        if ($c1 eq $c2) { $hits_i{$id1}->[$si] = 1; } ## iden aa
        if ($c2 ne '-') { $hits_a{$id1}->[$si] = 1; } ## aligned aa
      }
      $si++ if ($c2 ne '-');
    }
  }

  my $id1;
  foreach $id1 (keys %hits_i) {
    my $len = $len_hit{$id1};

    my $iden = 0;
    my $coverage1 = 0;
    for ($j=0; $j<$len; $j++) {
      $coverage1++ if ( $hits_a{$id1}->[$j] );
      $iden++      if ( $hits_i{$id1}->[$j] );
    }
    $iden      = ($g_iden) ? ($iden/$len) : ($iden/$coverage1);
    #passed
    my $exp1 = $hits_2_exp{$id1};
    if (($iden > $NR_clstr or $exp1<$NR_clstre) 
      and ($coverage1/$len >= $coverage) ) {
      ($NR_clstre == -1) ? 
      push(@blhits, [$id1, int($iden * 100), int($coverage1/$len * 100), $coverage1]) :
      push(@blhits, [$id1, int($iden * 100), int($coverage1/$len * 100), $coverage1, $exp1]);
    }
  }

  ## longest identical segment threshold use following 
  if (0) {
    my $id1;
    foreach $id1 (keys %hits_i) {
      my $len = $len_hit{$id1};
   
      my $longest_iden = 0;
      for ($j=0, $k=0; $j<$len; $j++) {
        if ( $hits_i{$id1}->[$j] ) {
          $k++; $longest_iden = $k if ($k > $longest_iden);
        }
        else {
          $k=0;
        }
      }
   
      #passed
      if ( $longest_iden >= $NR_clstr ) {
        push(@blhits, [$id1, $longest_ide, 0, 0]);
      }
    }
  }
  ## end longest identical segment threshold

  return @blhits;
}
########## END process_blout


sub process_blout_blastn {
  my ($i, $j, $k, $i0, $j0, $k0);
  my $blout = shift;
  my @blhits = ();

  my $bl = readblast("", $blout, $blast_type);
  my %hits_i = ();
  my %hits_a = ();
  my $hit_no = 0;
  my $len_hit = ();
  my %hits_2_exp = ();
  for ($i=0; $i<$bl->{no}; $i++) {
    my $p = $bl->{sbj}->[$i];
    my $id1 = (split(/\s+/, $p->{sbj} ))[0]; $id1 = substr($id1,1);

    my $frame = $p->{frame}; $frame =~ s/\s//g; 
    $id1 .= ($frame =~ /Minus/i) ? "-" : "+";

    if (not defined( $hits_i{$id1} )) {
      $hits_i{$id1} = [];
      $hits_a{$id1} = [];
      for ($j=0; $j<$p->{ln}; $j++) {
        push(@{$hits_i{$id1}}, 0); 
        push(@{$hits_a{$id1}}, 0);
      }
      $len_hit{$id1} = $p->{ln};
      $hit_no++;
    }

    if (not defined(   $hits_2_exp{$id1})) { $hits_2_exp{$id1} = 99999; }
    if ($p->{expect} < $hits_2_exp{$id1}) { $hits_2_exp{$id1} = $p->{expect}; }

    my $qaln  = $p->{qaln};
    my $saln  = $p->{saln};
    my $alnln = length($qaln);
    my $si    = $p->{sfrom}-1;
    for ($j=0; $j<$alnln; $j++) {
      my $c1 = substr($qaln,$j,1);
      my $c2 = substr($saln,$j,1);
      if (not $hits_a{$id1}->[$si] ) { # if not covered by pervious HSP
        if ($c1 eq $c2) { $hits_i{$id1}->[$si] = 1; } ## iden aa
        if ($c2 ne '-') { $hits_a{$id1}->[$si] = 1; } ## aligned aa
      }
      $si++ if ($c2 ne '-');
    }
  }

  my $id1;
  my %id2_added = ();
  foreach $id1 (keys %hits_i) {
    my $len = $len_hit{$id1};
    my $id2 = $id1;
    my $strand = chop($id2);

    next if ($id2_added{$id2}); # if hit added in a strand already

    my $iden = 0;
    my $coverage1 = 0;
    for ($j=0; $j<$len; $j++) {
      $coverage1++ if ( $hits_a{$id1}->[$j] );
      $iden++      if ( $hits_i{$id1}->[$j] );
    }
    $iden      = ($g_iden) ? ($iden/$len) : ($iden/$coverage1);
    #passed
    my $exp1 = $hits_2_exp{$id1};

    if (($iden > $NR_clstr or $exp1<$NR_clstre) 
      and ($coverage1/$len >= $coverage) ) {
      $id2_added{$id2} = 1;
      ($NR_clstre == -1) ? 
      push(@blhits, [$id2, int($iden * 100), int($coverage1/$len * 100), $coverage1, $strand]) :
      push(@blhits, [$id2, int($iden * 100), int($coverage1/$len * 100), $coverage1, $exp1, $strand]);
    }
  }

  return @blhits;
}
########## END process_blout_blastn


## from xblast/blast-parse.pl
## -rwxr-xr-x    1 wli      nfsusers    21077 Apr  9 14:12 blast-parse.pl
##                                                   year 2003
#call readblast($seq, $blout_file, $bl_prog);
#$bl_prog BLAST program:
#                       blastpgp
#                       blastp
#                       tblastn
#                       blastn
#                       blastx
#                       tblastx
sub readblast {
  my ($i, $j, $k, $no);
  my ($q_seq, $filename, $bl_prog) = @_;
  $bl_prog = "blastp" unless ( defined($bl_prog) and $bl_prog);

  my $self = {
    'seq'    => $q_seq,
    'seqln'  => length($q_seq),
    'no'  => 0,
    'round' => 0,
    'sbj' => [],
  };

  ### define some regular expression
  my ($re_gi, $re_ln, $re_score, $re_query, $re_end);
  $re_gi = qr/^>/;
  $re_ln = qr/^\s+Length\s+=\s+\d+\s*$/;
  $re_score = qr/^\s+Score\s+=\s+/;
  $re_query = qr/^Query:\s*\d+/;
  $re_end = qr/^\s+Database:/;


  my $fh  = "BL" ;
  if ( $filename =~ /gz$/ ) { open($fh,"gunzip < $filename |") || return; }
  else {                      open($fh, "$filename") || return; }

  my $round = 0;
  my ($line, $ll, $oline);

  # for blastpgp , we may need found round in case of mulitply iteration
  if ( $bl_prog eq 'blastpgp' or $bl_prog eq 'pdbblast') {

    while( $line = <$fh> ) {
      $line =~ s/\n//;
      if ( $line =~ /^Results from round\s+(\d+)\s*$/ ) { $round = $1; }
    }
    close($fh);

    if ( $filename =~ /gz$/ ) { open($fh,"gunzip < $filename |") || return; }
    else {                      open($fh, "$filename") || return; }

    ### assign self->round
    $self->{round} = $round if ($round > 0);

    ### rewind filehandle or go to the point of last round
    if ( $round > 0) {
      while( $line = <$fh> ) {
        $line =~ s/\n//;
        if ( $line =~ /^Results from round (\d+)\s*$/ ) {
          last if ( $round == $1 );
        }
      }
    }
  }
  # END if ( $bl_prog eq 'blastpgp' or $bl_prog eq 'pdbblast' )


  ### get the first line contain like ">gi|1234567..."
  $oline = "";
  while( $line = ( $oline || <$fh> ) ) {
    $line =~ s/\n//;
    return $self if ( $line =~ /$re_end/); # no hits
    last if ( $line =~ /$re_gi/ );
  }
  $oline = $line;

  ### begin parse the blast output
  @this_sbj = ();
  $no = 0;
  while( $line = ( $oline || <$fh> ) ) {
    $line =~ s/\n//;
    if ( $line =~ /$re_gi/ ) {
      $oline = $line;
      $ll = "";
      while( $line = ( $oline || <$fh> ) ) {
        $line =~ s/\n//;
        if ( $line =~ /$re_ln/) { last; }
        $ll .= $line;
        $oline = "";
      }
      my @gi = ();
      my @pdb = ();
      my @pdbchn = ();
      my $this_id = "";
      my $this_ln = 0;
      my $this_sbj = $ll;
      my $this_cmpd = "";
      if ( $this_sbj =~ /\S+ (.+)/ ) {
        $this_cmpd = $1;
        $this_cmpd=~ s/>.+$//; $this_cmpd =~ s/^\s+//;
      }

      my @tmp_gi_line = split(/gi\|/, $ll);
      for ($i=1; $i<= $#tmp_gi_line; $i++) {
        push(@gi, (split(/\|| /,$tmp_gi_line[$i])) [0]);
      }

      my @tmp_pdb_line = split(/pdb\|/, $ll);
      for ($i=1; $i<= $#tmp_pdb_line; $i++) {
        my $pdbid = substr($tmp_pdb_line[$i],0,6);
        if ( substr($pdbid,5,1) eq " " ) {
          $pdbid = lc(substr($pdbid,0,4)). "_";
        }
        else { $pdbid = lc(substr($pdbid,0,4)) . uc(substr($pdbid,5,1)); }
        push(@pdbchn, $pdbid);
        push(@pdb, substr($pdbid,0,4));
      }

      $this_id = substr($ll,1);
      if ( $this_id =~ /^([\w|]+\|[\d\w]+)[ |]/ ) { $this_id = $1; }
      else { $this_id = ""; }

      if ( $line =~ /^\s+Length\s+=\s+(\d+)\s*$/ ) { $this_ln = $1; }
      else { $this_ln = 0;}
      ### now the @gi, @pdb, @pdbchn, $this_id, $this_ln got defined

      $oline = $line;
      while( $line = ( $oline || <$fh> ) ) {
        $line =~ s/\n//;
        if ( $line =~ /$re_score/ ) { last; }
        $oline = "";
      }
      ### begin a new hsp entry

      $oline = $line;
      while( $line = ( $oline || <$fh> ) ) {
        $line =~ s/\n//;
        if ( $line =~ /$re_score/) {
          ### or later sbj_obj become a package?
          $this_sbj[$no] = {
            'gi'      => [@gi],
            'sbj'     => $this_sbj,
            'cmpd'    => $this_cmpd,
            'id'      => $this_id,
            'ln'      => $this_ln,
            'score'   => 0,
            'bscore'  => 0,
            'expect'  => 0,
            'iden'    => 0,
            'posi'    => 0,
            'gaps'    => 0,
            'alnln'   => 0,
            'qaln'    => '',
            'qfrom'   => 0,
            'qend'    => 0,
            'saln'    => '',
            'sfrom'   => 0,
            'send'    => 0,
            'prog'    => $bl_prog,
            'frame'   => '',
            'alnseq'  => '',
            'gp_alnseq'  => '',
          };

          my $score = $line;
          if ( $score =~ /Score\s*=\s*([+-]?\d+\.?\d?)/ ) {$score = $1;}
          else { $score = 0; }

          my $bscore = $line;
          if ( $bscore =~ /bits\s*\((.+)\)/ )  {$bscore = $1;}
          else { $bscore = 0; }

          my $expect = $line;
          if    ( $expect =~ /Expect\s*=\s*(.+)\s*$/ )        {$expect=$1;}
          elsif ( $expect =~ /Expect\(\d+\)\s*=\s*(.+)\s*$/ ) {$expect=$1;}
          else                                                {$expect= 0;}
          $expect =~ s/^[eE]/1e/; # change e-100 to 1e-100

          $line = <$fh>; $line =~ s/\n//;
          my ($alnln,$iden,$posi,$gaps) = (0, 0, 0, 0);
          if ( $line=~ /Identities\s*=\s*\d+\/(\d+)\s*\((\d+)%\)/) {
            ($alnln,$iden) = ($1, $2);
          }
          if ( $line =~ /Positives\s*=\s*\d+\/\d+\s*\((\d+)%\)/) {
            $posi = $1;
          }
          if ( $line =~ /Gaps\s*=\s*\d+\/\d+\s*\((\d+)%\)/) {
            $gaps = $1;
          }

          $this_sbj[$no]->{score}  = $score;
          $this_sbj[$no]->{bscore} = $bscore;
          $this_sbj[$no]->{expect} = $expect;
          $this_sbj[$no]->{alnln}  = $alnln;
          $this_sbj[$no]->{iden}   = $iden;
          $this_sbj[$no]->{posi}   = $posi;
          $this_sbj[$no]->{gaps}   = $gaps;

          if ($bl_prog =~ /tblastn|blastx|tblastx|blastn/) {
            $line = <$fh>; $line =~ s/\n//;
            my $fr = 0;
            if ( $line =~ /Frame\s*=\s*(.+)/ ) {
              $fr = $1;
            }
            elsif ( $line =~ /Strand\s*=\s*(.+)/ ) {# like Strand = Plus / Plus
              $fr = $1;
            }
            $this_sbj[$no]->{frame} = $fr;
          }


          $oline = $line;
          while( $line = ( $oline || <$fh> ) ) {
            $line =~ s/\n//;
            if ( $line =~ /$re_query/ ){ last; }
            $oline = "";
          }

          my $qaln  = ""; my $saln = "";
          my $qfrom = ""; my $qend = "";
          my $sfrom = ""; my $send = "";
          $oline = $line;
          while( $line = ( $oline || <$fh> ) ) {
            $line =~ s/\n//;
            if ( $line =~ /$re_query/ ){
               if ( $line =~ /^Query:\s*(\d+)\s*([A-Za-z-*]+)\s*(\d+)\s*$/ ) {
                 $qaln .= $2;
                 $qfrom = $1 if ( $qfrom eq "");
                 $qend = $3;
               }
               $line = <$fh>; $line=<$fh>; $line=~ s/\n//;
               if ( $line =~ /^Sbjct:\s*(\d+)\s*([A-Za-z-*]+)\s*(\d+)\s*$/ ) {
                 $saln .= $2;
                 $sfrom = $1 if ( $sfrom eq "");
                 $send = $3;
               }
            }
            if ( ($line =~ /$re_score/) or ($line =~ /$re_gi/)
                                        or ($line =~ /$re_end/) ) {
              $oline=$line;last;
            }
            else { $oline = ""; }
          }

          $this_sbj[$no]->{qfrom} = $qfrom;
          $this_sbj[$no]->{qend}  = $qend;
          $this_sbj[$no]->{sfrom} = $sfrom;
          $this_sbj[$no]->{send}  = $send;
          $this_sbj[$no]->{qaln}  = $qaln;
          $this_sbj[$no]->{saln}  = $saln;

          $no++;
        }
        ### end if ( $line =~ /$re_score/) {
        if ( $line =~ /$re_score/ ) {
          $oline = $line; next;
        }
        if ( ($line =~ /$re_gi/) or ($line=~ /$re_end/) ) {
          $oline = $line; last;
        }
        else { $oline = "";}
      }
      ### end       while( $line = ( $oline || <$fh> ) ) {
    }
    ### end if ( $line =~ /$re_gi/ ) {
    if ( $line =~ /$re_end/ ) { $oline = $line; last; }
  }
  ### end   while( $line = ( $oline || <$fh> ) ) {

  close($fh);

  $self->{no} = $no;
  $self->{sbj} = [@this_sbj];
  return $self;
}
########## END readblast



# but no multi iteration please
# through file handle 
###### sub imported from LIWZblast.pm
#call readblast($seq, $fh, $bl_prog);
#$bl_prog BLAST program:
#                       blastpgp
#                       blastp
#                       tblastn
#                       blastn
#                       blastx
#                       tblastx
sub readblast_fh {
  my ($i, $j, $k, $no);
  my ($q_seq, $fh, $bl_prog) = @_;
  $bl_prog = "blastp" unless ( defined($bl_prog) and $bl_prog);

  my $self = {
    'seq'    => $q_seq,
    'seqln'  => length($q_seq),
    'no'  => 0,
    'round' => 0,
    'sbj' => [],
  };

  ### define some regular expression 
  my ($re_gi, $re_ln, $re_score, $re_query, $re_end);
  $re_gi = qr/^>/;
  $re_ln = qr/^\s+Length\s+=\s+\d+\s*$/;
  $re_score = qr/^\s+Score\s+=\s+/;
  $re_query = qr/^Query:\s*\d+/;
  $re_end = qr/^\s+Database:/;

  my $round = 0;
  my ($line, $ll, $oline);

  ### get the first line contain like ">gi|1234567..."
  $oline = "";
  while( $line = ( $oline || <$fh> ) ) {
    $line =~ s/\n//;
    return $self if ( $line =~ /$re_end/); # no hits
    last if ( $line =~ /$re_gi/ );
  }
  $oline = $line;




  ### begin parse the blast output
  @this_sbj = ();
  $no = 0;
  while( $line = ( $oline || <$fh> ) ) {
    $line =~ s/\n//;
    if ( $line =~ /$re_gi/ ) {
      $oline = $line;
      $ll = "";
      while( $line = ( $oline || <$fh> ) ) {
        $line =~ s/\n//;   
        if ( $line =~ /$re_ln/) {  
          if ($line =~ /Length = (\d+)/ ) {
            $slen=$1;
          } 
          last;
        }
        $ll .= $line;
        $oline = "";
      }
      my @gi = ();
      my @pdb = ();
      my @pdbchn = ();
      my $this_id = "";
      my $this_ln = 0;
      my $this_sbj = $ll;
      my $this_cmpd = "";
      if ( $this_sbj =~ /\S+ (.+)/ ) { 
        $this_cmpd = $1; 
        $this_cmpd=~ s/>.+$//; $this_cmpd =~ s/^\s+//;
      }

      my @tmp_gi_line = split(/gi\|/, $ll);
      for ($i=1; $i<= $#tmp_gi_line; $i++) {
        push(@gi, (split(/\|| /,$tmp_gi_line[$i])) [0]);
      }

      my @tmp_pdb_line = split(/pdb\|/, $ll);
      for ($i=1; $i<= $#tmp_pdb_line; $i++) {
        my $pdbid = substr($tmp_pdb_line[$i],0,6);
        if ( substr($pdbid,5,1) eq " " ) {
          $pdbid = lc(substr($pdbid,0,4)). "_";
        }
        else { $pdbid = lc(substr($pdbid,0,4)) . uc(substr($pdbid,5,1)); }
        push(@pdbchn, $pdbid);
        push(@pdb, substr($pdbid,0,4));
      }

      $this_id = substr($ll,1);
      if ( $this_id =~ /^([\w|]+\|[\d\w]+)[ |]/ ) { $this_id = $1; }
      else { $this_id = ""; }

      if ( $line =~ /^\s+Length\s+=\s+(\d+)\s*$/ ) { $this_ln = $1; }
      else { $this_ln = 0;}
      ### now the @gi, @pdb, @pdbchn, $this_id, $this_ln got defined

      $oline = $line;
      while( $line = ( $oline || <$fh> ) ) {
        $line =~ s/\n//;
        if ( $line =~ /$re_score/ ) { last; }
        $oline = "";
      }
      ### begin a new hsp entry

      $oline = $line;
      while( $line = ( $oline || <$fh> ) ) {
        $line =~ s/\n//;
        if ( $line =~ /$re_score/) {
          ### or later sbj_obj become a package?
          $this_sbj[$no] = {
            'gi'      => [@gi],
            'sbj'     => $this_sbj,
            'cmpd'    => $this_cmpd,
            'id'      => $this_id,
            'ln'      => $this_ln,
            'score'   => 0,
            'bscore'  => 0,
            'expect'  => 0,
            'iden'    => 0,
            'posi'    => 0,
            'gaps'    => 0,
            'alnln'   => 0,
            'qaln'    => '',
            'qfrom'   => 0,
            'qend'    => 0,
            'saln'    => '',
            'sfrom'   => 0,
            'send'    => 0,
            'prog'    => $bl_prog,
            'frame'   => '',
            'alnseq'  => '',
            'gp_alnseq'  => '',
            'slen'    => $slen,
          };

          my $score = $line;
          if ( $score =~ /Score\s*=\s*([+-]?\d+\.?\d?)/ ) {$score = $1;}
          else { $score = 0; }

          my $bscore = $line;
          if ( $bscore =~ /bits\s*\((.+)\)/ )  {$bscore = $1;}
          else { $bscore = 0; }

          my $expect = $line;
          if    ( $expect =~ /Expect\s*=\s*(.+)\s*$/ )        {$expect=$1;}
          elsif ( $expect =~ /Expect\(\d+\)\s*=\s*(.+)\s*$/ ) {$expect=$1;}
          else                                                {$expect= 0;}
          $expect =~ s/^[eE]/1e/; # change e-100 to 1e-100

          $line = <$fh>; $line =~ s/\n//;
          my ($alnln,$iden,$posi,$gaps) = (0, 0, 0, 0);
          if ( $line=~ /Identities\s*=\s*\d+\/(\d+)\s*\((\d+)%\)/) {
            ($alnln,$iden) = ($1, $2);
          }
          if ( $line =~ /Positives\s*=\s*\d+\/\d+\s*\((\d+)%\)/) {
            $posi = $1;
          }
          if ( $line =~ /Gaps\s*=\s*\d+\/\d+\s*\((\d+)%\)/) {
            $gaps = $1;
          }

          $this_sbj[$no]->{score}  = $score;
          $this_sbj[$no]->{bscore} = $bscore;
	  $this_sbj[$no]->{expect} = $expect;
          $this_sbj[$no]->{alnln}  = $alnln;
	  $this_sbj[$no]->{iden}   = $iden;
          $this_sbj[$no]->{posi}   = $posi;
          $this_sbj[$no]->{gaps}   = $gaps; 

          if ($bl_prog =~ /tblastn|blastx|tblastx|blastn/) {
            $line = <$fh>; $line =~ s/\n//;
            my $fr = 0;
            if ( $line =~ /Frame\s*=\s*(.+)/ ) {
              $fr = $1;
            }
            elsif ( $line =~ /Strand\s*=\s*(.+)/ ) {# like Strand = Plus / Plus
              $fr = $1;
            }
            $this_sbj[$no]->{frame} = $fr;
          }

          $oline = $line;
          while( $line = ( $oline || <$fh> ) ) {
            $line =~ s/\n//;
            if ( $line =~ /$re_query/ ){ last; }
            $oline = "";
          }

          my $qaln  = ""; my $saln = "";
          my $qfrom = ""; my $qend = "";
          my $sfrom = ""; my $send = "";
          $oline = $line;
          while( $line = ( $oline || <$fh> ) ) {
            $line =~ s/\n//;
            if ( $line =~ /$re_query/ ){ 
               if ( $line =~ /^Query:\s*(\d+)\s*([A-Za-z-*]+)\s*(\d+)\s*$/ ) {
                 $qaln .= $2;
                 $qfrom = $1 if ( $qfrom eq "");
                 $qend = $3;
               }
               $line = <$fh>; $line=<$fh>; $line=~ s/\n//;
               if ( $line =~ /^Sbjct:\s*(\d+)\s*([A-Za-z-*]+)\s*(\d+)\s*$/ ) {
                 $saln .= $2;
                 $sfrom = $1 if ( $sfrom eq "");
                 $send = $3;
               }
            }
            if ( ($line =~ /$re_score/) or ($line =~ /$re_gi/)
                                        or ($line =~ /$re_end/) ) { 
              $oline=$line;last; 
            }
            else { $oline = ""; }
          }

          $this_sbj[$no]->{qfrom} = $qfrom;
          $this_sbj[$no]->{qend}  = $qend;
          $this_sbj[$no]->{sfrom} = $sfrom;
          $this_sbj[$no]->{send}  = $send;
          $this_sbj[$no]->{qaln}  = $qaln;
          $this_sbj[$no]->{saln}  = $saln;

          $no++;
        }
        ### end if ( $line =~ /$re_score/) {
        if ( $line =~ /$re_score/ ) {
          $oline = $line; next;
        }
        if ( ($line =~ /$re_gi/) or ($line=~ /$re_end/) ) { 
          $oline = $line; last; 
        }
        else { $oline = "";}
      }
      ### end       while( $line = ( $oline || <$fh> ) ) {
    }
    ### end if ( $line =~ /$re_gi/ ) {
    if ( $line =~ /$re_end/ ) { $oline = $line; last; }
  }
  ### end   while( $line = ( $oline || <$fh> ) ) {
  
  $self->{no} = $no;
  $self->{sbj} = [@this_sbj];
  return $self;
}
########## END readblast_fh


#used by psi-cd-hit-domain, just keep here, not used so far
## if $pq is partly represented by other repsentatives
## by $pq->{hit_me_no}, {hit_me_from} {hit_me_end}
## link the X regin below
## protein pq =========XXXXXXXXXXXX================
## then I should use only two "=" region to blast other seqs
## however, for simple calc, I still use the full length of $pq to blast others
## then I should delete the hits found by "X" region
## act, I also delete hits found by part "X" region
sub mask_bl_hits {
  my ($bl, $pq) = @_;
  my ($i, $j, $k);
  my @new_sbj = ();
  my $new_no = 0;
  for ($i=0; $i<$bl->{no}; $i++) {
    my $p = $bl->{sbj}->[$i];
    my $qfrom = $p->{qfrom};
    my $qend  = $p->{qend};

    my $overlap = 0;
    for ($j=0; $j<$pq->{hit_me_no}; $j++) {
      my $from1 = $pq->{hit_me_from}->[$j];
      my $end1  = $pq->{hit_me_end}->[$j];
      my $max_from = ($qfrom > $from1) ? $qfrom : $from1;
      my $min_end  = ($qend  < $end1)  ? $qend  : $end1;
      $overlap = 1 if ($min_end > $max_from);
    }
    next if ($overlap);
    push(@new_sbj, $bl->{sbj}->[$i]);
    $new_no++;
  }
  $bl->{no} = $new_no;
  $bl->{sbj} = [@new_sbj];
}
########## END mask_bl_hits


#used by psi-cd-hit-domain, just keep here, not used so far
sub remove_overlap_hits {
  my $bl = shift;
  my ($i, $j, $k);
  my @new_sbj = ();
  my $new_no = 0;
  my @froms = ();
  my @ends = ();
  my $old_id = "";

  for ($i=0; $i<$bl->{no}; $i++) {
    my $p = $bl->{sbj}->[$i];
    my $id1 = (split(/\s+/, $p->{sbj} ))[0]; $id1 = substr($id1,1);
    my $qfrom = $p->{qfrom};
    my $qend  = $p->{qend};

    if ($id1 ne $old_id) {
      @froms = (); @ends = ();
      $old_id = $id1;
    }
    else {
      my $overlap = 0;
      for ($j=0; $j<=$#froms; $j++) {
        my $from1 = $froms[$j];
        my $end1  = $ends[$j];
        my $max_from = ($qfrom > $from1) ? $qfrom : $from1;
        my $min_end  = ($qend  < $end1)  ? $qend  : $end1;
        $overlap = 1 if ($min_end > $max_from);
      }
      next if ($overlap);
    }
    push(@froms, $qfrom);
    push(@ends,  $qend);
    push(@new_sbj, $bl->{sbj}->[$i]);
    $new_no++;
  }

  $bl->{no} = $new_no;
  $bl->{sbj} = [@new_sbj];

}
########## END remove_overlap_hits

sub blast_formatdb_raw {
  my $cmd = `formatdb -i $db_in`;
}
########## END blast_formatdb_raw


sub blast_formatdb {
  my ($i0, $i, $j, $k);
  my $which_db = shift;
     $which_db = "db1" unless (defined($which_db));

  open(FDB, "> $tmp_db") || die;
  $j=0;
  if ($which_db eq "db2") { # for psi-cd-hit-2d
    $i0 = shift; # the one currently processing

    # if seq in db 2 is longer than $len_to_skip
    # it won't be grouped into db1 cluster
    # therefore be skipped
    my $len_to_skip = $lens[ $NR_idx[$i0] ] + $db2_len_over;

    for ($i=0; $i<$NR2_no; $i++) {
      next if ($passeds_db2[$i]);
      next if ( ($db2_len_over >= 0) and ($lens_db2[$i]>$len_to_skip));

      my $seq = $seqs_db2[$i];
      $seq =~ s/(.{70})/$1\n/g;
      $seq =~ s/\n$//;
      print FDB ">$i $dess_db2[$i]\n$seq\n";
      $j++;
    }
  }
  else { # for psi-cd-hit
    for ($i=0; $i<$NR_no; $i++) {
      next if ($passeds[$i]);
      my $seq = $seqs[$i];
      $seq =~ s/(.{70})/$1\n/g;
      $seq =~ s/\n$//;
      print FDB ">$i $dess[$i]\n$seq\n";
      $j++;
    }
  }
  close(FDB);
  $formatdb_no = $j;

  while(1) {
    opendir(SEQDB, $seq_dir) || next;
    my @leftseqs = grep {/lock/} readdir(SEQDB);
    closedir(SEQDB);
                                                                                       
    last unless @leftseqs;
    sleep(3);
  }

  my $cmd = `$formatdb -i $tmp_db`;
  ((-e "$tmp_db.phr") and (-e "$tmp_db.pin") and (-e "$tmp_db.psq")) ||
  ((-e "$tmp_db.nhr") and (-e "$tmp_db.nin") and (-e "$tmp_db.nsq"))
     || die "Can not formatdb";

  if ($host_no>0 and $local_db and (not $pbs)) {
    my %host_write = ();
    foreach $i (@hosts) {
      next if ($host_write{$i}); # prevent same host copy db mulitple times
      print LOG "copy $tmp_db to node $i`\n";
      $cmd = `ssh -x $i cp -fp $pwd/$tmp_db.p* $local_db`;
      $host_write{$i} = 1;
    }
    # sleep(120);
    print LOG "\n\n";
  }
}
########## END blast_formatdb


sub remove_blast_db {
  my ($i, $j, $k);
  $cmd = `rm -f $tmp_db`;
  $cmd = `rm -f $tmp_db.p*`;

  if ($host_no>0 and $local_db and (not $pbs)) {
    my %host_write = ();
    foreach $i (@hosts) {
      next if ($host_write{$i}); # prevent same host copy db mulitple times
      print LOG "remove $tmp_db no node $i\n";
      $cmd = `ssh -x $i rm -f $local_db/$tmp_db.p* >/dev/null 2>&1 &`;
      $host_write{$i} = 1;
    }
  }
}
########## END remove_blast_db


my $common_usage = <<EOD;

Options 
         -i  in_dbname, required
         -o  out_dbname, required
         -c  clustering threshold (sequence identity), default 0.3
         -ce clustering threshold (blast expect), default -1, 
             it means by default it doesn't use expect threshold, 
             but with positive value, the program cluster seqs if similarities
             meet either identity threshold or expect threshold 
         -L  coverage of shorter sequence ( aligned / full), default 0.0
         -M  coverage of longer sequence ( aligned / full), default 0.0
         -R  (1/0) use psi-blast profile? default 0
             perform psi-blast / pdb-blast type search
         -G  (1/0) use global identity? default 1
             sequence identity calculated as 
               total identical residues of local alignments / 
               length of shorter seq
             if you prefer to use -G 0, it is suggested that you also
             use -L, such as -L 0.8, to prevent very short matches.
         -d  length of description line in the .clstr file, default 30
             if set to 0, it takes the fasta defline and stops at first space
         -l  length_of_throw_away_sequences, default 10
         -p  profile search para, default
               "-a 2 -d nr80 -j 3 -F F -e 0.001 -b 500 -v 500"
         -bfdb profile database, default nr80
         -s  blast search para, default
               "-F F -e 0.000001 -b 100000 -v 100000"
         -be blast expect cutoff, default 0.000001
         -b  filename of list of hosts
             to run this program in parallel with ssh calls, you need provide
             a list of hosts
         -pbs No of jobs to send each time by PBS querying system
             you can not use both ssh and pbs at same time
         -k (1/0) keep blast raw output file, default 1
         -rs steps of save restart file and clustering output, default 5000
             everytime after process 5000 sequences, program write a 
             restart file and current clustering information
         -restart restart file, readin a restart file
             if program crash, stoped, termitated, you can restart it by
             add a option "-restart sth.restart"
         -rf steps of re format blast database, default 200,000
             if program clustered 200,000 seqs, it remove them from seq
             pool, and re format blast db to save time
         -local dir of local blast db, 
             when run in parallel with ssh (not pbs), I can copy blast dbs
             to local drives on each node to save blast db reading time
             BUT, IT MAY NOT FASTER
         -J  job, job_file, exe specific jobs like parse blast outonly
             DON'T use it, it is only used by this program itself
         -single files of ids those you known that they are singletons
            so I won't run them as queries
EOD


sub print_usage {
  print <<EOD;
Usage psi-cd-hit [Options] 
$common_usage

                                    ==============================
                                    by Weizhong Li, liwz\@sdsc.edu
                                    ==============================
    If you find cd-hit useful, please kindly cite:

    "Clustering of highly homologous sequences to reduce thesize of large protein database", Weizhong Li, Lukasz Jaroszewski & Adam GodzikBioinformatics, (2001) 17:282-283
    "Cd-hit: a fast program for clustering and comparing large sets of protein or nucleotide sequences", Weizhong Li & Adam Godzik Bioinformatics, (2006) 22:1658-1659

EOD
}
########## END print_usage


sub print_usage_2d {
  print <<EOD;
Usage psi-cd-hit-2d [Options] 
$common_usage
         -i2 second input database
         -blastn run blastn, default 0
         -lo how long can seq in db2 > db1 in a cluster, default 0
             means, that seq in db2 should <= seqs in db1 in a cluster

                                    ==============================
                                    by Weizhong Li, liwz\@sdsc.edu
                                    ==============================

    If you find cd-hit useful, please kindly cite:

    "Clustering of highly homologous sequences to reduce thesize of large protein database", Weizhong Li, Lukasz Jaroszewski & Adam GodzikBioinformatics, (2001) 17:282-283
    "Cd-hit: a fast program for clustering and comparing large sets of protein or nucleotide sequences", Weizhong Li & Adam Godzik Bioinformatics, (2006) 22:1658-1659


EOD
}
########## END print_usage_2d


## like above, but don't assign seqs to specific node
## while let nodes run them autoly
sub  run_batch_blast3 {
  my $i0 = shift;
  my ($id, $i, $j, $k);
 
  my $total_jobs = $batch_no_per_node * $host_no;

  $k=0;
  for ($j=0; $i0<$NR_no; $i0++, $j++) {
    $j=0 if ($j>=$host_no);

    $id = $NR_idx[$i0];
    next if ($passeds[$id]);
    next if ($in_bg[$id]);
    next if ($known_single and $known_singles[$id]);
    $in_bg[$id] = 1;

    my $seq = $seqs[$id];
    open(SEQ, "> $seq_dir/$id") || die "Can not write";
    print SEQ "$dess[$id]\n$seq\n";
    close(SEQ);

    $k++; last if ($k >= $total_jobs);
  }

  for ($j=0; $j<$host_no; $j++) {
    if ($pbs) { # using PBS querying system
      my $t = "psi-cd-hit-$j";
      print LOG "PBS querying $j\n";

      open(QUEUE,"| qsub -N $t");
      print QUEUE <<EOD;
#!/bin/sh
#\$ -S /bin/bash
#\$ -v PATH

cd $pwd
./$remote_perl_script $j

EOD
      close(QUEUE);


      #open(QUEUE,"| qsub -N $t -o $t.log -e $t.err");
      #print QUEUE "cd $pwd; ./$remote_perl_script $j";
      #close(QUEUE);
    }
    else {
      my $host1 = $hosts[$j];
      print LOG "running blast on $host1\n";
      my $cmd = `ssh -xqf $host1 'cd $pwd; ./$remote_perl_script $j >/dev/null 2>&1 &'`;
    }
  }
} 
########## END run_batch_blast3


sub run_this_blast {
  my $id = shift;
  my $cmd;

  open(SEQ, "> $seq_dir/$id") || die "Can not write";
  print SEQ "$dess[$id]\n$seqs[$id]\n";
  close(SEQ);

  #calculate profile
  if ((not (-e "$prof_dir/$id.prof")) and $use_prof) {
    $cmd = `$psi_blast $prof_para -i $seq_dir/$id -C $prof_dir/$id.prof`;
  }

  if ( not (-e "$bl_dir/$id") ) {
    if ((-e "$prof_dir/$id.prof") and $use_prof) {
      $cmd = `$psi_blast -d ./$tmp_db $bl_para -i $seq_dir/$id -R $prof_dir/$id.prof -o $bl_dir/$id`;
    }
    else {
      $cmd = `$blastp -d ./$tmp_db $bl_para -i $seq_dir/$id -o $bl_dir/$id`;
    }
  }

  $cmd = `rm -f $seq_dir/$id`;
}
########## END run_this_blast


sub write_remote_perl_script {
  my $dir1 = ($local_db) ? $local_db : ".";
  my $bl2  = ($use_prof) ?
     "$psi_blast -d $dir1/$tmp_db $bl_para -R $prof_dir/\$id.prof":
        "$blastp -d $dir1/$tmp_db $bl_para";

  my $blastn_yes = ($blast_type eq "blastn") ? "-blastn" : "";

  open(REPERL, "> $remote_perl_script") || die;
  print REPERL <<EOD;
#!/usr/bin/perl
\$host = shift;
\$arg = shift;

if (\$arg) {
  \@ids = split(/,/, \$arg);
}
else {
  while(1) {
    if (opendir(DDIR, "$seq_dir")) { 
      \@ids = grep {/^\\d+\$/} readdir(DDIR);
      last;
    }
    else {
      sleep(1);
    }
  }
}

foreach \$id (\@ids) {

  next unless (-e "$seq_dir/\$id");
  next if (-e "$seq_dir/\$id.lock");
  \$cmd = `touch $seq_dir/\$id.lock`;

  if ($use_prof) {
    if (not (-e "$prof_dir/\$id.prof")) {
      \$cmd = `$psi_blast $prof_para -i $seq_dir/\$id -C $prof_dir/\$id.prof`;
    }
  }

  if ( not (-e "$bl_dir/\$id") ) {
    \$cmd = `$bl2 -i $seq_dir/\$id -o $bl_dir/\$id`;
  }

  \$cmd = `$script_name -J parse_blout $bl_dir/\$id -c $NR_clstr -ce $NR_clstre -L $coverage -M $coverageR -G $g_iden $blastn_yes`;

  \$cmd = `rm -f  $seq_dir/\$id`;
  \$cmd = `rm -f  $seq_dir/\$id.lock`;
}

(\$tu, \$ts, \$cu, \$cs) = times();
\$tt = \$tu + \$ts + \$cu + \$cs;
\$cmd = `echo \$tt >> $seq_dir/host.\$host.cpu`;

EOD
  close(REPERL);
  my $cmd = `chmod 755 $remote_perl_script`;

} 
########## END write_remote_perl_script


sub wait_blast_out {
  my $out = shift;
  print LOG "waiting for $out";
  while(1) {
    if (-e $out) {
      my $last = `tail -1 $out`;
      chop($last);
      last if ($last =~ /^#$/);
    }
    sleep(1);
    print LOG ".";
  }
  print LOG "\n";
}
########## END wait_blast_out



1;

