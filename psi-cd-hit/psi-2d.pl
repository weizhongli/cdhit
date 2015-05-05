#!/usr/bin/perl


my $script_name = $0;
my $script_dir = $0;
   $script_dir =~ s/[^\/]+$//;
   chop($script_dir);
   $script_dir = "./" unless ($script_dir);


use Getopt::Std;
getopts("i:j:b:o:c:l:d:q:",\%opts);
die usage() unless ($opts{i} and $opts{j} and $opts{b} and $opts{o});

my ($i, $j, $k, $cmd, $ll);
my $db1        = $opts{i};
my $db2        = $opts{j};
my $bl         = $opts{b}; # $bl is output file that blast db against db2
                           # formatdb -i db2 -p F
                           # blastall -n blastn -i db1 -d db2 -m 8
my $output     = $opts{o};
my $circle     = $opts{c}; 
my $cutoff_len = $opts{l}; $cutoff_len = 0 unless ($cutoff_len);
my $cutoff_idens= $opts{d}; $cutoff_idens= 0 unless ($cutoff_idens);
my $cutoff_idenq= $opts{q}; $cutoff_idenq= 0 unless ($cutoff_idenq);

my @id_db1 = (); my %id_db1_2_len = (); %id_db1_2_len_nN = (); my $db1_no = 0;
my @id_db2 = (); my %id_db2_2_len = (); %id_db2_2_len_nN = (); my $db2_no = 0;
read_db1();
read_db2();

open(OUT, ">$output") || die "can not open $output\n";

open(BL, $bl) || die "can not open $bl\n";
my $last_qid = "";
my @bl_sbj;
my $bl_data_no=0;
while($ll = <BL>){
  next if ($ll =~ /^#/);
  chop($ll);
  my @lls = split(/\t/,$ll);

  if ($lls[0] ne $last_qid) {
    if ($bl_data_no>0) {
      my $bl_data = {
       'no'  => $bl_data_no,
       'sbj' => [@bl_sbj],
      };
      process_this_bl( $bl_data );
    }
    @bl_sbj = ();
    $bl_data_no = 0;
  }
  my $frame = "";
     $frame .= ($lls[6] < $lls[7]) ? "+" : "-";
     $frame .= ($lls[8] < $lls[9]) ? "+" : "-";
  $bl_sbj[$bl_data_no] = {
    'qid'     => $lls[0],
    'id'      => $lls[1],
    'iden'    => $lls[2],
    'alnln'   => $lls[3],
    'ms'      => $lls[4],
    'gap'     => $lls[5],
    'qfrom'   => $lls[6],
    'qend'    => $lls[7],
    'sfrom'   => $lls[8],
    'send'    => $lls[9],
    'expect'  => $lls[10],
    'score'   => $lls[11],
    'frame'   => $frame,
  };
  $bl_data_no++;
  $last_qid = $lls[0];
}
close(BL);
    if ($bl_data_no>0) {
      my $bl_data = {
       'no'  => $bl_data_no,
       'sbj' => [@bl_sbj],
      };
      process_this_bl( $bl_data );
    }

close(OUT);
################################################################################

sub process_this_bl {
  my $bl = shift;
  my ($i, $j, $k, $ll);
  keep_strand_with_top_hsp($bl);
  if ($circle) { reset_alignment_coor_for_circle_seq($bl); }

  filter_short_hits($bl); #### remove hits based on total lengths
                          #### still need to remove based on co-linner hits
  return unless ($bl->{no});

  my $qid = $bl->{sbj}->[0]->{qid};
  my %hit_ids = map { $_->{id}, 1} @{$bl->{sbj}};
  my $hit_no = 0;
  my $text_data = "";

  my @hsp = (); #### [id, len, qfrom, qend, sbegin, send, expect]
  my $hsp_no = 0;
  my $hit_len = 0;

  for ($i=0; $i<$bl->{no}; $i++) {
    my $p = $bl->{sbj}->[$i];
    my $id1 = $p->{id};

    if ($hsp_no) {
      if ($id1 ne $hsp[0]->[0]) {
        if (($hit_len >= $cutoff_len) and 
            ($hit_len >= $cutoff_idens * $id_db2_2_len{$hsp[0]->[0]}) and 
            ($hit_len >= $cutoff_idenq * $id_db1_2_len{$qid})) {
          @hsp = sort {$a->[2] <=> $b->[2]} @hsp;
          $text_data .= print_hsp($qid, @hsp);
          $hit_no++;
        }

        @hsp = ();
        @aln_lens = ();
        $hsp_no = 0;
        $hit_len = 0;
      }
    }

    #check whether overlap with previous high score HSPs
    my $overlap_flag = 0;
    for ($j=0; $j<$hsp_no; $j++) {
      if (overlap1($p->{qfrom}, $p->{qend}, $hsp[$j]->[2], $hsp[$j]->[3])) { $overlap_flag = 1; last; }
      if (overlap1($p->{sfrom}, $p->{send}, $hsp[$j]->[4], $hsp[$j]->[5])) { $overlap_flag = 1; last; }
    }
    next if ($overlap_flag);

    #check whether this HSP cross with previous high score HSPs
    my $cross_flag = 0;
    for ($j=0; $j<$hsp_no; $j++) {
      if (cross1($p->{qfrom}, $p->{qend}, $hsp[$j]->[2], $hsp[$j]->[3],
                 $p->{sfrom}, $p->{send}, $hsp[$j]->[4], $hsp[$j]->[5])) {
        $cross_flag = 1; last;
      }
    }
    next if ($cross_flag);

    push(@hsp, [$id1, $len_sub, $p->{qfrom}, $p->{qend}, $p->{sfrom}, $p->{send}, $p->{expect}, $p->{frame}]);
    $hit_len += abs($p->{sfrom} - $p->{send});
    $hsp_no++;
  }
  if ($hsp_no) {
    if (($hit_len >= $cutoff_len) and 
        ($hit_len >= $cutoff_idens * $id_db2_2_len{$hsp[0]->[0]}) and 
        ($hit_len >= $cutoff_idenq * $id_db1_2_len{$qid})) {
      @hsp = sort {$a->[2] <=> $b->[2]} @hsp;
      $text_data .= print_hsp($qid, @hsp);
      $hit_no++;
    }
  }

  print OUT ">$qid\t$id_db1_2_len{$qid}\t$hit_no\n";
  print OUT $text_data;
}
########## END process_this_bl

sub print_hsp {
  my ($qid, @hsp) = @_;
  my $hsp_no = $#hsp+1;
  my ($i, $j, $k);
  my $print_hsp_return = "";

  my @str1 = (); my $aln1 = 0;
  my @str2 = (); my $aln2 = 0;
  for ($i=0; $i<$hsp_no; $i++){
    push(@str1, "$hsp[$i]->[2]-$hsp[$i]->[3]"); $aln1 += $hsp[$i]->[3]-$hsp[$i]->[2];
    push(@str2, "$hsp[$i]->[4]-$hsp[$i]->[5]"); $aln2 += $hsp[$i]->[5]-$hsp[$i]->[4];
  }
  my $sid = $hsp[0]->[0];

  $print_hsp_return .= "\t$qid\t$id_db1_2_len{$qid}\t$id_db1_2_len_nN{$qid}\t$aln1\t". join("|", @str1). "\n";
  $print_hsp_return .= "\t$sid\t$id_db2_2_len{$sid}\t$id_db1_2_len_nN{$qid}\t$aln2\t". join("|", @str2). "\n";

  return $print_hsp_return;
}
########## END print_hsp


sub overlap1 {
  my ($b1, $e1, $b2, $e2) = @_;

  my $t; ### 
  if ($e1 < $b1) { $t  = $e1; $e1 = $b1; $b1 = $t; }
  if ($e2 < $b2) { $t  = $e2; $e2 = $b2; $b2 = $t; }

  return 0 if ($e2 < $b1);
  return 0 if ($b2 > $e1);
  return ( ($e1<$e2)? $e1:$e2 )-( ($b1>$b2)? $b1:$b2);
}
########## END overlap1


sub cross1 {
  my ($q_b1, $q_e1, $q_b2, $q_e2,
      $s_b1, $s_e1, $s_b2, $s_e2) = @_;

  my $fr_q1 = ($q_b1 < $q_e1) ? 1 : -1;
  my $fr_q2 = ($q_b2 < $q_e2) ? 1 : -1;
  my $fr_s1 = ($s_b1 < $s_e1) ? 1 : -1;
  my $fr_s2 = ($s_b2 < $s_e2) ? 1 : -1;

  my $fr1 = $fr_q1 * $fr_s1;
  my $fr2 = $fr_q2 * $fr_s2;
  return 1 if (($fr1 * $fr2) < 0); # one ++ and one +-


  my $t;
  if ($q_e1 < $q_b1) { $t    = $q_e1; $q_e1 = $q_b1; $q_b1 = $t; }
  if ($q_e2 < $q_b2) { $t    = $q_e2; $q_e2 = $q_b2; $q_b2 = $t; }
  if ($s_e1 < $s_b1) { $t    = $s_e1; $s_e1 = $s_b1; $s_b1 = $t; }
  if ($s_e2 < $s_b2) { $t    = $s_e2; $s_e2 = $s_b2; $s_b2 = $t; }

# after above transformation
#                               0           q_b1           q_e1             q_b2        q_e2       qlen
# query                      5' ====================================================================
# match                                     ||||||||||||||||                |||||||||||||
# subject                         5' ========================================================================>>>>>> frame +
#                                    0      s_b1           s_e1             s_b2        s_e2                 slen

# match                                     ||||||||||||||||                |||||||||||||
# subject                         3' ========================================================================>>>>>> frame -
#                                    slen   s_e1           s_b1             s_e2        s_b2                 0   

  if (($fr1 > 0) and ($fr2>0)) { # both ++
    return ( (($q_b2-$q_b1)*($s_b2-$s_b1) <0) ? 1 : 0);
  }
  else { # both --
    return ( (($q_b2-$q_b1)*($s_e1-$s_e2) <0) ? 1 : 0);
  }

}
########## END cross1


########## let the top hsp to start at 0 for both query and subject
########## i.e. the begining of HSP to be new original - coordinate 0
########## then reset all other HSPs' alignment coordinates
sub reset_alignment_coor_for_circle_seq {
  my $self = shift;
  my ($i,$j,$k);

  my $last_id = "";
  $j = 0;
  my $hsp_count = 0; # number of HSPs for a subject
  for ($i=0; $i<$self->{no}; $i++) {
    my $p = $self->{sbj}->[$i];
    my ($id1, $len_sub) = split(/\./, $p->{id});

    if ($id1 ne $last_id) {
      if ($hsp_count > 1) {  # it is necessary to reset coordinate when at least 2 HSP
        my $p_top_hsp = $self->{sbj}->[$j];
        my $len_q = (split(/\./, $p_top_hsp->{qid}))[1];
        my $len_s = (split(/\./, $p_top_hsp->{id}))[1];
        my $ref_q = ($p_top_hsp->{qfrom} < $p_top_hsp->{qend}) ? $p_top_hsp->{qfrom} : $p_top_hsp->{qend};
        my $ref_s = ($p_top_hsp->{sfrom} < $p_top_hsp->{send}) ? $p_top_hsp->{sfrom} : $p_top_hsp->{send};
        for ($k = $j; $k<$j+$hsp_count; $k++) {
          $self->{sbj}->[$k]->{qfrom} -= $ref_q; if ($self->{sbj}->[$k]->{qfrom} < 0) {$self->{sbj}->[$k]->{qfrom} += $len_q;}
          $self->{sbj}->[$k]->{qend}  -= $ref_q; if ($self->{sbj}->[$k]->{qend}  < 0) {$self->{sbj}->[$k]->{qend}  += $len_q;}
          $self->{sbj}->[$k]->{sfrom} -= $ref_s; if ($self->{sbj}->[$k]->{sfrom} < 0) {$self->{sbj}->[$k]->{sfrom} += $len_s;}
          $self->{sbj}->[$k]->{send}  -= $ref_s; if ($self->{sbj}->[$k]->{send}  < 0) {$self->{sbj}->[$k]->{send}  += $len_s;}
        }
      }
      $j = $i;
      $hsp_count = 0;
    }
    $last_id = $id1;
    $hsp_count++;
  }

      #last subject
      if ($hsp_count > 1) {  # it is necessary to reset coordinate when at least 2 HSP
        my $p_top_hsp = $self->{sbj}->[$j];
        my $len_q = (split(/\./, $p_top_hsp->{qid}))[1];
        my $len_s = (split(/\./, $p_top_hsp->{id}))[1];
        my $ref_q = ($p_top_hsp->{qfrom} < $p_top_hsp->{qend}) ? $p_top_hsp->{qfrom} : $p_top_hsp->{qend};
        my $ref_s = ($p_top_hsp->{sfrom} < $p_top_hsp->{send}) ? $p_top_hsp->{sfrom} : $p_top_hsp->{send};
        for ($k = $j; $k<$j+$hsp_count; $k++) {
          $self->{sbj}->[$k]->{qfrom} -= $ref_q; if ($self->{sbj}->[$k]->{qfrom} < 0) {$self->{sbj}->[$k]->{qfrom} += $len_q;}
          $self->{sbj}->[$k]->{qend}  -= $ref_q; if ($self->{sbj}->[$k]->{qend}  < 0) {$self->{sbj}->[$k]->{qend}  += $len_q;}
          $self->{sbj}->[$k]->{sfrom} -= $ref_s; if ($self->{sbj}->[$k]->{sfrom} < 0) {$self->{sbj}->[$k]->{sfrom} += $len_s;}
          $self->{sbj}->[$k]->{send}  -= $ref_s; if ($self->{sbj}->[$k]->{send}  < 0) {$self->{sbj}->[$k]->{send}  += $len_s;}
        }
      }

  return;
}
########## reset_alignment_coor_for_circle_seq

########### if a hit has multiple HSPs on both + - strands
########### keep only the HSPs, whose strand is same as the top HSP
sub keep_strand_with_top_hsp {
  my $self = shift;
  my ($i,$j,$k);

  my %id_2_strand = ();
  my @new_sbj = ();
  my $new_no = 0;
  for ($i=0; $i<$self->{no}; $i++) {
    my $p = $self->{sbj}->[$i];
    my $id1 = $p->{id};
    if (not defined($id_2_strand{$id1})) {
      $id_2_strand{$id1} = $p->{frame};
    }
    if ($p->{frame} eq $id_2_strand{$id1}) { #### this stand is same as the top strand
      push(@new_sbj, $self->{sbj}->[$i]);
      $new_no++;
    }
  }
  $self->{no} = $new_no;
  $self->{sbj} = [@new_sbj];
}
########## END keep_strand_with_top_hsp

########## remove short hits
sub filter_short_hits {
  my $self = shift;
  my ($i,$j,$k);

  my $qid = $self->{sbj}->[0]->{qid};

  my %hit_len = ();
  for ($i=0; $i<$self->{no}; $i++) {
    my $p = $self->{sbj}->[$i];
    my $id1 = $p->{id};
    $hit_len{$id1} += abs($p->{sfrom} - $p->{send});
  }

  my @new_sbj = ();
  my $new_no = 0;
  for ($i=0; $i<$self->{no}; $i++) {
    my $p = $self->{sbj}->[$i];
    my $id1 = $p->{id};

    if (($hit_len{$id1} >= $cutoff_len) and 
        ($hit_len{$id1} >= $cutoff_idens * $id_db2_2_len{$id1}) and
        ($hit_len{$id1} >= $cutoff_idenq * $id_db1_2_len{$qid})) {
      push(@new_sbj, $self->{sbj}->[$i]);
      $new_no++;
    }
  }
  $self->{no} = $new_no;
  $self->{sbj} = [@new_sbj];
}
########## END filter_short_hits 


sub read_db1 {
  my ($i, $j, $k, $ll, $len, $id, $len2);
  open(TMP, $db1) || die "Can not open $db1\n";

  $len = 0;
  $len2 = 0;
  while($ll=<TMP>){
    if ($ll =~ /^>/) {
      if ($len > 0) {
        $id_db1_2_len{$id} = $len;
        $id_db1_2_len_nN{$id} = $len2;
      }
      chop($ll);
      $ll =~ s/\s.+$//;
      $id = substr($ll,1);
      $len = 0;
      $len2 = 0;
    }
    else {
      chop($ll);
      $ll =~ s/\s//;
      $len += length($ll);
      $ll =~ s/N//ig;
      $len2 += length($ll);
    }
  }
      if ($len > 0) {
        $id_db1_2_len{$id} = $len;
        $id_db1_2_len_nN{$id} = $len2;
      }
  close(TMP);

  @id_db1 = keys %id_db1_2_len;
  $db1_no = $#id_db1 + 1;
  return;
}
########## END read_db1


sub read_db2 {
  my ($i, $j, $k, $ll, $len, $id, $len2);
  open(TMP, $db2) || die "Can not open $db2\n";

  $len = 0;
  $len2 = 0;
  while($ll=<TMP>){
    if ($ll =~ /^>/) {
      if ($len > 0) {
        $id_db2_2_len{$id} = $len;
        $id_db2_2_len_nN{$id} = $len2;
      }
      chop($ll);
      $ll =~ s/\s.+$//;
      $id = substr($ll,1);
      $len = 0;
      $len2 = 0;
    }
    else {
      chop($ll);
      $ll =~ s/\s//;
      $len += length($ll);
      $ll =~ s/N//ig;
      $len2 += length($ll);
    }
  }
      if ($len > 0) {
        $id_db2_2_len{$id} = $len;
        $id_db2_2_len_nN{$id} = $len2;
      }
  close(TMP);

  @id_db2 = keys %id_db2_2_len;
  $db2_no = $#id_db2 + 1;
  return;
}
########## END read_db2


sub usage {
<<EOD;
    -i fasta file1
    -j fasta file2
    -b blast output file 
    -o output
EOD
}
########## END usage
