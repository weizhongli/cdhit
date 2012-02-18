#!/usr/bin/perl

#usage: clstr_xml.pl [-len|-size] level1.clstr [level2.clstr level3.clstr ...]
#purpose: to create xml file from cd-hit or hierarchical cd-hit(h-cd-hit) results
#example:
#  if nr90 from nr100 and nr80 from nr90, with nr90.clstr and nr80.clstr
#  use
#     clstr_xml.pl -len nr90.clstr nr80.clstr
#  to create an xml file 

my $n=@ARGV;
my $option="-len";
if($n==0||($n==1&&$ARGV[0]=~/^-/)){print "Usage:\n\tclstr2xml.pl [-len|-size] input1.clstr [input2.clstr input3.clstr ...]\n";exit(0);}
if($ARGV[0]=~/^-/){
  $option=$ARGV[0];shift @ARGV;
  $n--;
}

my $topsize=0;
my $totalsequence=0;
for($ilevel=$n;$ilevel>0;$ilevel--){

  open(TMP, $ARGV[$ilevel-1]) || die "Can not open file";
  $readin = 0;
  my $gi = "";
  my $clstr = "";
  my $this_no = 0;
  while(my $ll=<TMP>) {
    if ($ll =~ /^>/ ) {
      if ($readin and $key and $this_no>0 ) {
        $nextlevel=$ilevel-1;
        $childofgi{"$key##$ilevel"}=[$this_no,$key,@members];
        if($ilevel>1){
          foreach $k (@members){$parent{"$k##$nextlevel"}=$key;}
          $parent{"$key##$nextlevel"}=$key;
        }
        $size{"$key##$ilevel"}=$this_no;
        $totalsequence+=$this_no;
        if($ilevel<$n){
          $totalsequence--;
          my $j=$ilevel;
          while(exists($parent{"$key##$j"})){
            $key=$parent{"$key##$j"};$j++;
            $size{"$key##$j"}+=($this_no-1);
          }
        }
      }
      $key="";
      $clstr="";
      $this_no=0;
      @members=();
      if($ilevel==$n){$topsize++;}
    }
#>Cluster 10
#0       187aa, >BPHE_BURCE\1-187... at 1:187:9:195/96.79%
#1       188aa, >Q7B113_9PSED\1-188... at 1:188:8:195/96.81%
#2       195aa, >Q9AEY4_9PSED\19-213... *

    else {
      $readin = 1;
      #$clstr .= $ll;
      if ($ll =~ /\*/ and $ll =~ /\d+\t(\d+)[a-z]{2}, >(.+)\.\.\./ ) { $key = $2; $len{$key}=$1;}
      else{$ll =~ /\d+\t(\d+)[a-z]{2}, >(.+)\.\.\./; $gi=$2;push(@members,$gi);$len{$gi}=$1;
        $ll=~/(\d+%|\d+\.\d+%)$/;$identity{"$gi##$ilevel"}=$1;
      }
      $this_no++;
    }
  }
  close(TMP);
  
  #last group
  if ($readin and $key and $this_no>0 ) {
    $nextlevel=$ilevel-1;
    $childofgi{"$key##$ilevel"}=[$this_no,$key,@members];
    if($ilevel>1){
      foreach $k (@members){$parent{"$k##$nextlevel"}=$key;}
      $parent{"$key##$nextlevel"}=$key;
    }
    $size{"$key##$ilevel"}=$this_no;
    $totalsequence+=$this_no;
    if($ilevel<$n){
      $totalsequence--;
      my $j=$ilevel;
      while(exists($parent{"$key##$j"})){
        $key=$parent{"$key##$j"};$j++;
        $size{"$key##$j"}+=($this_no-1);
      }
    }
  }#end of last group
}

print <<head;
<?xml version="1.0" encoding="ISO-8859-1"?>
<CDHIT_Clusters GroupNo="$topsize" SequenceNo="$totalsequence" xmlns="weizhong-lab.ucsd.edu">
head

my @allkeys;
if($option eq "-len"){
  @allkeys=sort {$len{substr($b,0,rindex($b,"##"))} <=> $len{substr($a,0,rindex($a,"##"))}} keys %childofgi;
}
elsif($option eq "-size"){
  @allkeys=sort {$size{$b} <=> $size{$a}} keys %childofgi;
}

foreach $k (@allkeys){
  printnode($k,$n);
}
print "</CDHIT_Clusters>\n";
exit(0);

sub printnode{
  my $ckey=shift;
  my $level=shift;
  my ($key,$keylevel)=split(/##/,$ckey);
  if($keylevel != $level){return;}
  my $space=" "x(3*($n-$level));
  my $mspace="$space   ";
  
  print "$space<representativeMember ";
  print "Level=\"$level\" ";
  if($level>1){print "GroupNo=\"$childofgi{$ckey}->[0]\" ";}
  print "SequenceNo=\"$size{$ckey}\" Length=\"$len{$key}\"";
  my $nextlevel=$level-1;
  my $upperlevel=$level+1;
  my $identitykey="$key##$upperlevel";
  if(exists($identity{$identitykey})){print " Identity=\"$identity{$identitykey}\"";}
  else{if($keylevel<$n){print " Identity=\"100%\"";}}
  print ">$key\n";
  
  my @ids=@{$childofgi{$ckey}};shift @ids;
  if($option eq "-len"||$nextlevel==0){
    @ids=sort {$len{$b} <=> $len{$a}} @ids;
  }
  elsif($option eq "-size"){
    @ids=sort {$size{"$b##$nextlevel"} <=> $size{"$a##$nextlevel"}} @ids;
  }
  for(my $i=0;$i<$childofgi{$ckey}->[0];$i++){
    if(!exists($childofgi{"$ids[$i]##$nextlevel"})) {
      print "$mspace<member Length=\"$len{$ids[$i]}\"";
      $identitykey="$ids[$i]##$level";
      if(exists($identity{$identitykey})){print " Identity=\"$identity{$identitykey}\"";}
      else{if($keylevel<$n||$n==1){print " Identity=\"100%\"";}}
      print ">$ids[$i]</member>\n";
    }
    else{
      $newkey="$ids[$i]##$nextlevel";
      printnode($newkey,$nextlevel);
    }
  }
  print "$space</representativeMember>\n";
}
