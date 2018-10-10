#!/usr/bin/perl

if(@ARGV==0){
   print "Usage:\n\tclstr_sql_tbl_sort.pl table_file level\n";
   exit(1);
}

my $clstr_file = shift;
my $level = shift;
my ($i, $j, $k, $ll);

my @lls = ();
  open(TMP, $clstr_file) || die;
  while($ll=<TMP>){
    chop($ll);
    push(@lls, $ll);
  }
  close(TMP);
  print STDERR "done reading $clstr_file\n";

   @sp=split(/\t/,$lls[0]);
   if(@sp < $level*2+2){
      print "error level\n";
      exit(1);
   }

  if ($level == 1) {
    @lls = sort {@a=split(/\t/,$a); @b=split(/\t/,$b);
                 $a[-2] <=> $b[-2] or $b[-1] <=> $a[-1] } @lls;
  }
  elsif ($level == 2) {
    @lls = sort {@a=split(/\t/,$a); @b=split(/\t/,$b);
                 $a[-2] <=> $b[-2] or $b[-1] <=> $a[-1] or 
                 $a[-4] <=> $b[-4] or $b[-3] <=> $a[-3]} @lls;
  }
  elsif ($level == 3) {
    @lls = sort {@a=split(/\t/,$a); @b=split(/\t/,$b);
                 $a[-2] <=> $b[-2] or $b[-1] <=> $a[-1] or
                 $a[-4] <=> $b[-4] or $b[-3] <=> $a[-3] or
                 $a[-6] <=> $b[-6] or $b[-5] <=> $a[-5]} @lls;
  }

foreach $ll (@lls) {
  print $ll , "\n";
}
