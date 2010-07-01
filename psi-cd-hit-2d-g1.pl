#!/usr/bin/perl -w

our $script_name = $0;
our $script_dir = $0;
   $script_dir =~ s/[^\/]+$//;
   $script_dir = "./" unless ($script_dir);
require "$script_dir/psi-cd-hit-local.pl";

parse_para_etc(@ARGV);
open_LOG();

our @dess    = ();
our @lens    = ();
our @NR_idx = ();
our $NR_no  = 0;
our %redundant_in_db2 = ();

my $cd_hit_div_pl     = "$script_dir/cd-hit-div.pl";
my ($i, $j, $k, $i0, $j0, $k0, $ll, $cmd);

read_db_no_seq();

@NR_idx = (0..($NR_no-1));
@NR_idx = sort { $lens[$b] <=> $lens[$a] or $a <=> $b } @NR_idx;

if (not -e "$db_in2-0") {
  $cmd = `$cd_hit_div_pl $db_in2 $db_in2 $div_no`;
}
if (not -e "$db_in.pin") {
  blast_formatdb_raw();
}

if (not -e "$db_in2-0.sh") {
  my @commands = ();
  print "Run following commands\n";
  for ($i=0; $i<$div_no; $i++) {
    $cmd = "$blastp -d $db_in -i $db_in2-$i -o $db_in2-$i-bl $bl_para";
    my $pwd = `pwd`;
    my $tsh = "$db_in2-$i.sh";
    open(TSH, "> $tsh") || die;
    print TSH <<EOD;
#!/bin/sh
#\$ -S /bin/bash
#\$ -v PATH,BLASTMAT,BLASTDB

cd $pwd
$cmd
      
EOD
    close(TSH);
    print "qsub -N run-$i $tsh\n";
  }

  die "re-run this script when jobs finished\n\n";
}


for ($i0=0; $i0<$NR_no; $i0++) {
  $i = $NR_idx[$i0];
  my $des = substr($dess[$i],1);
  $redundant_in_db2{$des} = [];
}

for ($i=0; $i<$div_no; $i++) {
  my $bl = "$db_in2-$i-bl";
  process_multi_bl("$bl");
}

write_db_clstr_db3();
close_LOG();

