#!/usr/bin/perl
################################################################################
# NGS workflow by Weizhong Li, http://weizhongli-lab.org
################################################################################

########## local variables etc. Please edit
$NGS_root     = "/home/oasis/data/NGS-ann-project";
$CD_HIT_dir   = "/home/oasis/data/etc/git/cdhit";

########## more local variables, do not edit next three lines
$NGS_tool_dir = "$NGS_root/NGS-tools";
$NGS_prog_dir = "$NGS_root/apps";
$NGS_bin_dir  = "$NGS_root/apps/bin";
$NGS_ref_dir  = "$NGS_root/refs";

########## computation resources for execution of jobs
%NGS_executions = ();
$NGS_executions{"qsub_1"} = {
  "type"                => "qsub-pe",
  "cores_per_node"      => 8,
  "number_nodes"        => 64,
  "user"                => "weizhong", #### I will use command such as qstat -u weizhong to query submitted jobs
  "command"             => "qsub",
  "command_name_opt"    => "-N",
  "command_err_opt"     => "-e",
  "command_out_opt"     => "-o",
  "template"            => <<EOD,
#!/bin/sh
#PBS -v PATH
#PBS -M liwz\@sdsc.edu
#PBS -q normal
#PBS -V
#PBS -l nodes=1:ppn=16,walltime=48:00:00,mem=60000mb

#\$ -v PATH
#\$ -V

EOD
};


$NGS_executions{"sh_1"} = {
  "type"                => "sh",
  "cores_per_node"      => 32,
  "number_nodes"        => 1,
};

$NGS_batch_jobs{"qc"} = {
  "CMD_opts"          => ["100"],
  "execution"         => "sh_1",               # where to execute
  "cores_per_cmd"     => 4,                     # number of threads used by command below
  "no_parallel"       => 1,                     # number of total jobs to run using command below
  "command"           => <<EOD,
java -jar $NGS_prog_dir/Trimmomatic/trimmomatic-0.32.jar PE -threads 4 -phred33 \\DATA.0 \\DATA.1 \\SELF/R1.fq \\SELF/R1-s.fq \\SELF/R2.fq \\SELF/R2-s.fq \\
    SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:\\CMDOPTS.0 MAXINFO:80:0.5 1>\\SELF/qc.stdout 2>\\SELF/qc.stderr

perl -e '\$i=0; while(<>){ if (/^\@/) {\$i++;  print ">Sample|\\SAMPLE|\$i ", substr(\$_,1); \$a=<>; print \$a; \$a=<>; \$a=<>;}}' < \\SELF/R1.fq > \\SELF/R1.fa &
perl -e '\$i=0; while(<>){ if (/^\@/) {\$i++;  print ">Sample|\\SAMPLE|\$i ", substr(\$_,1); \$a=<>; print \$a; \$a=<>; \$a=<>;}}' < \\SELF/R2.fq > \\SELF/R2.fa &

wait
rm -f \\SELF/R1.fq \\SELF/R2.fq \\SELF/R1-s.fq \\SELF/R2-s.fq
EOD
};

$NGS_batch_jobs{"otu"} = {
  "injobs"            => ["qc"],
  "CMD_opts"          => ["100","100", "0.97","0.0001","/home/oasis/data/projects/USDA/PE-test/plos-one-PRJEB4688/gg_13_5-PE99.100-R1",
                                                       "/home/oasis/data/projects/USDA/PE-test/plos-one-PRJEB4688/gg_13_5-PE99.100-R2", "75"],
  "execution"         => "sh_1",               # where to execute
  "cores_per_cmd"     => 2,                    # number of threads used by command below
  "no_parallel"       => 1,                    # number of total jobs to run using command below
  "command"           => <<EOD,
$CD_HIT_dir/cd-hit-est -i \\INJOBS.0/R1.fa -j \\INJOBS.0/R2.fa -o \\SELF/seq.nr -op \\SELF/seq.nr.2 -sf 1 -sc 1 -P 1 -r 0 \\
    -cx \\CMDOPTS.0 -cy \\CMDOPTS.1 -c 1.0  -n 10 -G 1 -b 1  -T 1 -M 8000  -d 0 -p 1 > \\SELF/seq.nr.log
$CD_HIT_dir/cd-hit-est -i \\SELF/seq.nr   -o \\SELF/seq.chimeric-clstr.R1 -r 0 -cx \\CMDOPTS.6 -c 0.99 -n 10 -G 1 -b 1  -T 1 -M 8000  -d 0 -p 1 > \\SELF/seq.chimeric-clstr.R1.log
$CD_HIT_dir/cd-hit-est -i \\SELF/seq.nr.2 -o \\SELF/seq.chimeric-clstr.R2 -r 0 -cx \\CMDOPTS.6 -c 0.99 -n 10 -G 1 -b 1  -T 1 -M 8000  -d 0 -p 1 > \\SELF/seq.chimeric-clstr.R2.log
$CD_HIT_dir/cd-hit-est -i \\SELF/seq.nr -j \\SELF/seq.nr.2 -o \\SELF/seq.99 -op \\SELF/seq.99.2 -P 1 -r 0 \\
    -cx \\CMDOPTS.0 -cy \\CMDOPTS.1 -c 0.99 -n 10 -G 1 -b 1  -T 1 -M 8000  -d 0 -p 1 > \\SELF/seq.99.log
$CD_HIT_dir/usecases/Miseq-16S/filter-chimeric-and-small.pl -c \\CMDOPTS.3 -k \\SELF/seq.nr.clstr \\
    -i \\SELF/seq.chimeric-clstr.R1.clstr -j \\SELF/seq.chimeric-clstr.R2.clstr \\
    -a \\SELF/seq.99.clstr -f \\SELF/seq.99 -g \\SELF/seq.99.2 -o \\SELF/seq.99f
cat \\CMDOPTS.4 \\SELF/seq.99f   > \\SELF/seq.99fwref
cat \\CMDOPTS.5 \\SELF/seq.99f.2 > \\SELF/seq.99fwref.2
$CD_HIT_dir/cd-hit-est -i \\SELF/seq.99fwref -j \\SELF/seq.99fwref.2 -o \\SELF/seq.97fwref -op \\SELF/seq.97fwref.2 -P 1 -r 0 \\
    -cx \\CMDOPTS.0 -cy \\CMDOPTS.1 -c 0.97 -n 10 -G 1 -b 10  -T 1 -M 8000  -d 0 -p 1 > \\SELF/seq.97fwref.log
$CD_HIT_dir/usecases/Miseq-16S/filter-refonly-cluster.pl < \\SELF/seq.97fwref.clstr > \\SELF/seq.97fwref.clstr.2
mv \\SELF/seq.97fwref.clstr.2 \\SELF/seq.97fwref.clstr
$CD_HIT_dir/clstr_rev.pl \\SELF/seq.nr.clstr       \\SELF/seq.99f.clstr     > \\SELF/seq.99f-full.clstr
$CD_HIT_dir/clstr_rev.pl \\SELF/seq.99f-full.clstr \\SELF/seq.97fwref.clstr > \\SELF/seq.97f-full.clstr
$CD_HIT_dir/clstr_sort_by.pl < \\SELF/seq.97f-full.clstr > \\SELF/seq.97f-full.clstr.s
mv \\SELF/seq.97f-full.clstr.s \\SELF/seq.97f-full.clstr
EOD
};



##############################################################################################
########## END
1;

