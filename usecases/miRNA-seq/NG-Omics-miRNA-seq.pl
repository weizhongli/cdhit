#!/usr/bin/perl
################################################################################
# NGS workflow by Weizhong Li, http://weizhongli-lab.org
################################################################################

########## local variables etc. Please edit
$CD_HIT_dir           = "/home/oasis/gordon-data/NGS-ann-project-new/apps/cd-hit-v4.6.8-2017-0621";
$NGS_prog_trimmomatic = "/home/oasis/gordon-data/NGS-ann-project-new/apps/Trimmomatic/trimmomatic-0.32.jar";

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
#PBS -V

#\$ -q RNA.q
#\$ -v PATH
#\$ -V

EOD
};


$NGS_executions{"sh_1"} = {
  "type"                => "sh",
  "cores_per_node"      => 8,
  "number_nodes"        => 1,
};


# $NGS_batch_jobs{"qc-pe"} = {
#   "CMD_opts"          => ["20"],
#   "execution"         => "qsub_1",               # where to execute
#   "cores_per_cmd"     => 4,                     # number of threads used by command below
#   "no_parallel"       => 1,                     # number of total jobs to run using command below
#   "command"           => <<EOD,
# java -jar $NGS_prog_trimmomatic PE -threads 4 -phred33 \\DATA.0 \\DATA.1 \\SELF/R1.fq \\SELF/R1-s.fq \\SELF/R2.fq \\SELF/R2-s.fq \\
#     SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:\\CMDOPTS.0 MAXINFO:80:0.5 1>\\SELF/qc.stdout 2>\\SELF/qc.stderr
# 
# perl -e '\$i=0; while(<>){ if (/^\@/) {\$i++;  print ">Sample|\\SAMPLE|\$i ", substr(\$_,1); \$a=<>; print \$a; \$a=<>; \$a=<>;}}' < \\SELF/R1.fq > \\SELF/R1.fa &
# perl -e '\$i=0; while(<>){ if (/^\@/) {\$i++;  print ">Sample|\\SAMPLE|\$i ", substr(\$_,1); \$a=<>; print \$a; \$a=<>; \$a=<>;}}' < \\SELF/R2.fq > \\SELF/R2.fa &
# 
# wait
# rm -f \\SELF/R1.fq \\SELF/R2.fq \\SELF/R1-s.fq \\SELF/R2-s.fq
# EOD
# };

 
$NGS_batch_jobs{"qc"} = {
  "CMD_opts"          => ["20"],
  "execution"         => "qsub_1",               # where to execute
  "cores_per_cmd"     => 4,                     # number of threads used by command below
  "no_parallel"       => 1,                     # number of total jobs to run using command below
  "command"           => <<EOD,
java -jar $NGS_prog_trimmomatic SE -threads 4 -phred33 \\DATA.0          \\SELF/R1.fq                                            \\
    SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:\\CMDOPTS.0 MAXINFO:80:0.5 1>\\SELF/qc.stdout 2>\\SELF/qc.stderr
perl -e '\$i=0; while(<>){ if (/^\@/) {\$i++;  print ">Sample|\\SAMPLE|\$i ", substr(\$_,1); \$a=<>; print \$a; \$a=<>; \$a=<>;}}' < \\SELF/R1.fq > \\SELF/R1.fa &
rm -f \\SELF/R1.fq
EOD
};


$NGS_batch_jobs{"clstr"} = {
  "injobs"            => ["qc"],
  "CMD_opts"          => ["0.95", "path_to_miRbase", "path_to_spike-in_db"],
  "execution"         => "qsub_1",               # where to execute
  "cores_per_cmd"     => 2,                    # number of threads used by command below
  "no_parallel"       => 1,                    # number of total jobs to run using command below
  "command"           => <<EOD,
#### cluster at 100%
$CD_HIT_dir/cd-hit-est -i \\INJOBS.0/R1.fa -o \\SELF/seq.nr -sf 1 -sc 1 -r 0 -c 1.00        -n 10 -p 1 -d 0 -G 1 -b 1 -T 4 -M 8000 > \\SELF/seq.nr.log
$CD_HIT_dir/usecases/miRNA-seq/filter-small-cluster.pl -i \\SELF/seq.nr.clstr -s \\SELF/seq.nr -o \\SELF/seq.nr-filtered.clstr -f \\SELF/seq.nr-filtered -c 1

$CD_HIT_dir/cd-hit-est -i \\SELF/seq.nr-filtered  -o \\SELF/seq.95  -r 0 -c \\CMDOPTS.0 -n 10 -p 1 -d 0 -G 1 -b 1 -T 4 -M 8000 > \\SELF/seq.95.log
$CD_HIT_dir/clstr_rev.pl \\SELF/seq.nr-filtered.clstr \\SELF/seq.95.clstr     > \\SELF/seq.95-full.clstr

$CD_HIT_dir/cd-hit-est-2d -i \\SELF/seq.95  -i2 \\CMDOPTS.1 -o \\SELF/seq.95.ref -r 1 -c \\CMDOPTS.0 -n 10 -p 1 -d 0 -G 0 -A 20 -b 1 -T 4 -M 8000 > \\SELF/seq.95.ref.log
$CD_HIT_dir/cd-hit-est-2d -i \\SELF/seq.95  -i2 \\CMDOPTS.2 -o \\SELF/seq.95.spk -r 1 -c \\CMDOPTS.0 -n 10 -p 1 -d 0 -G 0 -A 20 -b 1 -T 4 -M 8000 > \\SELF/seq.95.spk.log

$CD_HIT_dir/usecases/Miseq-16S/filter-nontop-ref.pl < \\SELF/seq.95.ref.clstr > \\SELF/seq.95.reftop.clstr
$CD_HIT_dir/usecases/Miseq-16S/filter-nontop-ref.pl < \\SELF/seq.95.spk.clstr > \\SELF/seq.95.spktop.clstr
$CD_HIT_dir/clstr_merge.pl \\SELF/seq.95-full.clstr \\SELF/seq.95.reftop.clstr > \\SELF/tmp.clstr
$CD_HIT_dir/clstr_merge.pl \\SELF/tmp.clstr         \\SELF/seq.95.spktop.clstr > \\SELF/miRNA.clstr

$CD_HIT_dir/clstr_sort_by.pl < \\SELF/miRNA.clstr > \\SELF/miRNA.clstr.s
mv \\SELF/miRNA.clstr.s \\SELF/miRNA.clstr
rm -f \\SELF/tmp.clstr
EOD
};


$NGS_batch_jobs{"clstr-pooled"} = {
  "CMD_opts"          => ["0.95", "path_to_miRbase", "path_to_spike-in_db"],
  "execution"         => "qsub_1",               # where to execute
  "cores_per_cmd"     => 2,                    # number of threads used by command below
  "no_parallel"       => 1,                    # number of total jobs to run using command below
  "command"           => <<EOD,
#### concat seq.nr-filtered seq.nr-filtered.clstr 

$CD_HIT_dir/cd-hit-est -i seq.nr-filtered    -o seq.95             -r 0 -c \\CMDOPTS.0 -n 10 -p 1 -d 0 -G 1 -b 1 -T 4 -M 8000 > seq.95.log
$CD_HIT_dir/clstr_rev.pl seq.nr-filtered.clstr       seq.95.clstr     > seq.95-full.clstr

$CD_HIT_dir/cd-hit-est-2d -i seq.95  -i2 \\CMDOPTS.1 -o seq.95.ref -r 1 -c \\CMDOPTS.0 -n 10 -p 1 -d 0 -G 0 -A 20 -b 1 -T 4 -M 8000 > seq.95.ref.log
$CD_HIT_dir/cd-hit-est-2d -i seq.95  -i2 \\CMDOPTS.2 -o seq.95.spk -r 1 -c \\CMDOPTS.0 -n 10 -p 1 -d 0 -G 0 -A 20 -b 1 -T 4 -M 8000 > seq.95.spk.log

$CD_HIT_dir/usecases/Miseq-16S/filter-nontop-ref.pl < seq.95.ref.clstr > seq.95.reftop.clstr
$CD_HIT_dir/usecases/Miseq-16S/filter-nontop-ref.pl < seq.95.spk.clstr > seq.95.spktop.clstr
$CD_HIT_dir/clstr_merge.pl seq.95-full.clstr seq.95.reftop.clstr > tmp.clstr
$CD_HIT_dir/clstr_merge.pl tmp.clstr         seq.95.spktop.clstr > miRNA.clstr

$CD_HIT_dir/clstr_sort_by.pl < miRNA.clstr > miRNA.clstr.s
mv miRNA.clstr.s miRNA.clstr

$CD_HIT_dir/usecases/miRNA-seq/clstr_2_miRNA-table.pl -i miRNA.clstr -o miRNA.txt

EOD
};

##############################################################################################
########## END
1;

