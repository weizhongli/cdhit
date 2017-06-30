#!/usr/bin/perl
################################################################################
# NGS workflow by Weizhong Li, http://weizhongli-lab.org
################################################################################

########## local variables etc. Please edit
$CD_HIT_dir           = "/home/oasis/data/etc/git/cdhit";
$NGS_prog_trimmomatic = "/home/oasis/data/NGS-ann-project/apps/Trimmomatic/trimmomatic-0.32.jar";


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

#\$ -v PATH
#\$ -V

EOD
};


$NGS_executions{"sh_1"} = {
  "type"                => "sh",
  "cores_per_node"      => 8,
  "number_nodes"        => 1,
};

$NGS_batch_jobs{"qc"} = {
  "CMD_opts"          => ["100"],
  "execution"         => "sh_1",               # where to execute
  "cores_per_cmd"     => 4,                     # number of threads used by command below
  "no_parallel"       => 1,                     # number of total jobs to run using command below
  "command"           => <<EOD,
java -jar $NGS_prog_trimmomatic PE -threads 4 -phred33 \\DATA.0 \\DATA.1 \\SELF/R1.fq \\SELF/R1-s.fq \\SELF/R2.fq \\SELF/R2-s.fq \\
    SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:\\CMDOPTS.0 MAXINFO:80:0.5 1>\\SELF/qc.stdout 2>\\SELF/qc.stderr

perl -e '\$i=0; while(<>){ if (/^\@/) {\$i++;  print ">Sample|\\SAMPLE|\$i ", substr(\$_,1); \$a=<>; print \$a; \$a=<>; \$a=<>;}}' < \\SELF/R1.fq > \\SELF/R1.fa &
perl -e '\$i=0; while(<>){ if (/^\@/) {\$i++;  print ">Sample|\\SAMPLE|\$i ", substr(\$_,1); \$a=<>; print \$a; \$a=<>; \$a=<>;}}' < \\SELF/R2.fq > \\SELF/R2.fa &

wait
rm -f \\SELF/R1.fq \\SELF/R2.fq \\SELF/R1-s.fq \\SELF/R2-s.fq
EOD
};


$NGS_batch_jobs{"otu"} = {
  "injobs"            => ["qc"],
  "CMD_opts"          => ["150", "100", "0.97", "0.0001", "path_to_spliced_ref_db-R1", "path_to_spliced_ref_db-R1", "75"],
  "execution"         => "sh_1",               # where to execute
  "cores_per_cmd"     => 2,                    # number of threads used by command below
  "no_parallel"       => 1,                    # number of total jobs to run using command below
  "command"           => <<EOD,
#### cluster at 100% PE
$CD_HIT_dir/cd-hit-est -i \\INJOBS.0/R1.fa -j \\INJOBS.0/R2.fa -o \\SELF/seq.nr -op \\SELF/seq.nr.2 -sf 1 -sc 1 -P 1 -r 0 \\
    -cx \\CMDOPTS.0 -cy \\CMDOPTS.1 -c 1.0  -n 10 -G 1 -b 1  -T 1 -M 8000  -d 0 -p 1 > \\SELF/seq.nr.log
#### cluster at 99% PE and SE for R1,R2 
$CD_HIT_dir/cd-hit-est -i \\SELF/seq.nr   -o \\SELF/seq.chimeric-clstr.R1 -r 0 -cx \\CMDOPTS.6 -c 0.99 -n 10 -G 0 -b 1 -A 50 -T 1 -M 8000  -d 0 -p 1 > \\SELF/seq.chimeric-clstr.R1.log
$CD_HIT_dir/cd-hit-est -i \\SELF/seq.nr.2 -o \\SELF/seq.chimeric-clstr.R2 -r 0 -cx \\CMDOPTS.6 -c 0.99 -n 10 -G 0 -b 1 -A 50 -T 1 -M 8000  -d 0 -p 1 > \\SELF/seq.chimeric-clstr.R2.log
$CD_HIT_dir/cd-hit-est -i \\SELF/seq.nr -j \\SELF/seq.nr.2 -o \\SELF/seq.99 -op \\SELF/seq.99.2 -P 1 -r 0 \\
    -cx \\CMDOPTS.0 -cy \\CMDOPTS.1 -c 0.99 -n 10 -G 1 -b 1  -T 1 -M 8000  -d 0 -p 1 > \\SELF/seq.99.log
$CD_HIT_dir/usecases/Miseq-16S/filter-chimeric-and-small.pl -c \\CMDOPTS.3 -k \\SELF/seq.nr.clstr \\
    -i \\SELF/seq.chimeric-clstr.R1.clstr -j \\SELF/seq.chimeric-clstr.R2.clstr \\
    -a \\SELF/seq.99.clstr -f \\SELF/seq.99 -g \\SELF/seq.99.2 -o \\SELF/seq.99f
$CD_HIT_dir/clstr_rev.pl \\SELF/seq.nr.clstr       \\SELF/seq.99f.clstr     > \\SELF/seq.99f-all.clstr
$CD_HIT_dir/cd-hit-est -i \\SELF/seq.99f -j \\SELF/seq.99f.2 -o \\SELF/seq.97 -op \\SELF/seq.97.2 -P 1 -r 0 \\
    -cx \\CMDOPTS.0 -cy \\CMDOPTS.1 -c 0.97 -n 10 -G 1 -b 10  -T 1 -M 8000  -d 0 -p 1 > \\SELF/seq.97.log
$CD_HIT_dir/cd-hit-est-2d -i \\SELF/seq.97 -j \\SELF/seq.97.2 -i2 \\CMDOPTS.4 -j2 \\CMDOPTS.5 -o \\SELF/seq.97.ref -op \\SELF/seq.97.ref.2 -P 1 -r 0 \\
    -cx \\CMDOPTS.0 -cy \\CMDOPTS.1 -c 0.97 -n 10 -G 1 -b 10  -T 1 -M 8000  -d 0 -p 1 > \\SELF/seq.97.ref.log
$CD_HIT_dir/clstr_rev.pl \\SELF/seq.99f-all.clstr \\SELF/seq.97.clstr > \\SELF/seq.97-all.clstr
$CD_HIT_dir/usecases/Miseq-16S/filter-nontop-ref.pl < \\SELF/seq.97.ref.clstr > \\SELF/seq.97.reftop.clstr
$CD_HIT_dir/clstr_merge.pl \\SELF/seq.97-all.clstr \\SELF/seq.97.reftop.clstr > \\SELF/OTU.clstr

rm -f \\SELF/seq.chimeric-clstr.R1    \\SELF/seq.chimeric-clstr.R1.log    \\SELF/seq.chimeric-clstr.R2    \\SELF/seq.chimeric-clstr.R2.log 
rm -f \\SELF/seq.97.ref \\SELF/seq.97.ref.2 \\SELF/seq.97.ref.log
mv \\SELF/seq.99f.log \\SELF/chimeric-small-clusters-list.txt

EOD
};


$NGS_batch_jobs{"otu-pooled"} = {
  "CMD_opts"          => ["150", "100", "0.97", "0.0001", "path_to_spliced_ref_db-R1", "path_to_spliced_ref_db-R1", "75"],
  "execution"         => "sh_1",               # where to execute
  "cores_per_cmd"     => 2,                    # number of threads used by command below
  "no_parallel"       => 1,                    # number of total jobs to run using command below
  "command"           => <<EOD,
#### before running
#### concat seq.99f seq.99f.2 seq.99f-all.clstr chimeric-small-clusters-list.txt
$CD_HIT_dir/cd-hit-est -i seq.99f -j seq.99f.2 -o seq.97 -op seq.97.2 -P 1 -r 0 \\
    -cx \\CMDOPTS.0 -cy \\CMDOPTS.1 -c 0.97 -n 10 -G 1 -b 10  -T 1 -M 8000  -d 0 -p 1 > seq.97.log
$CD_HIT_dir/cd-hit-est-2d -i seq.97 -j seq.97.2 -i2 \\CMDOPTS.4 -j2 \\CMDOPTS.5 -o seq.97.ref -op seq.97.ref.2 -P 1 -r 0 \\
    -cx \\CMDOPTS.0 -cy \\CMDOPTS.1 -c 0.97 -n 10 -G 1 -b 10  -T 1 -M 8000  -d 0 -p 1 > seq.97.ref.log
$CD_HIT_dir/clstr_rev.pl seq.99f-all.clstr seq.97.clstr > seq.97-all.clstr
$CD_HIT_dir/usecases/Miseq-16S/filter-nontop-ref.pl < seq.97.ref.clstr > seq.97.reftop.clstr
$CD_HIT_dir/clstr_merge.pl seq.97-all.clstr seq.97.reftop.clstr > OTU.clstr
$CD_HIT_dir/usecases/Miseq-16S/clstr_2_OTU_table.pl -i OTU.clstr -o OTU.txt
rm -f seq.97.ref seq.97.ref.2 seq.97.ref.log

EOD
};

##############################################################################################
########## END
1;

