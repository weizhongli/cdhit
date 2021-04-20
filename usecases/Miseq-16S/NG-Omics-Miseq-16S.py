#!/usr/bin/python
################################################################################
# NGS workflow by Weizhong Li, http://weizhongli-lab.org
################################################################################

queue_system = 'SGE'

########## local variables etc. Please edit
ENV={
  'CD_HIT_dir'           : '/home/oasis/data/NGS-ann-project/apps/cdhit.git',
  'NGS_prog_trimmomatic' : '/home/oasis/data/NGS-ann-project/apps/Trimmomatic/trimmomatic-0.32.jar',
}

########## computation resources for execution of jobs
NGS_executions = {}
NGS_executions['qsub_1'] = {
  'type'                : 'qsub-pe',
  'cores_per_node'      : 8,
  'number_nodes'        : 64,
  'user'                : 'weizhong', #### I will use command such as qstat -u weizhong to query submitted jobs
  'command'             : 'qsub',
  'command_name_opt'    : '-N',
  'command_err_opt'     : '-e',
  'command_out_opt'     : '-o',
  'template'            : '''#!/bin/bash
##$ -q RNA.q
#$ -v PATH
#$ -V

'''
}

NGS_executions['sh_1'] = {
  'type'                : 'sh',
  'cores_per_node'      : 8,
  'number_nodes'        : 1,
  'template'            : '''#!/bin/bash

'''
}


NGS_batch_jobs = {}
NGS_batch_jobs['qc'] = {
  'CMD_opts'         : ['100'],
  'non_zero_files'   : ['R1.fa.gz','R2.fa.gz'],
  'execution'        : 'qsub_1',               # where to execute
  'cores_per_cmd'    : 4,                    # number of threads used by command below
  'no_parallel'      : 1,                    # number of total jobs to run using command below
  'command'          : '''
java -jar $ENV.NGS_prog_trimmomatic PE -threads 4 -phred33 $DATA.0 $DATA.1 $SELF/R1.fq $SELF/R1-s.fq $SELF/R2.fq $SELF/R2-s.fq \\
    SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:$CMDOPTS.0 MAXINFO:80:0.5 1>$SELF/qc.stdout 2>$SELF/qc.stderr

perl -e '$i=0; while(<>){ if (/^@/) {$i++;  print ">Sample|$SAMPLE|$i ", substr($_,1); $a=<>; print $a; $a=<>; $a=<>;}}' < $SELF/R1.fq > $SELF/R1.fa &
perl -e '$i=0; while(<>){ if (/^@/) {$i++;  print ">Sample|$SAMPLE|$i ", substr($_,1); $a=<>; print $a; $a=<>; $a=<>;}}' < $SELF/R2.fq > $SELF/R2.fa &
wait
gzip $SELF/R1.fa &
gzip $SELF/R2.fa &
wait

rm -f $SELF/R1.fq $SELF/R2.fq $SELF/R1-s.fq $SELF/R2-s.fq
'''
}


NGS_batch_jobs['otu'] = {
  'injobs'            : ['qc'],
  'non_zero_files'    : ['seq.99f','seq.99f.2','seq.99f-all.clstr','pool.ok'],
  'CMD_opts'          : ['150', '100', '0.0005', '75', 'path_to_pooled_sample_dir'],
  'execution'         : 'qsub_1',               # where to execute
  'cores_per_cmd'     : 2,                    # number of threads used by command below
  'no_parallel'       : 1,                    # number of total jobs to run using command below
  'command'           : '''

#### 1. cluster at 100% PE
$ENV.CD_HIT_dir/cd-hit-est -i $INJOBS.0/R1.fa.gz -j $INJOBS.0/R2.fa.gz -o $SELF/seq.nr -op $SELF/seq.nr.2 -sf 1 -sc 1 -P 1 -r 0 \\
    -cx $CMDOPTS.0 -cy $CMDOPTS.1 -c 1.0  -n 10 -G 1 -b 1  -T 1 -M 8000  -d 0 -p 1 > $SELF/seq.nr.log

#### 2. cluster at 99% PE 
$ENV.CD_HIT_dir/cd-hit-est -i $SELF/seq.nr -j $SELF/seq.nr.2 -o $SELF/seq.99 -op $SELF/seq.99.2 -P 1 -r 0 \\
    -cx $CMDOPTS.0 -cy $CMDOPTS.1 -c 0.99 -n 10 -G 1 -b 1  -T 1 -M 8000  -d 0 -p 1 > $SELF/seq.99.log

#### 3. cluster at 99% SE for R1, R2 
$ENV.CD_HIT_dir/cd-hit-est -i $SELF/seq.nr   -o $SELF/seq.chimeric-clstr.R1 -r 0 -cx $CMDOPTS.3 -c 0.99 -n 10 -G 0 -b 1 -A 50 -T 1 -M 8000  -d 0 -p 1 > $SELF/seq.chimeric-clstr.R1.log
$ENV.CD_HIT_dir/cd-hit-est -i $SELF/seq.nr.2 -o $SELF/seq.chimeric-clstr.R2 -r 0 -cx $CMDOPTS.3 -c 0.99 -n 10 -G 0 -b 1 -A 50 -T 1 -M 8000  -d 0 -p 1 > $SELF/seq.chimeric-clstr.R2.log
rm -f $SELF/seq.chimeric-clstr.R1    $SELF/seq.chimeric-clstr.R1.log \\
      $SELF/seq.chimeric-clstr.R2    $SELF/seq.chimeric-clstr.R2.log 

#### 4. 5. filter chimeric sequences and sequences in small clusters
$ENV.CD_HIT_dir/usecases/Miseq-16S/filter-chimeric-and-small.pl -c $CMDOPTS.2 -k $SELF/seq.nr.clstr \\
    -i $SELF/seq.chimeric-clstr.R1.clstr -j $SELF/seq.chimeric-clstr.R2.clstr \\
    -a $SELF/seq.99.clstr -f $SELF/seq.99 -g $SELF/seq.99.2 -o $SELF/seq.99f
$ENV.CD_HIT_dir/clstr_rev.pl $SELF/seq.nr.clstr       $SELF/seq.99f.clstr     > $SELF/seq.99f-all.clstr
mv $SELF/seq.99f.log $SELF/chimeric-small-clusters-list.txt


####
if [ ! -e "$CMDOPTS.4" ]; then
  mkdir -p $CMDOPTS.4
fi

i="0"
while [ 1 ]; do

  if [ -e "$CMDOPTS.4/lock" ]; then
    echo "wait $CMDOPTS.4/lock"
    sleep 5
  else
    date > $CMDOPTS.4/lock
    
    cat $SELF/seq.99f                          >> $CMDOPTS.4/seq.99f
    cat $SELF/seq.99f.2                        >> $CMDOPTS.4/seq.99f.2
    cat $SELF/seq.99f-all.clstr                >> $CMDOPTS.4/seq.99f-all.clstr
    cat $SELF/chimeric-small-clusters-list.txt >> $CMDOPTS.4/chimeric-small-clusters-list.txt
    date > $SELF/pool.ok
    sleep 1

    rm -f $CMDOPTS.4/lock
    break
  fi

  i=$[$i+1]
  if [ "$i" -gt "50" ]; then
    echo "wait $CMDOPTS.4/lock for too long"
    break
  fi
done

'''
}


NGS_batch_jobs['otu-pooled'] = {
  'CMD_opts'          : ['150', '100', '0.97', 'path_to_spliced_ref_db-R1', 'path_to_spliced_ref_db-R1'],
  'non_zero_files'    : ['OTU.txt'],
  'execution'         : 'qsub_1',               # where to execute
  'cores_per_cmd'     : 2,                    # number of threads used by command below
  'no_parallel'       : 1,                    # number of total jobs to run using command below
  'command'           : '''
#### before running
#### concat seq.99f seq.99f.2 seq.99f-all.clstr chimeric-small-clusters-list.txt
$ENV.CD_HIT_dir/cd-hit-est -i seq.99f -j seq.99f.2 -o seq.97 -op seq.97.2 -P 1 -r 0 \\
    -cx $CMDOPTS.0 -cy $CMDOPTS.1 -c $CMDOPTS.2  -n 10 -G 1 -b 10  -T 1 -M 8000  -d 0 -p 1 > seq.97.log
$ENV.CD_HIT_dir/cd-hit-est-2d -i seq.97 -j seq.97.2 -i2 $CMDOPTS.3 -j2 $CMDOPTS.4 -o seq.97.ref -op seq.97.ref.2 -P 1 -r 0 \\
    -cx $CMDOPTS.0 -cy $CMDOPTS.1 -c $CMDOPTS.2  -n 10 -G 1 -b 10  -T 1 -M 8000  -d 0 -p 1 > seq.97.ref.log
$ENV.CD_HIT_dir/clstr_rev.pl seq.99f-all.clstr seq.97.clstr > seq.97-all.clstr
$ENV.CD_HIT_dir/usecases/Miseq-16S/filter-nontop-ref.pl < seq.97.ref.clstr > seq.97.reftop.clstr
$ENV.CD_HIT_dir/clstr_merge.pl seq.97-all.clstr seq.97.reftop.clstr > OTU.clstr
$ENV.CD_HIT_dir/usecases/Miseq-16S/clstr_2_OTU_table.pl -i OTU.clstr -o OTU.txt
rm -f seq.97.ref seq.97.ref.2 seq.97.ref.log

'''
}

