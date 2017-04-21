#!/usr/bin/perl
# =============================== NG-Omics-WF ==================================
#  _   _  _____         ____            _              __          ________ 
# | \ | |/ ____|       / __ \          (_)             \ \        / /  ____|
# |  \| | |  __ ______| |  | |_ __ ___  _  ___ ___ _____\ \  /\  / /| |__   
# | . ` | | |_ |______| |  | | '_ ` _ \| |/ __/ __|______\ \/  \/ / |  __|  
# | |\  | |__| |      | |__| | | | | | | | (__\__ \       \  /\  /  | |     
# |_| \_|\_____|       \____/|_| |_| |_|_|\___|___/        \/  \/   |_|     
#                                                                           
# =========================== Next Generation Omics data workflow tools ========
#
# Workflow tools for next generation genomics, metagenomics, RNA-seq 
# and other type of omics data analyiss, 
#
# Software originally developed since 2010 by Weizhong Li at UCSD
#                                               currently at JCVI
#
# http://weizhongli-lab.org/ngomicswf           liwz@sdsc.edu
# ==============================================================================

use Getopt::Std;
use POSIX;

getopts("i:R:s:J:Q:r:j:Z:t:S:T:",\%opts);
die usage() unless ($opts{i} and ($opts{s} or $opts{S}));

my $sample_in    = $opts{s};
my $sample_command_in = $opts{S}; #### ';' delimited samples, ':' delimited entries, e.g. sample1:R1.fq:R2.fq;sample2:R1.fq:R2.fq   or sample1;sample2;sample3
my $input_conf   = $opts{i};
my $this_task    = $opts{J};
our $G_NGS_root   = $opts{r};
my $queue_system = $opts{Q}; $queue_system = "SGE" unless $queue_system;
my $subset_wfs   = $opts{R};
my $subset_jobs  = $opts{j};
my $second_opt   = $opts{Z};
my $opt_file     = $opts{t};
my $opt_command_in = $opts{T}; #### ';' delimited jobs, ":" delimited entries, e.g. JobID_A:opt0:opt1:opt2;JobID_B:opt0:opt1

my $pwd            = `pwd`; chop($pwd);
my $sleep_time_min = 15;
my $sleep_time_max = 120;
my $log_dir        = "$pwd/WF-LOG";
my $log_file       = "$log_dir/LOG";
my $log_fileq      = "$log_dir/LOGq";
my $sh_dir         = "$pwd/WF-sh";
my $sh_bundle_dir  = "$pwd/WF-sh-bundle";
my $subset_flag    = 0;   #### run only one job, subset of jobs, or jobs in sub workflows
my %subset_jobs    = ();
my %qstat_xml_data = ();
my ($i, $j, $k, $ll, $cmd);

######## scan through WF configration 
######## and generate job list
require $input_conf;
my %job_list = (); # as $job_list{$t_job_id}{$t_sample_id} = {};
my ($t_sample_id, $t_job_id, $t_execution_id);
my ($t_sample,    $t_job,    $t_execution);
task_level_jobs();
my @NGS_batch_jobs = sort {($NGS_batch_jobs{$a}->{'job_level'} <=> $NGS_batch_jobs{$b}->{'job_level'}) or ($a cmp $b)} keys %NGS_batch_jobs;

$cmd = `mkdir -p $log_dir`       unless (-e $log_dir);
$cmd = `mkdir -p $sh_dir`        unless (-e $sh_dir);
$cmd = `mkdir -p $sh_bundle_dir` unless (-e $sh_bundle_dir);
open(LOG, ">> $log_file") || die "can not write to $log_file";

######## parse NGS_samples
my %NGS_sample_data = ();
my @NGS_samples = ();
if (defined($sample_in)) {
  open(TMP, $sample_in) || die "can not open $sample_in";
  while($ll=<TMP>){
    next if ($ll =~ /^#/);
    next unless ($ll =~ /^\w/); chop($ll);
    my ($id, @data) = split(/\s+/,$ll);
    push(@NGS_samples, $id);
    $NGS_sample_data{$id} = [@data];
    if (not (-e $id)) { $cmd = `mkdir $id`;}
  }
  close(TMP);
}
elsif (defined($sample_command_in)) {
  my @lls = split(/,/, $sample_command_in); 
  foreach $ll (@lls) {
    my ($id, @data) = split(/:/, $ll);
    push(@NGS_samples, $id);
    $NGS_sample_data{$id} = [@data];
    if (not (-e $id)) { $cmd = `mkdir $id`;}
  }
}
else {
  die "no input samples";
}

my %CMD_opts = (); 
if (-e $opt_file) {
  ##format example
  ##CMDOPT JobID_A:opt0:opt1:opt2
  ##CMDOPT JobID_B:opt0:opt1
  ##CMDOPT JobID_C:opt0:opt1:opt2:opt3
  open(TMP, $opt_file) || die "can not open $opt_file";
  while($ll = <TMP>){
    next if ($ll =~ /^#/);
    next unless ($ll =~ /^CMDOPT/);
    chop($ll);
    my ($i, $opt1) = split(/\s+/, $ll);
    my ($job_id, @opts) = split(/:/, $opt1);
    $CMD_opts{$job_id} = [@opts];
  }
  close(TMP);
}
elsif ($opt_command_in) {
  my @lls = split(/,/, $opt_command_in);
  foreach $ll (@lls) {
    my ($job_id, @opts) = split(/:/, $ll);
    $CMD_opts{$job_id} = [@opts];
  }
}

########## processing subset of jobs
if ($subset_wfs) {
  my @wfs = split(/,/, $subset_wfs);
  $subset_flag = 1;
  foreach $i (@wfs) {
    my @jobs = @{ $NGS_batch_sets{$i}->{"jobs"} };
    foreach $j (@jobs) { $subset_jobs{$j} = 1; }
  }
}
if ($subset_jobs) {
  $subset_flag = 1;
  my @jobs = split(/,/, $subset_jobs);
  foreach $j (@jobs) { $subset_jobs{$j} = 1; }
  add_subset_jobs_by_dependency(); 
}
if ($subset_flag) {
  my $job_str = join(" ", keys %subset_jobs);
  write_log("Running subset of jobs: $job_str");
}

my $verify_flag = 0;
foreach $t_job_id (keys %NGS_batch_jobs) {
  if ($subset_flag) {next unless ($subset_jobs{$t_job_id});}
  $t_job = $NGS_batch_jobs{$t_job_id};
  $t_execution = $NGS_executions{ $t_job->{"execution"} };

  my $pe_parameter = ""; #### setup pe parameters
  if ($t_execution->{'type'} eq "qsub-pe") {
    my $t_cores_per_cmd  = $t_job->{"cores_per_cmd"};
       $t_cores_per_cmd  = 1 unless ($t_cores_per_cmd);
    $pe_parameter = "#\$ -pe orte $t_cores_per_cmd";
  }

  if ($t_job->{"cores_per_cmd"} > $t_execution->{"cores_per_node"} ) {
    $verify_flag = 1;
    write_log("$t_job_id needs $t_job->{\"cores_per_cmd\"} cores, but $t_job->{\"execution\"} only has $t_execution->{\"cores_per_node\"} cores");
  }

  my $cmds_per_node = POSIX::floor( $t_execution->{"cores_per_node"}  / $t_job->{"cores_per_cmd"});
  my $nodes_total  = POSIX::ceil($t_job->{"no_parallel"} / $cmds_per_node);
  $t_job->{"cmds_per_node"} = $cmds_per_node;
  $t_job->{"nodes_total"}   = $nodes_total;

  if ($t_job->{"nodes_total"} > $t_execution->{"number_nodes"}) {
    $verify_flag = 1;
    write_log("$t_job_id needs $t_job->{\"nodes_total\"} nodes, but $t_job->{\"execution\"} only has $t_execution->{\"number_nodes\"} nodes");
  }

  my @CMD_opts = ();
     @CMD_opts = @{$t_job->{CMD_opts}}   if (defined($t_job->{CMD_opts}  ));
     @CMD_opts = @{$CMD_opts{$t_job_id}} if (defined($CMD_opts{$t_job_id})); #### command line take over default

  foreach $t_sample_id (@NGS_samples) {
    my @t_commands = split(/\t/, $t_job->{"command"});
    my $t_command  = "";
    foreach my $c0 (@t_commands) {
      my $c1 = $c0; 
      $c1 =~ s/\\SAMPLE/$t_sample_id/g;
      $c1 =~ s/\\SELF/$t_job_id/g;
      # take it easy, assuming maxium 20 input files
      $c1 =~ s/\\INFILES\.0/$t_job->{"infiles"}->[0]/g;    $c1 =~ s/\\INFILES\.10/$t_job->{"infiles"}->[10]/g;
      $c1 =~ s/\\INFILES\.1/$t_job->{"infiles"}->[1]/g;    $c1 =~ s/\\INFILES\.11/$t_job->{"infiles"}->[11]/g;
      $c1 =~ s/\\INFILES\.2/$t_job->{"infiles"}->[2]/g;    $c1 =~ s/\\INFILES\.12/$t_job->{"infiles"}->[12]/g;
      $c1 =~ s/\\INFILES\.3/$t_job->{"infiles"}->[3]/g;    $c1 =~ s/\\INFILES\.13/$t_job->{"infiles"}->[13]/g;
      $c1 =~ s/\\INFILES\.4/$t_job->{"infiles"}->[4]/g;    $c1 =~ s/\\INFILES\.14/$t_job->{"infiles"}->[14]/g;
      $c1 =~ s/\\INFILES\.5/$t_job->{"infiles"}->[5]/g;    $c1 =~ s/\\INFILES\.15/$t_job->{"infiles"}->[15]/g;
      $c1 =~ s/\\INFILES\.6/$t_job->{"infiles"}->[6]/g;    $c1 =~ s/\\INFILES\.16/$t_job->{"infiles"}->[16]/g;
      $c1 =~ s/\\INFILES\.7/$t_job->{"infiles"}->[7]/g;    $c1 =~ s/\\INFILES\.17/$t_job->{"infiles"}->[17]/g;
      $c1 =~ s/\\INFILES\.8/$t_job->{"infiles"}->[8]/g;    $c1 =~ s/\\INFILES\.18/$t_job->{"infiles"}->[18]/g;
      $c1 =~ s/\\INFILES\.9/$t_job->{"infiles"}->[9]/g;    $c1 =~ s/\\INFILES\.19/$t_job->{"infiles"}->[19]/g;

      $c1 =~ s/\\DATA\.0/$NGS_sample_data{$t_sample_id}->[0]/g;    $c1 =~ s/\\DATA\.10/$NGS_sample_data{$t_sample_id}->[10]/g;
      $c1 =~ s/\\DATA\.1/$NGS_sample_data{$t_sample_id}->[1]/g;    $c1 =~ s/\\DATA\.11/$NGS_sample_data{$t_sample_id}->[11]/g;
      $c1 =~ s/\\DATA\.2/$NGS_sample_data{$t_sample_id}->[2]/g;    $c1 =~ s/\\DATA\.12/$NGS_sample_data{$t_sample_id}->[12]/g;
      $c1 =~ s/\\DATA\.3/$NGS_sample_data{$t_sample_id}->[3]/g;    $c1 =~ s/\\DATA\.13/$NGS_sample_data{$t_sample_id}->[13]/g;
      $c1 =~ s/\\DATA\.4/$NGS_sample_data{$t_sample_id}->[4]/g;    $c1 =~ s/\\DATA\.14/$NGS_sample_data{$t_sample_id}->[14]/g;
      $c1 =~ s/\\DATA\.5/$NGS_sample_data{$t_sample_id}->[5]/g;    $c1 =~ s/\\DATA\.15/$NGS_sample_data{$t_sample_id}->[15]/g;
      $c1 =~ s/\\DATA\.6/$NGS_sample_data{$t_sample_id}->[6]/g;    $c1 =~ s/\\DATA\.16/$NGS_sample_data{$t_sample_id}->[16]/g;
      $c1 =~ s/\\DATA\.7/$NGS_sample_data{$t_sample_id}->[7]/g;    $c1 =~ s/\\DATA\.17/$NGS_sample_data{$t_sample_id}->[17]/g;
      $c1 =~ s/\\DATA\.8/$NGS_sample_data{$t_sample_id}->[8]/g;    $c1 =~ s/\\DATA\.18/$NGS_sample_data{$t_sample_id}->[18]/g;
      $c1 =~ s/\\DATA\.9/$NGS_sample_data{$t_sample_id}->[9]/g;    $c1 =~ s/\\DATA\.19/$NGS_sample_data{$t_sample_id}->[19]/g;

      $c1 =~ s/\\INJOBS\.0/$t_job->{"injobs"}->[0]/g;      $c1 =~ s/\\INJOBS\.10/$t_job->{"injobs"}->[10]/g;
      $c1 =~ s/\\INJOBS\.1/$t_job->{"injobs"}->[1]/g;      $c1 =~ s/\\INJOBS\.11/$t_job->{"injobs"}->[11]/g;
      $c1 =~ s/\\INJOBS\.2/$t_job->{"injobs"}->[2]/g;      $c1 =~ s/\\INJOBS\.12/$t_job->{"injobs"}->[12]/g;
      $c1 =~ s/\\INJOBS\.3/$t_job->{"injobs"}->[3]/g;      $c1 =~ s/\\INJOBS\.13/$t_job->{"injobs"}->[13]/g;
      $c1 =~ s/\\INJOBS\.4/$t_job->{"injobs"}->[4]/g;      $c1 =~ s/\\INJOBS\.14/$t_job->{"injobs"}->[14]/g;
      $c1 =~ s/\\INJOBS\.5/$t_job->{"injobs"}->[5]/g;      $c1 =~ s/\\INJOBS\.15/$t_job->{"injobs"}->[15]/g;
      $c1 =~ s/\\INJOBS\.6/$t_job->{"injobs"}->[6]/g;      $c1 =~ s/\\INJOBS\.16/$t_job->{"injobs"}->[16]/g;
      $c1 =~ s/\\INJOBS\.7/$t_job->{"injobs"}->[7]/g;      $c1 =~ s/\\INJOBS\.17/$t_job->{"injobs"}->[17]/g;
      $c1 =~ s/\\INJOBS\.8/$t_job->{"injobs"}->[8]/g;      $c1 =~ s/\\INJOBS\.18/$t_job->{"injobs"}->[18]/g;
      $c1 =~ s/\\INJOBS\.9/$t_job->{"injobs"}->[9]/g;      $c1 =~ s/\\INJOBS\.19/$t_job->{"injobs"}->[19]/g;

      $c1 =~ s/\\CMDOPTS\.0/$CMD_opts[0]/g;    $c1 =~ s/\\CMDOPTS\.10/$CMD_opts[10]/g;
      $c1 =~ s/\\CMDOPTS\.1/$CMD_opts[1]/g;    $c1 =~ s/\\CMDOPTS\.11/$CMD_opts[11]/g;
      $c1 =~ s/\\CMDOPTS\.2/$CMD_opts[2]/g;    $c1 =~ s/\\CMDOPTS\.12/$CMD_opts[12]/g;
      $c1 =~ s/\\CMDOPTS\.3/$CMD_opts[3]/g;    $c1 =~ s/\\CMDOPTS\.13/$CMD_opts[13]/g;
      $c1 =~ s/\\CMDOPTS\.4/$CMD_opts[4]/g;    $c1 =~ s/\\CMDOPTS\.14/$CMD_opts[14]/g;
      $c1 =~ s/\\CMDOPTS\.5/$CMD_opts[5]/g;    $c1 =~ s/\\CMDOPTS\.15/$CMD_opts[15]/g;
      $c1 =~ s/\\CMDOPTS\.6/$CMD_opts[6]/g;    $c1 =~ s/\\CMDOPTS\.16/$CMD_opts[16]/g;
      $c1 =~ s/\\CMDOPTS\.7/$CMD_opts[7]/g;    $c1 =~ s/\\CMDOPTS\.17/$CMD_opts[17]/g;
      $c1 =~ s/\\CMDOPTS\.8/$CMD_opts[8]/g;    $c1 =~ s/\\CMDOPTS\.18/$CMD_opts[18]/g;
      $c1 =~ s/\\CMDOPTS\.9/$CMD_opts[9]/g;    $c1 =~ s/\\CMDOPTS\.19/$CMD_opts[19]/g;
      $t_command .= "$c1\n";
    }


    my @t_infiles = map { "$t_sample_id/$_" } @{$t_job->{"infiles"}};
    my @t_injobs  =                           @{$t_job->{"injobs"}};
    my $t_sh_file = "$sh_dir/$t_job_id.$t_sample_id.sh";
    my $f_start    = "$pwd/$t_sample_id/$t_job_id/WF.start.date";
    my $f_complete = "$pwd/$t_sample_id/$t_job_id/WF.complete.date";
    my $f_cpu      = "$pwd/$t_sample_id/$t_job_id/WF.cpu";
    $job_list{$t_job_id}{$t_sample_id} = {
      'sample_id'    => $t_sample_id,
      'job_id'       => $t_job_id,
      'status'       => 'wait',       #### status can be wait (input not ready), ready (input ready), submitted (submitted or running), completed
      'command'      => $t_command,
      'sh_file'      => $t_sh_file,
      'infiles'      => [@t_infiles],
      'injobs'       => [@t_injobs],
      'start_file'   => $f_start,
      'complete_file'=> $f_complete,
      'cpu_file'     => $f_cpu,
    };

    my $v_command = "";
    foreach my $vf (@{$t_job->{"non_zero_files"}}) {
      $v_command .= "if ! [ -s $t_job_id/$vf ]; then echo \"zero size $t_job_id/$vf\"; exit; fi\n";
    }


    if (not -e $t_sh_file) {
      write_log("Write sh file to $t_sh_file");
      open(TSH, "> $t_sh_file") || die "can not write to $t_sh_file\n";
      print TSH <<EOD;
$t_execution->{"template"}
$pe_parameter

my_host=`hostname`
my_pid=\$\$
my_core=$t_job->{"cores_per_cmd"}
my_queue=$t_job->{"execution"}
my_time_start=`date +%s`;

cd $pwd
cd $t_sample_id
mkdir $t_job_id
if ! [ -f $f_start ]; then date +\%s > $f_start;  fi
$t_command
$v_command
date +\%s > $f_complete
#times >> $f_cpu

my_time_end=`date +%s`;
my_time_spent=\$((my_time_end-my_time_start))
echo "sample=$t_sample_id job=$t_job_id host=\$my_host pid=\$my_pid queue=\$my_queue cores=\$my_core time_start=\$my_time_start time_end=\$my_time_end time_spent=\$my_time_spent" >> $f_cpu

EOD
      close(TSH);
      #validate_cmd_line($t_command, $t_sh_file, $t_sample_id);
    }
  } ########## foreach my $c0 (@t_commands)
} ########## foreach $t_job (keys %NGS_batch_jobs)

die if ($verify_flag);

if    ($this_task eq  "log-cpu"     ) { task_log_cpu();      exit 0;}
elsif ($this_task eq  "list-jobs"   ) { task_list_jobs();    exit 0;}
elsif ($this_task eq  "snapshot"    ) { task_snapshot();     exit 0;}
elsif ($this_task eq  "delete-jobs" ) { task_delete_jobs($second_opt);  exit 0;}
elsif ($this_task eq  "write-sh"    ) {                                 exit 0;}
elsif ($this_task                   ) { die "undefined task $this_task";}

################################################################################################
#  _____               _   _  _____  _____  _           _       _           _       _         
# |  __ \             | \ | |/ ____|/ ____|| |         | |     | |         (_)     | |        
# | |__) |   _ _ __   |  \| | |  __| (___  | |__   __ _| |_ ___| |__        _  ___ | |__  ___ 
# |  _  / | | | '_ \  | . ` | | |_ |\___ \ | '_ \ / _` | __/ __| '_ \      | |/ _ \| '_ \/ __|
# | | \ \ |_| | | | | | |\  | |__| |____) || |_) | (_| | || (__| | | |     | | (_) | |_) \__ \
# |_|  \_\__,_|_| |_| |_| \_|\_____|_____/ |_.__/ \__,_|\__\___|_| |_|     | |\___/|_.__/|___/
#                                      ______                      ______ _/ |                
#                                     |______|                    |______|__/                 
########## Run NGS_batch_jobs for each samples http://patorjk.com/software/taag
################################################################################################


my %execution_submitted = (); # number of submitted jobs (qsub) or threads (local sh)
my $sleep_time     = $sleep_time_min;
while(1) {
  my $flag_job_done = 1;

  ########## reset execution_submitted to 0
  foreach $i (keys %NGS_executions) { $execution_submitted{$i} = 0; }

  my $flag_qstat_xml_call = 0;
  foreach $t_job_id (keys %NGS_batch_jobs) {
    my $t_job = $NGS_batch_jobs{$t_job_id};
    my $t_execution = $NGS_executions{ $t_job->{"execution"} };
    my $exe_type = $t_execution->{type};
    $flag_qstat_xml_call = 1 if (($queue_system eq "SGE") and (($exe_type eq "qsub") or ($exe_type eq "qsub-pe")));
  }
  SGE_qstat_xml_query() if $flag_qstat_xml_call;

  ########## check and update job status for submitted jobs
  foreach $t_job_id (keys %NGS_batch_jobs) {
    if ($subset_flag) {next unless ($subset_jobs{$t_job_id});} 
    my $t_job = $NGS_batch_jobs{$t_job_id};
    foreach $t_sample_id (@NGS_samples) {
      my $t_sample_job = $job_list{$t_job_id}{$t_sample_id};
      my $status = $t_sample_job->{'status'};

      next if ($status eq "completed");
      ########## check file system to update job status
      ########## in case this is a restart run
      check_submitted_job($t_job_id, $t_sample_id);
      next if ($t_sample_job->{'status'} eq "completed");
      $flag_job_done = 0;
    }
  }

  if ($flag_job_done) { write_log("job completed!"); last; }

  ########## check and update job status based on dependance 
  foreach $t_job_id (keys %NGS_batch_jobs) {
    if ($subset_flag) {next unless ($subset_jobs{$t_job_id});} 
    my $t_job = $NGS_batch_jobs{$t_job_id};
    foreach $t_sample_id (@NGS_samples) {
      my $t_sample_job = $job_list{$t_job_id}{$t_sample_id};
      my $status = $t_sample_job->{'status'};

      next unless ($status eq "wait");
      my @t_infiles = @{ $t_sample_job->{'infiles'} };
      my @t_injobs  = @{ $t_sample_job->{'injobs'} };
      my $t_ready_flag = 1;

      foreach $i (@t_infiles) {
        next if (-s $i); ####  non-zero size file
        $t_ready_flag = 0;
        last;
      }

      foreach $i (@t_injobs) {
        next if ( $job_list{$i}{$t_sample_id}->{'status'} eq "completed"); #### injob completed
        $t_ready_flag = 0;
        last;
      }
      if ($t_ready_flag) {
        $t_sample_job->{"status"} = "ready";
        write_log("$t_job_id,$t_sample_id: change status to ready");
      }
    }
  }

  ########## submit local sh jobs
  my $has_submitted_some_jobs = 0;
  foreach $t_job_id (keys %NGS_batch_jobs) {
    if ($subset_flag) {next unless ($subset_jobs{$t_job_id});} 
    my $t_job = $NGS_batch_jobs{$t_job_id};
    my $t_execution = $NGS_executions{ $t_job->{"execution"} };
    my $t_execution_id = $t_job->{"execution"};

    if ($subset_flag) {next unless ($subset_jobs{$t_job_id});} 
    next unless ($t_execution->{'type'} eq "sh");
    next if ( $execution_submitted{$t_execution_id} >= $t_execution->{"cores_per_node"} ); #### all cores are used

    foreach $t_sample_id (@NGS_samples) {
      my $t_sample_job = $job_list{$t_job_id}{$t_sample_id};
      my $status = $t_sample_job->{'status'};
      next unless ($status eq "ready");
      next if ( ($execution_submitted{$t_execution_id} + $t_job->{"cores_per_cmd"} * $t_job->{"no_parallel"}) > $t_execution->{"cores_per_node"} ); #### no enough available cores
      #### now submitting 

      my $t_sh_file = $t_sample_job->{'sh_file'};
      my $t_sh_pid  = "$t_sh_file.pids";
      for ($i=0; $i<$t_job->{"no_parallel"}; $i++) {
        $cmd = `sh $t_sh_file >/dev/null 2>&1 &`;
      }
      $cmd = `touch $t_sh_pid`;
      $t_sample_job->{'status'} = "submitted";
      write_log("$t_job_id,$t_sample_id: change status to submitted");
      $execution_submitted{ $t_execution_id } += $t_job->{"cores_per_cmd"} * $t_job->{"no_parallel"}; 
      $has_submitted_some_jobs = 1;
    }
  }

  ########## submit qsub-pe jobs, multiple jobs may share same node
  foreach $t_job_id (keys %NGS_batch_jobs) {
    if ($subset_flag) {next unless ($subset_jobs{$t_job_id});} 
    my $t_job = $NGS_batch_jobs{$t_job_id};
    my $t_execution = $NGS_executions{ $t_job->{"execution"} };
    my $t_execution_id = $t_job->{"execution"};

    next unless ($t_execution->{'type'} eq "qsub-pe");
    next if ( $execution_submitted{$t_execution_id} >= $t_execution->{"number_nodes"} ); #### resource full

    my $t_cores_per_node = $t_execution->{"cores_per_node"};
    my $t_cores_per_cmd  = $t_job->{"cores_per_cmd"};
    my $t_cores_per_job  = $t_cores_per_cmd * $t_job->{"no_parallel"};
    my $t_nodes_per_job  = $t_cores_per_job / $t_cores_per_node;

    foreach $t_sample_id (@NGS_samples) {
      my $t_sample_job = $job_list{$t_job_id}{$t_sample_id};
      my $status = $t_sample_job->{'status'};
      next unless ($status eq "ready");

      my $t_sh_file = $t_sample_job->{'sh_file'};
      my $t_sh_pid  = "$t_sh_file.pids";
      open(TID, "> $t_sh_pid")    || die "can not write to $t_sh_pid";

      for ($i=0; $i<$t_job->{"no_parallel"}; $i++) {
        my $t_stderr = "$t_sh_file.$i.stderr";
        my $t_stdout = "$t_sh_file.$i.stdout";
        $cmd = `qsub $t_execution->{"command_name_opt"} $t_job_id $t_execution->{"command_err_opt"} $t_stderr $t_execution->{"command_out_opt"} $t_stdout $t_sh_file 2>$log_fileq`;
        my $qsub_id = 0;
        if ($cmd =~ /(\d+)/) { $qsub_id = $1;} else {die "can not submit qsub job and return a id\n";}
        print TID "$qsub_id\n"; #### $cmd returns qsub id, write these ids to pid file for future qstat 
        $execution_submitted{$t_execution_id} += $t_nodes_per_job;
        write_log("$t_sh_bundle submitted for sample $t_sample_id, qsubid $cmd");
      }

      close(TID);
      $has_submitted_some_jobs = 1;
      $t_sample_job->{'status'} = "submitted";
    }
  } ########## END foreach $t_job_id (keys %NGS_batch_jobs) 

  ########## submit qsub jobs
  foreach $t_job_id (keys %NGS_batch_jobs) {
    if ($subset_flag) {next unless ($subset_jobs{$t_job_id});} 
    my $t_job = $NGS_batch_jobs{$t_job_id};
    my $t_execution = $NGS_executions{ $t_job->{"execution"} };
    my $t_execution_id = $t_job->{"execution"};

    if ($subset_flag) {next unless ($subset_jobs{$t_job_id});} 
    next unless ($t_execution->{'type'} eq "qsub");
    next if ( $execution_submitted{$t_execution_id} >= $t_execution->{"number_nodes"} ); #### resource full

    my $t_cores_per_node = $t_execution->{"cores_per_node"};
    my $t_cores_per_cmd  = $t_job->{"cores_per_cmd"};
    my $t_cores_per_job  = $t_cores_per_cmd * $t_job->{"no_parallel"};
    my $t_nodes_per_job  = POSIX::ceil($t_cores_per_job / $t_cores_per_node);
    my $t_cmds_per_node  = int($t_cores_per_node / $t_cores_per_cmd);
    my $t_jobs_per_node  = int($t_cores_per_node / $t_cores_per_job);

    ########## 1. this loop process jobs need 1 or more nodes per sample, ie. bundle within a sample, e.g. blast against refseq
    foreach $t_sample_id (@NGS_samples) {
      my $t_sample_job = $job_list{$t_job_id}{$t_sample_id};
      my $status = $t_sample_job->{'status'};
      next unless ($status eq "ready");
      next unless ($t_jobs_per_node <= 1);                 #### unless need >= 1 node, including jobs use between (51%-100%) cores per node
      last if ( ($execution_submitted{$t_execution_id} + $t_nodes_per_job) > $t_execution->{"number_nodes"}); #### no enough available queues

      my $t_sh_file = $t_sample_job->{'sh_file'};
      my $t_sh_bundle = "$sh_bundle_dir/$t_job_id.$t_sample_id.$$.sh";
      my $t_stderr    = "$t_sh_bundle.stderr";
      my $t_stdout    = "$t_sh_bundle.stdout";
      my $t_sh_pid  = "$t_sh_file.pids";

      open(TID, "> $t_sh_pid")    || die "can not write to $t_sh_pid";
      open(BSH, "> $t_sh_bundle") || die "can not write to $t_sh_bundle";
      print BSH <<EOD;
$t_execution->{"template"}
cd $pwd
EOD
      for ($i=0; $i<$t_cmds_per_node; $i++) {
        print BSH "sh $t_sh_file &\n";
        print BSH "sleep 3\n";
      }
      print BSH "wait\n";
      close(BSH);

      for ($i=0; $i<$t_nodes_per_job; $i++) {
        $cmd = `qsub $t_execution->{"command_name_opt"} $t_job_id $t_execution->{"command_err_opt"} $t_stderr $t_execution->{"command_out_opt"} $t_stdout $t_sh_bundle 2>$log_fileq`;
        my $qsub_id = 0;
        if ($cmd =~ /(\d+)/) { $qsub_id = $1;} else {die "can not submit qsub job and return a id\n";}
        print TID "$qsub_id\n"; #### $cmd returns qsub id, write these ids to pid file for future qstat 
        $execution_submitted{$t_execution_id}++;
        write_log("$t_sh_bundle submitted for sample $t_sample_id, qsubid $cmd");
      }
      close(TID);
      $has_submitted_some_jobs = 1;
      $t_sample_job->{'status'} = "submitted";
    } ########## END foreach $t_sample_id (@NGS_samples) 


    ########## 2. this loop process jobs need less than 1 node per sample, ie. bundle jobs across samples, e.g. qc 
    my @t_bundle = ();
    my $available_nodes = $t_execution->{"number_nodes"} - $execution_submitted{$t_execution_id};
    my $no_sample_can_be_processed = $available_nodes * $t_jobs_per_node;
    my @t_samples = ();
    my $t_batch_no = 0;

    foreach $t_sample_id (@NGS_samples) { #### same loop as next, to find out @t_samples and last sample can run
      my $t_sample_job = $job_list{$t_job_id}{$t_sample_id};
      my $status = $t_sample_job->{'status'};
      next unless ($status eq "ready");
      next unless ($t_jobs_per_node > 1);              #### unless a node can host 2 or more jobs
      last if ( $t_execution->{"number_nodes"} - $execution_submitted{$t_execution_id} <=0);
      push(@t_samples, $t_sample_id);
    }
    my $last_sample_can_run = $t_samples[-1];
    @t_samples = ();

    foreach $t_sample_id (@NGS_samples) {
      my $t_sample_job = $job_list{$t_job_id}{$t_sample_id};
      my $status = $t_sample_job->{'status'};
      next unless ($status eq "ready");
      next unless ($t_jobs_per_node > 1);              #### unless a node can host 2 or more jobs
      last if ( $t_execution->{"number_nodes"} - $execution_submitted{$t_execution_id} <=0);
      push(@t_samples, $t_sample_id);

      #### bundle @t_samples to one qsub job
      if ((($#t_samples+1) == $t_jobs_per_node) or ($t_sample_id eq $last_sample_can_run)) {
        my $t_sh_bundle = "$sh_bundle_dir/$t_job_id.samples-$t_batch_no.$$.sh";
        my $t_stderr    = "$t_sh_bundle.stderr";
        my $t_stdout    = "$t_sh_bundle.stdout";

        open(BSH, "> $t_sh_bundle") || die "can not write to $t_sh_bundle";
      print BSH <<EOD;
$t_execution->{"template"}
cd $pwd
EOD
        foreach $i (@t_samples) {
          my $t_sh_file = $job_list{$t_job_id}{$i}->{'sh_file'};
          for ($j=0; $j<$t_job->{"no_parallel"}; $j++) {
            print BSH "sh $t_sh_file &\n";
            print BSH "sleep 3\n";
          }
        }
        print BSH "wait\n";
        close(BSH);

        $cmd = `qsub $t_execution->{"command_name_opt"} $t_job_id $t_execution->{"command_err_opt"} $t_stderr $t_execution->{"command_out_opt"} $t_stdout $t_sh_bundle 2>$log_fileq`;
        my $qsub_id = 0;
        if ($cmd =~ /(\d+)/) { $qsub_id = $1;} else {die "can not submit qsub job and return a id\n";}

        foreach $i (@t_samples) {
          my $t_sh_file = $job_list{$t_job_id}{$i}->{'sh_file'};
          my $t_sh_pid  = "$t_sh_file.pids";
          open(TID, "> $t_sh_pid")    || die "can not write to $t_sh_pid";
          print TID "$qsub_id\n"; #### $cmd returns qsub id, write these ids to pid file for future qstat 
          write_log("$t_sh_bundle submitted for sample $i, qsubid $cmd");
          close(TID);
          $job_list{$t_job_id}{$i}->{'status'} = "submitted";
        }

        $has_submitted_some_jobs = 1;
        $execution_submitted{$t_execution_id}++;
        @t_samples = (); #### clear
        $t_batch_no++;
      }
    } ########## END foreach $t_sample_id (@NGS_samples)
  } ########## END foreach $t_job_id (keys %NGS_batch_jobs) 


  #### if has submitted some jobs, reset waiting time, otherwise double waiting time
  print_job_status_summary();
  if ($has_submitted_some_jobs) {
    $sleep_time = $sleep_time_min;
  }
  else {
    $sleep_time = $sleep_time*2;
    $sleep_time = $sleep_time_max if ($sleep_time > $sleep_time_max);
  }
  write_log("sleep $sleep_time seconds");
  sleep($sleep_time);
} ########## END  while(1) 

task_log_cpu();
################################################################################
########## END Run NGS_batch_jobs for each samples
################################################################################

close(LOG);
##########


sub write_log {
  my @txt = @_;
  my $i;
  my $date = `date`; chop($date);
  foreach $i (@txt) {
    print LOG    "$date $i\n";
    print STDERR "$date $i\n";
  }
  print LOG    "\n";
  print STDERR "\n";
}
########## END write_log

sub SGE_qstat_xml_query {
  my ($i, $j, $k, $cmd, $ll);
  %qstat_xml_data = (); #### global
  $cmd = `qstat -f -xml`;
  if ($cmd =~ /<queue_info/) { #### dummy 
    $qstat_xml_data{"NULL"}= ["NULL","NULL"];
  }

  my @lls = split(/\n/, $cmd);
  $i = 2; #### skip first 2 lines
  for (;     $i<$#lls+1; $i++) {
    if ($lls[$i] =~ /<job_list/) {
      my ($id, $name, $state);
      for (; $i<$#lls+1; $i++) {
        last if ($lls[$i] =~ /<\/job_list/);
        if ($lls[$i] =~ /<JB_job_number>(\d+)/) {  $id = $1;}
        if ($lls[$i] =~ /<JB_name>([^<]+)/) { $name = $1;}
        if ($lls[$i] =~ /<state>([^<]+)/) {$state = $1;}
      }
      if (defined($id) and defined($name) and defined($state)) {
        $qstat_xml_data{$id} = [$name, $state];
      }
    }
  }
}

########## check submitted job by checking pids, or qsub ids
########## update job status from wait|ready -> submitted if pid file exit (in case of restart of this script)
########## update job status from wait|ready|submitted -> completed if sh calls or qsub calls finished
##########    these pids or qsub ids are done
sub check_submitted_job {
  my ($t_job_id, $t_sample_id) = @_;
  my $t_sample_job = $job_list{$t_job_id}{$t_sample_id};
  my $t_job = $NGS_batch_jobs{$t_job_id};
  my $t_execution = $NGS_executions{ $t_job->{"execution"} };

  my ($i, $j, $k, $flag, $ll, $cmd);

  my $t_sh_file = $t_sample_job->{'sh_file'};
  my $t_sh_pid  = "$t_sh_file.pids";

  # status won't change unless there is a pid file
  return unless (-e $t_sh_pid);

  my $status = $t_sample_job->{'status'};
  if (($status eq "wait") or ($status eq "ready")) {
    $t_sample_job->{'status'} = "submitted";
    write_log("$t_job_id,$t_sample_id: change status to submitted");
  }

  my $exe_type = $t_execution->{type};

  if ($exe_type eq "sh") {
    $cmd = `ps -ef | grep "$t_sh_file" | grep -v grep`;
    if ($cmd =~ /\w/) { # still running 
      $execution_submitted{ $t_job->{"execution"} } += $t_job->{"cores_per_cmd"} * $t_job->{"no_parallel"};
    }
    elsif (validate_job_files($t_job_id, $t_sample_id)) {
      $t_sample_job->{'status'} = "completed";
      write_log("$t_job_id,$t_sample_id: change status to completed");
    }
    else {
      $t_sample_job->{'status'} = "error";
      write_log("$t_job_id,$t_sample_id: change status to error");
    }
    return;
  }
  elsif (($exe_type eq "qsub") or ($exe_type eq "qsub-pe")) {
    my @pids = ();
    open(CHECK, $t_sh_pid) || die "Can not open $t_sh_pid\n";
    while($ll = <CHECK>) {
      chop($ll); next unless ($ll =~ /\w/);
      push(@pids, $ll);
    }
    close(CHECK);

    my $finish_flag = 1;
    foreach $i (@pids) {
      if (($queue_system eq "SGE") and %qstat_xml_data) {
        if (defined($qstat_xml_data{$i})) {
          $t_sample_job->{'status'} = "running" if (($qstat_xml_data{$i}->[1] eq "r") and ($t_sample_job->{'status'} eq "submitted"));
          $finish_flag = 0;
          $execution_submitted{ $t_job->{"execution"} } ++; 
        }
      }
      elsif ($queue_system eq "SGE") {
        $cmd = `qstat -j $i | grep job_number`;
        if ($cmd =~ /$i/) {
          $finish_flag = 0;
          $execution_submitted{ $t_job->{"execution"} } ++; 
        }
      }
      else {
        $cmd = `qstat -r $i | grep $i`;
        $j = (split(/\D/,$cmd))[0];
        if ($j == $i) { # this job is running
          $finish_flag = 0;
          $execution_submitted{ $t_job->{"execution"} } ++; 
        }
      }
    }
    if ($finish_flag == 1) {
      if (validate_job_files($t_job_id, $t_sample_id)) {
        $t_sample_job->{'status'} = "completed";
        write_log("$t_job_id,$t_sample_id: change status to completed");
      }
      else {
        $t_sample_job->{'status'} = "error";
        write_log("$t_job_id,$t_sample_id: change status to error");
      }
    }
    return;
  }
  else {
    die "unknown execution type: $exe_type\n";
  }
}
########## END sub check_submitted_job 


# WF.start.date and WF.complete.date need to have non-zero size
sub validate_job_files {
  my ($t_job_id, $t_sample_id) = @_;
  my ($i, $j, $k);
  my $t_sample_job = $job_list{$t_job_id}{$t_sample_id};

  return 0 unless (-s $t_sample_job->{'start_file'} );
  return 0 unless (-s $t_sample_job->{'complete_file'} );
  return 0 unless (-s $t_sample_job->{'cpu_file'} );

  return 1; #### pass
}
########## END validate_job_files


sub print_job_status_summary {
  my ($t_job_id, $t_sample_id);
  my ($i, $j, $k);

  my %job_status = ();
  my $job_total = 0;
  foreach $t_job_id (keys %NGS_batch_jobs) {
    if ($subset_flag) {next unless ($subset_jobs{$t_job_id});} 
    foreach $t_sample_id (@NGS_samples) {
      my $t_sample_job = $job_list{$t_job_id}{$t_sample_id};
      my $status = $t_sample_job->{'status'};
      $job_status{$status}++;
      $job_total++;
    }
  }

  print STDERR "total jobs: $job_total,";
  foreach $i (sort keys %job_status) {
    print STDERR "$i: $job_status{$i},";
  } 
  print STDERR "\n"; 
}
########## END print_job_status_summary


sub validate_cmd_line {
  my ($i, $j, $k);
  my ($t_command, $t_sh_file, $t_sample_id) = @_;
  my @cmds = split(/\n/,$t_command);

  my @warn_path = ();
  foreach $i (@cmds) {
    my ($key_cmd, @opts) = split(/\s+/, $i);
    if ($key_cmd =~ /\//) {
      if (not -e $key_cmd) { push(@warn_path, $key_cmd); } 
    }
    @opts = grep {/\//} @opts;
    foreach $j (@opts) {
      my @opts1 = split(/,|;|>|<|\|/,$j);
      foreach $k (@opts1) {
        $k = "$t_sample_id/$k" unless (($k =~ /^\//) or ($k =~ /^\./));
        if (not -e $k) { push(@warn_path, $k); }
      }
    }
  }

  if (@warn_path) {
    print STDERR "File or program doesn't exist in $t_sh_file: ", join(" ", @warn_path), "\n";
  }

}
########## END validate_cmd_line

sub add_subset_jobs_by_dependency {
  my ($i, $j, $k, $ll, $t_job_id, $t_sample_id, $t_job);

  while(1) {
    my $num_subset_jobs = scalar keys %subset_jobs;

    foreach $t_job_id (keys %subset_jobs) {
      $t_job = $NGS_batch_jobs{$t_job_id};
      my @t_injobs  = @{$t_job->{"injobs"}};

      for $j (@t_injobs) {
        $subset_jobs{$j} = 1;
      }
    }

    last if ($num_subset_jobs == scalar keys %subset_jobs);
  }
}
########## END add_subset_jobs_by_dependency


sub task_level_jobs {
  my ($i, $j, $k, $ll, $t_job_id, $t_sample_id, $t_job);
  my %job_level = ();

  while(1) {
    my $change_flag = 0;

    foreach $t_job_id (keys %NGS_batch_jobs) {
      $t_job = $NGS_batch_jobs{$t_job_id};
      my @t_injobs  = @{$t_job->{"injobs"}};

      if (@t_injobs) {
        my $max_level_injob;
        foreach $j (@t_injobs) {
          next unless defined ($job_level{$j});
          $max_level_injob = $job_level{$j} if ($job_level{$j} > $max_level_injob);          
        }

        next unless (defined($max_level_injob));
        $max_level_injob++; #### one more level 
        if (not defined ($job_level{$t_job_id})) {
          $job_level{$t_job_id}=$max_level_injob;
          $change_flag = 1;
        }
        elsif ($max_level_injob > $job_level{$t_job_id}) {
          $job_level{$t_job_id}=$max_level_injob;
          $change_flag = 1;
        }
      }
      else {
        if (not defined ($job_level{$t_job_id})) {
          $job_level{$t_job_id}=1;
          $change_flag = 1;
        }
      }
    }
    last unless ($change_flag);
  }

  foreach $t_job_id (sort keys %NGS_batch_jobs) {
    $NGS_batch_jobs{$t_job_id}->{"job_level"} = $job_level{$t_job_id};
  }
}
########## END task_list_jobs

sub task_snapshot {
  my ($t_job_id, $t_sample_id);
  my ($i, $j, $k);

  if ($this_task) {
    my $flag_qstat_xml_call = 0;
    foreach $t_job_id (keys %NGS_batch_jobs) {
      my $t_job = $NGS_batch_jobs{$t_job_id};
      my $t_execution = $NGS_executions{ $t_job->{"execution"} };
      my $exe_type = $t_execution->{type};
      $flag_qstat_xml_call = 1 if (($queue_system eq "SGE") and (($exe_type eq "qsub") or ($exe_type eq "qsub-pe")));
    }
    SGE_qstat_xml_query() if $flag_qstat_xml_call;

    foreach $t_sample_id (@NGS_samples) {
      foreach $t_job_id (keys %NGS_batch_jobs) {
        check_submitted_job($t_job_id, $t_sample_id);
      }
    }
  }

  my $max_len_sample = 0;
  foreach $t_sample_id (@NGS_samples) {
    $max_len_sample = length($t_sample_id) if (length($t_sample_id) > $max_len_sample);
  }
  my $max_len_job = 0;
  foreach $t_job_id (@NGS_batch_jobs) {
    $max_len_job = length($t_job_id) if (length($t_job_id) > $max_len_job);
  }

  print <<EOD;
Job status: 
.\twait
-\tsubmitted
r\trunning  
+\tcompleted
!\terror
EOD

  for ($i=$max_len_job-1; $i>=0; $i--) {
    print ' 'x$max_len_sample, "\t";
    foreach $t_job_id (@NGS_batch_jobs) {
      print " ", ($i<length($t_job_id) ? substr(reverse($t_job_id), $i, 1):" ");
    }
    print "\n";
  }

  foreach $t_sample_id (@NGS_samples) {
    print "$t_sample_id\t";
    foreach $t_job_id (@NGS_batch_jobs) {
      my $t_sample_job = $job_list{$t_job_id}{$t_sample_id};
      my $status = $t_sample_job->{'status'};
      if    ($status eq "completed") { print " +";}
      elsif ($status eq "submitted") { print " -";}
      elsif ($status eq "running"  ) { print " r";}
      elsif ($status eq "wait"     ) { print " .";}
      elsif ($status eq "error"    ) { print " !";}
      else                           { print " _";}
    }
    print "\n";
  }
}
########## END task_snapshot

sub task_list_jobs {
  my ($i, $j, $k, $ll, $t_job_id, $t_sample_id, $t_job);
  foreach $t_job_id (@NGS_batch_jobs) {
    $t_job = $NGS_batch_jobs{$t_job_id};
    #my @t_infiles = @{$t_job->{"infiles"}};
    my @t_injobs  = @{$t_job->{"injobs"}};

    #print "\tInput_files:", join(",", @t_infiles) if @t_infiles;
    print "$t_job_id\tIn_jobs:[" , join(",", @t_injobs), "]\tJob_level:$t_job->{'job_level'}\n";
  }
}
########## END task_list_jobs

sub file1_after_file2 {
  my ($file1, $file2) = @_;

  # if not exist file1, assume it is in future, so it is newer
  if (not -e ($file1)) {return 0;}
  if (not -e ($file2)) {return 0;}

  my $mtime1 = (stat($file1))[9];
  my $mtime2 = (stat($file2))[9];

  return ( ($mtime1 > $mtime2) ? 1 : 0);
}
######## END file1_after_file2

sub file1_same_or_after_file2 {
  my ($file1, $file2) = @_;

  # if not exist file1, assume it is in future, so it is newer
  if (not -e ($file1)) {return 0;}
  if (not -e ($file2)) {return 0;}

  my $mtime1 = (stat($file1))[9];
  my $mtime2 = (stat($file2))[9];

  return ( ($mtime1 >= $mtime2) ? 1 : 0);
}
######## END file1_after_file2


sub task_delete_jobs {
  my $opt = shift;
  my ($i, $j, $k, $ll, $t_job_id, $t_sample_id);
  my ($mode, $c) = split(/:/, $opt);
  my $tmp_sh = "NGS-$$.sh";

  open(TMPSH, "> $tmp_sh") || die "can not write to file $tmp_sh";
  print TMPSH "#Please execute the following commands\n";
  foreach $t_sample_id (@NGS_samples) {
    my %job_to_delete_ids = ();
    if ($mode eq "jobids") {
       %job_to_delete_ids = map {$_, 1} split(/,/,$c);
    }
    elsif ($mode eq "run_after") {
      die "file $c doesn't exist!" unless (-e $c);
      foreach $t_job_id (keys %NGS_batch_jobs) {
        my $t_sample_job = $job_list{$t_job_id}{$t_sample_id};
        my $t_sh_file = $t_sample_job->{'sh_file'};
        my $t_sh_pid  = "$t_sh_file.pids";
        next unless (-e $t_sh_pid);   #### unless the job is submitted
        #$job_to_delete_ids{$t_job_id} = 1 if (file1_same_or_after_file2( $t_sample_job->{'start_file'} , $c));
        $job_to_delete_ids{$t_job_id} = 1 if (file1_same_or_after_file2( $t_sh_pid , $c));

      }
    }
    else {
      die "unknown option for deleting jobs: $opt";
    }

    # now %job_to_delete_ids are jobs need to be deleted
    # next find all jobs that depends on them, recrusively
    my $no_jobs_to_delete = scalar keys %job_to_delete_ids;
    while(1) {
      foreach $t_job_id (keys %NGS_batch_jobs) {
        my $t_sample_job = $job_list{$t_job_id}{$t_sample_id};
        my $t_sh_file = $t_sample_job->{'sh_file'};
        my $t_sh_pid  = "$t_sh_file.pids";
        next unless (-e $t_sh_pid);   #### unless the job is submitted
        my @t_injobs  = @{ $t_sample_job->{'injobs'} };
        foreach my $t_job_id_2 (@t_injobs) {
          $job_to_delete_ids{$t_job_id} = 1 if ($job_to_delete_ids{$t_job_id_2});
        }
      }
      last if ($no_jobs_to_delete == (scalar keys %job_to_delete_ids)); #### no more depending jobs
      $no_jobs_to_delete = scalar keys %job_to_delete_ids;
    }

    if ($no_jobs_to_delete) {
      print TMPSH "#jobs to be deleted for $t_sample_id: ", join(",", keys %job_to_delete_ids), "\n";
      print       "#jobs to be deleted for $t_sample_id: ", join(",", keys %job_to_delete_ids), "\n";
      foreach $t_job_id (keys %job_to_delete_ids) {
        my $t_sample_job = $job_list{$t_job_id}{$t_sample_id};
        my $t_sh_file = $t_sample_job->{'sh_file'};
        my $t_sh_pid  = "$t_sh_file.pids";
        print TMPSH "\\rm -rf $pwd/$t_sample_id/$t_job_id\n";
        print TMPSH "\\rm $t_sh_pid\n";        
        print TMPSH "\\rm $t_sh_file.*.std*\n";

        #### find the qsub ids to be deleted 
        my $qids = `cat $t_sh_pid`; $qids =~ s/\n/ /g; $qids =~ s/\s+/ /g;
        print TMPSH "qdel $qids\n";
      }
    }
  }
  close(TMPSH);
  print "The script is not delete the file, please run $tmp_sh to delete files!!!\n\n";
}
########## END task_list_jobs

sub task_log_cpu {
  my ($i, $j, $k, $ll, $t_job_id, $t_sample_id);

  my %cpu_info;
  foreach $t_job_id (keys %NGS_batch_jobs) {
    if ($subset_flag) {next unless ($subset_jobs{$t_job_id});} 
    my $t_job = $NGS_batch_jobs{$t_job_id};
    foreach $t_sample_id (@NGS_samples) {

      $cpu_info{$t_job_id}{$t_sample_id} = [$t_wall, $t_cpu];
    }
  }

  foreach $t_sample_id (@NGS_samples) {
    my $f_cpu = "$pwd/$t_sample_id/WF.cpu";
    open(CPUOUT, "> $f_cpu") || die "Can not open $f_cpu";
    print CPUOUT "#job_name\tCores\tWall(s)\tWall_time\tCPU(s)\tCPU_time\n";
    my $min_start = 1402092131 * 999999;
    my $max_end   = 0;
    my $sum_cpu   = 0;
    foreach $t_job_id (keys %NGS_batch_jobs) {
      if ($subset_flag) {next unless ($subset_jobs{$t_job_id});} 
      my $t_job = $NGS_batch_jobs{$t_job_id};
      my $t_core     = $t_job->{"cores_per_cmd"} * $t_job->{"no_parallel"};

      my $t_sample_job = $job_list{$t_job_id}{$t_sample_id};
      my $f_start    = $t_sample_job->{'start_file'};
      my $f_complete = $t_sample_job->{'complete_file'};
      my $f_cpu      = $t_sample_job->{'cpu_file'};
      my $t_start    = `cat $f_start`;    $t_start =~ s/\s//g; $min_start = $t_start if ($t_start < $min_start);
      my $t_end      = `cat $f_complete`; $t_end   =~ s/\s//g; $max_end   = $t_end   if ($t_end   > $max_end);
      my $t_wall     = int($t_end - $t_start);
         $t_wall     = 0 unless ($t_wall>0);

      my $t_cpu = 0;
      if (open(TCPU, $f_cpu)) {
        while($ll = <TCPU>) {
          chop($ll);
          if ($ll =~ /^(\d+)m(\d+)/) {
            $t_cpu += $1 * 60;
          }
        }
        close(TCPU);
      }
      $sum_cpu += $t_cpu;

      my $t_walls = time_str1($t_wall);
      my $t_cpus  = time_str1($t_cpu);
      print CPUOUT "$t_job_id\t$t_core\t$t_wall\t$t_walls\t$t_cpu\t$t_cpus\n";
    }
    my $t_wall = ($max_end - $min_start); $t_wall     = 0 unless ($t_wall>0);
    my $t_walls = time_str1($t_wall);
    my $sum_cpus= time_str1($sum_cpu);
    print CPUOUT "total\t-\t$t_wall\t$t_walls\t$sum_cpu\t$sum_cpus\n";
    close(CPUOUT);
  }
}
######### END task_log_cpu

sub time_str1 {
  my $s = shift;
  my $str = "";

  $str .= int($s/3600); $str .= "h"; $s = $s % 3600;
  $str .= int($s/60);   $str .= "m"; $s = $s % 60;
  $str .= $s;           $str .= "s";

  return $str;
}
########## END time_str1;






sub usage {
<<EOD;

# =============================== NG-Omics-WF ==================================
#  _   _  _____         ____            _              __          ________ 
# | \\ | |/ ____|       / __ \\          (_)             \\ \\        / /  ____|
# |  \\| | |  __ ______| |  | |_ __ ___  _  ___ ___ _____\\ \\  /\\  / /| |__   
# | . ` | | |_ |______| |  | | '_ ` _ \\| |/ __/ __|______\\ \\/  \\/ / |  __|  
# | |\\  | |__| |      | |__| | | | | | | | (__\\__ \\       \\  /\\  /  | |     
# |_| \\_|\\_____|       \\____/|_| |_| |_|_|\\___|___/        \\/  \\/   |_|     
#                                                                           
# =========================== Next Generation Omics data workflow tools ========

To run workflow: 
    $0 -s sample_file -i workflow_file

Options:

    -i workflow configration file, required

    -s sample data file, required unless -S is present
       File format example
#Sample data file example, TAB or space delimited for following lines
Sample_ID1 sample_data_0 sample_data_1
Sample_ID2 sample_data_0 sample_data_1
Sample_ID3 sample_data_0 sample_data_1

    -S sample data from command line, required unless -s is present
       format: Sample_ID1:sample_data_0:sample_data_0:sample_data_1,Sample_ID2:sample_data_0:sample_data_1

    -j run sub sets of jobs, optional, the workflow will run all jobs by default
       e.g. -j qc or -j qc,fastqc

    -t parameter file, optional, replace default paramters in workflow configration file
       File format example
#parameter file example, TAB or space delimited for following lines
CMDOPT JobID_A:opt0:opt1:opt2
CMDOPT JobID_B:opt0:opt1

    -T parameter from command line
       format:  JobID_A:opt0:opt1:opt2,JobID_B:opt0:opt1

    -r root directory of NGS-tools

    -J optional tasks
        write-sh: write sh files and quite
        log-cpu: gathering cpu time for each run for each sample
        list-jobs: list jobs
        snapshot: snapshot current job status
        delete-jobs: delete jobs, must supply jobs delete syntax by option -Z
          e.g. -J delete-jobs -Z jobids:assembly,blast  ---delete assembly,blast and all jobs depends on them
               -J delete-jobs -Z run_after:filename     ---delete jobs that has start time (WF.start.date) after this file, and all depending jobs

    -Z secondary parameter used by other options, such as -J

    -Q queue system, default SGE
        can be PBS, SGE

Question and comments:
     http://weizhongli-lab.org/ngomicswf           liwz\@sdsc.edu

EOD
}



############################################################################################
# _______    ________  _________       ___________________   ________  .____       _________
# \      \  /  _____/ /   _____/       \__    ___/\_____  \  \_____  \ |    |     /   _____/
# /   |   \/   \  ___ \_____  \   ______ |    |    /   |   \  /   |   \|    |     \_____  \ 
#/    |    \    \_\  \/        \ /_____/ |    |   /    |    \/    |    \    |___  /        \
#\____|__  /\______  /_______  /         |____|   \_______  /\_______  /_______ \/_______  /
#        \/        \/        \/                           \/         \/        \/        \/ 
############################################################################################

