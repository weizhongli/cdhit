#!/usr/bin/env python
# =============================== NG-Omics-WF ==================================
#  _   _  _____         ____            _              __          ________ 
# | \ | |/ ____|       / __ \          (_)             \ \        / /  ____|
# |  \| | |  __ ______| |  | |_ __ ___  _  ___ ___ _____\ \  /\  / /| |__   
# | . ` | | |_ |______| |  | | '_ ` _ \| |/ __/ __|______\ \/  \/ / |  __|  
# | |\  | |__| |      | |__| | | | | | | | (__\__ \       \  /\  /  | |     
# |_| \_|\_____|       \____/|_| |_| |_|_|\___|___/        \/  \/   |_|     
#                                                                           
# Workflow tools for next generation genomics, metagenomics, RNA-seq 
# and other type of omics data analyiss, 
#
# Software originally developed since 2010 by Weizhong Li at UCSD
#                                               currently at JCVI
#
# https://github.com/weizhongli/ngomicswf       liwz@sdsc.edu
# ==============================================================================

import os
import sys
import re
import argparse 
from argparse import RawTextHelpFormatter
import math
import subprocess
import time
import logging
import textwrap
import imp
import collections
import xml.etree.ElementTree as ET

__author__ = 'Weizhong Li'

############## Global variables
banner = '''
            ==================================================================
            Workflow tools for next generation genomics, metagenomics, RNA-seq
            and other type of omics data analyiss,
        
            Software originally developed since 2010 by Weizhong Li at UCSD
                                                          currently at JCVI
        
            http://weizhongli-lab.org/ngomicswf           liwz@sdsc.edu
            ==================================================================
'''
NGS_config = None
NGS_samples = []
NGS_sample_data = {}
NGS_opts = {}
pwd = os.path.abspath('.')
subset_flag = False
subset_jobs = []
qstat_xml_data = collections.defaultdict(dict)
job_list = collections.defaultdict(dict)  # as job_list[$t_job_id][$t_sample_id] = {}
execution_submitted = {}                  # number of submitted jobs (qsub) or threads (local sh)
############## END Global variables


def fatal_error(message, exit_code=1):
  print message
  exit(exit_code)


def read_parameters(args):
  '''read option parameters from file or command line'''
  if args.parameter_file:
    try:
      ##format example
      ##JobID_A opt0 opt1 opt2
      ##JobID_B opt0 opt1
      ##JobID_C opt0 opt1 opt2 opt3
      f = open(args.parameter_file, 'r')
      for line in f:
        if line[0] == '#':
          continue
        if not re.match('^\w', line):
          continue
        ll = re.split('\s+', line.rstrip());
        NGS_opts[ ll[0] ] = ll[1:]
      f.close()
    except IOError:
      fatal_error('cannot open ' + args.parameter_file, exit_code=1)

  elif args.parameter_name:
    for line in re.split(',', args.parameter_name):
      ll = re.split(':', line);
      NGS_opts[ ll[0] ] = ll[1:]
  return


def read_samples(args):
  '''read sample and sample data from file or command line'''
  if args.sample_file:
    try:
      f = open(args.sample_file, 'r')
      for line in f:
        if line[0] == '#':
          continue
        if not re.match('^\w', line):
          continue
        ll = re.split('\s+', line.rstrip());
        NGS_samples.append(ll[0]);
        NGS_sample_data[ ll[0] ] = ll[1:]
      f.close()
    except IOError:
      fatal_error('cannot open ' + args.sample_file, exit_code=1)

  elif args.sample_name:
    for line in re.split(',', args.sample_name):
      ll = re.split(':', line);
      NGS_samples.append(ll[0]);
      NGS_sample_data[ ll[0] ] = ll[1:]
  else:
    fatal_error('no input sample', exit_code=1)

  for sample in NGS_samples:
    if os.path.exists(sample):
      if os.path.isdir(sample):
        pass
      else:
        fatal_error('file exist: ' + sample, exit_code=1)
    else:
      if os.system("mkdir " + sample):
        fatal_error('can not mkdir: ' + sample, exit_code=1)
  return


def task_level_jobs(NGS_config):
  '''according to dependancy, make level of jobs'''
  job_level = {}
  while True:
    change_flag = False
    for t_job_id in NGS_config.NGS_batch_jobs.keys():
      t_job = NGS_config.NGS_batch_jobs[t_job_id]
      t_injobs = []
      if 'injobs' in t_job.keys():
        t_injobs = t_job['injobs']

      if len(t_injobs) > 0:
        max_level_injob = 0
        for j in t_injobs:
          if not j in job_level.keys():
            continue
            if job_level[j] > max_level_injob:
              max_level_injob = job_level[j]
        if max_level_injob == 1:
          continue
        max_level_injob +=1  #### one more level 
        if (t_job_id in job_level.keys()) and (job_level[t_job_id] >= max_level_injob):
          continue
        job_level[t_job_id]=max_level_injob
        change_flag = 1
      else:
        if not t_job_id in job_level.keys():
          job_level[t_job_id]=1
          change_flag = True

    if not change_flag: break

  # print job_level
  for t_job_id in NGS_config.NGS_batch_jobs.keys():
    NGS_config.NGS_batch_jobs[t_job_id]['job_level'] = job_level[t_job_id]
      
  return
#### END task_level_jobs(NGS_config)


def add_subset_jobs_by_dependency(NGS_config):
  '''add dependant jobs'''
  while True:
    num_subset_jobs = len(subset_jobs)
    for t_job_id in subset_jobs:
      t_job = NGS_config.NGS_batch_jobs[t_job_id]
      if 'injobs' in t_job.keys():
        for j in t_job['injobs']:
          if not (j in subset_jobs):
            subset_jobs.append(j)
    if num_subset_jobs == len(subset_jobs): break
  return
#### END add_subset_jobs_by_dependency()

       
def make_job_list(NGS_config):
  '''make sh script for each job / sample'''

  verify_flag = False
  for t_job_id in NGS_config.NGS_batch_jobs:
    if subset_flag and not (t_job_id in subset_jobs):
      continue

    t_job = NGS_config.NGS_batch_jobs[ t_job_id ]
    t_execution = NGS_config.NGS_executions[ t_job["execution"] ]

    pe_parameter = ''
    if t_execution[ 'type' ] == 'qsub-pe':
      t_cores_per_cmd  = t_job[ 'cores_per_cmd' ]
      pe_parameter = "#$ -pe orte " + str(t_cores_per_cmd)
      if 'pe_para' in t_execution.keys():
        pe_parameter = "#$ " + t_execution[ 'pe_para' ] + " " +  str(t_cores_per_cmd)

    if t_job[ 'cores_per_cmd' ] > t_execution[ 'cores_per_node' ]:
      fatal_error('not enough cores ' + t_job_id, exit_code=1)

    t_job[ 'cmds_per_node' ] = t_execution[ 'cores_per_node' ] / t_job[ 'cores_per_cmd' ]
    t_job[ 'nodes_total' ] = math.ceil( t_job[ 'no_parallel' ] / float(t_job[ 'cmds_per_node' ]))
 
    if t_job[ 'nodes_total' ] > t_execution[ 'number_nodes' ]:
      fatal_error('not enough nodes ' + t_job, exit_code=1)

    CMD_opts = []
    if 'CMD_opts' in t_job.keys():  
      CMD_opts = t_job[ 'CMD_opts' ]
    if t_job_id in NGS_opts.keys():
      CMD_opts = NGS_opts[ t_job_id ]

    for t_sample_id in NGS_samples:
      t_command = t_job[ 'command' ]
      t_command = re.sub('\$SAMPLE', t_sample_id, t_command)
      t_command = re.sub('\$SELF'  , t_job_id, t_command)

      for i in NGS_config.ENV.keys():
        t_command = re.sub('\$ENV.'+i, NGS_config.ENV[i], t_command)

      for i_data in range(0, len(NGS_sample_data[ t_sample_id ])):
        t_data = NGS_sample_data[ t_sample_id ][i_data]
        t_re = '\$DATA\.' + str(i_data)
        t_command = re.sub(t_re, t_data, t_command)

      t_injobs = []
      if 'injobs' in t_job.keys():
        t_injobs = t_job[ 'injobs' ]
        for i_data in range(0, len(t_job[ 'injobs' ])):
          t_data = t_job[ 'injobs' ][i_data]
          t_re = '\$INJOBS\.' + str(i_data)
          t_command = re.sub(t_re, t_data, t_command)

      for i_data in range(0, len(CMD_opts)):
        t_data = CMD_opts[i_data]
        t_re = '\$CMDOPTS\.' + str(i_data)
        t_command = re.sub(t_re, t_data, t_command)

      v_command = ''
      if 'non_zero_files' in t_job.keys():
        for t_data in t_job[ 'non_zero_files' ]:
          v_command = v_command + \
            'if ! [ -s {0}/{1} ]; then echo "zero size {2}/{3}"; exit; fi\n'.format(t_job_id, t_data, t_job_id, t_data)

      f_start    = pwd + '/' + t_sample_id + '/' + t_job_id + '/WF.start.date'
      f_complete = pwd + '/' + t_sample_id + '/' + t_job_id + '/WF.complete.date'
      f_cpu      = pwd + '/' + t_sample_id + '/' + t_job_id + '/WF.cpu'
      t_sh_file  = '{0}/WF-sh/{1}.{2}.sh'.format(pwd, t_job_id, t_sample_id)
      t_infiles = []
      if 'infiles' in t_job.keys():
        t_infiles = map(lambda x: t_sample_id + "/" + x, t_job[ 'infiles' ])
      job_list[ t_job_id ][ t_sample_id ] = {
        'sample_id'    : t_sample_id,
        'job_id'       : t_job_id,
        'status'       : 'wait',       #### status can be wait (input not ready), ready (input ready), submitted (submitted or running), completed
        'command'      : t_command,
        'sh_file'      : t_sh_file, 
        'infiles'      : t_infiles,
        'injobs'       : t_injobs,
        'start_file'   : f_start,
        'complete_file': f_complete,
        'cpu_file'     : f_cpu }

      if not os.path.exists( t_sh_file ):
        try:
          tsh = open(t_sh_file, 'w')
          tsh.write('''{0}
{1}

my_host=`hostname`
my_pid=$$
my_core={2}
my_queue={3}
my_time_start=`date +%s`

cd {4}/{5}
mkdir {6}
if ! [ -f {7} ]; then date +%s > {7};  fi
{8}
{9}
date +%s > {10}
my_time_end=`date +%s`;
my_time_spent=$((my_time_end-my_time_start))
echo "sample={5} job={6} host=$my_host pid=$my_pid queue=$my_queue cores=$my_core time_start=$my_time_start time_end=$my_time_end time_spent=$my_time_spent" >> {11}

'''.format(t_execution['template'], pe_parameter, t_job['cores_per_cmd'], t_job['execution'], pwd, t_sample_id, t_job_id, f_start, t_command, v_command, f_complete, f_cpu ))
          tsh.close()
          os.system('chmod u+x '+ t_sh_file)
        except IOError:
          fatal_error('cannot write to ' + job_list[ 't_job_id' ][ 't_sample_id' ][ 'sh_file' ], exit_code=1)
  return
### END def make_job_list(NGS_config):


def time_str1(s):
  str1 = str(s/3600) + 'h'
  s = s % 3600
  str1 = str1 + str(s/60) + 'm'
  s = s % 60
  str1 = str1 + str(s) + 's'
  return str1


def task_list_jobs(NGS_config):
  for t_job_id in NGS_config.NGS_batch_jobs.keys():
    t_job = NGS_config.NGS_batch_jobs[t_job_id]

    t_injobs = []
    if 'injobs' in t_job.keys():
      t_injobs  = t_job['injobs']
    print '{0}\tIn_jobs:[ {1} ]\tJob_level:{2}\n'.format(t_job_id, ','.join(t_injobs), t_job['job_level'] )


def task_snapshot(NGS_config):
  '''print job status'''
  queue_system = NGS_config.queue_system   #### default "SGE"
  this_task = True
  if this_task:
    flag_qstat_xml_call = False
    for t_job_id in NGS_config.NGS_batch_jobs.keys():
      t_job = NGS_config.NGS_batch_jobs[t_job_id]
      t_execution = NGS_config.NGS_executions[ t_job['execution']]
      exe_type = t_execution['type']
      if (queue_system == 'SGE') and (exe_type in ['qsub','qsub-pe']):
        flag_qstat_xml_call = True

    if flag_qstat_xml_call:
      SGE_qstat_xml_query()

    for t_sample_id in NGS_samples:
      for t_job_id in NGS_config.NGS_batch_jobs.keys():
        if subset_flag:
          if not (t_job_id in subset_jobs):
            continue
        check_submitted_job(NGS_config, t_job_id, t_sample_id)

  max_len_sample = 0;
  for t_sample_id in NGS_samples:
    if len(t_sample_id) > max_len_sample:
      max_len_sample = len(t_sample_id)
  max_len_job = 0;
  for t_job_id in NGS_config.NGS_batch_jobs.keys():
    if len(t_job_id) > max_len_job:
      max_len_job = len(t_job_id)

  print '''
Job status: 
.\twait
-\tsubmitted
r\trunning  
+\tcompleted
!\terror
x\tunselected job
'''

  for i1 in range(max_len_job):
    i = max_len_job - i1 - 1
    print ' ' * max_len_sample + "\t", 
    for t_job_id in NGS_config.NGS_batch_jobs.keys():
      if i < len(t_job_id):
        print ' ' + t_job_id[-i-1],
      else:
        print '  ',
    print ''
  print 'Sample\t' + ('-' * 30)

  for t_sample_id in NGS_samples:
    print t_sample_id + '\t',
    for t_job_id in NGS_config.NGS_batch_jobs.keys():
      if subset_flag:
        if not (t_job_id in subset_jobs):
          print ' x',
          continue
      t_sample_job = job_list[t_job_id][t_sample_id]
      status = t_sample_job['status']
      if   status == 'completed': print ' +',
      elif status == 'submitted': print ' -',
      elif status == 'running':   print ' r',
      elif status == 'wait':      print ' .',
      elif status == 'error':     print ' !',
      else:                       print ' _',
    print ''

  print '\n\n'
  return
### def task_snapshot():


def task_delete_jobs(NGS_config, opt):
  '''manually delete jobs and its dependant jobs'''

  mode, c = re.split(':', opt)
  tmp_sh = 'NGS-{0}.sh'.format(os.getpid())

  try:
    tsh = open(tmp_sh, 'w')
  except IOError:
    fatal_error('cannot write to ' + tmp_sh, exit_code=1)
  
  tsh.write('#Please execute the following commands\n')

  for t_sample_id in NGS_samples:
    job_to_delete_ids =[]
    if mode == 'jobids':
      job_to_delete_ids = re.split(',', c)
    elif mode == 'run after':
      if not os.path.exists(c):
        fatel_error('File does not exist:' + c, exit_code=1)
      for t_job_id in NGS_config.NGS_batch_jobs.keys():
        if subset_flag:
          if not (t_job_id in subset_jobs): continue
        t_sample_job = job_list[t_job_id][t_sample_id]

        t_sh_pid  = t_sample_job['sh_file'] + '.pids'
        if not os.path.exists(t_sh_pid): continue
        if file1_same_or_after_file2(t_sh_pid, c):
          job_to_delete_ids.append(t_job_id)
    else:
      fatel_error('unknown option for deleteing jobs: ' + opt, exit_code=1)

    # now job_to_delete_ids are jobs need to be deleted
    # next find all jobs that depends on them, recrusively
    no_jobs_to_delete = len(job_to_delete_ids)
    while True:
      for t_job_id in NGS_config.NGS_batch_jobs.keys():
        if subset_flag:
          if not (t_job_id in subset_jobs): continue
        t_sample_job = job_list[t_job_id][t_sample_id]

        t_sh_pid  = t_sample_job['sh_file'] + '.pids'
        if not os.path.exists(t_sh_pid): continue

        for t_job_id_2 in t_sample_job['injobs']:
          if t_job_id_2 in job_to_delete_ids:
            if not t_job_id in job_to_delete_ids:
              job_to_delete_ids.append(t_job_id)

      if no_jobs_to_delete == len(job_to_delete_ids): break
      no_jobs_to_delete = len(job_to_delete_ids)

    if not no_jobs_to_delete: continue
    tsh.write('#jobs to be deleted for ' + t_sample_id + ': ' + ' '.join(job_to_delete_ids) + '\n'), 

    for t_job_id in job_to_delete_ids:
      t_sample_job = job_list[t_job_id][t_sample_id]
      t_sh_pid  = t_sample_job['sh_file'] + '.pids'
      tsh.write('\\rm -rf {0}/{1}\n'.format(t_sample_id, t_job_id))
      tsh.write('\\rm '+ t_sh_pid + '\n')

      #### find the qsub ids to be deleted 
      pids = [] #### either pids, or qsub ids
      try:
        f = open(t_sh_pid, 'r')
        pids = f.readlines()
        f.close()
        pids = [x.strip() for x in pids]
      except IOError:
        fatal_error('cannot open ' + t_sh_pid, exit_code=1)
      tsh.write('qdel ' + ' '.join([str(x) for x in pids]) + '\n')

    tsh.write('\n\n')
  tsh.close()
  print 'The script does not delete the files, please run ' + tmp_sh + ' to delete files!!!\n\n';
  return
#### END def task_delete_jobs()


def file1_same_or_after_file2(file1, file2):
  # if not exist file1, assume it is in future, so it is newer
  if not os.path.exists(file1): return True
  # otherwise file1 exist 
  # if file2 not exist
  if not os.path.exists(file2): return False
  if os.path.getmtime(file1) >= os.path.getmtime(file2): return True
  else:                                                  return False


def SGE_qstat_xml_query():
  '''run qstat -f -xml and get xml tree'''
  global qstat_xml_data
  qstat_xml_data = collections.defaultdict(dict)
  t_out = ''
  try:
    t_out  = subprocess.check_output(['qstat -f -xml'], shell=True)
  except:
    fatal_error("can not run qstat", exit_code=1)

  qstat_xml_root = ET.fromstring(t_out)
  #qstat_xml_root = qstat_xml.getroot()
  for job_list in qstat_xml_root.iter('job_list'):
    job_id    = job_list.find('JB_job_number').text
    job_name  = job_list.find('JB_name').text
    job_state = job_list.find('state').text
    qstat_xml_data[job_id] = [job_name, job_state]

  return
#### END def SGE_qstat_xml_query()


def print_job_status_summary(NGS_config):
  '''print jobs status'''
  job_status = {}
  job_total = 0

  for t_job_id in NGS_config.NGS_batch_jobs.keys():
    if subset_flag:
      if not (t_job_id in subset_jobs):
        continue
    for t_sample_id in NGS_samples:
      status = job_list[t_job_id][t_sample_id]['status']
      if status in job_status.keys():
        job_status[status] +=1
      else:
        job_status[status] =1
      job_total +=1

  print 'total jobs: {0},'.format(job_total) , 
  for i in job_status.keys():
    print '{0}: {1}, '.format(i, job_status[i]), 
  print '\n'


def run_workflow(NGS_config):
  '''major loop for workflow run'''
  queue_system = NGS_config.queue_system   #### default "SGE"
  sleep_time_min = 15
  sleep_time_max = 120
  sleep_time = sleep_time_min

  while 1:
    flag_job_done = True
    ########## reset execution_submitted to 0
    for i in NGS_config.NGS_executions.keys():
      execution_submitted[ i ] = False

    flag_qstat_xml_call = False
    for t_job_id in NGS_config.NGS_batch_jobs.keys():
      t_job = NGS_config.NGS_batch_jobs[t_job_id]
      t_execution = NGS_config.NGS_executions[ t_job['execution']]
      exe_type = t_execution['type']
      if (queue_system == 'SGE') and (exe_type in ['qsub','qsub-pe']):
        flag_qstat_xml_call = True

    if flag_qstat_xml_call:
      SGE_qstat_xml_query()

    ########## check and update job status for submitted jobs
    for t_job_id in NGS_config.NGS_batch_jobs.keys():
      t_job = NGS_config.NGS_batch_jobs[t_job_id]
      if subset_flag:
        if not (t_job_id in subset_jobs):
          continue
      for t_sample_id in NGS_samples:
        t_sample_job = job_list[t_job_id][t_sample_id]
        if t_sample_job['status'] == 'completed':
          continue

        check_submitted_job(NGS_config, t_job_id, t_sample_id)
        if t_sample_job['status'] == 'completed':
          continue
        flag_job_done = False

    if flag_job_done:
      print_job_status_summary(NGS_config)
      print 'job completed!'
      break

    ########## check and update job status based on dependance 
    for t_job_id in NGS_config.NGS_batch_jobs.keys():
      t_job = NGS_config.NGS_batch_jobs[t_job_id]
      if subset_flag:
        if not (t_job_id in subset_jobs):
          continue
      for t_sample_id in NGS_samples:
        t_sample_job = job_list[t_job_id][t_sample_id]
        if t_sample_job['status'] != 'wait':
          continue

        t_ready_flag = True
        for i in t_sample_job['infiles']:
          if os.path.exists(i) and os.path.getsize(i) > 0:
            continue
          t_ready_flag = False
          break

        for i in t_sample_job['injobs']:
          if job_list[i][t_sample_id]['status'] == 'completed':
            continue
          t_ready_flag = False
          break

        if t_ready_flag:
          t_sample_job['status'] = 'ready'
          print '{0},{1}: change status to ready\n'.format(t_job_id, t_sample_id)
    ########## END check and update job status based on dependance 

    ########## submit local sh jobs
    has_submitted_some_jobs = False
    for t_job_id in NGS_config.NGS_batch_jobs.keys():
      t_job = NGS_config.NGS_batch_jobs[t_job_id]
      if subset_flag:
        if not (t_job_id in subset_jobs):
          continue
      t_execution = NGS_config.NGS_executions[ t_job['execution']]
      t_execution_id = t_job['execution']
      if t_execution['type'] != 'sh': 
        continue
      if execution_submitted[t_execution_id] >= t_execution['cores_per_node']:
        continue
      for t_sample_id in NGS_samples:
        t_sample_job = job_list[t_job_id][t_sample_id]
        if t_sample_job['status'] != 'ready':
          continue
        if (execution_submitted[t_execution_id] + t_job['cores_per_cmd'] * t_job['no_parallel']) > \
            t_execution['cores_per_node']: #### no enough available cores
          continue

        #### now submitting 
        pid_file = open( t_sample_job['sh_file'] + '.pids', 'w')
        for i in range(0, t_job['no_parallel']):
          p = subprocess.Popen(['/bin/bash', t_sample_job['sh_file']], shell=True)
          pid_file.write(str(p.pid)+'\n')
        pid_file.close()
        t_sample_job['status'] = 'submitted'
        print '{0},{1}: change status to submitted\n'.format(t_job_id, t_sample_id)
        execution_submitted[ t_execution_id ] += t_job['cores_per_cmd'] * t_job['no_parallel'] 
        has_submitted_some_jobs = True
    ########## END submit local sh jobs

    ########## submit qsub-pe jobs, multiple jobs may share same node
    for t_job_id in NGS_config.NGS_batch_jobs.keys():
      t_job = NGS_config.NGS_batch_jobs[t_job_id]
      if subset_flag:
        if not (t_job_id in subset_jobs):
          continue
      t_execution = NGS_config.NGS_executions[ t_job['execution']]
      t_execution_id = t_job['execution']

      if t_execution['type'] != 'qsub-pe':
        continue
      if execution_submitted[t_execution_id] >= t_execution['number_nodes']:
        continue
      t_cores_per_node = t_execution['cores_per_node']
      t_cores_per_cmd  = t_job['cores_per_cmd']
      t_cores_per_job  = t_cores_per_cmd * t_job['no_parallel']
      t_nodes_per_job  = t_cores_per_job / t_cores_per_node

      for t_sample_id in NGS_samples:
        t_sample_job = job_list[t_job_id][t_sample_id]
        if t_sample_job['status'] != 'ready':
          continue

        #### now submitting 
        pid_file = open( t_sample_job['sh_file'] + '.pids', 'w')
        for i in range(0, t_job['no_parallel']):
          t_stderr = t_sample_job['sh_file'] + '.' + str(i) + '.stderr'
          t_stdout = t_sample_job['sh_file'] + '.' + str(i) + '.stdout'

          qsub_exe = 'qsub'
          if 'qsub_exe' in t_execution.keys(): qsub_exe = t_execution['qsub_exe']
          command_line = qsub_exe + ' {0} {1} {2} {3} {4} {5} {6}'.format(t_execution['command_name_opt'], t_job_id,
                                                                          t_execution['command_err_opt'], t_stderr, 
                                                                          t_execution['command_out_opt'], t_stdout, t_sample_job['sh_file'])
          cmd = subprocess.check_output([command_line], shell=True)
          if re.search('\d+', cmd):
            pid = re.search('\d+', cmd).group(0)
            pid_file.write(pid + '\n')
          else:
            fatal_error('error submitting jobs')
          execution_submitted[t_execution_id] += t_nodes_per_job
          print '{0} submitted for {1}\n'.format(t_sample_job['sh_file'], t_sample_id)

        pid_file.close()
        t_sample_job['status'] = 'submitted'
        has_submitted_some_jobs = True
    ########## END submit qsub-pe jobs, multiple jobs may share same node
   
    ########## submit qsub jobs, job bundles disabled here, if need, check the original Perl script

    #### if has submitted some jobs, reset waiting time, otherwise double waiting time
    print_job_status_summary(NGS_config)
    if has_submitted_some_jobs:
      sleep_time = sleep_time_min
    else:
      sleep_time  *= 2
      if sleep_time > sleep_time_max:
        sleep_time = sleep_time_max
    time.sleep(sleep_time);
  #### END while 1:
  return
#### END def run_workflow(NGS_config)


def check_pid(pid):        
  '''Check For the existence of a unix pid. '''
# opt 1, this doesn't print OSError to stderr
  return os.path.exists('/proc/' + str(pid))

# opt 2,  
#  try:
#    os.kill(pid, 0)
#  except OSError: return False
#  else:           return True


def check_any_pids(pids):
  '''Check For the existence of a list of unix pids. return True if any one exist'''
  for pid in pids:
    if check_pid(pid):
      return True
  return False


def check_any_qsub_pids(pids):
  '''Check For the existence of a list of qsub pids. return True if any one exist'''
  for pid in pids:
    if pid in qstat_xml_data.keys():
      return True
  return False


def validate_job_files(t_job_id, t_sample_id):
  '''return True if necessary file exist'''
  t_sample_job = job_list[t_job_id][t_sample_id]
  if not (os.path.exists(t_sample_job['start_file'])    and os.path.getsize(t_sample_job['start_file']) > 0):    return False
  if not (os.path.exists(t_sample_job['complete_file']) and os.path.getsize(t_sample_job['complete_file']) > 0): return False
  if not (os.path.exists(t_sample_job['cpu_file'])      and os.path.getsize(t_sample_job['cpu_file']) > 0):      return False
  return True


#### def check_submitted_job()
def check_submitted_job(NGS_config, t_job_id, t_sample_id):
  '''
  check submitted jobs by checking pids or qsub ids
  update job status from wait|ready -> submitted if pid file exit (in case of restart of this script)
  update job status from wait|ready|submitted -> completed if sh calls or qsub calls finished
  '''
  t_sample_job = job_list[t_job_id][t_sample_id]
  t_job = NGS_config.NGS_batch_jobs[t_job_id]
  t_execution = NGS_config.NGS_executions[ t_job['execution']]

  t_sh_pid = t_sample_job['sh_file'] + '.pids'
  if not os.path.exists(t_sh_pid): return

  status = t_sample_job['status']
  if ((status == 'wait') or (status == 'ready')):
    t_sample_job['status'] = 'submitted'
    print '{0},{1}: change status to submitted\n'.format(t_job_id, t_sample_id)

  pids = [] #### either pids, or qsub ids
  try:
    f = open(t_sh_pid, 'r')
    pids = f.readlines()
    f.close()
    pids = [x.strip() for x in pids]
  except IOError:
    fatal_error('cannot open ' + t_sh_pid, exit_code=1)

  if len(pids) == 0:
    fatal_error('empty file ' + t_sh_pid, exit_code=1)
  
  exe_type = t_execution['type']
  if (exe_type == 'sh'):
    if check_any_pids(pids):    #### still running
      execution_submitted[ t_job['execution'] ] += t_job['cores_per_cmd'] * t_job['no_parallel']
    elif validate_job_files(t_job_id, t_sample_id):                       #### job finished
      t_sample_job['status'] = 'completed'
      print '{0},{1}: change status to completed\n'.format(t_job_id, t_sample_id)
    else:
      t_sample_job['status'] = 'error'
      print '{0},{1}: change status to error\n'.format(t_job_id, t_sample_id)
    return

  elif ((exe_type == 'qsub') or (exe_type == 'qsub-pe')):
    if check_any_qsub_pids(pids):    #### still running
      pass
    elif validate_job_files(t_job_id, t_sample_id):                       #### job finished
      t_sample_job['status'] = 'completed'
      print '{0},{1}: change status to completed\n'.format(t_job_id, t_sample_id)
    else:
      t_sample_job['status'] = 'error'
      print '{0},{1}: change status to error\n'.format(t_job_id, t_sample_id)
  else:
    fatal_error('unknown execution type: '+ exe_type, exit_code=1)
  return
#### END def check_submitted_job()


if __name__ == "__main__":
  parser = argparse.ArgumentParser(formatter_class = RawTextHelpFormatter,
                                   description     = textwrap.dedent(banner))

  parser.add_argument('-i', '--input',       help="workflow configration file, required", required=True)
  parser.add_argument('-s', '--sample_file', help='''
sample data file, required unless -S is present
File format example:
#Sample data file example, TAB or space delimited for following lines
Sample_ID1 sample_data_0 sample_data_1
Sample_ID2 sample_data_0 sample_data_1
Sample_ID3 sample_data_0 sample_data_1
  ''')
  parser.add_argument('-S', '--sample_name', help='''
sample data from command line, required unless -s is present
format:
Sample_ID1:sample_data_0:sample_data_0:sample_data_1,Sample_ID2:sample_data_0:sample_data_1
  ''')
  parser.add_argument('-t', '--parameter_file', help='''
replace default paramters in workflow configration file
File format example:
#parameter file example, TAB or space delimited for following lines
JobID_A opt0 opt1 opt2
JobID_B opt0 opt1
  ''')
  parser.add_argument('-T', '--parameter_name', help='''
parameter from command line
format:
JobID_A:opt0:opt1:opt2,JobID_B:opt0:opt1
  ''')
  parser.add_argument('-j', '--jobs', help='''run sub set of jobs, optional
the workflow will run all jobs by default.
to run sub set of jobs: -j qc or -j qc,fastqc
  ''')
  parser.add_argument('-J', '--task', help='''optional tasks
write-sh: write sh files and quite
log-cpu: gathering cpu time for each run for each sample
list-jobs: list jobs
snapshot: snapshot current job status
delete-jobs: delete jobs, must supply jobs delete syntax by option -Z
  e.g. -J delete-jobs -Z jobids:assembly,blast  ---delete assembly,blast and all jobs depends on them
       -J delete-jobs -Z run_after:filename     ---delete jobs that has start time (WF.start.date) after this file, and all depending jobs
  ''')
  parser.add_argument('-Z', '--second_parameter', help='secondary parameter used by other options, such as -J')
  parser.add_argument('-Q', '--queye', help='queue system, e.g. PBS, SGE', default='SGE')

  args = parser.parse_args()

  if (args.sample_file is None) and (args.sample_name is None) :
    parser.error('No sample file or sample name')
  NGS_config = imp.load_source('NGS_config', args.input)

  print banner
  read_samples(args)
  read_parameters(args)

  if args.jobs:
    subset_flag = True
    subset_jobs = re.split(',', args.jobs)
    add_subset_jobs_by_dependency(NGS_config)
    print 'Running subset jobs:',
    print subset_jobs
    print '\n'

  if not os.path.exists('WF-sh'): os.system('mkdir WF-sh')

  task_level_jobs(NGS_config)
  make_job_list(NGS_config)

  if args.task:
    if args.task == 'list-jobs':
      task_list_jobs(NGS_config)
      exit(0)
    elif args.task == 'snapshot':
      task_snapshot(NGS_config)
      exit(0)
    elif args.task == 'delete-jobs':
      task_delete_jobs(NGS_config, args.second_parameter)
      exit(0)
    elif args.task == 'write-sh':
      exit(0)
    else:
      fatal_error('undefined task' + args.task, exit_code=1)

  run_workflow(NGS_config)
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
