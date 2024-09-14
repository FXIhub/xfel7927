# use active jobs on slurm to see if jobs are running
#   
#   <job-name>-<EXP_ID>-<run-number>
#   e.g. cxi-007076-35

# use log files at:
#   ${EXP_PREFIX}/scratch/log/<job-name>-<EXP_ID>-%A-<run-number>.out
#   e.g. /gpfs/exfel/exp/SPB/202421/p007076/scratch/log/cxi-007076-234234234-35.out
# to see if job completed successfully

import fnmatch
import subprocess
import glob
import os


from constants import PREFIX, EXP_ID

class SLURM_status():
    def __init__(self):
        self.update_slurm()
        
    def update_slurm(self):
        jobs = subprocess.run(['squeue', '--format="%.30j"'], stdout=subprocess.PIPE)
        jobs = jobs.stdout.decode('utf-8').split('\n')
        self.slurm_jobs = jobs

    def running(self, job_name, run):
        job_string = f'{job_name}-{EXP_ID}-{run}'
        
        # is job running?   
        s = f'*{job_string}*"'
        match = fnmatch.filter(self.slurm_jobs, s)
        running = len(match) > 0
        return running


class LOG_status():
    
    def __init__(self, log_dir, EXP_ID):

def get_job_status(job_name, slurm_jobs, run, EXP_ID, log_dir):
    # Check SLURM
    #------------
    job_string = f'{job_name}-{EXP_ID}-{run}'
    
    # is job running?   
    s = f'*{job_string}*"'
    match = fnmatch.filter(jobs, s)
    if len(match) > 0:
        running = True

    # Check logs
    #------------
    log_files  = glob.glob(f'{log_dir}/{job_name}-{EXP_ID}-*-{run}')
    log_files.sort(key = os.path.getmtime, reverse = True)
    
    is_log_file = len(log_file) > 0 
    log_success = False
    if is_log_file :
        log_file = log_file[0]
         
        # has the job finished successfully?
        # check if log file exist and has "{job_name} done"     
        cmd = subprocess.run(['grep', f'{job_name} done', log_file])
        log_success = len(cmd.stdout) > 0
    return running, is_log_file, log_success
        
    

