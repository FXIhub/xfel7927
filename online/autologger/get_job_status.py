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
import time


from constants import PREFIX, EXP_ID


# job_name --> files to check
# run = zero padded run name
job_files = {'vds'   : ['{PREFIX}/scratch/vds/rrun.cxi'], 
             'events': ['{PREFIX}/scratch/events/rrun_events.h5', '{PREFIX}/scratch/powder/rrun_powder.h5'], 
             'cxi'   : ['{PREFIX}/scratch/cxi/rrun_hits.cxi', '{PREFIX}/scratch/emc/rrun.emc']}


class SLURM_status():
    def __init__(self):
        self.update_every = 10
        self.last_updated = None
        
    def update_slurm(self):
        t = time.time()
        if self.last_updated is None or (t-self.last_updated).seconds < self.update_every :
            jobs = subprocess.run(['squeue', '--format="%.30j"'], stdout=subprocess.PIPE)
            jobs = jobs.stdout.decode('utf-8').split('\n')
            self.slurm_jobs = jobs
            self.last_updated = time.time()
    
    def is_running(self, job_name, run):
        self.update_slurm()
        
        job_string = f'{job_name}-{EXP_ID}-{run}'
        
        # is job running?   
        s = f'*{job_string}*"'
        match = fnmatch.filter(self.slurm_jobs, s)
        running = len(match) > 0
        return running


class LOG_status():
    
    def __init__(self):
        self.log_dir = f'{PREFIX}/scratch/log'
    
    def check_log(self, job_name, run):
        log_files  = glob.glob(f'{self.log_dir}/{job_name}-{EXP_ID}-*-{run}')
        log_files.sort(key = os.path.getmtime, reverse = True)
        
        is_log_file = len(log_file) > 0 
        log_success = False
        if is_log_file :
            log_file = log_file[0]
             
            # has the job finished successfully?
            # check if log file exist and has "{job_name} done"     
            cmd = subprocess.run(['grep', f'{job_name} done', log_file])
            log_success = len(cmd.stdout) > 0
        return is_log_file, log_success

class FILE_status():
    
    def __init__(self):
        pass

    def check_files(self, job_name, run):
        assert(job_name in job_files) 
            
        file_status = True
        for file in job_files[job_name]:
            s = file.replace('run', f'{run:04}')
            if not os.path.exists(s) :
                file_status = False
        return file_status
        


