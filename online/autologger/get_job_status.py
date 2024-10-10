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
import h5py


from constants import PREFIX, EXP_ID


# job_name --> files to check
# run = zero padded run name
job_files = {'vds'    : [f'{PREFIX}/scratch/vds/rrun.cxi'], 
             'events' : [f'{PREFIX}/scratch/events/rrun_events.h5'], 
             'cxi'    : [f'{PREFIX}/scratch/saved_hits/rrun_hits.cxi', f'{PREFIX}/scratch/emc/rrun.emc'],
             'sizing' : [f'{PREFIX}/scratch/saved_hits/rrun_hits.cxi', f'{PREFIX}/scratch/events/rrun_events.h5'],
             'intensity' : [f'{PREFIX}/scratch/log/peak_intensity_report.pdf'],
             'static_emc' : [f'{PREFIX}/scratch/static_emc/rrun/recon.pdf']}


class SLURM_status():
    def __init__(self):
        self.update_every = 10
        self.last_updated = None
        
    def update_slurm(self):
        t = time.time()
        if self.last_updated is None or (t-self.last_updated) > self.update_every :
            print('updating slurm list')
            jobs = subprocess.run('squeue --Format="Name:30,ArrayTaskID:30"', shell=True, check=True, text=True, stdout=subprocess.PIPE)
            jobs = jobs.stdout.split('\n')
            self.slurm_jobs = jobs
            self.last_updated = time.time()
    
    def is_running(self, job_name, run):
        self.update_slurm()
        
        job_string = f'*{job_name}-{EXP_ID}* * {run} *'
        
        # is job running?   
        match = fnmatch.filter(self.slurm_jobs, job_string)
        running = len(match) > 0
        
        if job_name == 'events':
            print('slurm:', job_name, run, running, job_string)

        return running


class LOG_status():
    
    def __init__(self):
        self.log_dir = f'{PREFIX}/scratch/log'
    
    def check_log(self, job_name, run):
        search = f'{self.log_dir}/{job_name}-{EXP_ID}-*-{run}.out'
        log_files  = glob.glob(search)
        log_files.sort(key = os.path.getmtime, reverse = True)
        
        is_log_file = len(log_files) > 0 
        log_success = False
        #print('log:', run, search, is_log_file, log_files, log_success)
        if is_log_file :
            log_file = log_files[0]
             
            # has the job finished successfully?
            # check if log file exist and has "{job_name} done"     
            grp = f'grep "{job_name} done" {log_file}'
            cmd = subprocess.run(grp, shell=True, text=True, stdout=subprocess.PIPE).stdout
            #log_success = len(cmd.stdout) > 0

            if cmd :
                log_success = True
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
            
            elif job_name == 'sizing' and 'cxi' in file:
                with h5py.File(s, 'r') as f:
                    key = 'entry_1/sizing/short_axis_diameter'
                    if key not in f :
                        file_status = False

            elif job_name == 'sizing' and 'events' in file:
                with h5py.File(s, 'r') as f:
                    key = 'sizing/short_axis_diameter'
                    if key not in f :
                        file_status = False
            
            #print(s, file_status)
        return file_status
        


