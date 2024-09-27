from functools import reduce
from collections import OrderedDict
import os
import json
import numpy as np
import datetime
import authenticate_mdc
from metadata_client.metadata_client import MetadataClient
import socket
import glob
import fnmatch
import h5py
import subprocess
import time

from constants import NPULSES_DATASET, NPULSES_DA_NUM, PREFIX, EXP_ID
import get_job_status

hostname = socket.gethostname()
running_on_maxwell = 'desy' in hostname

def deep_get(dictionary, keys, default=None):
    return reduce(lambda d, key: d.get(key, default) if isinstance(d, dict) else default, keys.split("."), dictionary)

def get_date(s):
    """
    2024-02-15T06:01:21.000+01:00 -> 2024-02-15
    """
    return str(datetime.datetime.fromisoformat(s).date())

def get_time(s):
    """
    2024-02-15T06:01:21.000+01:00 -> 06:01:21
    """
    return str(datetime.datetime.fromisoformat(s).time())

def get_duration(s0, s1):
    if s0 is not None and s1 is not None :
        d0 = datetime.datetime.fromisoformat(s0)
        d1 = datetime.datetime.fromisoformat(s1)
        out = str(d1-d0)
    else :
        out = None
    return out

def get_calibrated_data_status(run_json):
    if run_json['cal_pipeline_reply'] == "Calibration jobs succeeded" :
        out = 'ready'
    # sometimes this is zero when cal is finished
    #if run_json['cal_num_requests'] == 0 :
    #    out = 'not requested'
    else :
        if run_json['cal_last_begin_at'] is not None :
            out = 'running'
        
        else :
            out = 'not running'
    return out


def get_vds_file_status(run, slurm_status, log_status, file_status):
    job_name = 'vds'

    if not running_on_maxwell :
        return None
        
    # check if it is running
    is_running = slurm_status.is_running(job_name, run)
    
    if is_running :
        out = 'running'
    else :
        # check logs
        is_log_file, log_success = log_status.check_log(job_name, run)

        
        # check files
        files_ok = file_status.check_files(job_name, run)
        
        if is_log_file :
            if files_ok and log_success :
                out = 'ready'
            # there is log file, it's not running, but no success line
            else :
                out = 'error'
        else :
            # log files were deleted or lost
            if files_ok :
                out = 'ready'
            # no job was submitted
            else :
                out = ''
    
    #print(job_name, run, is_running, is_log_file, files_ok) 
    return out


def get_num_trains(f, l):
    if f is not None and l is not None :
        out = l - f
    else :
        out = None
    return out
    


def get_num_pulses(run):
    npulses = None
    if running_on_maxwell :
        s = PREFIX+'/raw/r%.4d/*DA%.2d*.h5'%(run, NPULSES_DA_NUM)
        flist = sorted(glob.glob(s))
        if len(flist) > 0 :
            with h5py.File(flist[0]) as f:
                if NPULSES_DATASET in f:
                    npulses = int(f[NPULSES_DATASET][0])
    return npulses
    

def get_num_hits(run):
    hits = None
    if running_on_maxwell :
        fnam = PREFIX + '/scratch/events/r%.4d_events.h5' % run
        if os.path.exists(fnam) :
            print('get_num_hits:', fnam, os.path.exists(fnam))
            with h5py.File(fnam) as f:
                if 'is_hit' in f :
                    hits = int(np.sum(f['is_hit'][()]))
                    print('hits = ', hits)
    return hits
    

def get_events_file_status(run, slurm_status, log_status, file_status):
    job_name = 'events'

    if not running_on_maxwell :
        return None
        
    # check if it is running
    is_running = slurm_status.is_running(job_name, run)
    
    if is_running :
        out = 'running'
    else :
        # check logs
        is_log_file, log_success = log_status.check_log(job_name, run)
        
        # check files
        files_ok = file_status.check_files(job_name, run)

        if is_log_file :
            if files_ok and log_success :
                out = 'ready'
            # there is log file, it's not running, but no success line or no file
            else :
                out = 'error'
        else :
            # log files were deleted or lost
            if files_ok :
                out = 'ready'
            # no job was submitted
            else :
                out = ''
    
    #print(job_name, run, is_running, is_log_file, files_ok, 'out', out)
    return out

def get_cxi_file_status(run, slurm_status, log_status, file_status):
    job_name = 'cxi'

    if not running_on_maxwell :
        return None
        
    # check if it is running
    is_running = slurm_status.is_running(job_name, run)
    
    if is_running :
        out = 'running'
    else :
        # check logs
        is_log_file, log_success = log_status.check_log(job_name, run)
        
        # check files
        files_ok = file_status.check_files(job_name, run)
        
        if is_log_file :
            if files_ok and log_success :
                out = 'ready'
            # there is log file, it's not running, but no success line
            else :
                out = 'error'
        else :
            # log files were deleted or lost
            if files_ok :
                out = 'ready'
            # no job was submitted
            else :
                out = ''
    
    #print(job_name, run, is_running, is_log_file, files_ok) 
    return out

def get_sizing_file_status(run, slurm_status, log_status, file_status):
    job_name = 'sizing'
    
    if not running_on_maxwell :
        return None
        
    # check if it is running
    is_running = slurm_status.is_running(job_name, run)
    
    if is_running :
        out = 'running'
    else :
        # check logs
        is_log_file, log_success = log_status.check_log(job_name, run)
        
        # check files
        files_ok = file_status.check_files('sizing', run)
        
        if is_log_file :
            if files_ok and log_success :
                out = 'ready'
            # there is log file, it's not running, but no success line
            else :
                out = 'error'
        else :
            # log files were deleted or lost
            if files_ok :
                out = 'ready'
            # no job was submitted
            else :
                out = ''
    
    #print(job_name, run, is_running, is_log_file, files_ok) 
    return out

def get_static_emc_file_status(run, slurm_status, log_status, file_status):
    job_name = 'static_emc'

    if not running_on_maxwell :
        return None
        
    # check if it is running
    is_running = slurm_status.is_running(job_name, run)
    
    if is_running :
        out = 'running'
    else :
        # check logs
        is_log_file, log_success = log_status.check_log(job_name, run)

        # check files
        files_ok = file_status.check_files(job_name, run)

        if is_log_file :
            if files_ok and log_success :
                out = 'ready'
            # there is log file, it's not running, but no success line
            else :
                out = 'error'
        else :
            # log files were deleted or lost
            if files_ok :
                out = 'ready'
            # no job was submitted
            else :
                out = ''
    
    #print(job_name, run, is_running, is_log_file, files_ok) 
    return out


class Run_table():
    """
    Call the EuXFEL metadata catalog to get run info
    This table can then be written to file or sent to a google sheet
    """
    
    def __init__(self, proposal_number, credentials = 'credentials_mdc.json'):
        # communication object for requests to the mdc api
        self.comm = authenticate_mdc.get(credentials)
        self.proposal_number = proposal_number

        msg = MetadataClient.get_proposal_runs(self.comm, self.proposal_number)
        self.proposal_id = msg['data']['proposal']['id']
        
        # get sample name by id for later
        samples = self.comm.get_all_samples_by_proposal_id_api(self.proposal_id)
        self.sample_names = {s['id']: s['name'] for s in samples.json()}

        # get run type by id for later
        experiments = self.comm.get_all_experiments_by_proposal_id_api(self.proposal_id)
        self.run_types = {s['id']: s['name'] for s in experiments.json()}
        
        # get running jobs on maxwell
        self.slurm_status = get_job_status.SLURM_status()
        
        # check log files
        self.log_status = get_job_status.LOG_status()
        
        # check files
        self.file_status = get_job_status.FILE_status()
        
        # heading run_stats maping
        headings = OrderedDict([('Run number', lambda x: x['run_number']),
                                ('Date',       lambda x: get_date(x['begin_at'])),
                                ('Start Time', lambda x: get_time(x['begin_at'])),
                                ('Duration',   lambda x: get_duration(x['begin_at'], x['end_at'])),
                                ('Run Type',   lambda x: deep_get(x, 'experiment.name')),
                                ('Run Type',   lambda x: self.run_types[x['experiment_id']]),
                                ('Sample',     lambda x: self.sample_names[x['sample_id']]),
                                ('Num Trains', lambda x: get_num_trains(x['first_train'], x['last_train'])),
                                ('Num Pulses', lambda x: get_num_pulses(x['run_number'])),
                                ('Num Hits',   lambda x: get_num_hits(x['run_number'])),
                                ('Hit Rate',   lambda x: None),
                                ('Calib',      lambda x: get_calibrated_data_status(x)),
                                ('VDS',        lambda x: get_vds_file_status(x['run_number'], self.slurm_status, self.log_status, self.file_status)),
                                ('Events',     lambda x: get_events_file_status(x['run_number'], self.slurm_status, self.log_status, self.file_status)),
                                ('CXI',        lambda x: get_cxi_file_status(   x['run_number'], self.slurm_status, self.log_status, self.file_status)),
                                ('EMC files',  lambda x: get_cxi_file_status(   x['run_number'], self.slurm_status, self.log_status, self.file_status)),
                                ('Sizing',     lambda x: get_sizing_file_status(  x['run_number'], self.slurm_status, self.log_status, self.file_status)),
                                ('Static EMC', lambda x: get_static_emc_file_status(   x['run_number'], self.slurm_status, self.log_status, self.file_status)),
                                ('Comments',   lambda x: None)])
        self.headings = headings

    def update(self): 
        # dictionary
        msg = MetadataClient.get_proposal_runs(self.comm, self.proposal_number, page_size=1000)
        
        assert(msg['success'] == True)
        
        run_stats = msg['data']['runs']
        
        # sort by run number
        r = [run['run_number'] for run in run_stats]
        runs_sorted = np.argsort(r)
        
        run_table = []
        for r in runs_sorted:
            run = run_stats[r]
            row = []
            for h, l in self.headings.items():
                row.append(l(run))
                #try :
                #    row.append(l(run))
                #except Exception as e:
                #    #print(run)
                #    print(e)
            run_table.append(row)
        
        # make run dict for log
        # run_dict[run_number] = {'VDS': 'ready', ...}
        run_dict = OrderedDict()
        for row in run_table:
            run_number = int(row[0])
            run_dict[run_number] = OrderedDict()
            for r, h in zip(row, self.headings.keys()):
                run_dict[run_number][h] = r
        
        run_dict['last_update'] = time.time()
            
        # save run table, json or pickle?
        if running_on_maxwell :
            run_log = f'{PREFIX}/scratch/log/run_table.json'
        else :
            run_log = f'run_table.json'
        
        with open(run_log, 'w') as f:
            json.dump(run_dict, f, indent=4)
        
        return [h for h in self.headings], run_table


