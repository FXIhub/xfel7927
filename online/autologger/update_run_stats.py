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

from constants import NPULSES_DATASET, NPULSES_DA_NUM, PREFIX, EXP_ID

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
    d0 = datetime.datetime.fromisoformat(s0)
    d1 = datetime.datetime.fromisoformat(s1)
    return str(d1-d0)

def get_calibrated_data_status(run_json):
    if run_json['cal_num_requests'] == 0 :
        out = 'not requested'
    
    else :
        if run_json['cal_last_end_at'] is not None :
            out = 'ready'
        
        elif run_json['cal_last_begin_at'] is not None :
            out = 'running'
        
        else :
            out = 'not running'
    return out

def get_vds_file_status(run):
    """
    check for the existance of the file, e.g:
        ${PREFIX}/scratch/vds/r0034.cxi
    """
    fnam = f'{PREFIX}/scratch/vds/r{run:04}.cxi'
    if running_on_maxwell :
        if os.path.exists(fnam):
            out = 'ready'
        else :
            out = 'not ready'
    else :
        out = None
    return out


def get_events_file_status(run, jobs):
    """
    check for the existance of the file, e.g:
        ${PREFIX}/scratch/events/events_r0034.h5
    """
    fnam = f'{PREFIX}/scratch/events/events_r{run:04}.h5'
    out = None
    if running_on_maxwell :
        if os.path.exists(fnam):
            out = 'ready'
        else :
            out = 'not ready'
        
        if jobs :
            s = f'*events*{EXP_ID}*-{run}"'
            match = fnmatch.filter(jobs, s)
            print(run, s, match)
            if len(match) > 0:
                out = 'running'
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
        print(fnam, os.path.exists(fnam))
        if os.path.exists(fnam) :
            with h5py.File(fnam) as f:
                if 'is_hit' in f :
                    hits = np.sum(f['is_hit'][()])
                    print(hits)
    return hits
    



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
        self.
        
        # heading run_stats maping
        headings = OrderedDict([('Run number', lambda x: x['run_number']),
                                ('Date',       lambda x: get_date(x['begin_at'])),
                                ('Start Time', lambda x: get_time(x['begin_at'])),
                                ('Duration',   lambda x: get_duration(x['begin_at'], x['end_at'])),
                                ('Run Type',   lambda x: deep_get(x, 'experiment.name')),
                                ('Run Type',   lambda x: self.run_types[x['experiment_id']]),
                                ('Sample',     lambda x: self.sample_names[x['sample_id']]),
                                ('Num Trains', lambda x: x['last_train'] - x['first_train']),
                                ('Num Pulses', lambda x: get_num_pulses(x['run_number'])),
                                ('Num Hits',   lambda x: get_num_hits(x['run_number'])),
                                ('Hit Rate',   lambda x: None),
                                ('Calib',      lambda x: get_calibrated_data_status(x)),
                                ('VDS',        lambda x: get_vds_file_status(x['run_number'])),
                                ('Events',     lambda x: get_events_file_status(x['run_number'], jobs)),
                                ('Comments',   lambda x: None)])
        self.headings = headings

    def get_slurm_jobs(self):
        jobs = None
        if running_on_maxwell :
            jobs = subprocess.run(['squeue', '--format="%.30j"'], stdout=subprocess.PIPE)
            jobs = jobs.stdout.decode('utf-8').split('\n')
        return jobs
        
    def update(self): 
        # dictionary
        msg = MetadataClient.get_proposal_runs(self.comm, self.proposal_number)
        
        print(msg.keys())
        
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
        
        return [h for h in self.headings], run_table


