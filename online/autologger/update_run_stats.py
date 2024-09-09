from functools import reduce
from collections import OrderedDict
import os
import json
import numpy as np
import datetime
import authenticate_mdc
from metadata_client.metadata_client import MetadataClient


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
            out = 'available'
        
        elif run_json['cal_last_begin_at'] is not None :
            out = 'running'
        
        else :
            out = 'not running'
    return out

# heading run_stats maping
headings = OrderedDict([('Run number', lambda x: x['run_number']),
                        ('Date',       lambda x: get_date(x['begin_at'])),
                        ('Start Time', lambda x: get_time(x['begin_at'])),
                        ('Duration',   lambda x: get_duration(x['begin_at'], x['end_at'])),
                        ('Run Type',   lambda x: deep_get(x, 'experiment.name')),
                        ('Sample',     lambda x: deep_get(x, 'sample.name')),
                        ('Num Trains', lambda x: x['last_train'] - x['first_train']),
                        ('Num Pulses', lambda x: None),
                        ('Num Hits',   lambda x: None),
                        ('Hit Rate',   lambda x: None),
                        ('Calib',      lambda x: get_calibrated_data_status(x)),
                        ('VDS',        lambda x: None),
                        ('Events',     lambda x: None),
                        ('Comments',   lambda x: None)])


class Run_table():
    """
    Call the EuXFEL metadata catalog to get run info
    This table can then be written to file or sent to a google sheet
    """

    def __init__(self, proposal_number, credentials = 'credentials_mdc.json'):
        # communication object for requests to the mdc api
        self.comm = authenticate_mdc.get(credentials)
        self.proposal_number = proposal_number
        
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
            for h, l in headings.items():
                try :
                    row.append(l(run))
                except Exception as e:
                    print(run)
                    print(e)
            
            run_table.append(row)
        
        return [h for h in headings], run_table


