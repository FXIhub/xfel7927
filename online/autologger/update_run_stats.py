from functools import reduce
from collections import OrderedDict
import os
import json
import numpy as np
import datetime
import authenticate_mdc


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
                        ('Hit Rate',   lambda x: None)])

# remember experiment ids
experiment_ids = []

def get_experiment_id_from_run_id(run_id):
    creds = authenticate_mdc.get()
    run_stats = os.popen(
    f"""
    curl -X 'GET' \
      'https://in.xfel.eu/metadata/api/runs/{run_id}' \
      -H 'accept: application/json; version=1' \
      -H 'X-USER-EMAIL: {creds['email']}' \
      -H 'Authorization: Bearer {creds['token']}'
    """).read()
    
    # list of dictionaries
    run_stats = json.loads(run_stats)
    
    experiment_id = run_stats['experiment_id']
    return experiment_id

def get_experiment_ids_from_proposal_number(prop_no, creds):
    print(f'getting proposal id from proposal number {prop_no}')
    prop_stats = os.popen(
    f"""
    curl -X 'GET' \
      'https://in.xfel.eu/metadata/api/proposals' \
      --data number={prop_no} \
      -H 'accept: application/json; version=1' \
      -H 'X-USER-EMAIL: {creds['email']}' \
      -H 'Authorization: Bearer {creds['token']}'
    """).read()

    # list of dictionaries
    prop_stats = json.loads(prop_stats)

    prop_id = prop_stats[0]['id']
    print(f'proposal id: {prop_id}')
    
    print('getting experiment ids from proposal number')
    exp_stats = os.popen(
    f"""
    curl -X 'GET' \
      'https://in.xfel.eu/metadata/api/experiments' \
      --data proposal_id={prop_id} \
      -H 'accept: application/json; version=1' \
      -H 'X-USER-EMAIL: {creds['email']}' \
      -H 'Authorization: Bearer {creds['token']}'
    """).read()

    # list of dictionaries
    exp_stats = json.loads(exp_stats)
    
    global experiment_ids 
    experiment_ids = [exp['id'] for exp in exp_stats]
    

# I cannot get this to return all runs! actually this is becuase there are two different experiment ids in this proposal
# it seems that every run type has its own experiment id... how to get run info based on run or proposal number?!
def update_run_table(proposal_no):
    creds = authenticate_mdc.get()

    global experiment_ids
    if len(experiment_ids) == 0 :
        get_experiment_ids_from_proposal_number(proposal_no, creds)
        
    print('calling EuXFEL metadata api')
    run_stats = [] 
    for experiment_id in experiment_ids :
        call = f"""
        curl -X 'GET' \
          'https://in.xfel.eu/metadata/api/runs?experiment_id={experiment_id}' \
          --data page_size=500 \
          -H 'accept: application/json; version=1' \
          -H 'X-USER-EMAIL: {creds['email']}' \
          -H 'Authorization: Bearer {creds['token']}'
        """
        print(call)
        run_json = os.popen(call).read()
        
        # list of dictionaries
        run_stats += json.loads(run_json)
        print('found', len(run_stats),'runs in experiment', experiment_id)

    print(json.dumps(run_stats, indent=2))
    
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
