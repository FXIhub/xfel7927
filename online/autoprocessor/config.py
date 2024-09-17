import pathlib

# get offline dir 
# .../xfel7927/
root = pathlib.Path(__file__).parent.parent.parent.resolve()
script_dir = root.joinpath('offline/').resolve()
print('script directory:', script_dir)

s = ['ready', 'error', 'running']

def vds(run, run_table):
    run_script = False
    if all(i not in run_table[run]['VDS'] for i in s):
    
        if run_table[run]['Calib'] == 'ready':
            if run_table[run]['Sample'] == 'Cubic DNA origami':
                run_script = True
    
    command = f'. {script_dir}/submit_vds.sh {run}'
    return run_script, command
            
def events(run, run_table):
    run_script = False
    if all(i not in run_table[run]['Events'] for i in s):
        if run_table[run]['VDS'] == 'ready':
            if run_table[run]['Sample'] == 'Cubic DNA origami':
                run_script = True
    
    command = f'. {script_dir}/submit_events.sh {run}'
    return run_script, command

def cxi(run, run_table):
    run_script = False
    if all(i not in run_table[run]['CXI'] for i in s):
        if run_table[run]['Events'] == 'ready':
            if run_table[run]['Sample'] == 'Cubic DNA origami':
                run_script = True
        
    command = f". {script_dir}/submit_cxi.sh {run} {run_table[run]['Sample']}"
    return run_script, command

scripts = [vds, events, cxi]
