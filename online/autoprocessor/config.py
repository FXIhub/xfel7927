import pathlib

# get offline dir 
# .../xfel7927/
root = pathlib.Path(__file__).parent.parent.parent.resolve()
script_dir = root.joinpath('offline/').resolve()
print('script directory:', script_dir)

s = ['ready', 'error', 'running']

samples = ['DNA Origami Cube', 'DNA-origami Pointer', 'Erythrocruorin']

def vds(run, run_table):
    run_script = False
    if all(i not in run_table[run]['VDS'] for i in s):
    
        if run_table[run]['Calib'] == 'ready':
            if run_table[run]['Sample'] in samples :
                run_script = True
            if run_table[run]['Run Type'] in 'Gas background' :
                run_script = True
    
    command = f'. {script_dir}/submit_vds.sh {run}'
    return run_script, command
            
def events(run, run_table):
    run_script = False
    if all(i not in run_table[run]['Events'] for i in s):
        if run_table[run]['VDS'] == 'ready':
            if run_table[run]['Sample'] in samples :
                run_script = True
            if run_table[run]['Run Type'] in 'Gas background' :
                run_script = True
    
    command = f'. {script_dir}/submit_events.sh {run}'
    return run_script, command

def cxi(run, run_table):
    run_script = False
    if all(i not in run_table[run]['CXI'] for i in s):
        if run_table[run]['Events'] == 'ready':
            if run_table[run]['Sample'] in samples :
                run_script = True
        
    command = f". {script_dir}/submit_cxi.sh {run} {run_table[run]['Sample']}"
    return run_script, command

def static_emc(run, run_table):
    run_script = False
    if all(i not in run_table[run]['Static EMC'] for i in s):
        if run_table[run]['CXI'] == 'ready':
            if run_table[run]['Sample'] in samples :
                run_script = True
        
    command = f". {script_dir}/submit_static_emc.sh {run}"
    return run_script, command

def sizing(run, run_table):
    run_script = False
    if all(i not in run_table[run]['Sizing'] for i in s):
        if run_table[run]['CXI'] == 'ready':
            if run_table[run]['Sample'] in samples :
                run_script = True
        
    command = f". {script_dir}/submit_sizing.sh {run}"
    return run_script, command

scripts = [vds, events, cxi, static_emc, sizing]
