import pathlib

# get offline dir 
# .../xfel7927/
root = pathlib.Path(__file__).parent.parent.parent.resolve()
script_dir = root.joinpath('offline/').resolve()
print('script directory:', script_dir)

s = ['ready', 'error', 'running']

samples = ['DNA Origami Cube', 'DNA Origami Twister', 'DNA-origami Pointer', 'Erythrocruorin', 'AuNP', 'RNA polymerase', 'De novo designer octahedral nanoparticle']

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
        
    command = f". {script_dir}/submit_cxi.sh {run} {run_table[run]['Sample'].replace(' ', '_')}"
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

def intensity(run, run_table):
    run_script = False
    if all(i not in run_table[run]['Peak Intensity'] for i in s):
        if run_table[run]['Sizing'] == 'ready':
            run_script = True
        
    command = f". {script_dir}/submit_intensity.sh {run}"
    return run_script, command

def mask(run, run_table):
    run_script = False
    if all(i not in run_table[run]['Mask'] for i in s):
        if run_table[run]['Powder'] == 'ready':
            run_script = True
        
    command = f". {script_dir}/submit_mask.sh {run}"
    return run_script, command

def powder(run, run_table):
    run_script = False
    if all(i not in run_table[run]['Powder'] for i in s):
        if run_table[run]['Calib'] == 'ready':
            run_script = True
        
    command = f". {script_dir}/submit_powder.sh {run}"
    return run_script, command

#scripts = [vds, events, cxi, sizing, mask, powder]
scripts = [vds, events, sizing, mask, powder]
