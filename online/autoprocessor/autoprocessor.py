import config
import subprocess
import json
import time

from constants import PREFIX

def run_pipeline(sleep = 10):
    while True :
        # load run_table
        fnam = f'{PREFIX}/scratch/log/run_table.json'
        
        with open(fnam, 'r') as f:
            run_table = json.load(f)
        
        for run in run_table:   
            for script in config.scripts :
                run_script, command = script(run, run_table)
                 
                if run_script :
                    print(f'running {command}')
                    subprocess.run(command, shell=True, check=True, text=True)
        
        time.sleep(sleep)
                
        

    
    
    
