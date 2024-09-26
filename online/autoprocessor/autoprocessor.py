import config
import subprocess
import json
import time

from constants import PREFIX

# how to prevent re-running scripts at next loop?

def run_pipeline(sleep = 10):
    t = time.time()
    
    while True :
        # load run_table
        fnam = f'{PREFIX}/scratch/log/run_table.json'
        
        with open(fnam, 'r') as f:
            run_table = json.load(f)
        
        print('last update:', (run_table['last_update'] - t),' seconds ago')
        
        if (run_table['last_update'] - t) > 0 :
            for run in run_table:   
                if run == 'last_update' : 
                    continue 
                
                for script in config.scripts :
                    run_script, command = script(run, run_table)
                     
                    if run_script :
                        print(f'running {command}')
                        subprocess.run(command, shell=True, check=True, text=True)
        
        t = time.time()
        time.sleep(sleep)
                
        
if __name__ == '__main__':
    run_pipeline()
    
    
    
