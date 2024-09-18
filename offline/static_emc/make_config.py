import argparse
import pathlib

parser = argparse.ArgumentParser(description='Save hits in photon units to a cxi file')
parser.add_argument('runs', type=int, help='Run/s number', nargs='+')
parser.add_argument('-c', '--config', type=str, help='config template to use', default = 'config_template.py')
parser.add_argument('--cxi_dir', type=str, help='location of directory containing cxi files', required = True)
parser.add_argument('--output_dir', type=str, help='location of directory to put outut folder', required = True)
parser.add_argument('--mask', type=str, help='mask file', required = True)
parser.add_argument('--geom', type=str, help='crystfel geometry file', required = True)

args = parser.parse_args()

args.config = pathlib.Path(__file__).parent.resolve().joinpath(args.config).resolve()

from constants import PREFIX
import numpy as np
import sys

print('template file:', args.config, file=sys.stderr)

cxi_files = [str(pathlib.Path(args.cxi_dir).joinpath(f'r{run:04}_hits.cxi').resolve()) for run in args.runs]

if len(args.runs) == 1 :
    run = args.runs[0]
    temp = f'r{run:04}'
else :
    s = '-'.join(str(r) for r in args.runs)
    temp = f'runs_{s}'

working_directory = pathlib.Path(args.output_dir).joinpath(temp).resolve()

# make working directory if it does not exist
pathlib.Path(working_directory).mkdir(exist_ok=True)

config_fnam = pathlib.Path(working_directory).joinpath('config.py').resolve()

# read template
config_template = open(args.config, 'r').readlines()

# write data source
i = np.where(['data = []' in line for line in config_template])[0][0]
config_template[i] = f'data = {cxi_files}\n'

# write working directory
i = np.where(["working_dir = ''" in line for line in config_template])[0][0]
config_template[i] = f"working_dir = '{working_directory}'\n"

# write geom 
i = np.where(["geom_fnam = ''" in line for line in config_template])[0][0]
config_template[i] = f"geom_fnam = '{args.geom}'\n"

# write mask 
i = np.where(["static_emc_mask_fnam = ''" in line for line in config_template])[0][0]
config_template[i] = f"static_emc_mask_fnam = '{args.mask}'\n"

# write config
with open(config_fnam, 'w') as f:
    f.writelines(config_template)


# write config file name to stdout for calling bash script 
print(config_fnam)
