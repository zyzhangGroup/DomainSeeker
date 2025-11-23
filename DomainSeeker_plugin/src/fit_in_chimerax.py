import os,sys
import numpy as np
from chimerax.core.commands import run

domain_path=sys.argv[1]
map_path=sys.argv[2]
output_subdir=sys.argv[3]
ref_map_threshold=float(sys.argv[4])
resolution=float(sys.argv[5])
n_search=int(sys.argv[6])

# open files and fit
run(session,f'open {map_path}')
run(session,f"volume threshold #1 minimum {ref_map_threshold} set 0")
run(session,f'volume #2 level {ref_map_threshold}')
run(session,f'open {domain_path}')
fits=run(session,f"fitmap #3 inmap #2 resolution {resolution} search {n_search}")
log_data=[]
for fit in fits:
    transform_string_list=[]
    for element in fit.model_transforms()[0][1].matrix.reshape(-1):
        transform_string_list.append(f"{element:10.5f}")
    if fit.correlation() >= 0.00001:
        log_data.append([f"{fit.correlation():8.5f}",f"{fit.hits():8d}"]+transform_string_list)

# write transformation info into log
if len(log_data) >=1:
    os.makedirs(output_subdir,exist_ok=True)
    fitlog_subdir=os.path.join(output_subdir,"fitlogs")
    os.makedirs(fitlog_subdir,exist_ok=True)
    log_path=os.path.join(fitlog_subdir,os.path.basename(domain_path).replace('pdb','log'))
    print(f"{output_subdir=}")
    print(f"{log_path=}")
    np.savetxt(log_path,log_data,fmt="%s")

