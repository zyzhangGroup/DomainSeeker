import os,sys
import numpy as np
from chimerax.core.commands import run

domain_path=sys.argv[1]
map_path=sys.argv[2]
output_dir=sys.argv[3]
ref_map_threshold=float(sys.argv[4])
resolution=float(sys.argv[5])
n_search=int(sys.argv[6])

#open files and fit
run(session,f'open {map_path}')
run(session,f"volume threshold #1 minimum {ref_map_threshold} set 0")
run(session,f'volume #2 level {ref_map_threshold}')
run(session,f'open {domain_path}')
fits=run(session,f"fitmap #3 inmap #2 resolution {resolution} search {n_search}")
fit_message=fits[0].fit_message()

#write transform info and scores into log
os.makedirs(output_dir,exist_ok=True)
log_path=os.path.join(output_dir,os.path.basename(domain_path).replace('pdb','log'))
log=open(log_path,'w')
log.write(fit_message)
log.close()

