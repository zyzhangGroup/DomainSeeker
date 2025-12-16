# 全局异常处理
import domainseeker_errorLog
#----------------------------------------------------------------------------------------------------
import os,sys
import wget
import numpy as np
from tqdm import tqdm

UniprotID_list_path=sys.argv[2]
output_dir=sys.argv[3]

if len(sys.argv)>4:
    output_pdb_dir=sys.argv[4]
    output_pae_dir=sys.argv[5]
else:
    output_pdb_dir=os.path.join(output_dir,'pdb_files')
    output_pae_dir=os.path.join(output_dir,'pae_files')


os.makedirs(output_pdb_dir,exist_ok=True)
missing_pdbs=[]


os.makedirs(output_pae_dir,exist_ok=True)
missing_paes=[]


UniprotID_list=np.loadtxt(UniprotID_list_path,dtype=str,ndmin=1)
n=len(UniprotID_list)
for i in tqdm(range(n),desc='Downloading',file=sys.stdout):
    ID=UniprotID_list[i]
    pdb_path=os.path.join(output_pdb_dir,ID+'.pdb')
    if not os.path.exists(pdb_path):
        pdb_url=f"https://alphafold.ebi.ac.uk/files/AF-{ID}-F1-model_v6.pdb"
        try:
            wget.download(pdb_url,pdb_path,bar='')
        except Exception as e:
            missing_pdbs.append(ID)

    pae_path=os.path.join(output_pae_dir,ID+'.json')
    if not os.path.exists(pae_path):
        pae_url=f"https://alphafold.ebi.ac.uk/files/AF-{ID}-F1-predicted_aligned_error_v6.json"
        try:
            wget.download(pae_url,pae_path,bar='')
        except Exception as e:
            missing_paes.append(ID)


missing_pdb_log_path=os.path.join(output_dir,'missing_pdb.log')
missing_pdb_log=open(missing_pdb_log_path,'w')
if len(missing_pdbs)>0:
    for ID in missing_pdbs:
        missing_pdb_log.write(ID+'\n')
missing_pdb_log.close()

missing_pae_log_path=os.path.join(output_dir,'missing_pae.log')
missing_pae_log=open(missing_pae_log_path,'w')
if len(missing_paes)>0:
    for ID in missing_paes:
        missing_pae_log.write(ID+'\n')
missing_pae_log.close()



missing_pdb_log.close()
missing_pae_log.close()

print("Done fetching PDB and PAE files.")