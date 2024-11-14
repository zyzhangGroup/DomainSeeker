import os,sys
import numpy as np
import wget

UniprotID_list_path=sys.argv[1]
output_dir=sys.argv[2]

output_pdb_dir=os.path.join(output_dir,'pdb_files')
os.makedirs(output_pdb_dir,exist_ok=True)
missing_pdb_log_path=os.path.join(output_dir,'missing_pdb.log')
missing_pdb_log=open(missing_pdb_log_path,'w')

output_pae_dir=os.path.join(output_dir,'pae_files')
os.makedirs(output_pae_dir,exist_ok=True)
missing_pae_log_path=os.path.join(output_dir,'missing_pae.log')
missing_pae_log=open(missing_pae_log_path,'w')

UniprotID_list=np.loadtxt(UniprotID_list_path,dtype=str)
n=len(UniprotID_list)
for i in range(n):
    ID=UniprotID_list[i]
    pdb_path=os.path.join(output_pdb_dir,ID+'.pdb')
    if not os.path.exists(pdb_path):
        pdb_url=f"https://alphafold.ebi.ac.uk/files/AF-{ID}-F1-model_v4.pdb"
        try:
            wget.download(pdb_url,pdb_path,bar='')
        except Exception as e:
            missing_pdb_log.write(ID+'\n')
    
    pae_path=os.path.join(output_pae_dir,ID+'.json')
    if not os.path.exists(pae_path):
        pae_url=f"https://alphafold.ebi.ac.uk/files/AF-{ID}-F1-predicted_aligned_error_v4.json"
        try:
            wget.download(pae_url,pae_path,bar='')
        except Exception as e:
            missing_pae_log.write(ID+'\n')
    
    print(f"{i}/{n}",end="\r")


missing_pdb_log.close()
missing_pae_log.close()