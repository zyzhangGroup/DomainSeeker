import MDAnalysis as mda
import numpy as np
import sys,os

argv=sys.argv
domain_dir=argv[1]
fitout_subdir=argv[2]
domain_list=argv[3:]    #without suffix


def get_transformation_matrix_list(log_path):
    log_data=np.loadtxt(log_path,dtype=float)
    transform_matrix_list=log_data[:,2:].reshape((-1,3,4))
    return transform_matrix_list

def transform_positions(original_xyzs,transformation_matrix):
    rotation_matrix=transformation_matrix[:,:3]
    translation_vector=transformation_matrix[:,-1]
    return original_xyzs.dot(rotation_matrix.T)+translation_vector

def transform_save_single_fit(u,original_xyzs,transformation_matrix,save_path):
    transformed_xyzs=transform_positions(original_xyzs,transformation_matrix)
    u.atoms.positions=transformed_xyzs
    u.atoms.write(save_path)

os.makedirs(os.path.join(fitout_subdir,"fitted_domains"),exist_ok=True)

for domain in domain_list:
    domain_path=os.path.join(domain_dir,domain+'.pdb')
    u=mda.Universe(domain_path)
    original_xyzs=u.atoms.positions

    log_path=os.path.join(fitout_subdir,"fitlogs",domain+'.log')
    transformation_matrix_list=get_transformation_matrix_list(log_path)

    n=len(transformation_matrix_list)
    for i in range(n):
        transformation_matrix=transformation_matrix_list[i]
        save_path=os.path.join(fitout_subdir,"fitted_domains",f"{domain}_fit_{i:03d}.pdb")
        transform_save_single_fit(u,original_xyzs,transformation_matrix,save_path)