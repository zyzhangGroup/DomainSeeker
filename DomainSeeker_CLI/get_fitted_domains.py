import os,argparse
import MDAnalysis as mda
import numpy as np

parser = argparse.ArgumentParser(description="Create .pdb files of fitted domains")
parser.add_argument("domain_dir")
parser.add_argument("fitout_subdir", help="fitout_subdir is the subdirectory (e.g., 'A.mrc') of the fitout_dir")
parser.add_argument("domain_list", nargs="+")
argument = parser.parse_args()


domain_dir=argument.domain_dir
fitout_subdir=argument.fitout_subdir
domain_list=argument.domain_list


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
    domain_base = "_".join(domain.split("_")[0:2])
    fold = int(domain.split("_")[2])
    domain_path=os.path.join(domain_dir,domain_base+'.pdb')
    u=mda.Universe(domain_path)
    original_xyzs=u.atoms.positions

    log_path=os.path.join(fitout_subdir,"fitlogs",domain_base+'.log')
    transformation_matrix_list=get_transformation_matrix_list(log_path)

    transformation_matrix=transformation_matrix_list[fold]
    save_path=os.path.join(fitout_subdir,"fitted_domains",f"{domain_base}_fit_{fold:03d}.pdb")
    transform_save_single_fit(u,original_xyzs,transformation_matrix,save_path)