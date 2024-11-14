import MDAnalysis as mda
import numpy as np
import sys,os

argv=sys.argv
domain_dir=argv[1]
fitout_dir=argv[2]
domain_list=argv[3:]    #无后缀

for domain in domain_list:
    domain_path=os.path.join(domain_dir,domain+'.pdb')
    log_path=os.path.join(fitout_dir,domain+'.log')
    u=mda.Universe(domain_path)
    xyzs=u.atoms.positions
    transformation_matrix=np.loadtxt(log_path,skiprows=5,max_rows=3,dtype=float)
    rotation_matrix=transformation_matrix[:,:3]
    translation_vector=transformation_matrix[:,-1]
    xyzs=xyzs.dot(rotation_matrix.T)+translation_vector
    u.atoms.positions=xyzs
    fitted_domain_path=os.path.join(fitout_dir,domain+'_fitted.pdb')
    u.atoms.write(fitted_domain_path)