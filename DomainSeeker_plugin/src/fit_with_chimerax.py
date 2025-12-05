# 全局异常处理
import domainseeker_errorLog
#----------------------------------------------------------------------------------------------------

import os,sys
import subprocess,multiprocessing
import MDAnalysis as mda
import numpy as np
import scipy
import mrcfile
from tqdm import tqdm
import warnings


# 不同系统的ChimeraX
if sys.platform == 'win32':
    ChimeraX="ChimeraX-console"     # windows
else:
    ChimeraX="chimerax"     # Linux

def fit(params):
    domain_filename, ref_map_filename = params
    domain_path=os.path.join(domain_dir,domain_filename).replace("\\","/")
    map_path=os.path.join(map_dir,ref_map_filename).replace("\\","/")
    output_subdir=os.path.join(fitout_dir,ref_map_filename).replace("\\","/")
    cmd=f'{ChimeraX} --nogui --script "{script_dir}/fit_in_chimerax.py {domain_path} {map_path} {output_subdir} {ref_map_threshold} {resolution} {n_search}" --exit'
    # 舍弃输出
    subprocess.run(cmd,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    
# 全局变量
domain_dir=os.path.realpath(sys.argv[2]).replace("\\","/")
map_dir=os.path.realpath(sys.argv[3]).replace("\\","/")
fitout_dir=os.path.realpath(sys.argv[4]).replace("\\","/")
ref_map_threshold=float(sys.argv[5])
resolution=float(sys.argv[6])
n_search=int(sys.argv[7])

negtive_laplacian_cutoff=float(sys.argv[8])
positive_laplacian_cutoff=float(sys.argv[9])
fit_map_laplacian_cutoff_low=-2
fit_map_laplacian_cutoff_high=15

n_process=int(sys.argv[10])

script_dir=os.path.dirname(os.path.realpath(__file__)).replace("\\","/")

## 局部打分
def get_ref_map_mat_laplacian_and_params(ref_map_path):
    ref_map=mrcfile.open(ref_map_path)
    ref_map_mat=ref_map.data.T
    ref_map_mat=ref_map_mat*(ref_map_mat>=ref_map_threshold)
    origin=np.array(ref_map.header.origin.tolist())
    lengths=np.array(ref_map.header.cella.tolist())
    grid_shape=np.array([ref_map.header.nx,ref_map.header.ny,ref_map.header.nz])
    voxel=lengths[0]/grid_shape[0]
    density_param_dict={'origin':origin,'lengths':lengths,'grid_shape':grid_shape,'voxel':voxel}
    ref_map_mat_laplacian=scipy.ndimage.convolve(ref_map_mat,laplacian_kernel,mode='constant')
    return [ref_map_mat_laplacian,density_param_dict]

# generate density matrix

weight_dict={'C':12.011,'N':14.007,'O':16.000,'S':32.060}

# Laplacian kernel
laplacian_kernel=np.zeros([3,3,3])
laplacian_kernel[1,1,1]=-6
laplacian_kernel[1,1,0]=1
laplacian_kernel[1,1,2]=1
laplacian_kernel[1,0,1]=1
laplacian_kernel[1,2,1]=1
laplacian_kernel[0,1,1]=1
laplacian_kernel[2,1,1]=1


def add_atom(data,voxel,origin,xyz,weight):
    gxyz = (xyz-origin)/voxel
    xyz0 = np.floor(gxyz).astype(int)
    xyz1 = xyz0 + 1
    dxyz1 = xyz1 - gxyz
    dxyz0 = 1 - dxyz1
    data[xyz0[0], xyz0[1], xyz0[2]] += dxyz1[0]*dxyz1[1]*dxyz1[2]*weight
    data[xyz1[0], xyz0[1], xyz0[2]] += dxyz0[0]*dxyz1[1]*dxyz1[2]*weight
    data[xyz0[0], xyz1[1], xyz0[2]] += dxyz1[0]*dxyz0[1]*dxyz1[2]*weight
    data[xyz0[0], xyz0[1], xyz1[2]] += dxyz1[0]*dxyz1[1]*dxyz0[2]*weight
    data[xyz1[0], xyz1[1], xyz0[2]] += dxyz0[0]*dxyz0[1]*dxyz1[2]*weight
    data[xyz1[0], xyz0[1], xyz1[2]] += dxyz0[0]*dxyz1[1]*dxyz0[2]*weight
    data[xyz0[0], xyz1[1], xyz1[2]] += dxyz1[0]*dxyz0[1]*dxyz0[2]*weight
    data[xyz1[0], xyz1[1], xyz1[2]] += dxyz0[0]*dxyz0[1]*dxyz0[2]*weight

def get_gaussian_kernel(resolution,voxel):
    sigma_factor=3.0
    sigmap=resolution/(2.0*voxel*3.0**0.5)
    exth=int(np.ceil(sigma_factor*sigmap))
    kernel_size=int(2*exth - 1)
    bvalue=-1.0/(2.0*sigmap*sigmap)
    cvalue=sigma_factor*sigma_factor*sigmap*sigmap
    gaussian_kernel = np.zeros([kernel_size]*3)
    for i in range(kernel_size):
        for j in range(kernel_size):
            for k in range(kernel_size):
                d2 = (i-exth+1)**2 + (j-exth+1)**2 + (k-exth+1)**2
                if d2 <= cvalue:
                    gaussian_kernel[i, j, k] = np.exp(d2*bvalue)
    return gaussian_kernel

def gaussian_filter(fit_map_mat,resolution,voxel):
    return scipy.ndimage.convolve(fit_map_mat,get_gaussian_kernel(resolution,voxel),mode='constant')

# calculate scores

def get_transformation_matrix_list(log_path):
    log_data=np.loadtxt(log_path,dtype=float,ndmin=2)
    transform_matrix_list=log_data[:,2:].reshape((-1,3,4))
    return transform_matrix_list

def transform_positions(xyzs,transformation_matrix):
    rotation_matrix=transformation_matrix[:,:3]
    translation_vector=transformation_matrix[:,-1]
    return xyzs.dot(rotation_matrix.T)+translation_vector

def get_fit_map_mat(u,transformation_matrix,density_param_dict):
    grid_shape=density_param_dict['grid_shape']
    origin=density_param_dict['origin']
    lengths=density_param_dict['lengths']
    voxel=density_param_dict['voxel']
    xyzs=u.atoms.positions
    n_atoms=u.atoms.n_atoms
    # apply the transformation matrix
    transformed_xyzs=transform_positions(xyzs,transformation_matrix)
    weights=[weight_dict[e] for e in u.atoms.types]
    data=np.zeros(grid_shape)
    for i in range(n_atoms):
        xyz=transformed_xyzs[i]
        weight=weights[i]
        if np.sum(xyz>origin)+np.sum(xyz<(origin+lengths-voxel-0.001))==6:
            add_atom(data,voxel,origin,xyz,weight)
    data=scipy.ndimage.convolve(data,get_gaussian_kernel(resolution,voxel),mode='constant')
    return data

def get_scores_overlap(u,transformation_matrix,ref_map_mat_laplacian,density_param_dict):
    fit_map_mat=get_fit_map_mat(u,transformation_matrix,density_param_dict)
    fit_map_mat_laplacian=scipy.ndimage.convolve(fit_map_mat,laplacian_kernel,mode='constant')
    sel=(ref_map_mat_laplacian<negtive_laplacian_cutoff)*(fit_map_mat_laplacian<fit_map_laplacian_cutoff_low)
    sel_ref=(ref_map_mat_laplacian<negtive_laplacian_cutoff)+(ref_map_mat_laplacian>positive_laplacian_cutoff)
    sel_fit=(fit_map_mat_laplacian<fit_map_laplacian_cutoff_low)+(fit_map_mat_laplacian>fit_map_laplacian_cutoff_high)
    sel2=sel_ref*sel_fit
    overlap_volume=np.sum(sel)
    if overlap_volume==0:
        return [0,0]
    overlap_correlation=np.dot(ref_map_mat_laplacian[sel2],fit_map_mat_laplacian[sel2])/(np.linalg.norm(ref_map_mat_laplacian[sel2])*np.linalg.norm(fit_map_mat_laplacian[sel2]))
    return [overlap_volume,overlap_correlation]


def score(parms):
    density=parms[0]
    ref_map_mat_laplacian=parms[1]
    log=parms[2]
    density_param_dict=parms[3]
    domain_name=log[:-4]
    domain_path=os.path.join(domain_dir,domain_name+'.pdb')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UserWarning)
        # 抑制读取pdb时的冗余警告
        u=mda.Universe(domain_path)     # original domain
    log_path=os.path.join(fitout_dir,density,"fitlogs",log)
    transformation_matrix_list=get_transformation_matrix_list(log_path)
    score_list=[]
    for transformation_matrix in transformation_matrix_list:
        score_list.append(get_scores_overlap(u,transformation_matrix,ref_map_mat_laplacian,density_param_dict))
    return [log[:-4],score_list]


if __name__=="__main__":
    os.makedirs(fitout_dir,exist_ok=True)

    map_list=[file_name for file_name in os.listdir(map_dir) if file_name.endswith(".mrc")]
    domain_list=[file_name for file_name in os.listdir(domain_dir) if file_name.endswith(".pdb")]

    for i, ref_map_filename in enumerate(map_list):
        # fit
        output_subdir=os.path.join(fitout_dir,ref_map_filename).replace("\\","/")
        os.makedirs(output_subdir,exist_ok=True)
        fitlog_subdir=os.path.join(output_subdir,"fitlogs")
        os.makedirs(fitlog_subdir,exist_ok=True)

        config_log_path=os.path.join(output_subdir,"fit_config.txt").replace("\\","/")
        config_log=open(config_log_path,'w')
        config_log.write(f"{domain_dir=}\n"+
                         f"{map_dir=}\n"+
                         f"{ref_map_threshold=}\n"+
                         f"{resolution=:.2f}\n"+
                         f"{n_search=}\n"+
                         f"{n_process=}")
        config_log.close()

        params_list=[(domain_filename,ref_map_filename) for domain_filename in domain_list]

        with multiprocessing.Pool(n_process) as pool:
            # 强制迭代tqdm对象，更新进度条
            for result in tqdm(pool.imap_unordered(fit,params_list),total=len(domain_list),desc=f"{i+1}/{len(map_list)}--Fitting {ref_map_filename}",file=sys.stdout):
                # 处理结果
                pass

        # score
        ref_map_path=os.path.join(map_dir,ref_map_filename).replace("\\","/")
        ref_map_mat_laplacian,density_param_dict=get_ref_map_mat_laplacian_and_params(ref_map_path)

        log_list=[log for log in os.listdir(fitlog_subdir) if log.endswith('.log')]

        pool_params=[(ref_map_filename,ref_map_mat_laplacian,log,density_param_dict) for log in log_list]

        pool=multiprocessing.Pool(n_process)
        scores=list(tqdm(pool.imap_unordered(score,pool_params),total=len(log_list),desc=f"{i+1}/{len(map_list)}--Locally scoring {ref_map_filename}",file=sys.stdout)) #1.将结果存入列表。2.list强制迭代tqdm对象，更新进度条
        pool.close()
        pool.join()

        config_path=os.path.join(fitout_dir,ref_map_filename,"overlap_score_config.txt")
        f=open(config_path,'w')
        f.write(f"{domain_dir=}\n"
               +f"{ref_map_path=}\n"
               +f"{ref_map_threshold=}\n"
               +f"{negtive_laplacian_cutoff=}\n"
               +f"{positive_laplacian_cutoff=}\n"
               +f"{resolution=}\n"
               +f"{fit_map_laplacian_cutoff_low=}\n"
               +f"{fit_map_laplacian_cutoff_high=}\n"
               +f"{n_process=}\n")
        f.close()


        scores_dict={domain_name:score_list for domain_name,score_list in scores}

        scores_path=os.path.join(fitout_dir,ref_map_filename,"overlap_scores.npy")
        np.save(scores_path,scores_dict,allow_pickle=True)

    print("Done fitting and scoring.")



