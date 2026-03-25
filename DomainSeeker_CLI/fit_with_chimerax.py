import os,sys,shutil,argparse,math
from dataclasses import dataclass, field
import subprocess,multiprocessing
import MDAnalysis as mda
import numpy as np
import scipy
import mrcfile
from tqdm import tqdm
import warnings

script_dir=os.path.dirname(os.path.realpath(__file__)).replace("\\","/")

python_exec=os.path.realpath(sys.executable)
chimerax_dir = os.path.dirname(python_exec)
# mac
if sys.platform == 'darwin':
    chimerax_dir = chimerax_dir.replace("/Contents/bin", "/Contents/MacOS")

# ChimeraX
# windows
if sys.platform in ["win32","win64"]:
    ChimeraX=os.path.join(chimerax_dir,"ChimeraX-console.exe")
# mac
elif sys.platform=="darwin":
    ChimeraX = os.path.join(chimerax_dir,"ChimeraX")
# linux
else:
    ChimeraX = os.path.join(chimerax_dir,"chimerax")


    
# 

def parser() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Perform domain-density fitting. The EM density map must first be segmented (manually or automatically) into regions of interest, each saved as a separate .mrc file",
                                    epilog="Use get_fitted_domains.py to obtain .pdb files of the fitted domains")
    parser.add_argument("domain_dir")
    parser.add_argument("map_dir")
    parser.add_argument("fitout_dir")
    parser.add_argument("ref_map_threshold",help="Electron density map threshold value")
    parser.add_argument("resolution",help="Resolution of the electron density map, used to generate simulated densities for individual domains")
    parser.add_argument("n_search",help="Number of initial starting positions for gradient-based correlation optimization. A value of 200 is recommended for typical domain sizes, balancing computational cost and accuracy. Larger density regions may require more starting positions for optimal fitting")
    parser.add_argument("negtive_laplacian_cutoff",help="Parameters for evaluating overlap region matching. Set these to cover the main structural regions of the density map, while excluding noisy areas")
    parser.add_argument("positive_laplacian_cutoff",help="Parameters for evaluating overlap region matching. Set these to cover the main structural regions of the density map, while excluding noisy areas")
    parser.add_argument("ncores",help="Number of parallel processes to use")
    parser.add_argument("--batch-size", default = 10)
    parser.add_argument("--chimerax", default = False, help="Path to ChimeraX if it cannot be determined automatically", required=False)
    parser.add_argument("--calculation", choices=['default', 'chimerax-batch', 'direct-batch'])
    parser.add_argument("--max-domain", default = -1)

    return parser.parse_args()

# generate density matrix
WEIGHT_DICT = {
    'C':12.011,
    'N':14.007,
    'O':16.000,
    'S':32.060
    }

# Laplacian kernel
LAPLACIAN_KERNEL=np.zeros([3,3,3])
LAPLACIAN_KERNEL[1,1,1]=-6
LAPLACIAN_KERNEL[1,1,0]=1
LAPLACIAN_KERNEL[1,1,2]=1
LAPLACIAN_KERNEL[1,0,1]=1
LAPLACIAN_KERNEL[1,2,1]=1
LAPLACIAN_KERNEL[0,1,1]=1
LAPLACIAN_KERNEL[2,1,1]=1


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

def fit(
        ref_map_filename: str,
        map_dir: str,
        domain_dir: str,
        fitout_dir: str,
        ref_map_threshold: float,
        resolution: float,
        n_search: int,
        domain_filename_list: list[str]
):
    map_path=os.path.join(map_dir,ref_map_filename).replace("\\","/")
    output_subdir=os.path.join(fitout_dir,ref_map_filename).replace("\\","/")
    domain_filename_list = [os.path.join(domain_dir,domain_filename).replace("\\","/") for domain_filename in domain_filename_list]

    # Default command
    # cmd_list = [ChimeraX, "--nogui", "--script", f"{script_dir}/fit_in_chimerax.py {domain_path} {map_path} {output_subdir} {ref_map_threshold} {resolution} {n_search}", "--exit"]
    # subprocess.run(cmd_list,shell=False,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)

    # Use ChimeraX, but keep session open for batch
    # cmd_list = [ChimeraX, "--nogui", "--script", f"{script_dir}/fit_in_chimerax_pac.py {map_path} {output_subdir} {ref_map_threshold} {resolution} {n_search} {" ".join(domain_filename_list)}", "--exit"]
    # subprocess.run(cmd_list,shell=False,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)

    # Use ChimeraX in python script directly
    cmd_list = [
            '/usr/lib/ucsf-chimerax-daily/bin/python3.9', 
            '/cds/ban/ekummerant/domainseeker/DomainSeeker/DomainSeeker_CLI/fit_chimerax_directly.py', 
            map_path, output_subdir, 
            str(ref_map_threshold), 
            str(resolution), 
            str(n_search)
        ]
    [cmd_list.append(domain_path) for domain_path in domain_filename_list]
    subprocess.run(cmd_list,shell=False, stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)

def fitMultiProc(params):
    fit(*params)

# calculate scores

def get_transformation_matrix_list(log_path):
    log_data=np.loadtxt(log_path,dtype=float,ndmin=2)
    transform_matrix_list=log_data[:,2:].reshape((-1,3,4))
    return transform_matrix_list

def transform_positions(xyzs,transformation_matrix):
    rotation_matrix=transformation_matrix[:,:3]
    translation_vector=transformation_matrix[:,-1]
    return xyzs.dot(rotation_matrix.T)+translation_vector

def get_ref_map_mat_laplacian_and_params(ref_map_path, ref_map_threshold):
    ref_map=mrcfile.open(ref_map_path)
    ref_map_mat=ref_map.data.T
    ref_map_mat=ref_map_mat*(ref_map_mat>=ref_map_threshold)
    origin=np.array(ref_map.header.origin.tolist())
    lengths=np.array(ref_map.header.cella.tolist())
    grid_shape=np.array([ref_map.header.nx,ref_map.header.ny,ref_map.header.nz])
    voxel=lengths[0]/grid_shape[0]
    density_param_dict={'origin':origin,'lengths':lengths,'grid_shape':grid_shape,'voxel':voxel}
    ref_map_mat_laplacian=scipy.ndimage.convolve(ref_map_mat,LAPLACIAN_KERNEL,mode='constant')
    return [ref_map_mat_laplacian,density_param_dict]

def get_fit_map_mat(
        u,
        transformation_matrix,
        density_param_dict,
        resolution
    ):
    grid_shape=density_param_dict['grid_shape']
    origin=density_param_dict['origin']
    lengths=density_param_dict['lengths']
    voxel=density_param_dict['voxel']
    xyzs=u.atoms.positions
    n_atoms=u.atoms.n_atoms
    # apply the transformation matrix
    transformed_xyzs=transform_positions(xyzs,transformation_matrix)
    weights=[WEIGHT_DICT[e] for e in u.atoms.types]
    data=np.zeros(grid_shape)
    for i in range(n_atoms):
        xyz=transformed_xyzs[i]
        weight=weights[i]
        if np.sum(xyz>origin)+np.sum(xyz<(origin+lengths-voxel-0.001))==6:
            add_atom(data,voxel,origin,xyz,weight)
    data=scipy.ndimage.convolve(data,get_gaussian_kernel(resolution,voxel),mode='constant')
    return data

def get_scores_overlap(
        u,
        transformation_matrix,
        ref_map_mat_laplacian,
        density_param_dict,
        negtive_laplacian_cutoff,
        positive_laplacian_cutoff,
        fit_map_laplacian_cutoff_low,
        fit_map_laplacian_cutoff_high,
        resolution
    ):
    fit_map_mat=get_fit_map_mat(u,transformation_matrix,density_param_dict,resolution)
    fit_map_mat_laplacian=scipy.ndimage.convolve(fit_map_mat,LAPLACIAN_KERNEL,mode='constant')
    sel=(ref_map_mat_laplacian<negtive_laplacian_cutoff)*(fit_map_mat_laplacian<fit_map_laplacian_cutoff_low)
    sel_ref=(ref_map_mat_laplacian<negtive_laplacian_cutoff)+(ref_map_mat_laplacian>positive_laplacian_cutoff)
    sel_fit=(fit_map_mat_laplacian<fit_map_laplacian_cutoff_low)+(fit_map_mat_laplacian>fit_map_laplacian_cutoff_high)
    sel2=sel_ref*sel_fit
    overlap_volume=np.sum(sel)
    if overlap_volume==0:
        return [0,0]
    overlap_correlation=np.dot(ref_map_mat_laplacian[sel2],fit_map_mat_laplacian[sel2])/(np.linalg.norm(ref_map_mat_laplacian[sel2])*np.linalg.norm(fit_map_mat_laplacian[sel2]))
    return [overlap_volume,overlap_correlation]


def score(
    density,
    ref_map_mat_laplacian,
    log,
    density_param_dict,
    domain_name,
    domain_dir,
    fitout_dir,
    negtive_laplacian_cutoff,
    positive_laplacian_cutoff,
    fit_map_laplacian_cutoff_low,
    fit_map_laplacian_cutoff_high,
    resolution
):
    domain_path=os.path.join(domain_dir,domain_name+'.pdb')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UserWarning)
        # 
        u=mda.Universe(domain_path)     # original domain
    log_path=os.path.join(fitout_dir,density,"fitlogs",log)
    transformation_matrix_list=get_transformation_matrix_list(log_path)
    score_list=[]
    for transformation_matrix in transformation_matrix_list:
        score_list.append(
            get_scores_overlap(
                u,
                transformation_matrix,
                ref_map_mat_laplacian,
                density_param_dict,
                negtive_laplacian_cutoff,
                positive_laplacian_cutoff,
                fit_map_laplacian_cutoff_low,
                fit_map_laplacian_cutoff_high,
                resolution
            )
        )
    return [domain_name, score_list]

def scoreMultiProc(params):
    return score(*params)


def get_batch(
        domains: list[str], 
        batch_size: int, 
        n_procs: int, 
        max_domains: int = -1
    ):
    '''
    Returns batches of size batch_size, except if there are fewer domains
    available than are needed to fill all available processors. In that case,
    the batch_size is reduced to fill all processors
    '''
    if max_domains != -1:
        domains = domains[:max_domains]
    if len(domains) < (batch_size * n_procs):
        print('condition true')
        batch_size = math.ceil(len(domains)/n_procs)
    for i in range(0, len(domains), batch_size):
        yield domains[i:(i+batch_size)]

if __name__=="__main__":
    argument = parser()

    domain_dir = argument.domain_dir
    map_dir = argument.map_dir
    fitout_dir = argument.fitout_dir
    ref_map_threshold = float(argument.ref_map_threshold)
    resolution = float(argument.resolution)
    n_search = int(argument.n_search)
    negtive_laplacian_cutoff = float(argument.negtive_laplacian_cutoff)
    positive_laplacian_cutoff = float(argument.positive_laplacian_cutoff)
    n_process=int(argument.ncores)
    if argument.chimerax:
        ChimeraX = argument.chimerax

    fit_map_laplacian_cutoff_low=-2
    fit_map_laplacian_cutoff_high=15

    os.makedirs(fitout_dir,exist_ok=True)
    map_list=[file_name for file_name in os.listdir(map_dir) if file_name.endswith(".mrc")]
    domain_list=[file_name for file_name in os.listdir(domain_dir) if file_name.endswith(".pdb")]


    for i, ref_map_filename in enumerate(map_list):
        # fit
        output_subdir=os.path.join(fitout_dir,ref_map_filename).replace("\\","/")
        os.makedirs(output_subdir,exist_ok=True)
        fitlog_subdir=os.path.join(output_subdir,"fitlogs")
        os.makedirs(fitlog_subdir,exist_ok=True)

        seen_domains = [x.removesuffix('.log') for x in os.listdir(fitlog_subdir) if os.path.isfile(os.path.join(fitlog_subdir, x)) and x.endswith('.log')]
        filtered_domain_list = [x for x in domain_list if x.removesuffix('.pdb') not in seen_domains]

        print(f'Total number of domains: {len(domain_list)}')
        print(f'Number of domains seen: {len(seen_domains)} for map {ref_map_filename}')
        print(f'Remaining number of domains: {len(filtered_domain_list)}')


        config_log_path=os.path.join(output_subdir,"fit_config.txt").replace("\\","/")
        config_log=open(config_log_path,'w')
        config_log.write(f"{domain_dir=}\n"+
                         f"{map_dir=}\n"+
                         f"{ref_map_threshold=}\n"+
                         f"{resolution=:.2f}\n"+
                         f"{n_search=}\n"+
                         f"{n_process=}")
        config_log.close()

        params_list=[(ref_map_filename, map_dir, domain_dir, fitout_dir, ref_map_threshold, resolution, n_search, domain_filename_batch) for domain_filename_batch in get_batch(filtered_domain_list, int(argument.batch_size), n_process, int(argument.max_domain))]

        with multiprocessing.Pool(n_process) as pool:
            for result in tqdm(pool.imap_unordered(fitMultiProc,
                                                   params_list),
                               total=math.ceil(len(filtered_domain_list)/int(argument.batch_size)),
                               desc=f"{i+1}/{len(map_list)}--Fitting {ref_map_filename}",
                               file=sys.stdout,
                               smoothing=0.2):
                pass

        # score
        ref_map_path=os.path.join(map_dir,ref_map_filename).replace("\\","/")
        ref_map_mat_laplacian,density_param_dict=get_ref_map_mat_laplacian_and_params(ref_map_path, ref_map_threshold)

        log_list=[log for log in os.listdir(fitlog_subdir) if log.endswith('.log')]

        pool_params=[(ref_map_filename, ref_map_mat_laplacian, log, density_param_dict, log[:-4], domain_dir, fitout_dir, negtive_laplacian_cutoff, positive_laplacian_cutoff, fit_map_laplacian_cutoff_low, fit_map_laplacian_cutoff_high, resolution) for log in log_list]

        pool=multiprocessing.Pool(n_process)
        scores=list(tqdm(pool.imap_unordered(scoreMultiProc,pool_params),total=len(log_list),desc=f"{i+1}/{len(map_list)}--Locally scoring {ref_map_filename}",file=sys.stdout)) #1.将结果存入列表。2.list强制迭代tqdm对象，更新进度条
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



