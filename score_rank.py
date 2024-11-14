import MDAnalysis as mda
import numpy as np
import math
from scipy.interpolate import interp1d
import scipy
import mrcfile
import os,sys
import multiprocessing

domain_dir=os.path.realpath(sys.argv[1])
map_dir=os.path.realpath(sys.argv[2])
fit_out_dir=sys.argv[3]
ref_map_threshold=float(sys.argv[4])
ref_map_laplacian_cutoff_low=float(sys.argv[5])
ref_map_laplacian_cutoff_high=float(sys.argv[6])
resolution=float(sys.argv[7])
fit_map_laplacian_cutoff_low=-2
fit_map_laplacian_cutoff_high=15
n_thread=int(sys.argv[8])
box_num=int(sys.argv[9])
min_entries_per_box=int(sys.argv[10])

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

def get_fit_map_mat(domain_path,log_path,origin,lengths,grid_shape,voxel,resolution):
    u=mda.Universe(domain_path)
    n_atoms=u.atoms.n_atoms
    xyzs=u.atoms.positions
    # get and apply the transformation matrix
    transformation_matrix=np.loadtxt(log_path,skiprows=5,max_rows=3,dtype=float)
    rotation_matrix=transformation_matrix[:,:3]
    xyzs=xyzs.dot(rotation_matrix.T)
    translation_vector=transformation_matrix[:,-1]
    xyzs+=translation_vector
    weights=[weight_dict[e] for e in u.atoms.types]
    data=np.zeros(grid_shape)
    for i in range(n_atoms):
        xyz=xyzs[i]
        weight=weights[i]
        if np.sum(xyz>origin)+np.sum(xyz<(origin+lengths-voxel-0.001))==6:
            add_atom(data,voxel,origin,xyz,weight)
    data=scipy.ndimage.convolve(data,get_gaussian_kernel(resolution,voxel),mode='constant')
    return data

def get_scores_overlap(domain_path,log_path):
    fit_map_mat=get_fit_map_mat(domain_path,log_path,origin,lengths,grid_shape,voxel,resolution)
    fit_map_mat_laplacian=scipy.ndimage.convolve(fit_map_mat,laplacian_kernel,mode='constant')
    sel=(ref_map_mat_laplacian<ref_map_laplacian_cutoff_low)*(fit_map_mat_laplacian<fit_map_laplacian_cutoff_low)
    sel_ref=(ref_map_mat_laplacian<ref_map_laplacian_cutoff_low)+(ref_map_mat_laplacian>ref_map_laplacian_cutoff_high)
    sel_fit=(fit_map_mat_laplacian<fit_map_laplacian_cutoff_low)+(fit_map_mat_laplacian>fit_map_laplacian_cutoff_high)
    sel2=sel_ref*sel_fit
    overlap_volume=np.sum(sel)
    overlap_correlation=np.dot(ref_map_mat_laplacian[sel2],fit_map_mat_laplacian[sel2])/(np.linalg.norm(ref_map_mat_laplacian[sel2])*np.linalg.norm(fit_map_mat_laplacian[sel2]))
    return [overlap_volume,overlap_correlation]


def score(log):
    domain_name=log[:-4]
    domain_path=os.path.join(domain_dir,domain_name+'.pdb')
    log_path=os.path.join(fit_out_dir,density,log)
    return [domain_name]+get_scores_overlap(domain_path,log_path)
    
def ave_std(box,colume):
    data=[]
    for row in box:
        data.append(row[colume])
    return [np.mean(data),np.std(data)]
    
def get_zScores(scores):
    # sort and classity all scores by volume into boxes
    scores.sort(key=lambda row:row[1])
    min_vol=scores[0][1]
    max_vol=scores[-1][1]
    volume_len=max_vol-min_vol
    box_step=math.ceil(volume_len/box_num);
    boxes=[[] for i in range(box_num)]
    for i in range(len(scores)):
        if scores[i][1]!=max_vol:
            index=math.floor((scores[i][1]-min_vol)/box_step)
        else:
            index=box_num-1
        boxes[index].append(scores[i])
    # merge data in small boxes into their neighbours
    biggest_box_index=np.argsort(list(map(len,boxes)))[-1]
    for i in range(biggest_box_index-1):
        if len(boxes[i])<min_entries_per_box:
            boxes[i+1]=boxes[i]+boxes[i+1]
            boxes[i]=[]
    for i in range(box_num-1,biggest_box_index,-1):
        if len(boxes[i])<min_entries_per_box:
            boxes[i-1]=boxes[i-1]+boxes[i]
            boxes[i]=[]
    # calculate the ave and the std of each boxes
    ave_std_list=[ave_std(box,2) if len(box)>0 else [0,0] for box in boxes]
    for i in range(box_num):
        if len(boxes[i])==0:
            l=i-1
            while l>=0 and len(boxes[l])==0:
                l=l-1
            r=i+1
            while r<box_num and len(boxes[r])==0:
                r=r+1
            if l==-1:
                ave_std_list[i]=ave_std_list[r]
            elif r==box_num:
                ave_std_list[i]=ave_std_list[l]
            else:
                ave_std_list[i]=[np.interp(i,[l,r],[ave_std_list[l][0],ave_std_list[r][0]]),np.interp(i,[l,r],[ave_std_list[l][1],ave_std_list[r][1]])]
    ave_std_list=np.array([[2*ave_std_list[0][0]-ave_std_list[1][0],2*ave_std_list[0][1]-ave_std_list[1][1]]]
                          +ave_std_list
                          +[[2*ave_std_list[box_num-1][0]-ave_std_list[box_num-2][0],2*ave_std_list[box_num-1][1]-ave_std_list[box_num-2][1]]])
    # get the interpolation function of ave and std
    f_ave=interp1d(np.linspace(min_vol-0.5*box_step,max_vol+0.5*box_step,box_num+2),ave_std_list[:,0],3)
    f_std=interp1d(np.linspace(min_vol-0.5*box_step,max_vol+0.5*box_step,box_num+2),ave_std_list[:,1],3)
    # calculate z-scores
    zScores=list(map(lambda row:row+[(row[2]-f_ave(row[1]))/f_std(row[1])],scores))
    zScores.sort(key=lambda row:-row[-1])
    return zScores


if __name__=="__main__":
    density_list=os.listdir(fit_out_dir)
    for density in density_list:
        if density.endswith(".mrc"):
            ref_map_path=os.path.join(map_dir,density)
            ref_map=mrcfile.open(ref_map_path)
            ref_map_mat=ref_map.data.T
            ref_map_mat=ref_map_mat*(ref_map_mat>=ref_map_threshold)
            origin=np.array(ref_map.header.origin.tolist())
            lengths=np.array(ref_map.header.cella.tolist())
            grid_shape=np.array([ref_map.header.nx,ref_map.header.ny,ref_map.header.nz])
            voxel=lengths[0]/grid_shape[0]
            ref_map_mat_laplacian=scipy.ndimage.convolve(ref_map_mat,laplacian_kernel,mode='constant')
            
            log_list=[log for log in os.listdir(os.path.join(fit_out_dir,density)) if log.endswith('.log')]
            
            pool=multiprocessing.Pool(n_thread)
            scores=pool.map(score,log_list)
            pool.close()
            
            config_path=os.path.join(fit_out_dir,density,"score_rank_config.txt")
            f=open(config_path,'w')
            f.write(f"{domain_dir=}\n"
                   +f"{ref_map_path=}\n"
                   +f"{ref_map_threshold=}\n"
                   +f"{ref_map_laplacian_cutoff_low=}\n"
                   +f"{ref_map_laplacian_cutoff_high=}\n"
                   +f"{resolution=}\n"
                   +f"{fit_map_laplacian_cutoff_low=}\n"
                   +f"{fit_map_laplacian_cutoff_high=}\n"
                   +f"{box_num=}\n"
                   +f"{min_entries_per_box=}\n")
            f.close()
            
            zScores=get_zScores(scores)
            scores_path=os.path.join(fit_out_dir,density,"scores.txt")
            np.savetxt(scores_path,zScores,fmt="%s")