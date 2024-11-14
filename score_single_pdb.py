import MDAnalysis as mda
import numpy as np
import scipy
import mrcfile
import os,sys

domain_path=os.path.realpath(sys.argv[1])
ref_map_path=os.path.realpath(sys.argv[2])
ref_map_threshold=float(sys.argv[3])
ref_map_laplacian_cutoff_low=float(sys.argv[4])
ref_map_laplacian_cutoff_high=float(sys.argv[5])
resolution=float(sys.argv[6])
fit_map_laplacian_cutoff_low=-2
fit_map_laplacian_cutoff_high=15

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

def get_fit_map_mat(domain_path,origin,lengths,grid_shape,voxel,resolution):
    u=mda.Universe(domain_path)
    n_atoms=u.atoms.n_atoms
    xyzs=u.atoms.positions
    weights=[weight_dict[e] for e in u.atoms.types]
    data=np.zeros(grid_shape)
    for i in range(n_atoms):
        xyz=xyzs[i]
        weight=weights[i]
        if np.sum(xyz>origin)+np.sum(xyz<(origin+lengths-voxel-0.001))==6:
            add_atom(data,voxel,origin,xyz,weight)
    data=scipy.ndimage.convolve(data,get_gaussian_kernel(resolution,voxel),mode='constant')
    return data

def get_scores_overlap(domain_path):
    fit_map_mat=get_fit_map_mat(domain_path,origin,lengths,grid_shape,voxel,resolution)
    fit_map_mat_laplacian=scipy.ndimage.convolve(fit_map_mat,laplacian_kernel,mode='constant')
    sel=(ref_map_mat_laplacian<ref_map_laplacian_cutoff_low)*(fit_map_mat_laplacian<fit_map_laplacian_cutoff_low)
    sel_ref=(ref_map_mat_laplacian<ref_map_laplacian_cutoff_low)+(ref_map_mat_laplacian>ref_map_laplacian_cutoff_high)
    sel_fit=(fit_map_mat_laplacian<fit_map_laplacian_cutoff_low)+(fit_map_mat_laplacian>fit_map_laplacian_cutoff_high)
    sel2=sel_ref*sel_fit
    overlap_volume=np.sum(sel)
    overlap_correlation=np.dot(ref_map_mat_laplacian[sel2],fit_map_mat_laplacian[sel2])/(np.linalg.norm(ref_map_mat_laplacian[sel2])*np.linalg.norm(fit_map_mat_laplacian[sel2]))
    return [overlap_volume,overlap_correlation]



ref_map=mrcfile.open(ref_map_path)
ref_map_mat=ref_map.data.T
ref_map_mat=ref_map_mat*(ref_map_mat>=ref_map_threshold)
origin=np.array(ref_map.header.origin.tolist())
lengths=np.array(ref_map.header.cella.tolist())
grid_shape=np.array([ref_map.header.nx,ref_map.header.ny,ref_map.header.nz])
voxel=lengths[0]/grid_shape[0]
ref_map_mat_laplacian=scipy.ndimage.convolve(ref_map_mat,laplacian_kernel,mode='constant')

print(get_scores_overlap(domain_path))



