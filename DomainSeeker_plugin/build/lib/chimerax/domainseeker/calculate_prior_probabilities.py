# 全局异常处理
import domainseeker_errorLog
#----------------------------------------------------------------------------------------------------

import os, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import math
from scipy.interpolate import interp1d

map_dir = sys.argv[2]
fitout_dir = sys.argv[3]
box_num = int(sys.argv[4])
min_entries_per_box = int(sys.argv[5])
relative_density_cutoff = float(sys.argv[6])
z_score_offset = float(sys.argv[7])

grids_x=box_num
grids_y=box_num


def get_grid_index(point,min_vol,max_vol,min_cor,max_cor,grids_x,grids_y):
    step_x=(max_vol-min_vol)/grids_x
    step_y=(max_cor-min_cor)/grids_y
    vol=point[0]
    cor=point[1]
    index_x=math.floor((vol-min_vol)/step_x)
    if vol==max_vol:
            index_x=index_x-1
    index_y=math.floor((cor-min_cor)/step_y)
    if cor==max_cor:
            index_y=index_y-1
    index=[index_x,index_y]
    return index

def assign_points_to_grids(data,grids_x,grids_y):
    min_vol=np.min(data[:,0])
    max_vol=np.max(data[:,0])
    min_cor=np.min(data[:,1])
    max_cor=np.max(data[:,1])
    grids=[[[] for j in range(grids_y)] for i in range(grids_x)]
    for point in data:
        index=get_grid_index(point,min_vol,max_vol,min_cor,max_cor,grids_x,grids_y)
        grids[index[0]][index[1]].append(point)
    return grids

def get_filtered_grids(grids,relative_density_cutoff,number_of_data):
    grids_x=len(grids)
    grids_y=len(grids[0])
    filtered_grids=[[[] for j in range(grids_y)] for i in range(grids_x)]
    # remove points in grids with low relative density
    for i in range(grids_x):
          for j in range(grids_y):
            relative_density=len(grids[i][j])*grids_x*grids_y/number_of_data
            if relative_density>=relative_density_cutoff:
                 filtered_grids[i][j]=grids[i][j]
    return filtered_grids


def save_relative_grids_density_plot(grids,number_of_data, density_filename):
    colors = [(0,"white"),(0.1,"red"),(1,"red")]
    cmap = LinearSegmentedColormap.from_list('custom_cmap', colors)

    grids_x=len(grids)
    grids_y=len(grids[0])
    relative_densitys=[[0 for grid in row] for row in grids]
    for i in range(grids_x):
        for j in range(grids_y):
            relative_density=len(grids[i][j])*grids_x*grids_y/number_of_data
            relative_densitys[i][j]=relative_density
    relative_densitys=np.array(relative_densitys)
    plt.figure(figsize=(10,10))
    plt.matshow(relative_densitys.T[::-1], cmap=cmap)
    plt.colorbar()
    save_path=os.path.join(fitout_dir, density_filename, "local_assessing_relative_density.png")
    plt.savefig(save_path)

def merge_grids_to_box(grids):
     box=[]
     for grid in grids:
          box+=grid
     return box

def ave_std(box,colume):
    data=[]
    for row in box:
        data.append(row[colume])
    return [np.mean(data),np.std(data)]

def get_zScores(scores, density_filename):
    data=scores[:,1:].astype(float)
    max_vol=np.max(data[:,0])
    min_vol=np.min(data[:,0])
    box_step=(max_vol-min_vol)/box_num
    # assign points to grids
    grids=assign_points_to_grids(data,grids_x,grids_y)
    # remove anomalous points
    filtered_grids=get_filtered_grids(grids,relative_density_cutoff,len(data))
    save_relative_grids_density_plot(filtered_grids, len(data), density_filename)
    # assign all points by volume into boxes
    boxes=[merge_grids_to_box(filtered_grids[i]) for i in range(box_num)]
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
    ave_std_list=[ave_std(box,1) if len(box)>0 else [0,0] for box in boxes]
    for i in range(box_num):
        if len(boxes[i])==0:
            # find the left neighbour and right neighbour of the empty box
            l=i-1
            while l>=0 and len(boxes[l])==0:
                l=l-1
            r=i+1
            while r<box_num and len(boxes[r])==0:
                r=r+1
            # if box i the leftmost
            if l==-1:
                ave_std_list[i]=ave_std_list[r]
            # if box i the rightmost
            elif r==box_num:
                ave_std_list[i]=ave_std_list[l]
            else:
                ave_std_list[i]=[np.interp(i,[l,r],[ave_std_list[l][0],ave_std_list[r][0]]),np.interp(i,[l,r],[ave_std_list[l][1],ave_std_list[r][1]])]
    ave_std_list=np.array([[ave_std_list[0][0],ave_std_list[0][1]]]
                          +ave_std_list
                          +[[ave_std_list[box_num-1][0],ave_std_list[box_num-1][1]]])
    
    # get the interpolation function of ave and std
    f_ave=interp1d(np.linspace(min_vol-0.5*box_step,max_vol+0.5*box_step,box_num+2),ave_std_list[:,0],3)
    f_std=interp1d(np.linspace(min_vol-0.5*box_step,max_vol+0.5*box_step,box_num+2),ave_std_list[:,1],3)
    # show data
    plt.figure(figsize=(10,10))
    plt.scatter(data[:,0],data[:,1],s=0.5)
    plt.plot(np.linspace(min_vol-0.5*box_step,max_vol+0.5*box_step,box_num+2),ave_std_list[:,0],c='black')
    plt.plot(np.linspace(min_vol-0.5*box_step,max_vol+0.5*box_step,box_num+2),ave_std_list[:,0]+2*ave_std_list[:,1],c='black',linestyle='--')
    plot_save_path=os.path.join(fitout_dir, density_filename, "local_assessing.png")
    plt.savefig(plot_save_path)
    # calculate z-scores
    zScores=[]
    for row in scores:
            vol=int(row[1])
            cor=float(row[2])
            zScore=(cor-f_ave(vol))/f_std(vol)
            zScores.append(row.tolist()+[zScore])
    zScores.sort(key=lambda row:-row[-1])
    # save z-scores
    zScores_save_path=os.path.join(fitout_dir, density_filename, "zScores.txt")
    np.savetxt(zScores_save_path,zScores,fmt='%s')
    return zScores

def correlation_to_probability(c):
    return (c/(1-c))**3

def get_correct_fitting_probabilities(fit_log_data):
    fit_log_data=np.array(fit_log_data)
    correlations=fit_log_data[:,0]
    hits=fit_log_data[:,1]
    probabilities=correlation_to_probability(correlations)*hits
    total_probability=np.sum(probabilities)
    probabilities/=total_probability
    return probabilities

# sigmoid function
def zScore_to_probability(z, offset):
    return 1 / (1 + np.exp(offset - z))

def get_prior_probabilities(fitting_probabilities, zScores, offset):
    # 提取键并预分配数组以提高效率
    states=np.array([row[0] for row in zScores])
    factor_array = np.empty(len(states))
    
    # 计算因子值
    for i, name in enumerate(states):
        proteinId, domainId, fitId = name.split('_')
        z = float(zScores[i][-1])
        factor_array[i] = zScore_to_probability(z, offset) * fitting_probabilities[proteinId + "_" + domainId][int(fitId)]
    
    # 计算总和和排序索引
    partition_function = np.sum(factor_array)
    order_list = np.argsort(factor_array)[::-1]
    
    # 构造结果字典
    prior_probabilities = [[states[i], factor_array[i] / partition_function] for i in order_list]
    return prior_probabilities


density_filenames = os.listdir(map_dir)

for i, density_filename in enumerate(density_filenames):
    fitout_subdir = os.path.join(fitout_dir, density_filename)
    fitlog_subdir = os.path.join(fitout_subdir, "fitlogs")

    # 计算fitting_probability
    # 记录参数
    prior_config_path=os.path.join(fitout_subdir,"prior_config.txt")
    with open(prior_config_path,"w") as f:
        f.write("Parameter for calculating z-scores:\n")
        f.write(f"box_num={box_num}\n")
        f.write(f"min_entries_per_box={min_entries_per_box}\n")
        f.write(f"relative_density_cutoff={relative_density_cutoff}\n")
        f.write("\n")
        f.write("Parameter for calculating prior probabilities:\n")
        f.write(f"z_score_offset={z_score_offset}\n")

    fitting_probabilities_path=os.path.join(fitout_subdir,"fitting_probabilities.npy")
    fitting_probabilities={}
    print(f"{i+1}/{len(density_filenames)}--Calculating fitting probabilities for {density_filename}...")
    for file_name in os.listdir(fitlog_subdir):
        if file_name.endswith(".log"):
            fit_log_path=os.path.join(fitlog_subdir,file_name)
            if os.path.getsize(fit_log_path)>0:
                fit_log_data=np.loadtxt(fit_log_path,dtype=float,usecols=(0,1),ndmin=2)
                fitting_probabilities[file_name[:-4]]=get_correct_fitting_probabilities(fit_log_data)
    np.save(fitting_probabilities_path,fitting_probabilities)

    
    # 计算z-scores
    print(f"{i+1}/{len(density_filenames)}--Calculating z-scores for {density_filename}...")
    # 读取scores
    scores_path = os.path.join(fitout_subdir, "overlap_scores.npy")
    if os.path.exists(scores_path):
        data=np.load(scores_path,allow_pickle=True).item()
        scores=[]
        for domain in data.keys():
            for fit_id, score in enumerate(data[domain]):
                if score[0]!=0:
                    scores.append([f"{domain}_{fit_id}",score[0],score[1]])
        scores=np.array(scores)
        # 计算z-scores
        zScores=get_zScores(scores, density_filename)
    else:
        print(f"No scores found for {density_filename}",file=sys.stderr)

    # 计算prior_probabilities
    print(f"{i+1}/{len(density_filenames)}--Calculating prior probabilities for {density_filename}...")
    prior_probabilities = get_prior_probabilities(fitting_probabilities, zScores, z_score_offset)
    prior_probabilities_path = os.path.join(fitout_subdir, "prior_probabilities.txt")
    np.savetxt(prior_probabilities_path, prior_probabilities,fmt='%s')

print("Done calculating prior probabilities.")


    
    