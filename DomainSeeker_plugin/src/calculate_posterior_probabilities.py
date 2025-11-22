# 全局异常处理
import domainseeker_errorLog
#----------------------------------------------------------------------------------------------------

import numpy as np
import MDAnalysis as mda
import os,sys
import matplotlib.pyplot as plt
import networkx as nx
from scipy.spatial.transform import Rotation
import itertools
from collections import defaultdict
from tqdm import tqdm
import mrcfile
from scipy.spatial.distance import cdist

project_dir=sys.argv[2]
origin_domain_dir=sys.argv[3]
map_dir=sys.argv[4]
map_level=float(sys.argv[5])
fitout_dir=sys.argv[6]

acceptor_prior_probability_cutoff=float(sys.argv[7])
donor_prior_probability_cutoff=float(sys.argv[8])
evidence_strenth=float(sys.argv[9])

crosslink_files = sys.argv[10:]



def read_prior_probabilities(file_path):
    with open(file_path,"r") as f:
        lines = f.readlines()
        states=[0]*len(lines)
        prior_probabilities=[0]*len(lines)
        for i,line in enumerate(lines):
            state,prob=line.strip().split()
            states[i]=state
            prior_probabilities[i]=float(prob)
    return states,prior_probabilities

density_names=[file_name[:-4] for file_name in os.listdir(map_dir) if file_name.endswith('.mrc')]
density_names.sort()
density_name_to_index = {name:i for i,name in enumerate(density_names)}
n_densities=len(density_names)


states_of_densities=[]
state_to_id_of_densities=[]
prior_probabilities_of_densities=[]
for density_name in density_names:
    prior_probabilities_file_path=os.path.join(fitout_dir,density_name+".mrc","prior_probabilities.txt")
    states,prior_probabilities=read_prior_probabilities(prior_probabilities_file_path)
    states_of_densities.append(states)
    state_to_id_of_densities.append(dict(zip(states,range(1,len(states)+1))))
    prior_probabilities_of_densities.append(prior_probabilities)

# read crosslinks
crosslinked_residues_graph=nx.Graph()
crosslink_list=[]
for crosslink_file in crosslink_files:
    crosslink_list+=[item for item in np.loadtxt(crosslink_file,dtype=str,ndmin=2) if len(set(item))==2]
crosslinked_residues_graph.add_edges_from(crosslink_list)
print("number of unique crosslinks in files:",crosslinked_residues_graph.number_of_edges())

def is_residueId_in_selection(residueId,selection):
    for item in selection.split(','):
        start=int(item.split('-')[0])
        end=int(item.split('-')[1])
        if start<=int(residueId)<=end:
            return True
    return False

protein_to_crosslinked_residues_to_domains=defaultdict(dict)

for node in crosslinked_residues_graph.nodes:
    uniprot_id, residue_id = node.split(':')
    protein_to_crosslinked_residues_to_domains[uniprot_id][residue_id]=[]

for uniprot_id in protein_to_crosslinked_residues_to_domains.keys():
    domains_info_path=os.path.join(origin_domain_dir,uniprot_id+'.domains')
    if os.path.exists(domains_info_path):
        domains_info=np.loadtxt(domains_info_path,dtype=str).reshape(-1,2)
        for residue_id in protein_to_crosslinked_residues_to_domains[uniprot_id].keys():
            related_domains=protein_to_crosslinked_residues_to_domains[uniprot_id][residue_id]
            for domain_id,selection in domains_info:
                if is_residueId_in_selection(residue_id,selection):
                    related_domains.append(domain_id)

crosslink_nodes_without_domains=set()
# 删除没有对应结构域的残基
for uniprot_id in protein_to_crosslinked_residues_to_domains.keys():
    residues_without_domains = [residue_id for residue_id in protein_to_crosslinked_residues_to_domains[uniprot_id].keys() if len(protein_to_crosslinked_residues_to_domains[uniprot_id][residue_id])==0]
    for residue_id in residues_without_domains:
        protein_to_crosslinked_residues_to_domains[uniprot_id].pop(residue_id)
        crosslink_nodes_without_domains.add(uniprot_id+':'+residue_id)

# 删除没有对应结构域的蛋白质
crosslink_related_proteins_without_domains=[uniprot_id for uniprot_id in protein_to_crosslinked_residues_to_domains.keys() if len(protein_to_crosslinked_residues_to_domains[uniprot_id])==0]
for uniprot_id in crosslink_related_proteins_without_domains:
    protein_to_crosslinked_residues_to_domains.pop(uniprot_id)


crosslink_related_domains_dict=defaultdict(dict)
protein_to_domain_to_crosslink_residues={}
for uniprot_id in protein_to_crosslinked_residues_to_domains.keys():
    domain_to_crosslink_residues=defaultdict(dict)
    for residue_id in protein_to_crosslinked_residues_to_domains[uniprot_id].keys():
        for domain_id in protein_to_crosslinked_residues_to_domains[uniprot_id][residue_id]:
            domain_to_crosslink_residues[domain_id][residue_id]=np.zeros(3,dtype=np.float32)
    # 获取位置
    for domain_id in domain_to_crosslink_residues.keys():
        domain_path=os.path.join(origin_domain_dir,uniprot_id+'_'+domain_id+'.pdb')
        u=mda.Universe(domain_path)
        for residue_id in domain_to_crosslink_residues[domain_id].keys():
            position=u.select_atoms('resid '+residue_id+" and name CA").positions[0]
            crosslink_node=uniprot_id+':'+residue_id
            crosslink_related_domains_dict[crosslink_node][domain_id]=position
            domain_to_crosslink_residues[domain_id][residue_id]=position
    protein_to_domain_to_crosslink_residues[uniprot_id]=domain_to_crosslink_residues




# 判断两个mrc密度是否相邻

def get_mrc_data_and_params(mrc_filepath):
    with mrcfile.open(mrc_filepath, mode='r') as mrc:
        data = mrc.data.T.copy()
        origin=np.array(mrc.header.origin.tolist())
        lengths=np.array(mrc.header.cella.tolist())
        grid_shape=np.array([mrc.header.nx,mrc.header.ny,mrc.header.nz])
        voxel=lengths[0]/grid_shape[0]
        return data, origin, lengths, grid_shape, voxel
    
def get_cg_data(data, lengths, grid_shape, voxel, map_level=0, cg_voxel=10):
    # 低于level的数据点赋0
    data[data<map_level]=0
    # 将非0数据点划分进更大的盒子里
    cg_grid_shape=np.ceil(lengths/cg_voxel).astype(int)
    data_cg=np.zeros(cg_grid_shape)
    for i in range(grid_shape[0]):
        for j in range(grid_shape[1]):
            for k in range(grid_shape[2]):
                if data[i,j,k]>0:
                    cg_index=(((0.5+np.array([i,j,k]))*voxel)/cg_voxel).astype(int)
                    data_cg[*cg_index]+=1
    return data_cg, cg_voxel

def is_adjacent(mrc_path1, mrc_path2, map_level_1, map_level_2,cg_voxel=10,min_ratio=0.1, max_dist=60):
    data1, origin1, lengths1, grid_shape1, voxel1 = get_mrc_data_and_params(mrc_path1)
    data2, origin2, lengths2, grid_shape2, voxel2 = get_mrc_data_and_params(mrc_path2)
    cg_data1, cg_voxel1=get_cg_data(data1, lengths1, grid_shape1, voxel1, map_level_1, cg_voxel)
    cg_data2, cg_voxel2=get_cg_data(data2, lengths2, grid_shape2, voxel2, map_level_2, cg_voxel)
    # 所有超过一定ratio的cg格点的id
    cg_ids1=np.array(np.where(cg_data1>min_ratio*(cg_voxel1/voxel1)**3)).T
    cg_ids2=np.array(np.where(cg_data2>min_ratio*(cg_voxel2/voxel2)**3)).T
    # 坐标
    cg_coords1=(0.5+cg_ids1)*cg_voxel1+origin1
    cg_coords2=(0.5+cg_ids2)*cg_voxel2+origin2
    # 计算距离矩阵
    dist_matrix=cdist(cg_coords1,cg_coords2)
    # 最小距离是否小于阈值
    min_dist=np.min(dist_matrix)
    return min_dist<max_dist

# 将密度用编号表示为网络
Density_adj_graph=nx.Graph()
# 遍历两两密度，判断是否相邻
for density_id_1 in range(0,n_densities-1):
    for density_id_2 in range(density_id_1+1,n_densities):
        density_1=density_names[density_id_1]
        density_2=density_names[density_id_2]
        density_1_path=os.path.join(map_dir,density_1+".mrc")
        density_2_path=os.path.join(map_dir,density_2+".mrc")
        if is_adjacent(density_1_path,density_2_path,map_level,map_level):
            Density_adj_graph.add_edge(str(density_id_1),str(density_id_2))


# fit变换
def get_transformation_matrix_list(fit_log_path):
    log_data=np.loadtxt(fit_log_path,dtype=float)
    transform_matrix_list=log_data[:,2:].reshape((-1,3,4))
    return transform_matrix_list

def transform_positions(xyzs,transformation_matrix):
    xyzs=np.array(xyzs).reshape((-1,3))
    rotation_matrix=transformation_matrix[:,:3]
    translation_vector=transformation_matrix[:,-1]
    return xyzs.dot(rotation_matrix.T)+translation_vector


# 交联残基的距离分布
# 目前只考虑臂长和残基侧链长接近DSS的
mu=3      # 3
sigma=0.5   # 0.5

def lognormal_pdf(d,mu,sigma):
    max=np.exp((sigma**2)/2-mu)/(np.sqrt(2*np.pi)*sigma)
    return 1/(max*d*sigma*np.sqrt(2*np.pi))*np.exp(-0.5*((np.log(d)-mu)/sigma)**2)

def get_multiplier(distance,mu,sigma,evidence_strenth):
    return evidence_strenth*lognormal_pdf(distance,mu,sigma)





# 预先读取并计算每个密度中所有交联相关残基的位置，减少后续重复读取文件
# 暂时未考虑对称性
fitted_positions_of_crosslink_related_residues_in_Densities = {Density:{} for Density in Density_adj_graph.nodes}
for density_id in range(len(density_names)):
    for cl, data in crosslink_related_domains_dict.items():
        protein_id, residue_id = cl.split(':')
        domain_ids_and_fit_ids = []
        fitted_positions = []
        for domain_id, position in data.items():
            fit_log_path=os.path.join(fitout_dir,density_names[density_id]+".mrc", "fitlogs" ,f"{protein_id}_{domain_id}.log")
            if os.path.exists(fit_log_path):
                    transformation_matrix_list=get_transformation_matrix_list(fit_log_path)
                    for fit_id, transformation_matrix in enumerate(transformation_matrix_list):
                        # 只考虑先验概率高于阈值的状态
                        state=f"{protein_id}_{domain_id}_{fit_id}"
                        state_id=state_to_id_of_densities[density_id][state]
                        prior_probability=prior_probabilities_of_densities[density_id][state_id-1]
                        if prior_probability >= acceptor_prior_probability_cutoff:
                            fitted_position=transform_positions(position,transformation_matrix)[0]
                            domain_ids_and_fit_ids.append([str(domain_id), fit_id])
                            fitted_positions.append(fitted_position)
        if fitted_positions:
            fitted_positions_of_crosslink_related_residues_in_Densities[str(density_id)][cl]=[domain_ids_and_fit_ids,np.array(fitted_positions)]





# 构建状态类
crosslink_distance_cutoff=50    # DSS的参数
squared_crosslink_distance_cutoff=crosslink_distance_cutoff**2
# 非偶然交联概率/偶然交联概率


# 每个密度给出若干状态，每个状态都存在符合交联约束的残基
# 0表示others
crosslink_compliant_states_of_densities=[{0} for _ in range(len(density_names))]
crosslink_compliant_Density_states_graph=nx.Graph()
# 用于遍历交联和密度对的对应
crosslinks_di={tuple([str(row[0]),str(row[1])])
               for row in np.concatenate([np.array(crosslinked_residues_graph.edges),np.array(crosslinked_residues_graph.edges)[:,[1,0]]],axis=0)}
for Density_1,Density_2 in Density_adj_graph.edges:
    crosslink_related_residues_1=fitted_positions_of_crosslink_related_residues_in_Densities[Density_1].keys()
    crosslink_related_residues_2=fitted_positions_of_crosslink_related_residues_in_Densities[Density_2].keys()
    density_id_1 = int(Density_1)
    density_id_2 = int(Density_2)
    for cl_1,cl_2 in crosslinks_di:
        if cl_1 in crosslink_related_residues_1 and cl_2 in crosslink_related_residues_2:
            protein_id_1, residue_id_1 = cl_1.split(':')
            protein_id_2, residue_id_2 = cl_2.split(':')
            fitted_positions_1=fitted_positions_of_crosslink_related_residues_in_Densities[Density_1][cl_1][1]
            fitted_positions_2=fitted_positions_of_crosslink_related_residues_in_Densities[Density_2][cl_2][1]
            # 计算距离矩阵
            distace_matrix=cdist(fitted_positions_1, fitted_positions_2, metric="sqeuclidean")
            # 查找平方距离不大于阈值的索引
            filtered_index=np.array(np.where(distace_matrix<=squared_crosslink_distance_cutoff)).T
            # 记录每一对符合约束的密度对的状态对
            for i, j in filtered_index:
                domain_id_1, fit_id_1 = fitted_positions_of_crosslink_related_residues_in_Densities[Density_1][cl_1][0][i]
                domain_id_2, fit_id_2 = fitted_positions_of_crosslink_related_residues_in_Densities[Density_2][cl_2][0][j]
                # 判断自交联的可能性，有重复计算
                if protein_id_1==protein_id_2 and domain_id_1==domain_id_2:
                    position_1=crosslink_related_domains_dict[cl_1][domain_id_1]
                    position_2=crosslink_related_domains_dict[cl_2][domain_id_2]
                    squared_distance_self=np.sum((position_1-position_2)**2)
                    if squared_distance_self<=squared_crosslink_distance_cutoff:
                        continue
                # 更新交联相关状态字典
                state_1=f"{protein_id_1}_{domain_id_1}_{fit_id_1}"
                state_id_1=state_to_id_of_densities[density_id_1][state_1]
                prior_probability_1=prior_probabilities_of_densities[density_id_1][state_id_1-1]
                
                state_2=f"{protein_id_2}_{domain_id_2}_{fit_id_2}"
                state_id_2=state_to_id_of_densities[density_id_2][state_2]
                prior_probability_2=prior_probabilities_of_densities[density_id_2][state_id_2-1]

                # 至少有一个能作为供体
                if prior_probability_1<donor_prior_probability_cutoff and prior_probability_2<donor_prior_probability_cutoff:
                    continue
                crosslink_compliant_states_of_densities[density_id_1].add(state_id_1)
                crosslink_compliant_states_of_densities[density_id_2].add(state_id_2)
                # 更新状态图，节点为密度和状态，边表示符合交联约束
                node_1=f"{Density_1}:{state_id_1}"
                node_2=f"{Density_2}:{state_id_2}"
                if not crosslink_compliant_Density_states_graph.has_edge(node_1,node_2):
                    crosslink_compliant_Density_states_graph.add_edge(node_1,node_2,data_list=[])
                crosslink_compliant_Density_states_graph[node_1][node_2]["data_list"].append({"crosslink":frozenset({cl_1,cl_2}),
                                                                                            "crosslinked_residues":{node_1:cl_1,node_2:cl_2},
                                                                                            "distance":np.sqrt(distace_matrix[i,j]),
                                                                                            "multiplier":get_multiplier(np.sqrt(distace_matrix[i,j]),mu,sigma,evidence_strenth)})

# set转化为list
crosslink_compliant_states_of_densities=[list(item) for item in crosslink_compliant_states_of_densities]
# 计算others状态的先验概率
prior_probabilities_of_others_of_densities=[]
for density_id, state_ids in enumerate(crosslink_compliant_states_of_densities):
    p=1
    for state_id in state_ids[1:]:
        p-=prior_probabilities_of_densities[density_id][state_id-1]
    prior_probabilities_of_others_of_densities.append(p)


# 给所有用到的交联编号
sorted_crosslinks=set()
for node_1, node_2 in crosslink_compliant_Density_states_graph.edges:
    data_list=crosslink_compliant_Density_states_graph[node_1][node_2]["data_list"]
    for data in data_list:
        sorted_crosslinks.add(data["crosslink"])
sorted_crosslinks=sorted(sorted_crosslinks)
crosslink_to_index={cl:i for i,cl in enumerate(sorted_crosslinks)}

# 合并不同模块的相同密度
# 合并相同密度间的相同交联
crosslink_compliant_density_states_graph=nx.Graph()
for node_1,node_2 in crosslink_compliant_Density_states_graph.edges:
    Density_1,state_id_1=node_1.split(":")
    Density_2,state_id_2=node_2.split(":")
    density_id_1=int(Density_1.split(".")[0])
    density_id_2=int(Density_2.split(".")[0])
    new_node_1=f"{density_id_1}:{state_id_1}"
    new_node_2=f"{density_id_2}:{state_id_2}"
    if crosslink_compliant_density_states_graph.has_edge(new_node_1,new_node_2):
        crosslink_to_multiplier=crosslink_compliant_density_states_graph[new_node_1][new_node_2]["data_list"]
    else:
        crosslink_to_multiplier=defaultdict(float)
    for data in crosslink_compliant_Density_states_graph[node_1][node_2]["data_list"]:
        crosslink=crosslink_to_index[data["crosslink"]]
        multiplier=data["multiplier"]
        crosslink_to_multiplier[crosslink]+=float(multiplier)
    crosslink_compliant_density_states_graph.add_edge(new_node_1,new_node_2,data_list=crosslink_to_multiplier)
# 将data_list中的字典转化为列表
for node_1,node_2 in crosslink_compliant_density_states_graph.edges:
    crosslink_compliant_density_states_graph[node_1][node_2]["data_list"]=[item for item in crosslink_compliant_density_states_graph[node_1][node_2]["data_list"].items()]


# 保存符合交联约束的残基对
compliant_crosslinks = {} # density_1 : state_1 : [[density_2, state_2, residue_id_1, residue_id_2],... ]
for node_1, node_2 in crosslink_compliant_Density_states_graph.edges:
    Density_id_1, state_id_1 = node_1.split(':')
    density_id_1 = int(Density_id_1)
    state_id_1 = int(state_id_1)
    Density_id_2, state_id_2 = node_2.split(':')
    density_id_2 = int(Density_id_2)
    state_id_2 = int(state_id_2)
    if density_id_1 > density_id_2:
        node_1, node_2 = node_2, node_1
        density_id_1, density_id_2 = density_id_2, density_id_1
        state_id_1, state_id_2 = state_id_2, state_id_1
    density_1 = density_names[density_id_1]
    if compliant_crosslinks.get(density_1) is None:
        compliant_crosslinks[density_1] = {}
    density_2 = density_names[density_id_2]
    state_1 = states_of_densities[density_id_1][state_id_1-1]
    state_2 = states_of_densities[density_id_2][state_id_2-1]
    for data in crosslink_compliant_Density_states_graph[node_1][node_2]["data_list"]:
        residue_id_1 = int(data['crosslinked_residues'][node_1].split(":")[-1])
        residue_id_2 = int(data['crosslinked_residues'][node_2].split(":")[-1])
        if compliant_crosslinks[density_1].get(state_1) is None:
            compliant_crosslinks[density_1][state_1] = []
        compliant_crosslinks[density_1][state_1].append([density_2, state_2 ,residue_id_1, residue_id_2])
np.save(os.path.join(project_dir,"compliant_crosslinks.npy"), compliant_crosslinks)




possible_crosslink_related_densities_dict = defaultdict(set)
possible_crosslink_related_data_dict = defaultdict(list) 
compliant_states_of_densities_of_single_crosslink={}
compliant_crosslink_between_densities_graph=nx.Graph()
for node_1, node_2 in crosslink_compliant_density_states_graph.edges:
    density_id_1, state_id_1 = [int(item) for item in node_1.split(":")]
    density_id_2, state_id_2 = [int(item) for item in node_2.split(":")]
    for cl, multiplier in crosslink_compliant_density_states_graph[node_1][node_2]["data_list"]:
        if cl not in compliant_states_of_densities_of_single_crosslink.keys():
            compliant_states_of_densities_of_single_crosslink[cl]={}
        possible_crosslink_related_densities_dict[cl].add(density_id_1)
        possible_crosslink_related_densities_dict[cl].add(density_id_2)
        if density_id_1 not in compliant_states_of_densities_of_single_crosslink[cl].keys():
            compliant_states_of_densities_of_single_crosslink[cl][density_id_1]=[]
        compliant_states_of_densities_of_single_crosslink[cl][density_id_1].append(state_id_1)
        if density_id_2 not in compliant_states_of_densities_of_single_crosslink[cl].keys():
            compliant_states_of_densities_of_single_crosslink[cl][density_id_2]=[]
        compliant_states_of_densities_of_single_crosslink[cl][density_id_2].append(state_id_2)
        possible_crosslink_related_data_dict[cl].append((node_1, node_2, multiplier))
        if not compliant_crosslink_between_densities_graph.has_edge(density_id_1, density_id_2):
            compliant_crosslink_between_densities_graph.add_edge(density_id_1, density_id_2, compliant_crosslinks={})
        compliant_crosslink_between_densities_graph[density_id_1][density_id_2]["compliant_crosslinks"][cl]=multiplier

density_groups = []
for cl, densities in possible_crosslink_related_densities_dict.items():
    if len(densities)==2:
        continue
    for existing_group in density_groups:
        if densities.intersection(existing_group):
            existing_group.update(densities)
            break
    else:
        density_groups.append(densities)
# 每一个未分组的涉及交联的density，都单独作为一个group
for density in set(compliant_crosslink_between_densities_graph.nodes)-set().union(*density_groups):
    density_groups.append({density})
# 重写为列表并排序
density_groups = sorted([sorted(group) for group in density_groups], key=lambda x: min(x))
grouped_densities = set(density_id for group in density_groups for density_id in group)
# 逆索引
density_id_to_group_id = {density:i for i, group in enumerate(density_groups) for density in group}
# 组内索引
density_id_in_group_dict = {density:i for group in density_groups for i, density in enumerate(group)}




# 自定义网络结构
# 此处记录的density id都是组内id
# 组内数据
intra_group_multipliers_data=[{} for _ in density_groups] # group_id -> {density_id_in_group_1 : {state_id_1 : [(density_id_in_group_2, state_id_2,data_list)]}}
# 组间数据
# 组间无需保存交联编号
inter_group_multipliers_data=[] # (group_id_1, group_id_2, density_id_in_group_1, state_id_1, density_id_in_group_2, state_id_2, multiplier) 

for node_1, node_2 in crosslink_compliant_density_states_graph.edges:
    density_id_1, state_id_1 = [int(item) for item in node_1.split(':')]
    density_id_2, state_id_2 = [int(item) for item in node_2.split(':')]
    group_id_1 = density_id_to_group_id[density_id_1]
    group_id_2 = density_id_to_group_id[density_id_2]
    data_list=crosslink_compliant_density_states_graph[node_1][node_2]["data_list"]
    if group_id_1 == group_id_2:
        # 组内数据
        # 确保 density_id_1 < density_id_2
        if density_id_1 > density_id_2:
            density_id_1, density_id_2 = density_id_2, density_id_1
            state_id_1, state_id_2 = state_id_2, state_id_1
        # 组内id
        density_id_in_group_1 = density_id_in_group_dict[density_id_1]
        density_id_in_group_2 = density_id_in_group_dict[density_id_2]
        # 1为node，2为neighbor
        if density_id_in_group_1 not in intra_group_multipliers_data[group_id_1].keys():
            intra_group_multipliers_data[group_id_1][density_id_in_group_1]={}
        if state_id_1 not in intra_group_multipliers_data[group_id_1][density_id_in_group_1].keys():
            intra_group_multipliers_data[group_id_1][density_id_in_group_1][state_id_1]=[]
        intra_group_multipliers_data[group_id_1][density_id_in_group_1][state_id_1].append((density_id_in_group_2, state_id_2,data_list))
    else:
        # 组间数据
        multiplier = np.prod([row[1]+1 for row in data_list])
        # 确保 group_id_1 < group_id_2
        if group_id_1 > group_id_2:
            group_id_1, group_id_2 = group_id_2, group_id_1
            density_id_1, density_id_2 = density_id_2, density_id_1
            state_id_1, state_id_2 = state_id_2, state_id_1
        # 组内id
        density_id_in_group_1 = density_id_in_group_dict[density_id_1]
        density_id_in_group_2 = density_id_in_group_dict[density_id_2]

        inter_group_multipliers_data.append((group_id_1, group_id_2, density_id_in_group_1, state_id_1, density_id_in_group_2, state_id_2, multiplier))

# 将group分为不同super_group，不同super_group之间无交联
super_groups = []
for group_id_1, group_id_2, density_id_in_group_1, state_id_1, density_id_in_group_2, state_id_2, multiplier in inter_group_multipliers_data:
    for sg in super_groups:
        if group_id_1 in sg or group_id_2 in sg:
            sg.add(group_id_1)
            sg.add(group_id_2)
            break
    else:
        super_groups.append(set([group_id_1, group_id_2]))
# 所有剩余group都单独构成super_group
for group_id in range(len(density_groups)):
    if group_id not in [g for sg in super_groups for g in sg]:
        super_groups.append(set([group_id]))
# 集合转为列表并按最小group_id排序
super_groups = sorted([sorted(group_set) for group_set in super_groups], key=lambda x: min(x))

group_id_to_super_group_id = {group_id: super_group_id for super_group_id, group_set in enumerate(super_groups) for group_id in group_set}
group_id_in_super_group_dict = {group_id: group_id_in_super_group for group in super_groups for group_id_in_super_group, group_id in enumerate(group)}

inter_group_multipliers_data_of_super_groups = [[] for _ in super_groups]
for group_id_1, group_id_2, density_id_in_group_1, state_id_1, density_id_in_group_2, state_id_2, multiplier in inter_group_multipliers_data:
    for super_group_id, super_group in enumerate(super_groups):
        if group_id_1 in super_group:
            inter_group_multipliers_data_of_super_groups[super_group_id].append((group_id_1, group_id_2, density_id_in_group_1, state_id_1, density_id_in_group_2, state_id_2, multiplier))
            break



# 组间状态类
inter_group_crosslink_related_states_of_densities_of_groups_of_super_groups = [[{} for group_id_in_super_group in range(len(super_groups[super_group_id]))] for super_group_id in range(len(super_groups))]
for super_group_id, super_group in enumerate(super_groups):
    for group_id_1, group_id_2, density_id_in_group_1, state_id_1, density_id_in_group_2, state_id_2, multiplier in inter_group_multipliers_data_of_super_groups[super_group_id]:
        group_id_in_super_group_1 = group_id_in_super_group_dict[group_id_1]
        group_id_in_super_group_2 = group_id_in_super_group_dict[group_id_2]
        if inter_group_crosslink_related_states_of_densities_of_groups_of_super_groups[super_group_id][group_id_in_super_group_1].get(density_id_in_group_1) is None:
            inter_group_crosslink_related_states_of_densities_of_groups_of_super_groups[super_group_id][group_id_in_super_group_1][density_id_in_group_1]=set()
        inter_group_crosslink_related_states_of_densities_of_groups_of_super_groups[super_group_id][group_id_in_super_group_1][density_id_in_group_1].add(state_id_1)
        if inter_group_crosslink_related_states_of_densities_of_groups_of_super_groups[super_group_id][group_id_in_super_group_2].get(density_id_in_group_2) is None:
            inter_group_crosslink_related_states_of_densities_of_groups_of_super_groups[super_group_id][group_id_in_super_group_2][density_id_in_group_2]=set()
        inter_group_crosslink_related_states_of_densities_of_groups_of_super_groups[super_group_id][group_id_in_super_group_2][density_id_in_group_2].add(state_id_2)
# 重构
inter_group_crosslink_related_density_ids_in_group_of_group_of_super_group = [[[] for group_id_in_super_group in range(len(super_groups[super_group_id]))] for super_group_id in range(len(super_groups))]
inter_group_crosslink_related_partial_states_of_groups_of_super_groups = [[[] for group_id_in_super_group in range(len(super_groups[super_group_id]))] for super_group_id in range(len(super_groups))]
density_id_in_group_to_id_in_partial_group_of_group_of_super_group = [[{} for group_id_in_super_group in range(len(super_groups[super_group_id]))] for super_group_id in range(len(super_groups))]
for super_group_id, super_group in enumerate(super_groups):
    for group_id_in_super_group in range(len(super_group)):
        related_density_ids_in_group = sorted(inter_group_crosslink_related_states_of_densities_of_groups_of_super_groups[super_group_id][group_id_in_super_group].keys())
        related_partial_group_states = list(itertools.product(*[sorted(inter_group_crosslink_related_states_of_densities_of_groups_of_super_groups[super_group_id][group_id_in_super_group][density_id_in_group].union([-1]))
                                                                for density_id_in_group in related_density_ids_in_group]))
        inter_group_crosslink_related_density_ids_in_group_of_group_of_super_group[super_group_id][group_id_in_super_group]=related_density_ids_in_group
        density_id_in_group_to_id_in_partial_group_of_group_of_super_group[super_group_id][group_id_in_super_group]=dict(zip(related_density_ids_in_group, range(len(related_partial_group_states))))
        inter_group_crosslink_related_partial_states_of_groups_of_super_groups[super_group_id][group_id_in_super_group]=related_partial_group_states

        
# 按每个密度的状态进行状态分类
def get_intra_group_state_multipliers(group_state_class, group_id):
    crosslink_multipliers=[]
    for density_id_in_group_1, state_id_1 in enumerate(group_state_class):
        if intra_group_multipliers_data[group_id].get(density_id_in_group_1) is None:
            continue
        edges=intra_group_multipliers_data[group_id][density_id_in_group_1].get(state_id_1)
        if edges is None:
            continue
        if state_id_1 != 0:
            for density_id_in_group_2, state_id_2, data_list in edges:
                if state_id_2 == group_state_class[density_id_in_group_2]:
                    crosslink_multipliers.extend(data_list)
    return crosslink_multipliers

def get_multiplier_of_intro_group_state_class(group_state_class,group_id):
    crosslink_multipliers=np.array(get_intra_group_state_multipliers(group_state_class,group_id))
    if len(crosslink_multipliers)==0:
        return 1
    # 分组求和
    crosslink_ids=crosslink_multipliers[:,0]
    multipliers=crosslink_multipliers[:,1]
    unique_crosslink_ids, inverse_indices = np.unique(crosslink_ids, return_inverse=True, return_counts=False)
    group_sums = np.zeros_like(unique_crosslink_ids, dtype=multipliers.dtype)
    np.add.at(group_sums, inverse_indices, multipliers)
    group_sums = 1 + group_sums
    multiplier = np.prod(group_sums)
    return multiplier

def get_probabilities_of_group_state_class(group_state_class,group_id):
    p=1
    for density_id_in_group, state_id in enumerate(group_state_class):
        density_id = density_groups[group_id][density_id_in_group]
        if state_id != 0:
            p*=prior_probabilities_of_densities[density_id][state_id-1]
        else:
            p*=prior_probabilities_of_others_of_densities[density_id]
    return p

def get_group_state_class_info(group_state_class,group_id,group_state_class_digital):
    return group_state_class_digital,get_multiplier_of_intro_group_state_class(group_state_class,group_id),get_probabilities_of_group_state_class(group_state_class,group_id)


def get_multiplier_of_inter_group_state_class(inter_group_state_class,super_group_id):
    inter_group_multiplier = 1
    for group_id_1, group_id_2, density_id_in_group_1, state_id_1, density_id_in_group_2, state_id_2, multiplier in inter_group_multipliers_data_of_super_groups[super_group_id]:
        group_id_in_super_group_1 = group_id_in_super_group_dict[group_id_1]
        group_id_in_super_group_2 = group_id_in_super_group_dict[group_id_2]
        density_id_in_partial_group_1=density_id_in_group_to_id_in_partial_group_of_group_of_super_group[super_group_id][group_id_in_super_group_1][density_id_in_group_1]
        density_id_in_partial_group_2=density_id_in_group_to_id_in_partial_group_of_group_of_super_group[super_group_id][group_id_in_super_group_2][density_id_in_group_2]
        if inter_group_state_class[group_id_in_super_group_1][density_id_in_partial_group_1] == state_id_1 \
           and \
           inter_group_state_class[group_id_in_super_group_2][density_id_in_partial_group_2] == state_id_2:
            inter_group_multiplier *= multiplier
    return inter_group_multiplier


def get_inter_group_state_class_info(inter_group_state_class,super_group_id,inter_group_state_class_digital):
    return inter_group_state_class_digital,get_multiplier_of_inter_group_state_class(inter_group_state_class,super_group_id)





# 非标准进制转换

# 分组
n_digits_list_of_groups=[[len(crosslink_compliant_states_of_densities[density_id]) for density_id in group] for group in density_groups]
digital_weights_of_groups=[np.cumprod((n_digits_list_of_groups[group_id]+[1])[::-1]).tolist()[::-1] for group_id in range(len(density_groups))]

def group_state_class_code_to_digital(group_state_class_code,group_id):
    total = 0
    for i, digit in enumerate(group_state_class_code):
        total += digit * digital_weights_of_groups[group_id][i+1]
    return total

# 从数字计算state_class_code的第i位
def group_state_class_digital_to_code(group_state_class_digital,group_id,density_id_in_group):
    return (group_state_class_digital%digital_weights_of_groups[group_id][density_id_in_group])//digital_weights_of_groups[group_id][density_id_in_group+1]

# state_class_code转化为state_class
def group_state_class_code_to_state_class(group_state_class_code,group_id):
    return [crosslink_compliant_states_of_densities[density_id][group_state_class_code[density_id_in_group]] for density_id_in_group, density_id in enumerate(density_groups[group_id])]



n_digits_list_inter_groups = [[len(inter_group_crosslink_related_partial_states_of_groups_of_super_groups[super_group_id][group_id_in_super_group])
                                 for group_id_in_super_group in range(len(super_groups[super_group_id]))]
                                for super_group_id in range(len(super_groups))]
digital_weights_inter_groups = [np.cumprod((n_digits_list_inter_groups[super_group_id]+[1])[::-1]).tolist()[::-1] for super_group_id in range(len(super_groups))]

# 从数字计算inter_group_state_class_code的第i位
def inter_group_state_class_digital_to_code(inter_group_state_class_digital,super_group_id,group_id_in_super_group):
    return (inter_group_state_class_digital%digital_weights_inter_groups[super_group_id][group_id_in_super_group])//digital_weights_inter_groups[super_group_id][group_id_in_super_group+1]


# inter_group_state_class_code转化为partial_state_class
def inter_group_state_class_code_to_state_class(inter_group_state_class_code,super_group_id):
    return [inter_group_crosslink_related_partial_states_of_groups_of_super_groups[super_group_id][group_id_in_super_group][c] for group_id_in_super_group, c in enumerate(inter_group_state_class_code)]

# inter_group_state_class_digital转化为partial_state_class
def inter_group_state_class_digital_to_state_class(inter_group_state_class_digital,super_group_id):
    return inter_group_state_class_code_to_state_class([inter_group_state_class_digital_to_code(inter_group_state_class_digital,super_group_id,group_id_in_super_group)
                                                        for group_id_in_super_group in range(len(super_groups[super_group_id]))],super_group_id)




# 分组状态类遍历

# 按每个密度的状态进行状态分类
# 组内状态类
group_state_classes_info=[[] for _ in density_groups]

for group_id, group in enumerate(density_groups):
    number_of_group_classes=np.prod(n_digits_list_of_groups[group_id])
    # itertools.product计算若干列表的直积
    group_state_class_digital=0
    crosslink_compliant_states_of_densities_of_group=[crosslink_compliant_states_of_densities[density_id] for density_id in group]
    for group_state_class in tqdm(itertools.product(*crosslink_compliant_states_of_densities_of_group),total=number_of_group_classes):
        group_state_classes_info[group_id].append(get_group_state_class_info(group_state_class,group_id,group_state_class_digital))
        group_state_class_digital+=1


# 构建各super_group内的组间状态类
inter_group_state_classes_info = [[] for _ in super_groups]

for super_group_id, super_group in enumerate(super_groups):
    inter_group_state_class_digital=0
    for inter_group_state_class in itertools.product(*[inter_group_crosslink_related_partial_states_of_groups_of_super_groups[super_group_id][group_id_in_super_group] for group_id_in_super_group in range(len(super_group))]):
        inter_group_state_classes_info[super_group_id].append(get_inter_group_state_class_info(inter_group_state_class,super_group_id,inter_group_state_class_digital))
        inter_group_state_class_digital+=1




def get_intra_group_average_multiplier(group_id,super_group_id, restrictions = []):
    average_multiplier=0
    group_id_in_super_group=group_id_in_super_group_dict[group_id]
    
    for group_state_class_digital,group_multiplier,group_state_class_probability in group_state_classes_info[group_id]:
        # 验证是否满足所有约束
        flag=True
        for density_id_in_group_in_restriction, state_id_in_restriction in restrictions:
            density_id=density_groups[group_id][density_id_in_group_in_restriction]
            state_code=group_state_class_digital_to_code(group_state_class_digital,group_id,density_id_in_group_in_restriction)
            state_id=crosslink_compliant_states_of_densities[density_id][state_code]
            # 验证约束
            if state_id_in_restriction==-1:
                if state_id in inter_group_crosslink_related_states_of_densities_of_groups_of_super_groups[super_group_id][group_id_in_super_group][density_id_in_group_in_restriction]:
                    flag=False
                    break
            elif state_id_in_restriction!=state_id:
                flag=False
                break
        if not flag:
            continue
        # 计入平均
        average_multiplier+=group_multiplier*group_state_class_probability
    return average_multiplier

def inter_group_state_class_to_restrictions(inter_group_state_class,super_group_id):
    restrictions_of_groups=[[] for i in range(len(super_groups[super_group_id]))]
    for group_id_in_super_group in range(len(super_groups[super_group_id])):
        for density_id_in_partial_group, state_id in enumerate(inter_group_state_class[group_id_in_super_group]):
            density_id_in_group=inter_group_crosslink_related_density_ids_in_group_of_group_of_super_group[super_group_id][group_id_in_super_group][density_id_in_partial_group]
            restrictions_of_groups[group_id_in_super_group].append((density_id_in_group,state_id))
    return restrictions_of_groups

def get_average_multiplier(density_id_in_group,group_id_in_super_group,super_group_id,state_id):
    group_id=super_groups[super_group_id][group_id_in_super_group]
    density_id=density_groups[group_id][density_id_in_group]
    # 只在super_group内求平均
    average_multiplier=0
    # 设置分母，以便单独拿出“others”状态，减少重复计算
    if state_id==0:
        p_denominator=prior_probabilities_of_others_of_densities[density_id]
    else:
        p_denominator=prior_probabilities_of_densities[density_id][state_id-1]
    # 遍历所有组间状态类，获取约束
    for inter_group_state_class_digital, inter_group_multiplier in inter_group_state_classes_info[super_group_id]:
        average_multiplier_of_super_group=inter_group_multiplier
        inter_group_state_class=inter_group_state_class_digital_to_state_class(inter_group_state_class_digital,super_group_id)
        # 验证组间状态是否与给定状态一致
        density_id_in_partial_group=density_id_in_group_to_id_in_partial_group_of_group_of_super_group[super_group_id][group_id_in_super_group].get(density_id_in_group)
        if density_id_in_partial_group is None:
            restrictions_of_groups=inter_group_state_class_to_restrictions(inter_group_state_class,super_group_id)
            # 添加给定状态的约束
            restrictions_of_groups[group_id_in_super_group].append((density_id_in_group,state_id))
        else:
            state_id_in_inter_group_state_class=inter_group_state_class[group_id_in_super_group][density_id_in_partial_group]
            if state_id_in_inter_group_state_class==-1:
                if state_id not in inter_group_crosslink_related_states_of_densities_of_groups_of_super_groups[super_group_id][group_id_in_super_group][density_id_in_group]:
                    restrictions_of_groups=inter_group_state_class_to_restrictions(inter_group_state_class,super_group_id)
                    # 添加给定状态的约束
                    restrictions_of_groups[group_id_in_super_group].append((density_id_in_group,state_id))
                else:
                    continue
            else:
                if state_id==state_id_in_inter_group_state_class:
                    restrictions_of_groups=inter_group_state_class_to_restrictions(inter_group_state_class,super_group_id)
                else:
                    continue
        # 遍历组
        for group_id_in_super_group_i, group_id_i in enumerate(super_groups[super_group_id]):
            intra_group_average_multiplier_of_compliant_state_classes=get_intra_group_average_multiplier(group_id_i, super_group_id ,restrictions_of_groups[group_id_in_super_group_i])
            average_multiplier_of_super_group*=intra_group_average_multiplier_of_compliant_state_classes
        average_multiplier+=average_multiplier_of_super_group
    return average_multiplier/p_denominator


def normalization_and_sorting(state_ids,probabilities):
    order_list=np.argsort(probabilities)[::-1]
    state_ids=state_ids[order_list]
    probabilities=probabilities[order_list]/np.sum(probabilities)
    return state_ids,probabilities

def save_posterior_results_of_density(save_path,density_id,sorted_posterior_state_ids,sorted_posterior_probabilities):
    save_data=np.array([[states_of_densities[density_id][state_id-1] for state_id in sorted_posterior_state_ids],sorted_posterior_probabilities],dtype=str).T
    np.savetxt(save_path,save_data,fmt='%s')




# 调整概率
sorted_posterior_state_ids_of_densities=[0 for density_id in range(n_densities)]
sorted_posterior_probabilities_of_densities=[0 for density_id in range(n_densities)]

# 处理grouped_densities
print("Processing grouped_densities :")
for super_group_id in range(len(super_groups)):
    print("Super group", super_group_id)
    for group_id_in_super_group, group_id in enumerate(super_groups[super_group_id]):
        print("\t","Group", group_id)
        for density_id_in_group, density_id in enumerate(density_groups[group_id]):
            print("\t\t","Density", density_id,density_names[density_id])
            # 预处理
            posterior_state_ids=np.arange(1, len(states_of_densities[density_id])+1)
            posterior_probabilities=np.array(prior_probabilities_of_densities[density_id])

            # average multiplier for others
            average_multiplier_others=get_average_multiplier(density_id_in_group,group_id_in_super_group,super_group_id,state_id=0)
            
            # 遍历该density的所有状态
            for state_id in range(1, len(states_of_densities[density_id])+1):
                if state_id in crosslink_compliant_states_of_densities[density_id]:
                    average_multiplier=get_average_multiplier(density_id_in_group,group_id_in_super_group,super_group_id,state_id)
                else:
                    average_multiplier=average_multiplier_others
                posterior_probabilities[state_id - 1] *= average_multiplier
            # 归一化与排序
            sorted_posterior_state_ids, sorted_posterior_probabilities = normalization_and_sorting(posterior_state_ids, posterior_probabilities)
            sorted_posterior_state_ids_of_densities[density_id]=sorted_posterior_state_ids
            sorted_posterior_probabilities_of_densities[density_id]=sorted_posterior_probabilities
            # 保存结果
            save_path=os.path.join(fitout_dir,density_names[density_id]+".mrc","posterior_probabilities.txt")
            save_posterior_results_of_density(save_path,density_id,sorted_posterior_state_ids,sorted_posterior_probabilities)

print()

# 处理ungrouped_densities
print("Processing ungrouped_densities :")
for density_id in set(range(n_densities))-grouped_densities:
    print("Density", density_id,density_names[density_id])
    # 预处理
    posterior_state_ids=np.arange(1, len(states_of_densities[density_id])+1)
    posterior_probabilities=np.array(prior_probabilities_of_densities[density_id])
    # 直接归一化并写入结果
    sorted_posterior_state_ids, sorted_posterior_probabilities = normalization_and_sorting(posterior_state_ids, posterior_probabilities)
    sorted_posterior_state_ids_of_densities[density_id]=sorted_posterior_state_ids
    sorted_posterior_probabilities_of_densities[density_id]=sorted_posterior_probabilities
    save_path=os.path.join(fitout_dir,density_names[density_id]+".mrc","posterior_probabilities.txt")
    save_posterior_results_of_density(save_path,density_id,sorted_posterior_state_ids,sorted_posterior_probabilities)