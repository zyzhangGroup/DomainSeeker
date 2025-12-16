# 全局异常处理
import domainseeker_errorLog
#----------------------------------------------------------------------------------------------------

import numpy as np
import networkx as nx
import MDAnalysis as mda
import os,sys
import multiprocessing
from tqdm import tqdm
import warnings

pdb_dir=sys.argv[2]
pae_dir=sys.argv[3]
output_dir=sys.argv[4]
n_processors=int(sys.argv[5])

#parameters
if len(sys.argv)>6:
    plddt_cutoff=float(sys.argv[6])
    pae_cutoff=float(sys.argv[7])
    clique_cutoff=int(sys.argv[8])
    min_dege_ratio_between_cliques=float(sys.argv[9])
    min_common_nodes_ratio_between_cliques=float(sys.argv[10])
    minimum_domain_length=int(sys.argv[11])
    maximum_domain_length=int(sys.argv[12])
else:
    plddt_cutoff=70 
    pae_cutoff=5
    clique_cutoff=4
    min_dege_ratio_between_cliques=0.6
    min_common_nodes_ratio_between_cliques=0.5
    minimum_domain_length=40
    maximum_domain_length=1000


def parse_pae_file(pae_json_file):
    import json
    with open(pae_json_file, 'rt') as f:
        data=json.load(f)
        if type(data)==list:
            data=data[0]
    key=list({'pae','predicted_aligned_error'}&set(data.keys()))[0]
    pae_matrix=np.array(data[key],dtype=float)
    np.fill_diagonal(pae_matrix,0)
    return pae_matrix

def get_range_list(cluster):
    cluster=sorted(cluster)
    ranges=[]
    probe=cluster[0]
    range_=[probe]
    for number in cluster[1:]:
        if number==(probe+1):
            probe=number
        else:
            range_.append(probe)
            ranges.append(range_)
            probe=number
            range_=[probe]
    range_.append(probe)
    ranges.append(range_)
    return ranges


def get_select_string(cluster):
    range_list=get_range_list(cluster)
    select_string=''
    for range_ in range_list:
        select_string+=f"{range_[0]}-{range_[1]},"
    select_string=select_string.rstrip(',')
    return select_string 

def get_common_neighbors(G,clique):
    all_nodes=list(G.nodes)
    rest_nodes=list(set(all_nodes)-set(clique))
    tmp_neighbors=rest_nodes.copy()
    for i in rest_nodes:
        for j in clique:
            if not G.has_edge(i,j):
                tmp_neighbors.remove(i)
                break
    common_neighbors=tmp_neighbors
    return common_neighbors

def get_max_degree_node(G):
    degree_dict=dict(G.degree)
    max_degree_node=max(degree_dict,key=lambda x:degree_dict[x])
    return max_degree_node
        
def max_clique(G,clique=None):
    if clique is None:
        clique=[]
    common_neighbors=get_common_neighbors(G,clique)
    while len(common_neighbors)>0:
        g=G.subgraph(common_neighbors)
        max_degree_node=get_max_degree_node(g)
        clique.append(max_degree_node)
        common_neighbors=list(nx.neighbors(g,max_degree_node))
    clique=sorted(clique)
    return clique


def get_residue_graph(u,pae_matrix,plddt_cutoff,pae_cutoff):
    #get the adjacent matrix from the pae matrix
    averaged_pae_matrix=(pae_matrix+pae_matrix.T)/2
    adjacent_matrix=np.where((averaged_pae_matrix<pae_cutoff)*(averaged_pae_matrix>0),1,0)
    #constract the graph
    residue_graph=nx.from_numpy_array(adjacent_matrix)
    residue_graph=nx.relabel_nodes(residue_graph,{i:i+1 for i in range(len(adjacent_matrix))})
    # remove residues with low plddts
    residue_graph.remove_nodes_from(u.select_atoms(f'tempfactor 0:{plddt_cutoff}').residues.resids)
    return residue_graph


def cliques_covering_residue_graph(graph):
    # assert the nodes of G are integers
    queue=set(graph.nodes)
    cliques=[]
    while len(queue)>0:
        g=graph.subgraph(queue)
        clique=max_clique(g)
        clique=tuple(max_clique(graph,clique))
        cliques.append(clique)
        queue=queue-set(clique)
    return cliques


def get_clique_graph(residue_graph,cliques,clique_cutoff,min_dege_ratio_between_cliques,min_common_nodes_ratio_between_cliques):
    selected_cliques=[clique for clique in cliques if len(clique)>clique_cutoff]
    clique_graph=nx.Graph()
    clique_graph.add_nodes_from(selected_cliques)
    for i in range(len(selected_cliques)-1):
        c1=selected_cliques[i]
        for j in range(i+1,len(selected_cliques)):
            c2=selected_cliques[j]
            edge_ratio_between_cliques=len(list(nx.edge_boundary(residue_graph,c1,c2)))/(len(c1)*len(c2))
            common_nodes_ratio_between_cliques=len(set(c1).intersection(set(c2)))/min(len(c1),len(c2))
            if edge_ratio_between_cliques>=min_dege_ratio_between_cliques or common_nodes_ratio_between_cliques>=min_common_nodes_ratio_between_cliques:
                clique_graph.add_edge(c1,c2)
    return clique_graph


def clusters_from_clique_graph_component(clique_graph):
    components=list(nx.connected_components(clique_graph))
    components=[list(component) for component in components]
    clusters=[]
    for component in components:
        s=set(component[0])
        for clique in component[1:]:
            s=s.union(set(clique))
        cluster=tuple(sorted(s))
        clusters.append(cluster)
    clusters=sorted(clusters,key=lambda cluster:-len(cluster))
    return clusters


def get_clusters(u,pae_matrix,plddt_cutoff,pae_cutoff,min_dege_ratio_between_cliques,min_common_nodes_ratio_between_cliques,minimum_domain_length, maximum_domain_length):
    # constract the residue network
    residue_graph=get_residue_graph(u,pae_matrix,plddt_cutoff,pae_cutoff)
    # get cliques that cover the residue network
    cliques=cliques_covering_residue_graph(residue_graph)
    # constract the clique network
    clique_graph=get_clique_graph(residue_graph,cliques,clique_cutoff,min_dege_ratio_between_cliques,min_common_nodes_ratio_between_cliques)
    # get clusters
    clusters=clusters_from_clique_graph_component(clique_graph)
    selected_clusters=[cluster for cluster in clusters if minimum_domain_length<=len(cluster)<=maximum_domain_length]
    return selected_clusters

def save_pdbs(u,clusters,uniprot_id,output_dir):
    if len(clusters)>0:
        info_file=open(os.path.join(output_dir,f"{uniprot_id}.domains"),'w')
        for i in range(len(clusters)):
            cluster=clusters[i]
            sel=u.residues[np.array(cluster)-1]
            output_domain=os.path.join(output_dir,f"{uniprot_id}_D{i}.pdb")
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=UserWarning)
                # 抑制写入pdb时的冗余警告
                sel.atoms.write(output_domain)
            select_string=get_select_string(cluster)
            info_file.write(f"D{i}\t{select_string}\n")
        info_file.close()
    return 0

def run(pdb_file_name):
    uniprot_id=pdb_file_name[:-4]
    pae_path=os.path.join(pae_dir,f"{uniprot_id}.json")
    if os.path.exists(pae_path):
        # get the model
        pdb_path=os.path.join(pdb_dir,f"{uniprot_id}.pdb")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            # 抑制读取pdb时的冗余警告
            u=mda.Universe(pdb_path)
        # get pae matrix
        pae_matrix=parse_pae_file(pae_path)
        # clustering
        clusters=get_clusters(u,pae_matrix,plddt_cutoff,pae_cutoff,min_dege_ratio_between_cliques,min_common_nodes_ratio_between_cliques,minimum_domain_length,maximum_domain_length)
        # save pdbs
        save_pdbs(u,clusters,uniprot_id,output_dir)


# run in parallel
if __name__ == '__main__':
    os.makedirs(output_dir,exist_ok=True)
    file_list=[file_name for file_name in os.listdir(pdb_dir) if file_name.endswith(".pdb")]

    with multiprocessing.Pool(n_processors) as pool:
        # 强制迭代tqdm对象，更新进度条
        for result in tqdm(pool.imap_unordered(run, file_list),total=len(file_list),desc="Parsing",file=sys.stdout):
            # 处理结果
            pass
    
    print("Done domain parsing.")