import os,sys
import subprocess,multiprocessing

def fit(domain_filename):
    domain_path=os.path.join(domain_dir,domain_filename)
    subprocess.call(f'ChimeraX --nogui --offscreen --script \"{script_dir}/fit_in_chimerax.py {domain_path} {map_path} {output_subdir} {map_level} {resolution} {n_search}" --exit',shell=True)

if __name__=="__main__":
    domain_dir=os.path.realpath(sys.argv[1])
    map_dir=os.path.realpath(sys.argv[2])
    output_dir=os.path.realpath(sys.argv[3])
    map_level=float(sys.argv[4])
    resolution=float(sys.argv[5])
    n_search=int(sys.argv[6])
    n_thread=int(sys.argv[7])
    
    script_dir=os.path.dirname(os.path.realpath(__file__))

    map_list=os.listdir(map_dir)
    domain_list=os.listdir(domain_dir)

    os.makedirs(output_dir,exist_ok=True)

    for map_filename in map_list:
        output_subdir=os.path.join(output_dir,map_filename)
        os.makedirs(output_subdir,exist_ok=True)
        map_path=os.path.join(map_dir,map_filename)
        
        config_log_path=os.path.join(output_subdir,"fit_config.txt")
        config_log=open(config_log_path,'w')
        config_log.write(f"{domain_dir=}\n"+
                         f"{map_dir=}\n"+
                         f"{map_level=}\n"+
                         f"{resolution=:.2f}\n"+
                         f"{n_search=}\n"+
                         f"{n_thread=}")
        config_log.close()
        
        pool=multiprocessing.Pool(n_thread)
        pool.map(fit,domain_list)
        pool.close()
