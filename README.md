# DomainSeeker
DomainSeeker is a tool designed for de novo identification of protein domains in cryo-ET maps of protein complexes. DomainSeeker integrates AlphaFold2-predicted models to partition proteins into domains using predicted aligned errors and a clique-based graph cluttering algorithm, and subsequently fits these domains into segmented density maps. Finally, it uses a two-dimensional scoring metric to identify the most probable domain of each density map.

## Installation

### Install dependencies
UCSF ChimeraX 1.5 or higher [(https://www.cgl.ucsf.edu/chimerax/)](https://www.cgl.ucsf.edu/chimerax/)  
	Note: **Add ChimeraX to your system's PATH environment variable**.

#### Install Python dependencies
We recommand to use Python 3.9.
```
pip install mdanalysis
pip install networkx
pip install wget
```

## Workflow and example
We will use an example to demonstrate the workflow.  
We assume that your current working directory is the "example" folder for the exmaple commands.

### 1.Segment the target densities from the original density map.
We offered a target density in the "Example/em_data" folder. It corresponds to the Pacrg protein in the mouse sperm microtubule doublet.  

<img src="/Example/figures/density_Pacrg.jpg" width="400px">

### 2.Fetch pdb and pae files from AFDB
You need to provide a text file containing the Uniprot IDs of all candidate proteins, with each ID on a separate line. We offered a list of 400 proteins in the "protein_list.txt" file, in which Q9DAK2 is the one corresponds to the Pacrg target density.  

Usage:
```
fetch_pdb_pae.py protein_list output_dir
```
This script will download the proteins listed in the "protein_list" file, saving the pdb and pae files separately in the "output_dir/pdb_files" and "output_dir/pae_files" directories. Proteins that fail to download will be logged in the "missing_pdb.log" and "missing_pae.log" files.  

Example commands:
```
python ../fetch_pdb_pae.py protein_list.txt .
```
### 3.Domain parsing
Usage:
```
parse_with_pae.py pdb_dir pae_dir output_dir [plddt_cutoff=70]
```
> plddt_cutoff: Optional. The default value is 70.  
In the "output_dir", the "UniprotID.domains" files record the residue ranges of each domain for each protein. And each domain is named in the format "UniprotID_Dx.pdb", where x starts from 0.  

Example commands:
```
python ../parse_with_pae.py pdb_files pae_files domain_files
```
### 4.Fit domains into target densities.
Usage:
```
fit_with_chimerax.py domain_dir map_dir output_dir map_level resolution n_search n_thread
```
> map_level: Set the threshold of target densities to "map_level".  
> resolution: The resolution of the simulated densities of the input domains.  
> n_search: The number of initial placements of the fit model within the reference map.  
> n_thread: The number of domains fitted at the same time.  

This script will fit the domains in the "domain_dir" into each target density in the "map_dir" and save the results individually for each density in the "output_dir". The result files are named in the format "domain_.log" and contain information including the transformation matrix for the best match position. The parameters will be recorded in the "fit_config.txt" file.  

Example commands:
```
python ../fit_with_chimerax.py domain_files em_data fit_out 0.001 7 200 16
```

### 5.Score and rank.
Usage:
```
score_rank.py domain_dir map_dir fit_out_dir ref_map_threshold ref_map_laplacian_cutoff_low ref_map_laplacian_cutoff_high resolution n_thread box_num min_entries_per_box
```
> domain_dir: The folder containing domains not fitted into densities.  
> fit_out_dir: The "output_dir" used for the _"fit_with_chimerax.py"_ script.  
> ref_map_threshold: Voxel values below "ref_map_threshold" are ignored.  
> ref_map_laplacian_cutoff_low: After laplacian filtering, voxel values between "ref_map_laplacian_cutoff_low" and 0 are ignore.  
> ref_map_laplacian_cutoff_high: After laplacian filtering, voxel values  between 0 and "ref_map_laplacian_cutoff_high" are ignore.
> box_num: The number of bins into which the overlap volume is divided.
> min_entries_per_box: The minimum data points in each bin.

The "box_num" and "min_entries_per_box" parameters are used to caculate the local z-scores.  
The "ref_map_laplacian_cutoff_low" parameter is a negative value set to ignore some noise around 0. If the quality of the target density is quite good, you can also set this value to 0.  
The "ref_map_laplacian_cutoff_high" parameter is a positive value set to reduce the adverse effects at the density boundaries. It is better to set this parameter so that the positive-value surface surrounds the negative-value surface without exceeding it by too much.  
Pink: negative, yellow: positive.


<img src="/Example/figures/density_Pacrg_filtered.jpg" width="400px">

This script will create an "scores.txt" file in each density folder in the "fit_out_dir". And lines in the "scores.txt" file are sorted by local z-score and are in the following format:
> domain overlap_volume overlap_correlation local_z-score

Example commands:
```
python ../score_rank.py domain_files em_data fit_out 0.001 -0.0001 0.0001 7 16 10 50
```

### 6.Show results.
You can draw a 2D figure for the overlap volume and correlation.  
In this example, the results is plot below:

<img src="/Example/figures/score.jpg" width="400px">
In this figure, each point represents a domain. The cyan point standing out in the upper-right corner is a domain of Q9DAK2 (Pacrg).  

Besides, the local z-score of the cyan domain is 9.29, ranking first.  

Thus, DomainSeeker successfully identified the correct domain for the density of protein Pacrg.

## Additional Features
You can use the _"get_fitted_domain.py"_ script to generate the fitted domains from the output of fitting process.  
Usage:
```
get_fitted_domain.py domain_dir fitout_dir domain_list
```
> domain_dir: The directory of unfitted domains.
> fitout_dir: The output directory for the fitting results of a specific density.
> domain_list: A sequence of domain names, separated by spaces, without the '.pdb' suffix.

It will generate fitted domains named in the format "domainName_fitted.pdb", which will be saved in the "fitout_dir".

## Citation
If you use DomainSeeker in your work, please cite the following preprint:
```
Lu, Y., Chen, G., Sun, F., Zhu, Y., & Zhang, Z. (2024). De novo identification of protein domains in cryo-electron tomography maps from AlphaFold2 models. bioRxiv, 2024.2011.2021.623534. https://doi.org/10.1101/2024.11.21.623534 
```