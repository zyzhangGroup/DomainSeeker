import os

import tqdm

from powerfit_em.powerfit import setup_target, setup_template_structure, setup_rotational_matrix
from powerfit_em.powerfitter import PowerFitter
from powerfit_em.analyzer import Analyzer
from powerfit_em.helpers import write_fits_to_pdb

target_volume = 'em_data/A.mrc'
resolution = 6.5
no_resampling = True
resampling_rate = 2
no_trimming = False
trimming_cutoff = None

with open(target_volume, 'rb') as target_in:
    target = setup_target(
        target_in, 
        resolution, 
        no_resampling, 
        resampling_rate, 
        no_trimming, 
        trimming_cutoff
        )


template_structure = 'domain/Q6X6Z7_D0.pdb'
chain = None
# target = target_volume
# resolution =  
core_weighted = False

with open(template_structure, 'rb') as template_in:
    structure, template, mask, z_sigma = setup_template_structure(
        template_in, 
        chain, 
        target, 
        resolution, 
        core_weighted
    )

angle = 2

rotmat = setup_rotational_matrix(angle)

queue = None
nproc = 10
laplace = False

pf = PowerFitter(target, rotmat, template, mask, queue, nproc, laplace=laplace)
pf.scan(progress=tqdm.tqdm)


analyzer = Analyzer(
        pf.lcc,
        rotmat,
        pf.rot,
        voxelspacing=target.voxelspacing,
        origin=target.origin,
        z_sigma=z_sigma,
    )

analyzer.tofile('correct.solutions')

outputdir = 'powerfit_fits/'
if not os.path.exists(outputdir):
    os.mkdir(outputdir)
write_fits_to_pdb(structure, analyzer.solutions[:10], basename=os.path.join(outputdir, "out"))