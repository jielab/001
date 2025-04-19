#!/bin/bash
module load python/anaconda3/2020.7
source activate
conda activate gsMap

gsmap quick_mode --workdir '/work/sph-huangj/testjobs/0415/test2' --homolog_file '/work/sph-huangj/data/gsMap/gsMap_resource/homologs/mouse_human_homologs.txt' --sample_name 'E16.5_E1S1.MOSTA' 	--gsMap_resource_dir '/work/sph-huangj/data/gsMap/gsMap_resource' --hdf5_path '/work/sph-huangj/data/gsMap/gsMap_example_data/ST/E16.5_E1S1.MOSTA.h5ad' 	--annotation 'annotation' --data_layer 'count' --sumstats_file '/work/sph-huangj/data/gsMap/gsMap_example_data/GWAS/IQ_NG_2018.sumstats.gz' --trait_name 'IQ'


