#!/bin/bash

outdir=/work/sph-huangj/analysis/assoc

for dat in height; do # bmi bc_WBC bc_RBC bc_HGB bc_HCT bc_MCV bc_MCH bc_MCHC bc_RDW bc_PLT bc_PCT bc_MPV bc_PDW bc_LYMPH bc_MONO bc_NEUT bc_EO bc_BASO bc_NRBC bc_LYMPH_pt bc_MONO_pt bc_NEUT_pt bc_EO_pt bc_BASO_pt bc_NRBC_pt bc_RET_pt bc_RET bc_MRV bc_MSCV bc_IRF bc_HLR_pt bc_HLR bb_ALB bb_ALP bb_APOA bb_APOB bb_AST bb_BILD bb_BUN bb_CA bb_TC bb_CRE bb_CRP bb_CYS bb_GLU bb_HBA1C bb_HDL bb_IGF1 bb_LDL bb_LPA bb_oestradiol bb_PHOS bb_SHBG bb_TBIL bb_TES bb_TP bb_TG bb_UA bb_VITD; do
	
	echo "
	module load python/anaconda3/2020.7
	source activate
	conda activate R
	sed '29 s/M/$dat/2' /work/sph-huangj/scripts/assoc.bat.R > $dat.R
	R CMD BATCH $dat.R
	" > $outdir/$dat.cmd
	cd $outdir
	bsub -q smp -J $dat -o $dat.LOG -e $dat.ERR < $dat.cmd
done
