#!/bin/bash

dir=/work/sph-huangj/analysis

M="bc_WBC bc_RBC bc_HGB bc_HCT bc_MCV bc_MCH bc_MCHC bc_RDW bc_PLT bc_PCT bc_MPV bc_PDW bc_LYMPH bc_MONO bc_NEUT bc_EO bc_BASO bc_NRBC bc_LYMPH_pt bc_MONO_pt bc_NEUT_pt bc_EO_pt bc_BASO_pt bc_NRBC_pt bc_RET_pt bc_RET bc_MRV bc_MSCV bc_IRF bc_HLR_pt bc_HLR bb_ALB bb_ALP bb_APOA bb_APOB bb_AST bb_BILD bb_BUN bb_CA bb_TC bb_CRE bb_CRP bb_CYS bb_GLU bb_HBA1C bb_HDL bb_IGF1 bb_LDL bb_LPA bb_oestradiol bb_PHOS bb_SHBG bb_TBIL bb_TES bb_TP bb_TG bb_UA bb_VITD"
M_list=`echo $M | sed 's/^/c("/; s/ /","/g; s/$/")/'`

for Y in covid lung_cancer cad vte2 t2dm; do #asthma breast_cancer covid lung_cancer cad stroke t2dm vte2; do # asthma breast_cancer cad cancer circulatory colorectal_cancer copd covid dm dvt gallbladder_cancer gastric_cancer haemorrhoid hepatobiliary_cancer i81 leukemia_cancer lung_cancer lymphoma_cancer melanoma_cancer metabolic mi myeloma_cancer oesophageal_cancer oropharyngeal_cancer ovarian_cancer pancreatic_cancer pe prostate_cancer ra renal_cancer respiratory stroke t1dm t2dm varicose vte vte2; do

	outdir=$dir/assoc.ind/$Y; mkdir -p $outdir
	for X in height bmi walk_pace abo.a abo.o abo.a_b abo.se vte.ff apoe.e4 sp1.s sp1.z; do
        if [ $X == $Y ]; then continue; fi
	label="$X-$Y"
	if [ -f $outdir/$label.log ]; then echo $label already run; continue; fi

	echo "#!/bin/bash
	module load python/anaconda3/2020.7
	source activate
	conda activate R
	cat /work/sph-huangj/scripts/main/assoc.ind.R | sed -e '17 s/?/$X/' -e '18 s/?/$Y/' -e '19 s/?/$M_list/' -e '20 s/?/$label/' > $label.R
	R CMD BATCH $label.R
	" > $outdir/$label.cmd
	cd $outdir
	bsub -q short -J ind.$label -o $label.LOG -e $label.ERR < $label.cmd

done
done

# label=top; echo "X M Y N X2Y.b X2Y.p X2M.b X2M.p M2Y.b M2Y.p TOT.b TOT.p ACME.b ACME.p ADE.b ADE.p Prop.b Prop.p" > $label.txt; cat */*-*.log | sed 's/"//g' | cut -d ' ' -f 3-8,10-11,13-14,16-23 | awk 'NF==18 && $13*$15 >0 && $12 <0.05 && $14<0.05 && $18 <0.05 && $17>0.1' | sort -k 17,17gr >> $label.txt

