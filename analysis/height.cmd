
	module load gcc/11.2.0
	module load intel/2018.4

/work/sph-huangj/apps/regenie/regenie --bt --bsize 1000 --lowmem --lowmem-prefix tmp_preds --step 1 --strict --pgen /work/sph-huangj/data/ukb/gen/typ/ukb --phenoFile /work/sph-huangj/data/ukb/phe/common/ukb.pheno --phenoColList height --covarFile /work/sph-huangj/data/ukb/phe/common/ukb.pheno --covarColList age,sex,PC1,PC2 --out height.step1

