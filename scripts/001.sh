#!/bin/bash

dir=/restricted/projectnb/ukbiobank/jiehuang
datdir=$dir/data/gwas/pheweb
outdir=$dir/analysis/mr
mrlist=$dir/files/mr.list.txt


# download files
cut -f 3,4 $mrlist | sed '1d' | fgrep -vw "NA" | while IFS=$'\t' read -r pheno_code file_url; do
	file_name=$(basename $file_url)
	file_ext=${file_name##*.}
	if [[ -f $datdir/$pheno_code.gz ]]; then
		echo $pheno_code.gz already downloaded
		continue
	fi
	echo -e "#!/bin/bash
	wget $file_url
	mv $file_name $pheno_code.$file_ext
	if [[ $file_ext == zip ]]; then
		unzip -p $pheno_code.$file_ext | gzip -f > $pheno_code.gz
	fi
	if [[ $file_ext == txt || $file_ext == tsv ]]; then
		mv $pheno_code.$file_ext $pheno_code
		gzip -f $pheno_code
	fi
	" > $datdir/$pheno_code.cmd
done


cut 1-3 $mrlist | sed '1d' | while IFS=$'\t' read -r ex_name ex_ieu ex_gcat; do
	newdir=$outdir/$ex_name; mkdir -p $newdir
	top_snps_file=$dir/data/gwas/posCtrls/$ex_name.top.snps
	cut 1-3 $mrlist | sed '1d' | while IFS=$'\t' read -r ou_name ou_ieu ou_gcat; do
		if [[ $ex_name == $ou_name ]]; then
			echo $ex_name equal $ou_name
			continue
		fi
		if [[ -f $newdir/$ex_name.$ou_name.res.txt ]]; then
			echo $ex_name $ou_name already done
			continue
		fi
		echo -e "
		source('mr.f.R')
		mrfn(indir='$datdir', ex_name='$ex_name', ex_ieu='$ex_ieu', ex_gcat='$ex_gcat', ou_name='$ou_name', ou_ieu='$ou_ieu', ou_gcat='$ou_gcat', top_snps_file='$top_snps_file')
		" > $newdir/$ex_name.$ou_name.mr.R
		Rscript $ex_name.$ou_name.mr.R
	done
done
