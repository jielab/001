
Arr1=("snp" "chr" "pos" "ea" "nea" "eaf" "n" "beta" "se" "p")
Arr2=("snp|rsid|variant_id" "chr|chrom|chromosome" "pos|bp|base_pair" "ea|eff.allele|alt|a1|allele1" "nea|ref.allele|ref|a2|allele0" "eaf|a1freq" "n" "beta" "se|standard_error" "^p|pval|p_bolt_lmm")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 数据准备：GWAS download and format
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Linux三剑客示例：cat download.list.txt | sed -r 's/^/wget /; s/(\w+)$/\1\/\1_buildGRCh38.tsv.gz/' | awk '{cnt=int(NR/100); print $0 > "download"cnt".sh"}'
dir=/mnt/d
outdir=$dir/data/gwas/pheweb2
files=$dir/files/mr.list.txt
cut -f 1,4 $files | sed '1d' | fgrep -vw "NA" | while IFS=$'\t' read -r pheno_code file_url; do
	file_name=$(basename $file_url)
	file_ext=${file_name##*.}
	if [[ -f $pheno_code.gz ]]; then
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
	" > $pheno_code.cmd
	chmod 777 $pheno_code.cmd; ./$pheno_code.cmd 
done
for dat in `ls *gz | sed 's/\.gz$//g'`; do
	head_row=`zcat $dat.gz | head -1 | sed 's/\t/ /g'`;
	for i in ${!Arr1[@]}; do
		eval ${Arr1[$i]}=`echo $head_row | tr ' ' '\n' | grep -Einw ${Arr2[$i]} | sed 's/:.*//'`
	done
	echo $snp $chr $pos $ea $nea $eaf $n $beta $se $p
	echo "#!/bin/bash
	zcat $dat.gz | awk '{if ($p_col<1e-300) $p_col=1e-300; if (NR==1) print \"SNP CHR POS EA NEA EAF N BETA SE P\"; else print $snp,$chr,$pos, toupper($ea), toupper($nea),$eaf,$n, $beta,$se,$p}' | sort -k 2,2V -k 3,3n | gzip -f > $dat.gz
	" > $dat.cmd
	chmod 777 $dat.cmd; ./$dat.cmd 
done


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# C1: Correlation 主要指在基因水平上
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
conda activate ldsc # 采用 python2，而不是 python3
dir=/mnt/d/data/gwas/penguin
ldscdir=/mnt/d/data/ldsc
dats=`cd $dir; ls *.gz | sed 's/\.gz//g'`
for dat in $dats; do
	echo $dat
	datf=$dir/$dat.gz
	head_row=`zcat $datf | head -1 | sed 's/\t/ /g'`;
	munge_sumstats.py --chunksize 10000 --sumstats $dat.tmp.gz --merge-alleles $ldscdir/w_hm3.snplist --out $dat --snp SNP --a1 EA --a2 NEA --ignore SNPID,OR --signed-sumstats BETA,0 --p P --N 100000 
done
for dat in $dats; do
	echo $dat $dats | sed -e 's/ /.sumstats.gz,/g' -e 's/$/.sumstats.gz/' | xargs -n1 -I % ldsc.py --rg % --out $dat.rg --ref-ld-chr $ldscdir --w-ld-chr $ldscdir
	beginl=`awk '$1=="Summary" {printf NR}' $dat.rg.log` 
	awk -v s=$beginl 'FNR >s' $dat.rg.log | head -n -3 | sed 's/.sumstats.gz//g'  > $dat.rg.txt
done
awk 'NR==1 || FNR>1' *.rg.txt > all.rg.res


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# C2: Causation 主要是通过MR分析
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
exp_dir=/mnt/d/data/gwas/pheno
out_dir=/mnt/d/data/gwas/pheno
exp_files=`cd $exp_dir; ls *.gz | sed 's/\.gz//g'`
out_files=`cd $out_dir; ls *.gz | sed 's/\.gz//g'`
for exp in $exp_files; do
	exp_file0=$exp_dir/$exp.gz; exp_pheno=$exp
	head_row=`zcat $exp_file0 | head -1 | sed 's/\t/ /g'`
	for i in ${!Arr1[@]}; do
		eval ${Arr1[$i]}=`echo $head_row | tr ' ' '\n' | grep -Einw ${Arr2[$i]} | sed 's/:.*//'`
	done
	echo exp $exp, snp $snp, chr $chr, pos $pos, p $p
	exp_iv_file=$exp_dir/$exp.top.snp
	if [ ! -e $exp_iv_file ]; then
		exp_iv_file=$exp.top.snp # 这是本地的 $exp.top.snp 文件，不是上面的 $exp_dir/$exp.top.snp
		if [ ! -e $exp_iv_file ]; then
			zcat $exp_file0 | awk -v snp=$snp -v chr=$chr -v pos=$pos -v p=$p 'NR==1 || $p<=5e-8 {posm=int($pos/1e6); print $snp, $chr, $pos, $p, posm}' | sort -k 2,2n -k 5,5n -k 4,4g | awk '{if (arr[$NF] !="Y") print $1; arr[$NF] ="Y"}' > $exp.top.snp
		fi
	fi
	exp_iv_line=`wc -l $exp_iv_file | awk '{printf $1}'`
	if [ "$exp_iv_line" -le 1 ]; then # no IV SNP
		continue 
	fi
	exp_file=$exp.txt
	if [ ! -e $exp.txt ]; then
		zcat $exp_file0 | fgrep -wf $exp_iv_file | awk -v snp=$snp -v ea=$ea -v nea=$nea -v beta=$beta -v se=$se -v p=$p BEGIN'{print "SNP EA NEA BETA SE P"}{print $snp, toupper($ea), toupper($nea), $beta, $se, $p}' > $exp.txt	
	fi
	for out in $out_files; do
		if [ "$exp" = "$out" ]; then continue; fi
		if [ -e $exp.$out.Rout ]; then #already run
			continue
		fi
		echo Running: $exp $out 
		out_file0=$out_dir/$out.gz; out_pheno=$out
		head_row=`zcat $out_file0 | head -1 | sed 's/\t/ /g'`
		for i in ${!Arr1[@]}; do
			eval ${Arr1[$i]}=`echo $head_row | tr ' ' '\n' | grep -Einw ${Arr2[$i]} | sed 's/:.*//'`
		done
		echo out $out, snp $snp, ea $ea, nea $nea, beta $beta, se $se, p $p
		out_file=$exp.list.$out.txt
		if [ ! -e $out_file ]; then
			zcat $out_file0 | fgrep -wf $exp_iv_file | awk -v snp=$snp -v ea=$ea -v nea=$nea -v beta=$beta -v se=$se -v p=$p BEGIN'{print "SNP EA NEA BETA SE P"}{print $snp, toupper($ea), toupper($nea), $beta, $se, $p}' > $exp.list.$out.txt
		fi
		echo "
		source('/mnt/d/scripts/library/tsmr.f.R')
		tsmr_fn( exp_iv_file='$exp_iv_file',
			exp_name='$exp', exp_pheno='$exp_pheno', exp_ieu='NA', exp_file='$exp_file',  
			out_name='$out', out_pheno='$out_pheno', out_ieu='NA', out_file='$out_file'
		)
		" > $exp.$out.R
		R CMD BATCH $exp.$out.R
	done
done
cat *.tsmr.out | sed 's/|/\t/g' > ../pheno.res
# touch jie.txt; for f in *.tsmr.out; do cat $f >> jie.txt; done

	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# C3: Colocalization 从全局到局部local
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cut -f 1-3,10 x.education | awk '{if (NR==1) print $0,"ChrPos-m"; else {b=int($3/1e6); print $0,"chr"$2":"b"m"}}' > exp.chrpos
cut -f 1-3,10 x.education | awk 'NR==1 || $4<5e-8 {b=int($3/1e6); print $1,$2,$3,$4,b}' | sort -k 2,2n -k 5,5n -k 4,4g | awk '{if (arr[$NF] !="Y") print $1,"chr"$2":"$5"m"; arr[$NF] ="Y"}' > exp.top.loci
python /mnt/d/scripts/library/join_file.py -i "exp.chrpos,TAB,0 exp.top.loci,SPACE,0" -o exp.tmp
sort -k 2,2n -k 5,5V -k 4,4g exp.tmp | awk '{if ($7 !="NA") seg[$7]="Y"; if (NR==1 || seg[$5]=="Y") print $0}' | sort -k 2,2n -k 3,3n | sed 's/ /\t/g' | cut -f 1-6 > exp.txt
cut -f 1 exp.txt > exp.snp
python /mnt/d/scripts/library/join_file.py -i "exp.snp,TAB,0 x.education,TAB,0" -o exp.tmp2
paste exp.txt exp.tmp2 | sed '1 s/BETA/BETA.exp/; 1 s/SE/SE.exp/; s/ /\t/g' | cut -f 1,5,9-12,15-17 > exp.txt2 
for dat in `cd /mnt/d/data/gwas/penguin; ls y.*.gz | sed 's/\.gz//g'`; do
	echo $dat
	python /mnt/d/scripts/library/join_file.py -i "exp.snp,TAB,0 $dat,TAB,0" -o $dat.tmp
	cut -d " " -f 1,5,6,9,10 $dat.tmp | sed "1 s/ /\.$dat /g; 1 s/$/\.$dat/" > $dat.joined
done
paste exp.txt2 y.*.joined > all.coloc.tmp
awk '{print NF}' all.coloc.tmp | sort -nu
num=`ls -l y.*.joined | wc -l | awk '{printf $1}'`
awk -v num=$num '{for (i=1;i<=num;i++) {if (NR>1 && $(i*5+5) !=$1) print $1, $i, "ERR"; else if (NR>1 && $(i*5+6) !=$5 && $(i*5+6) !="NA") $(i*5+6) =-$(i*5+6); $(i*5+5)=""; $(i*5+6)=""; $(i*5+7)="" }; print $0}' all.coloc.tmp > all.coloc.txt; fgrep ERR all.coloc.txt


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# C4: Coevolution mRNA及蛋白质互作
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 参考 2023. On the Decoupling of Evolutionary Changes in mRNA and Protein Levels


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# C5: Community 共同体
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# A community with a shared future for mankind