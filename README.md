## 🧬1. 获取数据

### 这是初二生物学课本里面的一页。

![middle school](./images/middle.jpg)

* ### 📍1.1 [HAPMAP3 genotype 数据](https://www.broadinstitute.org/medical-and-population-genetics/hapmap-3), 1百多万个SNP，一般作为 LD 计算的 reference panel。

* ### 📍1.2 [千人基因组项目数据](https://www.internationalgenome.org/data)， 将近1亿个SNP，一般作为 imputation 的 reference panel。

* ### 📍1.3 UKB 数据，现在推荐用[UKB RAP](https://dnanexus.gitbook.io/uk-biobank-rap)。
```
1. 最新数据通知 https://community.ukbiobank.ac.uk/hc/en-gb/articles/26088595922333-New-and-Updated-Data
2. UKB RAP：https://ukbiobank.dnanexus.com/landing
```
<br/>


## 🧬2. GWAS

![GWAS](./images/GWAS.jpg)

* ### 📍2.1 GWAS数据获取，最经典的是[GWAS Catalog](https://www.ebi.ac.uk/gwas)。

* ### 📍2.2 GWAS数据QC示例
```
1. 确保文件每一行的列数目是一样的。将连续空格中插入NA，扣好第一粒纽扣。
   zcat GWAS.gz | awk '{print NF}' | sort -nu | wc -l 
   sed 's/^\t/NA\t/; s/\t\t/\tNA\t/g; s/\t\t/\tNA\t/g; s/\t$/\tNA/'。
2. 按照CHR和POS排序，否则pheweb 会报错
   sort -k 1,1V -k 2,2n # -V是为了把chrX和chrY排到最后，但是需要把第一行先写到新文件里。
3. GWAS数据本身的问题：
   (1) Allele 最好是大写，awk 和 R 都有 toupper()功能。
   (2) P值最好不要小于1e-312，awk 会把其当成0，有一些软件（比如LDSC）也会报错，这个时候要么用Z值，要么人为将这些P值设为1e-300。
   (3) BETA|SE|P出现“三缺一” 的情况： b = se * qnorm(p/2); se = abs(b/qnorm(p/2)); se = (CI_upper - CI_lower)/(1.96*2); p = 2*pnorm(-abs(b/se))

对检查没问题的GWAS，深加工示例：

1. 🚜liftOver 
   dat=XYZ; head -1 $dat.txt > $dat.sorted; tail -n +2 $dat.txt | sort -k 1,1V -k 2,2n > $dat.sorted
   python ~/scripts/f/add_rsid.py -i $dat.sorted --sep "\t" --chr CHR --pos POS --ref NEA --alt EA -d ~/data/dbsnp/rsids-v154-hg19.tsv.gz -o $dat.tmp1
   cat $dat.tmp1 | awk 'NR >1 {print "chr"$1, $2 -1, $2, $9}' | sed 's/^chr23/chrX/' > $dat.tolift
   liftOver $dat.tolift /work/sph-huangj/files/hg19ToHg38.over.chain.gz $dat.lifted $dat.unmapped
   cut -f 3,4 $dat.lifted > $dat.pos_snp
   python ~/scripts/f/join_file.py -i "$dat.tmp1,TAB,8 $dat.pos_snp,TAB,1" -o $dat.tmp2
   cut -d " " -f 1-10 $dat.tmp2 | sed '1s/POS/POS.37/; 1s/NA/POS/' | gzip -f > clean/$dat.gz
2. ⛄跟其他数据合并 
   python scripts/library/join_file.py -i "$dat,TAB,0 $dat.lifted.3col,TAB,2" -o $dat.NEW.tmp
   sed -i 's/  */\t/g' $dat.NEW.tmp; awk '$NF=="NA"' $dat.NEW.tmp | wc -l
   cut -f 1-10,12 $dat.NEW.tmp | sed '1 s/POS/POS.b38/' > $dat.NEW.txt
3. 🧢添加 rsID 
   用户也可以在pheweb网站下载 rsids-v154-hgXX.tsv.gz 文件（7亿多行）后，在本Github的 scripts文件夹下载本课题组修订的 add_rsid.py; dos2unix add_rsid.py，然后运行如下示例命令。
   python add_rsid.py -i test.tsv --sep "\t" --chr CHR --pos POS --ref NEA --alt EA -d data/dbsnp/rsids-v154-hg19.tsv.gz -o out.tsv
4. 🏃瘦身 ‍
   zcat $dat.gz | awk '{if (NR!=1) {$5=sprintf("%.4f",$5); $6=sprintf("%.4f",$6)} print $0}' | bgzip > $dat.lean.gz
5. 🔍索引 
   tabix -f -S 1 -s 1 -b 2 -e 2 GWAS.gz
6. 本课题组建议GWAS列名称： SNP CHR POS EA NEA EAF N BETA SE P。
```

* ### 📍2.3 GWAS数据可视化
> 密西根大学开发的[Pheweb](https://github.com/statgen/pheweb) ，UKB的几千个GWAS的数据放在[pheweb.org](https://pheweb.org/)上，[中国版本](https://pheweb.ckbiobank.org/)，[日本版本](https://pheweb.jp/)。
> Pheweb有一个强大的add_rsids.py 的功能，但是存在先天缺陷。根据该[聊天记录](https://github.com/statgen/pheweb/issues/217)，用户可以在安装pheweb 后找到 add_rsids.py 文件（find /home/ -name "add_rsid*" 或者 pip show --files pheweb），修改一行代码（第140行）。
> 密西根大学还开发了[locuszoom](http://locuszoom.org/) 实现基因组局部地区的可视化🔍。 
```
修改前：rsids = [rsid['rsid'] for rsid in rsid_group if cpra['ref'] == rsid['ref'] and are_match(cpra['alt'], rsid['alt'])]
修改后：rsids = [rsid['rsid'] for rsid in rsid_group if (cpra['ref'] == rsid['ref'] and are_match(cpra['alt'], rsid['alt'])) or (cpra['ref'] == rsid['alt'] and are_match(cpra['alt'], rsid['ref']))]
```
<br/>


## 🧬3. MR
> 如果有个体数据，可以用 [OneSampleMR包](https://cran.r-project.org/web/packages/OneSampleMR/index.html)。
> 如果只有已发表的summary数据，就可以使用Bristol大学开发的[TwoSampleMR R包](https://mrcieu.github.io/TwoSampleMR/index.html)或剑桥大学团队开发的[MendelianRandomization R包](https://wellcomeopenresearch.org/articles/8-449)。
> 工具变量，一般需要去掉 F_stats <10 或者位于 <b>[MHC区间]</b> 【chr6:28477897-33448354 [(GRCh37)](https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh37), chr6:28510120-33480577 [(GRCh38)](https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC)】 的SNP。
> <b>密西根大学</b>开发的 [imputation server](https://imputationserver.sph.umich.edu) 用的是： 从rs57232568 【29000554 (版本37), 29032777 (版本38)】 到 rs116206961【33999992 (版本37), 34032215 (版本38)】
> 10个注意事项示例：
```
> 1. 用 allele.qc，协调两组数据的 BETA和EAF，但是输出文件依然是原来的EA和 NEA。
> 2. 用 fread和fwrite 比 read.table 和 write.table 更快，但fwrite默认输出带quote。
> 3. EAF如果是0，fwrite 会将其视为NA，输出后TXT就是一个空格，导致数据少了一列。
> 4. gcta –cojo-slct 生成的 .ldr.cojo 文件最后多出一列TAB。
> 5. 添加代码，处理两个输入GWAS，其中一个或者两个都不存在 EAF 和 N 的问题。
> 6. 连着用几个 %>%，缩减代码，最后实际上跑偏了、失控了。
> 7. 用 group() 之后，如果后面没有跟上 ungroup()，后面会出问题。
> 8. 用 plink，.bim 文件中的染色体可以用 X表示。但是 gcta 必须用 23。
> 9. 注意HLA的GRCh37或GRCh37的确切 chr:start-end。
>10. 最后可能发现软件跑出来的结果有重大问题，不“鲁棒” https://github.com/ZhaotongL/cisMRcML/issues/6
```
<br/>


## 🧬4. 参考资料

基因注释信息🔍
```
> - dbSNP: https://www.ncbi.nlm.nih.gov/snp
> - UCSC genome browser: https://www.genome.ucsc.edu
> - 美国精准医学All of Us： https://databrowser.researchallofus.org 
> - 密西根大学公卫学院 TopMed browser: https://bravo.sph.umich.edu
> - 一天发了7篇 NATURE系列文章的Gnomad项目 browser: https://gnomad.broadinstitute.org
```

🐎GWAS-PRS-MR ”三驾马车“ 入门指南  
▸ GWAS入门： 2021. Nature RMP. [Genome-wide association studies](https://www.nature.com/articles/s43586-021-00056-9)  
▸ 🏮GWAS详解中文版：[gwaslab.org](https://gwaslab.org)，比如其中提到 [REGENIE](https://gwaslab.org/2021/03/28/regenie/)  
▸ PRS入门. Nature Protocols. [Tutorial: a guide to performing polygenic risk score analyses](https://www.nature.com/articles/s41596-020-0353-1)  
▸ MR入门： 2022. Nature RMP. [Mendelian randomization](https://www.nature.com/articles/s43586-021-00092-5)  

🛵R  
▸ 在“环境变量”里 将 R_LIBS_USER 设为 D:\software_win\R_lib  
▸ HPC上的环境变量位于 /work/sph-huangj/.conda/envs/R4.4.2  
▸ R画图集锦: https://r-graph-gallery.com/index.html  
▸ R新冠地图: https://statsandr.com/blog/top-r-resources-on-covid-19-coronavirus/  
▸ 供复现代码： https://globalenvhealth.org/code-data-download/  
▸ 🏮顾祖广炫酷生信图： [https://jokergoo.github.io/software/](https://jokergoo.github.io/software/)  
▸ 🏮梁志生R包荟萃🎇 [https://gitee.com/sheng0825/projects](https://gitee.com/sheng0825/projects)  



🤖Linux
```
> 安装: 遇到 press any key to continue，用管理员权限打开cmd, 运行 netsh winsock reset
> HPC后台运行： nohup ./assoc.sum.sh & 之后 ps aux | grep assoc.sum.sh 之后 kill
> HPC硬盘额度：mmlsquota -g sph-huangj --block-size auto
> 后台多线程下载: screen -dmS jack aria2c -x 4 -i url.txt --log-level=info --log=jack.log; screen -ls; screen -S jack -X quit 
> 三剑客🗡代码示例: awk '{cnt=int(NR/100); print $0 > "download"cnt".sh"}'
```

<br/>
<br/>
