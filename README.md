<br/>
<br/>

# #1. 下载和处理国际上公用公开的数据

## 这是初二生物学课本里面的一页。

![middle school](./images/middle.jpg)

## #1.1 HAPMAP3 genotype 数据, 一般作为 LD 计算的 reference panel

> 打开 <https://www.broadinstitute.org/medical-and-population-genetics/hapmap-3>， 
> 点击 How To Download This Release 下面的 A. SNP Genotype Data 段落的中间3个链接。
> 文件名字里面有 "b36"，现在一般都用 b37（比如 UK Biobank），甚至有的用 b38，
> 所以下载后解压后需要将那个 .map 文件先用 liftOver 转化为 b37 格式，然后用 PLINK 生成 bed/bim/fam 文件。
> 这个基因数据可供 LDSC 和 GSMR 等软件使用。
<br/>

## #1.2. 1000 genomes (千人基因组) genotype 数据， 一般作为 imputation 的 reference panel.

> 在[千人基因组官网](https://www.internationalgenome.org/data) 下载 Phase 3 对应的 VCF 链接，
> 有一个文件罗列了每一个样本的人群（pop）和人种 (super_pop)，以及性别，可以用PLINK --keep 选取特定人种的样本。
> 下载下来的数据，有将近一个亿的SNP，每个染色体都是单独的文件。后续跑 GWAS 或提取 PRS 的时候，也是每条染色体的数据分开来跑。
> PLINK的网站上也有“1000 genomes phase3” 数据。PLINK 不允许 SNP 名字有重复，可以用下面的命令来处理。
> > - awk '{if(array[$2]=="Y") {i++; $2=$2".DUP"i}; print $0; array[$2]="Y"}' chr1.bim.COPY > chr1.bim 
<br/>


# #2.  提取 UKB 表型数据
从 [ukbiobank](ukbiobank.ac.uk) 官网，点击 data showcase --> Essential information --> Accessing UK Biobank data，阅读 Data access guide 文件，里面会提到如何下载用来下载UKB数据的小软件（比如ukbunpack和unkconv）。
可以通过 Search 或 Browse（截图如下）去熟悉 UKB 里面的数据结构。具体的数据，需要申请得到批准后，从最上面的 Researcher log in 登录后获取，见百度网盘上的 ukb50136.enc。
![UKB](./images/ukb.png)


## #2.1 提取简单的表型数据（比如 age, sex, race, bmi, etc.）

> WINDOWS电脑建议安装系统自带的 Ubuntu Linux系统，然后用 cd /mnt/d/ （而不是 D:/）进入 D 盘。
> 1. 执行 ukbmd5 ukb50136.enc, 确认得到 981b47f85c6b2fb849320c7a3559ba23，确保数据完整。
> 2. 执行 ukbunpack 解压、解密数据：ukbunpack ukb50136.enc fd5b09bb6f6dbXe9e0774e91f81Xaf1a2d9e3aeaXddad42a2ff6ea60f2f23XX8
> 3. 执行 ukbconv，提取需要的列，并生成想要的格式：ukbconv ukb50136.enc_ukb r -iMY.fields.id -oMY
> 可以先写一个 MY.fields.txt 文件，列出想提取的变量和对应的 data-field，比如 21022 age ... ...
> 然后手动或者用下面的命令，提取出该文件的第一列（需要提取的表型的ID），并确认没有重复的ID。
>   - awk '{print $1}' MY.fields.txt > MY.fields.id; sort MY.fields.id | uniq -d
> 4. 参考 scripts/phe.R 代码，通过上面生成的 MY.r 读入上面生成的 MY.tab 数据，并且给每个变量赋予上述 MY.fields.txt给出的名字。需要注意 MY.r 第一行里面的路劲正确。
> 注：> 对于表型数据的提取，有一个 [ukbtools R软件包](https://kenhanscombe.github.io/ukbtools/articles/explore-ukb-data.html)。 但不是太好用，并且很慢，可供参考。

<br/>

## #2.2 提取跨越很多列的数据，比如 ICD (data field 42170）

> ICD 这样的指标，包含了很多不同时间的时间点，量很大，建议分开来处理。
> > ukbconv ukb42156.enc_ukb r -s42170 -oicd
> > sed -i 's/"//g icd.tab

> 将 icd.tab 文件整合为两列，便于读入R。
> > cat icd.tab | sed -e 's/\tNA//g' -e 's/\t/,/2g' | \
> > awk '{ if(NR==1) print "IID icd"; else if (NF==1) print $1 " NA"; else print $0"," }' > icd.2cols

<br/>

## #2.3. 对表型数据进行 GWAS 运行之前的处理

> 提取需要研究的表型数据和相关的covariates，比如 age, sex, PCs。
> 一般来说，quantitative的表型数据要 adjust for covariates 和转化成正态分布，这个可以在R里面用下面的命令来实现。
> trait_res = residuals(lm(trait ~ age+sex+PC1+PC2, na.action=na.exclude)
> trait_inv = qnorm((rank(trait_res,na.last="keep")-0.5) / length(na.omit(trait_res)))

> 对于疾病的binary 表型，只需要把需要 adjust 的covarites 和表型数据放在同一个表型数据文件里面，
> 然后在 GWAS里面的plink命令指明哪个是表型，哪些是 covariates。
> ![compareB](./images/T2D.Z.png) 
<br/>

## #2.4 关于UKB 基因数据
> UKB的基因数据很大，所有申请者都能得到一样的数据（样本的ID不一样），一般下载到服务器上去储存和使用。

<br/>


# #3. GWAS 运行和 结果 “挖掘”

![GWAS](./images/GWAS.jpg)

<br/>

## #3.1 专人在服务器上运行

> 目前GWAS 由专人负责运行，一般来说就是通过下面这样的PLINK命令来跑
> > for chr in {1..22}; do
> > - plink2 --memory 12000 --threads 16 --pfile chr$chr --extract ukb.chr$chr.good.snps --pheno cvd.EUR.pheno --no-psam-pheno --pheno-name XXX --1 --glm cols=+ax,+a1freq,+a1freqcc,+a1count,+a1countcc,+beta,+orbeta,+nobs hide-covar no-x-sex --covar pheno/ukb.cov --covar-name age,sex,PC1-PC10 --out chr$chr   
> > done

> 上述命令顺利跑完后，确认生成的文件没有问题后，可以把所有的染色体的数据串到一起，形成一个单一的 XXX.gwas.gz 文件。
最终合并成的 XXX.gwas.gz 文件用 TAB 分割，CHR:POS 排好序，要不然 LocusZoom 那样的软件不能处理。也可以用 tabix -f -S 1 -s 1 -b 2 -e 2 XXX.gwas.gz 对数据进行索引，便于 LocalZoom 那样的软件去处理。
<br/>

## #3.2 公开的GWAS数据进行练手，或对比

> 最经典的，起源于美国NIH 的 [GWAS Catalog](https://www.ebi.ac.uk/gwas). 这个页面也罗列了一些大型GWAS数据联盟。
> 欧洲版本，不需要下载就能通过 TwoSampleMR 远程读入。他们提倡 使用 VCF 格式的GWAS文件。
![Figure IEU](./images/ieu-open.png)

> UKB GWAS 完整的分析结果，网上发布
> > - 美国哈佛大学：http://www.nealelab.is/uk-biobank 
> > - 英国爱丁堡大学：geneatlas: http://geneatlas.roslin.ed.ac.uk
> > - 哈佛大学的 CVD knowlege portal: https://hugeamp.org/

<br/><br/>


## #3.3 GWAS的显示和注释，使用密西根大学开发的[Pheweb](https://github.com/statgen/pheweb) 流水线作业，日本人就用这个弄出了 [pheweb.jp](pheweb.jp)。密西根大学还开发了 [LocusZOOM](http://locuszoom.org), 具有类似和互补的功能。

> ### 如果不用上述的系统，也可以用 [PLINK](https://www.cog-genomics.org/plink/1.9/) 人工操作。点击左边菜单中的 Report postprocess 中的 3个命令（--annotate, --clump, --gene-report）

```
trait=MI

gunzip -c $trait.gwas.gz | sed '1 s/ POS/ BP/' > $trait.gwas.txt # 以后就不需要 sed 这一步了

plink --annotate $trait.gwas.txt NA ranges=glist-hg19 --border 10 --pfilter 5e-8 --out $trait.top

# 由于千人基因组 (g1k) 的基因数据过大（将近1亿个SNP），一般讲每一个染色体的GWAS数据分开来 clump
# plink clump 的结果，不包括那些 --bfile 里面没有的SNP，所以得要把那些SNP再添加到 clump 的结果里。
# 已沟通，可惜 PLINK的作者不想让 PLINK 来直接处理这个问题，
for chr in {1..22}; do
   plink1.9 --vcf g1k.chr$chr.vcf.gz --clump $trait.gwas.txt --clump-p1 5e-08 --clump-p2 5e-08 --clump-kb 1000 --clump-r2 0.2 --out $trait.chr$chr
   awk '$1 !="" {print $3,$1, $4,$5}' $trait.chr$chr.clumped > $trait.chr$chr.top
done

# 通过LD的计算来找到GWAS数据里面的independent top hits，也有一些问题（比如g1k的LD不是金标准，r2也不是最合理的筛选办法），并且计算量很大。 
# 如果不考虑 SNP之间的LD，只考虑距离，可以用下面这个简单的代码来寻找GWAS数据里面每1MB区间的top SNP。
# 假设GWAS的第1，2，3 列分别是 SNP, CHR, POS，最后一列是P。

zcat ABC.gwas.gz | awk 'NR==1 || $NF<5e-8 {b=sprintf("%.0f",$3/1e6); print $1,$2,$3,$NF,b}' | \
	sort -k 2,2n -k 5,5n -k 4,4g | awk '{if (arr[$NF] !="Y") print $0; arr[$NF] ="Y"}' 

# 要把上述得到的显著区域跟别人已经发表的 SNP进行比较，看是不是有重叠（1MB范围之内的重叠都算），可以用 bedtools 命令。
```
<br/>


## #3.5 GWAS 文件功能（functional）分析 

> ### 可阅读2020年的文章[From GWAS to Function: Using Functional Genomics to Identify the Mechanisms Underlying Complex Diseases](https://www.frontiersin.org/articles/10.3389/fgene.2020.00424/full)
> ### 可先尝试傻瓜相机式的[FUMA](https://fuma.ctglab.nl/) 网上解读系统

> ![FUMA](./images/fuma.png) 

<br/>

## #3.6 多基因风险评分PRS

> ### 相关的方法学，请参考经典版的 [PLINK0.9](http://zzz.bwh.harvard.edu/plink/profile.shtml) 和新版的 [PLINK1.9](https://www.cog-genomics.org/plink/1.9/score) 

<br/>
<br/>


## #3.7. 因果分析 Mendelian Randomization
> 如果有个体数据数据，可以采用 [GSMR](https://cnsgenomics.com/software/gcta/#GSMR)，参考[GSMR 文章](https://www.nature.com/articles/s41467-017-02317-2) 。
> 如果没有个体数据，只有别人报道的 exposure 和 outcome 的 BETA 和 SE，就可以使用 [MendelianRandomization R包](https://wellcomeopenresearch.org/articles/5-252/v2)，或 [TwoSampleMR R包](https://mrcieu.github.io/TwoSampleMR/index.html)。


<br/>
<br/>

# # 参考文献和网站

基因注释信息浏览器：
> - dbSNP: https://www.ncbi.nlm.nih.gov/snp/   
> - UCSC genome browser: https://www.genome.ucsc.edu/ 
> - 美国精准医学All of Us：https://www.researchallofus.org/ 和 https://databrowser.researchallofus.org/   
> - TopMed browser: https://bravo.sph.umich.edu/ 
> - Gnomad browser: https://gnomad.broadinstitute.org/ 
> - GlobalBiobankEngine：https://github.com/rivas-lab 


GWAS-PRS-MR 入门：
> GWAS:
> > - Y 2021. Nature Reviews Methods Primers. [Genome-wide association studies](https://www.nature.com/articles/s43586-021-00056-9)
> > - Y 2006. Nature Review Genetics. [A tutorial on statistical methods for population association studies](https://pubmed.ncbi.nlm.nih.gov/16983374/)
> > - Y 2020. NEJM. Genomewide Association Study of Severe Covid-19 with Respiratory Failure (https://www.nejm.org/doi/full/10.1056/NEJMoa2020283)
> > - 芬兰赫尔辛基大学 GWAS 课程：https://www.mv.helsinki.fi/home/mjxpirin/GWAS_course/

> PRS:
> > - Y 2020. Nature Protocols. [Tutorial: a guide to performing polygenic risk score analyses](https://www.nature.com/articles/s41596-020-0353-1)
> > - Y 2019. Lancet Respiratory Medicine. [Identification of risk loci and a polygenic risk score for lung cancer: a large-scale prospective cohort study in Chinese populations](pubmed.ncbi.nlm.nih.gov/31326317/)
> > - Y 2022. EHJ. [A polygenic risk score improves risk stratification of coronary artery disease: a large-scale prospective Chinese cohort study](pubmed.ncbi.nlm.nih.gov/35195259/)
> MR：
> > - Nature Reviews Methods Primers. [Mendelian randomization](https://www.nature.com/articles/s43586-021-00092-5)
> > - Y 2012. Lancet. [Plasma HDL cholesterol and risk of myocardial infarction: a mendelian randomisation study](https://pubmed.ncbi.nlm.nih.gov/22607825/)
> > - Y 2022. Diabetes Care. [Assessing the Causal Role of Sleep Traits on Glycated Hemoglobin: A Mendelian Randomization Study](https://pubmed.ncbi.nlm.nih.gov/35349659/)


<br/>
<br/>

# # 高阶理论与方法：
> - Y 2017. [Statistical methods to detect pleiotropy in human complex traits](https://pubmed.ncbi.nlm.nih.gov/29093210/)
> - Y 2019. [Meta-analysis and Mendelian randomization: A review](https://pubmed.ncbi.nlm.nih.gov/30861319/)
> - Y 2019. Nature Genetics. [A global overview of pleiotropy and genetic architecture in complex traits](https://www.nature.com/articles/s41588-019-0481-0)
> - Y 2021. [Genetic correlation and causal relationships between cardio-metabolic traits and lung function impairment](https://pubmed.ncbi.nlm.nih.gov/34154662/)
> - Y 2022. [Genetically predicted sex hormone levels and health outcomes: phenome-wide Mendelian randomization investigation](https://pubmed.ncbi.nlm.nih.gov/35218343/)

<br/>
<br/>

# # R分析和画图示例：
> - [The R Graph Gallery](https://r-graph-gallery.com/index.html)
> - [Doing and reporting your first mediation analysis in R](https://towardsdatascience.com/doing-and-reporting-your-first-mediation-analysis-in-r-2fe423b92171)
> - [Add P-values and Significance Levels to ggplots](https://www.r-bloggers.com/2017/06/add-p-values-and-significance-levels-to-ggplots/)
> - [Two-Way ANOVA Test in R](http://www.sthda.com/english/wiki/two-way-anova-test-in-r)
> - [ggfortify : Extension to ggplot2 to handle some popular packages](http://www.sthda.com/english/wiki/ggfortify-extension-to-ggplot2-to-handle-some-popular-packages-r-software-and-data-visualization)
> - [An Intro to Phylogenetic Tree Construction in R](https://fuzzyatelin.github.io/bioanth-stats/module-24/module-24.html)
> - [Top 100 R resources on COVID-19 Coronavirus](https://statsandr.com/blog/top-r-resources-on-covid-19-coronavirus/)
> - 以及 modelSummary, forplo，sankey diagram, CellChat, ComplexHeatmap，等等
