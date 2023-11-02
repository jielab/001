appID=66137 # 整个项目的ID，需要写到文章的致谢部分。

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. ukbunpack 打开包裹
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 上面两个 ukbunpack 命令可分开跑，跑完后分别跑下面的 unkbconvt 命令。
runID=50136 # 某年某日申请到的数据ID
	key=fd5b09bb6f6db1e9e0774e91f812af1a2d9e3aea3ddad42a2ff6ea60f2f23458
	echo -e "$appID\n$key" > $runID.key
	ukbunpack ukb$runID.enc $runID.key
runID=670287 # 某年某日申请到的新数据ID
	key=a93100fe54833270898b8e87f96646fdf52f509de2b87248ed1ea92b560c2ba1
	echo -e "$appID\n$key" > $runID.key
	ukbunpack ukb$runID.enc $runID.key


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. ukbconv 提取所需数据，并转换为所需格式
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 不用 ukbconv txt，而用 ukbconv r，因为前者生成的output文件第一列有多余的tab，最后一列有多余的^M
awk '{print $2}' /mnt/d/data/ukb/phe/common/ukb.vip.dat > ukb.vip.id
enc=ukbXXX.enc_ukb # 根据上面第1步生成的文件，填入XXX
ukbconv $enc docs # 生成 .html 和 fields.ukb 概括性文件
ukbconv $enc r -iukb.vip.id -ovip # 提取 ukb.vip.id 文件中罗列的表型，可以添加 –eencoding.dat 的 option 
ukbconv $enc r -s22182 -ohla # 提取 HLA
ukbconv $enc r -s41270 -oicd # 提取 ICD-10
	ukbconv $enc r -s41280 -oicdDate # 提取 ICD-10 date
	cat icd.tab | sed 's/\tNA//g; s/"//g; s/\t/,/2g; ' | awk '{if(NR==1) print "eid icd10"; else if (NF==1) print $1 "\tNA"; else print $0}' > icd.2cols
awk '$1>=130000 && $1 <=132604 && $1%2==0' fields.ukb > fod.id # 提取 fod(first occurence date) 数据
ukbconv $enc r -ifod.id -ofod


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. 提取 bulk 数据
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ukbconv $enc bulk -s22002 -oCEL 
ukbconv $enc bulk -s23151 -oexome_vcf # also 23152 
ukbconv $enc bulk -s23181 -owgs_cram_bgi # also 23182
ukbconv $enc bulk -s23183 -owgs_cram_broad # also 23184
