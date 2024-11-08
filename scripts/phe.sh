appID=66137 # 整个项目的ID，需要写到文章的致谢部分。

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. ukbunpack 打开包裹
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 上面两个 ukbunpack 命令可分开跑，跑完后分别跑下面的 unkbconvt 命令。
runID=48807 # Basket 2014862
    md5sum=eeff5373a57254ec5a0487e4ea2a943d
    key=e64a3f781*c2c8*00ed4ecfbe9a4f*2c0b5115aeb77e1cb1eba01f3d*e3463af
runID=50136 # Basket 2016432
	md5sum=981b47f85c6b2fb849320c7a3559ba23
	key=fd5b09bb6f6db*e9e0774e91f81*af1a2d9e*aea3ddad*2a2ff6ea60f2f23458
runID=670287 # Basket 4018893
	md5sum=6cfe4d40d12fdb1b8564dc67ce91ad20
	key=a93*00fe54833*70898b8e87f96646fdf52f509de2b87248ed1ea92b560c2ba1
echo -e "$appID\n$key" > .ukb.key
	ukbunpack ukb$runID.enc


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. ukbconv 提取所需数据，并转换为所需格式
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 不用 ukbconv txt，而用 ukbconv r，因为前者生成的output文件第一列有多余的tab，最后一列有多余的^M
awk '{print $1}' /mnt/d/data/ukb/phe/common/ukb.vip.dat > ukb.vip.id; sort ukb.vip.id | uniq -d
enc=ukb50136.enc_ukb
	ukbconv $enc docs
	ukbconv $enc r -iukb.vip.id -ovip # 提取 ukb.vip.id 文件中罗列的表型，可以添加 –eencoding.dat 的 option 
	ukbconv $enc r -s22182 -ohla
	ukbconv $enc r -s22009 -opc
enc=ukb670287.enc_ukb
	ukbconv $enc docs
	ukbconv $enc r -iukb.vip.id -ovip
	ukbconv $enc r -s20002 -osrd
	ukbconv $enc r -s87 -osrdTime
	ukbconv $enc r -s20003 -omed
	ukbconv $enc r -s40002 -odeath_icd10s # Contributory (secondary) causes of death: ICD10
	ukbconv $enc r -s41270 -oicd10 # Diagnoses - ICD10
	ukbconv $enc r -s41280 -oicd10Date # Date of first in-patient diagnosis - ICD10
	ukbconv $enc r -s41271 -oicd9 # Diagnoses - ICD9
	ukbconv $enc r -s41281 -oicd9Date # Date of first in-patient diagnosis - ICD9
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
