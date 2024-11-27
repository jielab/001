dir=/mnt/d/data/ukb/phe

# cat $dir/rap/raw/pheno.tsv | sed -e 's/^\t/NA\t/; s/\t\t/\tNA\t/g; s/\t\t/\tNA\t/g; s/\t$/\tNA/' > pheno.tab
zcat raw/pheno.tab | head -1 | tr '\t' '\n' > pheno.id
	awk '{print $1,$1}' pheno.id | sed 's/_\w\+//1' > pheno.id.2col
	sed 's/_.*//' pheno.id | sort | uniq > pheno.id.uniq; wc -l pheno.id.uniq


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 提取主要表型pheno和PCA数据
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat=dat # dat met 等
	awk 'BEGIN{print "eid"}{print $1}' $dir/common/ukb.vip.$dat > $dat.id; sort $dat.id | uniq -d
	fgrep -wf $dat.id pheno.id.uniq | wc -l
	fgrep -vwf pheno.id.uniq $dat.id | head # 显示UKB里面没有的 
	fgrep -nwf $dat.id pheno.id.2col | awk '$2 !~/_/ || $2 ~/_i0/' | awk -F ":" '{print $1}'  | tr '\n' ',' | sed 's/,$//' | xargs -n1 -I % cut -f % raw/pheno.tab > $dat.tab
# 一个 data field，多个column的数据
Arr1=(p87 p20002 p20003 p22009 p22182 p40002 p41270 p41271 p41280 p41281)
Arr2=(srdTime srd med pca hla death_icd10s icd10 icd9 icd10Date icd9Date)
for i in ${!Arr1[@]}; do
	field=${Arr1[$i]}; name=${Arr2[$i]}; echo RUN: $field $name
	fgrep -nw $field pheno.id.2col | awk -F ":" 'BEGIN{print 1}{print $1}' | tr '\n' ',' | sed 's/,$//' | xargs -n1 -I % cut -f % raw/pheno.tab > $name.tab
done
sed -i '1d; s/|/\t/g' icd10Date # 去掉第一行


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🈲 ukbunpack 打开包裹
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
appID=66137
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
