dir=/mnt/d/data/ukb/phe

# cat $dir/rap/raw/pheno.tsv | sed -e 's/^\t/NA\t/; s/\t\t/\tNA\t/g; s/\t\t/\tNA\t/g; s/\t$/\tNA/' > pheno.tab
head -1 raw/pheno.tab | tr '\t' '\n' > pheno.id
	awk '{print $1,$1}' pheno.id | sed 's/_\w\+//1' > pheno.id.2col
	sed 's/_.*//' pheno.id | sort | uniq > pheno.id.uniq; wc -l pheno.id.uniq


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 提取主要表型pheno和PCA数据
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
awk 'BEGIN{print "eid"}{print $1}' $dir/common/ukb.vip.dat > vip.id; sort vip.id | uniq -d
	fgrep -wf vip.id pheno.id.uniq | wc -l
	fgrep -vwf pheno.id.uniq vip.id | head # 显示UKB里面没有的 
fgrep -nwf vip.id pheno.id.2col | awk '$2 !~/_/ || $2 ~/_i0/' | awk -F ":" '{print $1}'  | tr '\n' ',' | sed 's/,$//' | xargs -n1 -I % cut -f % raw/pheno.tab > vip.tab
# 一个 data field，多个column的数据
Arr1=(p87 p20002 p20003 p22009 p40002 p41270 p41271 p41280 p41281)
Arr2=(srdTime srd med pca death_icd10s icd10 icd9 icd10Date icd9Date)
for i in ${!Arr1[@]}; do
	field=${Arr1[$i]}; name=${Arr2[$i]}; echo RUN: $field $name
	fgrep -nw $field pheno.id.2col | awk -F ":" 'BEGIN{print 1}{print $1}' | tr '\n' ',' | sed 's/,$//' | xargs -n1 -I % cut -f % raw/pheno.tab > $name.tab
done