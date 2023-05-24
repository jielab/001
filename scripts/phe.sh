appID=66137
runID=670287 #50136
enc=ukb$runID.enc_ukb
key=fd5b09bb6f6db1e9e0774e91f812af1a2d9e3aea3ddad42a2ff6ea60f2f23458
echo -e "$appID\n$key" > $runID.key
awk '{print $2}' /mnt/d/data/ukb/phe/common/ukb.vip.dat > ukb.vip.id


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# main PHE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# "ukbconv txt" creats leading tabs, windows ukbconv creates ^M
cd /mnt/d/data/ukb/pheno/raw-$runID/
ukbmd5 ukb$runID.enc # make sure it is 981b47f85c6b2fb849320c7a3559ba23
ukbunpack ukb$runID.enc $runID.key
ukbconv $enc docs # also creates fields.ukb
ukbconv $enc r -iukb.vip.id -ovip # –eencoding.dat
ukbconv $enc r -s22182 -ohla
ukbconv $enc r -s41270 -oicd
ukbconv $enc r -s41280 -oicdDate
cat icd.tab | sed 's/\tNA//g; s/"//g; s/\t/,/2g; ' |
	awk '{if(NR==1) print "eid icd10"; else if (NF==1) print $1 "\tNA"; else print $0}' > icd.2cols
awk '$1>=130000 && $1 <=132604 && $1%2==0' fields.ukb > fod.id 
ukbconv $enc r -ifod.id -ofod


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# bulk data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ukbconv $enc bulk -s22002 -oCEL 
ukbconv $enc bulk -s23151 -oexome_vcf # also 23152 
ukbconv $enc bulk -s23181 -owgs_cram_bgi # also 23182
ukbconv $enc bulk -s23183 -owgs_cram_broad # also 23184
