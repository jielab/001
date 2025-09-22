#%% 
import pandas as pd
import tqdm
import numpy as np

dir0 = "/mnt/d"
labels_file  = f"{dir0}/data/ukb/phe/delphi/labels.csv"
ukb_field_to_icd10_map_file = f"{dir0}/data/ukb/phe/delphi/icd10_codes_mod.tsv"
ubk_basket_tab_file = f"{dir0}/data/ukb/phe/rap/raw/pheno.tsv.gz" 
train_proportion = 0.8 # proportion of full data set to use for training (the rest will be used for validation)
output_prefix = 'real'


#%% 
# Read icd10 mapping file and defined index label link
icdict ={}
icdcodes = []
with open(ukb_field_to_icd10_map_file,'r') as f:
    for l in f:
        lvals=l.strip().split()
        icdict[lvals[0]]=lvals[5]
        icdcodes.append(lvals[5])

i = -1
label_dict = {}
with open(labels_file,'r') as f:
    for l in f:
        label_dict[l.strip().split(' ')[0]]=i
        i += 1

# hard coded sex and dob
icdict['p31'] = "sex"
icdict['p34'] = "YEAR"
icdict['p52'] = "MONTH"
icdict['p40000_i0'] = "Death" # ther is also p40000_i1

# cancer fields
for j in range(21):
    icdict['p40005_i'+str(j)] = "cancer_date_"+str(j)
    icdict['p40006_i'+str(j)] = "cancer_type_"+str(j)

# cancer hes fields 
#for j in range(213):
#    icdict['p41270_i'+str(j)] = "hicd_"+str(j)
#    icdict['p41280_i.'+str(j)] = "hicd_date_"+str(j)

icdict['p53_i0'] = "assessment_date"
icdict['p21001_i0']="BMI"
icdict['p1239_i0']="smoking"
icdict['p1558_i0']="alcohol"

len_icd = len(icdcodes)
#icdcodes.extend(['Death','assessment_date']+['cancer_date_'+str(j) for j in range(17)]+['hicd_date_'+str(j) for j in range(213)])
icdcodes.extend(['Death','assessment_date']+['cancer_date_'+str(j) for j in range(17)])


#%% 
# Read ukb basket file in chunks, select icd10 code occurance and dates, format for delphi
data_list = []
ukb_iterator = pd.read_csv(ubk_basket_tab_file, sep="\t", chunksize=1000, index_col=0)
first = next(ukb_iterator)                 # grab first chunk (DataFrame)
print(first.head(10))                  # head
print("cols:", first.columns.tolist()) # columns


#%% 
for _, dd in tqdm.tqdm(enumerate(ukb_iterator)):
    dd = dd.rename(columns=icdict)
    dd.dropna(subset=['sex'], inplace=True)
    dd['sex'] += 1
    dd = dd[[col for col in dd.columns if not col.startswith('p')]]
    dd['dob'] =  pd.to_datetime(dd[['YEAR', 'MONTH']].assign(DAY=1))
    dd[icdcodes] = dd[icdcodes].apply(pd.to_datetime, format="%Y-%m-%d")
    dd[icdcodes]=dd[icdcodes].sub(dd['dob'], axis=0)
    dd[icdcodes]=dd[icdcodes].apply(lambda x : x.dt.days)

    for col in icdcodes[:len_icd+1]:
        X = dd[col].dropna().reset_index().to_numpy().astype(int)
        data_list.append(np.hstack((X,label_dict[col]*np.ones([X.shape[0],1],X.dtype))))
    
    X = dd['sex'].reset_index().to_numpy().astype(int)
    data_list.append(np.c_[X[:,0],np.zeros(X.shape[0]),X[:,1]].astype(int))
     
    for j in range(21):
        dd_cancer = dd[['cancer_date_'+str(j),'cancer_type_'+str(j)]].dropna().reset_index()
        if not dd_cancer.empty:
            dd_cancer['cancer'] = dd_cancer['cancer_type_'+str(j)].str.slice(0,3)
            dd_cancer['cancer_label'] = dd_cancer["cancer"].map(label_dict)
            data_list.append(dd_cancer[['eid','cancer_date_'+str(j),'cancer_label']].dropna().astype(int).to_numpy())

    #for j in range(213):
    #    dd_hicd = dd[['hicd_date_'+str(j),'hicd_'+str(j)]].dropna().reset_index()
    #    if not dd_hicd.empty:
    #        dd_hicd['hicd'] = dd_hicd['hicd_'+str(j)].str.slice(0,3)
    #        dd_hicd['hicd_label'] = dd_hicd["hicd"].map(label_dict)
    #        data_list.append(dd_hicd[['eid','hicd_date_'+str(j),'hicd_label']].dropna().astype(int).to_numpy())
        
    dd_bmi = dd[['assessment_date','BMI']].dropna().reset_index()
    dd_bmi['bmi_status'] = np.where(dd_bmi['BMI']>28,5,np.where(dd_bmi.BMI>22,4,3))
    data_list.append(dd_bmi[['eid','assessment_date','bmi_status']].astype(int).to_numpy())
    
    dd_sm = dd[['assessment_date','smoking']].dropna().reset_index()
    dd_sm = dd_sm[dd_sm['smoking']!=-3]
    dd_sm['smoking_status'] = np.where(dd_sm['smoking']==1,8,np.where(dd_sm.smoking==2,7,6))
    data_list.append(dd_sm[['eid','assessment_date','smoking_status']].astype(int).to_numpy())
    
    dd_al = dd[['assessment_date','alcohol']].dropna().reset_index()
    dd_al = dd_al[dd_al['alcohol']!=-3]
    dd_al['alcohol_status'] = np.where(dd_al['alcohol']==1,11,np.where(dd_al.alcohol < 4,10,9))
    data_list.append(dd_al[['eid','assessment_date','alcohol_status']].astype(int).to_numpy())

  
#%% 
# reformat, split train and val and output to delphi format
data= np.vstack(data_list)
data = data[np.lexsort((data[:,1], data[:,2]==data[:,2].max(), data[:,0]))]
data = data[data[:,1]>=0]
data = pd.DataFrame(data).drop_duplicates([0,2]).values
data = data.astype(np.uint32)
data.tofile(output_prefix + '.bin')
ids = list(set(data[:,0]))
ids = ids.sort()
train_val_split = data[:,0] <= ids[int(len(ids)*train_proportion)]
data[train_val_split].tofile(output_prefix + '_train.bin')
data[~train_val_split].tofile(output_prefix + '_val.bin')