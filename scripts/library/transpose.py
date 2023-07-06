#!/usr/bin/python

## transpose.py 200 ' ' . "F12a.WGS.gz" transposed.gz

import sys,gzip,os,tempfile,time

if len(sys.argv) < 6 :
	print 'Usage: %s buf_limit_mb seperator work_dir in_file out_file'%(sys.argv[0])
	sys.exit(1)
buf_lim = int(sys.argv[1])*2**20 #in Bytes
sep=sys.argv[2]
work_dir=sys.argv[3]
in_file_name=sys.argv[4]
out_file_name=sys.argv[5]

fin = gzip.open(in_file_name)
# fout_tmp = gzip.GzipFile(fileobj=tempfile.TemporaryFile(dir=work_dir))
# fout_tmp.writelines(c+'\n' for c in fin.readline()[:-1].split(sep))
# fin_tmp = gzip.GzipFile(fileobj=fout_tmp.fileobj,mode='rb')
# fout_tmp.close()
# fin_tmp.fileobj.seek(0)
fout_tmp = tempfile.TemporaryFile(dir=work_dir)
fout_tmp.writelines(c+'\n' for c in fin.readline()[:-1].split(sep))
fout_tmp.flush()
fin_tmp = fout_tmp
fin_tmp.seek(0)

in_file_size = os.path.getsize(in_file_name)
t = time.time()
progress = 0

bytes_read = 0
matrix_buffer = []
for l in fin :
	if bytes_read + len(l)  >= buf_lim :
# 		fout_tmp = gzip.GzipFile(fileobj=tempfile.TemporaryFile(dir=work_dir))
# 		fout_tmp.writelines(sep.join(l2[:-1].split(sep) + [l3[i] for l3 in matrix_buffer])+'\n' for i,l2 in enumerate(fin_tmp))
# 		fin_tmp = gzip.GzipFile(fileobj=fout_tmp.fileobj,mode='rb')
# 		fout_tmp.close()
# 		fin_tmp.fileobj.seek(0)
		fout_tmp = tempfile.TemporaryFile(dir=work_dir)
		fout_tmp.writelines(sep.join(l2[:-1].split(sep) + [l3[i] for l3 in matrix_buffer])+'\n' for i,l2 in enumerate(fin_tmp))
		fout_tmp.flush()
		fin_tmp = fout_tmp
		fin_tmp.seek(0)
		
		bytes_read = 0
		matrix_buffer = []
	
	matrix_buffer.append(l[:-1].split(sep))
	bytes_read += len(l)
	p = int(fin.fileobj.tell()/float(in_file_size)*100)
	if p > progress :
		sys.stdout.write('\r%d%% Time elapsed:%ds'%(p,time.time()-t))
		sys.stdout.flush()
		progress = p
		
fout = gzip.open(out_file_name,'wb')
fout.writelines(sep.join(l1[:-1].split(sep) + [l2[i] for l2 in matrix_buffer])+'\n' for i,l1 in enumerate(fin_tmp))
sys.stdout.write('\nJob done in %ds\n'%(time.time()-t))
