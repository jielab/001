# coding=utf-8
import argparse
from collections import OrderedDict
	
sep = { 'COMMA' : ',' , "SPACE" : ' ','TAB' : '\t'}
	
# Parse all of our file
usage = 'python mergeFile.py -i "file1,SEP1,col1 file2,SEP,col2[ file3,SEP3,col3]" -o outFile'
parser = argparse.ArgumentParser(description=usage)
parser.add_argument("-i", dest = "input", help='file format')
parser.add_argument('-o', dest = 'out', help = 'output file name')
args = parser.parse_args()
	
try:
	# input parser
	input = (args.input).split(' ')
	input_list = []
	if len(input) < 2:
		print ('illegal -i option, please check again! No file processed!')
		exit()
	else:
		try:
			for file_par in input:
				file_par = file_par.split(',')
				line = []
				line.append(file_par[0]) # file name
				line.append(sep[file_par[1]]) # seperator
				line.append(int(file_par[2])) # col num
				input_list.append(line)
		except:
			print ('illegal -i option, please check again! No file processed!')
			exit()
	#output name
	out = args.out
except:
	print ("Wrong usage, run 'python mergeFile.py -h'")
	exit()


	
def build_hash(input):
	item = OrderedDict()
	c = 0
	line_count = 0;
	with open(input[0]) as fp2:
		for line in fp2.readlines():
			line = line.strip()
			line = line.strip(input[1])
			line_count = line_count + 1
			line_list = line.split(input[1])
			if c == 0:
				c = len(line_list)
			else:
				def test(list):
					# testing blank line
					for a in list:
						if len(a) == 0:
							return False
					return True
				if not test(line_list):
					print(line_list)
					return False,line_count
			tmp = ' '.join(line_list)
			item[line_list[input[2]]] = tmp
	cnt = ''
	for i in range(c - 1):
		cnt = cnt + 'NA '
	cnt  = cnt + 'NA'
	return item, cnt
	

hash_file1,cnt =  build_hash(input_list[0])
#print(hash_file1)
#print(cnt)
if (isinstance(hash_file1, bool)):
		print('file %s has inconsistent field numbers in line %d' % (input_list[0][0], cnt))
		exit()
del hash_file1	
	
hash_file_list = []	
for i in range(len(input_list) - 1):
	hash_file,cont  =  build_hash(input_list[i + 1])
	if isinstance(hash_file, bool):
		print('file %s has inconsistent field numbers in line %d' % (input_list[i + 1][0], cont))
		exit()
	hash_file_list.append((hash_file,cont))

	
with open(input_list[0][0]) as fp_in:
	sep = input_list[0][1]
	col = input_list[0][2]
	with open(out, 'w') as fp_out:
		for line in fp_in.readlines():
			line =  line.strip()
			line = line.strip(sep)
			line_list = line.split(sep)
			for hash_file in hash_file_list:
				if line_list[col] in hash_file[0]:
					line =  line + ' ' + hash_file[0][line_list[col]]
				else:
					line = line + ' ' + hash_file[1]
			line = line + '\n'
			fp_out.write(line)
