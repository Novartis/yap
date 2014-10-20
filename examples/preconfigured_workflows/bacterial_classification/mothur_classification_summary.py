#!/usr/bin/env python

from sys import *
from operator import add
import glob
import os 
if len(argv)==5:
        input_path = argv[2]
        output_path = argv[4]
else:
        print "Please issue the command in the following way:\n\n\tpython mothur_classification_summary -i <input_path> -o <output_path>"
	exit()
## Different file extensions/suffixes to be searched for in the input folder
count_type=['wang.pick.tax.summary','wang.tax.summary']

output_dir=os.path.split(output_path)[0]
if output_dir == '':
        output_dir = os.getcwd()

def merge_count(input_path,output_path,pattern):
	file_list=glob.glob(input_path+'/*'+pattern)
	if len(file_list) != 0:
		merge_dict = {}
		header=''
		for file in file_list:
			file_arr = []
			fh = open(file,'rb')
			file_arr = fh.readlines()
			for i in range(0,len(file_arr)):	
				temp_arr =[]
				temp_arr = file_arr[i].strip().split()
				key=''
				value=[]
				if len(temp_arr) != 5:
					print '''Input format error: Looking for mothur taxonomy columns
                        		[taxlevel,rankID,taxon,daughterlevels,total]'''
					exit()
				else:
					if i == 0:
						header = file_arr[i]
						pass
					else:	
						key = temp_arr[0] + "\t" + temp_arr[1] + "\t" + temp_arr[2]
						#key = temp_arr[0] + "\t" + temp_arr[1] + "\t" + temp_arr[2] \
							#+ "\t" + temp_arr[3]
						#value = [int(temp_arr[3]),int(temp_arr[4])]
						value = [int(temp_arr[4])]
		        			if key in merge_dict.keys():
							merge_dict[key] = map(add, merge_dict[key],value)
						else:
							merge_dict[key] = value 
		with open(output_path+ "-" + pattern,'w') as output_file:
			output_file.writelines(header)
			for key in merge_dict.keys():
				#output_file.writelines(key+'\t'+str(merge_dict[key][0]) + "\t" + str(merge_dict[key][1])+'\n')
				output_file.writelines(key+'\t'+str(merge_dict[key][0])+'\n')
			output_file.close()
			#sort output based on column 2 (rankID)		
                        output_filename = output_path+ "-" + pattern
			sort_cmd =  "sort -k 2 " + output_filename + " -o " + output_filename 
			os.system(sort_cmd)
if os.path.exists(input_path):
	if os.path.exists(output_dir):
		for i in count_type:
			merge_count(input_path,output_path,i)
	else:
		print "Error: In Merge counts - Output path does not exist"
else:
	print "Error: In Merge counts - Input path does not exist"

