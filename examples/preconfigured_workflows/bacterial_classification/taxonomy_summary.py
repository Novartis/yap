#!/usr/bin/env python

from sys import *
from operator import add
import glob
import os 
if len(argv)==5:
        input_path = argv[2]
        output_path = argv[4]
else:
        print "Please issue the command in the following way:\n\n\tpython taxonomy_summary -i <input_path> -o <output_path>"
	exit()
## Different file extensions/suffixes to be searched for in the input folder
count_type=['wang.pick.taxonomy','wang.taxonomy']

output_dir=os.path.split(output_path)[0]
if output_dir == '':
        output_dir = os.getcwd()

def merge_files(input_path,output_path,pattern):
	file_list=glob.glob(input_path+'/*'+pattern)
	if len(file_list) != 0:
		merge_dict = {}
		header = ''
		input_str = ''
		for file in file_list:
			input_str += file + " "
		os.system("cat " + input_str + " >" + output_path+ "-" + pattern)
if os.path.exists(input_path):
	if os.path.exists(output_dir):
		for i in count_type:
			merge_files(input_path,output_path,i)
	else:
		print "Error: In Merge counts - Output path does not exist"
else:
	print "Error: In Merge counts - Input path does not exist"

