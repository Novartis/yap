#!/usr/bin/env python
"""
Copyright 2014 Novartis Institutes for Biomedical Research

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import sys
import re
import glob
import os
import math
import shlex
import random
import numpy
import itertools
from multiprocessing import *
from subprocess import PIPE, Popen
from pytee import *
from string import Template
import yap_cmd_checks
import yap_file_io 
import yap_log
import yap_workflow_dict

## global variables
global basename_warnings
global missing_path_errors

## global variable intialization
basename_warnings = []
missing_path_errors = []

def create_dictionary(config_file):
    """ configuration file reader /
	creates dictionary of parameters """
    try:
        result = numpy.loadtxt(config_file, dtype = "a100",comments = '#',delimiter = ':=')
	workflow_config_dict = {}
    	for i in range(0,len(result)):
        	k = re.search('(.*\")(.*)(\".*)',result[i][0])
        	v = re.search('(.*\")(.*)(\".*)',result[i][1])
        	key = k.group(2).strip(" ")
        	val = v.group(2).strip(" ")
        	key = key.strip("\t")
        	val = val.strip("\t")
        	workflow_config_dict.setdefault(key,val)
    except IOError as (errno ,strerror):
	print "Error:  while Opening the file : " , config_file 
	print "	       I/O error({0}): {1}".format(errno,strerror)
	exit()
    except :
	print "Error : encountered format error while creating the dictionary for file : ",config_file
        print """	Please check the format the configuration file
	enclose variable and corresponding values in double quotes("...") 
       	use (:=) in variable value assignment (eg "variable" := "value")
        to add comments,start the line with symbol(#)
        please refer the documentation for the configuration file formats
              """
        exit()
    return workflow_config_dict

def split_array(totnum,nprocs,lpp,myrank):
    """ Returns the indices of the chunk(of the input file) to be split based on the 
        rank of the processor."""
    totnum = totnum/lpp
    n2=int(math.floor(totnum/nprocs))
    remain=int(math.fmod(totnum,nprocs))
    if (myrank < remain) :
        n1=n2+1
        ib1=(myrank)*n1
        ie1=ib1+(n1-1)
    else :
        n1=n2
        ib1=(myrank)*n1+remain
        ie1=ib1+n1-1
    ib=ib1*lpp
    ie=((ie1+1)*lpp)-1
    ret=[]
    ret.append(ib)
    ret.append(ie)
    return ret

def get_file_size(inputfile_name):
    """ Returns the size of a given file using file seek"""
    try:
	file_size_on_disk=os.path.getsize(inputfile_name)
    	cs=file_size_on_disk
    except Exception as e:
	print e,"Error in finding the file size"
    input_file  = yap_file_io.create_openfile_handler(inputfile_name)
    file_size=0
    input_file.seek(0)
    file_pos=0
    try:
            file_pos = 0
            end_file = 'False'
            file_pos1 = 0
            while end_file == 'False':
                file_pos=file_pos+cs
                input_file.seek(file_pos)
                read_buffer = input_file.read(1)
                length_read_buffer = len(read_buffer)
                if length_read_buffer == 0 :
                        file_size=input_file.tell()
                        end_file = 'True'
    finally:
           return file_size

def get_file_split_position_symbolbased(inputfile_name,nchunks,file_size ,input_file_format):
    """ Returns the split file indices based on the input file format and specified
        chunk size."""
    input_file = create_openfile_handler(inputfile_name)
    chunk_size = file_size / nchunks
    n_parts = 1
    chunk_size_arr = []
    search_symbol = ''
    j = 0
    input_file.seek(0)
    check_first_element = input_file.read(1)
    # file format checking 
    if input_file_format  == 'fastq':
        search_symbol = '@'
        if check_first_element != search_symbol:
                print " Error : Please specify the correct input file format in workflow configuration file"
                print " 	Input Data  is not in the given " , input_file_format , "format"
                exit()
    if input_file_format == 'fasta':
        search_symbol = '>'
        if check_first_element != search_symbol:
                print " Error : Please specify the correct input file format in workflow configuration file"
                print " 	Input Data  is not in the given " , input_file_format , "format"
                exit()
    if input_file_format == 'qseq' or input_file_format == 'tab':
        search_symbol = '\n'
    try:
            end_pos = 0
            i = 0
            buf = 0
            end_pos_arr = []
	    file_end = 'False'
	    file_pos = 0
            input_file.seek(0)
            for i in range(1,nchunks+1):
                hit_pos = -1
		input_file.seek(chunk_size,1)
		current_seek = input_file.tell()
                if i == nchunks or current_seek >= file_size:
                        end_pos = file_size
			end_pos_arr.append(end_pos)
                        break
                while hit_pos == -1:
                        buf += 1024
                        read_buffer = input_file.read(buf)
			length_read_buffer = len(read_buffer)
                        if length_read_buffer == buf:
                                hit_pos = read_buffer.rfind(search_symbol)
			else:
				hit_pos = length_read_buffer
		if search_symbol == "\n":
			end_pos = current_seek + (hit_pos+1)
		else:
			end_pos = current_seek + hit_pos
                end_pos_arr.append(end_pos)
		if end_pos >= file_size:
			break
            begin_end_index = numpy.zeros([len(end_pos_arr),3],int)
            for kk in range(0,len(end_pos_arr)):
                if kk == 0:
                        begin_end_index[kk][0] = 0
                else:
                        begin_end_index[kk][0] = end_pos_arr[kk-1]
                begin_end_index[kk][1] = end_pos_arr[kk]
                begin_end_index[kk][2] = kk
            file_split_info = []
            file_split_info.append(inputfile_name)
            file_split_info.append(begin_end_index)
    finally:
            input_file.close()
    return file_split_info

def convert_format(seqs_str):
    	""" Converts formats between qseq, fastq, fasta & tab if 
        the specified input and output formats are different. """
	seqs_arr = seqs_str.splitlines(1)
        seqs_arr_len = len(seqs_arr)
        if seqs_arr_len <= 0:
                print "Empty data : no sequences read or empty input file"
                exit()
        input_file_format = yap_workflow_dict.input_file_format
	output_file_format = yap_workflow_dict.preprocess_output_file_format
	format_seqs_arr = []
	if output_file_format == 'qseq':
        	print "please specify the output file format as fastq or fasta). output file format given := ", output_file_format
                exit()
	if input_file_format == 'qseq':
                try:
                        for i in range (0,len(seqs_arr)):
                                record = seqs_arr[i].strip("\n").split('\t')
				machine_name = record[0]
			        run_number = record[1]
            			lane_number = record[2]
            			tile = record[3]
            			x = record[4]
            			y = record[5]
				read = record[7]
                                sequence = record[8]
                                quality = record[9]
				pass_qc_msg = record[10]
				if pass_qc_msg == str(1):
                                	seq_id = '%s_%s:%s:%s:%s:%s#%s/%s' % (machine_name,run_number,lane_number,tile,x,y,read,pass_qc_msg)
					if output_file_format == 'fasta':
						format_seqs_arr.append(">" + seq_id + "\n" )
                                        	format_seqs_arr.append(sequence + "\n")
					elif output_file_format == 'tab':
						format_seqs_arr.append(seq_id + "\t" + sequence + "\n" )
					else:
						format_seqs_arr.append("@" + seq_id + "\n" )
						format_seqs_arr.append(sequence+ "\n")
						format_seqs_arr.append("+" + "\n")
						format_seqs_arr.append(quality + "\n")
                except Exception as e:
                        print " Error : Please specify the correct input file format in workflow configuration file"
                        print " 	Input Data  is not in the given " , input_file_format , "format"
                        exit()

        if input_file_format == 'fastq':
                if seqs_arr[0][0] != "@":
                        print " Error : Please specify the correct input file format in workflow configuration file"
                        print " 	Input Data  is not in the given " , input_file_format , "format"
                        exit()
                for i in range(0,seqs_arr_len,4):
				seq_id = seqs_arr[i].strip("\n")
                                sequence = seqs_arr[i+1].strip("\n")
				desc = seqs_arr[i+2].strip("\n")
				quality = seqs_arr[i+3].strip("\n")
                                if output_file_format == 'fasta':
                                        format_seqs_arr.append(">" + seq_id.lstrip('@') + "\n")
                                        format_seqs_arr.append(sequence + "\n")
                                elif output_file_format == 'tab':
                                        format_seqs_arr.append(seq_id.lstrip('@') + "\t" + sequence + "\n" )
				else:
					format_seqs_arr.append(seq_id  + "\n")
					format_seqs_arr.append(sequence + "\n")
					format_seqs_arr.append(desc + "\n")
					format_seqs_arr.append(quality + "\n")
	if input_file_format == 'fasta':
                if seqs_arr[0][0] != ">":
                        print " Error : Please specify the correct input file format in workflow configuration file"
                        print " 	Input Data  is not in the given " , input_file_format , "format"
                        exit()
                for  i in range(0,seqs_arr_len,2):
			seq_id = seqs_arr[i].strip("\n")
                	sequence = seqs_arr[i+1].strip("\n")
			if output_file_format  == 'tab' :
				format_seqs_arr.append(seq_id.lstrip('>') + "\t" + sequence + "\n" )
			elif output_file_format == 'fasta':
				format_seqs_arr.append(seq_id + "\n" ) 
				format_seqs_arr.append(sequence + "\n" )			
			else:
				format_seqs_arr.append("@" + seq_id.lstrip('>') + "\n")
                        	format_seqs_arr.append(sequence + "\n")
                        	format_seqs_arr.append("+" + "\n")
                        	format_seqs_arr.append(" "+ "\n")
	if input_file_format == 'tab':
	 	try:
                        for line in seqs_arr:
                                record = line.strip('\n').split('\t')
				if len(record) != 2:
                        		print " Error : Please specify the correct input file format in workflow configuration file"
                        		print " 	Input Data  is not in the given " , input_file_format , "format"
                			exit()
				seq_id = record[0] 
                                sequence = record[1] 
		                quality = " "
                        	if output_file_format  == 'tab' :
                                	format_seqs_arr.append(seq_id + "\t" + sequence + "\n")
                        	elif output_file_format == 'fasta':               
                                	format_seqs_arr.append(">"+seq_id + "\n")
                                	format_seqs_arr.append(sequence + "\n" )
				else:
	        	                format_seqs_arr.append("@" + seq_id + "\n")
        	                        format_seqs_arr.append(sequence + "\n")
                                	format_seqs_arr.append("+" + "\n")
                                	format_seqs_arr.append(" "+ "\n")
                except Exception as e:
                        print " Error : Please specify the correct input file format in workflow configuration file"
                        print " 	Input Data  is not in the given " , input_file_format , "format"
                        exit()
	return ''.join(format_seqs_arr)

def format_sequences(seqs_arr,output_dict):
    """Print the Sequence array into the corresponding format given in the Main Config File"""
    output_format = yap_workflow_dict.preprocess_output_file_format
    out_arr= []
    if output_format == "fasta":
        for i in range(0,len(seqs_arr)):
                j = i
                out_arr.append('>' + seqs_arr[j][0] + "\n" + seqs_arr[j][1] + "\n")
    if output_format == "fastq":
        for i in range(0,len(seqs_arr)):
                j = i
                out_arr.append('@'+ seqs_arr[j][0]+ "\n" + seqs_arr[j][1] + "\n" + '+' + "\n" + seqs_arr[j][2] + "\n")
    if output_format == "tab":
        for i in range(0,len(seqs_arr)):
                j = i
                out_arr.append(seqs_arr[j][0] + "\t" + seqs_arr[j][1] + "\n")
    return out_arr

def qc_basecount(seqs_str, workflow_prov):
    """ Returns the base count per read location """
    bases_dict = {'A':0,'a':0,'C':1,'c':1,'T':2,'t':2,'G':3,'g':3,'N':4,'n':4}
    output_file_format = yap_workflow_dict.preprocess_output_file_format
    max_read_length = int(yap_workflow_dict.max_read_length)
    if output_file_format == "fasta":
        loop_increment = 2 
    elif output_file_format == "tab":
        loop_increment = 1
    else:
        loop_increment = 4
    alphabet_size = 5
    base_count_per_read_location =numpy.zeros((max_read_length,alphabet_size),dtype = numpy.int)
    if seqs_str != '':
        seqs_arr = seqs_str.splitlines(1)
        for jj in range(0, len(seqs_arr), loop_increment):
            if output_file_format == 'tab':
                record = seqs_arr[jj].strip("\n").split("\t")
                str1 = record[1]
            else:
                str1 = seqs_arr[jj+1].strip("\n")
                read_length = len(str1)
            for i in range(read_length):
                ii = bases_dict[str1[i]]
                base_count_per_read_location[i,ii] =  base_count_per_read_location[i,ii]+1
    return base_count_per_read_location,workflow_prov

def plot_base_counts(x,a,c,t,g,n,output_fig):
        """ Generates plots of frequency of basecounts per read location. """
	import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
	max_y=max(max(a),max(c),max(t),max(g),max(n))
        min_y=min(min(a),min(c),min(t),min(g),min(n))
        min_x=min(x)
        max_x=max(x)
        fig = plt.figure(figsize=(9,6))
        ax = fig.add_subplot(111)
        frame  = ax.get_frame()
        frame.set_facecolor('1.00')    # set the frame face color to light gray
        ax.plot(x,a,'k-',x,c,'b-',x,t,'g-',x,g,'m-',x,n,'r-')
        leg=plt.legend(('A', 'C', 'T' , 'G' , 'N'),loc=(1.001,.01))
        ax.set_ylim([min_y,max_y])
        ax.set_xlim([min_x,max_x])
        ax.grid(False)
        ax.set_xlabel('Read Position')
        ax.set_ylabel('Base Count')
        ax.set_title('Position vs. Base Count')
	# set some legend properties.  All the code below is optional.  The
        # defaults are usually sensible but if you need more control, this
        # shows you how
        # the matplotlib.patches.Rectangle instance surrounding the legend
        frame  = leg.get_frame()
        frame.set_facecolor('1.00')    # set the frame face color to white
        # matplotlib.text.Text instances
        for t in leg.get_texts():
            t.set_fontsize('small')    # the legend text fontsize
        # matplotlib.lines.Line2D instances
        for l in leg.get_lines():
            l.set_linewidth(1.5)  # the legend line width
        pdf_file=output_fig+'.pdf'
        eps_file=output_fig+'.eps'
        plt.savefig(pdf_file, facecolor='w', edgecolor='w',
                orientation='landscape', papertype='letter', format='pdf',
                transparent=False, bbox_inches=None, pad_inches=0.1)
        plt.savefig(eps_file, facecolor='w', edgecolor='w',
                orientation='landscape', papertype='letter', format='eps',
                transparent=False, bbox_inches=None, pad_inches=0.1)

def create_dir(directoryname):
	""" Checks to see if a given directory exists, else creates it. """
        if os.path.exists(directoryname)==True:
                pass
        else:
                os.system ("mkdir -p " + " " + directoryname)
        return directoryname

def find_pair_end_files(file_list):
	""" Finds the corresponding paired end files in a given input list. 
            Uses '_1' and '_2' to pair them. """
        pairs=file_list
        pairs.sort()
        found_pairs=[]
        found_set=[]
        for i in range(len(pairs)):
                if pairs[i] in found_set:
                        continue
                else:
                        for j in range(len(pairs)):
                                if i==j or pairs[j] in found_set:
                                        continue
                                else:
                                        if pairs[i] not in found_set:
                                                if len(pairs[i]) == len(pairs[j]):
                                                        diff = sum(ch1 != ch2 for ch1, ch2 in zip(pairs[i],pairs[j]) )
                                                        if diff==1:
                                                                for ch1, ch2 in zip(pairs[i],pairs[j]):
                                                                        if ch1 != ch2:
                                                                                if ch1=='1' and ch2=='2' or ch1 == '2' and ch1 =='1':
                                                                                        found=[]
                                                                                        found.append(pairs[i])
                                                                                        found.append(pairs[j])
                                                                                        found_pairs.append(found)
                                                                                        found_set.append(pairs[i])
                                                                                        found_set.append(pairs[j])
                                                                                else:
                                                                                        pass
                                                        else :
                                                                if j == len(pairs)-1:
                                                                        pass
                                                                else :
                                                                        continue
                                                else:
                                                        continue
                                        else:
                                                break
        not_found_set= list(set(pairs)-set(found_set))
        return found_pairs,not_found_set

def find_unique_set(s1,s2):
	"""Reads sequence dictionaries for two pairs after preprocess step, 
	retains only the reads which belong to the same read"""
        only_headers_s1list=[]
        only_headers_s1dict={}
        only_headers_s2list=[]
        only_headers_s2dict={}
        unique_S1=[]
        unique_S2=[]
        for i in range(0,len(s1),4):
                only_headers_s1=re.match('^@\S+',str(s1[i]))
                if only_headers_s1:
			header_1 = only_headers_s1.group().strip("\n").strip("1")
                        only_headers_s1dict[header_1]=[]
                        only_headers_s1dict[header_1].append(s1[i+1])
                        only_headers_s1dict[header_1].append(s1[i+2])
                        only_headers_s1dict[header_1].append(s1[i+3])
        set_s1=set(only_headers_s1dict.keys())
        for i in range(0,len(s2),4):
                only_headers_s2=re.match('^@\S+',str(s2[i]))
                if only_headers_s2:
			header_2 = only_headers_s2.group().strip("\n").strip("2")
                        only_headers_s2dict[header_2]=[]
                        only_headers_s2dict[header_2].append(s2[i+1])
                        only_headers_s2dict[header_2].append(s2[i+2])
                        only_headers_s2dict[header_2].append(s2[i+3])
        set_s2=set(only_headers_s2dict.keys())
        unique_set= list(set_s1 & set_s2)
        for i in range(len(unique_set)):
                header_4_s1=unique_set[i]+"1"+"\n"
                header_4_s2=unique_set[i]+"2"+"\n"
                unique_S1.append(header_4_s1)
                unique_S2.append(header_4_s2)
                temp_S1=only_headers_s1dict[unique_set[i]]
                temp_S2=only_headers_s2dict[unique_set[i]]
                for j in range(len(temp_S1)):
                        unique_S1.append(temp_S1[j])
                for k in range(len(temp_S2)):
                        unique_S2.append(temp_S2[k])
        return ''.join(unique_S1),''.join(unique_S2)

def read_file_chunks(input_file,input_file2,chunk_number,nchunks,chunk_size,file_size ,format_specific_lines):
        """ Reads the indices of a file rather than physically chunking it.
            Uses this mechanism to split the input table."""
        end_pos_arr = []
        extra_lines = 0
        input_file.seek(0,1)
        current_pos = input_file.tell()
        inp1 = ''
        inptemp2 =''
        if input_file != '' :
            if chunk_number == nchunks-1:
                    inp1 = input_file.read()
                    tempcount = inp1.count("\n")
                    current_seek = input_file.tell()
            else:
                    inp1 = input_file.read(chunk_size)
                    tempcount = inp1.count("\n")

                    if tempcount < format_specific_lines :
                            extra_lines = format_specific_lines - tempcount
                            for k in range(0,extra_lines):
                                    inp1 += input_file.readline()
                            tempcount = tempcount + extra_lines
                            current_seek = input_file.tell()
                    else :
                            if tempcount % format_specific_lines == 0:
                                    if inp1[-1] != '\n':
                                            last_index = inp1.rindex('\n')
                                            back_steps = len(inp1) - last_index
                                            inp1 = inp1[0:(last_index+1)]
                                            input_file.seek(-(back_steps-1),1)
                                    current_seek = input_file.tell()
                            else :
                                    extra_lines = ((tempcount-(tempcount % format_specific_lines))+format_specific_lines) - tempcount
                                    for k in range(0,extra_lines):
                                            inp1 += input_file.readline()
                                    tempcount = tempcount + extra_lines
                                    current_seek = input_file.tell()
            if input_file2 != '':
                inptemp2 = ''.join(list(itertools.islice(input_file2,tempcount)))
                tempcount2 = inptemp2.count("\n")
                seekpos2 = input_file2.tell()
                back_steps_temp2 = 0
                diff = (tempcount2 - tempcount)
                if diff != 0 :
                    print "Error: Number of lines in paired chunks do not match"
        end_pos = current_seek
        end_pos_arr.append(end_pos)
        return end_pos,inp1,inptemp2

def find_variable(var,cmd):
	"""Given a string and varible prefix, extracts complete variable name,
	until while space is encountered""" 
	var_found = ''
	for jj in range(cmd.find(var),len(cmd)): 
        	if cmd[jj] != ' ':
                	var_found += cmd[jj]
                else:
                	break
	return var_found

def find_nth(str1, mystr, n):
    """ Finds a pattern in an input string and returns the starting index. """

    start = str1.find(mystr)
    while start >= 0 and n > 1:
        start = str1.find(mystr, start+len(mystr))
        n -=1
    return start

def run_function(cmd_string,inp_str,kk,output_dict,fh):
	""" Executes a function using the subprocess module.
            The funtion is called in the multiproc_function()."""

	try:
	        P1=Popen(cmd_string,stdin=PIPE,stdout=PIPE,stderr=fh,shell=True)
        	output_dict[str(kk)] = P1.communicate(inp_str)[0]
	except Exception as e:
		print cmd_string, " Failed!!"
		fh.write(str(e)+'\n')
	fh.close()

def multiproc_function(cmd_string,inp_str,lpp,header,err_log,stat_log):
	""" Parallelizes the commands across cores per node using the 
            multiprocessing module in python. """

        out_str=''
        err_str=''
        procs=[]
        st = 0
        en = 0
        try:
            lock = Lock()
            manager=Manager()
            output_dict=manager.dict()
            nprocs=cpu_count()
            tot_nlines=inp_str.count("\n")
            for i in range(0,nprocs):
		    fh=open(err_log+"_multiproc_"+str(i).zfill(4),'a')
                    ret=split_array(tot_nlines,nprocs,lpp,i)
                    ib=ret[0]
                    ie=ret[1]
                    nlines=ie-ib+1
                    n1=find_nth(inp_str[st:len(inp_str)], "\n", nlines)
                    n1=n1+1
                    en=st+n1
                    if header != '':
                            cmd_string_final = cmd_string + "_" + str(i)
                            procs.append(Process(target=run_function,args=(cmd_string_final,header + inp_str[st:en],i,output_dict,fh)))
                    else:
                            if inp_str[st:en] != '':
                                    procs.append(Process(target=run_function,args=(cmd_string,inp_str[st:en],i,output_dict,fh)))
                    st=en
            for i in range(0,len(procs)):
                    procs[i].start()
            for i in range(0,len(procs)):
                    procs[i].join()
 	    for i in range(0,len(procs)):
                    exit_code=procs[i].exitcode
                    yap_file_io.write_data("EXIT_CODE: "+str(exit_code)+"\n",err_log+"_multiproc_"+str(i).zfill(4))
            for i in range(0,len(procs)):
                    out_str += output_dict[str(i)]
        except Exception as e:
            print e
	del output_dict
        return out_str

def begin_end_checker(data,filename):
	""" Checks for corresponding begin and end terms to parse command sections of the configuration. """
        config_list=[]
        begin_list = []
        COMMENT_CHAR='#'
        for i in range(len(data)):
                line = data[i].strip()
                if COMMENT_CHAR in line:
                        line,comment=line.split(COMMENT_CHAR, 1)
                config_list.append(line)
        config_string= '\n'.join(config_list)
        command_begin_list = [":begin"]
        command_end_list = [":end"]
        # Syntax checker/Prilimnary checks
        dict_flag={}
        begin_status = 'off'
        end_status = 'off'
        results = []
        for i in range(len(config_list)):
                if config_list[i] in command_begin_list:
                                if begin_status == 'on':
                                        results.append(" Missing :end corresponding to :begin at line " + str(j+1))
                                        j = i
                                else:
                                        begin_status = 'on'
                                        j = i
                elif config_list[i] in command_end_list:
                        end_status = 'on'
                        j = i
                        if begin_status == 'off':
                                results.append("Missing :begin corresponding to :end at line " + str(j+1))
                        if begin_status == 'on' and end_status == 'on':
                                begin_status = 'off'
                                end_status = 'off'

        if begin_status == 'on':
                        results.append("Missing :end corresponding to :begin at line " + str(len(config_list)))

        if len(results) > 0:
                results.append("Note:Use symbol(:begin) and (:end) to define command sections,enclose variable and corresponding values in double quotes(eg \"..\")")
                begin_list = []
        else:
                begin_list=config_string.split(':begin')
        return begin_list,results

def yap_tee(initial_pipe_commands,commands2,input_file,err_log,stat_log):
	""" Emulates a Unix tee like function in python. """
	dir_path , input_file_name = os.path.split(input_file)
        input_file_name,format=os.path.splitext(input_file_name)
	random_number = str(random.random())
        if 'JOB_ID' in os.environ.keys():
                unique_jobid = os.environ['JOB_ID']
	else:
		unique_jobid = ''
        def func0(fifos,initial_pipe_commands,input_file):
                try:
                    format=os.path.splitext(input_file)[1]
                    tee1=create_tee(fifos,mode='w')
                    icommand=[]
                    num_pipe_commands=len(initial_pipe_commands)+1
                    if num_pipe_commands > 1 :
                            if format==".gz":
                                    p1=subprocess.Popen(["zcat", input_file], stdout=subprocess.PIPE)
                                    icommand.append(p1)
                            elif format==".bz2":
                                    p1=subprocess.Popen(["bzcat", input_file], stdout=subprocess.PIPE)
                                    icommand.append(p1)
                            else    :
                                    p1=subprocess.Popen(["cat", input_file], stdout=subprocess.PIPE)
                                    icommand.append(p1)
                            for i in range(1,num_pipe_commands):
                                    if( i == num_pipe_commands-1):
                                            p1=subprocess.Popen(shlex.split(initial_pipe_commands[i-1]),stdin=icommand[i-1].stdout,stdout=subprocess.PIPE)
                                            icommand.append(p1)
                                    else:
                                            p1=subprocess.Popen(shlex.split(initial_pipe_commands[i-1]),stdin=icommand[i-1].stdout,stdout=subprocess.PIPE)
                                            icommand.append(p1)
                    else :
                            if format==".gz":
                                    p1=subprocess.Popen(["zcat", input_file], stdout=subprocess.PIPE)
                                    icommand.append(p1)
                            elif format==".bz2":
                                    p1=subprocess.Popen(["bzcat", input_file], stdout=subprocess.PIPE)
                                    icommand.append(p1)
                            else    :
                                    p1=subprocess.Popen(["cat", input_file], stdout=subprocess.PIPE)
                                    icommand.append(p1)
                    tee1.write(icommand[num_pipe_commands-1].communicate()[0])
                    tee1.flush()
                    tee1.close()
                    rc=icommand[num_pipe_commands-1].poll()
                    if rc==0 :
                            pass
                except Exception as e:
                    print "\nError: yap_tee execution failed for initial pipe commands",e
	    	    yap_log.write_log(str(initial_pipe_commands).lstrip("[").rstrip("]"),input_file,'',str(e),err_log,stat_log)
        def func1(cmd1,fifoin,t_err_log):
                icommand=[]
                try:
		    fh=open(t_err_log,'a')
                    fh1=open(fifoin,'r')
                    num_cmd=len(cmd1)
                    cmd2=[]
                    for i in range(num_cmd):
                            tt=Template(cmd1[i])
                            cmd2.append(tt.substitute(os.environ))
                    for i in range(num_cmd):
                            if (num_cmd == 1) :
                               	p2=subprocess.Popen(cmd2[i],stdin=fh1,stdout=subprocess.PIPE,stderr=fh,shell=True)
                                icommand.append(p2)
                            elif (num_cmd > 1 and i == 0  ) :
                               	p2=subprocess.Popen(cmd2[i],stdin=fh1,stdout=subprocess.PIPE,stderr=fh,shell=True)
                                icommand.append(p2)
                            elif (num_cmd > 1 and i==num_cmd-1) :
                            	p2=subprocess.Popen(cmd2[i],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=fh,shell=True)
                                icommand.append(p2)
                            else :
                                p2=subprocess.Popen(cmd2[i],stdin=icommand[i-1].stdout,stdout=subprocess.PIPE,stderr=fh,shell=True)
                                icommand.append(p2)
                    icommand[num_cmd-1].communicate(icommand[num_cmd-2].stdout.read())
                    fh1.close()
                    icommand[num_cmd-2].stdout.close()
                    icommand[num_cmd-2].wait()
                    rc1=icommand[num_cmd-1].poll()
                    if rc1==0:
                            print cmd1, " Finished Successfully" , "return_code=", rc1
		    fh.close()
                except Exception as e:
                    print "\nError: yap_tee execution failed for commands",cmd1 ,e
	    	    yap_log.write_log(str(cmd1),input_file,'',str(e),err_log,stat_log)
        commands1=[]
        try:
            for i in range(len(commands2)):
                    t1=commands2[i].split('|')
                    commands1.append(t1)
            ncommands=len(commands1)
            nprocs=ncommands+1
            commands=[]
            fifos=[]
	    tmp_err_log=[]
            for i in range(0,ncommands):
                    commands.append(commands1[i])
                    temp_dir = yap_workflow_dict.yap_temp_user_dir
		    tt=Template(temp_dir)
                    fifo_name = tt.substitute(os.environ)+"/"+  input_file_name + unique_jobid + "_" + random_number + "_fifo"+ str(i)
                    fifos.append(fifo_name)
                    if not os.path.exists(fifos[i]):
                            os.mkfifo(fifos[i])      
                    tmp_err_log.append(err_log+"_yap_tee_"+str(i).zfill(4))
            lock=Lock()
            manager=Manager()
            procs=[]
	    fh=[]
            procs.append(Process(target=func0,args=(fifos,initial_pipe_commands,input_file)))
            for i in range(1,nprocs):
                    procs.append(Process(target=func1,args=(commands[i-1],fifos[i-1],tmp_err_log[i-1])))
            for i in range(nprocs):
                    procs[i].start()
            for i in range(nprocs):
                    procs[i].join()
	    for i in range(nprocs):
                    exit_code=procs[i].exitcode
                    yap_file_io.write_data("EXIT_CODE: "+str(exit_code)+"\n",tmp_err_log[i-1])
            for i in (fifos):
                    os.remove(i)
	    yap_log.merge_tee_files(str(initial_pipe_commands).lstrip("[").rstrip("]")+","+str(commands2).lstrip("[").rstrip("]"),input_file,err_log,stat_log)
        except Exception as e:
            print "\nError: yap_tee execution failed at multiprocess call ",e
	    yap_log.write_log(str(initial_pipe_commands).lstrip("[").rstrip("]")+","+str(commands2).lstrip("[").rstrip("]"),input_file,'',str(e),err_log,stat_log)
	    
def command_section_parser(data,filename,count):
	""" Parses command sections of the configuration and returns command arrays. """
	config_string= '$$'.join(data)
	begin_list=config_string.split(':begin')
	error_list=[]
        # command validator object 
        command_validator = yap_cmd_checks.yap_cmd_checks()
	if len(data) != 0:
		if len(begin_list[0])!=0:
        		count=len(begin_list[0])
        cmd_arr=[]
        for i in range(1,len(begin_list)):
                temp_list=[]
                temp_list=begin_list[i].split('$$')
                cmd=''
                cmd_name = ''
                blank = ' '
                execute_flag  = ''
                for j in range(0,(len(temp_list)-1)):
                	count+=1
                        line =temp_list[j].strip().strip('\t')
                        match_config=''
                        match_config=re.match('\"(.*)\"[ \t]*\:\=[ \t]*\"(.*)\"',line)
                        if match_config:
                                key=match_config.group(1).strip().strip('\t')
                                val=match_config.group(2).strip().strip('\t')
                                if key == "execute_command":
                                        if val == "yes":
                                                execute_flag = 'yes'
                                                # enable command validator
                                                command_validator.enable = "yes"
                                        else:
                                                # disable command validator
                                                command_validator.enable = "no"
                                elif key == "command_name":
                                        cmd += val + ' '
                                        path , cmd_name = os.path.split(val)
                                        cmd_name,ext = os.path.splitext(cmd_name)
                                        # validate command name
                                        validation_output = command_validator.is_valid_cmd_name(val)
                                        error_msg ="Error: In file: "+ filename + " Invalid command name at line " +\
                                                    str(count) +\
                                                    ". Please check if command executable exists, executable "+\
                                                    "symbolic link is not broken, if any, and it is set as executable. ==> " + line
                                        # check and print error if required
                                        check_cmd_validation_output(validation_output, error_msg, str(count),
                                                                    filename, cmd_name)
                                elif val == "no":
                                        pass
                                elif val == "yes":
                                        cmd += key + blank
                                        # validate command argument
                                        validation_output = command_validator.is_valid_cmd_args(key)
                                        error_msg ="Error: In file: "+ filename + " Invalid command argument at line " + \
                                                   str(count) + \
                                                    ". Please check if all file name and paths exist and any symbolic "+\
                                                    "links are not broken. ==> " + line
                                        # check and print error if required
                                        check_cmd_validation_output(validation_output, error_msg, str(count),
                                                                    filename, cmd_name)
                                else:
                                        cmd += key + blank + val + blank
                                        # validate command argument
                                        validation_output = command_validator.is_valid_cmd_args(key + blank + val)
                                        error_msg ="Error: In file: "+ filename + " Invalid command argument at line " + \
                                                   str(count) + \
                                                    ". Please check if all file name and paths exist and any symbolic "+\
                                                    "links are not broken. ==> " + line
                                        # check and print error if required
                                        check_cmd_validation_output(validation_output, error_msg, str(count),
                                                                    filename, cmd_name)                                        
                                        #extra check for empty invalid cuffdiff sample file name
                                        if key == 'list_of_samples_to_compare' and execute_flag == 'yes':
                                            if len(val.strip()) == 0:
                                                cuff_diff_error = "Error: In file: "+ filename + " Invalid command argument at line " +\
                                                    str(count) +\
                                                    ". Please provide a valid file name for argument: list_of_samples_to_compare. ==> " + line
                                                error_list.append(cuff_diff_error)#append error
                                        #extra check for empty cuff merge and cuff compare file name
                                        if key == 'list_of_samples' and execute_flag == 'yes':
                                            if len(val.strip()) == 0:
                                                cuff_diff_error = "Error: In file: "+ filename + " Invalid command argument at line " +\
                                                    str(count) +\
                                                    ". Please provide a valid file name for argument: list_of_samples. ==> " + line
                                                error_list.append(cuff_diff_error)#append error                                                                         
                        else:
                                if line ==':end':
                                        pass
                                elif re.match('^[\s\t]+$',line) or len(line)==0:
                                        pass
                                else:
                                        error_msg=''
                                        error_msg="Error: In file: "+ filename +\
                                                   " Check syntax at line " + str(count) + " ==> " + line
                                        error_list.append(error_msg)
                if execute_flag == 'yes':
                        cmd_arr.append([cmd_name,cmd])
        if len(error_list)!=0:
                cmd_arr = []
        return cmd_arr,error_list

def begin_end_checker1(data,filename,keywords):
	""" Checks for corresponding begin and end terms to parse command sections of the configuration. """
        command_begin_sym = keywords[0]
	command_end_sym = keywords[1]
        config_list=[]
        begin_list = []
        COMMENT_CHAR='#'
        for i in range(len(data)):
                line = data[i].strip()
                if COMMENT_CHAR in line:
                        line,comment=line.split(COMMENT_CHAR, 1)
                config_list.append(line)
        config_string= '\n'.join(config_list)
        # Syntax checker/Prilimnary checks
        dict_flag={}
        begin_status = 'off'
        end_status = 'off'
        results = []
        for i in range(len(config_list)):
                if config_list[i] == command_begin_sym:
                                if begin_status == 'on':
                                        results.append(" Missing " + command_end_sym + " corresponding to " + command_begin_sym + " at line " + str(j+1))
                                        j = i
                                else:
                                        begin_status = 'on'
                                        j = i
                elif config_list[i] == command_end_sym:
                        end_status = 'on'
                        j = i
                        if begin_status == 'off':
                                results.append("Missing " + command_begin_sym + " corresponding to " + command_end_sym + " at line " + str(j+1))
                        if begin_status == 'on' and end_status == 'on':
                                begin_status = 'off'
                                end_status = 'off'
        if begin_status == 'on':
                        results.append("Missing " + command_end_sym +" corresponding to " + command_begin_sym + " at line " + str(len(config_list)))
        return results

def command_parser(file_arr,filename):
	""" Parses the configuration files based on begin..end segments. 
            Returns syntactical and path errors if any. Else, returns the
            command list. """
        syntax_arr  = [[':begin',':end'],[':begin_tee',':end_tee']]
	meta_terms = ['input_file_type','input_directory']
        error_arr = []
	cmd_arr = []
        for i in syntax_arr:
                keywords = i
                results = begin_end_checker1(file_arr,filename,keywords)
                error_arr += results
        if len(error_arr) > 0:
                error_arr.append("Note:Use symbol(:begin) and (:end) to define command sections,enclose variable and corresponding values in double quotes(eg \"..\")")
        else:
                count =0
                tee_arr=[]
                tee_error_arr=[]
                cmd_arr=[]
                tmp_arr=[]
                tee_status='off'
                cmd_status='off'
                for i in range(len(file_arr)):
			cmd_meta_data = []
                        count=count+1
                        line = file_arr[i]
                        line = line.partition('#')[0]
                        line = line.lstrip().rstrip().rstrip('\n')
                        if line == ":begin_tee":
                                tee_status="on"
                                tee_arr.append(":begin_tee")
                        elif line == ":begin":
                                tmp_arr.append(line)
                                if tee_status=="off":
                                        tee_arr.append(":begin")
                                        cmd_status="on"
                        elif line == ":end":
                                tmp_arr.append(line)
                                if tee_status=="off":
					cmd_metadata= []
                                        tmp_arr,tee_error_arr=command_section_parser(tmp_arr,filename,(count-len(tmp_arr)))
                                        error_arr=error_arr+tee_error_arr
					if len(tmp_arr) >0:
						cmd_name = tmp_arr[0][0]
						cmd = tmp_arr[0][1]
		                                input_directory_obj = re.match( r'(.*) input_directory[\s\t]*([\S\T]*)[\s\t]*', cmd,re.M|re.I)
        		                        input_file_ext_matchobj = re.match(r'(.*) input_file_type [\s\t]*([\S\T]*)[\s\t]*', cmd,re.M|re.I)
						input_file_extension = ''
                        		        if input_file_ext_matchobj :
                                        		input_file_extension = input_file_ext_matchobj.group(2)
							cmd_metadata.append('input_file_type' + ' ' + input_file_extension  )
                               			matchobj = re.match( r'(.*) input_directory (\w*)', cmd,re.M|re.I)
                                		input_directory  = ''
                                		if matchobj:
                                        		input_directory = matchobj.group(2)
                                		else:
                                        		input_directory = "aligner_output"
						cmd_metadata.append('input_directory' + ' ' + input_directory )
						cmd= cmd.replace('input_directory' + ' ' + input_directory + ' ' ,'')
                                                cmd = cmd.replace('input_file_type' + ' ' + input_file_extension + ' ' ,'')
						tmp_arr= [[cmd_name,cmd]]
                                        	tee_arr.append(cmd_metadata)
                                        	tee_arr.append(tmp_arr)
                                        	cmd_arr.append(tee_arr)
                                        tmp_arr=[]
                                        tee_arr=[]
                                        cmd_status="off"
                        elif line==":end_tee":
                                ind=tmp_arr.index(':begin')
                                inp_arr=[]
				cmd_meta_data = []
                                inp_arr=tmp_arr[:ind]
				meta_count=0
				meta_key = ''
				meta_val = ''
                                for i in range (0, len(inp_arr)):
                                	meta_match=re.match('\"(.*)\"[ \t]*\:\=[ \t]*\"(.*)\"',inp_arr[i])
                                	if meta_match:
						meta_key , meta_val = inp_arr[i].strip("\n").replace('"','').replace(' ','').split(":=")
                                        	cmd_meta_data.append(meta_key+ " " + meta_val)
                                        	if meta_key in meta_terms:
                                                	meta_count+=1
                                	elif inp_arr[i] =='':
                                        	pass
                                	else:
                                        	error_arr.append("Error: Check syntax at line "+str(count-(len(tmp_arr))+i)+" ==> "+inp_arr[i]+"  Please follow the convention \"parameter\" := \"value\"\n")
                        	if meta_count<len(meta_terms):
                                	error_arr.append("Error: Missing terms at line " +  str(count-(len(tmp_arr))) + " ==>  Missing " + str(meta_terms) + "  in tee command section" )
                        	elif meta_count>len(meta_terms):
                                	error_arr.append("Error: Replication of terms at line " + str(count-(len(tmp_arr))) + " ==> Possible replication of terms "+str(meta_terms) )
                                tmp_arr,tee_error_arr=command_section_parser(tmp_arr[ind:],filename,(count-len(tmp_arr[ind:]))-1)
                                error_arr=error_arr+tee_error_arr
                                tee_arr.append(cmd_meta_data)
                                tee_arr.append(tmp_arr)
				if len(tmp_arr) >0:
                                	cmd_arr.append(tee_arr)
                                tmp_arr=[]
                                tee_arr=[]
                                tee_status="off"
                        else:
                                if tee_status=="on" or cmd_status=="on":
                                        tmp_arr.append(line)
                                else:
                                        pass
	if len(error_arr) > 0:
                error_arr.insert(0,"Errors found in " +filename)
                error_arr.append("Please follow the convention \"parameter\" := \"value\" for parameter\n")
		cmd_arr = []
	return cmd_arr,error_arr

def workflow_parser(data,filename,nprocs):
        """ Parses the key-value pairs in workflow configuration.
            Checks if paths are valid and for syntax errors. 
            Returns a list of dictionaries of all the workflows 
            contained in the configuration."""
        begin_list,error_list = begin_end_checker(data,filename)
        error_list2 = []
        count=0
        workflow_struct = []
        error_msg=''
        if len(begin_list) > 1 :
                expt_data=begin_list[0]
                ## BEFORE THE FIRST :begin :end BLOCK
                expt_data=begin_list[0]
                if len(expt_data)!=0:
                        workflow_dict={}
                        expt_list=[]
                        expt_list=expt_data.split('\n')
                        for i in range ( 0, (len(expt_list)-1)):
                                line = expt_list[i]
                                count+=1
                                match_expt=''
                                match_expt=re.match('\"(.*)\"[\s\t]*\:\=[\s\t]*\"(.*)\"',line)
                                if match_expt:
                                        key=match_expt.group(1).strip().strip('\t')
                                        val=match_expt.group(2).strip().strip('\t')
                                        workflow_dict[key]=val
                                else:
                                        if re.match('^[\s\t]+$',line) or len(line)==0:
                                                pass
                                        else:
                                                error_msg="Error at line " + str(count) + " ==> " + line
                                                error_list.append(error_msg)
                        workflow_struct.append(workflow_dict)
        #WORKFLOW CONFIGURATION CHECKER
        ## SCANS THROUGH BLOCK BY BLOCk ie :begin-:end
        # checks for errors
        # generates commands
        for i in range(1,len(begin_list)):
                count+=1
                temp_list=[]
                temp_list=begin_list[i].split('\n')
                workflow_dict_obj=yap_workflow_dict.workflow_dictionary()
                workflow_dict = workflow_dict_obj.create_default_wf_dict()
		workflow_dict["nprocs"]=str(nprocs)
                tmp_key_list = []
                for j in range(1,(len(temp_list)-1)):
                        line=temp_list[j].strip().strip('\t')
                        count+=1
                        match_config=''
                        match_config=re.match('\"(.*)\"[ \t]*\:\=[ \t]*\"(.*)\"',line)
                        if match_config:
                                key=match_config.group(1).strip().strip('\t')
                                val=match_config.group(2).strip().strip('\t')
                                if key in tmp_key_list:
                                        error_list.append("Error: Replication of terms at line " + str(count) + " ==> Possible replication of terms: "+ key )
                                else:
                                        workflow_dict[key]=val
                                        tmp_key_list.append(key)
                        else:
                                if line==':end':
                                        pass
                                elif re.match('^[\s\t]+$',line) or len(line)==0:
                                        pass
                                else:
                                        error_msg="Error at line " + str(count) + " ==> " + line
                                        error_list.append(error_msg)
                workflow_struct.append(workflow_dict)
                error_list2.extend(workflow_dict_obj.validate_wf_dict(workflow_dict))
        if len(error_list)!=0:
                error_list.insert(0,"Errors found in " +filename)
                error_list.append("Please follow the convention \"parameter\" := \"value\"\n")
                workflow_struct = []
        elif len(error_list2) != 0:
                error_list2.insert(0,"Errors found in " +filename)
                workflow_struct = []
                error_list = error_list2
        return workflow_struct,error_list

def rename_barcode(barcode):
    """ Returns the barcode information or "no_barcode_specified" 
        if there are no barcodes"""
    if barcode == "no_barcode_specified":
      	barcode_value = ''
    else:
       	barcode_value = barcode
    return barcode_value

def check_open_file_desc():
    """ Checks if a particular file descriptor is still open. """
    pid=os.getpid()
    ppid=os.getppid()
    procs=check_output("lsof -w -Ff -p "+str(ppid))

def check_cmd_validation_output(output, error_msg, count, file, command_name):
    """ Prints all the Warnings associated with a particular command. """
    # print values in global list
    global missing_path_errors
    global basename_warnings
    if output == False:
        missing_path_errors.append(error_msg)
    elif str(output).find("basename") >=0:
        # split the basename warnings
        basename_warn_data = str(output).split(":")
        for warning in basename_warn_data:
            if len(warning) > 0:
                warning_message = "Warning: At Line: " + count + " in file: "+file+". Files were found using " + warning + \
                                  ". Please make sure that command: " + command_name + " can work with basenames."
                # add warnings to the global list
                basename_warnings.append(warning_message)

def split_array_old(totnum,nprocs,myrank):
    """ Returns a split array based on the number of processes and ranks. """
    n2=math.floor(totnum/nprocs)
    remain=math.fmod(totnum,nprocs)
    if (myrank < remain) :
        n1=n2+1
        ib=(myrank)*n1
        ie=ib+(n1-1)
    else :
        n1=n2
        ib=(myrank)*n1+remain
        ie=ib+n1-1
    return int(ib),int(ie)

def split_files_each_proc(file_arr,nprocs):
	""" Returns array that distributes samples across all processors. """
        ntot = len(file_arr)
        post_proc_file_arr = []
        for i in range(0,nprocs):
                each_proc_arr = []
                ib,ie = split_array_old(ntot,nprocs,i)
                if i == 0:
                        max_no = (ie-ib)+1
                for j in range(ib,ie+1):
                        each_proc_arr.append(j)
                if len(each_proc_arr) > max_no:
                        max_no = len(each_proc_arr)
                elif len(each_proc_arr) < max_no :
                        for k in range(0,max_no-(len(each_proc_arr))):
                                each_proc_arr.append("no file")
                max_no = len(each_proc_arr)
                post_proc_file_arr.append(each_proc_arr)
        return post_proc_file_arr

def split_array(totnum,nprocs,lpp,myrank):
    """ Determines how the data must be split based on available resources. """
    totnum = totnum/lpp
    n2=int(math.floor(totnum/nprocs))
    remain=int(math.fmod(totnum,nprocs))
    if (myrank < remain) :
        n1=n2+1
        ib1=(myrank)*n1
        ie1=ib1+(n1-1)
    else :
        n1=n2
        ib1=(myrank)*n1+remain
        ie1=ib1+n1-1
    ib=ib1*lpp
    ie=((ie1+1)*lpp)-1
    ret=[]
    ret.append(ib)
    ret.append(ie)
    return ret

def get_filesplit_info(file_struct,inp_files_list,file_chunk_size,nprocs,format_specific_lines):
	""" Returns a list containing information specific to how the file is 
            split: file size, chunk size, number of chunks, number of lines 
            skipped according to format (eq: fastq = 4). """
        file_split_struct=[]
	for k in range(0,len(file_struct)):
        	if  file_struct[k] != 'no file':
			file_split_info=[]
                        each_file = inp_files_list[file_struct[k]][0]
                        chunk_size = (int(file_chunk_size))*1024*1024
                        file_size = get_file_size(each_file)
                        nchunks = file_size / chunk_size
			nprocs=int(nprocs)
                        if nchunks < nprocs:
                                nchunks = nprocs
                                chunk_size = file_size / nchunks
                        if nchunks > nprocs:
                                cpp=long(math.ceil(float(nchunks)/float(nprocs)))
                                nchunks=(nprocs*cpp)
                                chunk_size = file_size / nchunks
                        file_split_info.append(file_size)
                        file_split_info.append(chunk_size)
                        file_split_info.append(nchunks)
                        file_split_info.append(format_specific_lines)
                        if nchunks < nprocs:
                                print " Number of processors provided= ",nprocs
                                print " File contains less then ",nprocs," sequences"
                                print  " Please set the number of processors to smaller of equal to number of sequences"
                                print  "each processor should get atlest one sequence to process "
                                exit()
                	file_split_struct.append(file_split_info)
	return file_split_struct

def file_cleanup():
    """ Removes the temporary log directories after they have been merged. """
    os.system("rm -r " + yap_workflow_dict.temp_dir_path)
    os.system("rm -r " + yap_workflow_dict.err_log_path + "/*_log_temp")
    os.system("rm -r " + yap_workflow_dict.stat_log_path + "/*_log_temp")
