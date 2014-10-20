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

import os
import time
import glob
import re
import yap_tools
import yap_file_io
import yap_log
import yap_workflow_dict as wd
from subprocess import Popen,PIPE

def run_postprocess(
        postprocess_cmd_arr,
        file_basecount_dict,
        workflow_prov,
        err_log,
        stat_log):
    ''' 
    Prepare postprocess command with input/output paths according to sample name, 
    pass commands to yap_tee or subprocess for execution.
    '''
    if wd.regroup_output =='yes':
    	workflow_output_path = wd.workflow_output_path + "/regroup_output"
    else:
	workflow_output_path= wd.workflow_output_path
    for zz in range(0, len(postprocess_cmd_arr)):
        postprocess_tee_arr = []
        postprocess_nontee_arr = []
        initial_pipe_commands = []
        postprocess_temp_arr = []
        cmd_type = postprocess_cmd_arr[zz][0]
        cmd_meta_data = postprocess_cmd_arr[zz][1]
        postprocess_temp_arr = postprocess_cmd_arr[zz][2]
        input_file_extension = ''
	pipe1=''
	pipe2=''
        #set default input directory for postprocess stage as aligner_output
        #user can specify "postprocess_output" through  configuration file
        input_directory = "aligner_output" 
        for kk in range(0, len(cmd_meta_data)):
            input_meta_data = cmd_meta_data[kk].split(" ")
            if input_meta_data:
                if re.search('input_file_type', input_meta_data[0]) is not None:
                    input_file_extension = input_meta_data[1]
                if re.search('input_directory', input_meta_data[0]) is not None:
                    input_directory = input_meta_data[1]
        '''iterate over filename and barcode, fetch files from the source directory,
        file extensions and python glob module'''
        for filename_key in file_basecount_dict.iterkeys():
	    #fetch original input file pair for this sample
            for tmp_arr in wd.inp_files_list:
            	if tmp_arr[2] == filename_key:
			#store it in variables pipe1 and pipe2
                	pipe1 = tmp_arr[0]
                        pipe2 = tmp_arr[1]
                	break
            postprocess_input_file_arr = []
            path, file_name = os.path.split(filename_key)
            if wd.run_preprocess_analysis == "no" and wd.run_reference_alignment == "no" and wd.run_postprocess_analysis == "yes":
            	file_name, extension = os.path.splitext(file_name)
            	file_name, extension = os.path.splitext(file_name)
                file_basecount_dict[filename_key] = {
                    'no_barcode_specified': []}
            for barcode in file_basecount_dict[filename_key]:
                barcode_value = yap_tools.rename_barcode(barcode)
                aligner_dir_path = ''
                postprocess_dir_path = ''
                aligner_dir_path = workflow_output_path + "/" + file_name + "/" + barcode + "/" + input_directory
                postprocess_input = aligner_dir_path + "/" + "*" + input_file_extension + "*"
                postprocess_input_file_arr = glob.glob(postprocess_input)
            	if wd.run_preprocess_analysis == "no" and wd.run_reference_alignment == "no" and wd.run_postprocess_analysis == "yes":
            		if input_directory == "aligner_output":
                        	aligner_dir_path = path
                                postprocess_input = filename_key
                                temp_arr = glob.glob(aligner_dir_path + "/" + "*" + input_file_extension + "*")
                                if temp_arr > 0:
                                	for k in temp_arr:
                                        	if k == postprocess_input:
                                                	postprocess_input_file_arr = [postprocess_input]
                if input_file_extension == '':
                    postprocess_input_file_arr = []
                postprocess_dir_path = workflow_output_path + "/" + file_name + "/" + barcode + "/" + "postprocess_output"
                postprocess_input_file_arr.sort()
                if (len(postprocess_input_file_arr) > 0):
                    if wd.run_preprocess_analysis == "no" and wd.run_reference_alignment == "no" and wd.run_postprocess_analysis == "yes":
                    	if input_directory == "aligner_output":
                        	if postprocess_input_file_arr[0] == filename_key:
                            		pass
                        	else:
                            		break
                    for k in range(0, len(postprocess_temp_arr)):
                        postprocess_cmd = postprocess_temp_arr[k][1]
                        postprocess_cmd_name = postprocess_temp_arr[k][0]
                        if file_name == '':
                            if barcode_value != '':
                                postprocess_outfile = postprocess_dir_path + "/" + \
                                    barcode_value + "_" + postprocess_cmd_name
                            else:
                                postprocess_outfile = postprocess_dir_path + \
                                    "/" + postprocess_cmd_name
                        else:
                            if barcode_value != '':
                                postprocess_outfile = postprocess_dir_path + "/" + \
                                    file_name + "_" + barcode_value + \
                                    "_" + postprocess_cmd_name
                            else:
                                postprocess_outfile = postprocess_dir_path + \
                                    "/" + file_name + "_" + \
                                    postprocess_cmd_name
                        #replace generic keywords with appropriate file path variables
                        postprocess_cmd = postprocess_cmd.replace(
                            'input_directory', '')
                        postprocess_cmd = postprocess_cmd.replace(
                            'input_file_type' + ' ' + input_file_extension, '')
                        postprocess_cmd = postprocess_cmd.replace(
                            "aligner_output", '')
                        postprocess_cmd = postprocess_cmd.replace(
                            "postprocess_output", '')
                        postprocess_cmd = postprocess_cmd.replace(
                            'output_file', postprocess_outfile)
                        postprocess_cmd = postprocess_cmd.replace(
                            'output_directory', postprocess_dir_path)
                        postprocess_cmd = postprocess_cmd.replace(' =', '=')
                        postprocess_cmd = postprocess_cmd.replace(
                            'sample_name', file_name)
                        if re.search("file_based_input", postprocess_cmd) is not None:
                            postprocess_cmd = postprocess_cmd.replace(
                                'file_based_input',
                                postprocess_input_file_arr[0])
                            postprocess_nontee_arr = [
                                postprocess_cmd_name, postprocess_cmd]
                        elif re.search("directory_based_input", postprocess_cmd) is not None:
                            postprocess_cmd = postprocess_cmd.replace(
                                'directory_based_input', aligner_dir_path)
                            postprocess_nontee_arr = [
                                postprocess_cmd_name, postprocess_cmd]
                        else:
                            postprocess_tee_arr.append(postprocess_cmd)
                    workflow_prov.append(
                        "INPUT: " + postprocess_input_file_arr[0])
                    for kk in postprocess_tee_arr:
                        if kk != '':
                            workflow_prov.append(kk)
                    if len(postprocess_tee_arr) != 0:
		        #pass commands to yap_tee function
                        yap_tools.yap_tee(
                            initial_pipe_commands,
                            postprocess_tee_arr,
                            postprocess_input_file_arr[0],
                            err_log,
                            stat_log)
                    if len(postprocess_nontee_arr) != 0:
		        #pass commands to non_tee function which uses subproces
                        run_postprocess_nontee(
                            postprocess_nontee_arr,
                            workflow_prov,
                            err_log,
                            stat_log)
                else:
                    if file_name == '':
                        print "Warning: No aligner output for barcode = ", barcode, " ...skipping the postprocess step for command : \n", postprocess_temp_arr, "\n"
                    else:
                        print "Warning: No aligner output for filename = ", file_name, "  barcode = ", barcode, " ...skipping the postprocess step for command: \n", postprocess_temp_arr, "\n"
    return workflow_prov

def run_postprocess_nontee(
        postprocess_compare_arr,
        workflow_prov,
        err_log,
        stat_log):
    '''Executes command through subrocess, writes log information'''
    cmd = postprocess_compare_arr[1]
    try:
        std_out = ''
        std_err = ''
        prun = Popen(cmd, stdout=PIPE, stderr=PIPE, shell='True')
        std_out, std_err = prun.communicate()
        exit_code = prun.returncode
        yap_log.write_log(str(cmd), '', str(exit_code), std_err, err_log, stat_log)
        if cmd != '':
            workflow_prov.append(cmd)
    except Exception as e:
        print "Error : while running postprocess command ", cmd, "\n"
        print e
        yap_log.write_log(str(cmd), '', '', str(e), err_log, stat_log)
    return workflow_prov


def separate_postprocess_cmds(postprocess_cmd_arr):
    '''
    Separate commands into two categories, samples which needs to be
    analyzed individually or in groups of multiple samples, based on 
    key words present in the command string.
    '''
    postprocess_cmd_arr_new = []
    postprocess_compare_file_cmd_arr = []
    for j in range(0, len(postprocess_cmd_arr)):
        file_list_comp_matchobj = ''
        sample_file_list_matchobj = ''
        matchobj = ''
        cmd_type = postprocess_cmd_arr[j][0]
        cmd_meta_data = postprocess_cmd_arr[j][1]
        postprocess_temp_arr = postprocess_cmd_arr[j][2]
        postprocess_cmd_name = postprocess_temp_arr[0][0]
        postprocess_cmd = postprocess_temp_arr[0][1]
        file_list_comp_matchobj = re.match(
            r'(.*) list_of_samples_to_compare[\s\t]*([\S\T]*)[\s\t]*',
            postprocess_cmd,
            re.M | re.I)
        sample_file_list_matchobj = re.match(
            r'(.*) list_of_samples[\s\t]*([\S\T]*)[\s\t]*',
            postprocess_cmd,
            re.M | re.I)
        if file_list_comp_matchobj or sample_file_list_matchobj:
            postprocess_compare_file_cmd_arr.append(postprocess_cmd_arr[j])
        else:
            postprocess_cmd_arr_new.append(postprocess_cmd_arr[j])
    return postprocess_cmd_arr_new, postprocess_compare_file_cmd_arr

def get_input_sets(
        compare_file_set,
        compare_file_name,
        input_directory,
        input_file_extension):
    '''fetch the input file sets based on source directory and input file extension provided,
    returns list files matching the criteria'''
    compare_final_set = []
    if wd.regroup_output == 'yes':
	workflow_output_path = wd.workflow_output_path+ "/regroup_output"
    else:
	workflow_output_path = wd.workflow_output_path
    for j in range(0, len(compare_file_set)):
        compare_postprocess_dir_path = ''
        compare_file = ''
        if compare_file_set[j] != '' or (re.search('(\w+)', compare_file_set[j])) is not None:
            compare_file_temp = compare_file_set[j]
            path, compare_file = os.path.split(compare_file_temp)
            compare_file, extension = os.path.splitext(compare_file)
            compare_file, extension = os.path.splitext(compare_file)
            if wd.run_preprocess_analysis == "no" and wd.run_postprocess_analysis == "yes" and wd.run_reference_alignment == "no":
                if input_directory == 'postprocess_output': 
                    compare_postprocess_dir_path = workflow_output_path + "/" + "*" + compare_file + "*" + "/" + "no_barcode_specified" + "/" + input_directory + "/*" + input_file_extension + "*"
                else:
                    compare_postprocess_dir_path = compare_file_temp
            else:
                compare_postprocess_dir_path = workflow_output_path + "/" + "*" + compare_file + "*" + "/" + "no_barcode_specified" + "/" + input_directory + "/*" + input_file_extension + "*"
        if glob.glob(compare_postprocess_dir_path):
            temp_file = ''
            temp_file = glob.glob(compare_postprocess_dir_path)[0]
            if temp_file != '':
                compare_final_set.append(temp_file)
    return compare_final_set

def get_postprocess_file_compare_cmd_arr(
        postprocess_compare_file_cmd_arr,
        inp_file_list):
    '''Polish the postprocess commands which require multiple samples,
    for input/output paths'''
    postprocess_cmd_arr = postprocess_compare_file_cmd_arr
    postprocess_compare_arr = []
    temp_sample_file_list = ''
    list_samples_to_compare_dict = wd.list_of_samples_to_compare
    list_samples_dict = wd.list_of_samples
    if wd.regroup_output == 'yes':
	workflow_output_path = wd.workflow_output_path+ "/regroup_output"
    else:
	workflow_output_path = wd.workflow_output_path
    for j in range(0, len(postprocess_cmd_arr)):
        input_file_extension = ''
        file_list_comp_matchobj = ''
        sample_file_list_matchobj = ''
        input_directory_matchobj = ''
        compare_file_name = ''
        file_compare_list = []
        postprocess_cmd = postprocess_cmd_arr[j][1]
        postprocess_cmd_name = postprocess_cmd_arr[j][0]
        cmd_type = postprocess_cmd_arr[j][0]
        cmd_meta_data = postprocess_cmd_arr[j][1]
        postprocess_temp_arr = postprocess_cmd_arr[j][2]
        postprocess_cmd_name = postprocess_temp_arr[0][0]
        postprocess_cmd = postprocess_temp_arr[0][1]
        input_file_extension = ''
        list_delimiter = ''
        sample_name = ''
        #set default input directory for postprocess stage as aligner_output
        #user can specify "postprocess_output" through  configuration
        input_directory = "aligner_output" 
        list_delimiter_obj = re.match(
            r'(.*) list_delimiter[\s\t]*([\S\T]*)[\s\t]*', 
            postprocess_cmd,
            re.M | re.I)
        if list_delimiter_obj:
            list_delimiter = list_delimiter_obj.group(2).strip("\n")
            #postprocess_cmd = postprocess_cmd.replace(list_delimiter, '')
        #check for command input/output keywords from configuration variables 
        for kk in range(0, len(cmd_meta_data)):
            input_meta_data = cmd_meta_data[kk].split(" ")
            if input_meta_data:
                if re.search('input_file_type', input_meta_data[0]) is not None:
                    input_file_extension = input_meta_data[1] #fetch user provided input file type
                if re.search('input_directory', input_meta_data[0]) is not None:
                    input_directory = input_meta_data[1] # fetch user provided input directory
        if postprocess_cmd_name in list_samples_to_compare_dict:
            cmd_list_samples_to_compare = list_samples_to_compare_dict[
                postprocess_cmd_name]
            compare_file_name = cmd_list_samples_to_compare[0]
            file_compare_list = cmd_list_samples_to_compare[1]
            file_list_comp_matchobj = 'True'
        if postprocess_cmd_name in list_samples_dict:
            cmd_list_samples = list_samples_dict[postprocess_cmd_name]
            compare_file_name = cmd_list_samples[0]
            file_compare_list = cmd_list_samples[1]
            sample_file_list_matchobj = 'True'
        if wd.run_preprocess_analysis == "no" and wd.run_postprocess_analysis == "yes" and wd.run_reference_alignment == "no":
            inp_file_list1 = []
            for i in inp_file_list:
                compare_file_temp = i[0]
                path, compare_file = os.path.split(compare_file_temp)
                compare_file, extension = os.path.splitext(compare_file)
                compare_postprocess_dir_path = compare_file_temp
                for jj in glob.glob(compare_postprocess_dir_path):
                    inp_file_list1.append([jj, '', jj, ''])
            inp_file_list = inp_file_list1
        if compare_file_name == "all": #if all the sample space to be analyzed one-to-all
            files_temp_list = []
            if file_list_comp_matchobj == 'True':
                for j in range(0, len(inp_file_list)):
                    files_temp_list.append(inp_file_list[j][2])
                #generate the list of file sets one-to-all
                file_compare_list = (
                    generate_file_comparisons(files_temp_list))
            if sample_file_list_matchobj == 'True':
                for j in range(0, len(inp_file_list)):
                    files_temp_list.append(inp_file_list[j][2])
                file_compare_list.append([files_temp_list])
        '''Iterate over file sets to be analyzed together, Check if the files 
        exist according to source directory and input type'''
        for i in range(0, len(file_compare_list)):
            compare_file_set = file_compare_list[i]
            input_string = ''
            iter_var = 0
            cmd_compare_dir_path = workflow_output_path + "/" + postprocess_cmd_name + "_output"
            #create directory structure by command name and group number
            if os.path.exists(cmd_compare_dir_path) != True:
                os.system("mkdir" + " " + cmd_compare_dir_path)
            output_compare_dir_path = cmd_compare_dir_path + \
                "/" + postprocess_cmd_name + "_group_" + str(i + 1)
            err_log = wd.err_log_path+ '/' + postprocess_cmd_name + "_group_" + str(i + 1) + "_err.log"  
            stat_log = wd.stat_log_path + '/' + postprocess_cmd_name + "_group_" + str(i + 1) + "_stat.log" 
            if os.path.exists(output_compare_dir_path) != True:
                os.system("mkdir" + " " + output_compare_dir_path)
            temp_file_compare_list = cmd_compare_dir_path + \
                "/" + postprocess_cmd_name + "_input_file_list"
            temp_sample_file_list = output_compare_dir_path + "/" + \
                postprocess_cmd_name + "group" + \
                str(i + 1) + "_input_file_list"
            postprocess_outfile = output_compare_dir_path + \
                "/" + postprocess_cmd_name
            for jk in range(0, len(compare_file_set)):
                input_set = []
                iter_var = iter_var + 1
                compare_file = ''
                command_out_dir = ''
                postprocess_cmd_new = ''
                input_set = get_input_sets(
                    compare_file_set[jk],
                    compare_file_name,
                    input_directory,
                    input_file_extension)
                if file_list_comp_matchobj == 'True':
                    if len(input_set) > 0:
                        if iter_var == 1:
                            sample_name = postprocess_cmd_name + \
                                "_Group_" + str(i + 1)
                            yap_file_io.write_data(postprocess_cmd_name +
                                       "_Group_" +
                                       str(i +
                                           1) +
                                       ":" +
                                       "\n", temp_file_compare_list)
                        yap_file_io.write_data(
                            "Set" + str(jk + 1) + "=", temp_file_compare_list)
                        for zz in range(0, len(input_set)):
                            yap_file_io.write_data(
                                "\t" +
                                input_set[zz] +
                                "\n",
                                temp_file_compare_list)
                            input_string += input_set[zz] + ","
                        yap_file_io.write_data("\n", temp_file_compare_list)
                        if list_delimiter != '':
                            input_string = input_string.strip(
                                ",") + " " + list_delimiter + " "
                        else:
                            input_string = input_string.strip(",") + " "
                if sample_file_list_matchobj == 'True':
                    if len(input_set) > 0:
                        sample_name = postprocess_cmd_name + "_Group_" + str(1)
                        for kk in(input_set):
                            yap_file_io.write_data(kk + "\n", temp_sample_file_list)
            #polish commands according to type of command file comparison sets 
            if file_list_comp_matchobj == 'True':
                if input_string != '':
                    input_string = input_string.strip(list_delimiter + " ")
                    postprocess_cmd = postprocess_cmd.replace(
                        compare_file_name, '')
                    postprocess_cmd_new = postprocess_cmd.replace(
                        'input_directory', '')
                    postprocess_cmd_new = postprocess_cmd_new.replace(
                        input_directory, '')
                    postprocess_cmd_new = postprocess_cmd_new.replace(
                        'list_of_samples_to_compare', '')
                    postprocess_cmd_new = postprocess_cmd_new + input_string
                    postprocess_cmd_new = postprocess_cmd_new.replace(
                        'output_file', postprocess_outfile)
                    postprocess_cmd_new = postprocess_cmd_new.replace(
                        'output_directory', output_compare_dir_path)
                    postprocess_cmd_new = postprocess_cmd_new.replace(
                        'list_delimiter' + ' ' + list_delimiter, '')
                    postprocess_cmd_new = postprocess_cmd_new.replace(
                        'sample_name', sample_name)
            #polish commands according to type of command where all files together passed as a list  
            if sample_file_list_matchobj == 'True':
                if glob.glob(temp_sample_file_list):
                    postprocess_cmd = postprocess_cmd.replace(
                        compare_file_name, '')
                    postprocess_cmd_new = postprocess_cmd.replace(
                        'input_directory', '')
                    postprocess_cmd_new = postprocess_cmd_new.replace(
                        input_directory, '')
                    postprocess_cmd_new = postprocess_cmd_new.replace(
                        'output_file', postprocess_outfile)
                    postprocess_cmd_new = postprocess_cmd_new.replace(
                        'output_directory', output_compare_dir_path)
                    postprocess_cmd_new = postprocess_cmd_new.replace(
                        'list_of_samples', temp_sample_file_list)
                    postprocess_cmd_new = postprocess_cmd_new.replace(
                        'list_delimiter' + ' ' + list_delimiter , '')
                    postprocess_cmd_new = postprocess_cmd_new.replace(
                        'sample_name', sample_name)
            postprocess_compare_arr.append(
                [postprocess_cmd_name, postprocess_cmd_new, err_log, stat_log])
    return postprocess_compare_arr

def generate_file_comparisons(files_list):
    '''Generate unique one-to-all file sets from input file list'''
    file_compare_list = []
    for i in files_list:
        for j in files_list:
            if i < j:
                file_compare_list.append([[i], [j]])
    return file_compare_list
