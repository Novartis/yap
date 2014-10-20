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
import sys
import re
import glob
import time
import thread
import yap_tools 
import yap_file_io 
import yap_log 
import random
import yap_workflow_dict as wd
from subprocess import PIPE, Popen
def create_thread(thread_name, ss, fifo_name, err_log, stat_log):
        '''
	function to write data to a file or fifo
	'''

	op1 = open(fifo_name, 'a')
	try:
        	op1.writelines(ss)
        except Exception as e:
        	yap_log.write_log("Threading " + thread_name, fifo_name,'EXCEPTION', str(e), err_log, stat_log)
def align_func(aligner_cmd,aligner_out_basename,err_log,stat_log):
	'''
	function to execute aligner command through subprocess
	'''

	try:
            pbow = Popen(aligner_cmd,stdout=PIPE,stderr=PIPE,shell=True,close_fds=True)
            std_out, std_err = pbow.communicate()
            exit_code = pbow.returncode
            seqs_str1 = ''
            seqs_str2 = ''
            yap_log.write_log(aligner_cmd, aligner_out_basename,exit_code, std_err, err_log, stat_log)
            return std_out
        except Exception as e:
            print e
def run_aligner(
        seqs_str1,
        seqs_str2,
        fname,
        chunk_number,
	myrank,
        workflow_prov,
        err_log,
        stat_log):
    '''
    Runs alignment for chunk data or file,
    polishes commands for input/output paths and creates pipes.
    '''
    aligner_out_str = ''
    p1 = []
    n_cmd = len(wd.aligner_cmd_arr)
    for i in range(0, n_cmd):
        scratch_temp_dir = wd.yap_temp_user_dir
        cmd_type = wd.aligner_cmd_arr[i][0]
        cmd_meta_data = wd.aligner_cmd_arr[i][1]
        temp_arr = wd.aligner_cmd_arr[i][2]
        aligner_cmd_name = temp_arr[0][0]
        aligner_cmd = temp_arr[0][1]
        aligner_dir_path, file_name = os.path.split(fname)
        aligner_cmd = aligner_cmd.replace("output_directory", aligner_dir_path)
        aligner_cmd = aligner_cmd.replace("output_file", fname)
        aligner_cmd = aligner_cmd.replace("sample_name", file_name)
        aligner_cmd = aligner_cmd.replace("input_files_path", wd.input_files_path)
        aligner_cmd_tmp = ''
        pipe_var1 = ''
        pipe_var2 = ''
        pipe1_basename = ''
        pipe2_basename = ''
        if aligner_cmd.find("pipe1") != -1:
            pipe_var1 = yap_tools.find_variable("pipe1", aligner_cmd)
        if aligner_cmd.find("pipe2") != -1:
            pipe_var2 = yap_tools.find_variable("pipe2", aligner_cmd)
        if wd.data_distribution_method == "file_based":
            if pipe_var1 != '':
                aligner_cmd = aligner_cmd.replace(
                    pipe_var1, " " + seqs_str1 + " ", 1)
            if pipe_var2 != '':
                aligner_cmd = aligner_cmd.replace(
                    pipe_var2, " " + seqs_str2 + " ", 1)
            aligner_out_str = align_func(aligner_cmd,fname,err_log,stat_log)
        else:
            pipe1_basename = pipe_var1.replace('pipe1',scratch_temp_dir + "/" + file_name + wd.job_id + "_" + wd.random_id + '_pipe_' + str(chunk_number)+ "_" + str(i) + "_1")
            pipe2_basename = pipe_var2.replace('pipe2',scratch_temp_dir + "/" + file_name + wd.job_id + "_" + wd.random_id + '_pipe_' + str(chunk_number)+ "_" + str(i) + "_2")
            if pipe_var1 != '':
                aligner_cmd = aligner_cmd.replace(pipe_var1, " " + pipe1_basename + " ", 1)
            if pipe_var2 != '':
                aligner_cmd = aligner_cmd.replace(pipe_var2, " " + pipe2_basename + " ")
            if pipe_var1 != '' and pipe_var2 != '':
                if os.path.exists(pipe1_basename) != True:
                	os.mkfifo(pipe1_basename)
                try:
                    thread.start_new_thread(create_thread,("thread1",seqs_str1,pipe1_basename,err_log,stat_log))
                except:
                    print "Error: unable to start thread1"
                
		if os.path.exists(pipe2_basename) != True:
                	os.mkfifo(pipe2_basename)
                try:
                    thread.start_new_thread(create_thread,("thread2",seqs_str2,pipe2_basename,err_log,stat_log))
                except:
                    print "Error: unable to start thread2"
                aligner_out_str = align_func(aligner_cmd,fname,err_log,stat_log)
                os.unlink(pipe1_basename)
                os.unlink(pipe2_basename)
            elif pipe_var1 != '' and pipe_var2 == '':
                if os.path.exists(pipe1_basename) != True:
                	os.mkfifo(pipe1_basename)
                try:
                    thread.start_new_thread(create_thread,("thread1",seqs_str1,pipe1_basename,err_log,stat_log))
                except:
                    print "Error: unable to start thread1"
                aligner_out_str = align_func(aligner_cmd,fname,err_log,stat_log)
                os.unlink(pipe1_basename)

            elif pipe_var2 != '' and pipe_var1 == '':
                if os.path.exists(pipe2_basename) != True:
                	os.mkfifo(pipe2_basename)
                try:
                    thread.start_new_thread(create_thread,("thread2",seqs_str2,pipe2_basename,err_log,stat_log))
                except:
                    print "Error: unable to start thread"
                aligner_out_str = align_func(aligner_cmd,fname,err_log,stat_log)
                os.unlink(pipe2_basename)
            else:
                pbow = Popen(aligner_cmd,stdout=PIPE,stderr=PIPE,shell='True',close_fds='True')
                aligner_out_str, std_err = pbow.communicate()
                exit_code = pbow.returncode
                yap_log.write_log( aligner_cmd, fname, exit_code, std_err, err_log, stat_log)
        if aligner_cmd_name != '':
            workflow_prov.append(aligner_cmd)
        alignment_outfile_pos = 0
        format_ext = ''
        alignment_file_ext = ''
        while alignment_outfile_pos != -1:
            aligner_output_filename = ''
            alignment_outfile_pos = aligner_cmd.rfind(fname)
            for jj in range(alignment_outfile_pos, len(aligner_cmd)):
                if aligner_cmd[jj] != ' ':
                    aligner_output_filename += aligner_cmd[jj]
                else:
                    break
            aligner_output_filename_base, alignment_file_ext = os.path.splitext(
                aligner_output_filename)
            if alignment_file_ext == '.gz' or alignment_file_ext == 'bz2':
                aligner_output_filename_base, format_ext = os.path.splitext(
                    aligner_output_filename_base)
            if format_ext == '.sam' or alignment_file_ext == '.sam':
                alignment_outfile_pos = -1
            elif format_ext == '.bam' or alignment_file_ext == '.bam':
                alignment_outfile_pos = -1
            else:
                aligner_cmd_tmp = aligner_cmd[0:alignment_outfile_pos]
                aligner_cmd = aligner_cmd_tmp
	#treat tophat as exception; search if the file has been created and pass for sorting
        #this is because tophat's output filename cannot be customized according YAP output structure
        if re.search('tophat', aligner_cmd_name) is not None:
            if os.path.exists(aligner_dir_path + "/accepted_hits.bam"):
                aligner_output_filename = aligner_dir_path + "/accepted_hits.bam"
            if os.path.exists(aligner_dir_path + "/accepted_hits.sam"):
                aligner_output_filename = aligner_dir_path + "/accepted_hits.sam"
	#After alignment pass data for sorting 
        sort_alignment_output(chunk_number,aligner_cmd_name,aligner_cmd,aligner_output_filename,workflow_prov,err_log,stat_log)
    return aligner_out_str, workflow_prov

def execute_merge_alignment(
        final_output_name,
        sort_input_files_arr,
        file_type,
        file_name,
        barcode,
        sort_files_cmd,
        workflow_prov,
        err_log,
        stat_log):
    '''
    Executes merge data commands for alignment output data.
    '''
    sort_cmd_input = ''
    sort_input_files_new_arr = []
    if file_type != "sam":
        if len(sort_input_files_arr) > 0:
            if len(sort_input_files_arr) == 1:
                os.rename(sort_input_files_arr[0], final_output_name)
                workflow_prov.append(
                    'RENAMED FILE ' +
                    sort_input_files_arr[0] +
                    ' TO ' +
                    final_output_name)
            else:
                for z in range(0, len(sort_input_files_arr)):
                    sort_cmd_input += sort_input_files_arr[z].strip("\n") + " "
                if wd.alignment_sort_order == "unsorted":
                    sort_files_cmd = "samtools cat -o " + \
                        final_output_name + ' ' + sort_cmd_input
                else:
                    sort_files_cmd = sort_files_cmd + ' ' + \
                        final_output_name + ' ' + sort_cmd_input
                str_out = "*" * 50 + "MERGE ALIGNMENT STARTED" + "\t" + str(time.strftime("%Y/%m/%d %H:%M:%S", time.localtime())) + "*" * 50 + "\n"
                yap_file_io.write_data(str_out, err_log)
                yap_file_io.write_data(str_out, stat_log)
                pmerge = Popen(sort_files_cmd, stdout=PIPE, stderr=PIPE, shell='True')
                std_out, std_err = pmerge.communicate()
                exit_code = pmerge.returncode
                yap_log.write_log(sort_files_cmd, str(sort_input_files_arr).lstrip(
                    '[').rstrip(']'), exit_code, std_err, err_log, stat_log)
                str_out = "*" * 50 + "MERGE ALIGNMENT FINISHED" + "\t" + str(time.strftime("%Y/%m/%d %H:%M:%S", time.localtime())) + "*" * 50 + "\n"
                yap_file_io.write_data(str_out, err_log)
                yap_file_io.write_data(str_out, stat_log)
                if sort_files_cmd != '':
                    workflow_prov.append(sort_files_cmd)
                if exit_code != 0:
                    if file_name == '':
                        print "Error: chunk merge sort failed for barcode=", barcode, "\n"
                    else:
                        print "Error: chunks  merge sort failed for Filename=", file_name, "barcode=", barcode, "\n"
                for z in range(0, len(sort_input_files_arr)):
                    os.remove(sort_input_files_arr[z])
    else:
        if len(sort_input_files_arr) > 0:
            if len(sort_input_files_arr) == 1:
                os.rename(sort_input_files_arr[0], final_output_name)
                workflow_prov.append(
                    'RENAMED FILE ' +
                    sort_input_files_arr[0] +
                    ' TO ' +
                    final_output_name)
            else:
                str_out = "*" * 50 + "MERGE ALIGNMENT STARTED" + "\t" + str(time.strftime("%Y/%m/%d %H:%M:%S", time.localtime())) + "*" * 50 + "\n"
                yap_file_io.write_data(str_out, err_log)
                yap_file_io.write_data(str_out, stat_log)
                for z in range(0, len(sort_input_files_arr)):
                    sam_file_name = sort_input_files_arr[z]
                    sam_file_name_base, ext = os.path.splitext(sam_file_name)
                    sam_to_bam_cmd = "samtools view -bhS " + \
                        sam_file_name + " -o " + sam_file_name_base + ".bam"
                    pconv = Popen(
                        sam_to_bam_cmd, stdout=PIPE, stderr=PIPE, shell='True')
                    std_out, std_err = pconv.communicate()
                    exit_code = pconv.returncode
                    yap_log.write_log(
                        sam_to_bam_cmd,
                        final_output_name,
                        exit_code,
                        std_err,
                        err_log,
                        stat_log)
                    std_out = ""
                    std_err = ""
                    exit_code = 0
                    if exit_code != 0:
                        print " Sam to bam conversion failed"
                    sort_input_files_new_arr.append(
                        sam_file_name_base + '.bam')
                    os.remove(sam_file_name)
                for z in range(0, len(sort_input_files_new_arr)):
                    sort_cmd_input += sort_input_files_new_arr[
                        z].strip("\n") + " "
                if  wd.alignment_sort_order == "unsorted":
                    sort_files_cmd = "samtools cat -o - " + sort_cmd_input + \
                        " | samtools view -h - -o " + final_output_name
                else:
                    sort_files_cmd = sort_files_cmd + ' - ' + ' ' + sort_cmd_input + \
                        " | samtools view -h - -o " + final_output_name
                std_out = ''
                std_err = ''
                pmerge = Popen(
                    sort_files_cmd, stdout=PIPE, stderr=PIPE, shell='False')
                std_out, std_err = pmerge.communicate()
                exit_code = pmerge.returncode
                if sort_files_cmd != '':
                    workflow_prov.append(sort_files_cmd)
                if exit_code != 0:
                    if file_name == '':
                        print "Error: chunk merge sort failed for barcode=", barcode, "\n"
                    else:
                        print "Error: chunks  merge sort failed for Filename=", file_name, "barcode=", barcode, "\n"
                yap_log.write_log(sort_files_cmd, str(sort_input_files_arr).lstrip(
                    '[').rstrip(']'), exit_code, std_err, err_log, stat_log)
                str_out = "*" * 50 + "MERGE ALIGNMENT FINISHED" + "\t" + str(time.strftime("%Y/%m/%d %H:%M:%S", time.localtime())) + "*" * 50 + "\n"
                yap_file_io.write_data(str_out, err_log)
                yap_file_io.write_data(str_out, stat_log)
                for z in range(0, len(sort_input_files_new_arr)):
                    os.remove(sort_input_files_new_arr[z])

def merge_alignment_output(
        file_basecount_dict,
        workflow_prov,
        err_log,
        stat_log):
    '''
    Prepares for merging of chunk alignment output into single file,
    generates chunk list and commands based on filename and format.
    '''
    sort_order = wd.alignment_sort_order
    for filename_key in file_basecount_dict.iterkeys():
        path_name, file_name = os.path.split(filename_key)
        barcode_basecount_dict = file_basecount_dict[filename_key]
        file_type = ''
        sort_files_cmd = ''
        suffix_ext = ''
        for barcode in barcode_basecount_dict.iterkeys():
            barcode_value = yap_tools.rename_barcode(barcode)
            aligner_dir_path = wd.workflow_output_path + "/" + file_name + "/" + barcode + "/" + "aligner_output"
            aligner_output_base = aligner_dir_path + "/" + barcode_value + "_" + file_name
            if barcode_value != '':
                aligner_output_base = aligner_dir_path + "/" + "aligner_" + file_name + "_" + barcode_value
                aligner_final_output_base = aligner_dir_path + "/" + file_name + "_" + barcode_value
            else:
                aligner_output_base = aligner_dir_path + "/" + "aligner_" + file_name
                aligner_final_output_base = aligner_dir_path + "/" + file_name
            for aligner_output_key in wd.aligner_output_key_arr:
                sort_input_files_arr = []
                aligner_output_suffix = ''
                if len(aligner_output_key.split("output_file")) > 1:
	                aligner_output_suffix = aligner_output_key.split("output_file")[1]
                else:
                	if len(aligner_output_key.split("accepted_hits")) > 1:
                        	aligner_output_suffix = aligner_output_key.split("accepted_hits")[1]
                                aligner_output_base = aligner_dir_path + "/" + "accepted_hits"
                aligner_output_suffix = aligner_output_key.split("output_file")[1]
                aligner_output_suffix, suffix_ext = os.path.splitext(aligner_output_suffix)
                if aligner_output_suffix == '.sam' or aligner_output_suffix == '.bam':
                    aligner_output_suffix = ''
                else:
                    aligner_output_suffix, suffix_ext = os.path.splitext(aligner_output_suffix)
                output_with_suffix = aligner_output_base + "*" + aligner_output_suffix
                final_output_with_suffix = aligner_final_output_base + aligner_output_suffix
                if len(glob.glob(output_with_suffix + "*.sam")) > 0:
                    file_type = "sam"
                if len(glob.glob(output_with_suffix + "*.bam")) > 0:
                    file_type = "bam"
                if wd.alignment_sort_order == 'both':
                    sort_input_files_arr = glob.glob(
                        output_with_suffix + "*queryname*")
                    final_output_name = final_output_with_suffix + "_" + "queryname" + "." + file_type
                    sort_files_cmd = 'samtools merge -n '
                    execute_merge_alignment(
                        final_output_name,
                        sort_input_files_arr,
                        file_type,
                        file_name,
                        barcode,
                        sort_files_cmd,
                        workflow_prov,
                        err_log,
                        stat_log)
                    sort_input_files_arr = glob.glob(output_with_suffix + "*coordinate*")
                    final_output_name = final_output_with_suffix + "_" + "coordinate" + "." + file_type
                    sort_files_cmd = 'samtools merge '
                    execute_merge_alignment(
                        final_output_name,
                        sort_input_files_arr,
                        file_type,
                        file_name,
                        barcode,
                        sort_files_cmd,
                        workflow_prov,
                        err_log,
                        stat_log)
                else:
                    sort_input_files_arr = glob.glob(
                        output_with_suffix + "*" + file_type)
                    if wd.alignment_sort_order == 'unsorted':
                        final_output_name = final_output_with_suffix
                    else:
                        final_output_name = final_output_with_suffix + "_" + wd.alignment_sort_order
                    if file_type == "sam":
                        final_output_name += ".sam"
                    if file_type == "bam":
                        final_output_name += ".bam"
                    if file_type == 'sam' or file_type == "bam":
                        sort_cmd_input = ''
                        if wd.alignment_sort_order == 'queryname':
                            sort_files_cmd = 'samtools merge -n '
                        if wd.alignment_sort_order == 'coordinate':
                            sort_files_cmd = 'samtools merge '
                        execute_merge_alignment(
                            final_output_name,
                            sort_input_files_arr,
                            file_type,
                            file_name,
                            barcode,
                            sort_files_cmd,
                            workflow_prov,
                            err_log,
                            stat_log)
                if file_type == "sam" or file_type == "bam":
                    rm_cmd = "rm " + aligner_dir_path + "/" + barcode + \
                        "*" + aligner_output_suffix + "*." + file_type
                    if len(glob.glob(aligner_dir_path + "/" + barcode + "*" + aligner_output_suffix + "*." + file_type)) > 0:
                        prm = Popen(rm_cmd, shell='False')
                        std_out, std_err = prm.communicate()
                        exit_code = prm.returncode
                        if exit_code != 0:
                            if file_name == '':
                                print "Error: chunk merge clean up after sort failed for barcode=", barcode, "\n"
                            else:
                                print "Error: chunks clean up after merge sort failed for filename =", file_name, "barcode=", barcode, "\n"
    return workflow_prov

def sort_alignment_output(
        chunk_number,
        aligner_cmd_name,
        aligner_cmd,
        aligner_output_filename,
        workflow_prov,
        err_log,
        stat_log):
    '''
    sorts alignment output based on coordinate or queryname or both
    '''
    initial_pipe_commands = []
    sort_flag = 'False'
    format_ext = ''
    alignment_file_ext = ''
    sort_file_ext = ''
    after_sort_cmd = ''
    name_sort_cmd = ''
    coordinate_sort_cmd = ''
    if wd.alignment_sort_order != 'unsorted':
	#prepare sort commands, input/output filenames based on file extension
        aligner_dir_path = os.path.split(aligner_output_filename)[0]
        aligner_dir_path_tmp, filename_base = os.path.split(aligner_output_filename)
        if filename_base.find('queryname') == -1 and filename_base.find('coordinate') == -1:
            if os.path.exists(aligner_output_filename):
                aligner_output_filename_base, alignment_file_ext = os.path.splitext(
                    aligner_output_filename)
                if alignment_file_ext == '.gz' or alignment_file_ext == 'bz2':
                    aligner_output_filename_base, format_ext = os.path.splitext(
                        aligner_output_filename_base)
                sort_aligner_output_filename = aligner_output_filename_base + \
                    "_" + wd.alignment_sort_order
                if format_ext == '.sam' or alignment_file_ext == '.sam':
                    initial_pipe_commands = ['samtools view -bhS -']
                    sort_file_ext = ".sam"
                    name_sort_cmd = " samtools sort -on -m 100000000 - " + \
                        aligner_output_filename_base + "_" + "queryname"
                    coordinate_sort_cmd = " samtools sort -o -m 100000000 - " + \
                        aligner_output_filename_base + "_" + "coordinate"
                    name_sort_cmd = " samtools sort -on -m 100000000 - " + \
                        aligner_output_filename_base + "_" + "queryname"
                    coordinate_sort_cmd = " samtools sort -o -m 100000000 - " + \
                        aligner_output_filename_base + "_" + "coordinate"
                    after_sort_cmd = ' | samtools view -h - -o '
                else:
                    name_sort_cmd = " samtools sort -n -m 100000000 - "
                    coordinate_sort_cmd = " samtools sort -m 100000000 - "
                sort_aligner_output_filename += sort_file_ext
                if format_ext == '.sam' or format_ext == '.bam' or alignment_file_ext == '.sam' or alignment_file_ext == '.bam':
                    sort_flag = 'True'
        if sort_flag == 'True':
            if wd.alignment_sort_order == 'queryname':
                sort_commands = [
                    name_sort_cmd +
                    after_sort_cmd +
                    sort_aligner_output_filename]
            if wd.alignment_sort_order == 'coordinate':
                sort_commands = [
                    coordinate_sort_cmd +
                    after_sort_cmd +
                    sort_aligner_output_filename]
            if wd.alignment_sort_order == 'queryname' or wd.alignment_sort_order == 'coordinate':
		#if sort order is both; then call yap tee functionality
                yap_tools.yap_tee(initial_pipe_commands, sort_commands,aligner_output_filename, err_log, stat_log)
                if len(initial_pipe_commands) > 0:
                    workflow_prov.append(initial_pipe_commands[0] + sort_commands[0])
                else:
                    workflow_prov.append(sort_commands[0])
                rm_cmd = "rm " + aligner_output_filename
                prm = Popen(rm_cmd, shell='True').wait()
                if prm != 0:
                    print "Error: chunks clean up after merge sort failed for Filename=", aligner_output_filename, " chunk number=", chunk_number, "\n"
            if wd.alignment_sort_order == 'both':
                sort_queryname_output = aligner_output_filename_base + "_queryname" + sort_file_ext
                sort_coordinate_output = aligner_output_filename_base + "_coordinate" + sort_file_ext
                sort_commands = [name_sort_cmd + after_sort_cmd + sort_queryname_output,coordinate_sort_cmd + after_sort_cmd + sort_coordinate_output]
		#if sort order is both; then call yap tee functionality
                yap_tools.yap_tee(initial_pipe_commands, sort_commands,aligner_output_filename, err_log, stat_log)
                if len(initial_pipe_commands) > 0:
                    workflow_prov.append(
                        initial_pipe_commands[0] + sort_commands[0])
                    workflow_prov.append(
                        initial_pipe_commands[0] + sort_commands[1])
                else:
                    workflow_prov.append(sort_commands[0])
                    workflow_prov.append(sort_commands[1])
                try:
                    os.remove(aligner_output_filename)
                except Exception as e:
                    print "Error: chunks clean up after merge sort failed for Filename=", aligner_output_filename, " chunk number=", chunk_number, "cmd", "\n"
                    print e
    return workflow_prov

def regroup_files(
	regroup_arr,
        workflow_prov):
    '''
    Merge the alignment data based on regroup sample information given in-
    workflow configuration. Takes care of sorting and merging different file outputs
    into a single sample output. 
    '''
    regroup_title = regroup_arr[0]
    regroup_files_arr = regroup_arr[1]
    err_string = ""
    stat_string = ""
    for i in range(len(regroup_files_arr)):
        if os.path.exists(wd.err_log_path + "/" + regroup_files_arr[i] + "_err.log"):
            err_string += " " + wd.err_log_path + "/" + regroup_files_arr[i] + "_err.log"
        if os.path.exists(wd.stat_log_path + "/" + regroup_files_arr[i] + "_stat.log"):
            stat_string += " " + wd.stat_log_path + "/" + regroup_files_arr[i] + "_stat.log"
        if os.path.exists(wd.err_log_path + "/" + regroup_files_arr[i] + "_log_temp"):
            os.system("rm -rf " + wd.err_log_path + "/" + regroup_files_arr[i] + "_log_temp")
        if os.path.exists(wd.stat_log_path + "/" + regroup_files_arr[i] + "_log_temp"):
            os.system("rm -rf " + wd.stat_log_path + "/" + regroup_files_arr[i] + "_log_temp")
    if err_string != " " + wd.err_log_path + "/" + regroup_title + "_err.log":
        os.system("cat" +err_string +"> " + wd.err_log_path + "/" + regroup_title + "_err.log")
        os.system("rm" + err_string)
    os.system("mkdir " + wd.err_log_path + "/" + regroup_title + "_log_temp")
    if stat_string != " " + wd.stat_log_path + "/" + regroup_title + "_stat.log":
        os.system("cat" + stat_string + "> " + wd.stat_log_path + "/" + regroup_title + "_stat.log")
        os.system("rm" + stat_string)
    os.system("mkdir " + wd.stat_log_path + "/" + regroup_title + "_log_temp")
    err_log = wd.err_log_path + "/" + regroup_title + "_err.log"
    stat_log = wd.stat_log_path + "/" + regroup_title + "_stat.log"
    #write regrouped file information to log file
    regroup_merge_log = wd.regroup_output_path + "/" + "regroup_files.log"
    fw = open(regroup_merge_log, 'a')
    fw.write(regroup_title + " ")
    for i in range(0, len(regroup_files_arr)):
        if i == len(regroup_files_arr) - 1:
            fw.write(regroup_files_arr[i])
        else:
            fw.write(regroup_files_arr[i] + ",")
    fw.write("\n")
    fw.close()
    for aligner_output_key in wd.aligner_output_key_arr:
        sort_input_files_arr = []
        sort_input_files_arr_queryname = []
        sort_input_files_arr_coordinate = []
        aligner_output_suffix = ''
        if len(aligner_output_key.split("output_file")) > 1:
		aligner_output_suffix = aligner_output_key.split("output_file")[1]
        else:
        	if len(aligner_output_key.split("accepted_hits")) > 1:
                	aligner_output_suffix = aligner_output_key.split("accepted_hits")[1]
        aligner_output_suffix, suffix_ext = os.path.splitext(aligner_output_suffix)
        if aligner_output_suffix == '.sam' or aligner_output_suffix == '.bam':
            aligner_output_suffix = ''
        else:
            aligner_output_suffix, suffix_ext = os.path.splitext(
                aligner_output_suffix)
        sort_input_files_arr = []
        sort_input_files_arr_queryname = []
        sort_input_files_arr_coordinate = []
        for filename_key in regroup_files_arr:
            file_name = filename_key
            file_type = ''
            sort_files_cmd = ''
            suffix_ext = ''
            for barcode in wd.barcode_dict.iterkeys():
                barcode_value = yap_tools.rename_barcode(barcode)
                aligner_dir_path = wd.workflow_output_path + "/" + file_name + "/" + barcode + "/" + "aligner_output"
                regroup_aligner_dir_path = wd.regroup_output_path + "/" + regroup_title + "/" + barcode + "/" + "aligner_output"
                aligner_output_base = aligner_dir_path + \
                    "/" + barcode_value + "_" + file_name
                if barcode_value != '':
                    aligner_output_base = aligner_dir_path + "/" + \
                        "aligner_" + file_name + "_" + barcode_value
                    aligner_final_output_base = regroup_aligner_dir_path + \
                        "/" + regroup_title + "_" + barcode_value
                else:
                    aligner_output_base = aligner_dir_path + "/" + "aligner_" + file_name
                    aligner_final_output_base = regroup_aligner_dir_path + "/" + regroup_title
                if re.search('accepted_hits',aligner_output_key) != None:
			aligner_output_base = aligner_dir_path + "/" + "accepted_hits"
                output_with_suffix = aligner_output_base + "*" + aligner_output_suffix
                final_output_with_suffix = aligner_final_output_base + aligner_output_suffix
                if len(glob.glob(output_with_suffix + "*.sam")) > 0:
                    file_type = "sam"
                if len(glob.glob(output_with_suffix + "*.bam")) > 0:
                    file_type = "bam"
                if wd.alignment_sort_order == 'both':
                    sort_input_files_arr_queryname.extend(
                        glob.glob(output_with_suffix + "*queryname*"))
                    sort_input_files_arr_coordinate.extend(
                        glob.glob(output_with_suffix + "*coordinate*"))
                else:
                    sort_input_files_arr.extend(
                        glob.glob(output_with_suffix + "*" + file_type))
        if wd.alignment_sort_order == 'both':
            final_output_name_queryname = final_output_with_suffix + "_" + "queryname" + "." + file_type
            final_output_name_coordinate = final_output_with_suffix + "_" + "coordinate" + "." + file_type
            sort_files_cmd = 'samtools merge -n '
            execute_merge_alignment(
                final_output_name_queryname,
                sort_input_files_arr_queryname,
                file_type,
                file_name,
                barcode,
                sort_files_cmd,
                workflow_prov,
                err_log,
                stat_log)
            sort_files_cmd = 'samtools merge '
            execute_merge_alignment(
                final_output_name_coordinate,
                sort_input_files_arr_coordinate,
                file_type,
                file_name,
                barcode,
                sort_files_cmd,
                workflow_prov,
                err_log,
                stat_log)
        else:
            if wd.alignment_sort_order == 'unsorted':
                final_output_name = final_output_with_suffix
            else:
                final_output_name = final_output_with_suffix + "_" + wd.alignment_sort_order
            if file_type == "sam":
                final_output_name += ".sam"
            if file_type == "bam":
                final_output_name += ".bam"
            if file_type == 'sam' or file_type == "bam":
                sort_cmd_input = ''
                if wd.alignment_sort_order == 'queryname':
                    sort_files_cmd = 'samtools merge -n '
                if wd.alignment_sort_order == 'coordinate':
                    sort_files_cmd = 'samtools merge '
                execute_merge_alignment(
                    final_output_name,
                    sort_input_files_arr,
                    file_type,
                    file_name,
                    barcode,
                    sort_files_cmd,
                    workflow_prov,
                    err_log,
                    stat_log)
        if file_type == "sam" or file_type == "bam":
            rm_cmd = "rm " + aligner_dir_path + "/" + barcode + \
                "*" + aligner_output_suffix + "*." + file_type
            if len(glob.glob(aligner_dir_path + "/" + barcode + "*" + aligner_output_suffix + "*." + file_type)) > 0:
                prm = Popen(rm_cmd, shell='False').wait()
                if prm != 0:
                    if file_name == '':
                        print "Error: chunk merge clean up after sort failed for barcode=", barcode, "\n"
                    else:
                        print "Error: chunks clean up after merge sort failed for filename =", file_name, "barcode=", barcode, "\n"
    return workflow_prov
