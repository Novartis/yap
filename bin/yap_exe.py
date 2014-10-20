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
import numpy
import re
import time
import glob
from subprocess import PIPE, Popen
import yap_tools 
import yap_preprocess
import yap_aligner
import yap_file_io 
import yap_log 
import yap_workflow_dict as wd
'''
This script executes preprocess commands and prepared data for alignment,
for both chunkbased and file based workflows. 
'''
def execute_chunk(
        input_file_list_local,
        inp1,
        inp2,
        chunk_number,
	myrank,
        workflow_prov,
	eqp_dict):
    '''
    Executes preprocess commands for chunked data and passes to the alignment stage
    Takes chunked input data, filename list, chunk number, rank of the processor     
    and provenance list to append log data.
    ''' 
    # variable declaration
    input_filename_local = input_file_list_local[0]
    input_filename_local_2 = input_file_list_local[1]
    file_name = input_file_list_local[2]
    err_chunk_file = wd.err_log_path + "/" + file_name + \
        "_log_temp/" + file_name + "_" + str(chunk_number).zfill(6)
    stat_chunk_file = wd.stat_log_path + "/" + file_name + \
        "_log_temp/" + file_name + "_" + str(chunk_number).zfill(6)
    myhost = os.getenv('HOSTNAME')
    yap_file_io.write_data("HOSTNAME: " + str(myhost) + "\n", err_chunk_file)
    yap_file_io.write_data("HOSTNAME: " + str(myhost) + "\n", stat_chunk_file)
    yap_file_io.write_data("CHUNK NUMBER: " + str(chunk_number) + "\n", err_chunk_file)
    yap_file_io.write_data("CHUNK NUMBER: " + str(chunk_number) + "\n", stat_chunk_file)
    seqs_arr1 = []
    seqs_arr2 = []
    read_length = wd.max_read_length
    barcode_seqstruct_dict1 = {}
    barcode_seqstruct_dict2 = {}
    barcode_output_dict = {}
    aligner_out_str = ''
    sort_order = ''
    barcode_flag = 'False'
    sort_order = wd.alignment_sort_order
    # convert the input data based on format given in workflow configuration
    if wd.input_file_format == "qseq" or wd.input_file_format != wd.preprocess_output_file_format:
        inp1 = yap_tools.convert_format(inp1)
        if wd.paired_end_data == 'yes':
            inp2 = yap_tools.convert_format(inp2)
    if wd.run_preprocess_analysis == 'yes':
	str_out = "-"*20 + "PREPROCESS STARTED" +"\t" + str(time.strftime("%Y/%m/%d %H:%M:%S", time.localtime())) + "-"*20 + "\n"
	yap_file_io.write_data(str_out,err_chunk_file)
	yap_file_io.write_data(str_out,stat_chunk_file)
        # Run barcode splitter as first preprocess step
        for jj in range(0, len(wd.preprocess_cmd_arr)):
            preprocess_cmd_name = wd.preprocess_cmd_arr[jj][2][0][0]
            preprocess_cmd = wd.preprocess_cmd_arr[jj][2][0][1]
            if re.search('fastx_barcode_splitter', preprocess_cmd_name) is not None:
                barcode_flag = 'True'
                print "Entering " + preprocess_cmd_name + " : Filename=", input_filename_local, " chunk number=", chunk_number, "\n"
		str_out= "YAP_COMMAND: " + preprocess_cmd + "\n" + "INPUT FILE: " + input_filename_local
		yap_file_io.write_data(str_out,err_chunk_file)
		yap_file_io.write_data(str_out,stat_chunk_file)
                barcode_seqstruct_dict1, workflow_prov = yap_preprocess.fastx_barcode_splitter(
                    inp1, wd.preprocess_output_file_format, preprocess_cmd, workflow_prov, err_chunk_file, stat_chunk_file)
                yap_file_io.write_data("_" * 30 + "\n", err_chunk_file)
                yap_file_io.write_data("_" * 30 + "\n", stat_chunk_file)
                barcode_seqstruct_dict1["no_barcode_specified"] = ''
                print "Exiting " + preprocess_cmd_name + " : Filename=", input_filename_local, " chunk number=", chunk_number, "\n"
                if wd.paired_end_data == 'yes':
                    print "Entering " + preprocess_cmd_name + " : Filename=", input_filename_local_2, " chunk number=", chunk_number, "\n"
		    str_out= "YAP_COMMAND: " + preprocess_cmd + "\n" + "INPUT FILE: " + input_filename_local_2
		    yap_file_io.write_data(str_out,err_chunk_file)
		    yap_file_io.write_data(str_out,stat_chunk_file)
                    barcode_seqstruct_dict2, workflow_prov = yap_preprocess.fastx_barcode_splitter(
                        inp2,wd.preprocess_output_file_format , preprocess_cmd, workflow_prov, err_chunk_file, stat_chunk_file)
                    yap_file_io.write_data("_" * 30 + "\n", err_chunk_file)
                    yap_file_io.write_data("_" * 30 + "\n", stat_chunk_file)
                    barcode_seqstruct_dict2["no_barcode_specified"] = ''
                    print "Exiting " + preprocess_cmd_name + " : Filename=", input_filename_local, " chunk number=", chunk_number, "\n"
                break
        if barcode_flag == 'False':
            #if no barcode command; then create dictionary with one barcode tag
            barcode_seqstruct_dict1["no_barcode_specified"] = inp1
            barcode_seqstruct_dict2["no_barcode_specified"] = inp2
    else:
        #if no preprocess stage specified; then create dictionary with one barcode tag
        barcode_seqstruct_dict1["no_barcode_specified"] = inp1
        barcode_seqstruct_dict2["no_barcode_specified"] = inp2
    #iterate over the barcode dictionary 
    for barcode, inp1 in barcode_seqstruct_dict1.iteritems():
        run_unique_reads = 'False'
        barcode_value = yap_tools.rename_barcode(barcode)
        if wd.paired_end_data == "yes":
            inp2 = barcode_seqstruct_dict2[barcode]
        preprocessed_data_dict = {}
	#intialize matrix for basecount analysis
        aligner_output_str_local = ''
        basecount_matrix_local1 = numpy.zeros(
            (int(read_length), 5), dtype=numpy.int)
        basecount_matrix_local2 = numpy.zeros(
            (int(read_length), 5), dtype=numpy.int)
        barcode_output_dict.setdefault(barcode, [basecount_matrix_local1, basecount_matrix_local2])
        #set output file paths
        barcode_dir_path = wd.workflow_output_path + "/" + file_name + "/" + barcode
        preprocess_dir_path = barcode_dir_path + "/" + "preprocess_output"
        if wd.data_distribution_method != "file_based":
            if barcode_value != '':
                preprocess_out_filename1 = preprocess_dir_path + "/" + barcode_value + "_" + file_name + \
                    "_" + str(chunk_number).zfill(6) + "_" + \
                    str(myrank) + "_preprocessed_data_1.txt"
                preprocess_out_filename2 = preprocess_dir_path + "/" + barcode_value + "_" + file_name + \
                    "_" + str(chunk_number).zfill(6) + "_" + \
                    str(myrank) + "_preprocessed_data_2.txt"
            else:
                preprocess_out_filename1 = preprocess_dir_path + "/" + file_name + "_" + \
                    str(chunk_number).zfill(6) + "_" + \
                    str(myrank) + "_preprocessed_data_1.txt"
                preprocess_out_filename2 = preprocess_dir_path + "/" + file_name + "_" + \
                    str(chunk_number).zfill(6) + "_" + \
                    str(myrank) + "_preprocessed_data_2.txt"
        else:
            if barcode_value != '':
                preprocess_out_filename1 = preprocess_dir_path + "/" + \
                    "preprocess_data" + "_" + file_name + \
                    "_" + barcode_value + "_1.txt"
                preprocess_out_filename2 = preprocess_dir_path + "/" + \
                    "preprocess_data" + "_" + file_name + \
                    "_" + barcode_value + "_2.txt"
            else:
                preprocess_out_filename1 = preprocess_dir_path + "/" + \
                    "preprocess_data" + "_" + file_name + "_1.txt"
                preprocess_out_filename2 = preprocess_dir_path + "/" + \
                    "preprocess_data" + "_" + file_name + "_2.txt"
        aligner_dir_path = barcode_dir_path + "/" + "aligner_output"
        if barcode_value != '':
            aligner_output_filename = aligner_dir_path + "/" + "aligner_" + \
                file_name + "_" + barcode_value + \
                "_" + str(chunk_number).zfill(6)
        else:
            aligner_output_filename = aligner_dir_path + "/" + \
                "aligner_" + file_name + "_" + str(chunk_number).zfill(6)

        for jj in range(0, len(wd.preprocess_cmd_arr)):
            preprocess_cmd_name = wd.preprocess_cmd_arr[jj][2][0][1]
            preprocess_cmd = wd.preprocess_cmd_arr[jj][2][0][1]
            # skip fastqc and fastq screen and barcode splitter as they are
            # already executed
            if (re.search('fastqc', preprocess_cmd_name) is not None) or (re.search('fastq_screen', preprocess_cmd_name) is not None)or(re.search('fastx_barcode_splitter',
                                                                                                                                                  preprocess_cmd_name) is not None):
                pass
            else:
                if re.search('calculate_basecount_metrics', preprocess_cmd_name) is not None:
		    #excecute basecount calculation
                    basecount_matrix_local1, workflow_prov = yap_tools.qc_basecount(
                        inp1, workflow_prov)
                    basecount_matrix_local2, workflow_prov = yap_tools.qc_basecount(
                        inp2, workflow_prov)
                elif re.search('fastx_clipper', preprocess_cmd_name) is not None:
		    """
		    Check for fastx clipper as special case and execute.
		    This is because fastx clipper execution has been optimized by providing contaminants for every file,
		    instead of just applying contaminants universally. 
		    """ 
                    run_unique_reads = 'True'
                    if input_filename_local in wd.contaminant_dict.keys():
                        contaminants_arr1 = wd.contaminant_dict[
                            input_filename_local]
                        print "Entering " + preprocess_cmd_name + " : Filename=", input_filename_local, " chunk number=", chunk_number, "\n"
                        index = 0
                        for index in range(0, len(contaminants_arr1)):
			    #iterate over all the contaminants for this file 
                            fastx_clipper_cmd = preprocess_cmd
                            contaminant1 = contaminants_arr1[index].strip("\n")
                            if inp1 != '':
                                cont_replace = " -a " + contaminant1
                                fastx_clipper_cmd = fastx_clipper_cmd.replace(
                                    'pipe1', " - ") + " -a " + contaminant1
                                inp1 = yap_tools.multiproc_function(
                                    fastx_clipper_cmd, inp1, int(
                                        wd.format_specific_lines), '', err_chunk_file, stat_chunk_file)
                                yap_log.merge_multiproc_files(
                                    fastx_clipper_cmd,
                                    input_filename_local,
                                    barcode,
                                    err_chunk_file,
                                    stat_chunk_file)
                            if inp1 == '':
                                break
                        print "Exiting " + preprocess_cmd_name + " : Filename=", input_filename_local, " chunk number=", chunk_number, "\n"
                    if wd.paired_end_data == 'yes':
                        if input_filename_local_2 in wd.contaminant_dict.keys():
			    #repeat fastx clipper for the paired end
                            contaminants_arr2 = wd.contaminant_dict[
                                input_filename_local_2]
                            print "Entering " + preprocess_cmd_name + " : Filename=", input_filename_local_2, " chunk number=", chunk_number, "\n"
                            index = 0
                            for index in range(0, len(contaminants_arr2)):
                                fastx_clipper_cmd = preprocess_cmd
                                contaminant2 = contaminants_arr2[
                                    index].strip("\n")
                                if inp2 != '':
                                    cont_replace = " -a " + contaminant2
                                    fastx_clipper_cmd = fastx_clipper_cmd.replace(
                                        'pipe1',
                                        " - ") + " -a " + contaminant2
                                    inp2 = yap_tools.multiproc_function(
                                        fastx_clipper_cmd, inp2, int(
                                            wd.format_specific_lines), '', err_chunk_file, stat_chunk_file)
                                    yap_log.merge_multiproc_files(
                                        fastx_clipper_cmd,
                                        input_filename_local_2,
                                        barcode,
                                        err_chunk_file,
                                        stat_chunk_file)
                                if inp2 == '':
                                    break
                            print "Exiting " + preprocess_cmd_name + " : Filename=", input_filename_local_2, " chunk number=", chunk_number, "\n"
                elif re.search('eqp_rename_reads',preprocess_cmd_name) != None:
                        # this section renames reads according to specific format, applies to in-house use, neglect otherwise
                        inp1_arr = inp1.splitlines(1)
                        inp1=''
                        inp2_arr = inp2.splitlines(1)
                        inp2=''
                        read_count=1
                        if wd.data_distribution_method == "file_based":
                                if eqp_dict.has_key("eqp_read_counter"):
                                        if len(eqp_dict["eqp_read_counter"]) > 0:
                                                file_name, read_count = eqp_dict["eqp_read_counter"]
                                                if file_name !=  input_filename_local:
                                                        read_count = 1
                        format_lines = int(wd.format_specific_lines)
                        for i in range(0,len(inp1_arr),format_lines):
                                if wd.paired_end_data == 'yes':
                                        if (len(inp1_arr[i+1].strip("\n").replace('A','')) >= 5) and (len(inp2_arr[i+1].strip("\n").replace('A','')) >= 5) and (len(inp1_arr[i+1].strip("\n").replace('T','')) >= 5) and (len(inp2_arr[i+1].strip("\n").replace('T','')) >= 5) :
                                                inp1 += '@F'+str(read_count).zfill(9)+'/1'+'\n'
                                                inp2 += '@F'+str(read_count).zfill(9)+'/2'+'\n'
                                                for jj in range (1,format_lines):
                                                        inp1 += inp1_arr[i+jj]
                                                        inp2 += inp2_arr[i+jj]
                                                read_count += 1
                                else:
                                        if (len(inp1_arr[i+1].strip("\n").replace('A','')) >= 5) and (len(inp1_arr[i+1].strip("\n").replace('T','')) >= 5):
                                                inp1_arr[i] = '@F'+str(read_count).zfill(9)+'/1'+'\n'
                                                for jj in range (1,format_lines):
                                                        inp1 += inp1_arr[i+jj]
                                                read_count += 1
                        eqp_dict["eqp_read_counter"] = [ input_filename_local, read_count]
                        inp1_arr = []
                        inp2_arr = []
                else:
		    #set the flag to remove umatched pair after preprocesing 
                    run_unique_reads = 'True'
                    print "Entering " + preprocess_cmd_name + " : Filename=", input_filename_local, " chunk number=", chunk_number, "\n"
		    #for all other preprocess commands execute this section
                    if inp1 != '':
                        preprocess_cmd = preprocess_cmd.replace('pipe1', ' - ')
                        inp1 = yap_tools.multiproc_function(
                            preprocess_cmd, inp1, int(
                                wd.format_specific_lines), '', err_chunk_file, stat_chunk_file)
                        yap_log.merge_multiproc_files(
                            preprocess_cmd,
                            input_filename_local,
                            barcode,
                            err_chunk_file,
                            stat_chunk_file)
                    print "Exiting " + preprocess_cmd_name + " : Filename=", input_filename_local, " chunk number=", chunk_number, "\n"
                    if wd.paired_end_data == 'yes':
                        preprocess_cmd = preprocess_cmd.replace('pipe1', ' - ')
                        print "Entering " + preprocess_cmd_name + " : Filename=", input_filename_local_2, " chunk number=", chunk_number, "\n"
                        if inp2 != '':
                            inp2 = yap_tools.multiproc_function(
                                preprocess_cmd, inp2, int(
                                    wd.format_specific_lines), '', err_chunk_file, stat_chunk_file)
                            yap_log.merge_multiproc_files(
                                preprocess_cmd,
                                input_filename_local_2,
                                barcode,
                                err_chunk_file,
                                stat_chunk_file)
                        print "Exiting " + preprocess_cmd_name + " : Filename=", input_filename_local_2, " chunk number=", chunk_number, "\n"
        if wd.paired_end_data == 'yes':
            if run_unique_reads == 'True':
		#remove all the umatched pairs from two chunks belonging to the same sample
		#this is because each chunk goes through command separately, not as a pair.
                if inp1 != '' and inp2 != '':
                    inp1, inp2 = yap_tools.find_unique_set(
                        inp1.splitlines(1), inp2.splitlines(1))
	if wd.run_preprocess_analysis  == 'yes':
		#write log data
		str_out="-"*20 + "PREPROCESS FINISHED" + "\t" + str(time.strftime("%Y/%m/%d %H:%M:%S", time.localtime())) + "-"*20 + "\n"
		yap_file_io.write_data(str_out, err_chunk_file)
		yap_file_io.write_data(str_out, stat_chunk_file)
        if wd.data_distribution_method != "file_based":
	    #if the workflow is not filebased; then pass the chunks for alignment.
            if wd.run_reference_alignment == 'yes':
		str_out="-"*20 + "ALIGNMENT STARTED" + "\t" + str(time.strftime("%Y/%m/%d %H:%M:%S", time.localtime())) + "-"*20 + "\n"
		yap_file_io.write_data(str_out, err_chunk_file)
		yap_file_io.write_data(str_out, stat_chunk_file)
                if (wd.paired_end_data == 'yes' and inp1 != '' and inp2 != '') or (wd.paired_end_data != 'yes' and inp1 != ''):
                    print "Entering Alignment: Filename=", input_filename_local, "barcode=", barcode, " chunk number=", chunk_number, "\n"
                    if wd.paired_end_data == 'yes':
                        workflow_prov.append(
                            "INPUT: " +
                            input_filename_local +
                            " and " +
                            input_filename_local_2 +
                            " chunk number= " +
                            str(chunk_number))
                        aligner_out_str, workflow_prov = yap_aligner.run_aligner(
                            inp1, inp2,aligner_output_filename, chunk_number, myrank,workflow_prov, err_chunk_file, stat_chunk_file)
                    else:
                        workflow_prov.append(
                            "INPUT: " +
                            input_filename_local +
                            " chunk number= " +
                            str(chunk_number))
                        aligner_out_str, workflow_prov = yap_aligner.run_aligner(
                            inp1, '', aligner_output_filename, chunk_number,myrank,workflow_prov, err_chunk_file, stat_chunk_file)
                    rm_cmd = "rm " + aligner_output_filename + "*.sai"
                    if len(glob.glob(aligner_output_filename + "*.sai")) > 0:
                        prm = Popen(rm_cmd, shell='False').wait()
                    if len(glob.glob(aligner_output_filename + "*.head")) > 0:
                        prm = Popen(rm_cmd, shell='False').wait()

                else:
                    	print "Exiting Alignment: Filename=", input_filename_local, "barcode=", barcode, " chunk number=", chunk_number, "\n"
		str_out="-"*20 + "ALIGNMENT FINISHED" + "\t" + str(time.strftime("%Y/%m/%d %H:%M:%S", time.localtime())) + "-"*20+ "\n"
		yap_file_io.write_data(str_out, err_chunk_file)
		yap_file_io.write_data(str_out, stat_chunk_file)
            if wd.run_preprocess_analysis == 'yes':
                if wd.write_preprocessed_data == 'yes':
		    #write preprocessed data to file
                    yap_file_io.write_data(inp1, preprocess_out_filename1)
                    if wd.paired_end_data == "yes":
                        yap_file_io.write_data(inp2, preprocess_out_filename2)
                else:
		    #else empty input data chunk
                    inp1 = ''
                    inp2 = ''
            else:
		#else empty input data chunk
                inp1 = ''
                inp2 = ''

        else:
	    #if workflow is filebased; then write preprocessed data to file
            if wd.run_preprocess_analysis == "yes":
                if wd.write_preprocessed_data == 'yes' or wd.run_reference_alignment == "yes":
                    yap_file_io.write_data(inp1, preprocess_out_filename1)
                    if wd.paired_end_data == "yes":
                        yap_file_io.write_data(inp2, preprocess_out_filename2)
        barcode_output_dict[barcode][0] = basecount_matrix_local1
        barcode_output_dict[barcode][1] = basecount_matrix_local2
    return barcode_output_dict, workflow_prov

def execute_file(input_filename_local,input_filename_local_2,file_name,chunk_number,myrank,ii,file_basecount_dict):
        workflow_prov = []
        err_chunk_file = wd.err_log_path + "/" + file_name + \
                     "_log_temp/" + file_name + "_" + str(ii).zfill(6)
        stat_chunk_file = wd.stat_log_path + "/" + file_name + \
                      "_log_temp/" + file_name + "_" + str(ii).zfill(6)
        str_out="*" * 50 + "ALIGNMENT STARTED" + "\t" + str(time.strftime("%Y/%m/%d %H:%M:%S", time.localtime())) + "*" * 50 + "\n"
        yap_file_io.write_data(str_out,err_chunk_file)
        yap_file_io.write_data(str_out,stat_chunk_file)
        for filename_key in file_basecount_dict.iterkeys():
            if filename_key == file_name:
                for barcode in wd.barcode_dict.iterkeys():
                    barcode_value = yap_tools.rename_barcode(barcode)
                    barcode_dir_path = wd.workflow_output_path + "/" + file_name + "/" + barcode
                    aligner_dir_path = barcode_dir_path + "/" + "aligner_output"
                    if wd.alignment_sort_order != 'unsorted':
                        if barcode_value != '':
                            aligner_output_filename = aligner_dir_path + "/" + \
                                "aligner_" + file_name + \
                                "_" + barcode_value
                        else:
                            aligner_output_filename = aligner_dir_path + \
                                "/" + "aligner_" + file_name
                    else:
                        if barcode_value != '':
                            aligner_output_filename = aligner_dir_path + \
                                "/" + file_name + \
                                "_" + barcode_value
                        else:
                            aligner_output_filename = aligner_dir_path + \
                                "/" + file_name

                    if wd.run_preprocess_analysis == 'yes':
                        preprocessed_file_inp1 = ['pipe1']
                        preprocessed_file_inp2 = ['pipe2']
                        preprocess_dir_path = barcode_dir_path + \
                            "/" + "preprocess_output"
                        preprocessed_inp1 = preprocess_dir_path + \
                            "/" + "*preprocess_data*_1.txt"
                        preprocessed_inp2 = preprocess_dir_path + \
                            "/" + "*preprocess_data_*2.txt"
                        preprocessed_file_inp1 = glob.glob(
                            preprocessed_inp1)
                        if wd.paired_end_data == "yes":
                            preprocessed_file_inp2 = glob.glob(
                                preprocessed_inp2)
                        if (wd.paired_end_data== "yes" and preprocessed_file_inp1 and preprocessed_file_inp2) or (wd.paired_end_data != "yes" and preprocessed_file_inp1):
                            print "Entering Alignment section: Filename=", input_filename_local, "barcode=", barcode, "\n"
                            if wd.paired_end_data == 'yes':
                                workflow_prov.append(
                                    "INPUT: " +
                                    preprocessed_file_inp1[0] +
                                    " and " +
                                    preprocessed_file_inp2[0])
                                aligner_out_str, workflow_prov = yap_aligner.run_aligner(preprocessed_file_inp1[0], preprocessed_file_inp2[
                                                                                    0],aligner_output_filename, chunk_number,myrank,workflow_prov, err_chunk_file, stat_chunk_file)
                            else:
                                workflow_prov.append(
                                    "INPUT: " +
                                    preprocessed_file_inp1[0])
                                aligner_out_str, workflow_prov = yap_aligner.run_aligner(preprocessed_file_inp1[
                                                                                    0], '', aligner_output_filename, chunk_number,myrank, workflow_prov, err_chunk_file, stat_chunk_file)

                            if wd.write_preprocessed_data != 'yes':
                                prm1 = Popen(
                                    "rm " +
                                    preprocess_dir_path +
                                    "/" +
                                    "*preprocess_data*_1.txt",
                                    shell='False').wait()
                                if paired_end_data == "yes":
                                    if preprocessed_file_inp2:
                                        prm2 = Popen(
                                            "rm " +
                                            preprocess_dir_path +
                                            "/" +
                                            "*preprocess_data*_2.txt",
                                            shell='False').wait()
                        else:
                            print "Skipping Alignment for : Filename=", input_filename_local, "barcode=", barcode, "........", "No preprocessed data found"
                    else:
                        if wd.paired_end_data == 'yes':
                            workflow_prov.append(
                                "INPUT: " +
                                input_filename_local +
                                " and " +
                                input_filename_local_2)
                            aligner_out_str, workflow_prov = yap_aligner.run_aligner(
                                input_filename_local, input_filename_local_2, aligner_output_filename, 0, workflow_prov, err_chunk_file, stat_chunk_file)
                        else:
                            workflow_prov.append("INPUT: " + input_filename_local)
                            aligner_out_str, workflow_prov = yap_aligner.run_aligner(
                                input_filename_local, '', aligner_cmd_arr, aligner_output_filename, 0, workflow_prov, err_chunk_file, stat_chunk_file)
		    #remove temporary files created by aligners
                    rm_cmd = "rm " + \
                        aligner_output_filename + "*.sai"
                    if len(glob.glob(aligner_output_filename + "*.sai")) > 0:
                        prm = Popen(
                            rm_cmd, shell='False').wait()
                    if barcode in file_basecount_dict[filename_key]:
                        pass
                    else:
                        file_basecount_dict[
                            filename_key][barcode] = []
	#write to log
        str_out="*" * 50 + "ALIGNMENT FINISHED" + "\t" + str(time.strftime("%Y/%m/%d %H:%M:%S", time.localtime())) + "*" * 50 + "\n"
        yap_file_io.write_data(str_out,err_chunk_file)
        yap_file_io.write_data(str_out,stat_chunk_file)
        return workflow_prov, file_basecount_dict
