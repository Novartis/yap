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

import re
import os
import glob
from multiprocessing import *
import time
from subprocess import PIPE, Popen
import numpy
import yap_tools
import yap_file_io 
import yap_log 
import yap_workflow_dict as wd

def read_barcodes(bar_file):
    '''Reads barcode file and generates dictionary'''
    bar_arr = yap_file_io.read_file(bar_file)
    bar1 = []
    barcode_dict = {}
    for b in range(len(bar_arr)):
        matchobj = re.search(
            '\s*(\w*)(\s*).*', bar_arr[b].strip('\n'), re.M | re.I)
        if matchobj:
            split_by = matchobj.group(2)
        barcode = bar_arr[b].split(split_by)[1]
        bar_id = bar_arr[b].split(split_by)[0]
        barcode_dict[bar_id] = barcode
    barcode_dict['unmatched'] = 'barcode_unmatched'
    return barcode_dict

def write_basecount_matrix(file_basecount_dict):
    '''
    Prepare output filename information and write base count matrix 
    results to a file
    '''
    def write_final_matrix(
            file_basecount_matrix,
            basecount_metrics_filename,
            basecount_file_basename,
            ):
        '''Writes basecount matrix to file'''
        file_w_handler = open(basecount_metrics_filename, 'a+')
        len_matrix = len(file_basecount_matrix)
        len_row = len(file_basecount_matrix[0])
        file_w_handler.write("# A C T G N \n")
        for i in range(len_matrix):
            row_sum = numpy.sum(file_basecount_matrix[i])
            if row_sum != 0:
                file_w_handler.write(str(i) + " ")
                for j in range(len_row):
                    file_w_handler.write(
                        str(file_basecount_matrix[i][j]) + " ")
                file_w_handler.write("\n")
        file_w_handler.close()
        #create plot 
        os.system("gplot.sh " + basecount_metrics_filename +
                  " " + basecount_file_basename + '.eps')
    #For each sample, prepare basecount output filename information 
    for filename_key in file_basecount_dict.iterkeys():
        path_name, file_name = os.path.split(filename_key)
        file_name, extension = os.path.splitext(file_name)
        barcode_basecount_dict = file_basecount_dict[filename_key]
        for barcode in barcode_basecount_dict.iterkeys():
            barcode_value = yap_tools.rename_barcode(barcode)
            barcode_dir_path = wd.workflow_output_path + "/" + file_name + "/" + barcode
            preprocess_output_dir = barcode_dir_path + "/" + "preprocess_output"
            file_basecount_matrix1 = barcode_basecount_dict[barcode][0]
            file_basecount_matrix2 = barcode_basecount_dict[barcode][1]
            if barcode_value != '':
                basecount_file_basename1 = preprocess_output_dir + '/' + \
                    file_name + "_" + barcode_value + "_basecountmetrics_1"
                basecount_metrics_filename1 = preprocess_output_dir + '/' + \
                    file_name + "_" + barcode_value + \
                    "_basecountmetrics_1" + '.txt'
                basecount_file_basename2 = preprocess_output_dir + '/' + \
                    file_name + "_" + barcode_value + "_basecountmetrics_2"
                basecount_metrics_filename2 = preprocess_output_dir + '/' + \
                    file_name + "_" + barcode_value + \
                    "_basecountmetrics_2" + '.txt'
            else:
                basecount_file_basename1 = preprocess_output_dir + \
                    '/' + file_name + "_basecountmetrics_1"
                basecount_metrics_filename1 = preprocess_output_dir + \
                    '/' + file_name + "_basecountmetrics_1" + '.txt'
                basecount_file_basename2 = preprocess_output_dir + \
                    '/' + file_name + "_basecountmetrics_2"
                basecount_metrics_filename2 = preprocess_output_dir + \
                    '/' + file_name + "_basecountmetrics_2" + '.txt'

            if (numpy.sum(file_basecount_matrix1) > 0):
                print  "For filename = ", file_name, "barcode = ", barcode, " ...producing the basecount metrics results for first paired file", "\n"
                #pass output filename informaton to write function 
                write_final_matrix(
                    file_basecount_matrix1,
                    basecount_metrics_filename1,
                    basecount_file_basename1,
                    )
            else:
                print "No data for finename = ", file_name, " barcode = ", barcode, " ...skipping the basecount metrics output for first paired file", "\n"
            if (numpy.sum(file_basecount_matrix2) > 0):

                print "For filename = ", file_name, "barcode = ", barcode, " ...producing the basecount metrics results second paired file", "\n"
                #pass output filename informaton to write function 
                write_final_matrix(
                    file_basecount_matrix2,
                    basecount_metrics_filename2,
                    basecount_file_basename2,
                    )
            else:
                print "No data for finename = ", file_name, " barcode = ", barcode, " ...skipping the basecount metrics output for second paired file", "\n"

def merge_preprocessed_output(file_basecount_dict):
    '''For Every sample, merge preprocessed chunk data into a single file'''
    for filename_key in file_basecount_dict.iterkeys():
        #iterate over filename and barcode, get list of files to be merged
        file_name = filename_key  
        barcode_basecount_dict = file_basecount_dict[filename_key]
        for barcode in barcode_basecount_dict.iterkeys():
            barcode_value = yap_tools.rename_barcode(barcode)
            barcode_dir_path = wd.workflow_output_path + "/" + file_name + "/" + barcode + "/" + "preprocess_output"
            if barcode_value != '':
                preprocessed_output_basename = barcode_dir_path + "/" + "preprocess_data" + "_" + file_name + "_" + barcode_value
            else:
                preprocessed_output_basename = barcode_dir_path + "/" + "preprocess_data" + "_" + file_name
	    barcode_value = yap_tools.rename_barcode(barcode)
            preprocess_content1 = glob.glob(
            barcode_dir_path + "/" + "*preprocessed_data*_1.txt")
	    # cat files together 
            if len(preprocess_content1) > 0:
            	os.system("cat " + barcode_dir_path +"/" +"*preprocessed_data*_1.txt" +">" + preprocessed_output_basename + "_1.txt")
            	os.system("rm " + barcode_dir_path + "/" + barcode_value + "*preprocessed_data*_1.txt")
            if wd.paired_end_data == "yes":
            	preprocess_content2 = glob.glob(barcode_dir_path + "/" + "*preprocessed_data*_2.txt")
            	if len(preprocess_content2) > 0:
                	os.system("cat " + barcode_dir_path + "/" + "*preprocessed_data*_2.txt" + ">" + preprocessed_output_basename + "_2.txt")
                	os.system("rm " + barcode_dir_path + "/" + barcode_value + "*preprocessed_data*_2.txt")

def fastx_barcode_splitter(
        seqs_str,
        output_file_format,
        fastx_barcode_splitter_cmd,
        preprocess_prov,
        err_log,
        stat_log):
    '''
    Runs barcode splitter, returns output data into a dictionary,
    where key represents a barcdoe and sequence string as value
    '''
    bar_seq = ''
    preprocess_prov.append(fastx_barcode_splitter_cmd)
    globaldict = {}
    P1 = Popen(fastx_barcode_splitter_cmd, stdin=PIPE,
               stdout=PIPE, stderr=PIPE, shell=True)
    try:
        std_out, std_err = P1.communicate(seqs_str)
        exit_code = P1.returncode
        yap_log.write_log(fastx_barcode_splitter_cmd, "",
                  exit_code, std_err, err_log, stat_log)
        bar_seq_split = std_out.replace(" ", "").split("|")
        for i in range(0, len(bar_seq_split)):
            if bar_seq_split[i] != '':
                splited_S = bar_seq_split[i].split("=>")
                globaldict[splited_S[0]] = splited_S[1]
                del splited_S
    except Exception as e:
        write_data(str(e), err_log)
    yap_file_io.write_data("\n", err_log)
    yap_file_io.write_data("\n", stat_log)
    return globaldict, preprocess_prov


def run_fastqc(inp_files_list,fastqc_cmd):
    ''' 
    Runs fastqc command, writes log information to the files and returns log data list
    '''
    prov = []
    file_base_name = inp_files_list[2]
    fastqc_cmd += inp_files_list[0] + " " + inp_files_list[1] + " "
    err_log = wd.err_log_path + "/" + file_base_name + "_fastqc_err.log"
    stat_log = wd.stat_log_path + "/" + file_base_name + "_fastqc_stat.log"
    fastqc_cmd = fastqc_cmd.replace('output_directory',wd.workflow_output_path + "/" + file_base_name + "/" + "no_barcode_specified" + "/" + "preprocess_output")
    fastqc_cmd = fastqc_cmd.replace('pipe1', '')
    prov.append(fastqc_cmd)
    str_out="*" * 50 + "FASTQC STARTED" + "\t" + str(time.strftime("%Y/%m/%d %H:%M:%S", time.localtime())) + "*" * 50 + "\n"
    yap_file_io.write_data(str_out,err_log)
    yap_file_io.write_data(str_out,stat_log)
    prm = Popen(fastqc_cmd, stderr=PIPE, shell='False')
    std_out, std_err = prm.communicate()
    exit_code = prm.returncode
    if exit_code != 0:
        print "Error encountered during FastQC analysis", "\n"
    yap_log.write_log(
        fastqc_cmd, file_base_name, exit_code, std_err, err_log, stat_log)
    str_out="*" * 50 + "FASTQC FINISHED" + "\t" + str(time.strftime("%Y/%m/%d %H:%M:%S", time.localtime())) + "*" * 50 + "\n"
    yap_file_io.write_data(str_out,err_log)
    yap_file_io.write_data(str_out,stat_log)
    return prov

def run_fastq_screen(inp_files_list, fastq_screen_cmd):
    ''' 
    Runs fastq screen command, writes log information to the files and returns log data list
    '''
    prov = []
    file_base_name = inp_files_list[2]
    err_log = wd.err_log_path + "/" + file_base_name + "_fastqscreen_err.log"
    stat_log =wd.stat_log_path + "/" + file_base_name + "_fastqscreen_stat.log"
    fastq_screen_cmd += inp_files_list[0] + " " + inp_files_list[1] + " "
    fastq_screen_cmd = fastq_screen_cmd.replace('output_directory',wd.workflow_output_path + "/" + file_base_name + "/" + "no_barcode_specified" + "/" + "preprocess_output")
    fastq_screen_cmd = fastq_screen_cmd.replace('pipe1', '')
    fastq_screen_cmd += " "
    str_out="*" * 50 + "FASTQSCREEN STARTED" + "\t" + str(time.strftime("%Y/%m/%d %H:%M:%S", time.localtime())) + "*" * 50 + "\n"
    yap_file_io.write_data(str_out,err_log)
    yap_file_io.write_data(str_out,stat_log)
    prm = Popen(fastq_screen_cmd, stderr=PIPE, shell='False')
    std_out, std_err = prm.communicate()
    exit_code = prm.returncode
    prov.append(fastq_screen_cmd)
    yap_log.write_log(fastq_screen_cmd, file_base_name,
              exit_code, std_err, err_log, stat_log)
    str_out="*" * 50 + "FASTQSCREEN FINISHED" + "\t" + str(time.strftime("%Y/%m/%d %H:%M:%S", time.localtime())) + "*" * 50 + "\n"
    yap_file_io.write_data(str_out,err_log)
    yap_file_io.write_data(str_out,stat_log)
    return prov
