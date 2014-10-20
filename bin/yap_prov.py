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
from yap_tools import *
import yap_workflow_dict as wd

def preprocess_alignment_prov(
        file_split_struct,
        workflow_prov):
    if wd.run_preprocess_analysis == "yes" or wd.run_reference_alignment == "yes":
        for ii in range(0, len(wd.inp_files_list)):
            fh1 = ''
            fh2 = ''
            input_filename_local = wd.inp_files_list[ii][0]
            input_filename_local_2 = wd.inp_files_list[ii][1]
            path_name, file_name = os.path.split(input_filename_local)
            file_name, extension = os.path.splitext(file_name)
            if wd.run_preprocess_analysis == 'yes':
                for jj in range(0, len(wd.preprocess_cmd_arr)):
                    preprocess_cmd_name = wd.preprocess_cmd_arr[jj][2][0][0]
                    preprocess_cmd = wd.preprocess_cmd_arr[jj][2][0][1]
                    if re.search("fastx_barcode_splitter", preprocess_cmd_name) is not None:
                        fastx_barcode_splitter_cmd = "COMMAND: " + \
                            preprocess_cmd + '; INPUT: ' + input_filename_local
                        workflow_prov.append(fastx_barcode_splitter_cmd)
                    	if wd.paired_end_data == 'yes':
                        	fastx_barcode_splitter_cmd = "COMMAND: " + \
                            	preprocess_cmd + '; INPUT: ' + \
                            	input_filename_local_2
                        	workflow_prov.append(fastx_barcode_splitter_cmd)
                for barcode in wd.barcode_dict.iterkeys():
                    barcode_dir_path = wd.workflow_output_path + "/" + file_name + "/" + barcode + "/" + "preprocess_output"
                    if wd.write_preprocessed_data == 'yes':
                        preprocessed_output_basename_1 = barcode_dir_path + "/" + \
                            "preprocess_data" + "_" + file_name + \
                            "_" + barcode + "_1.txt"
                        preprocessed_output_basename_2 = barcode_dir_path + "/" + \
                            "preprocess_data" + "_" + file_name + \
                            "_" + barcode + "_2.txt"
                    else:
                        preprocessed_output_basename_1 = 'No output file as \"write_preprocessed_data\" was set as no'
                        preprocessed_output_basename_2 = 'No output file as \"write_preprocessed_data\" was set as no'

                    basecount_output_text_1 = barcode_dir_path + '/' + \
                        file_name + "_" + barcode + \
                        "basecountmetrics" + "_1.txt"
                    basecount_output_eps_1 = barcode_dir_path + '/' + \
                        file_name + "_" + barcode + \
                        "basecountmetrics" + "_1.eps"
                    basecount_output_text_2 = barcode_dir_path + '/' + \
                        file_name + "_" + barcode + \
                        "basecountmetrics" + "_2.txt"
                    basecount_output_eps_2 = barcode_dir_path + '/' + \
                        file_name + "_" + barcode + \
                        "basecountmetrics" + "_2.eps"
                    for jj in range(0, len(wd.preprocess_cmd_arr)):
                        preprocess_cmd_name = wd.preprocess_cmd_arr[jj][2][0][0]
                        preprocess_cmd = wd.preprocess_cmd_arr[jj][2][0][1]
                        if re.search("fastx_barcode_splitter",preprocess_cmd_name) or \
                        	re.search("fastqc",preprocess_cmd_name) or \
                                re.search("fastq_screen",preprocess_cmd_name):
                        	pass
                        else:
                            if re.search("fastx_clipper", preprocess_cmd_name) is not None:
                                if input_filename_local in wd.contaminant_dict.keys():
                                    contaminants_arr1 = wd.contaminant_dict[
                                        input_filename_local]
                                    index = 0
                                    for index in range(0, len(contaminants_arr1)):
                                        fastx_clipper_cmd = preprocess_cmd
                                        contaminant1 = contaminants_arr1[
                                            index].strip("\n")
                                        fastx_clipper_cmd = 'COMMAND: ' + fastx_clipper_cmd.replace(
                                            '-i -',
                                            " -a " + contaminant1 + " -i - ") + '\tINPUT: ' + input_filename_local + '; BARCODE: ' + barcode + '; OUTPUT: ' + preprocessed_output_basename_1
                                        workflow_prov.append(fastx_clipper_cmd)
                                        if wd.paired_end_data == 'yes':
                                            if os.path.isfile(preprocessed_output_basename_2):
                                                if input_filename_local_2 in wd.contaminant_dict.keys():
                                                    contaminants_arr2 = wd.contaminant_dict[
                                                        input_filename_local_2]
                                                    index = 0
                                                    for index in range(0, len(contaminants_arr2)):
                                                        fastx_clipper_cmd = preprocess_cmd
                                                        contaminant2 = contaminants_arr2[
                                                            index].strip("\n")
                                                        fastx_clipper_cmd = 'COMMAND: ' + fastx_clipper_cmd.replace(
                                                            '-i -',
                                                            " -a " + contaminant2 + " -i - ") + input_filename_local_2 + '; BARCODE: ' + barcode + '; OUTPUT: ' + preprocessed_output_basename_2
                                                        workflow_prov.append(
                                                            fastx_clipper_cmd)
                            elif re.search("fastq_quality_trimmer", preprocess_cmd_name) is not None:
                                fastq_quality_trimmer_cmd = 'COMMAND: ' + preprocess_cmd + '\tINPUT: ' + \
                                    input_filename_local + '; BARCODE: ' + barcode + \
                                    '; OUTPUT: ' + preprocessed_output_basename_1
                                workflow_prov.append(fastq_quality_trimmer_cmd)
                                if wd.paired_end_data == 'yes':
                                    fastq_quality_trimmer_cmd = 'COMMAND: ' + preprocess_cmd + '\tINPUT: ' + \
                                        input_filename_local_2 + '; BARCODE: ' + barcode + \
                                        '; OUTPUT: ' + \
                                        preprocessed_output_basename_2
                                    workflow_prov.append(fastq_quality_trimmer_cmd)
                            elif re.search("fastq_quality_filter", preprocess_cmd_name) is not None:
                                fastq_quality_filter_cmd = 'COMMAND: ' + preprocess_cmd + '\tINPUT: ' + \
                                    input_filename_local + '; BARCODE: ' + barcode + \
                                    '; OUTPUT: ' + preprocessed_output_basename_1
                                workflow_prov.append(fastq_quality_filter_cmd)
                                if wd.paired_end_data == 'yes':
                                    fastq_quality_filter_cmd = 'COMMAND: ' + preprocess_cmd + '\tINPUT: ' + \
                                        input_filename_local_2 + '; BARCODE: ' + barcode + \
                                        '; OUTPUT: ' + \
                                        preprocessed_output_basename_2
                                    workflow_prov.append(fastq_quality_filter_cmd)
                            elif re.search("calculate_basecount_metrics", preprocess_cmd_name) is not None:
                                base_count_metrics_cmd = 'BASE METRICS CALCULATED \tINPUT: ' + input_filename_local + \
                                    '; BARCODE: ' + barcode + '; OUTPUT: ' + \
                                    basecount_output_text_1 + \
                                    ', ' + basecount_output_eps_1
                                workflow_prov.append(base_count_metrics_cmd)
                                if wd.paired_end_data == 'yes':
                                    base_count_metrics_cmd = 'BASE METRICS CALCULATED\tINPUT: ' + input_filename_local_2 + \
                                        '; BARCODE: ' + barcode + '; OUTPUT: ' + \
                                        basecount_output_text_2 + \
                                        ', ' + basecount_output_eps_2
                                    workflow_prov.append(base_count_metrics_cmd)
                            else:
                                prov_cmd = 'COMMAND: ' + preprocess_cmd + '\tINPUT: ' + input_filename_local + \
                                    '; BARCODE: ' + barcode + '; OUTPUT: ' + \
                                    preprocessed_output_basename_1
                                workflow_prov.append(preprocess_cmd)
                                if wd.paired_end_data == 'yes':
                                    prov_cmd = 'COMMAND: ' + preprocess_cmd + '\tINPUT: ' + input_filename_local_2 + \
                                        '; BARCODE: ' + barcode + '; OUTPUT: ' + \
                                        preprocessed_output_basename_2
                                    workflow_prov.append(preprocess_cmd)
    return workflow_prov
