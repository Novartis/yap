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
import numpy
import time
import glob
import sys
import math
import re
import random
from yap_tools import *
from yap_file_io import *
from yap_preprocess import *
from yap_conflict_check import *
from yap_regroup_check import *
from yap_cuffdiff_check import *
from yap_cuff_merge_compare_check import *
from yap_printer import *
from yap_path_checks import *

def workflow_validator(workflow_struct, workflow_config_file, run_mode):
    """
    Given 3 parameters : Workflow Struct, configuration file and run mode.
    If given workflow_struct is empty, which means no correct workflow is provided, give error information and return.
    :param workflow_struct:
    :param workflow_config_file:
    :param run_mode:
    :return:
    """
    if not workflow_struct:
        no_workflow_error = "Error: No correct workflow is provided."
        print no_workflow_error
        return [no_workflow_error]

    # output variable
    results = []

    # Iterate over all workflows provided
    for wk in range(1, len(workflow_struct)):
        # Parameters Declaration
        inp_files = []
        aligner_arr = []
        output_dict = {}
        file_basecount_dict = {}
        file_split_struct = []
        cmd_str_arr = []
        postprocess_cmd_arr = []
        preprocess_cmd_arr = []
        fastqc_job_status = ''
        contaminant_dict = {}
        barcode_dict = {}
        fastqc_file_split_arr = []
        status_check = "healthy"
        status_arr = []
        preprocess_errorlist = []
        aligner_errorlist = []
        postprocess_errorlist = []
        parser_errorlist = []
        config_mismatch = []
        warning_list = []
        incompatibility_error_list = []
        cont_check_error_list = []
        conflict_validator = None
        barcode_dict = {'no_barcode_specified': ''}
        inp_files_list = []
        aligner_cmd_arr = []
        eqp_warning = ''
        # Set default values for workflow
        workflow_config_dict = workflow_struct[wk]
	YAP_HOME= str(os.getenv("YAP_HOME"))
        yap_src = YAP_HOME + "/bin"
        yap_version = "YAP-1.0"
        workflow_config_dict["yap_version"]=yap_version
        workflow_config_dict["yap_src"]=yap_src
        workflow_config_dict["usage_log"]="/usr/prog/yap/usage/usage.log"
        if ("max_read_length" in workflow_config_dict) != True:
            workflow_config_dict["max_read_length"] = "150"
        if ("file_chunk_size" in workflow_config_dict) != True:
            workflow_config_dict["file_chunk_size"] = "1024"
        workflow_config_dict["output_files_path"] = \
            os.path.abspath(workflow_config_dict["output_files_path"])
        username = os.getenv("USER")
        workflow_config_dict["username"] = username
        unique_jobid = ''
        if 'JOB_ID' in os.environ.keys():
            unique_jobid = os.environ['JOB_ID']
        workflow_config_dict["job_id"] = unique_jobid
        workflow_config_dict["random_id"] = str(random.random())
        workflow_config_dict["workflow_output_path"] = workflow_config_dict["output_files_path"] + \
                                                       "/" + workflow_config_dict["workflow_name"]
        if os.path.exists(workflow_config_dict["workflow_output_path"]):
            date_time = time.strftime("%Y.%m.%d-%H:%M:%S", time.localtime())
            workflow_config_dict["workflow_output_path"] = workflow_config_dict["workflow_output_path"] + "_" \
                                                           +date_time
            workflow_config_dict["workflow_name"] = workflow_config_dict["workflow_name"] + "_" + date_time
        workflow_config_dict["consolidated_output"] = workflow_config_dict["workflow_output_path"] + "/" + "consolidated_output"
        #turn off regroup output if alignment stage not in the workflow
	#regroup of output is valid only for alignment output
    	if  workflow_config_dict["run_reference_alignment"] == "no":
        	workflow_config_dict["regroup_output"] = "no"
        workflow_config_dict["consolidated_output"] = workflow_config_dict["workflow_output_path"] + "/" + "consolidated_output"
        temp_dir_path = workflow_config_dict["workflow_output_path"]+ "/" + "temp"
        workflow_config_dict["temp_dir_path"] = temp_dir_path
        workflow_config_dict["sample_track_log"] = workflow_config_dict["workflow_output_path"] + "/" + "sample_name_track.log"
        workflow_config_dict["log_path"] = workflow_config_dict["workflow_output_path"] + "/" + "yap_log_files"
        err_log_path = workflow_config_dict["log_path"] + "/error_log"
        stat_log_path = workflow_config_dict["log_path"] + "/status_log"
        workflow_config_dict["err_log_path"] = err_log_path
        workflow_config_dict["stat_log_path"] = stat_log_path
        workflow_config_dict["aligner_output_key_arr"] = []
        # Set chunk size
        # read main configuration file  and create dictionary
        # create preprocess configuration dictionary based on user input
        if workflow_config_dict["run_preprocess_analysis"] == "yes":
            preprocess_config_file = \
                workflow_config_dict["preprocess_configuration_file"]
            # proceed only if a valid pre-process config file is provided
            if os.path.exists(preprocess_config_file):
                preprocess_file_arr = yap_file_io.read_file(preprocess_config_file)
                preprocess_cmd_arr, preprocess_errorlist = yap_tools.command_parser(
                    preprocess_file_arr, preprocess_config_file)
            else:  # record error in output if preprocess config file
                # is not given
                status_check = 'unhealthy'
                preprocess_errorlist.append(
                    "Error:  while Opening the file : " +
                    preprocess_config_file +
                    "\nI/O error: No such file or directory")
        # setting up the variables for eqp ( support for in-house functionality, neglect otherwise)
        eqp_rename_flag= ''
        for jj in range(0,len(preprocess_cmd_arr)):
                preprocess_cmd_name = preprocess_cmd_arr[jj][2][0][0]
                preprocess_cmd = preprocess_cmd_arr[jj][2][0][1]
                if re.search('eqp_rename_reads',preprocess_cmd_name) == None:
                        eqp_rename_flag= 'False'
        if workflow_config_dict["workflow_type"].lower().find("eqp") == 0:
               	#if workflow_config_dict["data_distribution_method"] == "file_based":
                #	workflow_config_dict["data_distribution_method"] = "chunk_based"
                #      	eqp_warning = 'Warning : As workflow type is set to EQP, chunk_based data distribution method' +\
                #       	'is used for this workflow.'
                eqp_home=''
                eqp_home= os.getenv("EQP_HOME")
                if os.path.exists(str(eqp_home) + "/java") != True:
                       	status_check='unhealthy'
                       	preprocess_errorlist.append("Error: EQP java source does not exist. "  + "Please configure EQP environment")
		if workflow_config_dict["run_reference_alignment"] != "no" and workflow_config_dict["run_postprocess_analysis"] != "yes":
                	workflow_config_dict["run_preprocess_analysis"] = "yes"
                	if len(preprocess_cmd_arr) == 0:
                        	if len(eqp_warning) > 0:
                                	eqp_warning += '\n'

                        	eqp_warning = eqp_warning + 'Warning : As workflow type is set to EQP, YAP turned on the '+\
                        	'preprocess analysis and set the command : eqp_rename_reads to yes.'

                	#turn eqp on
                	if eqp_rename_flag == 'False' or len(preprocess_cmd_arr) == 0:
                    		preprocess_cmd_arr.append([':begin', ['input_directory aligner_output'], [['eqp_rename_reads','eqp_rename_reads']]])
                    		eqp_rename_flag == 'True'
			
        	if len(preprocess_cmd_arr) != 0:
            		if eqp_rename_flag == "True":
                		if workflow_config_dict["data_distribution_method"] == "chunk_based":
                    			if workflow_config_dict["merge_alignment_output"] == "yes":
                        			workflow_config_dict["merge_alignment_output"] = "no"
                        		eqp_warning = eqp_warning + 'Warning : As eqp_rename_reads preprocess command is turned on in chunk based mode, YAP is turning off the the '+\
                                		'merge_alignment_output'
        # check if reference alignment is required
        if workflow_config_dict["run_reference_alignment"] == "yes":
            # read aligner config file
            aligner_config_file = workflow_config_dict[
                "aligner_configuration_file"]
            if os.path.exists(aligner_config_file):
                # reading aligner configuration file into array
                aligner_arr = yap_file_io.read_file(aligner_config_file)
                # DOUBT: Code for command parser??
                # entering command parser to generate alignment commands
                aligner_cmd_arr, aligner_errorlist = yap_tools.command_parser(
                    aligner_arr, aligner_config_file)
                alignment_sort_options = [
                    'coordinate', 'queryname', 'both', 'unsorted', 'no']
                if workflow_config_dict["alignment_sort_order"] not in alignment_sort_options:
                    status_check = "unhealthy"
                    aligner_errorlist.append(
                        "Error: Invalid sort order for variable \"alignment_sort_order\". Entered value: " +
                        workflow_config_dict["alignment_sort_order"])
            else:  # record error if aligner config is not found
                status_check = 'unhealthy'
                aligner_errorlist.append(
                    "Error:  while Opening the file : " +
                    aligner_config_file +
                    "\nI/O error: No such file or directory")
        # read postprocess configuration if it is configured to run
        if workflow_config_dict["run_postprocess_analysis"] == "yes":
            postprocess_config_file = workflow_config_dict[
                "postprocess_configuration_file"]
            if os.path.exists(postprocess_config_file):
                postprocess_file_arr = yap_file_io.read_file(postprocess_config_file)
                postprocess_cmd_arr, postprocess_errorlist = yap_tools.command_parser(
                    postprocess_file_arr, postprocess_config_file)
                # create a new object of yap_cmd_check
                postprocess_cmd_check = yap_cmd_checks.yap_cmd_checks()
                if workflow_config_dict["run_reference_alignment"] == "yes":
                    # if run alignment, then check if needed sort
                    # in postprocess matches the alignment_sort_order
                    # which is set in workflow config file.
                    postprocess_cmd_check.check_sorted(
                        postprocess_cmd_arr,
                        workflow_config_dict["alignment_sort_order"])
                    if postprocess_cmd_check.warning_list:
                        # if postprocess_cmd_check.warning_list is not empty
                        # extend it to the global warning_list
                        warning_list.extend(postprocess_cmd_check.warning_list)
                else:
                    # if not run alignment, check if needed sorted
                    # files exist in 'aligner_output'
                    postprocess_cmd_check.check_files_when_postproc_only(
                        postprocess_cmd_arr,
                        workflow_config_dict['input_files_path'])
                    if postprocess_cmd_check.warning_list:
                        # if postprocess_cmd_check.warning_list is not empty
                        # extend it to the global warning_list
                        warning_list.extend(postprocess_cmd_check.warning_list)
                incompatibility_error_list = postprocess_cmd_check.error_list
                if incompatibility_error_list:
                    incompatibility_error_list.insert(
                        0, "workflow config file: " + workflow_config_file)
                    incompatibility_error_list.insert(
                        0,
                        "postprocess config file: " +
                        postprocess_config_file)
            else:
                status_check = 'unhealthy'
                postprocess_errorlist.append(
                    "Error:  while Opening the file : " +
                    postprocess_config_file +
                    "\nI/O error: No such file or directory")
        parser_errorlist = preprocess_errorlist + \
            aligner_errorlist + postprocess_errorlist
        # check for environment variable
        local_temp_dir = ''
        local_temp_dir = os.getenv("YAP_LOCAL_TEMPDIR")
        workflow_config_dict["yap_temp_source"]=local_temp_dir
        workflow_config_dict["yap_temp_user_dir"] = local_temp_dir + '/' + username + '/yap_temp'
        if local_temp_dir is None:
            status_check = "unhealthy"
            parser_errorlist.append(
                "Error: Variable YAP_LOCAL_TEMPDIR is not defined." +
                " Please define the variable for YAP to work properly.")
        else:
            if os.access(local_temp_dir, os.W_OK) == False:
                status_check = "unhealthy"
                parser_errorlist.append(
                    "Error: Write permission is not granted to YAP_LOCAL_TEMPDIR =" +
                    local_temp_dir +
                    ". Please define the variable for YAP to work properly.")
        # determining the files to be processed
        # get input files path and format
        input_file_path = workflow_config_dict["input_files_path"].split(";")
        input_files = workflow_config_dict["input_files"].split(";")
        file_format = workflow_config_dict["input_file_format"]
        inp_files_list = []
        # convert each input file to absolute path and add it back to
        # the workflow configuration variable
        for i in range(len(input_file_path)):
            if input_file_path[i] != '':
                input_file_path[i] = input_file_path[i].strip()
                input_file_path[i] = os.path.abspath(input_file_path[i])
        workflow_config_dict["input_files_path"] = ";".join(input_file_path)
        ''' user needs to specify same number of input paths and input
		files . Its an error if its a mismatch '''
        if len(input_file_path) == len(input_files):
            for i in range(0, len(input_file_path)):
                input_file_path[i] = input_file_path[i].strip()
                input_files[i] = input_files[i].strip()
                input_files_arr = input_files[i].split(",")
                if len(input_files_arr) > 1:
                    for inf in (input_files_arr):
                        if inf == '':
                            # DOUBT : Why is this
                            # not an error??
                            pass
                        else:
                            inf = inf.strip()
                            inf = input_file_path[i] + "/" + inf
                            if glob.glob(inf):
                                inp_files += glob.glob(inf)
                            else:
                                status_check = "unhealthy"
                                parser_errorlist.append(
                                    "Error: No input file " + inf)
                else:  # DOUBT: what's going on here?
                    if glob.glob(input_file_path[i] + "/" + input_files[i]):
                        inp_files = inp_files + \
                            glob.glob(
                                input_file_path[i] + "/" + input_files[i])
                    else:
                        status_check = "unhealthy"
                        parser_errorlist.append(
                            "Error: No input file: " +
                            input_file_path[i] +
                            "/" +
                            input_files[i])
            # DOUBT: why is sorting required??
            inp_files.sort()
            nfiles = len(inp_files)
            if nfiles == 0:
                status_check = "unhealthy"
                parser_errorlist.append("Error: No input files found")
                parser_errorlist.append(
                    "Please verify the file path and file name")
            # check for symbolic links in input files
            sym_link_checker = yap_path_checks()
            sym_link_check_results = filter(
                sym_link_checker.invalid_path_or_link, inp_files)
            if len(sym_link_check_results) > 0:
                status_check = "unhealthy"
                msg_pre = "Error: input file with broken symbolic link found: "
                msg_post = ". Please make sure that actual file exists and " +\
                           "required permissions are granted to the files."
                map((lambda x: parser_errorlist.append(
                    msg_pre + x + msg_post)), sym_link_check_results)
        else:
            status_check = "unhealthy"
            parser_errorlist.append(
                "Please see to it that each input path corresponds to an input file and its respective input format "
                "and separate them with a ';'\nFor eg:\n\t\t\"input_file_path\" := \"input_path1;input_path2\"\n\t"
                "\t\"input_files\" := \"file1.1,file1.2;file2\"")
        # Read configuration related to file formats
        if file_format == 'fasta':
            workflow_config_dict["format_specific_lines"] = 2
        elif file_format == 'fastq':
            workflow_config_dict["format_specific_lines"] = 4
        elif file_format == 'qseq':
            if workflow_config_dict["preprocess_output_file_format"] == 'fastq':
                workflow_config_dict["format_specific_lines"] = 4
            if workflow_config_dict["preprocess_output_file_format"] == 'fasta':
                workflow_config_dict["format_specific_lines"] = 2
            if workflow_config_dict["preprocess_output_file_format"] == 'sam':
                workflow_config_dict["format_specific_lines"] = 1
            if workflow_config_dict["preprocess_output_file_format"] == 'bam':
                workflow_config_dict["format_specific_lines"] = 1
            if workflow_config_dict["preprocess_output_file_format"] == 'tab':
                workflow_config_dict["format_specific_lines"] = 1
        elif file_format == 'tab':
            workflow_config_dict["format_specific_lines"] = 1
        elif file_format == 'sam':
            workflow_config_dict["format_specific_lines"] = 1
        elif file_format == 'bam':
            workflow_config_dict["format_specific_lines"] = 1
        elif file_format == 'bed':
            workflow_config_dict["format_specific_lines"] = 1
        else:
            status_check = "unhealthy"
            parser_errorlist.append ("Error: Unknown input file + \
                        			format mentioned :  " + file_format)
            parser_errorlist .append(
                "please specify correct file format [input file format supported are:  qseq ,fastq ,fasta,tab]")
        # read paired end data related configuration
        if workflow_config_dict["paired_end_data"] == "yes":
            nfiles = len(inp_files)
            if workflow_config_dict["run_preprocess_analysis"] == "yes" or workflow_config_dict["run_reference_alignment"] == "yes":
                if nfiles % 2 != 0:
                    status_check = "unhealthy"
                    parser_errorlist.append(
                        "Error: The files are lacking paired files")
                    parser_errorlist.append(
                        "please check the file names. The paired file names should only differ by \'1\' and \'2\' "
                        "in their naming convention. eg fasta_1.txt fasta_2.txt")
                # get paried end files list
                paired_list = []
                unpaired_list = []
                paired_list, unpaired_list = yap_tools.find_pair_end_files(inp_files)
                if len(unpaired_list) > 0:
                    status_check = "unhealthy"
                    parser_errorlist.append(
                        "Error: Files are lacking paired files")
                    parser_errorlist.append("Unpaired files list ")
                    parser_errorlist += unpaired_list
                inp_files_list = paired_list
            else:
                for i in range(0, nfiles):
                    inp_files_list.append([inp_files[i], ""])
        else:
            for i in range(0, nfiles):
                inp_files_list.append([inp_files[i], ""])
        # Check for replicate file names and apply the tags
        file_count_dict = {}
        inp_files_dict = {}
        n_inp_files_list = []
        extension_arr = [
            '.fq', '.fastq', '.bam', '.sam', '.fasta', '.fa', '.gtf', '.bed']
        for i in range(0, len(inp_files_list)):
            file1 = ''
            ext1 = ''
            file2 = ''
            ext2 = ''
            final1_ext = ''
            final2_ext = ''
            path, file1 = os.path.split(inp_files_list[i][0])
            file1, ext1 = os.path.splitext(file1)
            final1_ext += ext1
            file_n, ext1 = os.path.splitext(file1)
            if ext1 in extension_arr:
                file1, ext1 = os.path.splitext(file1)
                final1_ext = ext1 + final1_ext

            path, file2 = os.path.split(inp_files_list[i][1])
            file2, ext2 = os.path.splitext(file2)
            final2_ext += ext2
            file_n, ext2 = os.path.splitext(file2)
            if ext2 in extension_arr:
                file2, ext2 = os.path.splitext(file2)
                final2_ext = ext2 + final2_ext

            if (file1 in file_count_dict):
                tag = "_file" + str(file_count_dict[file1]).zfill(6)
                file1_tag = file1 + tag
                file2_tag = file2 + tag
                inp_files_dict[file1_tag] = [
                    inp_files_list[i][0],
                    inp_files_list[i][1],
                    file1_tag,
                    file2_tag]
                n_inp_files_list.append(
                    [inp_files_list[i][0], inp_files_list[i][1], file1_tag, file2_tag])
            else:
                file_count_dict[file1] = 1
                inp_files_dict[file1] = [
                    inp_files_list[i][0], inp_files_list[i][1], '', '']
                n_inp_files_list.append(
                    [inp_files_list[i][0], inp_files_list[i][1], file1, file2])

        # reading barcode file and creating barcode dictionary
        # DOUBT : what needs to be done with bar code??
        if workflow_config_dict["run_preprocess_analysis"] == "yes":
            barcode_splitter_index = ''
            fastqc_index = ''
            fastq_screen_index = ''
            for jj in range(0, len(preprocess_cmd_arr)):
                preprocess_cmd_type = preprocess_cmd_arr[jj][0]
                preprocess_meta_data = preprocess_cmd_arr[jj][1]
                preprocess_temp_arr = preprocess_cmd_arr[jj][2]
                preprocess_cmd_name = preprocess_temp_arr[0][0]
                preprocess_cmd = preprocess_temp_arr[0][1]
                if re.search('fastqc', preprocess_cmd_name) is not None:
                    fastqc_index = jj
                if re.search('fastq_screen', preprocess_cmd_name) is not None:
                    fastq_screen_index = jj
                if re.search('fastx_clipper', preprocess_cmd_name) is not None:
                    default_contaminants = ''
                    contaminants_file_name = ''
                    # read contaminants file
                    matchobj = re.match(
                        r'(.*) contaminants_file[\s\t]*([\S]*)[\s\t]*',
                        preprocess_cmd,
                        re.M | re.I)
                    if matchobj:
                        contaminants_file = matchobj.group(2)
                        contaminants_file_name = contaminants_file
                        contaminants_file_arr = read_file(contaminants_file)
                        if len(contaminants_file_arr) == 0:
                            status_check = "unhealthy"
                            parser_errorlist.append(
                                " ERROR  : No sequences in the contaminant file , specify 'contaminats_file' == "
                                "'no' or 'filename' in preprocess configuration file" +
                                contaminants_file)
                        preprocess_cmd = preprocess_cmd.replace(
                            'contaminants_file', '')
                        preprocess_cmd = preprocess_cmd.replace(
                            contaminants_file, '')
                    else:
                        default_contaminants = 'True'
                        contaminants_default_file = yap_src + \
                            "/" + "contaminants_default_file"
                        contaminants_file_name = contaminants_default_file
                        contaminants_default_file_arr = read_file(
                            contaminants_default_file)
                        contaminant_default_list = ''
                        contaminants_file_arr = []
                        for i in range(0, len(contaminants_default_file_arr)):
                            contaminant_default_list = contaminant_default_list + \
                                "," + \
                                contaminants_default_file_arr[i].strip("\n")
                        #contaminant_default_list = contaminant_default_list.lstrip(",") + "\n"
                        for i in range(0, len(inp_files_list)):
                            cont_file1 = inp_files_list[i][0]
                            contaminants_file_arr.append(
                                cont_file1 + "\t" + contaminant_default_list)
                            if workflow_config_dict["paired_end_data"] == "yes":
                                cont_file2 = inp_files_list[i][1]
                                contaminants_file_arr.append(
                                    cont_file2 +
                                    "\t" +
                                    contaminant_default_list)
                    cont_file_list = []
                    for i in range(0, len(contaminants_file_arr)):
                        matchobj = re.match(
                            r'\s*(\S*)(\s*).*',
                            contaminants_file_arr[i].strip("\n"),
                            re.M | re.I)
                        if matchobj:
                            split_by = matchobj.group(2)
                        cont_filename = contaminants_file_arr[
                            i].split(split_by)[0]
                        cont_file_list.append(cont_filename)
                    # print inp_files_list
                    # print inp_files
                    conflict_validator = yap_conflict_check(inp_files)
                    match_list, cont_not_exist_list, duplicate_dict = conflict_validator.validate_names_and_find_duplicates(
                        cont_file_list)
                    if cont_not_exist_list:
                        # all names which occur in contaminant file but not
                        # exist in input file list.
                        for ne in cont_not_exist_list:
                            error_str = "Error: In " + contaminants_file_name + ", file: '" + ne + \
                                "' is not one of the input files specified for this workflow. Please specify one of the input filenames for contaminance check."
                            cont_check_error_list.append(error_str)
                    if duplicate_dict:
                        # e.x. input file list contains /d1/f1, /d2/f1,
                        # then if f1 occurs in contaminants file,
                        # this duplicate happens
                        for d in duplicate_dict:
                            duplicate_str = "Error: In " + contaminants_file_name + ", name:'" + d + "' was found in " + \
                                conflict_validator.list_to_sentence(duplicate_dict[
                                                                    d]) + ". Please specify the " + "fullpath for this file to remove ambiguity."
                            cont_check_error_list.append(duplicate_str)
                    duplicate_name_dict = conflict_validator.find_duplicates_in_list(
                        match_list)
                    if duplicate_name_dict:
                        # in contaminants file, if one name occurs more than one time,
                        # this duplicate happens.
                        for dn in duplicate_name_dict:
                            duplicate_name_str = "Error: In " + contaminants_file_name + \
                                ", name:'" + dn + "' occurs several times. " + \
                                "Please give unique name."
                            cont_check_error_list.append(duplicate_name_str)

                    try:
                        for i in range(0, len(contaminants_file_arr)):
                            matchobj = re.match(
                                r'\s*(\S*)(\s*).*',
                                contaminants_file_arr[i].strip("\n"),
                                re.M | re.I)
                            if matchobj:
                                split_by = matchobj.group(2)
                            cont_filename = contaminants_file_arr[
                                i].split(split_by)[0]
                            contaminant_list = contaminants_file_arr[
                                i].split(split_by)[1]
                            # get input file names from contaminant file
                            path_name, file_name1 = os.path.split(
                                inp_files_list[0][0])
                            file_name1, extension = os.path.splitext(
                                file_name1)
                            file_ext = extension
                            while extension != '':
                                file_name1, extension = os.path.splitext(
                                    file_name1)
                                file_ext = file_ext + extension
                            tmp_path, cont_filename = os.path.split(
                                cont_filename)
                            cont_filename, tmp_ext = os.path.splitext(
                                cont_filename)
                            while tmp_ext != '':
                                cont_filename, tmp_ext = os.path.splitext(
                                    cont_filename)

                            cont_filename = path_name + \
                                "/" + cont_filename + file_ext
                            contaminant_list_arr = contaminant_list.split(",")
                            contaminant_dict[
                                cont_filename] = contaminant_list_arr
                    except Exception as e:
                        status_check = "unhealthy"
                        parser_errorlist.append(
                            "Format Error: in contaminants file")
                    #        parser_errorlist.append(e)
                if re.search('fastx_barcode_splitter', preprocess_cmd_name) is not None:
                    matchobj = re.match(
                        r'(.*) --bcfile[\s\t]*([\S]*)[\s\t]*',
                        preprocess_cmd,
                        re.M | re.I)
                    cmd_delete_index = jj
                    if matchobj:
                        barcode_file = matchobj.group(2)
                        barcode_dict = read_barcodes(barcode_file)
                        barcode_dict[
                            "no_barcode_specified"] = "no_barcode_specified"
            # order the preprocess commands: fastqc, fastq screen , barcode
            # splitter and then rest of the preprocess commands in cfg file
            # order
            if barcode_splitter_index != '':
                barcode_splitter_cmd = preprocess_cmd_arr.pop(
                    barcode_splitter_index)
                preprocess_cmd_arr.insert(0, barcode_splitter_cmd)
            if fastq_screen_index != '':
                fastq_screen_cmd = preprocess_cmd_arr.pop(fastq_screen_index)
                preprocess_cmd_arr.insert(0, fastq_screen_cmd)
            if fastqc_index != '':
                fastqc_cmd = preprocess_cmd_arr.pop(fastqc_index)
                preprocess_cmd_arr.insert(0, fastqc_cmd)
        # set paired end data for regroup and cuffdiff calculations
        paired_end_data = []
        if workflow_config_dict["paired_end_data"] == "yes":
            paired_end_data = inp_files_list
        # Re-group check
        # validate regroup file if regroup is specified in configuration
        regroup_file = ''  # regroup file path
        regroup_samples = {}  # regroup samples dictionary
        regroup_errors = {}  # regroup errors list
        regroup_done = False  # regroup execution status
        regroup_checker = None  # declare regroup checker object
        regroup_file, regroup_errors = \
            regroup_pre_validation_checks(workflow_config_dict)
        if len(regroup_errors) == 0 and len(regroup_file) > 0:
            # prepare regroup dictionary
            regroup_dict = get_regroup_dict(regroup_file,
                                            conflict_validator,
                                            inp_files, paired_end_data)
            # create regroup checker object
            regroup_checker = yap_regroup_check(regroup_dict)
            # validate regroup file using validator
            regroup_samples, regroup_errors  = \
                regroup_checker.validate_regroup_file()
            # set regroup execution flag
            regroup_done = True
        # print regroup errors if any
        if len(regroup_errors) > 0:
            # reset samples data if errors
            regroup_samples = {}
            # set status and extend error list
            status_check = "unhealthy"
            parser_errorlist.extend(regroup_errors)
        # add all types of samples comparisons
        add_sample_comparisons(workflow_config_dict, postprocess_cmd_arr)
        # cuffmerge and cuffcompare check
        compare_status = ''
        # cuffdiff check
        compare_status = execute_comparison_checks(
            workflow_config_dict,
            'list_of_samples_to_compare',
            yap_cuffdiff_check,
            regroup_samples,
            conflict_validator,
            inp_files,
            paired_end_data,
            parser_errorlist)
        if len(compare_status) > 0:
            status_check = "unhealthy"
        # End of cuffdiff processing loop
        compare_status = execute_comparison_checks(
            workflow_config_dict,
            'list_of_samples',
            yap_cuff_merge_compare_check,
            regroup_samples,
            conflict_validator,
            inp_files,
            paired_end_data,
            parser_errorlist)
        if len(compare_status) > 0:
            status_check = "unhealthy"
        # End of cuffmerge and cuffcompare processing loop
        # Print errors related to contaminants
        # extend contaminantls file errors
        if len(cont_check_error_list) > 0:
            # set status to unhealthy
            status_check = "unhealthy"
            parser_errorlist.extend(cont_check_error_list)
        # reassign the processed inp_files_list
        inp_files_list = n_inp_files_list
        # extend regroup samples
        if (regroup_checker is not None):
            regroup_checker.extend_regroup_samples(regroup_samples,
                                                   inp_files_list)
        # transform regroup samples dictionary for rest of yap
        regroup_samples_list = []
        for sample in regroup_samples:
            # get the tag name for each file in the sample
            sample_file_tags = []
            for file in regroup_samples[sample]:
                for inp_file in inp_files_list:
                    if file == inp_file[0]:
                        sample_file_tags.append(inp_file[2])
            regroup_samples_list.append([sample, sample_file_tags])
        # validate if output path exists
        if os.path.exists(workflow_config_dict["output_files_path"]):
            pass
        else:
            status_check = "unhealthy"
            parser_errorlist.append(
                "Error: The output files path does not exist")
            parser_errorlist.append(
                "Please make sure that user has the write access to output files directory ")
        if os.access(workflow_config_dict["output_files_path"], os.W_OK) == False:
            status_check = "unhealthy"
            parser_errorlist.append(
                "Error: User does not have write permission on output files path")
            parser_errorlist.append(
                "Please make sure that user has the write access to output files directory ")
        # print current workflow summary
        # add workflow information to the results
        results.append(
            "::new_section_title: YAP ANALYSIS SUMMARY FOR WORKFLOW = " +
            workflow_config_dict["workflow_name"] +
            " ")
        # append run mode specific data
        results.append("Operating System Information= " + ' '.join(os.uname()))
        results.append("USER= " + os.getenv("USER"))
        results.append("YAP SOURCE= " + os.getenv("YAP_HOME"))
        results.append("Python Source= " + (sys.version).replace("\n",' '))
        results.append(
            "Analysis Start Time For Workflow : " +
            workflow_config_dict["workflow_name"] +
            " " +
            time.strftime(
                "%Y/%m/%d %H:%M:%S",
                time.localtime()))
        results.append("YAP analysis general metadata: ")
        key_count = 0
        # print all keys in workflow struct and count of each key
        for key in workflow_struct[0].keys():
            key_count = key_count + 1
            results.append(
                str(key_count) + "." + key + ":" + workflow_struct[0][key])
        if run_mode == "--check":
            results.append(
                "Number of workflows provided= " + str(len(workflow_struct) - 1))
        # print workflow specific information
        results.append(
            "Instrument Type= " + workflow_config_dict["instrument_type"])
        results.append(
            "Specimen Information= " + workflow_config_dict["specimen_info"])
        results.append(
            "Workflow type= " + workflow_config_dict["workflow_type"])
        if workflow_config_dict["paired_end_data"] == "yes":
            results.append(
                "Number of input files= " + str(len(inp_files_list) * 2))
        else:
            results.append(
                "Number of input files= " + str(len(inp_files_list)))
        if "nprocs" in workflow_config_dict:
            results.append(
                "Number of processors= " + workflow_config_dict["nprocs"])
        results.append("Input files path for the workflow= " +
                       workflow_config_dict["input_files_path"])
        results.append("Input file provided: ")
        for i in range(0, len(inp_files_list)):
            if inp_files_list[i][2] == '' and inp_files_list[i][3] == '':
                if inp_files_list[i][0] != '':
                    results.append(str(i + 1) + "." + inp_files_list[i][0])
                if inp_files_list[i][1] != '':
                    results.append(" " + inp_files_list[i][1])
            else:
                results.append(str(i +
                                   1) +
                               "." +
                               inp_files_list[i][2] +
                               " => " +
                               inp_files_list[i][0])
                results.append("\t\t" + inp_files_list[i][1])
                results.append("\t")
        results.append("Output file path for the workflow= " +
                       workflow_config_dict["workflow_output_path"])
        if workflow_config_dict["paired_end_data"] == "yes":
            results.append("Sequence data type= paired end")
        else:
            results.append("Sequence data type= single end")
        results.append(
            "Input file format= " + workflow_config_dict["input_file_format"])
        results.append(
            "Maximum read length= " + workflow_config_dict["max_read_length"])
        results.append(
            "File chunk size (in megabytes)= " +
            workflow_config_dict["file_chunk_size"])
        results.append(
            "Data distribution method=" +
            workflow_config_dict["data_distribution_method"])
        results.append(
            "Output file path= " +
            workflow_config_dict["workflow_output_path"])
        results.append("::new_section")
        results.append("Analysis stages :")
        results.append(
            "Preprocess analysis= " +
            workflow_config_dict["run_preprocess_analysis"])
        results.append("Reference Sequence Alignment=" +
                       workflow_config_dict["run_reference_alignment"])
        results.append(
            "Postprocess Analysis= " +
            workflow_config_dict["run_postprocess_analysis"])
        results.append("::new_section")
        if workflow_config_dict["run_preprocess_analysis"] == "yes":
            if len(preprocess_cmd_arr) > 0:
                results.append("Preprocess Analysis commands:")
                results.append("Barcodes information: "),
                for i in barcode_dict:
                    results.append(i + " : " + barcode_dict[i])
                key_count = 0
                for i in range(0, len(preprocess_cmd_arr)):
                    cmd_type = preprocess_cmd_arr[i][0]
                    cmd_meta_data = preprocess_cmd_arr[i][1]
                    temp_arr = preprocess_cmd_arr[i][2]
                    preprocess_cmd_name = temp_arr[0][0]
                    preprocess_cmd = temp_arr[0][1]
                    results.append(str(i +
                                       1) +
                                   "." +
                                   " command name= " +
                                   preprocess_cmd_name +
                                   "," +
                                   "command line= " +
                                   preprocess_cmd)
                    if re.search('fastx_clipper', preprocess_cmd_name) is not None:
                        results.append("::new_section")
                        if len(contaminants_file) > 0:
                            results.append(
                                "contaminants_file= " + contaminants_file)
                        results.append(
                            "Contaminats used in fastx clipper analysis:")
                        contaminants_arr = []
                        for key in contaminant_dict.keys():
                            contaminants_arr = contaminant_dict[key]
                            if default_contaminants == 'True':
                                for ff in range(0, len(contaminants_arr)):
                                    # if contaminants_arr[ff] != '':
                                    results.append(contaminants_arr[ff])
                                break
                            else:
                                results.append(key)
                                for ff in range(0, len(contaminants_arr)):
                                    results.append(contaminants_arr[ff])
        results.append("::new_section")
        if workflow_config_dict["run_reference_alignment"] == "yes":
            results.append("Aligner commands:")
            for i in range(0, len(aligner_cmd_arr)):
                cmd_type = aligner_cmd_arr[i][0]
                cmd_meta_data = aligner_cmd_arr[i][1]
                temp_arr = aligner_cmd_arr[i][2]
                aligner_cmd_name = temp_arr[0][0]
                aligner_cmd = temp_arr[0][1]
                results.append(str(i +
                                   1) +
                               "." +
                               " command name= " +
                               aligner_cmd_name +
                               "," +
                               "command line= " +
                               aligner_cmd)
            results.append("Alignment output data sort order= " +
                           workflow_config_dict["alignment_sort_order"])
        results.append("::new_section")
        # print regroup samples information
        results.append("Samples re-grouped in this workflow:")
        if len(regroup_samples_list) == 0:
            results.append("None.")
        else:
            for item in regroup_samples_list:
                sample = item[0]
                results.append(
                    "Tags regrouped under name:'" + sample + "' are: ")
                for file in item[1]:
                    results.append("\t" + file)
        results.append("::new_section")
        if workflow_config_dict["run_postprocess_analysis"] == "yes":
            results.append("Potprocess analysis commands:")
            for i in range(0, len(postprocess_cmd_arr)):
                cmd_type = postprocess_cmd_arr[i][0]
                cmd_meta_data = postprocess_cmd_arr[i][1]
                cmd_arr = postprocess_cmd_arr[i][2]
                results.append(str(i + 1) + "." + " command type= " + cmd_type)
                results.append("\tcommand input : " + str(cmd_meta_data))
                for jj in range(0, len(cmd_arr)):
                    postprocess_cmd_name = cmd_arr[jj][0]
                    postprocess_cmd = cmd_arr[jj][1]
                    results.append("\t" +
                                   str(jj +
                                       1) +
                                   "." +
                                   " command name= " +
                                   postprocess_cmd_name +
                                   "," +
                                   "command line= " +
                                   postprocess_cmd)
            # print comparison information
            print_comparison_data(
                workflow_config_dict, 'list_of_samples_to_compare', results)
            print_comparison_data(
                workflow_config_dict, 'list_of_samples', results)
            results.append("::new_section")
        # append eqp warning to warning list (only meant for in-house use, neglect otherwise)
        if len(eqp_warning) > 0:
            warning_list.append(eqp_warning)
        # get status of all checks
        syn_stat = get_check_status(parser_errorlist, [])
        com_stat = get_check_status(incompatibility_error_list, warning_list)
        path_stat = get_check_status(missing_path_errors, basename_warnings)
        # get overall check status
        if "Failed" in [syn_stat, com_stat, path_stat]:
            stat = "Failed"
        elif "Passed With Warnings" in [syn_stat, com_stat, path_stat]:
            stat = "Passed With Warnings"
        else:
            stat = "Passed"
        # yap check summary
        results.append(
            "\n******************* YAP CHECK SUMMARY *******************")
        str1 = "* --Syntax check          : " + syn_stat
        str_format_in_box(str1, results, 57)
        str1 = "* --Compatibility check   : " + com_stat
        str_format_in_box(str1, results, 57)
        str1 = "* --File paths check      : " + path_stat
        str_format_in_box(str1, results, 57)
        results.append("*" * 57)
        str1 = "* YAP Configuration overall check status: " + stat
        results.append(str1)
        results.append("")
        results.append("::new_section_title: YAP Check Error/Warning Info ")
        # yap check detail info
        get_detail_error_info(parser_errorlist, [], "Syntax check", results)
        get_detail_error_info(
            incompatibility_error_list,
            warning_list,
            "Compatibility check",
            results)
        get_detail_error_info(
            missing_path_errors,
            basename_warnings,
            "File paths check",
            results)
        results.append("\n")
        results.append("::new_section")
        results.append(
            "::new_section_title: YAP configurations check end for Workflow = " +
            workflow_config_dict["workflow_name"] +
            " ")
        #----------------------------------------------------------------------
        # extend main error list with errors from all sections
        parser_errorlist.extend(incompatibility_error_list)
        parser_errorlist.extend(missing_path_errors)
        # update workflow config dict with required data
        workflow_config_dict["inp_files_list"] = inp_files_list
        workflow_config_dict["preprocess_cmd_arr"] = preprocess_cmd_arr
        workflow_config_dict["aligner_cmd_arr"] = aligner_cmd_arr
        workflow_config_dict["postprocess_cmd_arr"] = postprocess_cmd_arr
        workflow_config_dict["barcode_dict"] = barcode_dict
        workflow_config_dict["contaminant_dict"] = contaminant_dict
        workflow_config_dict["regroup_samples"] = regroup_samples_list
        #if len(regroup_samples_list) > 0:
        workflow_config_dict["regroup_output_path"] = workflow_config_dict["workflow_output_path"] + "/" + "regroup_output"
        # get split struct for paired end files
        workflow_config_dict["paired_files_split_arr"] = split_files_each_proc(inp_files_list, int(workflow_config_dict["nprocs"]))
        workflow_config_dict["check_results"] = results
        # creating aligner and postprocess command
        if workflow_config_dict["run_reference_alignment"] == "yes":
        # reading aligner configuration file into array
            aligner_output_key_arr = []
            for jj in range(0, len(aligner_cmd_arr)):
                cmd_type = aligner_cmd_arr[jj][0]
                cmd_meta_data = aligner_cmd_arr[jj][1]
                temp_arr = aligner_cmd_arr[jj][2]
                aligner_cmd_name = temp_arr[0][0]
                aligner_cmd = temp_arr[0][1]
                aligner_output_key = ''
                aligner_outfile_pos = aligner_cmd.rfind("output_file")
                for kk in range(aligner_outfile_pos, len(aligner_cmd)):
                    if aligner_cmd[kk] != ' ':
                        aligner_output_key += aligner_cmd[kk]
                    else:
                        break
                if aligner_output_key != '':
                    aligner_output_key_arr.append(aligner_output_key)
                if re.search('tophat', aligner_cmd_name) is not None:
                    aligner_output_key = 'accepted_hits.bam'
                    aligner_output_key_arr.append(aligner_output_key)
            aligner_output_key_arr.sort(key=len)
            aligner_output_key_arr = aligner_output_key_arr[::-1]
            workflow_config_dict["aligner_output_key_arr"] = aligner_output_key_arr
        # get file chunk information for each file in files list.  and create
        # output_directory structure for each file
        # converting megabytes into bytes. chunk_size is always in Megabytes(user
        # input)
        bytes_chunk_size = int(workflow_config_dict["file_chunk_size"]) * 1024 * 1024
        file_split_struct = []
        for i in range(0, len(inp_files_list)):
            file_split_info = []
            chunk_size = bytes_chunk_size
            each_file = inp_files_list[i][2]
            file_basecount_dict.setdefault(each_file, {})
    #update workflow config with filebase count dict
        workflow_config_dict["file_basecount_dict"] = file_basecount_dict
        workflow_config_dict.setdefault('list_of_samples', {})
        workflow_config_dict.setdefault('list_of_samples_to_compare', {})
        if run_mode == "--check" or len(parser_errorlist) > 0:
            print_list(results)
        return parser_errorlist


def get_check_status(error_list, warning_list):
    """
    Return check status.
            --"Passed":               no error, no warning
            --"Passed With Warnings": no error, but have warnings
            --"Failed":               have errors

    """
    if error_list:
        return "Failed"
    elif warning_list:
        return "Passed With Warnings"
    else:
        return "Passed"


def get_detail_error_info(error_list, warning_list, section, results):
    """
    Append detailed error info and warning info to results.
    """
    if error_list:
        results.append("::new_section")
        results.append("--YAP Configuration " + section + " status: Failed")
        for i in error_list:
            results.append(i)
        for i in warning_list:
            results.append(i)
    elif warning_list:
        results.append("::new_section")
        results.append("--YAP Configuration " + section +
                       " status: Passed With Warnings")
        for i in warning_list:
            results.append(i)


def str_format_in_box(str1, list, box_width):
    space_num = get_character_num(str1, box_width) - 1
    if space_num == 0:
        str1 += "*"
    else:
        str1 += ' ' * space_num
        str1 += "*"
    list.append(str1)


def add_sample_comparisons(workflow_config_dict, postprocess_cmd_arr):
    """
    Given workflow dict and postprocess command array,
    Updates dictionary with list of commands in which sample comparisons are
    required
    """
    file = ''  # output variable
    key = ''  # sample comparison key
    if workflow_config_dict["run_postprocess_analysis"] == "yes":
        for command in postprocess_cmd_arr:
            command_text = command[2]  # get command text
            for line in command_text:
                cmd_args = line[1]  # get line arguments
                matchobj = \
                    re.match(r'(.*) list_of_samples_to_compare' +
                             '[\s\t]*([\S]*)[\s\t]*', cmd_args, re.M | re.I)
                if matchobj is not None:
                    key = 'list_of_samples_to_compare'
                    file = matchobj.groups()[1]
                else:
                    matchobj = re.match(
                        r'(.*) list_of_samples' +
                        '[\s\t]*([\S]*)[\s\t]*',
                        cmd_args,
                        re.M | re.I)
                    if matchobj is not None:
                        key = 'list_of_samples'
                        file = matchobj.groups()[1]
                # add key to dict
                if len(key) > 0:
                    if key in workflow_config_dict:
                        # line 0 contains unique names
                        workflow_config_dict[key][line[0]] =\
                            [file, []]
                    else:
                        workflow_config_dict[key] =\
                            {line[0]: [file, []]}

def comparison_pre_validation_checks(
        workflow_config_dict,
        wf_key,
        command_name):
    """
    Checks if provided workflow is correct for
    commands which need to compare samples
    """

    file = ''
    file = workflow_config_dict[wf_key][command_name][0]
    if len(file) > 0:
        if file.strip() == "all":
            workflow_config_dict[wf_key][command_name] = ['all', []]
            return ''  # return empty if all
        else:
            return file
    else:
        workflow_config_dict[wf_key][command_name] = {}
        return ''  # return empty o/p as default

def regroup_pre_validation_checks(workflow_config_dict):
    """
    Checks if provided configuration is correct for regroup
    validation.
    Returns valid regroup filename if found, error list otherwise
    """
    if (workflow_config_dict["run_reference_alignment"] != "yes"):
        return '', []  # return default output
    if "regroup_output" in workflow_config_dict:
        if workflow_config_dict["regroup_output"] == "yes":
            regroup_file = workflow_config_dict["regroup_file"]
            if regroup_file != '':
                return regroup_file, []
            else:
                return '',\
                    ["Error: Invalid workflow configuration. " +
                     " regroup_output is set to yes, but no reroup " +
                     " file name is specified. Please specify the file" +
                     " name or set regroup_output to no."]
    return '', []  # return default output

def get_comparison_dict(input_file, regroup_samples, conflict_validator,
                        inp_files, paired_end_data):
    """
    Given all required arguments for validating comparison data
    Returns a comparison dictionary object with all key value pairs required by
    comparison checker object
    """
    compare_dict = {}  # output variable
    compare_dict["input_file"] = input_file
    compare_dict["samples"] = regroup_samples
    compare_dict["file_reader"] = yap_file_io.read_file
    compare_dict["paired_end_data"] = paired_end_data
    # set conflict validator
    if conflict_validator is None:
        conflict_validator = yap_conflict_check(inp_files)
    compare_dict["conflict_validator"] = conflict_validator
    return compare_dict  # return output

def get_regroup_dict(
        regroup_file,
        conflict_validator,
        inp_files,
        paired_end_data):
    """
    Given a regroup file, conflict validator and list of paired end and input
    files, Returns a dict object with all key value pairs required by regroup
    checker object
    """
    regroup_dict = {}  # output variable
    regroup_dict["regroup_file"] = regroup_file
    regroup_dict["file_reader"] = yap_file_io.read_file
    # set conflict validator
    if conflict_validator is None:
        conflict_validator = yap_conflict_check(inp_files)
    regroup_dict["conflict_validator"] = conflict_validator
    regroup_dict["paired_end_data"] = paired_end_data
    return regroup_dict  # return output

def print_comparison_data(workflow_config_dict, wf_key, results):
    if wf_key in workflow_config_dict:
        for key in workflow_config_dict[wf_key]:
            results.append("::new_section")
            results.append(key + " group(s) in this workflow: ")
            results.append("\t")
            file = workflow_config_dict[wf_key][key][0]
            if file == 'all':
                results.append("\tAll samples")
                continue
            group_counter = 0
            groups = workflow_config_dict[wf_key][key][1]
            for group in groups:
                group_counter += 1
                results.append("Below sets of files/samples " +
                               "in " + key + " group: " + str(group_counter))

                for i in range(0, len(group)):
                    set = group[i]
                    results.append("\tSet " + str(i + 1) + ": ")
                    for file in set:
                        results.append("\t\t" + file)
                results.append("\t")

def execute_comparison_checks(workflow_config_dict, wf_key,
                              compare_mod,
                              regroup_samples,
                              conflict_validator,
                              inp_files, paired_end_data,
                              parser_errorlist):
    """
    performs comparison checks and return status if unhealthy
    """
    compare_file = ''
    comparison_groups = {}  # groups from each compare file
    comparison_errors = {}
    status_check = ''
    if wf_key in workflow_config_dict:
        for key in workflow_config_dict[wf_key]:
            compare_file = \
                comparison_pre_validation_checks(workflow_config_dict,
                                                 wf_key, key)
            if len(compare_file) > 0:
                # prepare comparison dictionary
                compare_dict = get_comparison_dict(compare_file,
                                                   regroup_samples,
                                                   conflict_validator,
                                                   inp_files, paired_end_data)
                compare_dict["cmd_name"] = key
                # create comparison checker object
                compare_checker = compare_mod(compare_dict)
                # validate comparison file using validator
                comparison_groups, comparison_errors  = \
                    compare_checker.validate_compare_file()
                # extend main errorlist with  errors if any
                if len(comparison_errors) > 0:
                    # reset groups data if errors
                    comparison_groups = {}
                    # set status and extend error list
                    if len(status_check) == 0:
                        status_check = "unhealthy"
                    parser_errorlist.extend(comparison_errors)
                # transform for rest of yap
                groups_list = []
                for group in comparison_groups:
                    groups_list.append(comparison_groups[group])
                # append the data
                workflow_config_dict[wf_key][key][1] =\
                    groups_list
    return status_check
