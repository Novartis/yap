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
import glob
import re
import yap_file_io
from subprocess import Popen, PIPE
import yap_workflow_dict as wd


def merge_multiproc_files(command, filename, barcode, err_log, stat_log):
    """ Merges the temporary multiproc error files created. """
    exit_str = Popen("cat " + err_log + "_multiproc_* | grep EXIT_CODE",
                     stdout=PIPE, shell=True).communicate()[0]
    exit_code = 0
    for i in exit_str.split('\n'):
        m = re.match("EXIT_CODE: (.*)\n", i)
        if m:
            exit_code = exit_code + int(m.group(1))
    if exit_code == 0:
        yap_file_io.write_data(
            "YAP_COMMAND: " +
            command +
            "\nINPUT FILE: " +
            filename +
            "\n",
            stat_log)
        os.system("cat " + err_log + "_multiproc_* >>" + stat_log)
        yap_file_io.write_data("_" * 30 + "\n", stat_log)
        os.system("rm " + err_log + "_multiproc_*")
    else:
        yap_file_io.write_data(
            "YAP_COMMAND: " +
            command +
            "\nINPUT FILE: " +
            filename +
            "\n",
            err_log)
        os.system("cat " + err_log + "_multiproc_* >>" + err_log)
        yap_file_io.write_data("_" * 30 + "\n", err_log)
        os.system("rm " + err_log + "_multiproc_*")

def combine_temp_log(err_path, stat_path):
    """ Combines the temporary log file directory and 
        merges them with the sample error and status logs. """
    err_temp_list = glob.glob(err_path + "/*log_temp")
    stat_temp_list = glob.glob(stat_path + "/*log_temp")
    if len(err_temp_list) > 0:
        for i in err_temp_list:
            if os.path.isdir(i):
                file_name = i.replace("log_temp", "err.log")
                if len(glob.glob(i + "/*")):
                    os.system('cat ' + i + '/* >> ' + file_name)
                    os.system('rm ' + i + "/*")
    if len(stat_temp_list) > 0:
        for i in stat_temp_list:
            if os.path.isdir(i):
                file_name = i.replace("log_temp", "stat.log")
                if len(glob.glob(i + "/*")):
                    os.system('cat ' + i + '/* >> ' + file_name)
                    os.system('rm ' + i + "/*")

def merge_tee_files(command, filename, err_log, stat_log):
    """ Merges the temporary log files produced as a result of 
        the multiprocessing module. """
    exit_str = Popen("cat " + err_log + "_yap_tee_* | grep EXIT_CODE",
                     stdout=PIPE, shell=True).communicate()[0]
    exit_code = 0
    for i in exit_str.split('\n'):
        m = re.match("EXIT_CODE: (.*)\n", i)
        if m:
            exit_code = exit_code + int(m.group(1))
    if exit_code == 0:
        yap_file_io.write_data("YAP_COMMAND: " + command + "\nINPUT_FILES: " +
                   filename + "\nYAP_STATUS_MSG: ", stat_log)
        os.system("cat " + err_log + "_yap_tee_* >>" + stat_log)
        yap_file_io.write_data("_" * 30 + "\n", stat_log)
        os.system("rm " + err_log + "_yap_tee_*")
    else:
        yap_file_io.write_data("YAP_COMMAND: " + command + "\nINPUT_FILES: " +
                   filename + "\nYAP_ERROR_MSG:", err_log)
        os.system("cat " + err_log + "_yap_tee_* >>" + err_log)
        yap_file_io.write_data("_" * 30 + "\n", err_log)
        os.system("rm " + err_log + "_yap_tee_*")

def write_log(command, file, exit_code, std_err, err_log, stat_log):
    """ Checks if the command has succeeded/failed and
        logs it to the status/error log respectively. """
    cmd_sep = "_" * 30 + "\n"
    if str(exit_code) == '0':
        err_str = "YAP_COMMAND: %s\nINPUT_FILES: %s\nEXIT_CODE: %s\nYAP_STATUS_MSG: %s\n" % (
            command, file, exit_code, std_err)
        yap_file_io.write_data(err_str + cmd_sep, stat_log)
    else:
        err_str = "YAP_COMMAND: %s\nINPUT_FILES: %s\nEXIT_CODE: %s\nYAP_ERROR_MSG: %s\n" % (
            command, file, exit_code, std_err)
        yap_file_io.write_data(err_str + cmd_sep, err_log)

def pass_fail_dict(input_arr, dict):
    """ Called within the pass_fail_matrix() function. Accepts 
        the contents of a logfile as a list as input along with the sample_dict
        Gets the EXIT_CODE, err_log doesn't contain exit codes for a
        successful YAP run. """
    key = ""
    stage = "False"
    stage_dict = {'PREPROCESS': 'PREPROCESS',
                  'FASTQC': 'PREPROCESS',
                  'FASTQSCREEN': 'PREPROCESS',
                  'ALIGNMENT': 'ALIGNMENT',
                  'MERGE ALIGNMENT': 'ALIGNMENT',
                  'REGROUP': 'ALIGNMENT',
                  'POSTPROCESS': 'POSTPROCESS',
                  'CUFFDIFF': 'CUFFDIFF',
                  'CUFFMERGE': 'CUFFMERGE',
                  'CUFFCOMPARE': 'CUFFCOMPARE',
                  'MACS2': 'MACS2'}
    for i in range(len(input_arr)):
        # Pattern to match the start of a command
        start = re.match('\*{50}(.*) STARTED(.*)', input_arr[i])
        # Pattern to match the end of a command              
        finish = re.match('\*{50}(.*) FINISHED(.*)', input_arr[i])               
        exit_code = re.match("EXIT_CODE: (.*)", input_arr[i])
        if start and stage == "False":
            key = stage_dict[start.group(1)]
            stage = "True"
        elif exit_code and stage == "True":
            dict[key] = 'FAIL'
        elif finish:
            if key in dict.keys():
                pass
            else:
                dict[key] = 'PASS'
            stage = "False"
        elif stage == "False":
            pass

def print_matrix(dict, file):
    """ Called within the pass_fail_matrix() function.
        Input dict with sample name/command name if applicable with stagewise
        fail/pass dict and write path as input.
        prints the stagewise matrix to file. """
    normal_stage_arr = ['PREPROCESS', 'ALIGNMENT', 'POSTPROCESS']               
    # Stages of the matrix/sample
    compare_stage_arr = ['CUFFDIFF', 'CUFFCOMPARE', 'CUFFMERGE', 'MACS2']                
    # Multi-sample hence printed after.
    head_str = '\t' + '\t'.join(normal_stage_arr) + '\n'
    normal_str = ""
    compare_str = ""
    for i in sorted(dict.keys()):
        if sorted(dict[i].keys()) == sorted(normal_stage_arr):
            normal_str += i + "\t"
            for j in range(len(normal_stage_arr)):
                if j != len(normal_stage_arr) - 1:
                    normal_str += dict[i][normal_stage_arr[j]] + "\t"
                else:
                    normal_str += dict[i][normal_stage_arr[j]] + "\n"
        elif set(dict[i].keys()).issubset(set(compare_stage_arr)):
            for j in dict[i].keys():
                compare_str += i + ": " + dict[i][j] + "\n"
    yap_file_io.write_data(head_str + normal_str, file)
    if compare_str != '':
        yap_file_io.write_data("\n\n" + compare_str, file)

def pass_fail_matrix():
    """ Constructs the stagewise pass/fail matrix from the error_log. """
    err_log = wd.err_log_path
    matrix_path, junk = os.path.split(err_log)
    matrix_path += "/yap_pass_fail_matrix.log"
    sample_log_dict = {}
    pass_fail = {}
    for i in glob.glob(err_log + "/*"):
        path, file = os.path.split(i)
        file, ext = os.path.splitext(file)
        sample_log_dict[file.rstrip("_err")] = yap_file_io.read_file(i)  # len-1 chunks
    for i in sample_log_dict.keys():
        pass_fail_dict(sample_log_dict[i], pass_fail)
        if not set(pass_fail.keys()).issubset(set(['CUFFDIFF', 'CUFFCOMPARE', 'CUFFMERGE', 'MACS2'])):
            if wd.run_preprocess_analysis == "yes" and pass_fail.get("PREPROCESS") is None:
                pass_fail["PREPROCESS"] = "FAIL"
            elif wd.run_preprocess_analysis == "no":
                pass_fail["PREPROCESS"] = "N/A"
            if wd.run_reference_alignment == "yes" and pass_fail.get("ALIGNMENT") is None:
                pass_fail["ALIGNMENT"] = "FAIL"
            elif wd.run_reference_alignment == "no":
                pass_fail["ALIGNMENT"] = "N/A"
            if wd.run_postprocess_analysis == "yes" and pass_fail.get("POSTPROCESS") is None:
                pass_fail["POSTPROCESS"] = "FAIL"
            elif wd.run_postprocess_analysis == "no":
                pass_fail["POSTPROCESS"] = "N/A"
            if pass_fail["PREPROCESS"] == "FAIL":
                pass_fail["ALIGNMENT"] = "FAIL"
                pass_fail["POSTPROCESS"] = "FAIL"
        sample_log_dict[i] = pass_fail
        pass_fail = {}
    print_matrix(sample_log_dict, matrix_path)

def sample_filter():
    """ Filters failed samples after the alignment step. """
    ignore_list = []
    inp_files_list=[]
    list_of_samples=[]
    list_of_samples_to_compare=[]
    for i in wd.inp_files_list:
        pass_fail = {}
        err_log = wd.err_log_path + "/" + i[2] + "_err.log"
        if os.path.exists(err_log):
            pass_fail_dict(yap_file_io.read_file(err_log), pass_fail)
        if 'FAIL' in pass_fail.itervalues():
            ignore_list.append(i)
    if len(ignore_list) != 0:
        list_of_samples_to_compare = remove_corrupted_samples(wd.list_of_samples_to_compare,ignore_list)
        list_of_samples = remove_corrupted_samples(wd.list_of_samples, ignore_list)
        inp_files_list = [i for i in wd.inp_files_list if i not in ignore_list]
    	return inp_files_list,list_of_samples,list_of_samples_to_compare,ignore_list	
    else:
	return wd.inp_files_list, wd.list_of_samples, wd.list_of_samples_to_compare,ignore_list	

def remove_corrupted_samples(command_dict, ignore_list):
    """ Called withing the sample_filter() function
        Removes samples containing errors. Takes the 
        command_dict and the list of samples to be ignored as input. """
    if len(ignore_list) == 0:
        return command_dict
    new_dict = command_dict
    ignore_file_list = map((lambda x: x[0]), ignore_list)
    for key in command_dict:
        current_groups = command_dict[key][1]  # list of lines in file
        new_groups = []
        for group in current_groups:
            merged_sets = merge_group_sets(group)
            corrupt_samples = filter(
                (lambda x: x in ignore_file_list), merged_sets)
            if len(corrupt_samples) > 0:
                continue
            else:
                new_groups.append(group)
        new_dict[key][1] = new_groups
    return new_dict

def merge_group_sets(group):
    """ Called in the remove_corrupted_samples() function.
        'extends' group to merged_sets. """
    merged_set = []
    for i in range(0, len(group)):
        merged_set.extend(group[i])
    return merged_set
