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
import sys

""" Splits a single yap configuration into multiple yap jobs of a specified number samples each. """

def split_inp_file_path(inp_list, new_inp_list):
    """ Splits the input files list into a sets 
        containing only a specified number of samples. """

    inp_path_list = []
    inp_file_list = []
    for i in inp_list:
        for file in i:
            input_path, input_file = file.rsplit('/', 1)
            inp_path_list.append(input_path)
            inp_file_list.append(input_file)
    new_inp_list.append([';'.join(inp_path_list), ';'.join(inp_file_list)])


def find_pair_end_files(file_list):
    """ Finds the corresponding pairs of paired end files 
        given a file list. """
    pairs = sorted(file_list)
    found_pairs = []
    found_set = []

    for i in range(len(pairs)):
        if pairs[i] in found_set:
            continue
        else:
            for j in range(len(pairs)):
                if i == j or pairs[j] in found_set:
                    continue
                else:
                    if pairs[i] not in found_set:
                        if len(pairs[i]) == len(pairs[j]):
                            diff = sum(
                                ch1 != ch2 for ch1,
                                ch2 in zip(
                                    pairs[i],
                                    pairs[j]))
                            if diff == 1:
                                for ch1, ch2 in zip(pairs[i], pairs[j]):
                                    if ch1 != ch2:
                                        if ch1 == '1' and ch2 == '2' or ch1 == '2' and ch1 == '1':
                                            found = []
                                            found.append(pairs[i])
                                            found.append(pairs[j])
                                            found_pairs.append(found)
                                            found_set.append(pairs[i])
                                            found_set.append(pairs[j])
                                        else:
                                            pass
                            else:
                                if j == len(pairs) - 1:
                                    pass
                                else:
                                    continue
                        else:
                            continue
                    else:
                        break
    not_found_set = list(set(pairs) - set(found_set))
    return found_pairs, not_found_set
if len(sys.argv) < 7:
    print "command line : ", "yap_batch_submit -i full_path_to_config_file  -o output_path[where the different batch " \
                             "directories will be created] -n number_of_samples_per_batch"
    exit()
else:
    workflow_config_file = os.path.abspath(sys.argv[2])
    output_path = os.path.abspath(sys.argv[4])
    sample_number_per_batch = int(sys.argv[6])

workflow_file_data = open(workflow_config_file).readlines(1)
input_path, workflow_file_name = os.path.split(workflow_config_file)
input_files_path = []
input_files = []
split_output = []
metadata_index = []
config = []
workflow_name = ''
for i in range(len(workflow_file_data)):
    workflow_file_data[i] = workflow_file_data[i].strip('\n')
    key_val = workflow_file_data[i].split(":=", 1)
    if len(key_val) == 2 and not workflow_file_data[i].startswith("#"):
        key, value = key_val
        key = key.strip().strip('"')
        value = value.split("#")[0].strip().strip('"')
        config.append([key, value])
        if key.strip() == "input_files_path":
            input_files_path = value.split(';')
        elif key == "input_files":
            input_files = value.split(';')
        elif key == "workflow_name":
            workflow_name = value
    else:
        config.append([key_val[0], ''])
total_inp_files = []
for i in range(len(input_files_path)):
    if input_files_path[i] != '':
        input_files_path[i] = input_files_path[i].strip()
        if input_files_path[i] == "./" or input_files_path[i] == ".":
            input_files_path[i] = input_path
        else:
            input_files_path[i] = os.path.abspath(input_files_path[i])

        input_files_arr = input_files[i].split(',')
        for j in input_files_arr:
            total_inp_files += glob.glob(input_files_path[i] + "/" + j)
paired_files = []
unpaired_files = []
new_input_list = []
nfiles = len(total_inp_files)
if nfiles % 2 == 0:
    paired_files, unpaired_files = find_pair_end_files(total_inp_files)
    split_inp_file_list = [
        paired_files[
            i:i +
            sample_number_per_batch] for i in range(
            0,
            len(paired_files),
            sample_number_per_batch)]
    for inp in split_inp_file_list:
        split_inp_file_path(inp, new_input_list)
for i in range(len(new_input_list)):
    base_dir = os.path.dirname(workflow_config_file)
    path, base_dir = os.path.split(base_dir)
    batch_path = output_path + "/" + base_dir + "_batch_" + str(i)
    sge = open(output_path + "/" + base_dir + "_batch_submit.sh", "a")
    if os.path.exists(batch_path):
        print "Directory already exists: ", batch_path
        print "Exiting..."
        exit()
    os.system("mkdir " + batch_path)
    os.system("cp " + input_path + "/* " + batch_path)
    out = open(batch_path + "/" + workflow_file_name, 'w+')
    for key in config:
        if key[0] == "input_files_path":
            key[1] = new_input_list[i][0]
        if key[0] == "input_files":
            key[1] = new_input_list[i][1]
        if key[0] == "workflow_name":
            out.write(
                "\"" +
                key[0] +
                "\"" +
                " := " +
                "\"" +
                key[1] +
                "_batch_" +
                str(i) +
                "\"\n")
        elif key[1] == '':
            out.write(key[0] + "\n")
        else:
            out.write("\"" + str(key[0]) + "\" := \"" + key[1] + "\"\n")
    out.close()
    sge.write("cd " + batch_path + "\n")
    sge.write("qsub yap_sge.sh" + "\n")
    sge.close()
