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

import glob
import re
import sys
input_path = sys.argv[1]
output_file = "workflow_summary.txt"

def update_index(temp):
    count = 0
    for i in range(len(temp)):
        if temp[i] != '':
            temp_match = re.match(r'(\d)\..*', temp[i])
            if temp_match:
                count += 1
                temp[i] = re.sub(r'^\d(\.)', str(count) + r'\1', temp[i])
    return '\n'.join(temp)


def update_input_file(metadata):
    template = metadata[1]
    total_input_path_set = set()
    total_samples = []
    for i in metadata:
        path_match = re.search(
            r'Input files path for the workflow=\n?\s?(.*)\n', i)
        total_input_path_set.update(path_match.group(1).split(';'))
        sample_match = re.search(
            r'Input file provided:\n(.*)\nSequence data type=', i, re.DOTALL)
        total_samples.extend(sample_match.group(1).split('\n'))
    template = re.sub(r'(Input files path for the workflow=\n?\s?).*\n\t*\n?',
                      r'\1' + ';'.join(total_input_path_set) + r'\n', template)
    template = re.sub(
        r'(Input file provided:).*(Sequence data type=)',
        r'Input file provided: \n' +
        update_index(total_samples) +
        r'\nSequence data type=',
        template,
        0,
        re.DOTALL)
    return template


def multiple_split(in_string, pattern_list):
    out_list = []
    for i in pattern_list:
        temp_split = []
        temp_split = in_string.split(i)
        out_list.append(temp_split[0])
        if len(temp_split) > 1:
            in_string = temp_split[1]
        else:
            break
    return out_list


def combine_provenance(prov):
    temp_lis = ['PREPROCESS: \n', 'ALIGNMENT : \n',
                'MERGE ALIGNMENT OUTPUT : \n', 'POSTPROCESS : \n']
    pre_cmd = ''
    align_cmd = ''
    merge_cmd = ''
    post_cmd = ''
    head = ''
    for i in prov:
        split_list = []
        split_list = multiple_split(i, temp_lis)
        if head == '':
            head = split_list[0]
        try:
            pre_cmd += "PREPROCESS :\n" + split_list[1]
            align_cmd += "ALIGNMENT :\n" + split_list[2]
            merge_cmd += "MERGE ALIGNMENT OUTPUT :\n" + split_list[3]
            post_cmd += "POSTPROCESS :\n" + split_list[4]
        except Exception:
            pass
    return head + pre_cmd + align_cmd + merge_cmd + post_cmd


def combine_workflow_summary(input_path):
    file_dict = {
        'metadata': '',
        'wf_param': [],
        'stage': '',
        'preprocess': '',
        'barcode': '',
        'alignment': '',
        'regroup': '',
        'postprocess': '',
        'check_summary': '',
        'error_warning': '',
        'provenance': [],
        'end_time': ''}
    file_list = glob.glob(input_path + "/*/*workflow_summary.txt")
    for i in file_list:
        split_file = []
        split_file = open(i).read().split('-' * 109)
        if file_dict['metadata'] == '':
            file_dict['metadata'] = split_file[0]
        file_dict['wf_param'] += [split_file[1]]
        if file_dict['stage'] == '':
            file_dict['stage'] = split_file[2]
        if file_dict['preprocess'] == '':
            file_dict['preprocess'] = split_file[3]
        if file_dict['barcode'] == '':
            file_dict['barcode'] = split_file[4]
        if file_dict['alignment'] == '':
            file_dict['alignment'] = split_file[5]
        if file_dict['regroup'] == '':
            file_dict['regroup'] = split_file[6]
        if file_dict['postprocess'] == '':
            file_dict['postprocess'] = split_file[7]
        if file_dict['check_summary'] == '':
            file_dict['check_summary'] = split_file[8]
        if file_dict['error_warning'] == '':
            file_dict['error_warning'] = split_file[9]
        file_dict['provenance'] += [split_file[10]]
        if file_dict['end_time'] == '':
            file_dict['end_time'] = split_file[11]
    return file_dict['metadata'] + "-" * 109 + update_input_file(file_dict['wf_param']) + "-" * 109 + file_dict['preprocess'] + "-" * 109 + file_dict['barcode'] + "-" * 109 + file_dict['alignment'] + "-" * 109 + file_dict[
        'regroup'] + "-" * 109 + file_dict['postprocess'] + "-" * 109 + file_dict['check_summary'] + "-" * 109 + file_dict['error_warning'] + "-" * 109 + combine_provenance(file_dict['provenance']) + "-" * 109 + file_dict['end_time']
output_path = glob.glob(input_path)[0]
output_path = re.sub(r'_batch_[\d]$', '', output_path)
output_file = output_path + output_file
print output_file
fh = open(output_file, 'w')
fh.write(combine_workflow_summary(input_path))
combine_workflow_summary(input_path)
fh.close()
