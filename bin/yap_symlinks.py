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
import glob
from subprocess import *

""" 
Given a yap output path, creates symlinks of the aligner_output across all samples at a given output directory 
"""

if len(sys.argv) != 5:

    print 'Please enter the respective input and output path directories in the format: \n\t\t yap_symlinks.py -i <input_path> -o <output_path>'

else:
    yap_symlinks_input = ''
    yap_symlinks_output = ''
    yap_symlinks_input = sys.argv[2]
    yap_symlinks_output = sys.argv[4]
    yap_symlinks_output = yap_symlinks_output.rstrip('\n')

if os.path.exists(yap_symlinks_input) and os.path.exists(yap_symlinks_output) == 'False':
    print 'Please see to it that the specified paths exist.'
else:
    yap_symlinks_aligner_output = glob.glob(
        yap_symlinks_input + '/*/*/aligner_output/*')

    def unique_file_name(file_path_list):
        unq = set()
        for i in file_path_list:
            if i not in unq:
                unq.add(i)
                yield i

    for i in yap_symlinks_aligner_output:
        temp_arr = []
        alignment_filename = ''
        temp_arr = i.split('/')
        if temp_arr[-1].find(temp_arr[-4]) != -1:
            alignment_filename = temp_arr[-1]
        else:
            alignment_filename = temp_arr[-4] + \
                "_" + temp_arr[-3] + "_" + temp_arr[-1]
        alignment_filename = alignment_filename.replace(
            '_no_barcode_specified', '')
        alignment_filename = alignment_filename.replace('__', '_')
        alignment_filename = alignment_filename.strip('_')
        psymlink = Popen('ln -s ' + i + ' ' + yap_symlinks_output +
                         '/' + alignment_filename, shell='True').wait()
