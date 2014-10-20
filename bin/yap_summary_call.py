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

import sys
import os
import glob
from subprocess import PIPE, Popen

"""
-> Generates summary file which is a cumulation of the log file and all reports and graphs generated in the postprocess.
-> Generates special summary tables for htseq, fastqc, fastqscreen
"""
sample_dict = {}

if len(sys.argv) < 2:
    print " command line : ", "yap_summary [file with yap_output directory path]"
    exit()

summary_input_file = sys.argv[1]
if os.path.exists(summary_input_file):
    directory_info = open(summary_input_file, 'rb').readlines()
else:
    print "Yap summary input file does not exists"
    exit()


for input_string in directory_info:
    summary_cmd = 'yap_summary' + ' ' + input_string
    psum = Popen(summary_cmd, stdout=PIPE, shell='False', close_fds='True')
    psum.communicate()
