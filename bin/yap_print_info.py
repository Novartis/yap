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

import time
import os
from yap_printer import *
from yap_check import *
import yap_workflow_dict as wd 
def print_info():
    file_split_struct = []
    # Set default values for workflow
    f_summary_file =  wd.workflow_output_path + "/" + wd.workflow_name + "_workflow_summary.txt"
    f_summary = open(f_summary_file, 'wb')
    number_of_samples = len(wd.inp_files_list)
    if os.path.exists(wd.usage_log):
        f_usage = open(wd.usage_log, 'a')
        if wd.job_id != '':
            f_usage.write(
                wd.yap_version +
                "\t" +
                wd.username +
                "\t" +
                time.strftime(
                    "%Y/%m/%d %H:%M:%S",
                    time.localtime()) +
                "\t" +
                str(number_of_samples) +
                "\t" +
                str(wd.nprocs) +
                "\t" +
                str(wd.job_id) +
                "\t" +
                "START" +
                "\n")
        else:
            f_usage.write(
                wd.yap_version +
                "\t" +
                wd.username +
                "\t" +
                time.strftime(
                    "%Y/%m/%d %H:%M:%S",
                    time.localtime()) +
                "\t" +
                str(number_of_samples) +
                "\t" +
                str(wd.nprocs) +
                "\t" +
                str(wd.random_id) +
                "\t" +
                "START" +
                "\n")

        f_usage.close()
    else:
        "Warning: YAP usage log not found"
    # print information on screen and file
    print_list(wd.check_results, f_summary)
    # close the summary file
    f_summary.close()
#-------------------------------------------------------------------------


def print_usage():
    print "Options : "
    print " To print help : ", " yap --help or yap -h"
    print " To do configuration file format check : ", " yap --check [workflow_configuration_filename], ", "\"eg: yap --check workflow_configuration.cfg \""
    print " To run YAP: ", " yap -n [number of processors] [workflow_configuration_filename], ", "\"eg: yap -n 2 workflow_configuration.cfg \"", "(yap run with 2 processors)"
    print "exiting the program"
