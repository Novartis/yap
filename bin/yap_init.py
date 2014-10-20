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
import yap_workflow_dict as wd
import yap_tools  
def initialize_dir_struct():
    print "Initializing YAP..."
    print "Creating output directory structure"
    #create temporary workspace directory
    yap_tools.create_dir(wd.yap_temp_user_dir)
    os.system("mkdir" + " " + wd.workflow_output_path)
    os.system("mkdir" + " " + wd.temp_dir_path)
    os.system("mkdir " + wd.consolidated_output)
    # create output directory for regrouped samples
    if wd.regroup_output =='yes':
        os.system("mkdir" + " " + wd.regroup_output_path)
    for i in wd.regroup_arr:  # iterate over regroup arr
        regroup_title = i[0]
        file_dir_path = wd.regroup_output_path + "/" + regroup_title
        os.system("mkdir" + " " + file_dir_path)
        for barcode in wd.barcode_dict.iterkeys():
            barcode_dir_path = file_dir_path + "/" + barcode
            os.system("mkdir" + " " + barcode_dir_path)
            preprocess_dir_path = barcode_dir_path + \
                "/" + "preprocess_output"
            aligner_output_dir_path = barcode_dir_path + \
                "/" + "aligner_output"
            barcode_basecount_dir = preprocess_dir_path + \
                "/" + "basecount_analysis"
            postprocess_dir_path = barcode_dir_path + \
                "/" + "postprocess_output"
            os.system("mkdir" + " " + postprocess_dir_path)
            os.system("mkdir" + " " + preprocess_dir_path)
            os.system("mkdir" + " " + aligner_output_dir_path)
    # create log file to track the regroup samples
    fw = open(wd.sample_track_log, 'wb')
    for group in wd.inp_files_list:
        each_file = group[2]
        file_dir_path = wd.workflow_output_path + "/" + each_file
        fw.write(each_file + "\t=>\t" + group[0] + "," + group[1] + "\n")
        os.system("mkdir" + " " + file_dir_path)
        for barcode in wd.barcode_dict.iterkeys():
            barcode_dir_path = file_dir_path + "/" + barcode
            os.system("mkdir" + " " + barcode_dir_path)
            preprocess_dir_path = barcode_dir_path + "/" + "preprocess_output"
            aligner_output_dir_path = barcode_dir_path + "/" + "aligner_output"
            barcode_basecount_dir = preprocess_dir_path + \
                "/" + "basecount_analysis"
            postprocess_dir_path = barcode_dir_path + \
                "/" + "postprocess_output"
            os.system("mkdir" + " " + postprocess_dir_path)
            os.system("mkdir" + " " + preprocess_dir_path)
            os.system("mkdir" + " " + aligner_output_dir_path)
    fw.close()
    os.system("mkdir " + wd.log_path)
    os.system("mkdir " + wd.err_log_path)
    os.system("mkdir " + wd.stat_log_path)
    for group in wd.inp_files_list:
        each_file = group[2]
        temp_err_path = ""
        temp_stat_path = ""
        temp_err_path = wd.err_log_path + "/" + each_file + "_log_temp"
        temp_stat_path = wd.stat_log_path + "/" + each_file + "_log_temp"
        os.system("mkdir" + " " + temp_err_path)
        os.system("mkdir" + " " + temp_stat_path)
