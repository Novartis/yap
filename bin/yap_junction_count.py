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
from subprocess import PIPE, Popen

""" Counts the number of junctions.
    Usage: -exon_coordinates_file [ucsc coordinate bed file] 
           -exon_CoordToNumber_file [filename] 
           -i [inputfilename or for stdin '-'] 
           -o [output filename]   """
def print_usage():
    print "DESCRIPTION: Counts number of junctions .\n"
    print "USAGE:"
    print sys.argv[0] + "-exon_coordinates_file [ucsc coordinate bed file] -exon_CoordToNumber_file [filename] -i [inputfilename or for stdin '-'] -o [output filename]\n"
    exit()
if len(sys.argv) != 9:
    print_usage()
for i in range(0, len(sys.argv)):
    if sys.argv[i] == "-input_type":
        file_type = sys.argv[i + 1]
    if sys.argv[i] == "-exon_coordinates_file":
        exon_coordinates_file = sys.argv[i + 1]
    if sys.argv[i] == "-exon_CoordToNumber_file":
        exon_coordinate_number_file = sys.argv[i + 1]
    if sys.argv[i] == '-i':
        input_file = sys.argv[i + 1]
    if sys.argv[i] == '-o':
        output_file = sys.argv[i + 1]
final_output_file = output_file + "_junctioncount.txt"
unknown_junction_file = output_file + "_unknown_junction.txt"
temp_sorted_junctions = output_file + "_temp_sorted_junction.bed"
temp_sorted_junctionscount = output_file + "_temp_sorted_junctionscount.bed"
if input_file == "-":
    try:
        junction_cmd = " intersectBed -a - -b " + exon_coordinates_file + \
            " -wa -wb | cut -f1-6,16 > " + final_output_file
        prun = Popen(junction_cmd, stdin=PIPE, stdout=PIPE, shell='False')
        input_str = sys.stdin.read()
        out, err = prun.communicate(input_str)

    except:
        print "Error: Junction counts failed", "\ncommand=", junction_cmd
else:
    try:
        junction_cmd = " intersectBed -a " + input_file + " -b " + \
            exon_coordinates_file + \
            " -wa -wb | cut -f1-6,16 > " + final_output_file
        prun = Popen(junction_cmd, stdin=PIPE, stdout=PIPE, shell='False')
        out, err = prun.communicate()

    except:
        print "Error: Junction counts failed", "\ncommand=", junction_cmd
if input_file == "-":
    try:
        sort_junction_cmd = " cut -f4 " + " | " + \
            "sort -u " + "> " + temp_sorted_junctions
        prun = Popen(sort_junction_cmd, stdin=PIPE, stdout=PIPE, shell='False')
        out, err = prun.communicate(input_str)
    except:
        print "Error: Junction file sort  failed", "\ncommand=", sort_junction_cmd
else:
    try:
        sort_junction_cmd = " cut -f4 " + input_file + \
            " | " + "sort -u " + "> " + temp_sorted_junctions
        prun = Popen(sort_junction_cmd, stdin=PIPE, stdout=PIPE, shell='False')
    except:
        print "Error: Junction file sort  failed", "\ncommand=", sort_junction_cmd
try:
    sort_junctioncount_cmd = " cut -f4  " + final_output_file + \
        " | " + "sort -u " + "> " + temp_sorted_junctionscount
    prun = Popen(
        sort_junctioncount_cmd, stdin=PIPE, stdout=PIPE, shell='False')
    out, err = prun.communicate()
    unknown_junction_cmd = "diff " + temp_sorted_junctionscount + \
        " " + temp_sorted_junctions + " | " + " grep '^>' "
    prun = Popen(unknown_junction_cmd, stdin=PIPE, stdout=PIPE, shell='False')
    out, err = prun.communicate()
    unknwn_junc = []
    input_list = []
    unknwn_junc = out.replace('> ', '').split("\n")
    if input_file == "-":
        input_list = input_str.splitlines()
    else:
        with open(input_file, 'r') as fi:
            input_list = fi.readlines()
    uk = open(unknown_junction_file, 'w')
    for i in range(len(input_list)):
        temp_arr = []
        temp_arr = input_list[i].strip('\n').split('\t')
        if len(temp_arr) == 12:
            if str(temp_arr[3]) in unknwn_junc:
                output_str = temp_arr[
                    0] + '_' + temp_arr[1] + '_' + temp_arr[2] + '\t' + temp_arr[4] + '\n'
                uk.write(output_str)
    uk.close()
    os.remove(temp_sorted_junctionscount)
    os.remove(temp_sorted_junctions)
except:
    print "Error: In creating unknown junction file "
try:
    fh = open(exon_coordinate_number_file, 'rb')
except IOError as e:
    print "Error while opening the ", exon_coordinate_number_file
    print e
exon_number_arr = fh.readlines()
fh.close()
line = ''
exon_coord_num_dict = {}
for line in exon_number_arr:
    temp_arr = []
    temp_arr = line.strip("\n").split("\t")
    if len(temp_arr) < 5:
        print "Not enough fields in exon number file"
        print "please refer to the documentation"
        exit()
    else:
        exon_coordinate = temp_arr[3]
        exon_number = temp_arr[4]
        exon_coord_num_dict[exon_coordinate] = exon_number
try:
    fh = open(output_file + '_junctioncount.txt', 'rb')
except IOError as e:
    print "Error while opening the ", output_file + '_junctioncount.txt'
    print e
    exit()
junction_dict = {}
junction_range_dict = {}
for line in fh:
    temp_arr = []
    temp_arr = line.strip("\n").split("\t")
    key = ''
    if len(temp_arr) != 7:
        print "Not enough fields in junction counts file"
        exit()
    gene_name = (temp_arr[6].split("|"))[1]
    key = temp_arr[3] + "\t" + temp_arr[4] + \
        "\t" + temp_arr[5] + "\t" + gene_name
    begin_range = temp_arr[1]
    end_range = temp_arr[2]
    if temp_arr[6] in exon_coord_num_dict:
        exon_info = exon_coord_num_dict[temp_arr[6]]
        if key in junction_dict:
            junction_dict[key].append(exon_info)
        else:
            junction_dict[key] = [exon_info]
            junction_range_dict[key] = [begin_range, end_range, temp_arr[6]]
fh.close()
fh_out = open(output_file + "_junctioncount_summary.txt", 'wb')
for key in junction_dict.iterkeys():
    val = []
    exon_info = ''
    val = junction_dict[key]
    range_info = junction_range_dict[key]
    junction_info = key.split("\t")
    junction_name = junction_info[0]
    junction_length = junction_info[1]
    junction_direction = junction_info[2]
    gene_name = junction_info[3]
    old_exon_number = 0
    genename_direction_list = range_info[2].split("|")
    junction_direction = genename_direction_list[2]
    ref_range = genename_direction_list[0].split(".")[1].split("-")
    if len(val) == 1:
        exon_info1 = ''
        exon_info2 = ''
        if int(range_info[0]) in range(int(ref_range[0]), int(ref_range[1])):
            exon_number = val[0].split("-exon-")[1]
            exon_info1 = "exon" + exon_number
        if int(range_info[1]) in range(int(ref_range[0]), int(ref_range[1])):
            exon_number = val[0].split("-exon-")[1]
            exon_info2 = "exon" + exon_number
        if exon_info1 != '' or exon_info2 != '':
            if exon_info1 == '':
                exon_info1 = range_info[0]
            if exon_info2 == '':
                exon_info2 = range_info[1]
            exon_info = exon_info1 + "_" + exon_info2
    else:
        for i in val:
            exon_number = i.split("-exon-")[1]
            if junction_direction == "-":
                if old_exon_number <= int(exon_number):
                    exon_info = "exon" + exon_number + "_" + exon_info
                else:
                    exon_info += "exon" + exon_number + "_"
            if junction_direction == "+":
                if old_exon_number <= int(exon_number):
                    exon_info += "exon" + exon_number + "_"
                else:
                    exon_info = "exon" + exon_number + "_" + exon_info
            old_exon_number = int(exon_number)
        exon_info = exon_info.strip("_")
    if exon_info != '':
        fh_out.write(junction_name + "\t" + junction_length +
                     "\t" + gene_name + "_" + exon_info + "\n")
fh_out.close()
