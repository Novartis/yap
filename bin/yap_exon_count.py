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
    Options: -f [PercentOverlap]
           -exon_coordinates_file [ucsc coordinate bed file] 
           -exon_CoordToNumber_file [filename] 
           -i [inputfilename or for stdin '-'] 
           -o [output filename]   """

def print_usage():
    print "DESCRIPTION: Counts number of exons .\n"
    print "USAGE:"
    print sys.argv[0] + " -f [PercentOverlap]  -exon_coordinates_file [ucsc coordinate bed file] -exon_CoordToNumber_file [filename] -i [inputfilename or for stdin '-'] -o [output filename]\n"
    exit()

if len(sys.argv) < 11:
    print_usage()
for i in range(0, len(sys.argv)):
    if sys.argv[i] == "-input_type":
        file_type = sys.argv[i + 1]
    if sys.argv[i] == "-f":
        percent_overlap = sys.argv[i + 1]
    if sys.argv[i] == "-exon_coordinates_file":
        exon_coordinates_file = sys.argv[i + 1]
    if sys.argv[i] == "-exon_CoordToNumber_file":
        exon_coordinate_number_file = sys.argv[i + 1]
    if sys.argv[i] == '-i':
        input_file = sys.argv[i + 1]
    if sys.argv[i] == '-o':
        output_file = sys.argv[i + 1]
final_output_file = output_file + "_exoncount.txt"
file_type = ''


def check_file_type(input_str, input_file):
    file_type = ''
    def get_file_type(input_str):
        if '\0' in input_str:
            # file is in binary fomat
            file_type = 'bam'
        else:
            file_type = "sam"
        return file_type
    if input_file != '':
        input_str = open(input_file, 'r').read(200)
        file_type = get_file_type(input_str)
    if input_file == '':
        file_type = get_file_type(input_str)
    return file_type

input_str = ''
if input_file == "-":
    input_str = sys.stdin.read()
    file_type = check_file_type(input_str, '')
    if file_type == "sam":
        exon_cmd = " awk '$6!~/[N]/' | samtools view -bS - | bamToBed -i - | intersectBed -a - -b " + exon_coordinates_file + " -f " + percent_overlap + \
            " | sort +0 -1 +1n -2 +2n -3 +3 -4 - | uniq -c | sed 's/^      //' | sed 's/ /\t/' | cut -f2-7 | coverageBed -a - -b " + \
            exon_coordinates_file + " > " + final_output_file

    if file_type == "bam":
        exon_cmd = "samtools view -h - | awk '$6!~/[N]/' | samtools view -bS - | bamToBed -i - | intersectBed -a - -b " + exon_coordinates_file + " -f " + percent_overlap + \
            " | sort +0 -1 +1n -2 +2n -3 +3 -4 - | uniq -c | sed 's/^      //' | sed 's/ /\t/' | cut -f2-7 | coverageBed -a - -b " + \
            exon_coordinates_file + " > " + final_output_file

    if file_type == 'sam' or file_type == 'bam':
        try:
            prun = Popen(exon_cmd, stdin=PIPE, stdout=PIPE, shell='False')
            out, err = prun.communicate(input_str)
            input_str = ''
        except:
            print "Error: Exon counts failed", "\ncommand=", exon_cmd
            print "Please check if you provided the correct the alignment file type[sam or bam] "

else:
    file_type = check_file_type('', input_file)
    if file_type == 'sam':
        exon_cmd = "cat " + input_file + " | awk '$6!~/[N]/' | samtools view -bS - | bamToBed -i - | intersectBed -a - -b " + exon_coordinates_file + " -f " + percent_overlap + \
            " | sort +0 -1 +1n -2 +2n -3 +3 -4 - | uniq -c | sed 's/^      //' | sed 's/ /\t/' | cut -f2-7 | coverageBed -a - -b " + \
            exon_coordinates_file + " > " + final_output_file

    if file_type == 'bam':
        exon_cmd = "cat " + input_file + "|samtools view -h - | awk '$6!~/[N]/' | samtools view -bS - | bamToBed -i - | intersectBed -a - -b " + exon_coordinates_file + " -f " + percent_overlap + \
            "  | sort +0 -1 +1n -2 +2n -3 +3 -4 - | uniq -c | sed 's/^      //' | sed 's/ /\t/' | cut -f2-7 | coverageBed -a - -b " + \
            exon_coordinates_file + " > " + final_output_file

    if file_type == 'sam' or file_type == 'bam':
        try:
            prun = Popen(exon_cmd, stdin=PIPE, stdout=PIPE, shell='False')
            prun.communicate()
        except:
            print "Error: Exon counts failed", "\ncommand=", exon_cmd
            print "Please check if you provided the correct the alignment file type[sam or bam] "

# convert exon coordinate information to exon number
# created just for easy readability, the file with exon coordinate will
# still be kept intact
try:
    fh = open(exon_coordinate_number_file, 'rb')
except IOError as e:
    print "Error while opening the ", exon_coordinate_number_file
    print e
exon_number_arr = fh.readlines()
fh.close()

line = ''
temp_arr = []
exon_coord_num_dict = {}
for line in exon_number_arr:
    temp_arr = line.strip("\n").split("\t")
    if len(temp_arr) < 5:
        print "Not enough fields in exon number file"
        print "please refer to the documentation"
        exit()
    else:
        exon_coordinate = temp_arr[3]
        exon_number = temp_arr[4]
        if len(temp_arr) > 5:
            gene_coordinate = temp_arr[5]
        else:
            gene_coordinate = ''
        exon_coord_num_dict[exon_coordinate] = [exon_number, gene_coordinate]
try:
    fh = open(final_output_file, 'rb')
except IOError as e:
    print "Error while opening the ", final_output_file
    print e
    exit()
line = ''
temp_arr = []
fh_out = open(output_file + "_exoncount_summary.txt", 'wb')

for line in fh:
    temp_arr = line.strip("\n").split("\t")
    if len(temp_arr) != 8:
        print "Not enough fields in exon counts file"
        exit()
    else:
        exon_coordinate_temp = temp_arr[3]
        counts = temp_arr[4]
        if (exon_coordinate_temp in exon_coord_num_dict):
            exon_number = exon_coord_num_dict[exon_coordinate_temp][0]
            gene_coordinate = exon_coord_num_dict[exon_coordinate_temp][1]
            if gene_coordinate == '':
                fh_out.write(exon_number + "\t" + counts + "\n")
            else:
                fh_out.write(
                    exon_number +
                    "\t" +
                    gene_coordinate +
                    "\t" +
                    counts +
                    "\n")
fh.close()
fh_out.close()
