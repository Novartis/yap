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
import os
import gzip
import bz2
import sys
"""
functions : create_openfile_handler, read_file, read_file_chunk
"""


def create_openfile_handler(input_file):
    """ Creates an open file handler given an input file.
        Checks if file exists first. Uncompresses file if necessary. """

    path_name, file_name = os.path.split(input_file)
    filebase_name, file_ext = os.path.splitext(file_name)
    file_handler = ''
    if file_ext == '.bz2':
        try:
            file_handler = bz2.BZ2File(input_file, 'rb')

        except IOError as xxx_todo_changeme:
            (errno, strerror) = xxx_todo_changeme.args
            print "Error:  while Opening the file : ", input_file
            print "        I/O error({0}): {1}".format(errno, strerror)

    elif file_ext == '.gz':
        try:
            file_handler = gzip.open(input_file, 'rb')

        except IOError as xxx_todo_changeme1:
            (errno, strerror) = xxx_todo_changeme1.args
            print "Error:  while Opening the file : ", input_file
            print "        I/O error({0}): {1}".format(errno, strerror)
    else:
        try:
            file_handler = open(input_file, 'r')

        except IOError as xxx_todo_changeme2:
            (errno, strerror) = xxx_todo_changeme2.args
            print "Error:  while Opening the file : ", input_file
            print "        I/O error({0}): {1}".format(errno, strerror)

    return file_handler


def read_file(input_file):
    """ Returns a list of file contents, else throws an exception. """
    
    seqs_arr = []
    file_handler = ''
    try:
        file_handler = create_openfile_handler(input_file)
        if file_handler != '':
            seqs_arr = file_handler.readlines()
            file_handler.close()
    except Exception as e:
        print "Error:  Reading the file : ", input_file, e
    return seqs_arr


def read_file_chunk(input_file, file_begin_pos, file_end_pos):
    """ Reads a file chunk, given start and end position of the file chunk. """

    seqs_arr = []
    file_handler = ''
    try:
        file_handler = create_openfile_handler(input_file)
        if file_handler != '':
            file_handler.seek(0)
            file_handler.seek(file_begin_pos)
            seqs_arr = file_handler.read(file_end_pos - file_begin_pos)
            file_handler.close()
    except Exception as e:
        print "Error:  while Reading from the file : ", input_file, e
    return seqs_arr


def write_data(data, output_filename):
    """ Takes a list and filename as input. Writes data to that file. """

    fw = open(output_filename, 'a')
    fw.writelines(data)
    fw.close()
