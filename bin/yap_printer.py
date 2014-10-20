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
'''
Contains functions used for printing data in yap
'''


def wrap_line_to_width(line, width=80, file=None):
    '''
    Given a line and width,
    Returns list of lines wrapped at spaces in given line
    '''
    line = str(line)  # convert line to string first
    line_length = len(line)  # get line_length
    tab_count = line.count("\t")  # count tabs in line
    line_length = line_length + (7 * tab_count)  # adjust for tab character
    if line_length <= width:
        print_wrapped_line(line, file)
        return  # return to the caller
    wrapped_line = line[0:width]  # get wrapped line
    # find the split point
    actual_width = 0
    last_space_index = -1
    for counter in range(0, len(wrapped_line)):  # iterate over each character
        letter = wrapped_line[counter]  # get current character
        if letter == ' ':
            last_space_index = counter  # save space index
        if letter == '\t':
            actual_width = actual_width + 7  # tab adjust
        else:
            actual_width = actual_width + 1
        if actual_width > width:
            break  # stop if beyond allowed width
    end_index = -1  # end index of current line
    if counter < actual_width:  # if line had tabs
        if last_space_index > 0:
            end_index = last_space_index
        else:
            end_index = counter
    else:
        end_index = counter
    print_wrapped_line(line[0:end_index], file)  # print this line
    wrap_line_to_width(line[end_index:], width, file)  # print rest of the line
#-------------------------------------------------------------------------


def print_list(data, file=None, wrap_width=110,
               section_marker="::new_section", section_splitter='-'):

    # get the current terminal window size
    size_coord = []
    size_coord = os.popen('stty size', 'r').read().split()
    if len(size_coord) == 0:
        columns = 110  # if terminal is not available
    else:
        columns = int(size_coord[1])  # actual size of terminal window
    if columns > 0:
        wrap_width = columns - 1  # set wrap width
    if wrap_width < 80:
        wrap_width = 80  # set minimum size

    # print data
    for line in data:
        line = str(line)  # convert all lines to string
        # if line.find("command name=") >=0:#skip wrapping any commands
        #	print_wrapped_line(line,file)
        # continue#skip rest in this case
        # print new section if required
        if line.find(section_marker) >= 0:
            print_section(line, section_marker, section_splitter,
                          file, wrap_width)
            continue  # skip rest of the checks in this case
        else:
            print_wrapped_line(line, file)

        # print line if starts with tab
        # if line.find("\t") == 0: print "\n" +  line
        # lines = line.split('\n')#split lines starting new lines
        #map((lambda x: wrap_line_to_width(x, wrap_width,file)), lines)
#-------------------------------------------------------------------------


def print_section(line, section_marker, section_splitter, file, wrap_width):
    if line.find(section_marker + "_title") >= 0:
        title = line.replace(section_marker + "_title:", "")
        pad_count = get_character_num(title, wrap_width) / 2
        line = ("-" * pad_count + title + "-" * pad_count)
    else:
        line = (section_splitter * wrap_width)
    # print line
    print_wrapped_line(line, file)
#-------------------------------------------------------------------------


def print_wrapped_line(line, file):
    # print line on screen and in file
    print line
    if file is not None:
        print >>file, line
#-------------------------------------------------------------------------


def get_character_num(str, max_width):
    char_num = 0
    if max_width >= len(str):
        char_num = max_width - len(str)
    return char_num
#-------------------------------------------------------------------------
