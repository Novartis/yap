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

import re
import string
import copy
import os
from yap_comparison import *


class yap_cuff_merge_compare_check(yap_comparison_parent):

    '''
    Provides methods required to validate comparison files used in
    yap post processing
    '''

    def __init__(self, dict):
        yap_comparison_parent.__init__(self, dict)

#-------------------------------------------------------------------------

    def get_groups(self, file_data):
        '''
        Given data from input file,
        Returns the dictionary of groups from this file
        '''
        if len(file_data) == 0:
            return {}  # return empty o/p
        groups = {}  # output variable
        for i in range(0, len(file_data)):
            line = file_data[i].strip('\n')  # remove new line character
            if len(line) > 0:  # continue only if line has some content
                # search current line for occurences of whitespace characters
                matchobj = re.match(r'\s*(\S*)(\s*).*', line, re.M | re.I)
                if matchobj is None:
                    return {}  # return empty if any invalid line
                delimiter = matchobj.groups()[1]  # get separator character
                if len(delimiter) <= 0:
                    return {}  # return empty if any invalid line
                all_sets = line.split(delimiter)  # get sets of group
                # remove empty elements
                all_sets = filter((lambda x: len(x) > 0), all_sets)
                # convert each elements to single element list to get common
                # formatting
                all_sets = map((lambda x: [x]), all_sets)
                # append list of sets to group only if non-empty
                if len(all_sets) > 0:
                    groups[i + 1] = all_sets
        return groups  # return output
#-------------------------------------------------------------------------
