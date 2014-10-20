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


class yap_cuffdiff_check(yap_comparison_parent):

    '''
    Provides methods required to validate cuffdiff file data
    and return appropriate error messages.
    '''

    def __init__(self, dict):
        yap_comparison_parent.__init__(self, dict)

#-------------------------------------------------------------------------

    def get_sets(self, cuff_diff_data):
        '''
        Given data from cuff diff file,
        Returns the dictionary of all sets from this file
        '''
        if len(cuff_diff_data) == 0:
            return {}  # return empty o/p
        sets = {}  # output variable
        for i in range(0, len(cuff_diff_data)):
            line = cuff_diff_data[i].strip('\n')  # remove new line character
            if len(line) > 0:  # continue only if line has some content
                # search current line for occurences of whitespace characters
                matchobj = re.match(r'\s*(\S*)(\s*).*', line, re.M | re.I)
                if matchobj is None:
                    return {}  # return empty if any invalid line
                delimiter = matchobj.groups()[1]  # get separator character
                if len(delimiter) <= 0:
                    return {}  # return empty if any invalid line
                # get sets of group in comma separated(or single member) format
                all_sets_comma = line.split(delimiter)
                # create an entry for each set in o/p
                for set in all_sets_comma:
                    set = set.split(',')  # convert to a list
                    print set
                    if set in sets:
                        sets[set] = sets[set].append(i)
                    else:
                        sets[set] = [i]
        return sets  # return output
#-------------------------------------------------------------------------

    def get_groups(self, cuff_diff_data):
        '''
        Given data from cuff diff file,
        Returns the dictionary of groups from this file
        '''
        if len(cuff_diff_data) == 0:
            return {}  # return empty o/p
        groups = {}  # output variable
        for i in range(0, len(cuff_diff_data)):
            line = cuff_diff_data[i].strip('\n')  # remove new line character
            line = line.replace(", ", ",")  # remove spaces from same set
            line = line.replace(",\t", ",")  # remove tabs from same set
            if len(line) > 0:  # continue only if line has some content
                # search current line for occurences of whitespace characters
                matchobj = re.match(r'\s*(\S*)(\s*).*', line, re.M | re.I)
                if matchobj is None:
                    return {}  # return empty if any invalid line
                delimiter = matchobj.groups()[1]  # get separator character
                if len(delimiter) <= 0:
                    return {}  # return empty if any invalid line
                # get sets of group in comma separated(or single member) format
                all_sets_comma = line.split(delimiter)
                # get all sets with each set as list
                all_sets = map(
                    (lambda x: string.split(x, ",")), all_sets_comma)
                # remove empty elements
                all_sets = map((lambda x: self.remove_empty_elements(x)),
                               all_sets)
                # remove empty sets
                all_sets = filter((lambda x: len(x) > 0), all_sets)
                # append list of sets to group only if non-empty
                if len(all_sets) > 0:
                    groups[i + 1] = all_sets
        return groups  # return output
#-------------------------------------------------------------------------

    def remove_empty_elements(self, input):
        '''
        Given a list of string,
        Returns a new list of non-empty elements
        after removing leading and trailing tabs and spaces
        '''
        return filter((lambda x: len(x) > 0), input)
#-------------------------------------------------------------------------
