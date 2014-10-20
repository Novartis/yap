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

class yap_conflict_check:

    """
    Provides methods to perform file-file, file-sample, file-group and
    sample-group comparisons and find conflicts.

    """

    def __init__(self, input_files):
        self.input_files = map(self.translate_path, input_files)
        self.filename_dict = \
            self.generate_filename_dict(self.input_files)

    def translate_path(self, path):
        """
        Given a path,
        Returns a path after expanding environment and user variables and
        relative paths to absolute path
        """
        path = os.path.expandvars(path)  # expand environment variables
        path = os.path.expanduser(path)  # expand user's home directory
        # don't convert to absolute if just filename
        if len(os.path.dirname(path)) == 0 and (path not in ['.', ".."]):
            return path
        path = os.path.abspath(path)  # convert relative path to absolute
        return path  # return output

    def translate_paired_end_paths(self, paired_end_files):
        '''
        Given a list of paired end files
        Returns a new list of paired end files with each file translated
        using translate path function
        '''
        if len(paired_end_files) <= 0:
            return []  # return empty o/p
        paired_end_files_out = []  # output variable
        for paired_list in paired_end_files:  # translate each paths
            paired_list_out = map(self.translate_path, paired_list)
            paired_end_files_out.append(paired_list)  # append to o/p
        return paired_end_files_out  # return output

    def get_paths(self, name):
        '''
        Given a name,
        Returns the list of paths matching to the key similar to the
        name
        '''
        if len(name) <= 0:
            return None  # return null for empty input
        # return if an exact match is found
        if name in self.filename_dict:
            return self.filename_dict[name]
        # return all values for a partial match
        matches = []
        for key in self.filename_dict:
            if key.find(name) == 0:
                new_paths = self.find_new_items(matches,
                                                self.filename_dict[key])
                # extend only if a unique match is found
                if len(new_paths) > 0:
                    matches.extend(new_paths)
        if len(matches) == 0:
            return None  # return null if no matches
        else:
            return matches  # return output

    def find_new_items(self, current_list, new_list):
        '''
        Given two lists,
        Returns items which are not available in current lists,
        Return empty list if no such items are found
        '''
        if len(current_list) == 0:
            return new_list  # all paths are new
        # select an items not in current list and return list
        return filter((lambda item: item not in current_list),
                      new_list)

    def validate_names_and_find_duplicates(self, names):
        '''
        Given list of filenames,
        Calls validate_names_and_find_duplicates_with_finder with
        get_paths as finder and returns the result
        '''
        return self.validate_names_and_find_duplicates_with_finder(
            names,
            self.get_paths)

    def validate_names_and_find_duplicates_with_finder(self, filenames,
                                                       finder):
        """
        Input:
                --filenames: a list of filenames occured in contaminant file
        Check if all filenames exist in input files name and
        there is no filename duplicate in filenames.
        Return values:
                --match_list:
                --error_list: all filenames which not exist in input files name
                --duplicate_dict: [key:value]
                        -key: filename which duplicate happens
                        -value: all path this filename occurs
        """
        match_list = []
        error_list = []
        duplicate_dict = {}
        # translate all filenames paths to complete paths
        filenames = map(self.strip_space_tab_newline, filenames)
        filenames = map(self.translate_path, filenames)
        for fn in filenames:
            if fn in self.input_files:
                # filename exist in self.input_files
                match_list.append(fn)
            else:
                # treat fn as basename
                paths = finder(fn)
                if paths is not None:
                    # basename exists
                    if len(paths) > 1:
                        # duplicate happens
                        duplicate_dict[fn] = paths
                    else:
                        # no duplicate
                        match_list.extend(paths)
                else:
                    # basename not exists
                    error_list.append(fn)
        return match_list, error_list, duplicate_dict

    def generate_filename_dict(self, paths):
        """
        Given a list of complete filepaths,
        Returns a dictionary, with keys as filenames and values as list of
        all paths that contain the corresponding key
        Invariant: Paths contain filenames complete with extension.

        """
        output = {}  # output variable
        if len(paths) <= 0:
            return output  # return empty output for empty input
        for path in paths:
            output[path] = [path]  # add each path as key also.
            basename = os.path.basename(path)  # get filename from path
            if len(basename) <= 0:
                continue  # skip if no filename in path
            # get name without extension
            basename_no_ext = os.path.splitext(basename)[0]
            # create a new entry if it does not exist, append otherwise
            if basename in output:
                output[basename].append(path)
            else:
                output[basename] = [path]
            # include a name with filename without extension also
            if len(basename_no_ext) <= 0:
                continue  # skip if name is exmpty
            if basename_no_ext != basename:  # add an entry for just filename
                if basename_no_ext in output:
                    output[basename_no_ext].append(path)
                else:
                    output[basename_no_ext] = [path]
        return output  # return dict

    def find_duplicates_in_list(self, input):
        """
        Given a list,
        Returns a dictionary of all duplicates in the list,
        Return empty dictionary if no duplicate entries are found.

        """
        output = {}  # output variable
        if len(input) <= 0:
            return output  # return empty output for empty input
        for item in input:
            if item not in output:  # check only if item not seen earlier
                item_count = input.count(item)  # count items
            # add to output if item occurs more than once in list
            if item_count > 1:
                output[item] = item_count
        return output

    def list_to_sentence(self, list):
        """
        Translate the given list to a string.

        """
        sentence = ""
        for i in range(0, len(list)):
            if i == len(list) - 1:
                sentence += "'" + list[i] + "'"
            else:
                sentence += "'" + list[i] + "' and "
        return sentence

    def strip_space_tab_newline(self, input):
        '''
        Given a string,
        Returns a string after removing starting and trailing spaces,
        tabs and new line character
        '''
        if len(input) <= 0:
            return ''  # empty o/p for empty i/p
        input = input.strip()
        input = input.strip('\n')
        input = input.strip('\t')
        return input
