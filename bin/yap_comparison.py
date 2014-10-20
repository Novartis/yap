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

class yap_comparison_parent:

    '''
    Provides methods required to validate cuffdiff, cuffcompare, cuffmerge
    file data
    '''

    def __init__(self, dict):
        self.file_reader = dict["file_reader"]
        self.cmd_file = dict["input_file"]
        self.samples = dict["samples"]
        self.validator = dict["conflict_validator"]
        self.paired_end_files =\
            self.validator.translate_paired_end_paths(dict["paired_end_data"])
        self.cmd_name = dict["cmd_name"]

    def validate_compare_file(self):
        '''
        validates comparison file data
        Returns formatted error messages if there are any discrepancies.
        '''
        file_data = []  # input file data
        if os.path.exists(self.cmd_file) == False:  # return if invalid file
            return {}, ["Error: File: " + self.cmd_file +
                        " does not exist. Please provide a valid filename."]
        file_data = self.file_reader(self.cmd_file)  # read input file
        groups = self.get_groups(file_data)  # get all groups in file
        # return error if no groups could be retereived
        if len(groups) == 0:
            return self.empty_file_error()
        # validate group length
        length_errors = self.validate_groups_length(groups)
        if len(length_errors) > 0:
            return groups, self.single_set_group_error(length_errors)
        # validate if file paths in each group can be found uniquely in
        # list of input files or in list of sample names
        groups, errors, name_conflicts =\
            self.validate_groups_and_expand_paths(groups, self.validator)
        # return if errors found in any of the groups
        if len(errors) > 0:
            return groups, self.errors_to_messages(errors)
        # return if duplicates found in any of the samples
        if len(name_conflicts):
            return groups, self.name_conflicts_to_messages(name_conflicts)
        # find duplicates - this iteration of duplicates will find file names or
        # sample names duplicates in each group
        duplicates = self.find_duplicates_in_groups(groups)
        # return error if any duplicates are found
        if len(duplicates) > 0:
            return groups, self.duplicates_to_messages(duplicates)
        # include paired end files and find duplicates in paired end file
        groups = self.include_paired_end_files_groups(
            groups, self.paired_end_files)
        # find duplicates - this iteration of duplicates will find only paired
        # end duplicates
        duplicates = self.find_duplicates_in_groups(groups)
        # return error if any duplicates are found
        if len(duplicates) > 0:
            return groups, self.paired_end_duplicates_to_messages(duplicates)
        # include sample names and find duplicates in sample names
        groups = self.include_sample_names_in_groups(groups, self.samples)
        # find duplicates - this iteration will find only sample name
        # duplicates
        duplicates = self.find_duplicates_in_groups(groups)
        # return groups and errors if any duplicates are found
        return groups, self.sample_name_duplicates_to_messages(duplicates)

#-------------------------------------------------------------------------
    def empty_file_error(self):
        '''
        Returns a formatted error to indicate an empty compare file
        '''
        return ["Error:" + self.cmd_name +
                " groups could not be retrieved from " +
                self.cmd_name + " comparison file: " + self.cmd_file +
                ". Please make sure that " + self.cmd_name +
                " comaprison  file is not empty and  is correctly formatted."]
#-------------------------------------------------------------------------

    def single_set_group_error(self, errors_in):
        '''
        Returns a formatted error to indicate single set groups
        '''
        errors_out = []
        for error in errors_in:
            errors_out.append(
                "Error: In " +
                self.cmd_name +
                " comparison file, group at line :" +
                str(error) +
                " contains only one set of file(s)." +
                "Please make sure that each line in " +
                self.cmd_name +
                " comparison file contains multiple file(s)/set of files(s)" +
                " separated by a space or tab.")
        return errors_out
#-------------------------------------------------------------------------

    def validate_groups_length(self, groups):
        '''
        Given a dictionary of groups,
        Returns the group names which contain only one set
        '''
        if len(groups) <= 0:
            return {}  # return empty o/p for empty i/p
        errors = []  # output variable
        for group in groups:  # validate each group
            if len(groups[group]) <= 1:
                errors.append(group)
        return errors  # return output

#-------------------------------------------------------------------------

    def validate_groups_and_expand_paths(self, groups_in, validator):
        '''
        Given a dictionary of groups and a validator object,
        Converts file paths in each group to complete paths,
        Returns the new groups dictionary and name conflict errors, if any
        '''
        if len(groups_in) <= 0:
            return {}  # return empty o/p for empty i/p
        groups_out = {}
        errors = {}
        name_conflicts = {}  # output variables
        for group in groups_in:  # validate each group
            group_out, group_errors, group_name_conflicts =\
                self.validate_group_and_expand_paths(groups_in[group],
                                                     validator)
            # add conflicts,errors and groups
            if len(group_name_conflicts) > 0:
                # add name conflicts
                name_conflicts[group] = group_name_conflicts
            elif len(group_errors) > 0:
                errors[group] = group_errors  # add errors
            else:
                groups_out[group] = group_out  # append expanded path to output
        return groups_out, errors, name_conflicts  # return output

#-------------------------------------------------------------------------
    def validate_group_and_expand_paths(self, group_in, validator):
        '''
        Given a list of all sides of a group and a validator object,
        Converts file paths in each side of group to complete paths,
        Returns new group list and name conflicts error list, if any
        '''
        if len(group_in) <= 0:
            return {}  # return empty o/p for empty i/p
        group_out = []
        errors = []
        name_conflicts = []  # output variables
        # validate each side using path and samples finder
        for side in group_in:
            side_out, side_errors, side_name_conflicts =\
                validator.validate_names_and_find_duplicates_with_finder(side,
                                                                         self.get_paths_or_samples)
            # add duplicates and sides
            if len(side_name_conflicts) > 0:
                name_conflicts.append(side_name_conflicts)  # add duplicates
            elif len(side_errors) > 0:
                errors.extend(side_errors)  # add errors
            else:
                group_out.append(side_out)  # append expanded path to output
        return group_out, errors, name_conflicts  # return output
#-------------------------------------------------------------------------

    def get_paths_or_samples(self, name):
        '''
        Given a name,
        Returns list of matching results from list of input files and
        list of sample names
        '''
        if len(name) <= 0:
            return []  # return empty o/p for empty i/p
        matches = []  # output variable
        matches = self.validator.get_paths(name)  # check if name is a filename
        if matches is not None:
            return matches  # return if a filename match is found
        # if no filename match is found, check if name matches any of samples
        matches = filter((lambda x: x.find(name) == 0), self.samples.keys())
        if len(matches) == 0:
            return None  # return none if no matches
        else:
            return matches  # return matches otherwise
#-------------------------------------------------------------------------

    def errors_to_messages(self, errors):
        '''
        Given a dictionary of errors,
        Returns list of formatted error messages for each error.
        '''
        if len(errors) <= 0:
            return []  # return empty output
        messages = []  # output variable
        for group in errors:
            messages.append(
                "Error: In " +
                self.cmd_name +
                " comparison file, " +
                "files: " +
                self.list_to_sentence(
                    errors[group]) +
                " specified in " +
                "line number: " +
                str(group) +
                " are not present in the list of" +
                " input files or regroup samples for this " +
                "workflow. Please provide a" +
                " valid list of files or samples for this group.")
        return messages  # return output

#-------------------------------------------------------------------------
    def name_conflicts_to_messages(self, dict):
        '''
        Given a dictionary of duplicates,
        Returns list of formatted error messages for each duplicate.
        '''
        if len(dict) <= 0:
            return []  # return empty output for empty input
        messages = []  # output variable
        for group in dict:
            for name_dict in dict[group]:
                for name in name_dict:
                    messages.append(
                        "Error: In " +
                        self.cmd_name +
                        " comparison " +
                        "file, basename : '" +
                        name +
                        "' specified in line number: " +
                        str(group) +
                        " was found in: " +
                        self.list_to_sentence(
                            name_dict[name]) +
                        ". Please provide a complete path or name" +
                        " to avoid ambiguity.")
        return messages
#-------------------------------------------------------------------------

    def find_duplicates_in_groups(self, groups):
        '''
        Given a dictionary of groups,
        Return dictionary of groups which have duplicates, empty dict otherwise
        '''
        if len(groups) <= 0:
            return {}  # empty o/p for empty i/p
        duplicates = {}  # output variable
        for group in groups:
            group_duplicates = self.find_duplicates_in_group(groups[group])
            if len(group_duplicates) > 0:
                # add duplicates to the o/p
                duplicates[group] = group_duplicates
        return duplicates  # return o/p
#-------------------------------------------------------------------------

    def find_duplicates_in_group(self, group):
        '''
        Given a group as a list,
        Returns a dictionary containing duplicate names and their count,
        Return empty dict if no duplicates found
        '''
        if len(group) <= 0:
            return {}  # empty o/p for empty i/p
        duplicates = {}  # output dictionary to store duplicates
        names_table = {}  # instantiate a names table to keep track of all
        # items in current group
        for set in group:  # each set can contain multiple names/files
            for name in set:
                if name in names_table:
                    names_table[name] = names_table[name] + 1
                    # add repetition count
                    duplicates[name] = names_table[name]
                else:
                    names_table[name] = 1  # add first occurrence of a name
        return duplicates  # return o/p
#-------------------------------------------------------------------------

    def duplicates_to_messages(self, duplicates):
        '''
        Given a dictionary of duplicates,
        Returns list of formatted error messages
        '''
        if len(duplicates) <= 0:
            return []  # return empty o/p
        messages = []  # output variable
        for group in duplicates:
            for name in duplicates[group]:
                messages.append("Error: In " +
                                self.cmd_name +
                                " comparison file, " +
                                self.sample_or_file(name) +
                                " has been specified " +
                                str(duplicates[group][name]) +
                                " times " +
                                "in line number :" +
                                str(group) +
                                " Please " +
                                "provide unique values to avoid ambiguity.")
        return messages  # return output
#-------------------------------------------------------------------------

    def sample_or_file(self, name):
        '''
        Given a name,
        Determines if it is a sample or file name and
        Returns formatted name accordingly
        '''
        if len(name) <= 0:
            return ''  # return empty o/p for empty i/p
        if len(os.path.dirname(name)) == 0:
            return "Sample : " + name
        else:
            return "File : " + name
#-------------------------------------------------------------------------

    def include_paired_end_files_groups(self, groups_in, paired_end_files):
        '''
        Given a groups dict and a list of paired end files,
        Returns an updated groups dictionary with all paired end files included
        '''
        if len(groups_in) <= 0:
            return {}  # empty o/p for empty i/p
        if len(paired_end_files) <= 0:
            return groups_in  # no updates required
        groups_out = {}  # output variable
        for group in groups_in:  # add paired end files for each group
            groups_out[group] = \
                self.include_paired_end_files_group(groups_in[group],
                                                    paired_end_files)
        return groups_out  # return output
#-------------------------------------------------------------------------

    def include_paired_end_files_group(self, group_in, paired_end_files):
        '''
        Given a group as list and a list of paired end files,
        Returns updated group  with all paired end files included
        '''
        if len(group_in) <= 0:
            return []  # empty o/p for empty i/p
        group_out = []  # output variable
        group_seen = []  # names in all the sets seen so far
        for set in group_in:  # each set can contain multiple file names
            # copy current set
            set_out = copy.deepcopy(set)
            for name_counter in range(0, len(set)):
                name = set[name_counter]  # loop over each name
                # check if name has a pair
                for paired_list in paired_end_files:
                    if name in paired_list:
                        # dont't consider if name is already seen in group
                        if name in group_seen:
                            break
                        # don't consider if name is already duplicate
                        if set_out.count(name) > 1:
                            break
                        # remove file from output if it has a pair
                        set_out.remove(name)
                        # add the pair to the output
                        set_out.append(paired_list[0])
                        break  # break current loop
            group_seen.extend(set_out)  # add current set elements as seen
            group_out.append(set_out)  # append updated set to the group
        return group_out  # return o/p
#-------------------------------------------------------------------------

    def paired_end_duplicates_to_messages(self, duplicates):
        '''
        Given a dictionary of paired end duplicates,
        Returns list of formatted error messages
        '''
        if len(duplicates) <= 0:
            return []  # return empty o/p
        messages = []  # output variable
        for group in duplicates:
            for name in duplicates[group]:
                messages.append("Error: In " + self.cmd_name + " comparison " +
                                "file, Paired-end " +
                                " duplicate found for file : " + name +
                                " in line number :" + str(group) +
                                ". Please provide the complete path for the" +
                                " first pair at most once in one line to " +
                                "avoid ambiguity.")
        return messages  # return output
#-------------------------------------------------------------------------

    def include_sample_names_in_groups(self, groups_in, samples):
        '''
        Given groups and samples dictionary,
        Returns updated groups dictionary with file names replaced with sample
        names
        '''
        if len(groups_in) <= 0:
            return {}  # empty o/p for empty i/p
        if len(samples) <= 0:
            return groups_in  # no updates required
        groups_out = {}  # output variable
        for group in groups_in:
            groups_out[group] = self.include_sample_names_in_group(
                groups_in[group],
                samples)
        return groups_out  # return o/p
#-------------------------------------------------------------------------

    def include_sample_names_in_group(self, group_in, samples):
        '''
        Given a group as list and samples dictionary,
        Returns updated group with all file names replaced with sample names
        '''
        if len(group_in) <= 0:
            return []  # empty o/p for empty i/p
        group_out = []  # output variable
        for set in group_in:  # each set can contain multiple file names
            set_out = copy.deepcopy(set)  # copy current set
            # will contain unique samples found in current set
            samples_found = []
            for name_counter in range(0, len(set)):
                name = set[name_counter]  # loop over each name
                for sample in samples:  # check if name is part of a sample
                    if name in samples[sample]:
                        # remove file from output as it will be replaced
                        # with sample name
                        set_out.remove(name)
                        # add to samples found list if not already there
                        if sample not in samples_found:
                            samples_found.append(sample)
            # add sample names found in above search to the set
            set_out.extend(samples_found)
            group_out.append(set_out)  # append updated set to the group
        return group_out  # return o/p
#-------------------------------------------------------------------------

    def sample_name_duplicates_to_messages(self, duplicates):
        '''
        Given a dictionary of sample name duplicates,
        Returns list of formatted error messages
        '''

        if len(duplicates) <= 0:
            return []  # return empty o/p
        messages = []  # output variable
        for group in duplicates:
            for name in duplicates[group]:
                messages.append(
                    "Error: In " +
                    self.cmd_name +
                    " comparison file," +
                    " Sample name: " +
                    name +
                    " or files from this sample" +
                    " have been specified multiple times in line " +
                    "number: " +
                    str(group) +
                    ". Please specify sample name at most once " +
                    "per line to avoid ambiguity.")
        return messages  # return output
#-------------------------------------------------------------------------

    def list_to_sentence(self, list):
        """
        Translate the given list to a string.
        """
        sentence = "\n"
        for i in range(0, len(list)):
            if i == len(list) - 1:
                sentence += "'" + list[i] + "'"
            else:
                sentence += "'" + list[i] + "'\n"
        return sentence
#-------------------------------------------------------------------------
