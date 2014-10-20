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

import string
import copy
import re
import os


class yap_regroup_check:

    '''
    Provides methods required to validate yap regroup file
    and return appropriate error messages
    Invariant1: An existing regroup file is provided
    Invariant2: A non-null conflict validator and file reader is provided
    Invariant3: A valid paired end data is provided
    '''

    def __init__(self, dict):
        # initialize all regroup global variables
        self.regroup_file = dict["regroup_file"]
        self.validator = dict["conflict_validator"]
        self.file_reader = dict["file_reader"]
        self.paired_end_files =\
            self.validator.translate_paired_end_paths(dict["paired_end_data"])
#-------------------------------------------------------------------------

    def validate_regroup_file(self):
        '''
        Validates the given regroup file using given validator object,
        Returns the samples dictionary and  error list if there are
        any discrepancies.
        '''
        samples = {}  # samples data in regroup file
        file_data = []  # input file data
        # return if invalid file
        if os.path.exists(self.regroup_file) == False:
            return {}, ["Error: Regroup file: " + self.regroup_file +
                        " does not exist. Please provide a valid filename."]
        file_data = self.file_reader(self.regroup_file)  # read regroup file
        samples, duplicate = self.get_sample_names(
            file_data)  # get all sample names
        # return if any errors or duplicates found
        if len(duplicate) > 0:
            return samples, self.duplicate_to_messages(duplicate)
        # return error if no samples could be retereived
        if len(samples) == 0:
            return {}, self.empty_regroup_file_error()
        # validate if file paths in each sample can be found uniquely in
        # list of input file
        samples, errors, name_conflicts =\
            self.validate_samples_and_expand_paths(samples, self.validator)
        # return if errors found in any of the samples
        if len(errors) > 0:
            return samples, self.errors_to_messages(errors)
        # return if duplicates found in any of the samples
        if len(name_conflicts):
            return samples, self.name_conflicts_to_messages(name_conflicts)
        # check if any sample has multiple files with same name
        # this check will find any duplicates except paired end duplicates
        duplicates = self.find_duplicates_in_samples(
            samples,
            self.validator.find_duplicates_in_list)
        if len(duplicates) > 0:  # return if any sample has same filenames
            return samples, self.duplicates_in_samples_to_messages(duplicates)
        # include paired end files in samples
        samples = self.include_paired_end_files(samples, self.paired_end_files)
        # check if any sample has multiple files with same name
        # this check will now find paired end duplicates
        duplicates = self.find_duplicates_in_samples(
            samples,
            self.validator.find_duplicates_in_list)
        if len(duplicates) > 0:  # return if any sample has same filenames
            return samples, self.paired_end_duplicates_to_messages(duplicates)
        # check if any file appears in more than one sample
        duplicates = self.find_duplicate_files_accross_samples(samples)
        # return if any errors or duplicates found
        if len(duplicates) > 0:
            return samples, self.duplicates_to_messages(duplicates)
        # check if any sample name is same as filename and does not contain that
        # file
        conflicts = self.conflicting_sample_file_names(
            samples,
            self.validator.filename_dict)
        # return output data
        return samples, self.conflicts_to_messages(conflicts)
#-------------------------------------------------------------------------

    def extend_regroup_samples(self, samples, input_files_list):
        '''
        Given samples dictionary, list of input files with tag names
        Returns updated samples dictionary with new samples added for the files
        which are not currently in regroup samples
        '''
        if len(samples) == 0:
            return {}  # empty o/p for empty i/p
        # get dictionary of files which are already in samples
        samples_files = self.samples_to_list(samples)
        for item in input_files_list:
            if item[0] not in samples_files:  # check only for first file
                # add file/pair tag as key and path(s) as list
                samples[item[2]] = [item[0]]  # add first file
                # add second file of pair if exists
                #if len(item[1]) > 0: samples[item[2]].append(item[1])
#-------------------------------------------------------------------------

    def samples_to_list(self, samples):
        '''
        Given a list of samples,
        Returns a dictionary of all files in samples with file path as key
        '''
        if len(samples) == 0:
            return {}  # empty o/p for empty i/p
        output = {}  # output variable
        for sample in samples:
            for file in samples[sample]:
                if file not in output:
                    output[file] = ""  # add paths only as keys
        return output  # return output dictionary

#-------------------------------------------------------------------------
    def empty_regroup_file_error(self):
        '''
        Returns a formatted error to indicate an empty regroup file
        '''
        return ["Error: Sample Information could not be retrieved from " +
                "regroup file: " + self.regroup_file +
                ". Please make sure that regroup file is not empty and is" +
                " correctly formatted."]

#-------------------------------------------------------------------------
    def find_duplicates_in_samples(self, samples, duplicate_finder):
        '''
        Given a list of samples,
        Returns a dictionary of duplicates
        '''
        if len(samples) <= 0:
            return {}  # return empty output
        duplicates_out = {}  # output duplicate dictionary
        for sample in samples:
            duplicates = duplicate_finder(samples[sample])
            if len(duplicates) > 0:
                duplicates_out[sample] = duplicates
        return duplicates_out  # return output
#-------------------------------------------------------------------------

    def name_conflicts_to_messages(self, dict):
        '''
        Given a dictionary of samples with files with duplicate basename ,
        Returns list of formatted error messages for each sample.
        '''
        if len(dict) <= 0:
            return []  # return empty output for empty input
        messages = []  # output variable
        for sample in dict:
            for name in dict[sample]:
                messages.append("Error: In regroup file, basename : '" +
                                name + "' specified in Sample: " + sample +
                                " was found in files: " +
                                self.list_to_sentence(dict[sample][name]) +
                                ". Please provide the complete path for the" +
                                " first pair to avoid ambiguity.")
        return messages
#-------------------------------------------------------------------------

    def paired_end_duplicates_to_messages(self, dict):
        '''
        Given a dictionary of samples with duplicate files,
        Returns list of formatted error messages for each sample.
        '''
        if len(dict) <= 0:
            return []  # return empty output for empty input
        messages = []  # output variable
        for sample in dict:
            for file in dict[sample]:
                messages.append(
                    "Error: In regroup file, Paired-end duplicate " +
                    "found for File:  " +
                    file +
                    " in Sample:" +
                    sample +
                    ". Please specify only one file of paired end " +
                    "data for regrouping. YAP will automatically " +
                    "include other files of the pair in this sample" +
                    " while regrouping.")
        return messages
#-------------------------------------------------------------------------

    def duplicates_in_samples_to_messages(self, dict):
        '''
        Given a dictionary of samples with duplicate files,
        Returns list of formatted error messages for each sample.
        '''
        if len(dict) <= 0:
            return []  # return empty output for empty input
        messages = []  # output variable
        for sample in dict:
            for file in dict[sample]:
                messages.append("Error: In regroup file, File:  " +
                                file +
                                " is mentioned " +
                                str(dict[sample][file]) +
                                " times in Sample:" +
                                sample +
                                ". Please remove duplicate entires from this " +
                                "sample.")
        return messages
#-------------------------------------------------------------------------

    def errors_to_messages(self, errors):
        '''
        Given a dictionary of errors,
        Returns list of formatted error messages for each error.
        '''
        if len(errors) <= 0:
            return []  # return empty output
        messages = []  # output variable
        for error in errors:
            messages.append(
                "Error: In regroup file, files: " +
                self.list_to_sentence(
                    errors[error]) +
                " specified in " +
                "Sample: " +
                error +
                " are not present in the list of" +
                " input files for this workflow. Please provide" +
                " valid list of files for this sample.")
        return messages  # return output

#-------------------------------------------------------------------------
    def duplicates_to_messages(self, dict):
        '''
        Given a dictionary of duplicates,
        Returns list of formatted error messages for each duplicate.
        '''
        if len(dict) <= 0:
            return []  # return empty output for empty input
        messages = []  # output variable
        for duplicate in dict:
            messages.append("Error: In regroup file, same files " +
                            "or different files from same pair " +
                            "are part of samples: \n'" + duplicate + "'" +
                            self.list_to_sentence(dict[duplicate]) +
                            ". Please remove duplicate files from these " +
                            "samples.")
        return messages
#-------------------------------------------------------------------------

    def duplicate_to_messages(self, duplicate):
        '''
        Given a duplicate sample name,
        Returns a list of error messages for this duplicate
        '''
        return [
            "Error: In regroup file, Sample name: " +
            duplicate +
            " is used" +
            " more than once. Please specify unique sample names."]
#-------------------------------------------------------------------------

    def conflicts_to_messages(self, conflicts):
        '''
        Given a list of conflicting sample names,
        Returns list of formatted error messages for each conflict
        '''
        if len(conflicts) <= 0:
            return[]  # return empty output
        return map(
            (lambda conflict: "Error: In regroup file: Sample name: " +
             conflict +
             " matches a filename in the input files. Please change the" +
             " sample name to avoid any ambiguity."),
            conflicts)

#-------------------------------------------------------------------------
    def conflicting_sample_file_names(self, samples, filenames_dict):
        '''
        Given dictionary of samples and filenames(without paths),
        Return list of sample names which match file names and do not contain
        these file names.
        '''
        if len(samples) <= 0 or len(filenames_dict) <= 0:
            return {}  # return empty dict
        conflicts = []  # store conficts output
        for sample in samples:
            if sample in filenames_dict:  # if sample name is also a filename
                    # add sample to output if matching file is not part of
                    # sample
                if filenames_dict[sample][0] not in samples[sample]:
                    conflicts.append(sample)
        return conflicts  # return output
#-------------------------------------------------------------------------

    def find_duplicate_files_accross_samples(self, samples):
        '''
        Given a samples dictionary,
        Returns a dictionary of samples which contain same files
        '''
        if len(samples) <= 0:
            return {}  # return empty for empty input
        duplicates = {}  # output dictionary
        # compare files in each sample with files of every other sample
        for first_counter in range(0, len(samples)):
            key1 = samples.keys()[first_counter]
            for second_counter in range(first_counter + 1, len(samples)):
                key2 = samples.keys()[second_counter]
                # find duplicates among current sample pairs
                duplicate = filter((lambda file: file in samples[key2]),
                                   samples[key1])
                if len(duplicate) > 0:  # add duplicate information to output
                    if key1 in duplicates:
                        duplicates[key1].append(key2)
                    else:
                        duplicates[key1] = [key2]
        return duplicates  # return output
#-------------------------------------------------------------------------

    def include_paired_end_files(self, samples_in, paired_end_files):
        '''
        Given a samples dictionary and list of paired end data,
        Returns sample dictionary with all files of a pair included in sample if
        one of the files was originally included
        '''
        if len(paired_end_files) <= 0:
            return samples_in  # return input as it is
        samples_out = {}  # ouptut samples dictionary
        for sample in samples_in:
            # add current sample to output
            samples_out[sample] = copy.deepcopy(samples_in[sample])
            for file_counter in range(0, len(samples_in[sample])):
                file = samples_in[sample][file_counter]  # loop over each file
                # check if file has a pair
                for paired_list in paired_end_files:
                    if file in paired_list:
                        # break if file is already a duplicate in sample
                        if samples_out[sample].count(file) > 1:
                            break
                        # remove file from output if it has a pair
                        samples_out[sample].remove(file)
                        # add the first file of pair to the output
                        samples_out[sample].append(paired_list[0])
                        break  # break current loop
        return samples_out  # return output

#-------------------------------------------------------------------------
    def validate_samples_and_expand_paths(self, samples_in, validator):
        '''
        Given a samples dictionary and a validator object,
        Validates if all files in all samples are there in input files
        Returns errors list if any error else returns empty list
        '''
        errors_out = {}  # output error dictionary
        samples_out = {}  # output samples dictionary
        duplicates_out = {}  # output duplicates dictionary
        for sample in samples_in:
            # validate each sample
            files, errors, duplicates = \
                validator.validate_names_and_find_duplicates(
                    samples_in[sample])
            # add duplicates,errors and samples
            if len(duplicates) > 0:
                duplicates_out[sample] = duplicates  # add duplicates
            elif len(errors) > 0:
                errors_out[sample] = errors  # add errors
            else:
                samples_out[sample] = files  # append expanded path to output
        return samples_out, errors_out, duplicates_out  # return output

#-------------------------------------------------------------------------
    def get_sample_names(self, regroup_file_data):
        '''
        Given regroup file data as list,
        Returns dictionary of all sample names and corresponding files
        and errors or duplicates, if any.
        '''
        samples = {}  # output samples dictionary
        duplicate = ''
        for i in range(0, len(regroup_file_data)):
            # remove new line character
            line = regroup_file_data[i].strip("\n")
            if len(line) > 0:  # continue only if line has content
                # search current line for occurence of whitespace characters
                matchobj = re.match(r'\s*(\S*)(\s*).*', regroup_file_data[i],
                                    re.M | re.I)
                delimiter = matchobj.groups()[1]  # find delimiter character
                if len(delimiter) == 0:
                    return {}, ''  # return if invalid line found
                sample_name, sample_files = regroup_file_data[
                    i].split(delimiter, 1)
                sample_files = sample_files.strip(" ").strip("\n").split(",")
                sample_files = map(string.strip, sample_files)  # strip spaces
                # add to output only if unique
                if sample_name not in samples:
                    samples[sample_name] = sample_files
                else:
                    return samples, sample_name  # return if duplicate found
        return samples, duplicate  # return samples and duplicates data
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
