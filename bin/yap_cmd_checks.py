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
import distutils.spawn
from yap_path_checks import *

class yap_cmd_checks:

    """
    Provides methods to validate yap commands and command arrays

    """
    # global variable
    # store command execution flag
    global enable
    global exceptions_list

    # Class Initialization
    def __init__(self):
        self.enable = "yes"
        self.error_list = []
        self.warning_list = []
        self.exceptions_list = ["calculate_basecount_metrics",
                                "cd","eqp_rename_reads"]

    def is_valid_cmd_name(self, cmd):
        """
        Given a yap command,
        if enable is set to yes, validates the command
        Returns true if command is valid, a warning message if command
        argument contains a basename, false otherwise.
        If enable is set to no, returns true
        A command is valid if an executable corresponding to the command
        name exists in system and all arguments are valid

        """
        if self.enable != "yes":
            return True  # True is default when disabled
        if len(cmd) <= 0:
            return False  # return false if empty input
        # split name and argumens, name is the command executable
        name, args = self.split_cmd_name_args(cmd)
        # check if command executable exists and is set as executable
        if self.check_exe(name) == False:
            return False
        # for valid command name, return the result of cmd arg validaton
        return self.is_valid_cmd_args(args)

#-------------------------------------------------------------------------
    def split_cmd_name_args(self, cmd):
        """
        Given a yap command,
        Splits command name and arguments into two strings
        Returns the strings back to caller.

        """
        if len(cmd) <= 0:
            return "", ""  # return empty strings if input was empty
        # get index of first space character
        spc_idx = cmd.find(" ")
        # no command arguments if there is no space in command
        if spc_idx < 0:
            return cmd, ""
        else:
            return cmd[0:spc_idx], cmd[spc_idx + 1:]

#-------------------------------------------------------------------------
    def is_valid_cmd_args(self, cmd_arg_string):
        """
        Given a yap command argument string,
        if enable is set to yes, validates the command arguments
        Returns true if it is valid,a warning message if a basename is found in
        one of the arguments, false otherwise.
        If enable is set to no, returns true
        A command argument is valid if all paths or basenames contained in it
        exist in system

        """
        if self.enable != "yes":
            return True  # True is default when disabled
        # strip starting and trailing spaces
        cmd_arg_string = cmd_arg_string.strip()
        if len(cmd_arg_string) <= 0:
            return True  # return True if empty string
        # find all file names and paths in argument string
        path_ops = yap_path_checks()
        cmd_arg_items = path_ops.extract_files_and_paths(cmd_arg_string)
        # return True if no files/paths in command arguments
        if len(cmd_arg_items) <= 0:
            return True
        # check if all files and paths exist in command argument
        return self.is_valid_cmd_paths(
            cmd_arg_items,
            path_ops.valid_path_or_link,
            path_ops.basename_exists)
#-------------------------------------------------------------------------

    def is_valid_cmd_paths(self, items, path_checker, basename_checker):
        """
        Given a list of files and paths and a checking function,
        Returns
                1. true if all items are valid as per path checking funciton,
                2. warning message if some items are valid as per only basename
                checking function
                3. false otherwise.

        Invariant 1: If an item contains output_file then it is always valid,
        regardless of checker output

        """
        warning = ""  # output warning message
        for item in items:
            # contiune only if a path containing output_file does not exist
            if (path_checker(item) == False and item.find("output_file") < 0):
                # check for basename, return false if basename doesn't exist
                if basename_checker(item) == False:
                    return False
                else:
                    warning = warning + "basename in " + item + ":"
        if len(warning) > 0:
            return warning  # return warning if there
        else:
            return True  # return True if all items are valid

#-------------------------------------------------------------------------

    def check_sorted(self, postprocess_cmd_arr, alignment_sort_order):
        """
        Check if input_file_type of all commands in postprocess match the
        alignment_sort_order which is set in workflow config file.
        Give Error infomation if not match.

        """

        for post_cmd in postprocess_cmd_arr:
            # each command in postprocess
            file_info = post_cmd[1]  # file_info e.x.:
            # ['input_file_type *queryname*.sam',
            #  'input_directory aligner_output']
            command_name = post_cmd[2][0][0]
            if len(file_info) >= 2:
                self.check_sorted_helper(file_info,
                                         alignment_sort_order,
                                         command_name)
            else:
                self.error_list.append("Error: you didn't provide" +
                                       " input_file_type or input_directory.")

#-------------------------------------------------------------------------
    def check_sorted_helper(
            self,
            file_info,
            alignment_sort_order,
            command_name):
        """
        file_info includes 'input_file_type'.
        Check if 'input_file_type' matches alignment_sort_order when it
        is 'queryname' or 'coordinate'.

        """
        for info in file_info:
            if info.find("input_file_type") != -1:
                if 'queryname' in info:
                    self.check_match(
                        'queryname', alignment_sort_order, command_name)
                elif 'coordinate' in info:
                    self.check_match(
                        'coordinate', alignment_sort_order, command_name)

#-------------------------------------------------------------------------
    def check_files_when_postproc_only(
            self,
            postprocess_cmd_arr,
            input_files_path):
        """
        When run postprocess alone, check if needed file of given
        input_file_type exists in given input_directory.

        """
        path_check = yap_path_checks()  # new object of yap_path_checks
        for post_cmd in postprocess_cmd_arr:
            # each command in postprocess
            file_info = post_cmd[1]
            command_name = post_cmd[2][0][0]
            found = 'not check'
            # file_info e.g.:
            # ['input_file_type *queryname*.sam',
            #  'input_directory aligner_output']
            for info in file_info:
                # get input_file_type and input_directory
                if info.startswith("input_file_type"):
                    l = len('input_file_type')
                    input_file_type = info[l:].strip()
                if info.startswith("input_directory"):
                    l = len("input_directory")
                    input_directory = info[l:].strip()
            if input_directory == 'aligner_output':
                # when the input_directory is 'aligner_output',
                # check if the required input_file_type
                # exists in the input_directory
                input_directory = input_files_path.split(';')
                found = \
                    path_check.file_or_basename_exist_in_any(
                        input_directory,
                        input_file_type)
            if found == 'files_found' or found == 'not check':
                pass
            elif found == 'basename found':
                warning_string = "Warning: basename found " + \
                    "instead of filename in " + \
                    " input file type. Please" + \
                    " make sure that command '" + \
                    command_name + "' can work with" + \
                    "basenames."
                self.warning_list.append(warning_string)
            elif found == 'not found':
                self.error_list.append("Error: required file " +
                                       input_file_type + " not exists" +
                                       " in aligner_output directory.")
#-------------------------------------------------------------------------

    def check_match(self, input_file_type, alignment_sort_order, command_name):
        """
        Check if given input_file_type matches alignment_sort_order


        """

        if alignment_sort_order == 'both':
            return 1  # match
        elif input_file_type != alignment_sort_order:
            # not match
            error_string_tail = " but your alignment_sort_order in " + \
                "workflow config file is set to " + \
                alignment_sort_order + '.'
            self.error_list.append("Error: in postprocess, command '"
                                   + command_name + "' requires " +
                                   input_file_type +
                                   " sort," + error_string_tail)
        return 0

#-------------------------------------------------------------------------

    def check_exe(self, command):
        """
        Given a command string,
        Retuns true if command is executable, false otherwise
        """
        # expand environment variables in command
        command = os.path.expandvars(command)
        if command in self.exceptions_list:
            return True  # true for all exceptions
        # find path of the given command exectuable
        find_exe_result = distutils.spawn.find_executable(command)
        if find_exe_result is None:
            return False  # not an executable
        # check if executable path exists (applies only when full path is
        # given)
        path_checker = yap_path_checks()
        if path_checker.valid_path_or_link(find_exe_result):
            # return true if execute is set on executable file
            return os.access(find_exe_result, os.X_OK)
        else:
            return False  # return false if executable does not exist

#-------------------------------------------------------------------------
