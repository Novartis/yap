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
import string
import glob

class yap_path_checks:

    """
    Provides methods to check if a string contains a path and if a given
    path exists or not

    """

#-------------------------------------------------------------------------

    def extract_files_and_paths(self, input, count=0):
        """
        Given input string and a count
        Returns, if found,  all(or upto the given count) unix like path or file
        patterns present in the given string,
        else an empty list if no file or folder is found
        Invariant 1: count=0 means return all unix like path patterns
        Invariant 2: File names and paths do not contain any spaces

        """
        # output variable
        result = []
        # perform checks only for non-empty paths
        if len(input) > 0:

            # split input with space as delimiter
            input = input.strip()
            input_data = input.split(" ")

            # count tracker
            counter = 0

            # loop through each word in input data
            for word in input_data:
                # perform check only if non empty word
                word = word.strip()
                if len(word) > 0:
                    # add to output if path pattern is found
                    if self.is_path_pattern(word):
                        result.append(word)
                        counter = counter + 1
                    ''' Excluding file pattern checks as users will be requested
                     to use path like patterns for files also. '''

                    # else add to output if file pattern is found
                    # elif self.is_file_pattern(word) == True:
                    #    result.append(word)
                    #    counter = counter + 1

                # return result if requested number of patterns have been found
                if (count > 0 and count == counter):
                    return result

        # end if len(input) > 0

        # return output
        return result
#-------------------------------------------------------------------------

    def is_path_pattern(self, input):
        """
        Given an input string ,
        Returns true if string contains a unix like path pattern,
        false otherwise
        """

        # Unix path can be of one of below type:
        # 1. ../../
        # 2. /<some string>
        # 3. ./

        # perform checks only if input is non empty
        if len(input) > 0:

            # Doing an easy check currently , assuming any string containing /
            # character is a path
            if input.find("/") >= 0:
                return True
            else:
                return False

            # Below code is currently being ignored
            # set loop counter to max 3 characters
            loop_counter = min(3, len(input))
            # Loop only upto first 3 characters of given string
            for i in range(0, loop_counter):
                # continue next iteration if character is .
                if input[i] == ".":
                    pass
                # return true if a forward slash is found
                elif input[i] == "/":
                    return True

            # end for
        # return false if loop did not return true
        return False
#-------------------------------------------------------------------------

    def is_file_pattern(self, input):
        """
        Given an input string ,
        Returns true if string contains a file name like pattern,
        false otherwise
        Invariant : any input with a . in it is considered a file name
        """
        # return true if a . is found in input
        if string.find(input, ".") < 0:
            return False
        else:
            return True
#-------------------------------------------------------------------------

    def find_paths_with_forward_slash(self, input, count=0):
        """
        Given input string and a count
        Returns, if found,  all(or upto the given count) unix like path patterns
        present in the given string, else an empty list if no path is found
        Invariant 1: count=0 means return all unix like path patterns
        Invariant 2: Paths do not contain any spaces

        """
        # output variable
        result = []

        # loop variable
        counter = 0

        while len(input) > 0:
            # return empty string if input does not contain a "/" character
            path_start_index = string.find(input, "/")

            if path_start_index < 0:
                return result

            # get the index of first space after the "/" character
            space_index = string.find(input, " ", path_start_index)

            # if no space is found, path must end till end of string
            if space_index < 0:
                result.append(input[path_start_index:])
                input = ""
            else:
                result.append(input[path_start_index:space_index])
                input = input[space_index + 1:]

            # increment counter
            counter = counter + 1

            # return result if requested number of paths have been found
            if (count > 0 and count == counter):
                return result

        # end while

        # return output
        return result
#-------------------------------------------------------------------------

    def path_exists(self, input):
        """
        Given a unix like path pattern,
        Returns true if path exists on current system, false otherwise

        """

        # expand environment variables in path
        input = os.path.expandvars(input)
        # get absolute path
        abs_path = os.path.abspath(input)

        return os.path.exists(abs_path)
#-------------------------------------------------------------------------

    def basename_exists(self, path, basename=""):
        """
        Given a path and/or basename,
        Returns True if any files exist in given path with the given basename,
        False otherwise

        """

        # get absolute path
        abs_path = os.path.abspath(path)
        # add wildcard character * to basename
        if len(basename) > 0:
            basename = "/*" + basename + "*"
        else:
            abs_path = abs_path + "*"
        # return the result if basename is found
        if glob.glob(abs_path + basename):
            return True
        else:
            return False
#-------------------------------------------------------------------------

    def files_exist(self, path, file_search_string):
        """
        Given a unix path and search string for files,
        Returns true if path contains any files based on the given pattern,
        false otherwise

        """

        # return false if invalid path is provided
        if self.path_exists(path) == False:
            return False
        # get absolute path
        abs_path = os.path.abspath(path)
        # In case of valid path, return result of file search
        if glob.glob(abs_path + "/" + file_search_string):
            return True
        else:
            return False
#-------------------------------------------------------------------------

    def file_or_basename_exist_in_any(self, paths, file_search_string):
        """
        Given a list of unix paths and search string for files or
        file basename,
        Returns
             1. 'files found.' if any of the given paths contain files based
             on the given pattern,
             2. 'basename found.' if any of the given paths contain files
             based on basename in the given pattern
             3. 'not found' if no files are found in above scenarios

        """

        # check for just file names first
        if self.files_exist_in_any(paths, file_search_string):
            return 'files found'
        # check for base names if file names were not found
        for path in paths:
            if self.basename_exists(path, file_search_string):
                return 'basename found'
        # return none if neither was found
        return 'not found'

#-------------------------------------------------------------------------
    def files_exist_in_any(self, paths, file_search_string):
        """
        Given a list of unix paths and search string for files,
        Returns true if any of the given paths contain file based on the given
        pattern, false otherwise

        """

        # check if files are found in any paths
        for path in paths:
            if self.files_exist(path, file_search_string):
                return True

        # return false if files were not found in any paths
        return False

#-------------------------------------------------------------------------
    def files_exist_in_all(self, paths, file_search_string):
        """
        Given a list of unix paths and search string for files,
        Returns true only if all of the given paths contain file based on the
        given pattern, false otherwise

        """

        # check if files are found in all paths
        for path in paths:
            if self.files_exist(path, file_search_string) == False:
                return False

        # return true if files were found in all paths
        return True

#-------------------------------------------------------------------------
    def invalid_path_or_link(self, path):
        '''
        Given a path,
        Returns false if it is an aboslute path or
        if it is a valid link to a symbolic path,
        true otherwise
        '''
        import re
        yap_meta_words = [
            "workflow_directory",
            "output_directory",
            "output_file",
            "file_based_input",
            "input_files_path",
            "directory_based_input"]
        path = self.translate_path(path)  # translate env variables
        path_new = path  # copy path
        for i in yap_meta_words:
            if re.search(i, path) is not None:
                return False
        if os.path.islink(path):
            path_new = os.readlink(path)  # check links
        if self.path_exists(path_new):
            if os.access(path_new, os.R_OK) == False:
                return True  # invalid if no access
            else:
                return False  # false if valid path and access
        else:
            if self.is_path_pattern(path_new) == False:
                # check if actual path exists in current directory
                return self.invalid_path_or_link(
                    os.path.dirname(path) +
                    "/" +
                    path_new)
            else:
                return True  # broken symbolic link

#-------------------------------------------------------------------------
    def valid_path_or_link(self, path):
        '''
        Given a path,
        Returns true if path is not found to be invalid, false otherwise
        '''
        return (not self.invalid_path_or_link(path))

#-------------------------------------------------------------------------
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

#------------------------------------------------------------------------------
