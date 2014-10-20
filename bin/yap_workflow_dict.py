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

class workflow_dictionary:

    def __init__(self):
        pass

    def create_default_wf_dict(self):
        """
        Return the default workflow dictionary

        """
        default_wf_dict = {
            "comment": "",
            "analyst_name": "",
            "organisation_name": "",
            "instrument_type": "",
            "speciment_info": "",
            "seq_type": "",
            "workflow_type": "",
            "workflow_name": "",
            "input_files_path": "",
            "input_files": "",
            "input_file_format": "",
            "regroup_output": "no",
            "regroup_file": "",
            "paired_end_data": "yes",
            "max_read_length": "150",
            "file_chunk_size": "1024",
            "data_distribution_method": "chunk_based",
            "output_files_path": "",
            "preprocess_output_file_format": "fastq",
            "run_preprocess_analysis": "no",
            "preprocess_configuration_file": "",
            "write_preprocessed_data": "no",
            "run_reference_alignment": "no",
            "aligner_configuration_file": "",
            "alignment_sort_order": "unsorted",
            "merge_alignment_output": "yes",
            "run_postprocess_analysis": "no",
            "postprocess_configuration_file": ""
	    }
        return default_wf_dict

    def validate_wf_dict(self, wf_dict):
        """
        Validate the given workflow dictionary.
        Return the error list.

        """
        error_list = []
        value_non_whitespace_keys = [
            "workflow_name",
            "input_files_path",
            "input_files",
            "input_file_format",
            "output_files_path"]
        value_yes_or_no_keys = [
            "regroup_output",
            "paired_end_data",
            "run_preprocess_analysis",
            "write_preprocessed_data",
            "run_reference_alignment",
            "merge_alignment_output",
            "run_postprocess_analysis"]

        value_is_int_keys = ["max_read_length", "file_chunk_size"]

        dependency_keys = [
            [
                "regroup_output", "regroup_file"], [
                "run_preprocess_analysis", "preprocess_configuration_file"], [
                "run_reference_alignment", "aligner_configuration_file"], [
                    "run_postprocess_analysis", "postprocess_configuration_file"]]

        for key in value_non_whitespace_keys:
            # Check if wf_dict[key] is a non_whitespace string or empty string
            if not wf_dict[key]:
                error_message = "Error: " + key + \
                    " = ''. Please give a valid value for " + key
                error_list.append(error_message)
            if wf_dict[key].isspace():
                error_message = "Error: The value of " + key + \
                    " should be non_whitespace string. Please check."
                error_list.append(error_message)

        for key in value_yes_or_no_keys:
            # Check if wf_dict[key] is either "yes" or "no"
            if wf_dict[key] not in ["no", "yes"]:
                error_message = "Error: The value of " + key + \
                    " should be 'yes' or 'no'. Please check."
                error_list.append(error_message)

        for key in value_is_int_keys:
            # Check if wf_dict[key] is positive int
            if not self.str_is_positive_int(wf_dict[key]):
                error_message = "Error: The value of " + key + \
                    " should be positive int number, such as '100'. Please check."
                error_list.append(error_message)

        for pair in dependency_keys:
            # If pair[0] is 'yes', Check if wf_dict[pair[1]] is
            # a non_whitespace string or empty string
            if wf_dict[pair[0]] == "yes":
                key = pair[1]
                if not wf_dict[key]:
                    error_message = "Error: " + key + \
                        " = ''. Please give a valid value for " + key
                    error_list.append(error_message)
                if wf_dict[key].isspace():
                    error_message = "Error: The value of " + key + \
                        " should be non_whitespace string. Please check."
                    error_list.append(error_message)

        if wf_dict["preprocess_output_file_format"] not in ["fastq", "fasta", "tab"]:
            error_message = "Error: The value of preprocess_output_file_format should be one of 'fastq', 'fasta' and 'tab'. Please check."
            error_list.append(error_message)
        if wf_dict["data_distribution_method"] not in ["chunk_based", "file_based"]:
            error_message = "Error: The value of data_distribution_method should be either 'chunk_based' or 'file_based'. Please check."
            error_list.append(error_message)
        if wf_dict["alignment_sort_order"] not in ["both", "coordinate", "queryname", "unsorted"]:
            error_message = "Error: The value of alignment_sort_order should be one of 'both', 'coordinate', 'queryname' and 'unsorted'. Please check."
            error_list.append(error_message)
        return error_list

    def str_is_positive_int(self, str):
        """
        Check if given string is int. e.x. "100" is positive int,
         "100abc", "-100", "100.00" are not.

        """
        try:
            a = int(str)
        except ValueError:
            return False
        else:
            if a > 0:
                return True
            else:
                return False

    def make_global(self,workflow_config_dict):
        """Unpacks configuration dictionary and make the variables available 
        globally"""
        #workflow general variables 
        global yap_version
        yap_version= workflow_config_dict["yap_version"] 
        global yap_src
	yap_src = workflow_config_dict["yap_src"]
        global username
        username= workflow_config_dict["username"] 
        global workflow_name
	workflow_name = workflow_config_dict["workflow_name"]
        global nprocs
	nprocs = int(workflow_config_dict["nprocs"])
        global job_id
        job_id = workflow_config_dict["job_id"]
	'''if no job id available, assign a random number, appended to 
        temporary files incombination with processors number to maintain unique name'''
        global random_id
        random_id = workflow_config_dict["random_id"]
        global data_distribution_method
	data_distribution_method = workflow_config_dict["data_distribution_method"]
        global file_chunk_size
	file_chunk_size = workflow_config_dict["file_chunk_size"]
        global regroup_output
	regroup_output = workflow_config_dict["regroup_output"]
	#configuration sanity check results 
        global check_results
        check_results = workflow_config_dict["check_results"]
        #output, log paths and temp dirs
        global yap_temp_source
        yap_temp_source = workflow_config_dict["yap_temp_source"]
        global yap_temp_user_dir
        yap_temp_user_dir = workflow_config_dict["yap_temp_user_dir"]
        global log_path
        log_path = workflow_config_dict["log_path"]
        global err_log_path
        err_log_path = workflow_config_dict["err_log_path"] 
        global stat_log_path
        stat_log_path = workflow_config_dict["stat_log_path"]
        global usage_log
        usage_log = workflow_config_dict["usage_log"]
        global temp_dir_path
        temp_dir_path = workflow_config_dict["temp_dir_path"]
        global workflow_output_path
        workflow_output_path = workflow_config_dict["workflow_output_path"]
        global regroup_output_path
        regroup_output_path = workflow_config_dict["regroup_output_path"]
        global sample_track_log
        sample_track_log = workflow_config_dict["sample_track_log"]
        global consolidated_output
	consolidated_output = workflow_config_dict["consolidated_output"]
        #input  
        global inp_files_list
        inp_files_list = workflow_config_dict["inp_files_list"]
        global input_files_path
        input_files_path = workflow_config_dict["input_files_path"]
        global file_basecount_dict
        global file_basecount_dict
        file_basecount_dict = workflow_config_dict["file_basecount_dict"]
        global regroup_arr
        regroup_arr = workflow_config_dict["regroup_samples"]
        global paired_files_split_arr
        paired_files_split_arr = workflow_config_dict["paired_files_split_arr"]
        global input_files
	input_files = workflow_config_dict["input_files"]
        global input_file_format
	input_file_format = workflow_config_dict["input_file_format"] 
        global max_read_length
	max_read_length = workflow_config_dict["max_read_length"]
        global paired_end_data
	paired_end_data = workflow_config_dict["paired_end_data"]
        global format_specific_lines
	format_specific_lines = workflow_config_dict["format_specific_lines"]
        #stages 
        global run_preprocess_analysis
	run_preprocess_analysis=workflow_config_dict["run_preprocess_analysis"]
        global run_reference_alignment
	run_reference_alignment = workflow_config_dict["run_reference_alignment"]
        global run_postprocess_analysis
	run_postprocess_analysis = workflow_config_dict["run_postprocess_analysis"]
	#commands 
        global preprocess_cmd_arr
        preprocess_cmd_arr = workflow_config_dict["preprocess_cmd_arr"]
        global aligner_cmd_arr
        aligner_cmd_arr = workflow_config_dict["aligner_cmd_arr"] 
        global postprocess_cmd_arr
        postprocess_cmd_arr = workflow_config_dict["postprocess_cmd_arr"]
        #stage specific variables
	#preprocess
        global barcode_dict
        barcode_dict = workflow_config_dict["barcode_dict"]
        global contaminant_dict
        contaminant_dict = workflow_config_dict["contaminant_dict"]
        global preprocess_output_file_format
	preprocess_output_file_format = workflow_config_dict["preprocess_output_file_format"]
        global write_preprocessed_data
	write_preprocessed_data = workflow_config_dict["write_preprocessed_data"]
	#alignment
	global aligner_output_key_arr
        aligner_output_key_arr = workflow_config_dict["aligner_output_key_arr"]
	global alignment_sort_order
	alignment_sort_order = workflow_config_dict["alignment_sort_order"]
	global merge_alignment_output
	merge_alignment_output = workflow_config_dict["merge_alignment_output"]
	#postprocess 
	global list_of_samples
	list_of_samples = workflow_config_dict["list_of_samples"]
	global list_of_samples_to_compare
	list_of_samples_to_compare = workflow_config_dict["list_of_samples_to_compare"]
