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
from mpi4py import MPI
import time
import glob
import sys
import random
import re
import yap_print_info 
import yap_exe
import yap_tools 
import yap_file_io 
import yap_preprocess
import yap_init 
import yap_postprocess 
import yap_check
import yap_prov
import yap_log
import yap_aligner
import yap_workflow_dict as wd
'''
This script intializes MPI communicator, configuration paramteres,
and controls parallel data flow through different analysis stages.
'''
#intializng mpi communicator
comm = MPI.COMM_WORLD
nprocs = comm.Get_size()
myrank = comm.Get_rank()
#intializing variables 
run_mode = ''
print_usage = 'false'
exit_status = ''
fastqc_job_status = ''
fastq_screen_job_status = ''
preqc_prov=[]
preproc_alignment_prov=[]
workflow_prov=[]
filebased_file_split_arr=[]
ignore_files_list=[]
command_count=0
workflow_struct=[]
inp_files_list=[]
file_split_strut=[]
postprocess_cmd_arr=[]
postprocess_compare_file_cmd_arr=[]
workflow_prov_merge=[]
workflow_prov_post=[]
regroup_arr=[]
regroup_split_index=[]
filebased_split_arr_index=[]
parser_errorlist=[]
file_split_struct=[]
#root processor checks the command options
if (myrank == 0):
    if len(sys.argv) == 2 or len(sys.argv) == 3:
        if len(sys.argv) == 2:
            if sys.argv[1] == '-h' or sys.argv[1] == '--help' or sys.argv[1] == "--check" or sys.argv[1] == "-check":
                print_usage = 'True'
            else:
                workflow_config_file = sys.argv[1]
        if len(sys.argv) == 3:
            if sys.argv[1] == "--check" or sys.argv[1] == "-check":
                workflow_config_file = sys.argv[2]
                run_mode = "--check"
            else:
                print_usage = 'True'
    else:
        print_usage = 'True'
    if print_usage == 'True':
        print "Running :", wd.yap_version
        print "Options : "
        print " To print help : ", " yap --help or yap -h"
        print " To do configuration file format check : ", " yap --check [workflow_configuration_filename], ", "\"eg: yap --check workflow_configuration.cfg \""
        print " To run YAP: ", " yap -n [number of processors] [workflow_configuration_filename], ", "\"eg: yap -n 2 workflow_configuration.cfg \"", "(yap run with 2 processors)"
        print "exiting the program"
        exit_status = 'True'
    if exit_status != 'True':
	#reads the main workflow configuration file
        workflow_file_data = yap_file_io.read_file(workflow_config_file)
        if len(workflow_file_data) > 0:
	    #pass the file data to workflow parser
            workflow_struct, workflow_errorlist = yap_tools.workflow_parser(workflow_file_data, workflow_config_file,nprocs)
            if len(workflow_struct) < 2 or len(workflow_errorlist) > 0:
                print"Format Error : while parsing the workflow configuration file= ", workflow_config_file, "\n"
                print"YAP analysis general metadata mising or syntax error.", "\n"
                print"Note:Use symbol(:begin) and (:end) to define command sections,enclose variable and corresponding values in double quotes(eg \"..\")"
                print"Use (:=) in variable value assignment (eg \"variable\" := \"value\")"
                print"To add comments,start the line with symbol(#)"
                for i in workflow_errorlist:
                    print i
                exit_status = 'True'
    if exit_status != 'True':
        if len(workflow_struct) != 0:
	    #check if the workflow configuration passes the sanity checks
            parser_errorlist = yap_check.workflow_validator(workflow_struct, workflow_config_file, run_mode)
        else:
            exit_status = 'True'
        if run_mode == "--check" or len(parser_errorlist) > 0:
	    #decide the exit status based on run mode and santiy check results
            exit_status = 'True'
#communicate exit status to every processor
exit_status = comm.bcast(exit_status, root=0)
if exit_status == 'True': #all the processors to execute this
    exit()
#communicate all the workflows to all processors 
workflow_struct = comm.bcast(workflow_struct, root=0)
#loop over each workflow for data processing
for wk in range(1, len(workflow_struct)):
    #store the current workflow details into a dictionary  
    workflow_config_dict = workflow_struct[wk]
    #make workflow configuration variables global 
    workflow_obj=wd.workflow_dictionary()
    workflow_obj.make_global(workflow_config_dict)
    #summary file to store all the workflow provenance details
    f_summary_file = wd.workflow_output_path + "/" + wd.workflow_name + "_workflow_summary.txt"
    basecount_metrics_flag = ''
    file_basecount_dict=wd.file_basecount_dict
    #create temp directory, this is local for everynode
    yap_tools.create_dir(wd.yap_temp_user_dir)
    if (myrank == 0):
	#create output directory structure 
        yap_init.initialize_dir_struct()
    	# printing analysis summary
        yap_print_info.print_info()
        str_out= "-"*20 +" PROVENANCE "+ "-"*20 +"\n\n"
        yap_file_io.write_data(str_out,f_summary_file)
    comm.barrier()
    if wd.run_preprocess_analysis == "yes":
	#if preprocess is set to 'yes', perform initial qc commands
        for i in range(0, len(wd.preprocess_cmd_arr)):
            preprocess_cmd_name = wd.preprocess_cmd_arr[i][2][0][0]
            preprocess_cmd = wd.preprocess_cmd_arr[i][2][0][1]
            if re.search('calculate_basecount_metrics', preprocess_cmd_name) is not None:
                basecount_metrics_flag = 'True'
            if re.search('fastqc', preprocess_cmd_name) is not None:
                fastqc_split_index_arr = wd.paired_files_split_arr[myrank]
                fastqc_inp_files = []
                for k in range(0, len(fastqc_split_index_arr)):
                    if fastqc_split_index_arr[k] == 'no file':
                        fastqc_job_status = "Done"
                        pass
                    else:
                        fastqc_inp_files = wd.inp_files_list[
                            fastqc_split_index_arr[k]]
                        if len(fastqc_inp_files) != 0:
                            try:
                                preqc_prov += yap_preprocess.run_fastqc(fastqc_inp_files, preprocess_cmd)
                                fastqc_job_status = "Done"
                            except Exception as err:
                                print " Error : Fastqc failed"
                                print err
                print "-"*20,"FastQC analysis End Time", time.strftime("%Y/%m/%d %H:%M:%S", time.localtime()), "-"*20
		comm.barrier()
            if re.search('fastq_screen', preprocess_cmd_name) is not None:
                fastq_screen_split_arr_index = wd.paired_files_split_arr[myrank]
                fastq_screen_inp_files = []
                for k in range(0, len(fastq_screen_split_arr_index)):
                    if fastq_screen_split_arr_index[k] == 'no file':
                        fastq_screen_job_status == "Done"
                        pass
                    else:
                        fastq_screen_inp_files = wd.inp_files_list[
                            fastq_screen_split_arr_index[k]]

                        if len(fastq_screen_inp_files) != 0:
                            try:
                                preqc_prov += yap_preprocess.run_fastq_screen(fastq_screen_inp_files, preprocess_cmd)
                                fastq_screen_job_status = "Done"
                            except Exception as err:
                                print " Error : Fastq screening failed"
                                print err
                print "-"*20,"FastqScreen analysis End Time", time.strftime("%Y/%m/%d %H:%M:%S", time.localtime()), "-"*20
		comm.barrier()
        comm.barrier()
        if myrank == 0:
	    #write log information
            if len(glob.glob(wd.err_log_path + "/*fastqc_err.log")):
                for i in glob.glob(wd.err_log_path + "/*fastqc_err.log"):
                    os.system("cat" + " " + i + " >> " + i.replace("fastqc_err.log","err.log"))
                    os.system("rm" + " " + i)
            if len(glob.glob(wd.stat_log_path + "/*fastqc_stat.log")):
                for i in glob.glob(wd.stat_log_path + "/*fastqc_stat.log"):
                    os.system("cat" + " " + i + " >> " + i.replace("fastqc_stat.log","stat.log"))
                    os.system("rm" + " " + i)

            if len(glob.glob(wd.err_log_path + "/*_fastqscreen_err.log")):
                for i in glob.glob(wd.err_log_path + "/*fastqscreen_err.log"):
                    os.system("cat " + i + " >> " + i.replace("fastqscreen_err.log","err.log"))
                    os.system("rm " + i)
            if len(glob.glob(wd.stat_log_path + "/*fastqscreen_stat.log")):
                for i in glob.glob(wd.stat_log_path + "/*_fastqscreen_stat.log"):
                    os.system("cat " + i + " >> " + i.replace("fastqscreen_stat.log","stat.log"))
                    os.system("rm " + i)
            #write preprocess heading to provenance file
            yap_file_io.write_data("PREPROCESS:" + "\n\n" ,f_summary_file)
        #collect preprocess commands from all the processors and write to provenance file
        for k in range(0, wd.nprocs):
            if myrank == k:
                for gg in (preqc_prov):
            	    yap_file_io.write_data("\t" + gg + "\n\n",f_summary_file)
	    comm.barrier()
    #get file chunk split information for files assigned to each processor
    if wd.run_preprocess_analysis == "yes" or wd.run_reference_alignment == "yes":
        paired_files_split_index_arr = wd.paired_files_split_arr[myrank]
	file_split_struct=yap_tools.get_filesplit_info(paired_files_split_index_arr,wd.inp_files_list,wd.file_chunk_size,wd.nprocs,wd.format_specific_lines)
    comm.barrier()
    #combine the file chunk split information from all processors 
    for k in range(1, wd.nprocs):
        if myrank == k:
            comm.send(file_split_struct, dest=0, tag=myrank + 1001)
        if myrank == 0:
            file_split_struct_temp = comm.recv(source=k, tag=k + 1001)
            file_split_struct = file_split_struct + file_split_struct_temp
        comm.barrier()
    if myrank == 0:
        print "-"*20,"File Decompression End Time ", time.strftime("%Y/%m/%d %H:%M:%S", time.localtime()), "-"*20
        if wd.data_distribution_method == "file_based":
            #if file based workflow; then assign file names to processors, instead of chunks
            filebased_file_split_arr = yap_tools.split_files_each_proc(file_split_struct, wd.nprocs)
    #broadcast the files assigned to all processors 
    file_split_struct = comm.bcast(file_split_struct, root=0)
    filebased_file_split_arr = comm.bcast(filebased_file_split_arr, root=0)
    #set dictionary variable for eqp
    eqp_dict = {'eqp_read_counter':[]} #applies to in-house use, neglect otherwise
    if wd.run_preprocess_analysis == "yes" or wd.run_reference_alignment == "yes":
        if wd.data_distribution_method != "file_based":
	    '''root reads chunks based  precalculated file offsets, 
            sends information to assigned processor '''
            for ii in range(0, len(file_split_struct)):
                fh1 = ''
                fh2 = ''
                input_filename_local = wd.inp_files_list[ii][0]
                input_filename_local_2 = wd.inp_files_list[ii][1]
                file_name = wd.inp_files_list[ii][2]
                fh1 = yap_file_io.create_openfile_handler(input_filename_local)
                if wd.paired_end_data == 'yes':
                    fh2 = yap_file_io.create_openfile_handler(input_filename_local_2)
                file_size = file_split_struct[ii][0]
                chunk_size = file_split_struct[ii][1]
                nchunks = file_split_struct[ii][2]
                format_specific_lines = file_split_struct[ii][3]
                barcode_output_dict = {}
                postprocess_file_spilt_index = []
                end_pos = 0
                cn = 0
                for jj in range(0, (nchunks / wd.nprocs)):
                    if myrank == 0:
                        chunk_number = cn
                        end_pos, inp1, inp2 = yap_tools.read_file_chunks(
                            fh1, fh2, cn, nchunks, chunk_size, file_size, format_specific_lines)
                        cn = cn + 1
                    if wd.nprocs > 1:
                        for p in range(1, wd.nprocs):
                            if(myrank == 0):
                                end_pos, inptemp1, inptemp2 = yap_tools.read_file_chunks(
                                    fh1, fh2, cn, nchunks, chunk_size, file_size, format_specific_lines)
                                comm.send(inptemp1, dest=p, tag=p + 10)
                                comm.send(inptemp2, dest=p, tag=p + 100)
                                comm.send(cn, dest=p, tag=p + 1000)
                                cn = cn + 1
                                inptemp1 = ''
                                inptemp2 = ''
                            if(myrank == p):
                                inp1 = ''
                                inp2 = ''
                                inp1 = comm.recv(source=0, tag=myrank + 10)
                                inp2 = comm.recv(source=0, tag=myrank + 100)
                                chunk_number = comm.recv(
                                    source=0, tag=myrank + 1000)
                            #wait for all the processors to finish work
                            comm.barrier()
                    comm.barrier()
                    if myrank == 0:
			#print status information
                        print "-"*20,"Preprocess and Alignment Start Time ", time.strftime("%Y/%m/%d %H:%M:%S", time.localtime()), " For ", input_filename_local, "-"*20 + "\n" 
                    #pass chunk information to yap_exe for preprocess and alignment execution 
                    print "Entering Preprocess and Alignment section: ", "rank=", myrank, " input file = ", input_filename_local, "chunk_number=", chunk_number, "\n\n"
                    barcode_output_dict, workflow_prov = yap_exe.execute_chunk(wd.inp_files_list[ii], inp1, inp2, chunk_number,myrank,workflow_prov,eqp_dict)
                    print "Exiting Preprocess and Alignment section: ", "rank=", myrank, " input file = ", input_filename_local, "chunk_number=", chunk_number, "\n\n"
                    comm.barrier()
                    #empty input strings to release memory
                    inp1 = ''
                    inp2 = ''
                    #merging chunk basecount matrix for each file, for a particular processor 
                    for filename_key in file_basecount_dict.iterkeys():
                        if filename_key == file_name:
                            for barcode in barcode_output_dict.iterkeys():
                                if barcode in file_basecount_dict[filename_key]:
                                    file_basecount_dict[filename_key][barcode][0] += barcode_output_dict[barcode][0]
                                    file_basecount_dict[filename_key][barcode][1] += barcode_output_dict[barcode][1]
                                else:
                                    file_basecount_dict[filename_key][barcode] = [barcode_output_dict[barcode][0],barcode_output_dict[barcode][1]]
                        comm.barrier()
                    if myrank == 0:
			#print status information
                        print "-"*20,"Preprocess and Alignment End Time ", time.strftime("%Y/%m/%d %H:%M:%S", time.localtime()), "For ", input_filename_local, "-"*20 + "\n"
                comm.barrier()

        else:
            #this section applies to file based mode
            #based on myrank of the processors, get the filename index list 
            filebased_split_arr_index = filebased_file_split_arr[myrank]
            filebased_split_inp_files = []
            filebased_split_struct = []
	    #assigns filename list for each processors based on the precalculated file index information 
            for k in range(0, len(filebased_split_arr_index)):
                if filebased_split_arr_index[k] == 'no file':
                    pass
                else:
                    filebased_split_inp_files.append(
                        wd.inp_files_list[filebased_split_arr_index[k]])
                    filebased_split_struct.append(
                        file_split_struct[filebased_split_arr_index[k]])
            #each processor to iterate over assigned list of files
            if len(filebased_split_inp_files) != 0:
                for ii in range(0, len(filebased_split_inp_files)):
                    fh1 = ''
                    fh2 = ''
                    input_filename_local = filebased_split_inp_files[ii][0]
                    input_filename_local_2 = filebased_split_inp_files[ii][1]
                    file_name = filebased_split_inp_files[ii][2]
                    err_chunk_file = wd.err_log_path + "/" + file_name + \
                        "_log_temp/" + file_name + "_" + str(ii).zfill(6)
                    stat_chunk_file = wd.stat_log_path + "/" + file_name + \
                        "_log_temp/" + file_name + "_" + str(ii).zfill(6)
                    if myrank == 0:
                        print "-"*20,"Preprocess and Alignment Start Time ",time.strftime("%Y/%m/%d %H:%M:%S", time.localtime()), " For ", input_filename_local,"-"*20
                    fh1 = yap_file_io.create_openfile_handler(input_filename_local)
                    if wd.paired_end_data == 'yes':
                        fh2 = yap_file_io.create_openfile_handler(input_filename_local_2)
                    file_size = filebased_split_struct[ii][0]
                    chunk_size = filebased_split_struct[ii][1]
                    nchunks = filebased_split_struct[ii][2]
                    format_specific_lines = filebased_split_struct[ii][3]
                    barcode_output_dict = {}
                    postprocess_file_spilt_index = []
                    end_pos = 0
                    cn = 0
                    if wd.run_preprocess_analysis == 'yes':
			'''
			For preprocessing, data is processed in chunks and merged together before alignment
                        '''
                        for jj in range(0, (nchunks)):
                            chunk_number = jj
                            end_pos, inp1, inp2 = yap_tools.read_file_chunks(
                                fh1, fh2, jj, nchunks, chunk_size, file_size, format_specific_lines)
                            print "Entering Preprocess  section: ", "rank=", myrank, " input file = ", input_filename_local, "chunk_number=", chunk_number, "\n\n"
                            str1 = "Entering Preprocess section: rank=" + \
                                str(myrank) + " input file = " + input_filename_local + \
                                "chunk_number=" + str(chunk_number) + "\n"
			    #calling yap_exe to execute preprocess section for all the chunks
                            barcode_output_dict, workflow_prov = yap_exe.execute_chunk(filebased_split_inp_files[
                                                                         ii], inp1, inp2, chunk_number,myrank, workflow_prov,eqp_dict)
                            print "Exiting Preprocess  section: ", "rank=", myrank, " input file = ", input_filename_local, "chunk_number=", chunk_number, "\n\n"
                            inp1 = ''
                            inp2 = ''
                            #merging chunk basecount matrix for each file, for a particular processor 
                            for filename_key in file_basecount_dict.iterkeys():
                                if filename_key == file_name:
                                    for barcode in barcode_output_dict.iterkeys():
                                        if barcode in file_basecount_dict[filename_key]:
                                            file_basecount_dict[filename_key][barcode][
                                                0] += barcode_output_dict[barcode][0]
                                            file_basecount_dict[filename_key][barcode][
                                                1] += barcode_output_dict[barcode][1]
                                        else:
                                            file_basecount_dict[filename_key][barcode] = [
                                                barcode_output_dict[barcode][0],
                                                barcode_output_dict[barcode][1]]
                    if wd.run_reference_alignment == 'yes':
				#calling alignment section for filebased approach
                    		workflow_prov,file_basecount_dict=yap_exe.execute_file(input_filename_local,input_filename_local_2,file_name,chunk_number,myrank,ii,file_basecount_dict,eqp_dict)
		    if myrank == 0:
                        print "-"*20 + " Preprocess and Alignment End Time " + time.strftime("%Y/%m/%d %H:%M:%S", time.localtime()) + " For " + input_filename_local, "-"*20
        comm.barrier()
    if myrank == 0:
        if wd.run_reference_alignment == "yes":
            #combining temporary logs
            yap_log.combine_temp_log(wd.err_log_path, wd.stat_log_path)
    if wd.run_preprocess_analysis == "yes":
        if basecount_metrics_flag == 'True':
            #merging chunk basecount matrix for each file, from all the processors 
            for p in range(1, wd.nprocs):
                if (p == myrank):
                    comm.send(file_basecount_dict, dest=0, tag=p + 999)
                if (myrank == 0):
                    file_basecount_dict1 = comm.recv(source=p, tag=p + 999)
                    for filename_key in file_basecount_dict.iterkeys():
                        barcode_output_dict1 = file_basecount_dict1[
                            filename_key]
                        for barcode in barcode_output_dict1.iterkeys():
                            if barcode in file_basecount_dict[filename_key]:
                                file_basecount_dict[filename_key][barcode][
                                    0] += barcode_output_dict1[barcode][0]
                                file_basecount_dict[filename_key][barcode][
                                    1] += barcode_output_dict1[barcode][1]
                            else:
                                file_basecount_dict[filename_key][barcode] = [
                                    barcode_output_dict1[barcode][0],
                                    barcode_output_dict1[barcode][1]]
                comm.barrier()
    comm.barrier()
    if myrank == 0:
        print "-"*20, "Postproces Start Time", time.strftime("%Y/%m/%d %H:%M:%S", time.localtime()), "-"*20
	#Based on error log, filter the samples which have failed
        inp_files_list,list_of_samples,list_of_samples_to_compare,ignore_files_list = yap_log.sample_filter()
	#separate commands into two categories, samples which needs to be analyzed individually or in groups of multiple samples 
        postprocess_cmd_arr, postprocess_compare_file_cmd_arr = yap_postprocess.separate_postprocess_cmds(wd.postprocess_cmd_arr)
        preproc_alignment_prov = yap_prov.preprocess_alignment_prov(file_split_struct,preproc_alignment_prov)
        #write the provenance list from preprocess and alignment section to provenance file
        for gg in (preproc_alignment_prov):
            yap_file_io.write_data("\t"+gg+"\n\n",f_summary_file)
        if wd.run_preprocess_analysis == "yes":
            if basecount_metrics_flag == 'True':
                print "Writing basecount metrics : rank= ", myrank
                # write the basecount matrix to file
                yap_preprocess.write_basecount_matrix(file_basecount_dict)
        if wd.regroup_output == "yes":
		regroup_arr=wd.regroup_arr
		#filter out error samples from regroup list, discard groups where even one of the sample has failed.
                if len(ignore_files_list) !=0:
                    regroup_arr=[i for i in wd.regroup_arr if not(set(i[1]) and set(ignore_files_list))]
                #calculate how many files needs to distrbiuted for each processor for postprocess 
                regroup_split_index = yap_tools.split_files_each_proc(regroup_arr, wd.nprocs)
    #broadcast filtered file list information to all processors 
    inp_files_list = comm.bcast(inp_files_list, root=0)
    ignore_files_list = comm.bcast(ignore_files_list, root=0)
    if wd.regroup_output == "yes":
    	workflow_output_path = wd.workflow_output_path + "/" + "regroup_output"
    else:
	workflow_output_path = wd.workflow_output_path 
    if inp_files_list != []:
        try:
            #regrouping of the output is only valid for alignment stage 
            if wd.regroup_output == "yes" :
                #broadcast regroup list to all processors 
                regroup_split_index = comm.bcast(regroup_split_index, root=0)
                regroup_arr = comm.bcast(regroup_arr, root=0)
                local_regroup_index_arr = regroup_split_index[myrank]
                for jj in range(len(local_regroup_index_arr)):
	            #if there is no file assigned; then pass
                    if local_regroup_index_arr[jj] == 'no file':
                        pass
                    else:
                        #assign regroup sample list based on precalculated index 
                        local_regroup_arr = regroup_arr[local_regroup_index_arr[jj]]
			#merge alignment output data for files based on groups from regroup list
                        workflow_prov_merge = yap_aligner.regroup_files(local_regroup_arr,workflow_prov_merge)
                    comm.barrier()
                comm.barrier()
            if myrank == 0:
                #regrouping of the output is only valid for alignment stage
                if wd.regroup_output == "yes" :
                    '''if regroup workflow; then change workflow input/output path for postprocess
                    based on regroup output structure '''
                    inp_files_list = []
                    file_base_count_dict = {}
                    for ii in range(0, len(regroup_arr)):
                        inp_files_list.append( ['', '', regroup_arr[ii][0], ''])
	                file_basecount_dict[regroup_arr[ii][0]] = wd.barcode_dict
            #calculate filename list for postprocessing for each processors 
            postprocess_file_split_index = yap_tools.split_files_each_proc(inp_files_list, wd.nprocs)
	    '''broadcast file basecount dictionary and input file list; if there are any changes
            because of regrouping, those will be passed to all processors'''
            file_basecount_dict = comm.bcast(file_basecount_dict, root=0)
            inp_files_list = comm.bcast(inp_files_list, root=0)
            #broadcast individual sample postprocess command list 
            postprocess_cmd_arr = comm.bcast(postprocess_cmd_arr, root=0)
            #broadcast command list requiring multiple samples for processing  
            postprocess_compare_file_cmd_arr = comm.bcast(postprocess_compare_file_cmd_arr, root=0)
            #broadcast the posprocess sample assignment list
            postprocess_file_split_index = comm.bcast(postprocess_file_split_index, root=0)
            #based on my rank, get the filename list assigned for postprocessig 
            file_split_index_arr = postprocess_file_split_index[myrank]
            #iterate over each filename
            for k in range(0, len(file_split_index_arr)):
                file_index = file_split_index_arr[k]
                if file_index == "no file":
                    pass
                else:
                    postprocess_file_name = inp_files_list[file_index][2]
                    temp_err_log = wd.err_log_path + "/" + postprocess_file_name + "_err.log"
                    temp_stat_log = wd.stat_log_path + "/" + postprocess_file_name + "_stat.log"
                    local_file_basecount_dict = {}
                    if wd.run_preprocess_analysis == "no" and wd.run_reference_alignment == "no" and wd.run_postprocess_analysis == "yes":
                        postprocess_file_name_tmp = inp_files_list[file_index][0]
                        local_file_basecount_dict[postprocess_file_name_tmp] = file_basecount_dict[postprocess_file_name]
                    else:
                        local_file_basecount_dict[postprocess_file_name] = file_basecount_dict[postprocess_file_name]
		    #merge all the chunk preprocessed data per file
                    if wd.run_preprocess_analysis == "yes":
                        if wd.write_preprocessed_data == "yes":
                            if wd.merge_alignment_output == "yes":
                                if wd.data_distribution_method != "file_based":
                                    print "Merging preprocessed output: ", "rank=", myrank, " input file = ", postprocess_file_name, "\n\n"
                                    yap_preprocess.merge_preprocessed_output(local_file_basecount_dict)
                                    print "Done Merging preprocessed output", "rank=", myrank, " input file = ", postprocess_file_name, "\n\n"
                    #merge all the chunk alignment data per file
                    if wd.regroup_output == "yes":
                    	wd.merge_alignment_output == "no"
                    if wd.run_reference_alignment == "yes":
                        if wd.merge_alignment_output == "yes":
                            print "Merging Alignment output: ", "rank=", myrank, " input file = ", postprocess_file_name, "\n"
                            workflow_prov_merge = yap_aligner.merge_alignment_output(local_file_basecount_dict,workflow_prov_merge,temp_err_log,temp_stat_log)
                            print "Done Merging Alignment output: ", "rank=", myrank, " input file = ", postprocess_file_name, "\n"
                    if wd.run_postprocess_analysis == "yes":
                        print "Entering Postprocess section: ", "rank=", myrank, " input file = ", postprocess_file_name, "\n"
                        if len(postprocess_cmd_arr) != 0:
                            #set postprocess log file paths 
                            post_err_log = wd.err_log_path + "/" + postprocess_file_name + "_log_temp/" + postprocess_file_name
                            post_stat_log = wd.stat_log_path + "/" + postprocess_file_name + "_log_temp/" + postprocess_file_name
			    str_out="-"*20 +"POSTPROCESS STARTED" + "\t" + str(time.strftime("%Y/%m/%d %H:%M:%S",time.localtime())) +"-"*20 +"\n"
	                    yap_file_io.write_data(str_out,post_err_log)
                            yap_file_io.write_data(str_out,post_stat_log)
                            #execute postprocess commands using yap_postprocess 
                            workflow_prov_post = yap_postprocess.run_postprocess(postprocess_cmd_arr,local_file_basecount_dict,workflow_prov_post,post_err_log,post_stat_log)
			    str_out="*" *10 +"POSTPROCESS FINISHED" + "\t" + str(time.strftime("%Y/%m/%d %H:%M:%S",time.localtime())) + "*" * 10 +"\n"
                            yap_file_io.write_data(str_out,post_err_log)
                            yap_file_io.write_data(str_out,post_stat_log)
                        print "Exiting Postprocess section: ", "rank=", myrank, " input file = ", postprocess_file_name, "\n"
                comm.barrier()
            comm.barrier()
	    #collect alignment provenance information from all processors and write to file
            if myrank == 0:
                preproc_alignment_prov = yap_prov.preprocess_alignment_prov(
                    file_split_struct,
                    preproc_alignment_prov)
                for gg in (preproc_alignment_prov):
            	    yap_file_io.write_data("\t"+gg+"\n\n",f_summary_file)
            	yap_file_io.write_data("ALIGNMENT :" + "\n\n",f_summary_file)
            #iterate over each processors, write provenance data 
            for k in range(0, wd.nprocs):
                if myrank == k:
                    for gg in (workflow_prov):
            	    	yap_file_io.write_data("\t"+gg+"\n\n",f_summary_file)
                comm.barrier()
	    #write provenance related to alignment merge commands
            if myrank == 0:
                yap_file_io.write_data("MERGE ALIGNMENT OUTPUT :" + "\n\n",f_summary_file)
                    
            for k in range(0, wd.nprocs):
                if myrank == k:
                    for gg in (workflow_prov_merge):
            	    	yap_file_io.write_data("\t"+gg+"\n\n",f_summary_file)
                comm.barrier()
	    #write provenance for postprocess commands which work on individual file (sample)
            if myrank == 0:
                command_count = 0
                yap_file_io.write_data("POSTPROCESS :" +"\n\n", f_summary_file)
            for k in range(0, wd.nprocs):
                if myrank == k:
                    for gg in (workflow_prov_post):
            	    	yap_file_io.write_data("\t"+gg+"\n\n",f_summary_file)
                comm.barrier()
            #this section is for postprocess commands requiring multiple files (samples) for processing 
            workflow_prov = []
            postprocess_compare_arr = []
            postprocess_file_compare_split_index = []
            for kk in range(0, len(postprocess_compare_file_cmd_arr)):
                if myrank == 0:
                    postprocess_compare_arr = yap_postprocess.get_postprocess_file_compare_cmd_arr(
                        [postprocess_compare_file_cmd_arr[kk]], inp_files_list)
                    postprocess_file_compare_split_index = yap_tools.split_files_each_proc(
                        postprocess_compare_arr,
                        wd.nprocs)
                postprocess_compare_arr = comm.bcast(
                    postprocess_compare_arr, root=0)
                postprocess_file_compare_split_index = comm.bcast(
                    postprocess_file_compare_split_index, root=0)
                file_compare_split_index_arr = postprocess_file_compare_split_index[
                    myrank]
                for k in range(0, len(file_compare_split_index_arr)):
                    file_index = file_compare_split_index_arr[k]
                    if file_index == "no file":
                        pass
                    else:
                        postprocess_compare_local_arr = postprocess_compare_arr[
                            file_index]
                        if len(postprocess_compare_local_arr) > 0:
                            post_err_log = postprocess_compare_local_arr[2]
                            post_stat_log = postprocess_compare_local_arr[3]
                            str_out= "*" * 10 + postprocess_compare_arr[file_index][0].upper() + " STARTED" + "\t" + str( time.strftime("%Y/%m/%d %H:%M:%S",time.localtime())) + "*" * 10 + "\n"
                            yap_file_io.write_data(str_out,post_err_log)
                            yap_file_io.write_data(str_out,post_stat_log)
                            workflow_prov = yap_postprocess.run_postprocess_nontee(
                                postprocess_compare_local_arr,
                                workflow_prov,
                                post_err_log,
                                post_stat_log)
                            str_out= "*" * 10 + postprocess_compare_arr[file_index][0].upper() + " FINISHED" + "\t" + str( time.strftime("%Y/%m/%d %H:%M:%S",time.localtime())) + "*" * 10 +"\n"
                            yap_file_io.write_data(str_out,post_err_log)
                            yap_file_io.write_data(str_out,post_stat_log)
                    comm.barrier()
                comm.barrier()

            if myrank == 0:
                print "-"*20, "Postprocess End Time ", time.strftime("%Y/%m/%d %H:%M:%S", time.localtime()), "-"*20
                if wd.run_postprocess_analysis == "yes":
                    yap_log.combine_temp_log(wd.err_log_path, wd.stat_log_path)

            comm.barrier()
        except Exception as e:
        	print "yap.py: Encountered error in postprocess flow"
		print "Error", e
    else:
        print "Postprocess not executed owing to errors in the PREPROCESS/ALIGNMENT step.\n"
    #write provenance data from postprocess multiple sample commands 
    for k in range(0, wd.nprocs):
        if myrank == k:
            for gg in (workflow_prov):
            	yap_file_io.write_data("\t"+gg+"\n",f_summary_file)
        comm.barrier()
    #create summary input file listing yap workflow output path, used by yap_summary_call script	
    if myrank == 0:
        if wd.job_id !='':
            summary_input = 'yap_summary_input.' + str(wd.job_id)
            fhr = open(summary_input, 'a')
            fhr.write("-i " + workflow_output_path + " -o " + workflow_output_path +"\n")
            fhr.close()
	str_out = "-"*20 + "Analysis End Time For Workflow : " + wd.workflow_name + " " + time.strftime("%Y/%m/%d %H:%M:%S", time.localtime()) + "-"*20
	print str_out
        yap_file_io.write_data(str_out +"\n",f_summary_file)
	'''create symlinks from aligner output into consolidated output directory,
		useful to trigger postprocess-only job in the next run'''	
        os.system(wd.yap_src + "/" + "yap_symlinks.py -i " + workflow_output_path + " -o " + wd.consolidated_output)
        # write usage log
        number_of_samples = len(wd.inp_files_list)
        if os.path.exists(wd.usage_log):
        	f_usage = open(wd.usage_log, 'a')
        	if wd.job_id == '':
                	job_id=wd.random_id
                f_usage.write(wd.yap_version +"\t" +wd.username +"\t" +time.strftime("%Y/%m/%d %H:%M:%S",time.localtime()) +"\t" + str(number_of_samples) +"\t" +str(wd.nprocs) +"\t" +str(wd.job_id) +"\t" +"END" +"\n")
                f_usage.close()
        else:
                "Warning: YAP usage log not found"
    comm.barrier()
comm.barrier()
if myrank == 0:
    #clean the temporary files
    yap_tools.file_cleanup()
    #write samples pass/fail log 
    yap_log.pass_fail_matrix()
    print "-"*20 + "YAP Analysis End Time : " + time.strftime("%Y/%m/%d %H:%M:%S", time.localtime()) + "-"*20
MPI.Finalize()
