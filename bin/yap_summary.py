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


import sys
import os
import glob
import re
from subprocess import *
from yap_file_io import *
import time
"""
Generates summary file which is a cumulation of the log file and all reports 
and graphs generated in the postprocess. Generates special summary tables for 
htseq, fastqc, fastqscreen
"""

sample_dict = {}

#Checks the number of system arguments.
if len(sys.argv) < 5:
    print "command line : ", "yap_summary -i input_directory_name[full qualified path to yap output directory]  \
           -o summary_output_path[where you want to create yap_summary output] "
    exit()


# Removes trailing '\n'. Checks for the validity of the file paths provided.
# Removes trailing '/' of file path
if os.path.exists(sys.argv[2]):
    outputdir = sys.argv[2].strip("\n")

    if outputdir[-1] == '/':
        outputdir = outputdir[0:-1]

    path, dir_name = os.path.split(outputdir)

else:
    print "ERROR: Input path specified doesn't exist."
    exit()


print "Beginning YAP Summary", time.strftime("%Y/%m/%d %H:%M:%S", time.localtime())
summary_dir = sys.argv[4] + "/" + dir_name + "_summary"

# Remove, previous summary_dir, make a new summary directory and a temp
# directory in it only if the path exists
if os.path.exists(summary_dir):
    prm = Popen("rm -r " + summary_dir, shell='True').wait()
pmkdir = Popen("mkdir " + summary_dir, shell='True').wait()

# print glob.glob(sys.argv[4]+'/*')
temp_dir = summary_dir + "/" + "temp"
ptemp = Popen("mkdir " + temp_dir, shell='True').wait()
summary_dir = glob.glob(summary_dir)[0]
# Searches output directory at one level
sample_list = glob.glob(outputdir + "/*")
regroup_list = []


for i in sample_list:
    if i.split('/')[-1] == "regroup_output":
        regroup_list = glob.glob(i + "/*")

sample_list = sample_list + regroup_list

# Delete all sample directories with no output.
for i in sample_list:
    if os.path.isdir(i) and os.path.split(i)[-1] != 'consolidated_output' and \
        os.path.split(i)[-1] != 'yap_log_files' and \
        os.path.split(i)[-1] != dir_name + '_summary':

        for j in glob.glob(i + '/*'):
            if os.path.isdir(j):

                if os.path.exists(j + "/preprocess_output") and \
                    os.path.exists(j + "/aligner_output") and \
                    os.path.exists(j + "/postprocess_output"):

                    if glob.glob(j + "/preprocess_output/*") == [] and \
                        glob.glob(j + "/aligner_output/*") == [] and \
                        glob.glob(j + "/postprocess_output/*") == []:
                        os.system('rm -r ' + i)

# If all sample directories are empty, script exits out
if len(sample_list) == 0:
    print "Error : Opening the directory,", outputdir, " No such input directory "
    print "please provide the YAP output directory path as input"
    print "command line : ", "yap_summary output_directory_name"
    exit()
print "Creating YAP Summary output...."


# Declares a list for each tool/software, for which, 
# summary files are going to be generated
fastqc_list = []
fastq_screen_list = []
summary_file = ''
htseq_list = []
rna_bias_list = []
gc_bias_list = []
insert_size_list = []
qs_cycle_list = []
qs_distribution_list = []
mark_dup_list = []
align_sum_list = []
exon_count_list = []
junc_count_list = []
uk_junc_list = []
cufflinks_list = []
hs_list = []
target_pcr_list = []
#define variables for eqp ( specific to in-house function, neglect otherwise)
eqp_gene=[]
eqp_junc=[]
eqp_exon=[]

# Fetches list of inputs, barcodes, htseq, fastqc and fastq screen files
# with full path and appends them to their respective lists.
# Else, seeks out the workflow_summary.txt to get the provenance.
for i in range(len(sample_list)):
    if glob.glob(sample_list[i] + '/*/postprocess_output/'):
        barcode_list = glob.glob(sample_list[i] + "/*")
        sample_dict[str(sample_list[i])] = barcode_list
        if glob.glob(sample_list[i] + '/*/postprocess_output/*htseq-count.out'):
            htseq_list.append(
                glob.glob(
                    sample_list[i] +
                    '/*/postprocess_output/*htseq-count.out'))
        if glob.glob(sample_list[i] + '/*/postprocess_output/*CollectRnaSeqMetrics.txt'):
            rna_bias_list.append(
                glob.glob(
                    sample_list[i] +
                    '/*/postprocess_output/*CollectRnaSeqMetrics.txt'))
        if glob.glob(sample_list[i] + '/*/postprocess_output/*CollectInsertSizeMetrics.txt'):
            insert_size_list.append(
                glob.glob(
                    sample_list[i] +
                    '/*/postprocess_output/*CollectInsertSizeMetrics.txt'))
        if glob.glob(sample_list[i] + '/*/postprocess_output/*GcBiasMetrics_summary.txt'):
            gc_bias_list.append(
                glob.glob(
                    sample_list[i] +
                    '/*/postprocess_output/*GcBiasMetrics_summary.txt'))
        if glob.glob(sample_list[i] + '/*/postprocess_output/*MeanQualityByCycle.txt'):
            qs_cycle_list.append(
                glob.glob(
                    sample_list[i] +
                    '/*/postprocess_output/*MeanQualityByCycle.txt'))
        if glob.glob(sample_list[i] + '/*/postprocess_output/*CollectTargetedPcrMetrics.txt'):
            target_pcr_list.append(
                glob.glob(
                    sample_list[i] +
                    '/*/postprocess_output/*CollectTargetedPcrMetrics.txt'))
        if glob.glob(sample_list[i] + '/*/postprocess_output/*CalculateHsMetrics.txt'):
            hs_list.append(
                glob.glob(
                    sample_list[i] +
                    '/*/postprocess_output/*CalculateHsMetrics.txt'))
        if glob.glob(sample_list[i] + '/*/postprocess_output/*QualityScoreDistribution.txt'):
            qs_distribution_list.append(
                glob.glob(
                    sample_list[i] +
                    '/*/postprocess_output/*QualityScoreDistribution.txt'))
        if glob.glob(sample_list[i] + '/*/postprocess_output/*MarkDuplicates.txt'):
            mark_dup_list.append(
                glob.glob(
                    sample_list[i] +
                    '/*/postprocess_output/*MarkDuplicates.txt'))
        if glob.glob(sample_list[i] + '/*/postprocess_output/*CollectAlignmentSummaryMetrics.txt'):
            align_sum_list.append(
                glob.glob(
                    sample_list[i] +
                    '/*/postprocess_output/*CollectAlignmentSummaryMetrics.txt'))
        if glob.glob(sample_list[i] + '/*/postprocess_output/*exoncount_summary.txt'):
            exon_count_list.append(
                glob.glob(
                    sample_list[i] +
                    '/*/postprocess_output/*exoncount_summary.txt'))
        if glob.glob(sample_list[i] + '/*/postprocess_output/*junctioncount_summary.txt'):
            junc_count_list.append(
                glob.glob(
                    sample_list[i] +
                    '/*/postprocess_output/*junctioncount_summary.txt'))
        if glob.glob(sample_list[i] + '/*/postprocess_output/*unknown_junction.txt'):
            uk_junc_list.append(
                glob.glob(
                    sample_list[i] +
                    '/*/postprocess_output/*unknown_junction.txt'))
        if glob.glob(sample_list[i] + '/*/postprocess_output/genes.fpkm_tracking'):
            cufflinks_list.append(
                glob.glob(
                    sample_list[i] +
                    '/*/postprocess_output/genes.fpkm_tracking'))
        #fetch files specific to eqp ( in-house functionality, neglect otherwise)
        if glob.glob(sample_list[i]+'/*/postprocess_output/*merge_eqp_counts-exon.cnt'):
                eqp_exon.append(glob.glob(sample_list[i]+'/*/postprocess_output/*merge_eqp_counts-exon.cnt'))

        if glob.glob(sample_list[i]+'/*/postprocess_output/*merge_eqp_counts-gene.cnt'):
                eqp_gene.append(glob.glob(sample_list[i]+'/*/postprocess_output/*merge_eqp_counts-gene.cnt'))

        if glob.glob(sample_list[i]+'/*/postprocess_output/*merge_eqp_counts-junction.cnt'):
                eqp_junc.append(glob.glob(sample_list[i]+'/*/postprocess_output/*merge_eqp_counts-junction.cnt'))

    if glob.glob(sample_list[i] + '/*/preprocess_output/'):
        barcode_list = glob.glob(sample_list[i] + "/*")
        sample_dict[str(sample_list[i])] = barcode_list
        if glob.glob(sample_list[i] + '/*/preprocess_output/*_fastqc/summary.txt'):
            fastqc_list.append(
                glob.glob(
                    sample_list[i] +
                    '/*/preprocess_output/*_fastqc/summary.txt'))
        if glob.glob(sample_list[i] + '/*/preprocess_output/*_screen.txt'):
            fastq_screen_list.append(
                glob.glob(
                    sample_list[i] +
                    '/*/preprocess_output/*_screen.txt'))
    else:
        log_file = re.match(r'(.*)_summary.txt', sample_list[i], re.M | re.I)
        if log_file:
            summary_file = sample_list[i]

def output_summary(sample_dict):
    """ Produces a pdf file containing the txt, png result files across all samples run.
        Accepts the sample_dict {sample:barcode} and merges into one pdf. """

    complete_dict = {}
    all_pdfs = []
    main_page = temp_dir + "/" + dir_name + "_" + "front"
    fw_main_page = open(main_page + ".txt", "wb")
    fw_main_page.write("\n\n\n")
    fw_main_page.write("\t\t\t\t\tYAP Postprocess Summary Results" + "\n\n")
    fw_main_page.write("Workflow Name= " + dir_name + "\n")
    fw_main_page.write("Input Directory Path Provided= " + outputdir + "\n")
    fw_main_page.write("Summary Ouput Directory Path= " + summary_dir + "\n")
    fw_main_page.close()
    for k in range(len(sorted(sample_dict.iterkeys()))):
        filename_key = ''
        filename_key = sorted(sample_dict.iterkeys())[k]
        path, file_name = os.path.split(filename_key)
        for j in range(len(sample_dict[filename_key])):
            barcode_path = ''
            barcode_path = sample_dict[filename_key][j]
            temp_pdfs = []
            pdf_files = []
            path, barcode_name = os.path.split(barcode_path)
            aligner_dir_path = ''
            postprocess_dir_path = barcode_path + "/" + "postprocess_output"
            pdf_files = glob.glob(postprocess_dir_path + "/*.pdf")
            text_files = glob.glob(postprocess_dir_path + "/*.txt")
            sample_page = temp_dir + "/" + file_name + \
                "_" + barcode_name + "_" + "main"
            fw_sample_page = open(sample_page + ".txt", "a+")
            fw_sample_page.write("\n\n\n\n\n")
            fw_sample_page.write("\t\t\t\t\tSample : " + file_name + "\n")
            fw_sample_page.write("\t\t\t\t\tBarcode : " + barcode_name)
            fw_sample_page.close()
            if len(text_files) > 0:
                for i in range(0, len(text_files)):
                    if not text_files[i].find('exoncount.txt'):
                        if not text_files[i].find('junctioncount.txt'):
                            if not text_files[i].find('unknown_junction.txt'):
                                post_path, post_file_name = os.path.split(
                                    text_files[i])
                                post_file_name, exten = os.path.splitext(
                                    post_file_name)
                                txt_pdf = Popen(
                                    "text2pdf.py " +
                                    text_files[i] +
                                    " -o " +
                                    temp_dir +
                                    "/" +
                                    str(k) +
                                    "_" +
                                    str(j) +
                                    "_" +
                                    str(i) +
                                    ".pdf",
                                    shell='False').wait()
                                temp_pdfs.append(
                                    temp_dir +
                                    "/" +
                                    str(k) +
                                    "_" +
                                    str(j) +
                                    "_" +
                                    str(i) +
                                    ".pdf")
            txt_pdf = Popen(
                "text2pdf.py " +
                sample_page +
                ".txt -o " +
                sample_page +
                ".pdf",
                shell='False').wait()
            all_pdfs.append(sample_page + ".pdf")
            all_pdfs.extend(temp_pdfs)
            all_pdfs.extend(pdf_files)
    final_pdf_files = []
    txt_pdf = Popen("text2pdf.py " + main_page + ".txt -o " +
                    main_page + ".pdf", shell='False').wait()
    if summary_file != '':
        path, sufile = os.path.split(summary_file)
        su_file, ext = os.path.splitext(sufile)
        su_file = temp_dir + '/' + su_file
        txt_pdf = Popen(
            "text2pdf.py " +
            summary_file +
            " -o " +
            su_file +
            ".pdf",
            shell='False').wait()

        final_pdf_files.append(su_file + ".pdf")
    final_pdf_files.insert(0, main_page + ".pdf")
    final_pdf_files.extend(all_pdfs)
    final_input_pdfs = ''
    summary_final_pdf = summary_dir + "/" + dir_name + "_summary" + ".pdf"
    for i in range(0, len(final_pdf_files)):
        if len(final_input_pdfs) + len(summary_final_pdf) > 8170:
            summary_temp_pdf = summary_dir + "/temp/" + \
                dir_name + "_temp_summary_" + str(i) + ".pdf"
            txt_pdf = Popen("mergepdf.py -i " + final_input_pdfs +
                            " -o " + summary_temp_pdf, shell='False').wait()
            final_input_pdfs = summary_temp_pdf + " "
        final_input_pdfs += final_pdf_files[i] + " "
    txt_pdf = Popen("mergepdf.py -i " + final_input_pdfs +
                    " -o " + summary_final_pdf, shell='False').wait()
    ptemp = Popen("rm -r " + temp_dir, shell='True').wait()


#
# 
# 
def count_summary(count_list, file_ext):
    """ Summarizes count files, namely - HTSeq and Exon counts
        Takes the list of file paths and the desired output file extension.
        Writes out a flat file listing results across all samples. """

    count_arr = []
    count_dict = {}
    count_gene_set = set()
    sample_set = set()
    for i in range(len(count_list)):
        for j in count_list[i]:
            if j and os.path.getsize(j) > 0:
                path = j.split('/')
                sample_set.update(['/'.join(path[:-3])])
                with open(j, 'rb') as foo:
                    count_arr = foo.readlines()
                    foo.close()
                for k in count_arr:
                    temp = k.strip('\n').split('\t')
                    count_gene_set.update([temp[0]])
                    count_dict[((path[-4], path[-3]), temp[0])] = temp[1]
    bar_str = ''
    file_str = '\t'
    for k in sorted(sample_set):
        temp1 = k.split('/')
        file_str = file_str + temp1[-1]
        for v in sample_dict[k]:
            temp = v.split('/')
            bar_str = bar_str + '\t' + temp[-1]
            file_str = file_str + '\t'
    bar_str = "BARCODE" + bar_str.rstrip('\t') + '\n'
    file_str = "SAMPLE" + file_str.rstrip('\t') + '\n'
    with open(summary_dir + '/' + dir_name + file_ext + '.txt', 'a+') as ct:
        ct.write(file_str + bar_str)
        for g in sorted(count_gene_set):
            ct.write(g)
            for v in sorted(sample_dict.values()):
                for l in v:
                    path = l.split('/')
                    try:
                        ct.write('\t' + count_dict[((path[-2], path[-1]), g)])
                    except KeyError:
                        pass
            ct.write('\n')
    ct.close()


def unknown_junc_summary(uk_list, file_ext):
    """ Summarizes unknown junctions across all samples.
        Takes the list of file paths and the desired output file extension. """
    with open(summary_dir + '/' + dir_name + file_ext, 'w') as uk:
        for i in range(len(uk_list)):
            for j in uk_list[i]:
                if j:
                    uk_file_arr = []
                    path = j.split('/')
                    with open(j, 'rb') as foo:
                        uk = foo.readlines()
                        foo.close()
                file_str = path[-4] + '\n' + path[-3] + '\n'
                uk.write(file_str)
                for k in uk_file_arr:
                    uk.write(k)
                uk.write('\n')
    uk.close()


def junc_count_summary(file_list, file_ext):
    """ Summarizes count files, namely - Junction counts
        Takes the list of file paths and the desired output file extension. """
    junc = []
    file_dict = {}
    gene_set = set()
    sample_set = set()
    junc_sample_dict = {}
    for j in file_list:
        for i in j:
            file_con = []
            file_con = read_file(i)
            path = []
            path = i.split('/')
            sample_set.update(['/'.join(path[:-3])])
            sample = path[-4]
            barcode = path[-3]
            for line in file_con:
                full_key = tuple()
                junction, count, key = line.strip().rstrip('\n').split('\t')
                full_key = ((sample, barcode), key)
                gene_set.update([key])
                try:
                    if full_key in file_dict:
                        file_dict[full_key] += int(count)
                    else:
                        file_dict[full_key] = int(count)
                except KeyError as e:
                    print e
    bar_str = '\t'
    file_str = '\t'
    for k in sorted(sample_set):
        temp1 = k.split('/')
        file_str = file_str + temp1[-1]
        for v in sample_dict[k]:
            temp = v.split('/')
            bar_str = bar_str + '\t' + temp[-1]
            file_str = file_str + '\t'
            junc_sample_dict[k] = v
    bar_str = 'BARCODE' + bar_str.rstrip('\t') + '\n'
    file_str = 'SAMPLE' + file_str.rstrip('\t') + '\n'
    with open(summary_dir + '/' + dir_name + file_ext + '.txt', 'w') as hq:
        hq.write(file_str + bar_str)
        for g in sorted(gene_set):
            hq.write(g)
            for k in sorted(sample_set):
                for v in sample_dict[k]:
                    path = v.split('/')
                    try:
                        hq.write(
                            '\t' + str(file_dict[((path[-2], path[-1]), g)]))
                    except KeyError:
                        hq.write('\tN/A')
            hq.write('\n')
    hq.close()


# Summarizes cufflinks generated files across samples
# Takes the list of file paths and the desired output file extension.
# Writes out a flat file listing results across all samples. 
def cufflinks_summary(cuff_list, file_ext):
    """ Summarizes cufflinks generated files across samples.
        Takes the list of file paths and the desired output file extension. """

    cuff = []
    cuff_dict = {}
    cuff_gene_set = set()
    sample_set = set()
    cuff_sample_dict = {}
    for i in range(len(cuff_list)):
        for j in cuff_list[i]:
            if j:
                path = j.split('/')
                sample_set.update(['/'.join(path[:-3])])
                with open(j, 'rb') as foo:
                    htseq = foo.readlines()
                    foo.close()
                for k in range(1, len(htseq)):
                    temp = []
                    temp = htseq[k].strip('\n').split('\t')
                    gene_string = temp[0] + '|' + temp[6] + '|'
                    cuff_gene_set.update([gene_string])
                    cuff_dict[
                        ((path[-4], path[-3]), gene_string)] = temp[-4] + "\t" + temp[-1]
    bar_str = ''
    file_str = '\t'
    fpkm_str = ''
    for k in sorted(sample_set):
        temp1 = k.split('/')
        file_str = file_str + temp1[-1]
        for v in sample_dict[k]:
            temp = v.split('/')
            bar_str = bar_str + '\t\t' + temp[-1]
            file_str = file_str + '\t\t'
            fpkm_str = fpkm_str + '\tFPKM\tFPKM_Status'
            cuff_sample_dict[k] = v
    bar_str = 'BARCODE\t' + bar_str.rstrip('\t').lstrip('\t') + '\n'
    file_str = 'SAMPLE' + file_str.rstrip('\t') + '\n'
    fpkm_str = 'TRACKING_ID' + fpkm_str + '\n'
    with open(summary_dir + '/' + dir_name + file_ext, 'w') as hq:
        hq.write(file_str + bar_str + fpkm_str)
        for g in sorted(cuff_gene_set):
            hq.write(g)
            for k in sorted(sample_set):
                for v in sample_dict[k]:
                    path = v.split('/')
                    try:
                        hq.write('\t' + cuff_dict[((path[-2], path[-1]), g)])
                    except KeyError:
                        hq.write('\tN/A')
            hq.write('\n')
    hq.close()

# Transposes the fastqc_summ list from the fastqc_summary function. fout is the file handle of the output file.
# if Flag == True then it is treated as a title string.
def fastqc_transpose(fin, fout, flag):
    """ Transposes the fastqc_summ list from the fastqc_summary function. fout is the file handle of the output file.
        if Flag == True then it is treated as a title string. """

    file_name = fin.split(
        '/')[-5] + "-" + fin.split('/')[-2].rstrip('.fq_fastqc')
    with open(fin) as foo:
        ls = foo.readlines()
    foo.close()
    entry = dict()
    tup = ()
    sample = set()
    vals = set()
    for i in range(len(ls)):
        ls[i] = ls[i].rstrip('\n').split('\t')
        ls[i][2], junk = os.path.split(ls[i][2])
        sample.update([file_name])
        vals.update([ls[i][1]])
        tup = (file_name, ls[i][1])
        entry[tup] = ls[i][0]
    with open(fout, 'a+') as f:
        if flag == 'True':
            outstring = ''
            outstring = outstring + 'Sample_ID\t'
            for v in vals:
                outstring = outstring + v + '\t'
            f.write(outstring.rstrip('\t'))
            f.write('\n')
        for i in sample:
            out_string = ''
            out_string = out_string + i + '\t'
            for j in vals:
                out_string = out_string + entry[(i, j)] + '\t'
            f.write(out_string.rstrip('\t'))
            f.write('\n')

def fastqc_summary(sample_dict):
    """ Generates the summary files across all samples for results produced by the FastQC package. """
    count = 0
    print "Generating Fastqc summary.."
    for j in sample_dict:
        fastqc_summ = j
        for k in range(len(fastqc_summ)):
            fout = summary_dir + "/" + dir_name + "_FastQC_summary_report.txt"
            if count == 0:
                fastqc_transpose(fastqc_summ[k], fout, 'True')
                count += 1
            else:
                fastqc_transpose(fastqc_summ[k], fout, 'False')

def summary_contaminants(fin, fout, flag):
    """ Transposes the fastqscreen_summ list from the fastqscreen_summary function. fout is the file handle of the output file.
        if Flag == True then it is treated as a title string. """

    with open(fin) as foo:
        line = foo.readlines()
    foo.close()
    file = fin.split('/')[-4]
    if fin.split('/')[-3] != 'no_barcode_specified':
        file = file + ',' + os.path.split(fin)[-3]
    with open(fout, 'a+') as f:
        out_string = ''
        header = len(line)
        if line[-1].startswith('%'):
            header -= 2
        for i in range(len(line)):
            if line[i].startswith('#'):
                header -= 1
            if not line[i].startswith('#'):
                if re.search('.+\t.+\t.+', line[i]):
                    line[i] = line[i].replace('\n', '\t')
                    if line[i].startswith('Library') and flag == 'True':
                        out_string = out_string + 'Sample_ID\t'
                        out_string = out_string + \
                            line[i] * (header - 1) + '\n' + file + '\t'
                    elif line[i].startswith('Library') and flag == 'False':
                        out_string = out_string + file + '\t'
                    else:
                        out_string = out_string + line[i]
        f.write(out_string.rstrip('\t'))
        f.write('\n')

def fastq_screen_summary(sample_dict):
    """ Generates the summary files across all samples for results produced by the Fastqscreen package. """
    count = 0
    for j in sample_dict:
        count = count + 1
        for k in range(len(j)):
            fout = summary_dir + "/" + dir_name + \
                "_FastqScreen_summary_report.txt"
            if count == 1:
                summary_contaminants(j[k], fout, 'True')
            else:
                summary_contaminants(j[k], fout, 'False')


def collect_metrics_summary(mark_dup_list, file_string):
    """ Generates the summary files across all samples for results produced by the Picardtools across all samples except QS Distribution and QS cycle. """

    dict_md = {}
    title_md = ''
    sample = ''
    barcode = ''
    mark_dup_summ = ''
    list_md = []
    length_title = []
    row_summ = []
    for i in range(len(mark_dup_list)):
        for k in range(len(mark_dup_list[i])):
            with open(mark_dup_list[i][k]) as mark_fh:
                file_md = mark_fh.read()
            mark_fh.close()
            list_md = file_md.split('#')
            for l in range(len(list_md)):
                match_md = re.match(
                    r'\s*METRICS CLASS\t.*\n(.*)\n', list_md[l], re.I | re.M)
                if match_md:
                    row_summ = list_md[l][1:].rstrip(
                        '\n').strip('\t ').split('\n')
                    mark_title = []
                    mark_value = []
                    for m in range(1, len(row_summ)):
                        if m == 1:
                            mark_title = row_summ[m].split('\t')
                        else:
                            row_summ[m] = row_summ[m].rstrip('\n')
                            row_summ[m] = row_summ[m].rstrip('\t')
                            mark_value.append(
                                row_summ[m].rstrip('\t').split('\t'))
                    for n in range(len(mark_value)):
                        for j in range(len(mark_value[n])):
                            dict_md[mark_dup_list[i][k].split(
                                '/')[-4], mark_dup_list[i][k].split('/')[-3], n, mark_title[j]] = mark_value[n][j]
    mark_dup_summ = summary_dir + "/" + dir_name + \
        '_' + file_string + '_summary.txt'
    title_md = '\t'.join(mark_title)
    with open(mark_dup_summ, 'a+') as mark_fout:
        mark_fout.write("SAMPLE\tBARCODE\t" + title_md + '\n')
        for i in sorted(sample_dict.keys()):
            path, sample = os.path.split(i)
            for j in sorted(sample_dict[i]):
                path2, barcode = os.path.split(j)
                out_string = ''
                for l in range(len(mark_value)):
                    out_string = out_string + sample + '\t' + barcode + '\t'
                    for k in mark_title:
                        try:
                            out_string = out_string + \
                                dict_md[sample, barcode, l, k]

                        except KeyError:
                            out_string = out_string + 'N/A'
                        out_string = out_string + '\t'
                    out_string.rstrip('\t')
                    out_string = out_string + '\n'
                mark_fout.write(out_string.rstrip('\t'))
                mark_fout.write('\n')



def quality_score_summary(qs_list, phrase):
    """ Generates the summary files across all samples for results produced by the Picardtools - QS Distribution and QS cycle. """

    match_qs = {}
    for i in range(len(qs_list)):
        for j in range(len(qs_list[i])):
            barcode = ''
            sample = ''
            barcode = qs_list[i][j].split('/')[-3]
            sample = qs_list[i][j].split('/')[-4]
            with open(qs_list[i][j]) as qs:
                qs_file = qs.read()
            qs.close()
            match_summary_qs = re.search(
                'HISTOGRAM\t[\w\.]+\n(.*)', qs_file, re.DOTALL)
            if match_summary_qs:
                match_qs[sample, barcode] = (
                    match_summary_qs.group(1).split('\n'))
    qs_file_string = summary_dir + '/' + dir_name + '_' + phrase + '.txt'
    with open(qs_file_string, 'a+') as qs_out:
        qs_array = []
        out_string = ''
        for k in sample_dict:
            path, temp_sample = os.path.split(k)
            for v in sample_dict[k]:
                path1, temp_bar = os.path.split(v)
                out_string = out_string + temp_sample + '\t' + temp_bar + '\t'
                if (temp_sample, temp_bar) in match_qs:
                    qs_array.append(match_qs[temp_sample, temp_bar])
        qs_out.write(out_string.rstrip('\t'))
        qs_out.write('\n')
        count = 0
        while count != len(qs_array[0]):
            out_string = ''
            for i in qs_array:
                out_string = out_string + i[count] + '\t'
            qs_out.write(out_string.rstrip('\t'))
            qs_out.write('\n')
            count = count + 1
    qs_out.close()

output_summary(sample_dict)
if len(htseq_list) > 0:
    print "Generating HTSeq summary..."
    count_summary(htseq_list, '_htseq_summary.raw')
if len(exon_count_list) > 0:
    print "Generating YAP Exon counts summary..."
    count_summary(exon_count_list, '_exon_count_summary')
if len(junc_count_list) > 0:
    print "Generating YAP Junction counts summary..."
    junc_count_summary(junc_count_list, '_junction_count_summary')
if len(uk_junc_list) > 0:
    print "Generating Unknown Junction summary..."
    unknown_junc_summary(uk_junc_list, '_unknown_junction_summary')
if len(fastqc_list) > 0:
    print "Generating FastQC summary..."
    fastqc_summary(fastqc_list)
if len(fastq_screen_list) > 0:
    print "Generating Fastq Screen summary..."
    fastq_screen_summary(fastq_screen_list)
if len(rna_bias_list) > 0:
    print "Generating RNA Bias summary..."
    collect_metrics_summary(rna_bias_list, 'RnaSeqMetrics')
if len(gc_bias_list) > 0:
    print "Generating GC Bias summary..."
    collect_metrics_summary(gc_bias_list, 'GcBiasMetrics')
if len(mark_dup_list) > 0:
    print "Generating Mark Duplicates summary..."
    collect_metrics_summary(mark_dup_list, 'MarkDuplicates')
if len(insert_size_list) > 0:
    print "Generating Insert Size summary..."
    collect_metrics_summary(insert_size_list, 'InsertSizeMetrics')
if len(target_pcr_list) > 0:
    print "Generating Targeted Pcr Metrics summary..."
    collect_metrics_summary(target_pcr_list, 'TargetedPcrMetrics')
if len(hs_list) > 0:
    print "Generating Calculated HS Metrics summary..."
    collect_metrics_summary(hs_list, 'CalculateHsMetrics')
if len(align_sum_list) > 0:
    print "Generating Aligner Metrics summary..."
    collect_metrics_summary(align_sum_list, 'AlignmentSummaryMetrics')
if len(qs_distribution_list) > 0:
    print "Generating Picard-Quality Score Distribution summary..."
    quality_score_summary(qs_distribution_list, 'QSdistribution')
if len(qs_cycle_list) > 0:
    print "Generating Picard-QS_Cycle summary..."
    quality_score_summary(qs_cycle_list, 'QScycle')
if len(cufflinks_list) > 0:
    print "Generating Cufflinks summary..."
    cufflinks_summary(cufflinks_list, '_cufflinks_summary.fpkm')

#Generating summary for EQP results (meant for in-house use only, neglect otherwise)
def mapped_read_count(input,output):
        sample=input.split("/")[-1]
        barcode_list=sample_dict[input]
        for barcode in barcode_list:
                total=(0,0)
                barcode_str=barcode.split("/")[-1]
                cmd_string='cat '+barcode+'/aligner_output/*_count.log | grep -A 3 combined.sam.gz |grep reads\ are\ mapped'
                cmd_out=Popen(cmd_string,stdout=PIPE,shell=True).communicate()[0]
                read_pattern=re.findall(r"(\d+) of (\d+) reads are mapped \([\d\.]+\%\)\.",cmd_out)
                for x,v in read_pattern:
                        total =  ((total[0] + int(x)) , (total[1] + int(v)))
                output.write(sample+"\t"+barcode_str+"\t"+str(total[0])+"\t"+str(total[1])+"\n")

for i in sample_list:
    if glob.glob(i+"/*/aligner_output/*_count.log"):
            print "Generating EQP Total Mapped Read counts summary..."
            eqp_map=open(summary_dir+"/eqp_mapped_read_count_summary.txt","a")
            eqp_map.write("SAMPLE\tBARCODE\tMAPPED_READS\tTOTAL_READS\n")
            mapped_read_count(i,eqp_map)
            eqp_map.close()
if len(eqp_gene)>0:
        print "Generating EQP Gene summary..."
        count_summary(eqp_gene,'_merge_eqp_counts-gene_summary')

if len(eqp_junc)>0:
        print "Generating EQP Junction summary..."
        count_summary(eqp_junc,'_merge_eqp_counts-junction_summary')

if len(eqp_exon)>0:
        print "Generating EQP Exon summary..."
        count_summary(eqp_exon,'_merge_eqp_counts-exon_summary')
  
print "YAP summary finished!", time.strftime("%Y/%m/%d %H:%M:%S", time.localtime())
