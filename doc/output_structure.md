# YAP Output structure

* Once you have run a YAP job, you need to know how your results have been saved. 
* When you list your output directory, you should see a list of directories with the names of the (first sequence in the case of of the paired end) samples you provided.
* Every sample directory contains preprocess, alignment and postprocess output directories.
* Apart from these, YAP has a few directories and files described in the table below.

```
>ls yap2.3_output
consolidated_output  human_sample1_1  human_sample2_1  human_sample3_1  human_sample4_1  sample_name_track.log  yap2.3_test_summary   yap2.3_test_workflow_summary.txt  yap_log_file
```

| Directory | Description |
|---|---|
| consolidated_ouput| When aligner is run, this directory contains symlinks to all the alignment files across samples run in that YAP job. This makes it easier to generate counts/metrics other postprocess packages even after your YAP job is complete. |
|workflow_summary.txt | Chronologically logs the different commands executed by yap, hence giving the user an easy way to reproduce their results at any given time.|
| yap2.3_test_summary| The summary folder contains summarized counts/metrics results across all the samples in one text file. Please refer to the yap_summary section for more.|
| yap_log_files| please refer to our Logging your errors section|
| sample_name_track.log| Incase of duplication of sample names YAP assigns a suffix of 001, 002 etc. This file specifies the metadata related to that.|

It must be noted that packages like cuffdiff, cuffmerge and cuffcompare will reside on this level as they generate statistics by comparing samples and hence can't be stored under the sample directories. They will be denoted as cuffdiff_output, cuffmerge_output and cuffcompare_output respectively.

Sample directories
--------------------

By default YAP names the sample directories after the first sequence of your paired end sample. Under each sample directory you have its barcode information, if you've run fastx_barcode_splitter package, or a no_barcode_specified by default.

```
> ls human_sample1_1   
no_barcode_specified
```
Under the barcode level, we have 3 directories, dedicated for storing results of the three YAP stages - preprocess, alignment and postprocess.

```
> ls 
aligner_output  postprocess_output  preprocess_output
```

```
yap_output_directory

consolidated_output 
    # Symlinks of aligner_output
    human_sample1_1_coordinate.sam  human_sample2_1_coordinate.sam
    human_sample1_1_queryname.sam   human_sample2_1_queryname.sam 

cuffdiff_output 
    cuffdiff_group_0 
         transcripts.gtf 
cuffmerge_output
 
human_sample1_1 
    aligner_output
         human_sample1_1_coordinate.sam  human_sample1_1_queryname.sam
    preprocess_output
         human_sample1_1.fq_fastqc      
         human_sample1_1.fq_screen.txt 
         human_sample1_1.fq_fastqc.zip            
         human_sample1_2.fq_fastqc
         human_sample1_1.fq_screen.png  
         human_sample1_2.fq_fastqc.zip
    postprocess_output
         genes.fpkm_tracking 
         human_sample1_1_CalculateHsMetrics.txt      human_sample1_1_CollectAlignmentSummaryMetrics.txt 
         human_sample1_1_CollectGcBiasMetrics.pdf    human_sample1_1_CollectGcBiasMetrics.txt 
         human_sample1_1_CollectGcBiasMetrics_summary.txt   human_sample1_1_CollectInsertSizeMetrics.pdf
         human_sample1_1_CollectInsertSizeMetrics.txt   human_sample1_1_CollectRnaSeqMetrics.pdf
         human_sample1_1_CollectRnaSeqMetrics.txt    human_sample1_1_CollectTargetedPcrMetrics.txt           
         human_sample1_1_MarkDuplicates.bam   human_sample1_1_MarkDuplicates.txt  
         human_sample1_1_MeanQualityByCycle.pdf   human_sample1_1_MeanQualityByCycle.txt  
         human_sample1_1_QualityScoreDistribution.pdf   human_sample1_1_QualityScoreDistribution.txt  
         human_sample1_1_htseq-count.out   human_sample1_1_yap_exon_count_exoncount.txt  
         human_sample1_1_yap_exon_count_exoncount_summary.txt isoforms.fpkm_tracking 
         skipped.gtf  transcripts.gtf 

human_sample2_1 
    aligner_output 
         human_sample2_1_coordinate.sam  human_sample2_1_queryname.sam 
    preprocess_output 
         human_sample2_1.fq_fastqc
         human_sample2_1.fq_screen.txt
         human_sample2_1.fq_fastqc.zip              
         human_sample2_2.fq_fastqc
         human_sample2_1.fq_screen.png
         human_sample2_2.fq_fastqc.zip  
    postprocess_output 
         genes.fpkm_tracking 
         human_sample2_1_CalculateHsMetrics.txt      human_sample2_1_CollectAlignmentSummaryMetrics.txt 
         human_sample2_1_CollectGcBiasMetrics.pdf    human_sample2_1_CollectGcBiasMetrics.txt 
         human_sample2_1_CollectGcBiasMetrics_summary.txt   human_sample2_1_CollectInsertSizeMetrics.pdf
         human_sample2_1_CollectInsertSizeMetrics.txt   human_sample2_1_CollectRnaSeqMetrics.pdf
         human_sample2_1_CollectRnaSeqMetrics.txt    human_sample2_1_CollectTargetedPcrMetrics.txt           
         human_sample2_1_MarkDuplicates.bam   human_sample2_1_MarkDuplicates.txt  
         human_sample2_1_MeanQualityByCycle.pdf   human_sample2_1_MeanQualityByCycle.txt  
         human_sample2_1_QualityScoreDistribution.pdf   human_sample2_1_QualityScoreDistribution.txt  
         human_sample2_1_htseq-count.out   human_sample2_1_yap_exon_count_exoncount.txt  
         human_sample2_1_yap_exon_count_exoncount_summary.txt isoforms.fpkm_tracking 
         skipped.gtf  transcripts.gtf 

sample_name_track.log 

human_sample1_1 => /usr/prog/yap/2.3/examples/sample_input/human_sample1_1.fq,/usr/prog/yap/2.3/examples/sample_input/human_sample1_2.fq
human_sample2_1 => /usr/prog/yap/2.3/examples/sample_input/human_sample2_1.fq,/usr/prog/yap/2.3/examples/sample_input/human_sample2_2.fq 

yap2.3_test_summary 
     yap2.3_test_AlignmentSummaryMetrics_summary.txt   yap2.3_test_CalculateHsMetrics_summary.txt
     yap2.3_test_FastQC_summary_report.txt   yap2.3_test_FastqScreen_summary_report.txt
     yap2.3_test_GcBiasMetrics_summary.txt   yap2.3_test_InsertSizeMetrics_summary.txt
     yap2.3_test_MarkDuplicates_summary.txt  yap2.3_test_QScycle.txt
     yap2.3_test_QSdistribution.txt   yap2.3_test_RnaSeqMetrics_summary.txt
     yap2.3_test_TargetedPcrMetrics_summary.txt   yap2.3_test_cufflinks_summary.fpkm
     yap2.3_test_exon_count_summary.txt   yap2.3_test_htseq_summary.raw.txt
     yap2.3_test_summary.pdf  yap2.3_test_workflow_summary.txt 

yap_log_files 
  error_log 
    human_sample1_1_err.log  human_sample2_1_err.log 
  status_log 
    human_sample1_1_stat.log  human_sample2_1_stat.log   yap_pass_fail_matrix.log
```
