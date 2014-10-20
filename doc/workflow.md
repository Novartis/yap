# Workflow configuration

You can think of the workflow_configuration as the head configuration file. You can enter as many workflows as you want within this file, separated by ":begin .. :end" segments. In this file, as will be illustrated below, you will be able to enter your general metadata and your workflow specific metadata, as well as workflow details and parameters necessary to run your workflow. YAP will use this information to accept your input and name your results in an intelligible manner.

### General outline of a workflow configuration

```
general metadata for your full YAP job to be entered in a "key" := "value" syntax
:begin
workflow1 metadata
workflow1 parameters
:end 
:begin
workflow2 metadata
workflow2 parameters
:end
```

## Cheatsheet

Now we encourage you to open a workflow_configuration file for reference. In the table below we highlight the keys required per workflow.

### General metadata

|Keys|Definition|
|-------------------------------------------------|----------------------------------------------------------|
|"comment" := "120 chars" |Are all general metadata variables. This section is flexible and hence you can edit and add your own metadata in a "key" := "value" format.|
|"analyst_name" := "120 chars"|Example: "comment" := "yap_test"|
|"organization_name" := "NIBR"|"analyst_name" := "analyst1" etc.|

### Workflow specific metadata

|Keys|Definition|
|-------------------------------------------------------|----------------------------------------------------------|
|_:begin_|marks the start of a workflow|
|"instrument_type" := "Illumina"|Enter workflow specific metadata. Flexible field and hence you can edit or add your own metadata Example: "instrument_type" := "Illumina"|
|"specimen_info" := "[tissue type]"|"specimen_info" := "Ovarian"|
|"seq_type" := "[DNA or RNA]"|"seq_type" := "RNA"|
||"dataset" := "Human body map"|

### Required workflow variables

|Keys|Definition|
|-------------------------------------------------------|----------------------------------------------------------|
|"workflow_type" := "rnaseq"|Type of analysis you are performing|
|"workflow_name" := "yap2.3_test"|Name of your YAP job. The output folder is going to be named based on this.|
|"input_files_path" := "/usr/prog/yap/2.3/examples/sample_input"|Specify the full directory path that contains your samples. If you have more than one input directory location, use a ';' to separate the two paths.|
|"input_files" := "human*.fq" |Enter a wildcard that captures your set of input files or list them out separated by commas. If you specify more than one input path, then you need to correspondingly separate this field by ';' s as well.|

```
Example:
"input_files_path" := "yap/2.3/examples/sample_input;alt_path/human_body_map_trunc/"
"input_files" := "human_sample1*.fq,human_sample2*.fq;s_3_*_20000000.fq*"
```
The above combination will capture all qualifying files specific to the wildcard expression, similar to an 'ls' in UNIX.

|Keys|Definition|
|-------------------------------------------------------|----------------------------------------------------------|
|"input_file_format" := "fastq"|Specify your input file format. fastq, qseq, fasta etc|
|"paired_end" := "yes"|"yes" for paired end, "no" for single end|
|"data_distribution_method" := "chunk_based"|There are two ways YAP processes your data, we either chunk your seqfiles trim, filter, align these chunks and then merge them into a sam/bam file or [In the case of tophat] keep them as it is. Hence, you have to choose between "file_based" or "chunk_based", depending on what you want to run.|
```
We recommend,
"data_distribution_method" := "file_based" for tophat
"data_distribution_method" := "chunk_based" for the rest of the yap recognized aligners. 

If you have an aligner of your own, we encourage you to see if splitting the file into smaller chunks will interfere with the integrity of your data.
```

|Keys|Definition|
|-------------------------------------------------------|----------------------------------------------------------|
|"output_files_path" := "./"|Specify your output path where the YAP job results can be written to. './' represents the current directory and is recognized in YAP.|

### Merging multiple lane samples.

Please read the section on [Merging multiple lane samples](regroup.md)

|Keys|Definition|
|-------------------------------------------------------|----------------------------------------------------------|
|"regroup_output" := "no"|If you want merge samples run on multiple lanes, set this option as yes.|
|"regroup_file" := "regroup_sample.cfg"|If you set regroup_output to yes, then you have to specify the (first sequence of the pair's name) samples that need to be merged and what to call the merged samples as in a configuration and enter the name of the configuration file here.|

### Preprocessing

|Keys|Definition|
|-------------------------------------------------------|----------------------------------------------------------|
|"preprocess_output_file_format" := "fastq"|Specify the format of the file after preprocessing steps.|
|"run_preprocess_analysis" := "yes"|If you want to run preprocessing at all.|
|"preprocess_configuration_file" := "preprocess_configuration.cfg"|The name of the configuration file containing your preprocessing command in the syntax specified above.|
|"write_preprocessed_data" := "no"|Do you want have a copy of your preprocessed data?(yes/no)|

### Alignment
|Keys|Definition|
|-------------------------------------------------------|----------------------------------------------------------|
|"run_reference_alignment" := "yes"|"yes" if you are running any alignment steps. _Tip: You cannot run only preprocessing and postprocessing steps without running an alignment for obvious reasons, YAP will throw an error in such a case._|
|"aligner_configuration_file" := "bowtie_1.0.0_configuration.cfg"|Specify the configuration file that contains your commands. Or choose from your preset [aligner configuration] files.|
|"alignment_sort_order" := "both"|What order do you want your alignment file to be sorted in. Choose from: "queryname", "coordinate", "both" or "unsorted". You would need to choose "both" if you are running HTSeq and Picard for instance. Tip: Tophat already gives a coordinate sorted output, you are better off leaving it unsorted|
|"merge_alignment_output" := "yes"|In the case of chunk_based data distribution, if you want a single merged alignment output file, choose "yes". Else, if you have a postprocess package that can work on chunked data, (like DMP's EQP Pipeline) you can change it to "no".|

### Postprocess

|Keys|Definition|
|-------------------------------------------------------|----------------------------------------------------------|
|"run_postprocess_analysis" := "yes"|"yes" to run postprocess.|
|"postprocess_configuration_file" := "postprocess_configuration.cfg"|Specify the configuration file containing the postprocessing commands.|
|:end|Marks the end of the workflow|

[Click](design_workflow.md) to go back to "Designing your own workflow"
