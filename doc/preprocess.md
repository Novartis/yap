# Preprocess configuration

The preprocess commands are executed first in YAP. Any script/package that performs QC, filters or massages the raw data before alignment can be inserted here. The YAP code internally matches the pair at the end of execution before alignment.

If you want to insert your own script.
* Make sure that the package supports data processing in parts/chunks - parallelization.
* If you split input into two parts and run them in tandem, the merged output should be same as the whole file output.
* Accepts UNIX standard input and can be written to UNIX standard output.
* Works on one file at a time, even for paired end reads.
* _pipe1_ in preprocess represents the sample in a particular command. YAP iterates over the samples and substitutes pipe1 appropriately. 
* FastQC, Fastqscreen: Are the only exceptions to this rule because they work on whole files.
If these conditions aren't met, please place your command in the aligner_configuration, before the actual aligner command.


Below is a schematic showing you the conversion of a fastqscreen command to a YAP configuration

```
/usr/prog/fastx-toolkit/0.0.13/fastq_quality_filter -q 20 -p 70 -i pipe1
```

```
:begin
"execute_command" := "yes"
"command_name" := " /usr/prog/fastx-toolkit/0.0.13/fastq_quality_filter"
"" := "-q 20 -p 70 -i pipe1"
:end
```

[Click](design_workflow.md) to go back to "Designing your own workflow"
