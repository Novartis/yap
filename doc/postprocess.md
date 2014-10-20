# Postprocess configuration

We execute the commands in the postprocess_configuration file after alignment. Owing to a number of different alignment output formats, we have decided to introduce 2 fields.

|field|description|
|--------------------------------------------------|---------------------------------------------------------|
|"input_file_type"|A wild card expression that obtains a list of the qualifying files. When you do choose queryname and coordinate sorting in the workflow_configuration, Enter "*queryname.bam" and "*coordinate.bam" respectively.|
|"input_directory"|Where do your input files reside?
||_aligner_output_ - corresponds to the aligner_output directory of yap for every possible sample|
||_postprocess_output_ - corresponds to the postprocess_output directory of yap for every possible sample|
||For example: Cuffdiff takes the "transcripts.gtf" from cufflinks, which is a postprocess package.|

While inserting your own command make sure:

* Command line should reside between a :begin .. :end segment.
* Every :begin marks the start of a command and :end marks the end of that command.
* Every option-parameter couplet is accepted in a "keyword" := "value" format.
* Every command needs a switch keyword i.e. "execute_command" := "yes"
* You enter the executable as the value to the "command_name" keyword.
* Insert appropriate values for "input_directory" and "input_file_type"
* Pipes do not work in postprocess_configuration.

##:begin_tee vs :begin

We know that entering the input_file_type and input_directory for every command can be tedious. Hence, we have a _:begin_tee_ separator that takes a single input_file_type and input_directory and applies them to multiple commands specified within it using the UNIX Tee mechanism within YAP.

The commands below namely the YAP exon count and Picard Collect Multiple Metrics will take coordinate sorted sam files from the aligner output directory of the workflow.

```
:begin_tee        
"input_directory" := "aligner_output"
"input_file_type" := "*coordinate.sam"
       :begin
           #YAP module for exon count
           "execute_command" := "yes"
           "command_name" := "yap_exon_count"
           "" := "-f 1.0 -exon_coordinates_file ucsc-final_exon_coord.bed -exon_CoordToNumber_file human-ucsc-final_exon_coord_number.bed -i - -o output_file"        
       :end
	
       ## Picard Collect Multiple Metrics
       :begin
           "execute_command" := "no"
           "command_name" := "java -Xmx4g -jar CollectMultipleMetrics.jar"
           "" := "VALIDATION_STRINGENCY= SILENT I= /dev/stdin O= output_file ASSUME_SORTED= True REFERENCE_SEQUENCE= hg19.fa QUIET= True"
       :end
:end_tee
```

_Note: Indentation is solely for readability and is not a requirement_

## Multiple sample comparisons

Commands involving multiple sample comparisons such as cuffdiff or cuffmerge require the user to specify the samples necessary for the comparison. If the user has to perform more than one comparison, there needs to be more than one cuffdiff command. 

Instead, in YAP, we require an accessory file that lists each combination in a new line. We have classified it into two types and have two different keywords.

|keyword|description|
|----|----|
|list_of_samples_to_compare| Used for cuffdiff. YAP takes care of preparing the comparison sets. Possible values: "all" or full path to the comparison file.|
|list_of_samples| Used for cuffmerge and cuffcompare like commands.Possible values: "all" or full path to the comparison file|
|list_delimiter| YAP uses spaces and tabs as default delemiters. Incase your command requires another delimiter, specify with this field. |

```
For example:
human_sample1 human_sample2,human_sample3
human_sample4 human_sample5
```

_The example above will correspond to two comparison groups. ',' separated values are experimental replicates._

[Click](design_workflow.md) to go back to "Designing your own workflow"
