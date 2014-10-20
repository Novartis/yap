# Design your own YAP workflow

## How to insert your command in YAP

Before you create a workflow, it is essential to understand that YAP accepts commands through its configuration files which follow a particular syntax.

* The command should reside between a _:begin .. :end_ segment.
* Every _:begin_ marks the start of a command and _:end_ marks the end of that command.
* Every option-parameter couplet is accepted in a _"keyword" := "value"_ format.
* Every command needs a switch keyword i.e. _"execute_command" := "yes" # or no_
* The path to the executable is provided as a value to the keyword _"command_name"_.

Lets take an example.

**Command on the command line**

`bowtie hg19 -q -v 2 -k 10 -m 10 --best -S -p 8 -1 human_sample1_1.fq -2 human_sample1_2.fq > human_sample1.sam`

**Command in the configuration**

`:begin`

`"execute_command" := "yes"`

`"command_name" := "bowtie"`

`"" := "hg19 -q -v 2 -k 10 -m 10 --best -S -p 8 -1 pipe1 -2 pipe2 >output_file.sam"`

`:end`

_Note: the terms pipe1 and pipe2 will be explained later._

## Where to insert your command in YAP?

YAP divides the processing into 3 main stages. preprocess, alignment and postprocess. The commands are accepted through stage specific configuration files. There is a main workflow configuration file which takes the metadata and controls the flow of data between the three stages.

* [workflow_configuration.cfg](workflow.md)       -       manage metadata, specify input files, paths and output directories.

* [preprocess_configuration.cfg](preprocess.md)   -       pre-alignment packages to massage your seqdata
    * Make sure that the package supports data processing in parts/chunks - parallelization.
    * Works on one file at a time (For eg: independently on both reads of a pair).

* [aligner_configuration.cfg](aligner.md)          -       bwa, bowtie, bowtie2, tophat or insert your own aligner
    * Can take both single end or paired end input at the same time.
    * Input files can be specified using _pipe1_ and in the case of paired end data, as _pipe1_ and _pipe2_.
    * Use key _output_file_ and _output_directory_ to specify the output filenames and output directories. The framework populates the full path for you based on the sample name.

* [postprocess_configuration.cfg](postprocess.md) -       postalignment packages, generate counts or metrics
    * Every command needs two additional keywords to specify the location and format of the input file.

    |field|description|
    |----------------------------------------------|---------------------------------------------------------|
    |"input_file_type"|A _wild card_ expression that obtains a list of the qualifying files. When you do choose     queryname and coordinate sorting in the workflow_configuration, Enter "*queryname.bam" and "*coordinate.bam" respectively.|
    |"input_directory"|Where do your input files reside?
    ||_aligner_output_ - corresponds to the aligner_output directory of yap for every possible sample|
    ||_postprocess_output_ - corresponds to the postprocess_output directory of yap for every possible sample|
    ||For example: Cuffdiff takes the "transcripts.gtf" from cufflinks, which is a postprocess package.|

    * Implements a UNIX tee like feature to process multiple commands on one sample.
    * Specify a set of commands within a _:begin_tee_ .. _:end_tee_ command block.

    ```
    :begin_tee
        "input_directory" := "aligner_output"
        "input_file_type" := "*coordinate*"

        :begin
        "execute_command" := "yes"
        "command_name" := "yap_exon_count -f 1.0 -exon_coordinates_file human-ucsc-final_exon_coord.bed -exon_CoordToNumber_file human-ucsc-final_exon_coord_number.bed -i - -o := output_file"
        :end

        :begin
        "execute_command" := "no"
        "command_name" := "java -Xmx4g -jar /usr/prog/picard-tools/1.89/CollectMultipleMetrics.jar VALIDATION_STRINGENCY=SILENT I=/dev/stdin O=output_file ASSUME_SORTED=True REFERENCE_SEQUENCE=hg19.fa QUIET=True"
        :end
    :end_tee
    ```

    * The _"input_file_type"_ and _"input_directory"_ specified under _:begin_tee_ is the input passed to all commands in that _tee_ code block.

**Order of execution**
```
workflow configuration -> preprocess_configuration -> aligner_configuration -> postprocess_configuration.
```
