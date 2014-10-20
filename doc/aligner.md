# Aligner configuration

While inserting your own command make sure
* Command line resides between a :begin .. :end segment.
* Every _:begin_ marks the start of a command and _:end_ marks the end of that command.
* Every option-parameter couplet is accepted in a "keyword" := "value" format.
* Every command needs a switch keyword i.e. _"execute_command"_ := _"yes"_
* You enter the executable as the value to the "command_name" keyword.
* Input files can be specified using pipe1 and in the case of paired end data, as pipe1 and pipe2.
* If you want the right output prefix, which is by default the basename, use the variable "_output_file_"
* _output_file_ will contain the yap output path and the appropriate basename prefix.
Eg: output_file.sam = yap_output/human_sample1_1/no_barcode_specified/
                      aligner_output/human_sample1.sam
* If the package writes out a directory like tophat does, use the variable "_output_directory_".
* _output_directory_ will point directly to the "aligner_output" directory. Eg: /clscratch/521/yap_output/human_sample1_1/no_barcode_specified/aligner_output

Lets take an example.

**Command on the command line**

`bowtie hg19 -q -v 2 -k 10 -m 10 --best -S -p 8 -1 human_sample1_1.fq -2 human_sample1_2.fq > human_sample1.sam`

**Command in the configuration**
```
:begin
"execute_command" := "yes"
"command_name" := "bowtie"
"" := "hg19 -q -v 2 -k 10 -m 10 --best -S -p 8 -1 pipe1 -2 pipe2 >output_file.sam"
:end
```
[Click](design_workflow.md) to go back to "Designing your own workflow"
