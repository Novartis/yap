# Merging multiple lane samples (Regrouping)

YAP always uses the first pair in a paired end sample to identify that sample. Keeping this in mind while merging multiple lane samples.

* Create a file called regroup.cfg if there isn't already one
* In the file follow a replacement sample name followed by a space followed by the old sample names separated with commas.

```
new_sample_name file1_1,file2_1,file3_1
```

* The above will merge the bam/sam files of file1, file2, file3 and call them _"new_sample_name"_
* Every line represents a distinct sample group.
* Once you have created the file, you need to update the appropriate workflow_configuration parameters. Check [workflow_configuration] page to do this.
* The "yap --check" command tells you if the syntax is correct and if the correct samples are being chosen to be regrouped.

Example:

**Contents of regroup_output.cfg**
```
sample_ABC human_sample1_1,human_sample3_1
XYZ human_sample2_1,human_sample4_1
```

### Regroup Output

There will be a change in the output structure of YAP.

* The alignment and postprocess output for all samples will be stored under the regroup_output directory.
* Samples that have no replicates will automatically be written to the regroup_output directory.
* All postprocess commands will consider the respective merged alignment output.
* As preprocess is run before merging they would follow the regular YAP directory structure with regular sample directories outside the yap_output directory.

If you have 3 files, human_sample1, human_sample2 to be merged as sample_ABC and human_sample3 without replicate data, then the regroup_output will contain:

regroup_ouput/
    sample_ABC
    human_sample3_1

