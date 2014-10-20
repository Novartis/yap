# Running a YAP job

Once you've set your environment, it is best to run a quick demo job to get the feel of running YAP. The following section is meant to be interactive and hence you would need Linux account access and access to the cluster.

Download the demo configuration files from [here]

There are 3 stages in YAP - Preprocess, Alignment and Postprocess. You can have command level control of these three stages in the namesake configuration files and a workflow level control in the workflow_configuration.

|configuration|purpose|
|----------|------------|
|aligner_configuration|bwa, bowtie, bowtie2, tophat or insert your own aligner|
|postprocess_configuration|postalignment packages, generate counts or metrics|
|preprocess_configuration|pre-alignment packages to massage your seqdata|
|workflow_configuration|manage metadata, specify input files, paths and output directories|
|yap_sge|submitting your job to the cluster|

The demo runs a RNASeq QC and counts workflow consisting of:
* Preprocess: FastQC, Fastqscreen
* Alignment: Bowtie, both queryname and coordinate sorted
* Postprocess: yap junction and exon counts, Picard tools (PostQC), HTSeq (Raw counts) and Cufflinks (normalized counts)

We run this workflow on 2 nodes on the UGE cluster.

To run the yap_demo job, we next need to check to see if our configuration files are correct using the command.

```
 cd <your_working directory>
yap --check workflow_configuration.cfg 
```

The "yap --check"  command checks to see

* If all paths specified are valid
* If YAP finds the appropriate input files
* Checks for syntax errors
* Lists commands to be executed.
* Gives section-wise error/warning report.

The next step is to run the YAP job.

```
mpirun -n <number_of_cores> yap workflow_configuration.cfg
```

If you have a SGE environment. Pass the number of slots variable $NSLOTS.
