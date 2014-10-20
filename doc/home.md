What is YAP?
-------------

YAP is an extensible parallel framework, written in python using openmpi libraries. It allows researchers to quickly build high throughput big data pipelines without extensive knowledge of parallel programming. The user interacts with the framework through simple configuration files to capture analysis parameters and user directed metadata, enabling reproducible research. Using YAP, analysts have been able to achieve a significant speed up of upto 36x in RNASeq workflow execution times.

## Data processing
The scalability in YAP is achieved  using data parallel paradigms. The input data is split into a number of chunks, where the chunk size in megabytes is user defined. The chunks are then distributed evenly to a number of processors which then process them in parallel. Synchronization points are identified and is taken care of automatically by the package. This has a two fold purpose. By splitting the data into smaller chunks it is possible to run it on hardware with small amounts of RAM. This is especially important as the design philosophy of YAP is to do as much work as possible once the sequences are read into memory.

The data is processed through three main stages Preprocess, Alignment and Postprocess, with appropriate checkpoints. The user can restart the analysis from any particular stage. Each stage can contain multiple packages/software/scripts. 

**_Fig A_: Shows the order of processing. _Fig B_: Depicts the algorithm that is used by YAP.**
![Image](flowchart.png "(A) Data flow. (B) Dataprocessing")

## [Setup and Requirements](setup.md)

## Preconfigured Workflows
 
* Try preconfigured workflows which are set-up to do specific analysis, these are designed with the help of Analysts - click [here](preconfigured_workflow.md)
* For a full list of packages that have been tested with YAP - click [here](packages.md) 
