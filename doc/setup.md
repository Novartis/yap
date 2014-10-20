Requirements
---------------
**YAP only runs on linux!**

Dependencies: 

The following dependencies have to be first installed in your environment. Once installed, make sure these dependencies are added to your path.

* Recent versions of gcc (gcc 4.8.x is well tested)
* Python 2.7.7
* Openmpi 1.6.5
* Python modules:
    * MPI4py - 1.3
    * PyPdf - 1.13
    * Numpy - 1.7.1
    * netsa-utils - 1.4.3
* bedtools - 2.15.0
* samtools - 0.1.18
	

System Configuration:

YAP provides a framework to run external tools and data, so the tools used in the workflows drive the system requirements. It can be installed on multicore linux workstation with a decent amount of memory for small data, or on large cluster systems to scale optimally for large data processing. The framework has been tested extensively for NGS data on clusters with minimum system configuration of 8-12 cores and 24-48Gb memory.


YAP Setup
----------

* Download the yap source from [here](https://github.com/Novartis/yap/archive/master.zip)
* Uncompress the source directory 

	for example: uncompress the directory as "/home/packages/YAP"

* Set YAP_HOME environment variable to the source directory.

	$ export $YAP_HOME=/home/packages/YAP

* Add bin directory to path

			$ export PATH=$PATH:$YAP_HOME/bin/
				
* Set YAP_LOCAL_TEMPDIR  environment variable for temporary computation. For optimum performance point this directory to a location which is  local to the machine. 

			$ export YAP_LOCAL_TEMPDIR=/scratch/username/yap_temp

**Verification**

	$ echo YAP_HOME
		/home/packages/YAP

	$ echo $YAP_LOCAL_TEMPDIR
		/scratch/username/yap_temp

	$ which yap
		/home/packages/YAP/bin/yap

Demo Run
---------

The following section is meant to be interactive and hence you would need Linux account access and install the necessary packages.

* [Running a YAP job](run_job.md)
* [YAP - Output Structure](output_structure.md)

Design your own Workflow
------------------------

* Learn to create your own workflow using YAP framework - click [here](design_workflow.md)

Preconfigured Workflows
-----------------------

* Try preconfigured workflows which are set-up to do specific analysis, these are designed with the help of NIBR Analysts - click [here](preconfigured_workflow.md)
* For a full list of packages that have been tested with YAP - click [here](packages.md)
