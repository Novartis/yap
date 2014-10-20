# Preconfigured workflows 

Once you have setup YAP feel free to download and run any of our preconfigured workflows listed below.

* RNASeq QC and Counts pipeline - [Download configuration](rnseq_qc_counts.zip)

	FastQC, FastqScreen- > Bowtie ->cooridnate and queryname sort -> Picard tools, Cufflinks, and HTSeq counts

* GATK Variant calling pipeline - [Download configuration](variant_call.zip)

	BWA -> Picard MarkDuplicates.bam -> GATk - RealignerTargetCreator, IndelRealigner, CountCovariates, TableRecalibration; Picard - QualityScoreDistribution, CalculateHsMetrics, CollectTargetedPcrMetrics

* ChIP-Seq peakcalling - [Download configuration](chipseq.zip)
       
	Fastq, fastq screen ->  alignment to hg19 with BWA -> filter low map quality reads -> coordinate sort -> peak calls using MACS2

* Bacterial classification - [Download configuration](bacterial_classification.zip)
	
	Fasta files -> classfication using mothur package -> summarize the results
