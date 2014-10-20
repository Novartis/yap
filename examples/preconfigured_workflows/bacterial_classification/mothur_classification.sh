#!/bin/bash
f1=$1
outdir=$2
NSLOTS=$3
temp=$4
src=$5
count_table_path=$6
sample_name=$7
var1=`echo $outdir | awk -F'/no_barcode_specified/' '{print $1}'`
var2=`basename $var1 | awk -F'.unique' '{print $1}'`
count_table=$count_table_path/$var2.count_table

## mothur cmds
../../packages/mothur/1.33.2/mothur <<EOF
set.dir(output=$outdir,tempdefault=$temp)
align.seqs(fasta=$f1, reference=$src/ncbi-phiX.fa, processors=$NSLOTS)
screen.seqs(fasta=current, minlength=100, maxambig=0)
list.seqs(fasta=current)
remove.seqs(accnos=current, fasta=$f1)
classify.seqs(fasta=current, count=$count_table, reference=$src/SSURef_NR99_115_tax_silva_IDsonly.fasta, taxonomy=$src/SSURef_NR99_115.tax, cutoff=80, iters=1000)
remove.lineage(taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)
summary.tax(taxonomy=current, reftaxonomy=$src/SSURef_NR99_115.tax)
EOF
