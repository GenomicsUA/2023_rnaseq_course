#!/bin/bash

panel=$1

SCRIPT=`readlink -f $0`
SCRIPT_PATH=`dirname $SCRIPT`

if [ ! -f $panel.all_exons.raw.bed ]
then
    Rscript ${SCRIPT_PATH}/generate_panel.R grch37 $panel.genes.csv $panel.all_exons.raw.bed 15 $panel.genes.raw.bed
fi

bedtools sort -i $panel.all_exons.raw.bed -g ${SCRIPT_PATH}/chr_n.txt | bedtools merge -c 4,5,6,7,8,9 -o distinct,collapse,distinct,distinct,collapse,collapse -delim ";" > $panel.all_exons.bed
cat $panel.all_exons.bed | awk '{print "chr"$0}' > $panel.all_exons.chr.bed

bedtools sort -i $panel.genes.raw.bed -g ${SCRIPT_PATH}/chr_n.txt > $panel.genes.bed
cat $panel.genes.bed | awk '{print "chr"$0}' > $panel.genes.chr.bed

echo "Total intervals: " `cat $panel.all_exons.chr.bed | wc -l` > $panel.summary.txt
echo "Total bp: " `cat $panel.all_exons.chr.bed | awk -F '\t' '{sum+=$3-$2+1}END{print sum}'` >> $panel.summary.txt

