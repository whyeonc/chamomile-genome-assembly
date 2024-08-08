#!/bin/zsh

projectdir="/filer-dg/agruppen/dg2/cho/chamomile"

#Index construction

fasta='/filer-dg/agruppen/dg2/cho/chamomile/diploid/phased_assembly/assembly/chamomile/chamomile2.hic.p_utg.fa'
gmap_build='/opt/Bio/gmap/2019-09-12/bin/gmap_build'
dir=$fasta:h'/gmap_index'
name=${fasta:t:r}
$gmap_build $fasta -D $dir -d $name > ${name}_gmap_build.out  2> ${name}_gmap_build.err

#Alignment
gmapl='/opt/Bio/gmap/2019-09-12/bin/gmapl'
gmap_build='/opt/Bio/gmap/2019-09-12/bin/gmap_build'
dir=$fasta:h'/gmap_index'
name=${fasta:t:r}
threads=50
query='/filer-dg/agruppen/dg2/cho/chamomile/diploid/phased_assembly/guide_map/pseudomolecules_v1.fa'
prefix="${name}"

$gmapl -d $name -D $dir -t $threads -f 2 $query > ${prefix}.gff 2> ${prefix}.err  
grep mRNA $prefix.gff | cut -f 1,4,5,9 | tr ';' '\t' | cut -f 1-3,5,7-8 \
   | sed -e 's/coverage=//' -e 's/identity=//' -e 's/Name=//' > ${prefix}_table.txt 
