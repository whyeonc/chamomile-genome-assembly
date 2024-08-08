#!/bin/zsh

projectdir="/filer-dg/agruppen/dg2/cho/chamomile"

reads=$projectdir/hifi_reads

prefix=${projectdir:r}/chamomile

find $reads | grep 'fastq.gz$' | xargs hifiasm -t 50 -o $prefix  2> ${prefix}_asm.err
