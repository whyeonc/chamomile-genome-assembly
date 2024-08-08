#https://bitbucket.org/tritexassembly/tritexassembly.bitbucket.io/src/master/R/
source('/filer-dg/agruppen/seq_shared/mascher/code_repositories/triticeae.bitbucket.io/R/pseudomolecule_construction.R')

#read chromosome lengths of guide map
fread('/filer-dg/agruppen/dg2/cho/chamomile/diploid/phased_assembly/guide_map/pseudomolecules_v1.fa.fai', sel=1:2, col.names=c("chr", "len"))->fai

#read centromere positions
fread('/filer-dg/agruppen/dg2/cho/chamomile/diploid/phased_assembly/pseudomolecules_v1_centromere_positions.tsv')->cen
setnames(cen, c("chr", "cen_pos"))

#read Hi-C mapping output (list of read pairs connecting RE fragments)
fread('/filer-dg/agruppen/dg2/cho/chamomile/diploid/phased_assembly/mapping/hic/hic_fragment_pairs.tsv.gz')->f
setnames(f, c("ctg1", "pos1", "ctg2", "pos2"))

#read contig lengths of assebmbly
faif <- '/filer-dg/agruppen/dg2/cho/chamomile/diploid/phased_assembly/assembly/chamomile/chamomile2.hic.p_utg.fa.fai'
fread(faif, sel=1:2, col.names=c("contig", "contig_length"))->fb

#read GMAP alignments of guide map genes
fread('/filer-dg/agruppen/dg2/cho/chamomile/diploid/phased_assembly/assembly/chamomile/gmap_index/chamomile2.hic.p_utg_chamomile_table.txt')->x
setnames(x, c("contig", "start", "end", "transcript", "alnlen", "id"))

#read position of genes in guide map assembly
fread('/filer-dg/agruppen/dg2/cho/chamomile/diploid/phased_assembly/guide_map/pseudomolecules_v1_singlecopy_100bp.bed')->p
setnames(p, c("chr", "chr_start", "chr_end", "transcript"))

#merge mapping with gene positions
p[x, on="transcript"] -> px
#filter for high-quality alignments
px[alnlen >= 90 & id >= 97]->a
#only consider genes with up to 2 alignments (one to each haplotype)
px[a[, .N, key=transcript][N <= 2]$transcript, on="transcript"]->aa

#keep contigs with at least 4 aligned genes, at least 75 % of which come from the same chromosome
aa[, .N, key=.(chr, contig)][, p := N/sum(N), by=contig][order(-p)][!duplicated(contig)] -> cc
cc[N >= 4 & p >= 0.75] -> cc

#merge with contig lengths
fb[cc, on="contig"] -> cc
cc[contig_length>=300000]->cc

#check proportion of anchored sequence
sum(cc$contig_length)/1e6

#calculate approximate chromosomal locations
aa[cc[, .(contig, chr)], on=c("chr", "contig")][, .(pos=median(as.numeric(chr_start)), pos_mad=mad(as.numeric(chr_start))), key=contig][cc, on="contig"] -> cc
cc[, mr := pos_mad / contig_length]
setorder(cc, chr, pos) 
cc[, chr_idx := 1:.N, by=chr]
cc[, agp_pos := c(0, cumsum(contig_length[-.N])), by=chr]

#read coverage information and convert to correct assembly coordinates
fread(cmd="grep '^S' /filer-dg/agruppen/dg2/cho/chamomile/diploid/phased_assembly/assembly/chamomile/chamomile2.hic.p_utg.noseq.gfa | cut -f 2,5 | tr ':' '\\t' | cut -f 1,4", head=F) -> cov
setnames(cov, c("contig", "cc"))
cov[, contig := sub('l$', '', sub("utg0*", "contig_", contig))]
cc[cov, on="contig"] ->  cc

#### Get the distribution of coverages, infer 1X and 2X boundaries

pdf("coverage_chr.pdf", height=10, width=10)
pcc1<-cc[chr!="<NA>"]
plot(density(pcc1$cc),xlim=c(0,100))
dev.off()

##
# Import data and create assembly object
##

#read guide map gene-based guide map ("pseudo-POPSEQ")
readRDS('/filer-dg/agruppen/dg2/cho/chamomile/phased_assembly/guide_map/chamomile_pseudopopseq.Rds') -> pg

#read unitig lengths
f <- '/filer-dg/agruppen/dg2/cho/chamomile/diploid/phased_assembly/assembly/chamomile/chamomile2.hic.p_utg.fa.fai'
fread(f, head=F, select=1:2, col.names=c("scaffold", "length"))->fai

#read guide map gene alignment (GMAP output) and merge with pseudo-POPSEQ table
fread('/filer-dg/agruppen/dg2/cho/chamomile/diploid/phased_assembly/assembly/chamomile/gmap_index_v3/chamomile2.hic.p_utg_pseudomolecule_v3_table.txt')->x
setnames(x, c("contig", "start", "end", "transcript", "alnlen", "id"))
fai[, .(scaffold, scaffold_length=length)][x[, .(css_contig=transcript, scaffold=contig, pos=start)], on="scaffold"] -> aln
pg[aln, on="css_contig"] -> aln

#read Hi-C map
dir <- '/filer-dg/agruppen/dg2/cho/chamomile/diploid/phased_assembly/mapping/hic'
fread(paste('find', dir, '| grep "fragment_pairs.tsv.gz$" | xargs zcat'),
      header=F, col.names=c("scaffold1", "pos1", "scaffold2", "pos2"))->fpairs

#initialize assembly object, anchor to guide map and calculate physical coverage
init_assembly(fai=fai, cssaln=aln, fpairs=fpairs) -> assembly
anchor_scaffolds(assembly = assembly, popseq=pg, species="chamomile") -> assembly
add_hic_cov(assembly, binsize=1e4, binsize2=1e6, minNbin=50, innerDist=3e5, cores=20)->assembly

#save uncorrected assembly object
saveRDS(assembly, file="chamomile_assembly.Rds")
