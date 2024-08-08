#https://bitbucket.org/tritexassembly/tritexassembly.bitbucket.io/src/master/R/
.libPaths(c("/home/mascher/tmp/3.5.1", "/filer-dg/agruppen/seq_shared/mascher/Rlib/3.5.1", "/opt/Bio/R_LIBS/3.5.1"))
source('/filer-dg/agruppen/dg2/cho/chamomile/diploid/phased_assembly/bitbucket/R/pseudomolecule_construction.R')

#read chromosome lengths of Morex V3
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

#read GMAP alignments of Morex V3 genes
fread('/filer-dg/agruppen/dg2/cho/chamomile/diploid/phased_assembly/assembly/chamomile/gmap_index/chamomile2.hic.p_utg_chamomile_table.txt')->x
setnames(x, c("contig", "start", "end", "transcript", "alnlen", "id"))

#read position of genes in Morex assembly
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
sum(cc$contig_length)/1e6 ##[1] 4343.794


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

pdf("coverage_chr_new.pdf", height=10, width=10)
pcc1<-cc[chr!="<NA>"]
plot(density(pcc1$cc),xlim=c(0,100))
abline(v=c(36), col="blue")
dev.off()

pdf("coverage_all.pdf", height=10, width=10)
plot(density(cov$cc),xlim=c(0,100))
abline(v=c(36,38, 40), col="blue")
dev.off()


#merge syntenic positions with Hi-C link list
f[, .N, key=.(ctg1, ctg2)] -> v
cc[, .(ctg1=contig, chr1=chr, pos1=pos,cc1=cc)][v, on="ctg1"] -> v
cc[, .(ctg2=contig, chr2=chr, pos2=pos,cc2=cc)][v, on="ctg2"] -> v


#run PCA of each chromosomes
rbindlist(mclapply(mc.cores=9, fai[1:9, chr], function(j){
 #fill empty values in link list (with 0)
 setnames(v[chr1 == chr2 & chr1 == j][, chr2 := NULL], "chr1", "chr") -> x
 dcast(x, ctg1 ~ ctg2, value.var="N", fill = 0) -> x
 melt(x, id="ctg1", measure=setdiff(names(x), "ctg1"), variable.factor=F, value.name="N", variable.name="ctg2") -> x
 #convert link counts to log
 x[, l := log10(0.1 + N)]
 #merge with positional information 
 cc[, .(ctg1=contig, chr, pos1=pos)][x, on="ctg1"] -> x
 cc[, .(ctg2=contig, pos2=pos)][x, on="ctg2"] -> x
 fai[x, on="chr"] -> x
 cen[x, on="chr"] -> x
 #specify linear model to remove DDD and Rabl effects 
 x[, end_dist1 := pmin(pos1, len - pos1)]
 x[, end_dist2 := pmin(pos2, len - pos2)]
 x[, cen_dist1 := abs(cen_pos - pos1)]
 x[, cen_dist2 := abs(cen_pos - pos2)]
 x[, ldist := abs(pos1 - pos2)]
 x[, res := lm(l ~ log(1+ldist) * log(cen_dist1) * log(cen_dist2) * log(end_dist1) * log(end_dist2))$res]

 #run PCA on residuals
 dcast(x, ctg1 ~ ctg2, value.var="res", fill=0) -> y
 y[, ctg1 := NULL]
 #PCA of correlation matrix, e.g. SVD of co-variance of correlation of Hi-C distance matrix (seems complicated, but gives good results)
 prcomp(cor(y), scale=T, center=T)->pca
 #extract first four PCs
 data.table(contig=rownames(pca$rotation), pca$rotation[, 1:4]) -> p
 #merge PC table with positional information
 setorder(cc[, .(contig, chr, pos, agp_pos, contig_length,cc)][p, on="contig"], chr, agp_pos) -> pp
 pp[, chr := j]
})) -> pp

#add index columns for plotting purposes
pp[, idx := 1:.N]

#plot the PCA results: PC scores along the chromosomes
pdf("chamomile_haplotype_separation_HiC_all.pdf", height=4, width=10)
par(mfrow=c(2,5))
lapply(c("PC1", "PC2", "PC3", "PC4"), function(j){
 lapply(fai[1:9, chr], function(i){
  pp[chr == i, plot(agp_pos/1e6, las=1, bty='l', get(j), type='n', main=sub("chr", "", chr[1]), ylab=j,
        xlab="Hv syntenic position [Mb]")]
  pp[chr == i, lines(lwd=3, c(agp_pos/1e6, (agp_pos + contig_length)/1e6), c(get(j), get(j))), by=idx]
 })
 plot(axes=F, xlab="", ylab="", 0, type='n')
})
dev.off()



pp[cc<=40, col := "red"]
pp[cc>40,col:="blue"]

pdf("chamomile_haplotype_separation_HiC_all_pc1_2.pdf", height=10, width=10)

lapply(fai[1:9, chr], function(i){
  pp[chr == i, plot(PC2, PC1, col=col, pch=20, las=1,bty='l', type="p", main=sub("chr", "", chr[1]), ylab="PC1",
            xlab="PC2")]
  legend("topleft", pch=19, bty='n', legend=c("40", ">40"), col=c("red", "blue"))
 })

dev.off()

pdf("chamomile_haplotype_separation_HiC_all_pc2_3.pdf", height=10, width=10)

lapply(fai[1:9, chr], function(i){
  pp[chr == i, plot(PC3, PC2, col=col, pch=20, las=1,bty='l', type="p", main=sub("chr", "", chr[1]), ylab="PC2",
            xlab="PC3")]
  legend("topleft", pch=19, bty='n', legend=c("40", ">40"), col=c("red", "blue"))
 })

dev.off()

pdf("chamomile_haplotype_separation_HiC_all_pc3_4.pdf", height=10, width=10)

lapply(fai[1:9, chr], function(i){
  pp[chr == i, plot(PC4, PC3, col=col, pch=20, las=1,bty='l', type="p", main=sub("chr", "", chr[1]), ylab="PC3",
            xlab="PC4")]
  legend("topleft", pch=19, bty='n', legend=c("40", ">40"), col=c("red", "blue"))
 })

dev.off()


#run PCA of each chromosomes for filter cc>40
rbindlist(mclapply(mc.cores=1, fai[1:9, chr], function(j){
 #fill empty values in link list (with 0)
 setnames(v[chr1 == chr2 & chr1 == j & cc1<=40 &cc2<=40][, chr2 := NULL], "chr1", "chr") -> x
 dcast(x, ctg1 ~ ctg2, value.var="N", fill = 0) -> x
 melt(x, id="ctg1", measure=setdiff(names(x), "ctg1"), variable.factor=F, value.name="N", variable.name="ctg2") -> x
 #convert link counts to log
 x[, l := log10(0.1 + N)]
 #merge with positional information 
 cc[, .(ctg1=contig, chr, pos1=pos)][x, on="ctg1"] -> x
 cc[, .(ctg2=contig, pos2=pos)][x, on="ctg2"] -> x
 fai[x, on="chr"] -> x
 cen[x, on="chr"] -> x
 #specify linear model to remove DDD and Rabl effects 
 x[, end_dist1 := pmin(pos1, len - pos1)]
 x[, end_dist2 := pmin(pos2, len - pos2)]
 x[, cen_dist1 := abs(cen_pos - pos1)]
 x[, cen_dist2 := abs(cen_pos - pos2)]
 x[, ldist := abs(pos1 - pos2)]
 x[, res := lm(l ~ log(1+ldist) * log(cen_dist1) * log(cen_dist2) * log(end_dist1) * log(end_dist2))$res]

 #run PCA on residuals
 dcast(x, ctg1 ~ ctg2, value.var="res", fill=0) -> y
 y[, ctg1 := NULL]
 #PCA of correlation matrix, e.g. SVD of co-variance of correlation of Hi-C distance matrix (seems complicated, but gives good results)
 prcomp(cor(y), scale=T, center=T)->pca
 #extract first four PCs
 data.table(contig=rownames(pca$rotation), pca$rotation[, 1:4]) -> p
 #merge PC table with positional information
 setorder(cc[, .(contig, chr, pos, agp_pos, contig_length,cc)][p, on="contig"], chr, agp_pos) -> pp
 pp[, chr := j]
})) -> pp1

#add index columns for plotting purposes
pp1[, idx := 1:.N]


pdf("chamomile_haplotype_separation_HiC_cc32_pc1_2.pdf", height=10, width=10)

lapply(fai[1:9, chr], function(i){
  pp1[chr == i, plot(PC2, PC1,  pch=20, las=1,bty='l', type="p", main=sub("chr", "", chr[1]), ylab="PC1",
            xlab="PC2")]
 })

dev.off()

pdf("chamomile_haplotype_separation_HiC_cc32_pc2_3.pdf", height=10, width=10)

lapply(fai[1:9, chr], function(i){
  pp1[chr == i, plot(PC3, PC2,  pch=20, las=1,bty='l', type="p", main=sub("chr", "", chr[1]), ylab="PC2",
            xlab="PC3")]
 })

dev.off()

pdf("chamomile_haplotype_separation_HiC_cc32_pc3_4.pdf", height=10, width=10)

lapply(fai[1:9, chr], function(i){
  pp1[chr == i, plot(PC4, PC3,  pch=20, las=1,bty='l', type="p", main=sub("chr", "", chr[1]), ylab="PC3",
            xlab="PC4")]
 })

dev.off()


#plot the PCA results: PC scores along the chromosomes
pdf("chamomile_haplotype_separation_HiC_cc32.pdf", height=4, width=10)
par(mfrow=c(2,5))
lapply(c("PC1", "PC2", "PC3", "PC4"), function(j){
 lapply(fai[1:9, chr], function(i){
  pp1[chr == i, plot(agp_pos/1e6, las=1, bty='l', get(j), type='n', main=sub("chr", "", chr[1]), ylab=j,
            xlab="Hv syntenic position [Mb]")]
  pp1[chr == i, lines(lwd=3, c(agp_pos/1e6, (agp_pos + contig_length)/1e6), c(get(j), get(j))), by=idx]
 })
 plot(axes=F, xlab="", ylab="", 0, type='n')
})
dev.off()


#####


##
# Import data and create assembly object
##

#create guide map
.libPaths(c("/home/mascher/tmp/3.5.1", "/filer-dg/agruppen/seq_shared/mascher/Rlib/3.5.1", "/opt/Bio/R_LIBS/3.5.1"))
source('/filer-dg/agruppen/dg2/cho/chamomile/phased_assembly/bitbucket/R/pseudomolecule_construction.R')

# Read single-copy regions generated in last step
f<-'/filer-dg/agruppen/dg2/cho/chamomile/diploid/phased_assembly/guide_map/pseudomolecules_v3_singlecopy_100bp.bed'
fread(f) -> p
setnames(p, c("agp_chr", "agp_start", "agp_end", "seq"))

# Change chromosome names to their respective chr numbers
# These are example names from maize

ids <- c("chr1",
         "chr2",
         "chr3",
         "chr4",
         "chr5",
         "chr6",
         "chr7",
         "chr8",
         "chr9")

for (i in 1:9){
  p[grepl(ids[i], agp_chr), agp_chr:=as.integer(i)]
}

# This is for unplaced scaffolds, assigning them as "NA"
p[grepl("NW", agp_chr), agp_chr:=as.integer(NA)]

p[, .(css_contig=seq, popseq_cM=(agp_start+agp_end)/2/1e6, #set pseudo-cM as Mb position
    sorted_genome = as.character(NA),
    css_contig_length = agp_end - agp_start,
    popseq_alphachr = as.character(agp_chr), # it needs to be character column
    sorted_arm = NA) ] -> pp

pp[,popseq_chr:=as.integer(popseq_alphachr)]
pp[, sorted_alphachr := popseq_alphachr]
pp[, sorted_chr := popseq_chr]
pp[, sorted_lib := popseq_alphachr]

# Save as RDS file for use in TRITEX pipeline
saveRDS(pp, "chamomile_v3_pseudopopseq.Rds")

#read MorexV3 gene-based guide map ("pseudo-POPSEQ")
readRDS('/filer-dg/agruppen/dg2/cho/chamomile/phased_assembly/guide_map/chamomile_pseudopopseq.Rds') -> pg

#read unitig lengths
f <- '/filer-dg/agruppen/dg2/cho/chamomile/diploid/phased_assembly/assembly/chamomile/chamomile2.hic.p_utg.fa.fai'
fread(f, head=F, select=1:2, col.names=c("scaffold", "length"))->fai

#read MorexV3 gene alignment (GMAP output) and merge with pseudo-POPSEQ table
fread('/filer-dg/agruppen/dg2/cho/chamomile/diploid/phased_assembly/assembly/chamomile/gmap_index_v3/chamomile2.hic.p_utg_pseudomolecule_v3_table.txt')->x
setnames(x, c("contig", "start", "end", "transcript", "alnlen", "id"))
fai[, .(scaffold, scaffold_length=length)][x[, .(css_contig=transcript, scaffold=contig, pos=start)], on="scaffold"] -> aln
pg[aln, on="css_contig"] -> aln

#read Hi-C map
dir <- '/filer-dg/agruppen/dg2/cho/chamomile/diploid/phased_assembly/mapping/hic'
fread(paste('find', dir, '| grep "hic_fragment_pairs.tsv.gz$" | xargs zcat'),
      header=F, col.names=c("scaffold1", "pos1", "scaffold2", "pos2"))->fpairs

#initialize assembly object, anchor to guide map and calculate physical coverage
init_assembly(fai=fai, cssaln=aln, fpairs=fpairs) -> assembly
anchor_scaffolds(assembly = assembly, popseq=pg, species="chamomile") -> assembly
add_hic_cov(assembly, binsize=1e4, binsize2=1e6, minNbin=50, innerDist=3e5, cores=20)->assembly

#save uncorrected assembly object
saveRDS(assembly, file="chamomile_assembly_v3.Rds")

##
# Check for and correct chimeras 
##

readRDS('chamomile_assembly.Rds') -> assembly

#make diagnostic plot for one contig

#plot all > 5e6 contig 
assembly$info[length >= 5e6, .(scaffold, length)][order(-length)] -> ss

i=1
ss[i]$scaffold -> s 
assembly$cov[s, on='scaffold'] -> b

for (ni in 2:length(ss$scaffold)) {
i=ni
ss[i]$scaffold -> s 
rbind(b, assembly$cov[s, on='scaffold'])->b
} 

plot_chimeras(assembly=assembly, scaffolds=b,  species="chamomile", refname="chamomile", autobreaks=F, mbscale=1, file="chamomile_assembly_all5e6_contig.pdf", cores=50)



#plot all > 1e6 contig & mri <= -2.5 
#not doing this step
assembly$info[length >= 1e6& mri <= -2.5, .(scaffold, length)][order(-length)] -> ss

i=1
ss[i]$scaffold -> s 
assembly$cov[s, on='scaffold'] -> b

for (ni in 2:length(ss$scaffold)) {
i=ni
ss[i]$scaffold -> s 
rbind(b, assembly$cov[s, on='scaffold'])->b
} 

plot_chimeras(assembly=assembly, scaffolds=b,  species="chamomile", refname="chamomile", autobreaks=F, mbscale=1,file="chamomile_assembly_all_contig.pdf", cores=50)



###############################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################


#set break point coordinates
#not break for chamomile
i=4
ss[i]$scaffold -> s 
assembly$cov[s, on='scaffold'][bin >= 10e6 & bin <= 12e6][order(r)][1, .(scaffold, bin)] -> b

i=44
ss[i]$scaffold -> s 
rbind(b, assembly$cov[s, on='scaffold'][bin >= 3e6 & bin <= 4e6][order(r)][1, .(scaffold, bin)])->b


setnames(b, "bin", "br")
plot_chimeras(assembly=assembly, scaffolds=b, br=b, species="chamomile", refname="chamomile",  mbscale=1,
         file="chamomile_assembly_chimeras_final.pdf", cores=30)

#implement the correction
break_scaffolds(b, assembly, prefix="contig_corrected_v1_", slop=1e4, cores=30, species="chamomile") -> assembly_v2

#save the object
saveRDS(assembly_v2, file="chamomile_assembly_v2.Rds")

#don't break chimera
assembly -> assembly_v2
saveRDS(assembly_v2, file="chamomile_assembly_v2.Rds")

##################################################################################################
##################################################################################################
##################################################################################################
#################################### Hap Split ###################################################
##################################################################################################
##################################################################################################
##################################################################################################


#read assembly object. This is the hifiasm unitig assembly after breaking chimeras based on coverage drops
readRDS('chamomile_assembly_v2.Rds') -> assembly_v2

##
# Positions contigs based on guide map
##

#read alignment of MorexV3 HC genes 
fread('/filer-dg/agruppen/dg2/cho/chamomile/diploid/phased_assembly/assembly/chamomile/gmap_index_v3/chamomile2.hic.p_utg_pseudomolecule_v3_table.txt')->x
setnames(x, c("contig", "start", "end", "transcript", "alnlen", "id"))

#read positions of MorexV3 genes
fread('/filer-dg/agruppen/dg2/cho/chamomile/diploid/phased_assembly/guide_map/pseudomolecules_v3_singlecopy_100bp.bed')->p
setnames(p, c("chr", "chr_start", "chr_end", "transcript"))

#keep only genes aligned for >= 80 % of their length with 90 % identity, allow up to two alignment (one per haplotype), modify for tetraploids
p[x, on="transcript"] -> px
px[alnlen >= 90 & id >= 97]->a
px[a[, .N, key=transcript][N <= 2]$transcript, on="transcript"]->aa

#convert original contig positions to corrected (assembly_v2) positions 
assembly_v2$info[, .(contig=orig_scaffold, start=orig_start, orig_start, scaffold)][aa, on=c("contig", "start"), roll=T] -> aa

#chromosome assignment
aa[, .N, key=.(chr, scaffold)][, p := N/sum(N), by=scaffold][order(-p)][!duplicated(scaffold)] -> cc

#keep contigs with at least 4 aligned genes, 75 % of aligned are from the major chromosome
cc[N >= 4 & p >= 0.75] -> cc

#check how much of the assembly can be assigned to chromosomes
assembly_v2$info[, .(scaffold, scaffold_length=length)][cc, on="scaffold"] -> cc

cc[scaffold_length>=300000]->cc 

sum(cc$scaffold_length)/1e6  ##[1] 4343.794



#get the approximate chromosome positions (median of alignment coordinates)
aa[cc[, .(scaffold, chr)], on=c("chr", "scaffold")][, .(pos=median(as.numeric(chr_start)), pos_mad=mad(as.numeric(chr_start))), key=scaffold][cc, on="scaffold"] -> cc
cc[, mr := pos_mad / scaffold_length]  
setorder(cc, chr, pos) 
cc[, chr_idx := 1:.N, by=chr]
cc[, agp_pos := c(0, cumsum(scaffold_length[-.N])), by=chr]

#save results
saveRDS(cc, file="chamomile_assembly_v2_Hv_guide.Rds")

##
# Positioning additional contigs by Hi-C
## 

#load guide map positions
readRDS('chamomile_assembly_v2_Hv_guide.Rds') -> cc
#get Hi-C links in the terminal 2 Mb of each scaffold
assembly_v2$fpairs[, .(scaffold=scaffold1, pos=pos1, link=scaffold2)] -> ff
assembly_v2$info[, .(scaffold, length)][ff, on="scaffold"] -> ff
ff[pos <= 2e6 | length - pos <= 2e6] -> ff
ff[scaffold != link] -> ff
ff[, .N, key=.(scaffold, link)] -> fa
#exclude scaffold-link paris with only a single Hi-C pair
fa[N > 1] -> fa
#add guide map positions
cc[, .(scaffold, scaffold_chr=chr, scaffold_pos=pos)][fa, on="scaffold"] -> fa
cc[, .(link=scaffold, link_chr=chr, link_pos=pos)][fa, on="link"] -> fa
#get approximate position based on Hi-C links
fa[!is.na(link_chr), .(n=sum(N), pos=weighted.mean(link_pos, N)), key=.(scaffold, scaffold_chr, scaffold_pos, link_chr)] -> fv
fv[, p := n/sum(n), by=scaffold]
assembly_v2$info[, .(scaffold, length)][fv, on="scaffold"] -> fv
#get chromosome assignment for scaffolds without position in guide map
fv[is.na(scaffold_chr)][order(-p)][!duplicated(scaffold)][order(-length)] -> cc_lift0

#exclude contigs shorter than 300 kb
cc_lift0[length >= 3e5] -> cc_lift
#merge guide map and Hi-C lift tables; write output
rbind(cc[, .(scaffold, chr, pos)], cc_lift[, .(scaffold, chr=link_chr, pos)]) -> a
a[assembly_v2$info[, .(scaffold, length)], on="scaffold"] -> a
saveRDS(a, file="chamomile_assembly_v2_Hv_guide+HiClift.Rds")

#check how much of the assembly is placed
 sum(a[!is.na(chr)]$length) 
#[1] 4572064073
 sum(a$length)
#[1] 5769542610
 4572064073/5769542610
#[1] 0.7924483

##
# Haplotype phasing
##

#read lengths of Morex chromosomes
fread('/filer-dg/agruppen/dg2/cho/chamomile/diploid/phased_assembly/guide_map/pseudomolecules_v3.fasta.fai', sel=1:2, col.names=c("chr", "len"))->chamofai

#read centromere position in Morex 
fread('/filer-dg/agruppen/dg2/cho/chamomile/diploid/phased_assembly/pseudomolecules_v3_centromere_positions.tsv')->cen
setnames(cen, c("chr", "cen_pos"))

#read scaffold posiitions and add them to the Hi-C link table
readRDS(file="chamomile_assembly_v2_Hv_guide.Rds") -> cc
assembly_v2$fpairs[, .N, key=.(ctg1=scaffold1, ctg2=scaffold2)] -> v
cc[, .(ctg1=scaffold, chr1=chr, pos1=pos)][v, on="ctg1"] -> v
cc[, .(ctg2=scaffold, chr2=chr, pos2=pos)][v, on="ctg2"] -> v

#run a PCA on the intra-chromosomal matrices
rbindlist(mclapply(mc.cores=7, chamofai[1:9, chr], function(j){
 #fill in empty scaffold pairs with 0
 setnames(v[chr1 == chr2 & chr1 == j][, chr2 := NULL], "chr1", "chr") -> x
 dcast(x, ctg1 ~ ctg2, value.var="N", fill = 0) -> x
 melt(x, id="ctg1", measure=setdiff(names(x), "ctg1"), variable.factor=F, value.name="N", variable.name="ctg2") -> x
 x[, l := log10(0.1 + N)]
 cc[, .(ctg1=scaffold, chr, pos1=pos)][x, on="ctg1"] -> x
 cc[, .(ctg2=scaffold, pos2=pos)][x, on="ctg2"] -> x
 chamofai[x, on="chr"] -> x
 cen[x, on="chr"] -> x
 x[, end_dist1 := pmin(pos1, len - pos1)]
 x[, end_dist2 := pmin(pos2, len - pos2)]
 x[, cen_dist1 := abs(cen_pos - pos1)]
 x[, cen_dist2 := abs(cen_pos - pos2)]
 x[, ldist := abs(pos1 - pos2)]
 #get "normalized" Hi-C ounts by removing the factor linear distance between loci, distance from centromere and distance from chromosome end (all log scaled)
 x[, res := lm(l ~ log(1+ldist) * log(cen_dist1) * log(cen_dist2) * log(end_dist1) * log(end_dist2))$res]

 #convert to matrix
 dcast(x, ctg1 ~ ctg2, value.var="res", fill=0) -> y
 y[, ctg1 := NULL]
 #run PCA on correlation matrix
 prcomp(cor(y), scale=T, center=T)->pca
 #get first four eigenvector
 data.table(contig=rownames(pca$rotation), pca$rotation[, 1:4]) -> p
 setorder(cc[, .(contig=scaffold, chr, pos)][p, on="contig"], chr, pos) -> pp
 pp[, chr := j]
})) -> pp
assembly_v2$info[, .(contig=scaffold, contig_length=length)][pp, on="contig"] -> pp
pp[, idx := 1:.N]
#save results
saveRDS(pp, file="chamomile_assembly_v2_haplotype_separation_HiC.Rds")

#read coverage information and convert to correct assembly coordinates
fread(cmd="grep '^S' /filer-dg/agruppen/dg2/cho/chamomile/diploid/phased_assembly/assembly/chamomile/chamomile2.hic.p_utg.gfa | cut -f 2,5 | tr ':' '\\t' | cut -f 1,4", head=F) -> cov
setnames(cov, c("contig", "cc"))
cov[, contig := sub('l$', '', sub("utg0*", "contig_", contig))]
assembly_v2$info[, .(contig=orig_scaffold, scaffold)][cov, on=c("contig")] -> cov
setorder(cc[, .(scaffold, scaffold_length, chr, pos)][cov, on="scaffold"], chr, pos) -> cov
cov[, idx := 1:.N]
saveRDS(file="chamomile_assembly_v2_cov.Rds", cov)

pp-> pq

#manually define cuts between haplotype for each chromosome
#2x contigs get "haplotype 3", i.e. present in both 1 and 2
cov[, .(contig=scaffold, cc)][pq, on="contig"] -> pq
pq[chr == "chr1" & PC1 > 0, hap := 1]
pq[chr == "chr1" & PC1 < 0, hap := 2]
pq[chr == "chr1" & cc >= 36, hap := 3]
pq[chr == "chr2" & PC1 > 0, hap := 1]
pq[chr == "chr2" & PC1 < 0, hap := 2]
pq[chr == "chr2" & cc >= 36, hap := 3]
pq[chr == "chr3" & PC1 > 0, hap := 1]
pq[chr == "chr3" & PC1 < 0, hap := 2]
pq[chr == "chr3" & cc >= 36, hap := 3]
pq[chr == "chr4" & PC1 > 0, hap := 1]
pq[chr == "chr4" & PC1 < 0, hap := 2]
pq[chr == "chr4" & cc >= 36, hap := 3]
pq[chr == "chr5" & PC1 > 0, hap := 1]
pq[chr == "chr5" & PC1 < 0, hap := 2]
pq[chr == "chr5" & cc >= 36, hap := 3]
pq[chr == "chr6" & PC1 > 0, hap := 1]
pq[chr == "chr6" & PC1 < 0, hap := 2]
pq[chr == "chr6" & cc >= 36, hap := 3]
pq[chr == "chr7" & PC1 > 0, hap := 1]
pq[chr == "chr7" & PC1 < 0, hap := 2]
pq[chr == "chr7" & cc >= 36, hap := 3]
pq[chr == "chr8" & PC1 > 0, hap := 1]
pq[chr == "chr8" & PC1 < 0, hap := 2]
pq[chr == "chr8" & cc >= 36, hap := 3]
pq[chr == "chr9" & PC1 > 0, hap := 1]
pq[chr == "chr9" & PC1 < 0, hap := 2]
pq[chr == "chr9" & cc >= 36, hap := 3]
pq[, col := "black"]
pq[hap == 1, col := "red"]
pq[hap == 2, col := "blue"]
pq[hap == 3, col := "purple"]


fread('/filer-dg/agruppen/dg2/cho/chamomile/diploid/phased_assembly/guide_map/pseudomolecules_v1.fa.fai', sel=1:2, col.names=c("chr", "len"))->chamofai
#plot results (PC1 score and coverage)
pdf("270623_assembly_v2_haplotype_separation_HiC.pdf", height=8, width=10)
par(mfrow=c(2, 1))
lapply(c("PC1"), function(j){
 lapply(chamofai[1:9, chr], function(i){
  pq[chr == i, plot(pos/1e6, las=1, bty='l', get(j), type='n', main=sub("chr", "", chr[1]), ylab="PC1",
        xlab="Hv syntenic position [Mb]", xlim=c(0, chamofai[i, len/1e6, on="chr"]))]
   pq[chr == i, lines(lwd=3, c(pos/1e6, (pos + contig_length)/1e6), c(PC1, PC1), col=col), by=idx]
  abline(v=c(0, chamofai[i, len/1e6, on="chr"]), col="blue")
  cov[chr == i, plot(pos/1e6, las=1, bty='l', cc, type='n', main=sub("chr", "", chr[1]), ylab="coverage",
        xlab="Hv syntenic position [Mb]", xlim=c(0, chamofai[i, len/1e6, on="chr"]))]
   cov[chr == i, lines(lwd=3, c(pos/1e6, (pos + scaffold_length)/1e6), c(cc, cc), col=1), by=idx]
 })
})

dev.off()

saveRDS(pq, file="chamomile_assembly_v2_haplotype_separation_v0.Rds")


##
# Assign contigs that are not place in the guide map to haplotype using Hi-C
##

#read haplotype separation and guide map table
readRDS(file="chamomile_assembly_v2_haplotype_separation_v0.Rds") -> pq
readRDS(file="chamomile_assembly_v2_Hv_guide+HiClift.Rds") -> cc

readRDS('chamomile_assembly_v2.Rds') -> assembly_v2
readRDS('chamomile_assembly_v2_cov.Rds') -> cov


#find and tabulate links between unplaced and placed con tigs
assembly_v2$fpairs[, .(scaffold=scaffold1, pos=pos1, link=scaffold2)] -> ff
assembly_v2$info[, .(scaffold, length)][ff, on="scaffold"] -> ff
ff[pos <= 2e6 | length - pos <= 2e6] -> ff
ff[scaffold != link] -> ff
ff[, .N, key=.(scaffold, link)] -> fa
fa[N > 1] -> fa
fa[scaffold %in% setdiff(cc$scaffold, pq$contig)] -> fa
fa[link %in% pq[hap %in% 1:2]$contig] -> fa
cc[, .(scaffold, scaffold_chr=chr, scaffold_pos=pos)][fa, on="scaffold"] -> fa
cc[, .(link=scaffold, link_chr=chr, link_pos=pos)][fa, on="link"] -> fa
pq[, .(link=contig, hap)][fa, on="link"] -> fa
fa[scaffold_chr == link_chr] -> fa
fa[, .(n=sum(N)), key=.(scaffold, chr=scaffold_chr, hap)] -> fv
fv[, p := n/sum(n), by=scaffold]
assembly_v2$info[, .(scaffold, length)][fv, on="scaffold"] -> fv
fv[length >= 0][order(-p)][!duplicated(scaffold)][order(-length)] -> hap_lift0
cov[, .(scaffold, cc)][hap_lift0, on="scaffold"] -> hap_lift0

saveRDS(hap_lift0, file="chamomile_assembly_v2_hap_lift0.Rds") 

#keep contigs unassigned to one haplotype OR contigs assigned to both haplotyped with double coverage
hap_lift0[p >= 0.8 | (p <= 0.55 & cc >= 36)] -> hap_lift

#keep onlt contigs >= 300 kb
hap_lift[length >= 3e5] -> hap_lift

#contigs assigned to both haplotype, get haplotype "3"
hap_lift[p <= 0.55, hap := 3]

#combine both haplotype tables
rbind(pq[, .(scaffold=contig, hap)],  hap_lift[, .(scaffold, hap)]) -> hh
hh[cc, on="scaffold"] -> hh
saveRDS(hh, file="chamomile_assembly_v2_Hv_guide+HiClift+hap.Rds") 


sum(hh[!is.na(chr)&hap==1]$length)
#[1] 2189264147
sum(hh[!is.na(chr)&hap==2]$length)
#[1] 2273679966
sum(hh[!is.na(chr)&hap==3]$length)
#[1] 80993772


##
# Construct first Hi-C map
##

#read Hi-C fragment data
f <- '/filer-dg/agruppen/dg2/cho/chamomile/diploid/phased_assembly/assembly/chamomile/chamomile2.hic.p_utg_DpnII_fragments_30bp.bed'
read_fragdata(info=assembly_v2$info, file=f)->frag_data

#read haplotype information
readRDS(file="chamomile_assembly_v2_Hv_guide+HiClift+hap.Rds") -> hh

#Hi-C map for haplotype 1
hh[hap %in% c(1,3), .(scaffold, length, chr=as.integer(substr(chr, 4, 4)), cM = pos/1e6)] -> hic_info
frag_data$info[, .(scaffold, nfrag)][hic_info, on="scaffold"] -> hic_info

hic_map(info=hic_info, assembly=assembly_v2, frags=frag_data$bed, species="chamomile", ncores=21,
  min_nfrag_scaffold=30, max_cM_dist = 1000,
  binsize=2e5, min_nfrag_bin=10, gap_size=100)->hic_map_v1_hap1

saveRDS(hic_map_v1_hap1, file="chamomile_hic_map_v1_hap1.Rds")

#Hi-C map for haplotype 2
hh[hap %in% c(2,3), .(scaffold, length, chr=as.integer(substr(chr, 4, 4)), cM = pos/1e6)] -> hic_info
frag_data$info[, .(scaffold, nfrag)][hic_info, on="scaffold"] -> hic_info

hic_map(info=hic_info, assembly=assembly_v2, frags=frag_data$bed, species="chamomile", ncores=21,
  min_nfrag_scaffold=30, max_cM_dist = 1000,
  binsize=2e5, min_nfrag_bin=10, gap_size=100)->hic_map_v1_hap2

saveRDS(hic_map_v1_hap2, file="chamomile_hic_map_v1_hap2.Rds")

#create Hi-C plots for haplotype 1
snuc <- '/filer-dg/agruppen/dg2/cho/chamomile/phased_assembly/assembly/chamomile/chamomile2.hic.p_utg_DpnII_fragments_30bp_split.nuc.txt'
hic_plots(rds="chamomile_hic_map_v1_hap1.Rds", cov=F, assembly=assembly_v2, cores=30, species="chamomile", nuc=snuc) -> hic_map_v1_hap1_l

#create Hi-C plots for haplotype 2
snuc <- '/filer-dg/agruppen/dg2/cho/chamomile/phased_assembly/assembly/chamomile/chamomile2.hic.p_utg_DpnII_fragments_30bp_split.nuc.txt'
hic_plots(rds="chamomile_hic_map_v1_hap2.Rds", cov=F, assembly=assembly_v2, cores=30, species="chamomile", nuc=snuc) -> hic_map_v1_hap2_l

#combine both Hi-C maps
combine_hic(hic_map_v1_hap1, hic_map_v1_hap2, assembly,species="chamomile_2X") -> ch
ch$hic_map -> hic_map_v1
ch$assembly -> assembly_v2_hap

#create inter-chromosomal plot for combined hap1+hap2 Hi-C haplotypes
nuc <- '/filer-dg/agruppen/dg2/cho/chamomile/diploid/phased_assembly/assembly/chamomile/chamomile2.hic.p_utg_DpnII_fragments_30bp_split.nuc.txt'
add_psmol_fpairs(assembly=assembly_v2_hap, hic_map=hic_map_v1, nucfile=nuc)->hic_map_v1_l

bin_hic_step(hic=hic_map_v1_l$links, frags=hic_map_v1_l$frags, binsize=1e6, chrlen=hic_map_v1_l$chrlen, cores=14)->hic_map_v1_l$hic_1Mb

f <- "chamomile_hic_map_v1_interchromosomal.png"
interchromosomal_matrix_plot(hic_map=hic_map_v1_l, file=f, species="chamomile_2X")

#write Hi-C map for editing
write_hic_map(rds="chamomile_hic_map_v1_hap1.Rds", file="chamomile_hic_map_v1_hap1.xlsx", species="chamomile")
write_hic_map(rds="chamomile_hic_map_v1_hap2.Rds", file="chamomile_hic_map_v1_hap2.xlsx", species="chamomile")

###### write fasta

readRDS('chamomile_assembly_v3.Rds') -> assembly_v3
readRDS('chamomile_hic_map_v6_hap1.Rds') -> hic_map_v6_hap1

readRDS('chamomile_hic_map_v5_hap2.Rds') -> hic_map_v5_hap2


fasta <- '/filer-dg/agruppen/dg2/cho/chamomile/diploid/phased_assembly/assembly/chamomile/chamomile2.hic.p_utg.fa'
sink("chamomile_pseudomolecules_v6_hap1.log")
compile_psmol(fasta=fasta, output="chamomile_pseudomolecules_v6_hap1", hic_map=hic_map_v6_hap1, assembly=assembly_v3, cores=30)
sink()

sink("chamomile_pseudomolecules_v5_hap2.log")
compile_psmol(fasta=fasta, output="chamomile_pseudomolecules_v5_hap2", hic_map=hic_map_v5_hap2, assembly=assembly_v3, cores=30)
sink()


####run alginment for 2 haplotypes 
/opt/Bio/minimap2/2.17/bin/minimap2 -t 20 -f 0.005 -2 -I 20G -K5G -x asm5 /filer-dg/agruppen/dg2/cho/chamomile/phased_assembly/phasing_assembly/chamomile_pseudomolecules_v2_hap1/chamomile_pseudomolecules_v2_hap1.fasta /filer-dg/agruppen/dg2/cho/chamomile/phased_assembly/phasing_assembly/chamomile_pseudomolecules_v2_hap2/chamomile_pseudomolecules_v2_hap2.fasta > chamomile_compare_hap_v2.paf 2> chamomile_compare_hap_v2.err &


####hap1


.libPaths(c("/home/mascher/tmp/3.5.1", "/filer-dg/agruppen/seq_shared/mascher/Rlib/3.5.1", "/opt/Bio/R_LIBS/3.5.1"))
source("/filer-dg/agruppen/seq_shared/mascher/code_repositories/triticeae.bitbucket.io/R/pseudomolecule_construction.R")

read_paf("/filer-dg/agruppen/dg2/cho/chamomile/diploid/phased_assembly/phasing_assembly/compare_hap/compare_pseudomolecule_v3_ref_hap1_v6.paf.gz", save=F, primary_only=T) -> paf #### A_B qry=B ref=A 
#paf[, reference := sub("_1", "", reference)]
paf[query == reference, ] -> paf
paf[mapq >= 30] -> paf
paf[alnlen >= 2000] -> paf 
paf[, chr := sub("chr", "", query)]
#paf[, chr := sub("H", "", chr)]
paf[, idx := 1:.N]

pdf("hap1_chamomile_v6_correlation_plots.pdf") 
par(mar=c(5,5,3,1))
par(cex.main=1)
par(cex.lab=1)
par(cex.axis=1)
lapply(1:9, function(i){
 paf[chr == i] -> pafL
 pafL[chr == i, plot(las=1, bty='l',type='n', 0, cex.lab=1.2, xlab="chamomile hap1 CCS (Mb)", ylab="chamomile integrated CCS (Mb)", xlim=c(0,max(query_end))/1e6, ylim=c(0,max(reference_end))/1e6, main=paste0(i, "H"))] ### xlab = qry = B and ylab = ref = A  
 pafL[orientation == 1, lines(c(query_start/1e6, query_end/1e6), c(reference_start/1e6, reference_end/1e6), col="#000000", lwd=2), by=idx]
 pafL[orientation == -1, lines(c(query_start/1e6, query_end/1e6), c(reference_end/1e6, reference_start/1e6), col="#000000", lwd=2), by=idx]
})
dev.off()



####hap2


.libPaths(c("/home/mascher/tmp/3.5.1", "/filer-dg/agruppen/seq_shared/mascher/Rlib/3.5.1", "/opt/Bio/R_LIBS/3.5.1"))
source("/filer-dg/agruppen/seq_shared/mascher/code_repositories/triticeae.bitbucket.io/R/pseudomolecule_construction.R")

read_paf("/filer-dg/agruppen/dg2/cho/chamomile/diploid/phased_assembly/phasing_assembly/compare_hap/compare_pseudomolecule_v3_ref_hap2_v5.paf.gz", save=F, primary_only=T) -> paf #### A_B qry=B ref=A 
#paf[, reference := sub("_1", "", reference)]
paf[query == reference, ] -> paf
paf[mapq >= 30] -> paf
paf[alnlen >= 2000] -> paf 
paf[, chr := sub("chr", "", query)]
#paf[, chr := sub("H", "", chr)]
paf[, idx := 1:.N]

pdf("hap2_chamomile_v5_correlation_plots.pdf") 
par(mar=c(5,5,3,1))
par(cex.main=1)
par(cex.lab=1)
par(cex.axis=1)
lapply(1:9, function(i){
 paf[chr == i] -> pafL
 pafL[chr == i, plot(las=1, bty='l',type='n', 0, cex.lab=1.2, xlab="chamomile hap2 CCS (Mb)", ylab="chamomile integrated CCS (Mb)", xlim=c(0,max(query_end))/1e6, ylim=c(0,max(reference_end))/1e6, main=paste0(i, "H"))] ### xlab = qry = B and ylab = ref = A  
 pafL[orientation == 1, lines(c(query_start/1e6, query_end/1e6), c(reference_start/1e6, reference_end/1e6), col="#000000", lwd=2), by=idx]
 pafL[orientation == -1, lines(c(query_start/1e6, query_end/1e6), c(reference_end/1e6, reference_start/1e6), col="#000000", lwd=2), by=idx]
})
dev.off()

#####hap1 and hap2


read_paf("/filer-dg/agruppen/dg2/cho/chamomile/diploid/phased_assembly/phasing_assembly/compare_hap/compare_hap1_v6_hap2_v5.paf.gz", save=F, primary_only=T) -> paf #### A_B qry=B ref=A 

paf[query == reference, ] -> paf
paf[mapq >= 30] -> paf
paf[alnlen >= 2000] -> paf 
paf[, chr := sub("chr", "", query)]
#paf[, chr := sub("H", "", chr)]
paf[, idx := 1:.N]

pdf("hap1_v6_hap2_v5_chamomile_correlation_plots.pdf") 
par(mar=c(5,5,3,1))
par(cex.main=1)
par(cex.lab=1)
par(cex.axis=1)
lapply(1:9, function(i){
 paf[chr == i] -> pafL
 pafL[chr == i, plot(las=1, bty='l',type='n', 0, cex.lab=1.2, xlab="chamomile_hap2 CCS (Mb)", ylab="chamomile_hap1 CCS (Mb)", xlim=c(0,max(query_end))/1e6, ylim=c(0,max(reference_end))/1e6, main=paste0(i, "H"))] ### xlab = qry = B and ylab = ref = A  
 pafL[orientation == 1, lines(c(query_start/1e6, query_end/1e6), c(reference_start/1e6, reference_end/1e6), col="#000000", lwd=2), by=idx]
 pafL[orientation == -1, lines(c(query_start/1e6, query_end/1e6), c(reference_end/1e6, reference_start/1e6), col="#000000", lwd=2), by=idx]
})
dev.off()


#import Hi-C maps after editing, v1 -> v2
read_hic_map(rds="chamomile_hic_map_v1_hap1.Rds", file="chamomile_hic_map_v1_hap1_edit6.xlsx") -> nmap
hic_map(species="chamomile", agp_only=T, map=nmap)->hic_map_v6_hap1
saveRDS(hic_map_v6_hap1, file="chamomile_hic_map_v6_hap1.Rds")

read_hic_map(rds="chamomile_hic_map_v1_hap2.Rds", file="chamomile_hic_map_v1_hap2_edit5.xlsx") -> nmap
hic_map(species="chamomile", agp_only=T, map=nmap)->hic_map_v5_hap2
saveRDS(hic_map_v5_hap2, file="chamomile_hic_map_v5_hap2.Rds")

#proceed with Hi-C plots and repeat edit cycle if need be



#create Hi-C plots for haplotype 1
snuc <- '/filer-dg/agruppen/dg2/cho/chamomile/diploid/phased_assembly/assembly/chamomile/chamomile2.hic.p_utg_DpnII_fragments_30bp_split.nuc.txt'
hic_plots(rds="chamomile_hic_map_v6_hap1.Rds", cov=F, assembly=assembly_v2, cores=30, species="chamomile", nuc=snuc) -> hic_map_v6_hap1_l
saveRDS(hic_map_v6_hap1_l, file="hic_map_v6_hap1_l.Rds")

#create Hi-C plots for haplotype 2
snuc <- '/filer-dg/agruppen/dg2/cho/chamomile/diploid/phased_assembly/assembly/chamomile/chamomile2.hic.p_utg_DpnII_fragments_30bp_split.nuc.txt'
hic_plots(rds="chamomile_hic_map_v5_hap2.Rds", cov=F, assembly=assembly_v2, cores=30, species="chamomile", nuc=snuc) -> hic_map_v5_hap2_l
saveRDS(hic_map_v5_hap2_l, file="hic_map_v5_hap2_l.Rds")



#readRDS('FB20_005_1_assembly_v2.Rds') -> assembly_v2
readRDS('chamomile_hic_map_v6_hap1.Rds') -> hic_map_v6_hap1
readRDS('chamomile_hic_map_v4_hap2.Rds') -> hic_map_v4_hap2


fasta <- '/filer-dg/agruppen/dg2/cho/chamomile/diploid/phased_assembly/assembly/chamomile/chamomile2.hic.p_utg.fa'
sink("chamomile_pseudomolecules_v6_hap1.log")
compile_psmol(fasta=fasta, output="chamomile_pseudomolecules_v6_hap1", hic_map=hic_map_v6_hap1, assembly=assembly_v2, cores=30)
sink()

sink("chamomile_pseudomolecules_v5_hap2.log")
compile_psmol(fasta=fasta, output="chamomile_pseudomolecules_v5_hap2", hic_map=hic_map_v5_hap2, assembly=assembly_v2, cores=30)
sink()


combine_hic(hic_map_v6_hap1, hic_map_v5_hap2, assembly_v2) -> ch2
ch2$hic_map -> hic_map_v3
ch2$assembly -> assembly_v3_hap

nuc <- '/filer-dg/agruppen/dg2/cho/chamomile/diploid/phased_assembly/assembly/chamomile/chamomile2.hic.p_utg_DpnII_fragments_30bp_split.nuc.txt'
add_psmol_fpairs(assembly=assembly_v3_hap, hic_map=hic_map_v3, nucfile=nuc)->hic_map_v3_l

bin_hic_step(hic=hic_map_v3_l$links, frags=hic_map_v3_l$frags, binsize=1e6, chrlen=hic_map_v3_l$chrlen, cores=14)->hic_map_v3_l$hic_1Mb

f <- "chamomile_hic_map_v3_interchromosomal.png"
interchromosomal_matrix_plot(hic_map=hic_map_v3_l, file=f, species="chamomile_2X")
saveRDS(hic_map_v3_l, file="hic_map_v3_l.Rds")




####hap1


.libPaths(c("/home/mascher/tmp/3.5.1", "/filer-dg/agruppen/seq_shared/mascher/Rlib/3.5.1", "/opt/Bio/R_LIBS/3.5.1"))
source('/filer-dg/agruppen/dg2/cho/chamomile/diploid/phased_assembly/bitbucket/R/pseudomolecule_construction.R')

read_paf("hap1_chamomile_v2_hap1.paf.gz", save=F, primary_only=T) -> paf #### A_B qry=B ref=A 
paf[, reference := sub("_1", "", reference)]
paf[query == reference, ] -> paf
paf[mapq >= 30] -> paf
paf[alnlen >= 2000] -> paf 
paf[, chr := sub("chr", "", query)]
paf[, chr := sub("H", "", chr)]
paf[, idx := 1:.N]

pdf("hap1_FB20_005_1_v2_hap1_correlation_plots.pdf") 
par(mar=c(5,5,3,1))
par(cex.main=1)
par(cex.lab=1)
par(cex.axis=1)
lapply(1:7, function(i){
 paf[chr == i] -> pafL
 pafL[chr == i, plot(las=1, bty='l',type='n', 0, xlab="FB20_005_1 hap1 CCS (Mb)", ylab="FB19 Hap1 CCS (Mb)", xlim=c(0,max(query_end))/1e6, ylim=c(0,max(reference_end))/1e6, main=paste0(i, "H"))] ### xlab = qry = B and ylab = ref = A  
 pafL[orientation == 1, lines(c(query_start/1e6, query_end/1e6), c(reference_start/1e6, reference_end/1e6), col="#000000", lwd=2), by=idx]
 pafL[orientation == -1, lines(c(query_start/1e6, query_end/1e6), c(reference_end/1e6, reference_start/1e6), col="#000000", lwd=2), by=idx]
})
dev.off()



####hap2


.libPaths(c("/home/mascher/tmp/3.5.1", "/filer-dg/agruppen/seq_shared/mascher/Rlib/3.5.1", "/opt/Bio/R_LIBS/3.5.1"))
source("/filer-dg/agruppen/seq_shared/mascher/code_repositories/triticeae.bitbucket.io/R/pseudomolecule_construction.R")

read_paf("compare_hap/compare_pseudomolecule_v3_ref_hap2.paf.gz", save=F, primary_only=T) -> paf #### A_B qry=B ref=A 
paf[, reference := sub("_1", "", reference)]
paf[query == reference, ] -> paf
paf[mapq >= 30] -> paf
paf[alnlen >= 2000] -> paf 
paf[, chr := sub("chr", "", query)]
paf[, chr := sub("H", "", chr)]
paf[, idx := 1:.N]

pdf("hap1_FB20_005_1_v2_hap2_correlation_plots.pdf") 
par(mar=c(5,5,3,1))
par(cex.main=1)
par(cex.lab=1)
par(cex.axis=1)
lapply(1:9, function(i){
 paf[chr == i] -> pafL
 pafL[chr == i, plot(las=1, bty='l',type='n', 0, xlab="FB20_005_1 hap2 CCS (Mb)", ylab="FB19 Hap1 CCS (Mb)", xlim=c(0,max(query_end))/1e6, ylim=c(0,max(reference_end))/1e6, main=paste0(i, "H"))] ### xlab = qry = B and ylab = ref = A  
 pafL[orientation == 1, lines(c(query_start/1e6, query_end/1e6), c(reference_start/1e6, reference_end/1e6), col="#000000", lwd=2), by=idx]
 pafL[orientation == -1, lines(c(query_start/1e6, query_end/1e6), c(reference_end/1e6, reference_start/1e6), col="#000000", lwd=2), by=idx]
})
dev.off()




source('/filer-dg/agruppen/seq_shared/mascher/code_repositories/triticeae.bitbucket.io/R/pseudomolecule_construction.R')

readRDS('FB20_005_1_hic_map_v2_hap1.Rds') -> hic_map_v4_hap1
readRDS('FB20_005_1_hic_map_v2_hap2.Rds') -> hic_map_v2_hap2



hic_map_v4_hap1$agp[!is.na(chr)]$scaffold->h1
hic_map_v2_hap2$agp[!is.na(chr)]$scaffold->h2



hic_map_v4_hap1$agp[,.(scaffold,scaffold_length)][scaffold!="gap"] -> chrh

chrh[,hap:=0]
chrh[scaffold %in% h1, hap:=hap+10]
chrh[scaffold %in% h2, hap:=hap+2]


write.table(chrh,"FB20_005_1_contig_all.txt",quote=F,row.names=F,sep="\t")

sum(hic_map_v4_hap1$agp[!is.na(chr)]$scaffold_length)
sum(hic_map_v2_hap2$agp[!is.na(chr)]$scaffold_length)


sum(chrh[hap==0]$scaffold_length)


mkdir uncontig
cd uncontig

less ../FB20_005_1_contig_all.txt | awk '{if($3==0) print $1}' > uncontig

seqkit grep -f uncontig ../FB20_005_1_pseudomolecules_v2_hap1/*_assembly_v2.fasta > FB20_005_1_unanchor.fasta








###########################################################

source('/filer-dg/agruppen/dg2/cho/chamomile/phased_assembly/bitbucket/R/pseudomolecule_construction.R')


chrNames <- function(agp=F, species="wheat") {
 if(species == "wheat"){
  data.table(alphachr=apply(expand.grid(1:7, c("A", "B", "D"), stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:21)->z
 } else if (species == "barley"){
  data.table(alphachr=apply(expand.grid(1:7, "H", stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:7)->z
 }
 else if (species == "rye"){
  data.table(alphachr=apply(expand.grid(1:7, "R", stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:7)->z
 }
 else if (species == "lolium"){
  data.table(alphachr=as.character(1:7), chr=1:7)->z
 }
 else if (species == "maize"){
  data.table(alphachr=as.character(1:10), chr=1:10)->z
 }
 else if (species == "sharonensis"){
  data.table(alphachr=apply(expand.grid(1:7, "S", stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:7)->z
 }
 else if (species == "oats"){
  data.table(alphachr=sub(" ", "", apply(expand.grid(1:21, "M", stringsAsFactors=F), 1, function(x) paste(x, collapse=""))), chr=1:21)->z
 }
 else if (species == "oats_new"){
  data.table(alphachr=apply(expand.grid(1:7, c("A", "C", "D"), stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:21)->z
 }
 else if (species == "avena_barbata"){
  data.table(alphachr=apply(expand.grid(1:7, c("A", "B"), stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:14)->z
 }
 else if (species == "hordeum_bulbosum"){
  data.table(alphachr=sort(apply(expand.grid(1:7, c("H_"), c(1:2), stringsAsFactors=F), 1, function(x) paste(x, collapse=""))), chr=1:14)->z
 }
 else if (species == "faba_bean"){
  data.table(alphachr=apply(expand.grid(1:6, "", stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:6)->z
 }
 else if (species == "chamomile"){
  data.table(alphachr=apply(expand.grid(1:9, stringsAsFactors=F), 1,function(x) paste(x, collapse="")), chr=1:9)->z
 }
 else if (species == "chamomile_2X"){
  data.table(alphachr=sort(apply(expand.grid(1:9, c("_"), c(1:2), stringsAsFactors=F), 1, function(x) paste(x, collapse=""))), chr=1:18)->z
 }
 else if (species == "Pginseng"){
  data.table(alphachr=apply(expand.grid(1:24, "", stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:24)->z
 }
 if(agp){
  rbind(z, data.table(alphachr="Un", chr=0))[, agp_chr := paste0("chr", alphachr)]->z
 }
 z[]
}

wheatchr<-function(agp=F, species="wheat") {
 if(species == "wheat"){
  data.table(alphachr=apply(expand.grid(1:7, c("A", "B", "D"), stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:21)->z
 } else if (species == "barley"){
  data.table(alphachr=apply(expand.grid(1:7, "H", stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:7)->z
 }
 else if (species == "rye"){
  data.table(alphachr=apply(expand.grid(1:7, "R", stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:7)->z
 }
 else if (species == "lolium"){
  data.table(alphachr=as.character(1:7), chr=1:7)->z
 }
 else if (species == "maize"){
  data.table(alphachr=as.character(1:10), chr=1:10)->z
 }
 else if (species == "sharonensis"){
  data.table(alphachr=apply(expand.grid(1:7, "S", stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:7)->z
 }
 else if (species == "oats"){
  data.table(alphachr=sub(" ", "", apply(expand.grid(1:21, "M", stringsAsFactors=F), 1, function(x) paste(x, collapse=""))), chr=1:21)->z
 }
 else if (species == "oats_new"){
  data.table(alphachr=apply(expand.grid(1:7, c("A", "C", "D"), stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:21)->z
 }
 else if (species == "avena_barbata"){
  data.table(alphachr=apply(expand.grid(1:7, c("A", "B"), stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:14)->z
 }
 else if (species == "hordeum_bulbosum"){
  data.table(alphachr=sort(apply(expand.grid(1:7, c("H_"), c(1:2), stringsAsFactors=F), 1, function(x) paste(x, collapse=""))), chr=1:14)->z
 }
 else if (species == "faba_bean"){
  data.table(alphachr=apply(expand.grid(1:6, "", stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:6)->z
 }
 else if (species == "chamomile"){
  data.table(alphachr=apply(expand.grid(1:9, stringsAsFactors=F), 1,function(x) paste(x, collapse="")), chr=1:9)->z
 }
 else if (species == "chamomile_2X"){
  data.table(alphachr=sort(apply(expand.grid(1:9, c("_"), c(1:2), stringsAsFactors=F), 1, function(x) paste(x, collapse=""))), chr=1:18)->z
 }
 else if (species == "Pginseng"){
  data.table(alphachr=apply(expand.grid(1:24, "", stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:24)->z
 }
 if(agp){
  rbind(z, data.table(alphachr="Un", chr=0))[, agp_chr := paste0("chr", alphachr)]->z
 }
 z[]
}

#Function to combine the Hi-C maps from two haplotypes
combine_hic <- function(hap1, hap2, assembly, species="chamomile_2X"){
 assembly_v2 <- assembly
 hic_map_v1_hap1 <- hap1
 hic_map_v1_hap2 <- hap2
 hic_map_v1_hap1$agp[agp_chr != "chrUn"] -> a1
 hic_map_v1_hap2$agp[agp_chr != "chrUn"] -> a2
 a1[, agp_chr := paste0(agp_chr, "_1")]
 a2[, agp_chr := paste0(agp_chr, "_2")]
 a1[, chr := NULL]
 a2[, chr := NULL]
 chrNames(agp=T, species)[, .(chr, agp_chr)][a1, on="agp_chr"] -> a1
 chrNames(agp=T, species)[, .(chr, agp_chr)][a2, on="agp_chr"] -> a2

 c(a1[scaffold != "gap"]$scaffold, a2[scaffold != "gap"]$scaffold) -> s
 s[duplicated(s)] -> s

 a1[s, on="scaffold", scaffold := paste0(scaffold, "_hap1")]
 a2[s, on="scaffold", scaffold := paste0(scaffold, "_hap2")]
 rbind(a1, a2) -> a
 
 hic_map_v1_hap1$chrlen[!is.na(chr)][, .(agp_chr=paste0(agp_chr, "_1"), length, truechr)] -> l1
 hic_map_v1_hap2$chrlen[!is.na(chr)][, .(agp_chr=paste0(agp_chr, "_2"), length, truechr)] -> l2
 rbind(l1, l2) -> l
 l[, offset := cumsum(c(0, length[1:(.N-1)]))]
 l[, plot_offset := cumsum(c(0, length[1:(.N-1)]+1e8))]
 chrNames(agp=T, species)[l, on="agp_chr"] -> l

 copy(assembly_v2$info) -> ai
 ai[!s, on="scaffold"] -> u
 ai[s, on="scaffold"][, scaffold := paste0(scaffold, "_hap1")] -> i1
 ai[s, on="scaffold"][, scaffold := paste0(scaffold, "_hap2")] -> i2
 rbind(u, i1, i2) -> i

 assembly_v2$fpairs[, .(scaffold1, scaffold2, pos1, pos2)] -> f
 f[scaffold1 %in% s, scaffold1 := paste0(scaffold1, ifelse(runif(.N) > 0.5, "_hap1", "_hap2"))]
 f[scaffold2 %in% s, scaffold2 := paste0(scaffold2, ifelse(runif(.N) > 0.5, "_hap1", "_hap2"))]

 list(info=i, fpairs=f) -> assembly_hap
 list(agp=a, chrlen=l) -> hic_map
 list(assembly_hap=assembly_hap, hic_map=hic_map)
}




