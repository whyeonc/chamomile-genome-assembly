chamomile

/opt/Bio/R/3.5.1/bin/R
.libPaths(c("/home/mascher/tmp/3.5.1", "/filer-dg/agruppen/seq_shared/mascher/Rlib/3.5.1", "/opt/Bio/R_LIBS/3.5.1"))
source("/filer-dg/agruppen/seq_shared/mascher/code_repositories/triticeae.bitbucket.io/R/pseudomolecule_construction.R")

f<-'reference_assembly_singlecopy_100bp.bed'
fread(f) -> p
setnames(p, c("agp_chr", "agp_start", "agp_end", "seq"))
ids <- c("OU343217.1",
         "OU343218.1",
         "OU343219.1",
         "OU343220.1",
         "OU343221.1",
         "OU343222.1",
         "OU343223.1",
         "OU343224.1",
         "OU343225.1")
for (i in 1:9){
  p[grepl(ids[i], agp_chr), agp_chr:=as.integer(i)]
}

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
saveRDS(pp, "chamomile2_pseudopopseq.Rds")

source("/filer-dg/agruppen/dg2/cho/chamomile/diploid/hifi+hic/hap1_hic/bitbucket/R/pseudomolecule_construction.R")
readRDS('/filer-dg/agruppen/dg2/cho/chamomile/hifi+hic/hap1_hic/chamomile2_hap1_pseudopopseq.Rds')-> popseq
f<- '/filer-dg/agruppen/dg2/cho/chamomile/diploid/integrated_/hap1_hic/assembly/chamomile2/chamomile2.hic.hap1.p_ctg.fa.fai'
fread(f, head=F, select=1:2, col.names=c("scaffold", "length")) -> fai
f <- '/filer-dg/agruppen/dg2/cho/chamomile/hifi+hic/hap1_hic/guide_map/chamomile_hap1.paf.gz'
read_morexaln_minimap(paf=f, popseq=popseq, minqual=30, minlen=500, prefix=F) -> morexaln
dir <- '/filer-dg/agruppen/dg2/cho/chamomile/hifi+hic/hap1_hic/mapping/hic'
fread(paste('find', dir, '| grep "hic_fragment_pairs.tsv.gz$" | xargs zcat'),
       header=F, col.names=c("scaffold1", "pos1", "scaffold2", "pos2")) -> fpairs
init_assembly(fai=fai, cssaln=morexaln, fpairs=fpairs) -> assembly
anchor_scaffolds(assembly = assembly, popseq=popseq, species="chamomile") -> assembly
add_hic_cov(assembly, binsize=1e4, binsize2=1e6, minNbin=50, innerDist=3e5, cores=40) -> assembly
saveRDS(assembly, file="assembly_re.Rds")

assembly$info[length >= 1e6, .(scaffold, length)][order(-length)] -> s
plot_chimeras(assembly=assembly, scaffolds=s, species="chamomile", refname="makinoi", autobreaks=F, mbscale=1, file="assembly_1Mb.pdf", cores=40)

assembly$info[length >= 1e6 & mri <= -3, .(scaffold, length)][order(-length)] -> s

plot_chimeras(assembly=assembly, scaffolds=s, species="chamomile", refname="makinoi", autobreaks=F, mbscale=1, file="assembly_chimeras.pdf", cores=30)


assembly-> assembly_v2

saveRDS(assembly_v2, file="assembly_v2.Rds")

fbed <- '/filer-dg/agruppen/dg2/cho/chamomile/diploid/hifi+hic/hap1_hic/assembly/chamomile2/chamomile2.hic.hap1.p_ctg_DpnII_fragments_30bp.bed'
read_fragdata(info=assembly_v2$info, file=fbed) -> frag_data

# Consider only contigs >= 1 Mb first
frag_data$info[!is.na(hic_chr) & length >= 1e6, .(scaffold, nfrag, length, chr=hic_chr, cM=popseq_cM)] -> hic_info

hic_map(info=hic_info, assembly=assembly_v2, frags=frag_data$bed,
	species="chamomile", ncores=30, min_nfrag_scaffold=30, 
	max_cM_dist = 1000, binsize=2e5, min_nfrag_bin=10, gap_size=100) -> hic_map_v1

saveRDS(hic_map_v1, file="hic_map_v1.Rds")

# Make the Hi-C plots, files will be placed in your working directory.
snuc <- '/filer-dg/agruppen/dg2/cho/chamomile/diploid/hifi+hic/hap1_hic/assembly/chamomile2/chamomile2.hic.hap1.p_ctg_DpnII_fragments_30bp_split.nuc.txt'

hic_plots(rds="hic_map_v1.Rds", assembly=assembly_v2,
	cores=40, species="chamomile", nuc=snuc) -> hic_map_v1

nuc <- '/filer-dg/agruppen/dg2/cho/chamomile/diploid/hifi+hic/hap1_hic/assembly/chamomile2/chamomile2.hic.hap1.p_ctg_DpnII_fragments_30bp_split.nuc.txt'
add_psmol_fpairs(assembly=assembly_v2, hic_map=hic_map_v1, nucfile=nuc)->hic_map_v1_l

bin_hic_step(hic=hic_map_v1_l$links, frags=hic_map_v1_l$frags, binsize=1e6, chrlen=hic_map_v1_l$chrlen, cores=14)->hic_map_v1_l$hic_1Mb

f <- "chamomile2_interchromosomal.png"
interchromosomal_matrix_plot(hic_map=hic_map_v1_l, file=f, species="chamomile")

 write_hic_map(rds="hic_map_v1.Rds", file="hic_map_v1.xlsx", species="chamomile")

ssh cho@vm-123
cd /srv/shiny-server/Blocks/Theme

HiC inspector
user: dg_member
password: 3aeC1+Irt


read_hic_map(rds="hic_map_v1.Rds", file="hic_map_v2.xlsx") -> nmap
diff_hic_map(rds="hic_map_v1.Rds", nmap, species="chamomile")
hic_map(species="chamomile", agp_only=T, map=nmap) -> hic_map_v2

saveRDS(hic_map_v2, file="hic_map_v2.Rds")

snuc <- '/filer-dg/agruppen/dg2/cho/chamomile/hifi+hic/hap1_hic/assembly/chamomile2/chamomile2.hic.hap1.p_ctg_DpnII_fragments_30bp_split.nuc.txt'
hic_plots(rds="hic_map_v2.Rds", assembly=assembly_v2,
	cores=30, species="chamomile", nuc=snuc) -> hic_map_v2

hic_map_v2$agp[agp_chr != "chrUn" & gap == F][, .("Ncontig" = .N,
	"N50 (Mb)"=round(n50(scaffold_length/1e6), 1),
	"max_length (Mb)"=round(max(scaffold_length/1e6),1),
	"min_length (kb)"=round(min(scaffold_length/1e3),1)), key=agp_chr] -> res
hic_map_v2$chrlen[, .(agp_chr, "length (Mb)"=round(length/1e6, 1))][res, on="agp_chr"] -> res
setnames(res, "agp_chr", "chr")

hic_map_v2$agp[gap == F & agp_chr != "chrUn"][, .("Ncontig" = .N,
	"N50 (Mb)"=round(n50(scaffold_length/1e6), 1),
	"max_length (Mb)"=round(max(scaffold_length/1e6),1),
	"min_length (kb)"=round(min(scaffold_length/1e3),1))][, agp_chr := "1-10"] -> res2
hic_map_v2$chrlen[, .(agp_chr="1-10", "length (Mb)"=round(sum(length)/1e6, 1))][res2, on="agp_chr"] -> res2
setnames(res2, "agp_chr", "chr")

hic_map_v2$agp[gap == F & agp_chr == "chrUn"][, .(chr="un", "Ncontig" = .N,
	"N50 (Mb)"=round(n50(scaffold_length/1e6), 1),
	"max_length (Mb)"=round(max(scaffold_length/1e6),1),
	"min_length (kb)"=round(min(scaffold_length/1e3),1),
	"length (Mb)"=sum(scaffold_length/1e6))] -> res3

rbind(res, res2, res3) -> res

write.xlsx(res, file="hic_map_v2_pseudomolecule_stats.xlsx") ##write.xlsx 안돼서 write.csv2 하고 수정함

fasta <- '/filer-dg/agruppen/dg2/cho/chamomile/hifi+hic/hap1_hic/assembly/chamomile2/chamomile2.hic.hap1.p_ctg.fa'
sink("pseudomolecules_v1.log") 
compile_psmol(fasta=fasta, output="pseudomolecules_v1",
	hic_map=hic_map_v2, assembly=assembly_v2, cores=30)
sink() 

##hic plot with rearranged sequence
source("/filer-dg/agruppen/dg2/cho/chamomile/diploid/hifi+hic/hap1_hic/bitbucket/R/pseudomolecule_construction.R")
read_hic_map(rds="/filer-dg/agruppen/dg2/cho/chamomile/diploid/hifi+hic/hap1_hic/pseudomolecules/hic_map_v1.Rds", file="/filer-dg/agruppen/dg2/cho/chamomile/diploid/hifi+hic/hap1_hic/pseudomolecules/hic_map_v5.xlsx") -> nmap 
diff_hic_map(rds="/filer-dg/agruppen/dg2/cho/chamomile/diploid/hifi+hic/hap1_hic/pseudomolecules/hic_map_v1.Rds", nmap, species="chamomile")
hic_map(species="chamomile", agp_only=T, map=nmap) -> hic_map_v5
saveRDS(hic_map_v5, file="hic_map_v5.Rds")
snuc <- '/filer-dg/agruppen/dg2/cho/chamomile/diploid/hifi+hic/hap1_hic/assembly/chamomile2/chamomile2.hic.hap1.p_ctg_DpnII_fragments_30bp_split.nuc.txt'
readRDS('/filer-dg/agruppen/dg2/cho/chamomile/diploid/hifi+hic/hap1_hic/pseudomolecules/assembly_v2.Rds') -> assembly_v2
hic_plots(rds="hic_map_v5.Rds", assembly=assembly_v2, cores=30, species="chamomile", nuc=snuc) -> hic_map_v5
saveRDS(hic_map_v5, file="hic_map_v5.Rds")

.libPaths(c("/home/mascher/tmp/3.5.1", "/filer-dg/agruppen/seq_shared/mascher/Rlib/3.5.1", "/opt/Bio/R_LIBS/3.5.1"))
source('/filer-dg/agruppen/dg7/cho/imperata_cylindrica/WK14/pseudohaploid_assembly/bitbucket/R/pseudomolecule_construction.R')
read_hic_map(rds="hic_map_v1.Rds", file="hic_map_v6.xlsx") -> nmap 
hic_map(species="chamomile", agp_only=T, map=nmap) -> hic_map_v6
saveRDS(hic_map_v6, file="hic_map_v6.Rds")
snuc <- '/filer-dg/agruppen/dg2/cho/chamomile/diploid/integrated_assembly/hap1_hic/assembly/chamomile2/chamomile2.hic.hap1.p_ctg_DpnII_fragments_30bp_split.nuc.txt'
readRDS('assembly_v2.Rds') -> assembly_v2
hic_plots(rds="hic_map_v6.Rds", assembly=assembly_v2, cores=30, species="chamomile", nuc=snuc) -> hic_map_v6
saveRDS(hic_map_v6, file="hic_map_v6.Rds")

read_hic_map(rds="hic_map_v1.Rds", file="hic_map_v7.xlsx") -> nmap 
hic_map(species="chamomile", agp_only=T, map=nmap) -> hic_map_v7
saveRDS(hic_map_v7, file="hic_map_v7.Rds")
snuc <- '/filer-dg/agruppen/dg2/cho/chamomile/diploid/integrated_assembly/hap1_hic/assembly/chamomile2/chamomile2.hic.hap1.p_ctg_DpnII_fragments_30bp_split.nuc.txt'
readRDS('assembly_v2.Rds') -> assembly_v2
hic_plots(rds="hic_map_v7.Rds", assembly=assembly_v2, cores=30, species="chamomile", nuc=snuc) -> hic_map_v7
saveRDS(hic_map_v7, file="hic_map_v7.Rds")



fasta <- '/filer-dg/agruppen/dg2/cho/chamomile/diploid/integrated_assembly/hap1_hic/assembly/chamomile2/chamomile2.hic.hap1.p_ctg.fa'
sink("pseudomolecules_v5.log") 
compile_psmol(fasta=fasta, output="pseudomolecules_v5",
	hic_map=hic_map_v7, assembly=assembly_v2, cores=30)
sink() 

##centromere checking with 5Mb bin size
read_hic_map(rds="/filer-dg/agruppen/dg2/cho/chamomile/diploid/hifi+hic/hap1_hic/pseudomolecules/hic_map_v1.Rds", file="/filer-dg/agruppen/dg2/cho/chamomile/diploid/hifi+hic/hap1_hic/pseudomolecules/hic_map_v4.xlsx") -> nmap 
hic_map(species="chamomile", agp_only=T, map=nmap) -> hic_map_v4
readRDS('/filer-dg/agruppen/dg2/cho/chamomile/diploid/hifi+hic/hap1_hic/pseudomolecules/assembly_v2.Rds') -> assembly_v2
hic_plots(rds="hic_map_v4.Rds", assembly=assembly_v2, cores=30, species="chamomile", nuc=snuc) -> hic_map_v4
nuc <- '/filer-dg/agruppen/dg2/cho/chamomile/diploid/hifi+hic/hap1_hic/assembly/chamomile2/chamomile2.hic.hap1.p_ctg_DpnII_fragments_30bp_split.nuc.txt'
add_psmol_fpairs(assembly=assembly_v2, hic_map=hic_map_v4, nucfile=nuc)->hic_map_v2_l

bin_hic_step(hic=hic_map_v2_l$links, frags=hic_map_v2_l$frags, binsize=5e6, chrlen=hic_map_v2_l$chrlen, cores=14)->hic_map_v2_l$hic_5Mb
f <- "chamomile2_interchromosomal_5Mb.png"
interchromosomal_matrix_plot(hic_map=hic_map_v2_l, file=f, species="chamomile")

saveRDS(hic_map_v2_l, file="hic_map_v2_l.Rds")

###annotation 
####with helixer
singularity pull docker://gglyptodon/helixer-docker:helixer_v0.3.1_cuda_11.2.0-cudnn8
singularity run -B /filer-dg/agruppen --nv helixer-docker_helixer_v0.3.1_cuda_11.2.0-cudnn8.sif Helixer.py --fasta-path annotation_result/pseudomolecules_v3.fasta --lineage land_plant --gff-output-path annotation_result/chamomile_pseudomolecules_v3.gff3 &
sbatch --auks=yes -p gpu --cpus-per-task=1 --mem-per-cpu=1G --constraint=avx test.sh

#####################################################################################################################################################################################################################
###########################################################################tetraploid##(don't use this, use phased assembly pipeline)##################################################################################################
######################################################################################################################################################################################################
################################################################################################################################################################################
##genome size expectation with kmer(c 10 mem 400G)
jellyfish count -m 51 -s 1G -t 20 -C <(zcat ../raw_reads/hifi_reads/*.fastq.gz) &
jellyfish histo -t 10 mer_counts.jf > 51mer.histo

.libPaths(c("/home/mascher/tmp/3.5.1", "/filer-dg/agruppen/seq_shared/mascher/Rlib/3.5.1", "/opt/Bio/R_LIBS/3.5.1"))
library("findGSE")
findGSE(histo="51mer.histo",sizek=51, outdir="51_het", exp_hom=30)


#####guide_map
source("/filer-dg/agruppen/dg2/cho/chamomile/tetraploid/bitbucket/R/pseudomolecule_construction.R")
.libPaths(c("/home/mascher/tmp/3.5.1", "/filer-dg/agruppen/seq_shared/mascher/Rlib/3.5.1", "/opt/Bio/R_LIBS/3.5.1"))

f<-'/filer-dg/agruppen/dg2/cho/chamomile/tetraploid/guide_map/pseudomolecules_v3_singlecopy_100bp.bed'
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
pp[, sorted_alphachr := paste("chr",popseq_alphachr,sep="")]
pp[, sorted_chr := popseq_chr]
pp[, sorted_lib := sorted_alphachr]

# Save as RDS file for use in TRITEX pipeline
saveRDS(pp, "/filer-dg/agruppen/dg2/cho/chamomile/tetraploid/guide_map/tetra_chamomile_pseudopopseq.Rds")

/opt/Bio/minimap2/2.17/bin/minimap2 -t 20 -2 -I 20G -K5G -x asm5 /filer-dg/agruppen/dg2/cho/chamomile/tetraploid/assembly/tetra_chamomile/tetra_chamomile.hic.p_ctg.gfa.fa pseudomolecules_v3_singlecopy_100bp.fasta | gzip > guidemap.paf.gz 2> guidemap.err &

.libPaths(c("/home/mascher/tmp/3.5.1", "/filer-dg/agruppen/seq_shared/mascher/Rlib/3.5.1", "/opt/Bio/R_LIBS/3.5.1"))
source("/filer-dg/agruppen/dg2/cho/chamomile/tetraploid/bitbucket/R/pseudomolecule_construction.R")
readRDS('tetra_chamomile_pseudopopseq.Rds')-> popseq
f<- '/filer-dg/agruppen/dg2/cho/chamomile/tetraploid/assembly/tetra_chamomile/tetra_chamomile.hic.p_ctg.gfa.fa.fai'
fread(f, head=F, select=1:2, col.names=c("scaffold", "length")) -> fai
f <- '/filer-dg/agruppen/dg2/cho/chamomile/tetraploid/guide_map/guidemap.paf.gz'
read_morexaln_minimap(paf=f, popseq=popseq, minqual=30, minlen=500, prefix=F) -> morexaln
dir <- '/filer-dg/agruppen/dg2/cho/chamomile/tetraploid/mapping/hic'
fread(paste('find', dir, '| grep "hic_fragment_pairs.tsv.gz$" | xargs zcat'),
       header=F, col.names=c("scaffold1", "pos1", "scaffold2", "pos2")) -> fpairs
init_assembly(fai=fai, cssaln=morexaln, fpairs=fpairs) -> assembly
anchor_scaffolds(assembly = assembly, popseq=popseq, species="chamomile") -> assembly
add_hic_cov(assembly, binsize=1e4, binsize2=1e6, minNbin=50, innerDist=3e5, cores=40) -> assembly
saveRDS(assembly, file="assembly.Rds")

assembly$info[length >= 1e6, .(scaffold, length)][order(-length)] -> s
plot_chimeras(assembly=assembly, scaffolds=s, species="chamomile", refname="chamomile", autobreaks=F, mbscale=1, file="assembly_1Mb.pdf", cores=40)

assembly$info[length >= 1e6 & mri <= -3, .(scaffold, length)][order(-length)] -> s

plot_chimeras(assembly=assembly, scaffolds=s, species="chamomile", refname="chamomile_diploid", autobreaks=F, mbscale=1, file="assembly_chimeras.pdf", cores=30)


assembly-> assembly_v2

saveRDS(assembly_v2, file="assembly_v2.Rds")

fbed <- '/filer-dg/agruppen/dg2/cho/chamomile/tetraploid/assembly/tetra_chamomile/tetra_chamomile.hic.p_ctg.gfa_DpnII_fragments_30bp.bed'
read_fragdata(info=assembly_v2$info, file=fbed) -> frag_data

# Consider only contigs >= 1 Mb first
frag_data$info[!is.na(hic_chr) & length >= 1e6, .(scaffold, nfrag, length, chr=hic_chr, cM=popseq_cM)] -> hic_info

hic_map(info=hic_info, assembly=assembly_v2, frags=frag_data$bed,
	species="chamomile_tetraploid", ncores=30, min_nfrag_scaffold=30, 
	max_cM_dist = 1000, binsize=2e5, min_nfrag_bin=10, gap_size=100) -> hic_map_v1

saveRDS(hic_map_v1, file="hic_map_v1.Rds")

# Make the Hi-C plots, files will be placed in your working directory.
snuc <- '/filer-dg/agruppen/dg2/cho/chamomile/tetraploid/assembly/tetra_chamomile/tetra_chamomile.hic.p_ctg.gfa_DpnII_fragments_30bp_split.nuc.txt'

hic_plots(rds="hic_map_v1.Rds", assembly=assembly_v2,
	cores=40, species="chamomile_tetraploid", nuc=snuc) -> hic_map_v1

nuc <- '/filer-dg/agruppen/dg2/cho/chamomile/tetraploid/assembly/tetra_chamomile/tetra_chamomile.hic.p_ctg.gfa_DpnII_fragments_30bp_split.nuc.txt'
add_psmol_fpairs(assembly=assembly_v2, hic_map=hic_map_v1, nucfile=nuc)->hic_map_v1_l

bin_hic_step(hic=hic_map_v1_l$links, frags=hic_map_v1_l$frags, binsize=1e6, chrlen=hic_map_v1_l$chrlen, cores=14)->hic_map_v1_l$hic_1Mb

f <- "tetra_chamomile_interchromosomal.png"
interchromosomal_matrix_plot(hic_map=hic_map_v1_l, file=f, species="chamomile_tetraploid")

 write_hic_map(rds="hic_map_v1.Rds", file="hic_map_v1.xlsx", species="chamomile_tetraploid")
