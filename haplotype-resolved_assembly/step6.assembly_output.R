##
# Construct first Hi-C map
##

#read Hi-C fragment data
f <- '../FB20_005_1_hifiasm.p_utg_MboI_fragments_30bp.bed'
read_fragdata(info=assembly_v2$info, file=f)->frag_data

#read haplotype information
readRDS(file="FB20_005_1_assembly_v2_Hv_guide+HiClift+hap.Rds") -> hh

#Hi-C map for haplotype 1
hh[hap %in% c(1,3), .(scaffold, length, chr=as.integer(substr(chr, 4, 4)), cM = pos/1e6)] -> hic_info
frag_data$info[, .(scaffold, nfrag)][hic_info, on="scaffold"] -> hic_info

hic_map(info=hic_info, assembly=assembly_v2, frags=frag_data$bed, species="barley", ncores=21,
  min_nfrag_scaffold=30, max_cM_dist = 1000,
  binsize=2e5, min_nfrag_bin=10, gap_size=100)->hic_map_v1_hap1

saveRDS(hic_map_v1_hap1, file="FB20_005_1_hic_map_v1_hap1.Rds")

#Hi-C map for haplotype 2
hh[hap %in% c(2,3), .(scaffold, length, chr=as.integer(substr(chr, 4, 4)), cM = pos/1e6)] -> hic_info
frag_data$info[, .(scaffold, nfrag)][hic_info, on="scaffold"] -> hic_info

hic_map(info=hic_info, assembly=assembly_v2, frags=frag_data$bed, species="barley", ncores=21,
  min_nfrag_scaffold=30, max_cM_dist = 1000,
  binsize=2e5, min_nfrag_bin=10, gap_size=100)->hic_map_v1_hap2

saveRDS(hic_map_v1_hap2, file="FB20_005_1_hic_map_v1_hap2.Rds")

#create Hi-C plots for haplotype 1
snuc <- '/filer-dg/agruppen/dg7/fengj/hbulbosum_genome/FB20_005_1/FB20_005_1_hifiasm.p_utg_MboI_fragments_30bp_split.nuc.txt'
hic_plots(rds="FB20_005_1_hic_map_v1_hap1.Rds", cov=F, assembly=assembly_v2, cores=30, species="barley", nuc=snuc) -> hic_map_v1_hap1_l

#create Hi-C plots for haplotype 2
snuc <- '/filer-dg/agruppen/dg7/fengj/hbulbosum_genome/FB20_005_1/FB20_005_1_hifiasm.p_utg_MboI_fragments_30bp_split.nuc.txt'
hic_plots(rds="FB20_005_1_hic_map_v1_hap2.Rds", cov=F, assembly=assembly_v2, cores=30, species="barley", nuc=snuc) -> hic_map_v1_hap2_l

#combine both Hi-C maps
combine_hic(hic_map_v1_hap1, hic_map_v1_hap2, assembly_v2) -> ch
ch$hic_map -> hic_map_v1
ch$assembly -> assembly_v2_hap

#create inter-chromosomal plot for combined hap1+hap2 Hi-C haplotypes
nuc <- '/filer-dg/agruppen/dg7/fengj/hbulbosum_genome/FB20_005_1/FB20_005_1_hifiasm.p_utg_MboI_fragments_30bp_split.nuc.txt'
add_psmol_fpairs(assembly=assembly_v2_hap, hic_map=hic_map_v1, nucfile=nuc)->hic_map_v1_l

bin_hic_step(hic=hic_map_v1_l$links, frags=hic_map_v1_l$frags, binsize=1e6, chrlen=hic_map_v1_l$chrlen, cores=14)->hic_map_v1_l$hic_1Mb

f <- "FB20_005_1_hic_map_v1_interchromosomal.png"
interchromosomal_matrix_plot(hic_map=hic_map_v1_l, file=f, species="hordeum_bulbosum")

#write Hi-C map for editing
write_hic_map(rds="FB20_005_1_hic_map_v1_hap1.Rds", file="FB20_005_1_hic_map_v1_hap1.xlsx", species="barley")
write_hic_map(rds="FB20_005_1_hic_map_v1_hap2.Rds", file="FB20_005_1_hic_map_v1_hap2.xlsx", species="barley")

###### write fasta

readRDS('FB20_005_1_assembly_v2.Rds') -> assembly_v2
readRDS('FB20_005_1_hic_map_v1_hap1.Rds') -> hic_map_v1_hap1

readRDS('FB20_005_1_hic_map_v1_hap2.Rds') -> hic_map_v1_hap2


fasta <- '/filer-dg/agruppen/dg7/fengj/hbulbosum_genome/FB20_005_1/FB20_005_1_hifiasm.p_utg.fasta'
sink("FB20_005_1_pseudomolecules_v1_hap1.log")
compile_psmol(fasta=fasta, output="FB20_005_1_pseudomolecules_v1_hap1", hic_map=hic_map_v1_hap1, assembly=assembly_v2, cores=30)
sink()

sink("FB20_005_1_pseudomolecules_v1_hap2.log")
compile_psmol(fasta=fasta, output="FB20_005_1_pseudomolecules_v1_hap2", hic_map=hic_map_v1_hap2, assembly=assembly_v2, cores=30)
sink()


###Do mannel check for 2 haplotypes

#import Hi-C maps after editing, v1 -> v2
read_hic_map(rds="FB20_005_1_hic_map_v1_hap1.Rds", file="FB20_005_1_hic_map_v1_hap1_edit.xlsx") -> nmap
hic_map(species="barley", agp_only=T, map=nmap)->hic_map_v2_hap1
saveRDS(hic_map_v2_hap1, file="FB20_005_1_hic_map_v2_hap1.Rds")

read_hic_map(rds="FB20_005_1_hic_map_v1_hap2.Rds", file="FB20_005_1_hic_map_v1_hap2_edit.xlsx") -> nmap
hic_map(species="barley", agp_only=T, map=nmap)->hic_map_v2_hap2
saveRDS(hic_map_v2_hap2, file="FB20_005_1_hic_map_v2_hap2.Rds")

#proceed with Hi-C plots and repeat edit cycle if need be


#create Hi-C plots for haplotype 1
snuc <- '/filer-dg/agruppen/dg7/fengj/hbulbosum_genome/FB20_005_1/FB20_00#read Hi-C fragment data
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
