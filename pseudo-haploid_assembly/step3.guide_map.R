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
