library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)
library(stringr)
library(mgcv)


sgRNA_dir <- '/data/CARD_iNDI/CRISPRi/PIPseq/20240125/sgRNA_sensitivity_5'
sc <- CreateSeuratObject(Read10X(sgRNA_dir))

# Filtering
sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
lower_bound <- quantile(sc[['nFeature_RNA']]$nFeature_RNA, 0.025)
upper_bound <- quantile(sc[['nFeature_RNA']]$nFeature_RNA, 0.975)
sc <- subset(sc, subset = nFeature_RNA > lower_bound & nFeature_RNA < upper_bound & percent.mt <= 5)


guides <- as.data.table(LayerData(object = sc, layer = 'counts.SNT Capture'), keep.rownames=TRUE)



# GetAssayData(object = sc.sgrna, assay = "RNA", slot = "counts")
# LayerData(object=sc.sgrna, layer='counts.Gene Expression', features='counts')
# gRNAcounts <- GetAssayData(sc, slot='counts')[grepl(rownames(sc), pattern = 'P1P2'),]

sc <- NormalizeData(sc, layer='counts.Gene Expression')
genes <- rownames(LayerData(object=sc, layer='counts.Gene Expression'))
sc <- ScaleData(sc, features = genes)

# Assign cells to gRNA groups following Datlinger et. al Nat. Methods 2017

# Datligner et. al defines cells based on the ratio of the highest count to all 
# the other gRNA counts. The paper requires that the most expressed gRNA is 3x
# higher than all the other gRNA summed together in order to call it as a 
# singlet. That multiple can be changed here:
singlet_fc <- 3
rn <- guides[, rn]
guides[, rn := NULL]

cell_assignments <- as.data.table(t(apply(guides, 2, function(x){
  total <- sum(x)
  if(total == 0){
    type <- 'None'
    gene <- 'None'
  } else if(max(x) > singlet_fc * (total - max(x))){
    type <- 'Singlet'
    gene <- rn[which(x == max(x))]
  } else {
    type <- 'Multiplet'
    gene <- 'None'
  }
  return(c(type, gene))
})))

# Wrangling
cell_assignments[, 'barcode' := colnames(guides)]
setnames(cell_assignments, old=c('V1','V2'), new=c('type','guide'))
cell_assignments[guide %like% 'non-targeting', 'gene' := 'nt']
cell_assignments[! guide %like% 'non-targeting' & guide != 'None', gene := tstrsplit(guide, split='-sgRNA')[[1]]]
setcolorder(cell_assignments, c('barcode','type','guide','gene'))


# Define high count genes to plot
top_n <- 12
high_counts <- cell_assignments[type=='Singlet', .N, by=gene][order(-N)][1:top_n, gene]

# Calculate targeting and non-targeting counts per gene
o <- foreach(i=high_counts, .combine='rbind', .errorhandling='remove') %do% {
  cells_with_this_guide <- cell_assignments[gene==i, barcode]
  cells_non_targeting <- cell_assignments[gene=='nt', barcode]

  kd_expression <- sc@assays$RNA$scale.data[i, cells_with_this_guide]
  nt_expression <- sc@assays$RNA$scale.data[i, cells_non_targeting]
  
  dt1 <- data.table(cell=names(kd_expression), gene=i, expression=kd_expression, sgRNA='targeting')
  dt2 <- data.table(cell=names(nt_expression), gene=i, expression=nt_expression, sgRNA='non-targeting')
  rbindlist(list(dt1, dt2))
}

# plot
g <- ggplot(data=o, aes(x=gene, y=expression, color=sgRNA)) + geom_boxplot()
ggsave(g, file='test.png', width=40, height=15, units='cm')
