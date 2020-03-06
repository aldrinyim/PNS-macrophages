library(Seurat)
library(dplyr)
library(Matrix)
library(monocle)

neg.data <- Read10X(data.dir = "/Users/aldrinyim/Documents/WashU/Flt3-data/neg-clec7a/raw_gene_bc_matrices/grcm38/")
pos.data <- Read10X(data.dir = "/Users/aldrinyim/Documents/WashU/Flt3-data/pos-clec7a/raw_gene_bc_matrices/grcm38/")

neg <- CreateSeuratObject(raw.data = neg.data, project = "Flt3-neg", min.cells = 3, min.genes = 200)
neg@meta.data$flt3 <- "Flt3-neg"
neg.mito.genes <- grep(pattern = "^mt-", x = rownames(x = neg@data), value = TRUE)
neg.percent.mito <- Matrix::colSums(neg@raw.data[neg.mito.genes, ])/Matrix::colSums(neg@raw.data)
neg <- AddMetaData(object = neg, metadata = neg.percent.mito, col.name = "percent.mito")
VlnPlot(object = neg, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
par(mfrow = c(1, 2))
GenePlot(object = neg, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = neg, gene1 = "nUMI", gene2 = "nGene")
neg <- FilterCells(object = neg, subset.names = c("nGene", "percent.mito"), 
                   low.thresholds = c(500, -Inf), high.thresholds = c(4000, 0.1))
neg <- NormalizeData(object = neg)
neg <- ScaleData(object = neg)
neg <- FindVariableGenes(object = neg)

neg

pos <- CreateSeuratObject(raw.data = pos.data, project = "Flt3-pos", min.cells = 3, min.genes = 200)
pos@meta.data$flt3 <- "Flt3-pos"
pos.mito.genes <- grep(pattern = "^mt-", x = rownames(x = pos@data), value = TRUE)
pos.percent.mito <- Matrix::colSums(pos@raw.data[pos.mito.genes, ])/Matrix::colSums(pos@raw.data)
pos <- AddMetaData(object = pos, metadata = pos.percent.mito, col.name = "percent.mito")
VlnPlot(object = pos, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
par(mfrow = c(1, 2))
GenePlot(object = pos, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = pos, gene1 = "nUMI", gene2 = "nGene")
pos <- FilterCells(object = pos, subset.names = c("nGene", "percent.mito"), 
                   low.thresholds = c(500, -Inf), high.thresholds = c(4000, 0.1))
pos <- NormalizeData(object = pos)
pos <- ScaleData(object = pos)
pos <- FindVariableGenes(object = pos)

pos

neg@meta.data[, "protocol"] <- "Flt3-neg"
pos@meta.data[, "protocol"] <- "Flt3-pos"

second_merge <- MergeSeurat(object1=pos,object2=neg, add.cell.id1="Flt3-pos",add.cell.id2="Flt3-neg")
second_merge

second_merge <- NormalizeData(object = second_merge, normalization.method = "LogNormalize", 
                              scale.factor = 10000)

par(mfrow = c(1, 1))
second_merge <- FindVariableGenes(object = second_merge, mean.function = ExpMean, dispersion.function = LogVMR, 
                                  x.low.cutoff = 0.05, x.high.cutoff = 5, y.cutoff = 0.5)

second_merge <- ScaleData(object = second_merge, vars.to.regress = c("nUMI", "percent.mito"))

second_merge <- RunPCA(object = second_merge, pc.genes = second_merge@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 30)
VizPCA(object = second_merge, pcs.use = 1:2)
PCAPlot(object = second_merge, dim.1 = 1, dim.2 = 2)

PCHeatmap(object = second_merge, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
PCHeatmap(object = second_merge, pc.use = 1:10, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)

second_merge <- JackStraw(object = second_merge, num.replicate = 100)
JackStrawPlot(object = second_merge, PCs = 1:12)

second_merge <- FindClusters(object = second_merge, reduction.type = "pca", dims.use = 1:6, 
                             resolution = 0.3, print.output = 0, algorithm = 1, save.SNN = TRUE, force.recalc = TRUE)

second_merge <- RunTSNE(object = second_merge, dims.use = 1:6, do.fast = T)
TSNEPlot(object = second_merge, group.by = "protocol", do.return = TRUE, pt.size = 0.8)
TSNEPlot(object = second_merge, do.return = TRUE, pt.size = 0.8, do.label=F)

second_merge <- SubsetData(object = second_merge, cells.use = WhichCells(object = second_merge, ident=c(0,1,2,3,4)))
second_merge

second_merge <- readRDS("/Users/aldrinyim/Box Sync/10X genomics/2019/Flt3/merged_analysis/non-CCA-merged-flt3negpos.RData")
saveRDS(second_merge,file="/Users/aldrinyim/Box Sync/10X genomics/2019/Flt3/merged_analysis/non-CCA-merged-flt3negpos.RData")

# Store the cluster identities in a new column in 'object@meta.data'
second_merge <- Seurat::StashIdent(object = second_merge, save.name = "clusterID")

# Set the experimental condition as cell identity
second_merge <- Seurat::SetAllIdent(object = second_merge, id = "flt3")

# Subset data for 'CTRL' cells
second_merge.neg <- Seurat::SubsetData(object = second_merge, ident.use = "Flt3-neg")

# Subset data for 'STIM' cells
second_merge.pos <- Seurat::SubsetData(object = second_merge, ident.use = "Flt3-pos")

# Restore identities stored in 'object@meta.data$clusterID'
second_merge.neg <- Seurat::SetAllIdent(object = second_merge.neg, id = "clusterID")
second_merge.pos <- Seurat::SetAllIdent(object = second_merge.pos, id = "clusterID")

Seurat::TSNEPlot(object = second_merge.neg, group.by = "protocol",do.return = TRUE, pt.size = 0.8, do.label=F)
Seurat::TSNEPlot(object = second_merge.pos, group.by = "protocol", do.return = TRUE, pt.size = 0.8, do.label=F)

cluster_defining_genes = c("Mki67","Top2a","Ccr2","Fos","Lyve1","Ly6e")
FeaturePlot(object = second_merge, features.plot = cluster_defining_genes, cols.use = colorRampPalette(c("navyblue","Cyan","yellow","red"))(4), pt.size = 0.4, no.axes=T)

