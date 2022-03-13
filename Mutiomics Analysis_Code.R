#P. Augsornworawat
#Sample codes used for integration and analysis of multiome GEX and ATAC sequencing datasets
#this is a utilization of Signac package. Please refer to https://satijalab.org/signac/ for more details.
#Figure 1
#SC ISLETS


library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(S4Vectors)
library(ggplot2)
library(patchwork)
library(GenomicRanges)

library(BiocParallel)
library(chromVAR)
library(JASPAR2020)
library(motifmatchr)
library(SeuratWrappers)
library(monocle3)
library(Matrix)

set.seed(1234)


# Set working directory
setwd("...")

#Laod gene annotations for human hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "hg38"

#set processes for  Quantifying peaks in each dataset using FeatureMatrix
plan("multiprocess", workers = 2)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM

#READ RNA DATA
RNA.scislet1 <- Read10X_h5(".../scislet1/filtered_feature_bc_matrix.h5")
RNA.scislet2 <- Read10X_h5(".../sciselt2/filtered_feature_bc_matrix.h5")
RNA.scislet3 <- Read10X_h5(".../week2/filtered_feature_bc_matrix.h5")

#READ in peak sets
peaks.scislet1 <- read.table(file = ".../scislet2/atac_peaks.bed",col.names = c("chr", "start", "end"))
peaks.scislet2 <- read.table(file = ".../scislet3/atac_peaks.bed",col.names = c("chr", "start", "end"))
peaks.scislet3 <- read.table(file = ".../week2/atac_peaks.bed",col.names = c("chr", "start", "end")) 

# convert to genomic ranges
gr.scislet1 <- makeGRangesFromDataFrame(peaks.scislet1)
gr.scislet2 <- makeGRangesFromDataFrame(peaks.scislet2)
gr.scislet3 <- makeGRangesFromDataFrame(peaks.scislet3)

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(gr.scislet1, gr.scislet2, gr.scislet3))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

#CREATE FRAGMENT FILES
# load metadata
md.scislet1 <- read.table(file = ".../scislet1/per_barcode_metrics.csv",stringsAsFactors = FALSE,sep = ",", header = TRUE,row.names = 1)[-1, ] 
md.scislet2<- read.table(file = ".../scislet2/per_barcode_metrics.csv",stringsAsFactors = FALSE,sep = ",", header = TRUE,row.names = 1)[-1, ] 
md.scislet3 <- read.table(file = ".../week2/per_barcode_metrics.csv",stringsAsFactors = FALSE,sep = ",", header = TRUE,row.names = 1)[-1, ] 

# create fragment objects
frags.scislet1 <- CreateFragmentObject(path = ".../scislet1/atac_fragments.tsv.gz",cells = rownames(md.scislet1))
frags.scislet2 <- CreateFragmentObject(path = ".../scislet2/atac_fragments.tsv.gz",cells = rownames(md.scislet2))
frags.scislet3 <- CreateFragmentObject(path = ".../week2/atac_fragments.tsv.gz",cells = rownames(md.scislet3))

#Quantify peaks in each dataset
scislet1.counts <- FeatureMatrix(fragments = frags.scislet1, features = combined.peaks,cells = rownames(md.scislet1))
scislet2.counts <- FeatureMatrix(fragments = frags.scislet2, features = combined.peaks,cells = rownames(md.scislet2))
scislet3.counts <- FeatureMatrix(fragments = frags.scislet3, features = combined.peaks,cells = rownames(md.scislet3))

#Create Seurat object using RNA Data
scislet1 <- CreateSeuratObject(counts = RNA.scislet1$`Gene Expression`, assay = "RNA")
scislet2 <- CreateSeuratObject(counts = RNA.scislet2$`Gene Expression`, assay = "RNA")
scislet3 <- CreateSeuratObject(counts = RNA.scislet3$`Gene Expression`, assay = "RNA")

#remove non-common cells in multiome data
include_list1 <- colnames(scislet1)
include_list2 <- colnames(scislet2)
include_list3 <- colnames(scislet3)

scislet1.counts<-scislet1.counts[, include_list1]
scislet2.counts<-scislet2.counts[, include_list2]
scislet3.counts<-scislet3.counts[, include_list3]

#Add RNA assay to Seurat Object
scislet1[["ATAC"]] <- CreateChromatinAssay(counts = scislet1.counts, sep = c(":", "-"), fragments = frags.scislet1, annotation = annotation)
scislet2[["ATAC"]] <- CreateChromatinAssay(counts = scislet2.counts, sep = c(":", "-"), fragments = frags.scislet2, annotation = annotation)
scislet3[["ATAC"]] <- CreateChromatinAssay(counts = scislet3.counts, sep = c(":", "-"), fragments = frags.scislet3, annotation = annotation)

#QUALITY CONTROL
DefaultAssay(scislet1) <- "ATAC"
DefaultAssay(scislet2) <- "ATAC"
DefaultAssay(scislet3) <- "ATAC"

scislet1 <- NucleosomeSignal(scislet1)
scislet2 <- NucleosomeSignal(scislet2)
scislet3 <- NucleosomeSignal(scislet3)

scislet1 <- TSSEnrichment(scislet1)
scislet2 <- TSSEnrichment(scislet2)
scislet3 <- TSSEnrichment(scislet3)

#visualize QC features
VlnPlot(object = scislet1,features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),ncol = 4,pt.size = 0)
VlnPlot(object = scislet2,features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),ncol = 4,pt.size = 0)
VlnPlot(object = scislet3,features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),ncol = 4,pt.size = 0)

# filter out low quality cells
scislet1 <- subset( x = scislet1,subset = nCount_ATAC < 50000 & nCount_RNA < 50000 & nCount_ATAC > 1000 & nCount_RNA > 1000 & nucleosome_signal < 1.25 & TSS.enrichment > 2)
scislet2 <- subset( x = scislet2,subset = nCount_ATAC < 40000 & nCount_RNA < 40000 & nCount_ATAC > 1000 & nCount_RNA > 1000 & nucleosome_signal < 1.5 & TSS.enrichment > 2)
scislet3 <- subset( x = scislet3,subset = nCount_ATAC < 40000 & nCount_RNA < 40000 & nCount_ATAC > 1000 & nCount_RNA > 1000 & nucleosome_signal < 1.5 & TSS.enrichment > 2)


#Peak calling for each data
# call peaks using MACS2
peakscislet1 <- CallPeaks(scislet1 , macs2.path = ".../macs2")
peakscislet2 <- CallPeaks(scislet2, macs2.path = ".../macs2")
peakscislet3 <- CallPeaks(scislet3, macs2.path = ".../macs2")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peakscislet1 <- keepStandardChromosomes(peakscislet1, pruning.mode = "coarse")
peakscislet1 <- subsetByOverlaps(x = peakscislet1, ranges = blacklist_hg38_unified, invert = TRUE)
peakscislet2 <- keepStandardChromosomes(peakscislet2, pruning.mode = "coarse")
peakscislet2 <- subsetByOverlaps(x = peakscislet2, ranges = blacklist_hg38_unified, invert = TRUE)
peakscislet3 <- keepStandardChromosomes(peakscislet3, pruning.mode = "coarse")
peakscislet3 <- subsetByOverlaps(x = peakscislet3, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_countsscislet1 <- FeatureMatrix(fragments = Fragments(scislet1), features = peakscislet1, cells = colnames(scislet1))
macs2_countsscislet2 <- FeatureMatrix(fragments = Fragments(scislet2), features = peakscislet2, cells = colnames(scislet2))
macs2_countsscislet3 <- FeatureMatrix(fragments = Fragments(scislet3), features = peakscislet3, cells = colnames(scislet3))

# create a new assay using the MACS2 peak set and add it to the Seurat object
scislet1[["peaks"]] <- CreateChromatinAssay(counts = macs2_countsscislet1, fragments = frags.scislet1, annotation = annotation)
scislet2[["peaks"]] <- CreateChromatinAssay(counts = macs2_countsscislet2, fragments = frags.scislet2, annotation = annotation)
scislet3[["peaks"]] <- CreateChromatinAssay(counts = macs2_countsscislet3, fragments = frags.scislet3, annotation = annotation)


#MERGING of DATASETS
#save condition info
scislet1$dataset <- 'scislet1'
scislet2$dataset <- 'scislet2'
scislet3$dataset <- 'scislet3'

#rename cells with ID
scislet1 = RenameCells(scislet1, add.cell.id = "scislet1")
scislet2 = RenameCells(scislet2, add.cell.id = "scislet2")
scislet3 = RenameCells(scislet3, add.cell.id = "scislet3")

#RUN SCTransform
scislet1 <- SCTransform(scislet1,assay = "RNA")
scislet2 <- SCTransform(scislet2,assay = "RNA")
scislet3 <- SCTransform(scislet3,assay = "RNA")

#Fomd variable feature and Run PCA
scislet1 <- FindVariableFeatures(scislet1)
scislet2 <- FindVariableFeatures(scislet2)
scislet3 <- FindVariableFeatures(scislet3)

scislet1 <- RunPCA(scislet1)
scislet2 <- RunPCA(scislet2)
scislet3 <- RunPCA(scislet3)

#ANCHORING datasets and merge objects using integrate
# find integration anchors 
integratefeatures<- SelectIntegrationFeatures(object.list = list(scislet1, scislet2, scislet3),nfeatures = 500)
anchors <- PrepSCTIntegration(object.list = list(scislet1, scislet2, scislet3))
anchors <- FindIntegrationAnchors(object.list = anchors,reduction = "rpca",dims = 1:20,normalization.method = "SCT", anchor.features = integratefeatures,  assay = c('SCT', 'SCT','SCT'))

#INTEGRATE data and create a new merged object using SCT data
integratedRNA <- IntegrateData(anchorset = anchors, normalization.method = "SCT",dims = 1:20,new.assay.name = "integratedRNA")

#INTEGRATE data and create a new merged object using LSI 
integratedlsi<-integratedRNA
DefaultAssay(integratedlsi) <- "peaks"
integratedlsi <- FindTopFeatures(integratedlsi, min.cutoff = 5)
integratedlsi <- RunTFIDF(integratedlsi)
integratedlsi <- RunSVD(integratedlsi)
integratedlsi <- IntegrateEmbeddings(anchorset = anchors, new.reduction.name = "integratedLSI",reductions = integratedlsi@reductions[["lsi"]])

#Combined assays and reduction into one seurat object
integrated <- integratedlsi
integrated@assays[["integratedRNA"]] <- integratedRNA@assays[["integratedRNA"]]

#JOINT UMAP VISUALIZATION
DefaultAssay(integrated) <- "integratedRNA"
integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated, verbose = FALSE)

DefaultAssay(integrated) <- "peaks"
integrated <- FindTopFeatures(integrated, min.cutoff = 5)
integrated <- RunTFIDF(integrated)
integrated <- RunSVD(integrated)
# build a joint neighbor graph using both assays
integrated <- FindMultiModalNeighbors(object = integrated,reduction.list = list("pca", "integratedLSI"), dims.list = list(1:24, 2:24),modality.weight.name = "RNA.weight", verbose = TRUE)
# build a joint UMAP visualization
integrated <- RunUMAP(object = integrated, nn.name = "weighted.nn",assay = "SCT",verbose = TRUE)
integrated <- FindClusters(integrated,resolution = 0.06, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

#PLOT UMAP figures
DefaultAssay(integrated) <- "SCT"
DimPlot(integrated, label = TRUE, repel = TRUE, reduction = "umap") 
DimPlot(integrated, label = TRUE, repel = TRUE, reduction = "umap", group.by = 'dataset') 


#Find gene expression markers for each cluster
DefaultAssay(integrated) <- "SCT"
cluster0 <- FindMarkers(integrated, ident.1 = 0,ident.2 = NULL, logfc.threshold = 0.25, test.use = "wilcox", only.pos = TRUE) 
cluster1 <- FindMarkers(integrated, ident.1 = 1,ident.2 = NULL, logfc.threshold = 0.25, test.use = "wilcox", only.pos = TRUE) 
cluster2 <- FindMarkers(integrated, ident.1 = 2,ident.2 = NULL, logfc.threshold = 0.25, test.use = "wilcox", only.pos = TRUE) 
cluster3 <- FindMarkers(integrated, ident.1 = 3,ident.2 = NULL, logfc.threshold = 0.25, test.use = "wilcox", only.pos = TRUE) 
cluster4 <- FindMarkers(integrated, ident.1 = 4,ident.2 = NULL, logfc.threshold = 0.25, test.use = "wilcox", only.pos = TRUE) 
cluster5 <- FindMarkers(integrated, ident.1 = 5,ident.2 = NULL, logfc.threshold = 0.25, test.use = "wilcox", only.pos = TRUE) 
cluster6 <- FindMarkers(integrated, ident.1 = 6,ident.2 = NULL, logfc.threshold = 0.25, test.use = "wilcox", only.pos = TRUE) 
cluster7 <- FindMarkers(integrated, ident.1 = 7,ident.2 = NULL, logfc.threshold = 0.25, test.use = "wilcox", only.pos = TRUE) 

#RENAME clusters with identity based on expression markers
integrated <- RenameIdents(integrated,'0'='EC','1'='Beta','2'='Alpha','3'='Delta','4'='EB','5'='Mesenchyme','6'='Proliferating Alpha','7'='Exocrine')
integrated$celltype <- integrated@active.ident

#ANNOTATED UMAP COMBINED PLOT WITH COLORS
DefaultAssay(integrated) <- "SCT"
DimPlot(integrated, label = FALSE, repel = TRUE,group.by="celltype", reduction = "umap", cols= c("#FFA500", "#8E2C85", "#134D9C","#E0082D","#D877CF","#7D000E","#AAAAAA","#679327")) 
DimPlot(integrated, label = FALSE, repel = TRUE, reduction = "umap", group.by = 'dataset', cols = c("#2f4b7c", "#d45087", "#ffa600"), shuffle=TRUE) 

#ATAC PEAKS plots of key hormone genes
DefaultAssay(integrated) <- "peaks"
Idents(integrated) <- "celltype"
CoveragePlot(object = integrated,  region = "INS",  features = "INS", expression.assay = "SCT", extend.upstream = 5000,  extend.downstream = 2000) & scale_fill_manual(values = c("#FFA500", "#8E2C85", "#134D9C","#E0082D","#D877CF","#7D000E","#AAAAAA","#679327"))
CoveragePlot(object = integrated,  region = "GCG",  features = "GCG", expression.assay = "SCT", extend.upstream = 1000,  extend.downstream = 7000) & scale_fill_manual(values = c("#FFA500", "#8E2C85", "#134D9C","#E0082D","#D877CF","#7D000E","#AAAAAA","#679327"))
CoveragePlot(object = integrated,  region = "SST",  features = "SST", expression.assay = "SCT", extend.upstream = 5000,  extend.downstream = 2000) & scale_fill_manual(values = c("#FFA500", "#8E2C85", "#134D9C","#E0082D","#D877CF","#7D000E","#AAAAAA","#679327"))


#PERFORM ATAC or GEX only Analysis
#ATAC ONLY
DefaultAssay(integrated) <- "peaks"
integrated <- RunUMAP(object = integrated, reduction = 'integratedLSI', dims = 2:20, reduction.name = 'umap.atac')
integrated <- FindNeighbors(object = integrated, reduction = 'integratedLSI', dims = 2:20)
integrated <- FindClusters(object = integrated, verbose = FALSE, algorithm = 3,resolution = 0.27)
DimPlot(object = integrated, label = TRUE,reduction='umap.atac') 
gene.activities <- GeneActivity(integrated) 
#Add the gene activity (promoter accessibility) matrix to the Seurat object as a new assay and normalize it
integrated[['ActivityRNA']] <- CreateAssayObject(counts = gene.activities)
integrated <- NormalizeData(object = integrated, assay = 'ActivityRNA', normalization.method = 'LogNormalize', scale.factor = median(integrated$nCount_ActivityRNA))


#RNA ANALYSIS ONLY
DefaultAssay(integrated) <- "integratedRNA"
integrated <- FindVariableFeatures(integrated, nfeatures = 500)
integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated, npcs = 20)
integrated <- RunUMAP(integrated,reduction = 'pca', dims = 1:17, reduction.name = "umap.rna") 
integrated <- FindNeighbors(integrated, dims = 1:17)
integrated <- FindClusters(integrated, resolution = 0.2, algorithm = 3)

DefaultAssay(integrated) <- "SCT"
DimPlot(integrated, label = TRUE,reduction='umap.rna') 


#MOTIF activity analysis
#Adding motif information to the Seurat object using chromVAR package
register(SerialParam()) 
DefaultAssay(integrated) <- 'peaks'
integrated <- RunChromVAR(object = integrated, genome = BSgenome.Hsapiens.UCSC.hg38)
#Use this to find motif id from gene name
DefaultAssay(integrated) <- 'peaks'
motif.name <- ConvertMotifID(integrated, name = 'MAFA')
