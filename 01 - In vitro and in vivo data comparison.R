library(ggplot2)
library(Seurat)
library(dplyr)
library(Matrix)
library(tibble)
library(biomaRt)
library(reshape2)
library(cowplot)


################ LOAD THE DATA #######################
######################################################

##### Loading the raw files  for the human in vivo data to create the seurat merged object 

CS12_raw_counts<-read.table("GSM3993423_CS12_CH_UMI_raw.txt")
CS13_raw_counts<-read.table('GSM3993424_CS13_DA_UMI_raw.txt')
CS14_raw_counts<-read.table('GSM3993425_CS14_DA_UMI_raw.txt')


CS12 <- CreateSeuratObject(CS12_raw_counts, assay = "RNA", project = "AGM")
CS13 <- CreateSeuratObject(CS13_raw_counts, assay = "RNA", project = "AGM")
CS14 <- CreateSeuratObject(CS14_raw_counts, assay = "RNA", project = "AGM")

# set return.only.var.genes = F to avoid error in the merging
CS12 <- SCTransform(CS12, return.only.var.genes = F)
CS13 <- SCTransform(CS13, return.only.var.genes = F)
CS14 <- SCTransform(CS14, return.only.var.genes = F)

# store the carneige stages in meta.data$CS which I create here
CS12@meta.data$CS <- "CS12"
CS13@meta.data$CS <- "CS13"
CS14@meta.data$CS <- "CS14"


CS12_14 <- merge(CS12, y = c(CS13, CS14), 
                 add.cell.ids = c("CS12", "CS13", "CS14"),
                 project = "Hemogenic_endothelium", merge.data = T)

rm(CS12, CS14, CS13, CS12_raw_counts, CS13_raw_counts, CS14_raw_counts)

# 2_dimension reduction, clustering and marker

CS12_14 <- FindVariableFeatures(CS12_14)
CS12_14 <- RunPCA(CS12_14)
CS12_14 <- RunUMAP(object = CS12_14, dims = 1:20)
CS12_14 <- FindNeighbors(object = CS12_14, dims = 1:20, verbose = FALSE)
CS12_14 <- FindClusters(CS12_14, dims= 1:20)

# Rename identites to arterial, hemogenic arterial and venous

CS12_14 <- RenameIdents(CS12_14, "0" = "aHEC", "1" = "vEC", "2" = "aEC")

# load the in vitro data from Fidanza at al, Blood 2020 and normalise both dataset 

IVD_EndoHPC <- readRDS("Fidanza2020_Blood.RDS")
IVD_EndoHPC <- SCTransform(IVD_EndoHPC, return.only.var.genes = F)
CS12_14 <- SCTransform(CS12_14, return.only.var.genes = F)


#### FIND THE MARKERS OF THE AGM CLUSTERS (aHEC, vEC, aEC) ####
###########################################################

AGM.markers <- FindAllMarkers(CS12_14, min.pct = 0.5, only.pos = T)
# sort the marker of aHEC only
aHEC.AGMdata.markers <- AGM.markers$gene[AGM.markers$cluster=="aHEC"]


######## SORT TFs from aHEC markers ###########
###############################################

ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
go_id_TF = c('GO:0003700', 'GO:0000130', 'GO:0001071')
aHEC.TF <- getBM(attributes = c('hgnc_symbol'),
                 filters = c('go', 'hgnc_symbol'), 
                 values = list(go_id_TF, aHEC.AGMdata.markers), 
                 mart = ensembl, uniqueRows = T)

######### MERGE and SCT ##########

# merge the In vitro and in vivo data 
IVD_AGM <- merge(IVD_EndoHPC, y = CS12_14, add.cell.ids = c("IVD", "AGM"))
IVD_AGM <- FindVariableFeatures(IVD_AGM)
IVD_AGM <- RunPCA(IVD_AGM)
DimPlot(IVD_AGM)
IVD_AGM <- RenameIdents(IVD_AGM, "7" = "IVD_endo", "5" = "IVD_HPCs")
IVD_AGM@meta.data$Celltype <- IVD_AGM@active.ident

#run SCT again to normalise the data
IVD_AGM_SCT <- SCTransform(IVD_AGM, return.only.var.genes = F)
IVD_AGM_SCT <- FindVariableFeatures(IVD_AGM_SCT)
IVD_AGM_SCT <- RunPCA(IVD_AGM_SCT)
ElbowPlot(IVD_AGM_SCT)
IVD_AGM_SCT <- RunUMAP(IVD_AGM_SCT, dims = 1:20)
DimPlot(IVD_AGM_SCT)

# use Celltype as active identites
Idents(IVD_AGM_SCT) <- IVD_AGM_SCT@meta.data[["Celltype"]]
DimPlot(IVD_AGM_SCT)

############ SORTING TFs targets for cas9 activation ###############
####################################################################

# get expression values for TF genes from the object
TF_expression <- GetAssayData(IVD_AGM_SCT)[aHEC.TF$hgnc_symbol, Idents(IVD_AGM_SCT)=="aHEC"]
# calculate the average expression per gene which are on the rows
mean_TF_expression <- rowMeans(TF_expression)
# order the TF by expression average
aHEC.TF.ordered <- aHEC.TF$hgnc_symbol[order(mean_TF_expression, decreasing = T)]

# store the data plotted in the dotplot into dp 
# to get the % of cells that express each gene
# we also added some reference gene to the list
dp <- DotPlot(IVD_AGM_SCT, features = c(aHEC.TF.ordered, "GAPDH", "ACTB", "B2M"), dot.scale = 4) + 
  coord_flip() + 
  theme(axis.text.y = element_text(size = 8), axis.text.x = element_text(angle = 45, hjust = 1))

# sort the genes to be expressed in more than 50% of the aHEC cells
aHEC.TF.sorted <- dp$data[(dp$data$id=="aHEC" & 
                             dp$data$pct.exp>50),]

# sort the genes in aHEC.TF.sorted to be expressed in less than 25% of the IVD_endo cells
aHEC.TF.sorted <- dp$data[(dp$data$id=="IVD_endo" & 
                             dp$data$pct.exp<25  & 
                             dp$data$features.plot %in% aHEC.TF.sorted$features.plot),]




