library(Seurat)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(SeuratWrappers)
library(monocle3)

####### gRNA analysis ####
seuratObj <- readRDS("SeuratObject.rds")

# Find cells with more than a guide for AGM library cells
# and cells with a guide (NT) for NT cells
seuratObj.gRNA<- subset(seuratObj, 
                        subset = (nFeature_gRNA>1 & libraryID >=3)|
                          (nFeature_gRNA==1 & libraryID <=2))

DimPlot(seuratObj.gRNA, split.by = "libraryID", pt.size = 0.8)

##### COMBINE gRNA #####

targetedGenes <- c("ZNF124", "ZNF33A", "SOX6", "NR4A1", "NFAT5", "SMAD7", 
                   "GATA2", "TFDP2", "RUNX1T1", "NT")
featureNames <- rownames(seuratObj.gRNA[["gRNA"]])

# We go through our targeted genes one at a time, 
# and calculate the sum of the expression of all the guides for that gene
combinedExpression <- sapply(targetedGenes, function(gene){
  # We use grep to get all the features containing the gene name
  # GetAssayData will return a matrix with guides on rows and cells on columns
  # We can then sum over the columns to get our combined expression
  expr <- as.matrix(GetAssayData(seuratObj.gRNA, assay = "gRNA")[grep(gene, featureNames),])
  if (gene != "NT")
    colSums(expr)
  # The NT guide is an exception, since it will returns a single column
  # with cells on the rows. We sum rows in this case
  else
    rowSums(expr)
})

# combinedExpression now holds a matrix of counts, with cells on rows and
# guides on the columns. We transpose it to be able to add it as a new assay
# in the Seurat object
combinedExpression <- t(combinedExpression)
rownames(combinedExpression) <- paste0(rownames(combinedExpression), "_gRNA")
# Add a new assay combining multiple gRNA for the same target into one
seuratObj.gRNA[["combined_gRNA"]] <- CreateAssayObject(counts = combinedExpression)

DimPlot(seuratObj.gRNA, reduction = "umap", 
        split.by = "libraryID", ncol = 2)

# Save to file
saveRDS(seuratObj.gRNA, "SeuratObject_combined_guides.rds")


####### TARGET gene analysis #########

target <- c ("ZNF124", "SOX6", "NR4A1", "SMAD7",  "GATA2", 
             "ZNF33A", "TFDP2", "RUNX1T1", "NFAT5")

# Extract expression of target genes from Seurat object and put in a dataframe
target_expr <- data.frame(t(as.matrix(GetAssayData(seuratObj.gRNA)[target,])))

target_expr_long <- target_expr %>%
  # Put barcode into a column
  rownames_to_column("Barcode") %>%
  # Transform from wide data to long data
  pivot_longer(-Barcode, names_to = "Gene", values_to = "Expression") %>% 
  # Add a Library column 
  mutate(Library = substr(Barcode, 18, 18))

# Calculate mean by gene and library
mean_expr <- target_expr_long %>% 
  group_by(Gene, Library) %>% 
  summarise(MeanExpr = mean(Expression))

mean_expr$Gene <- factor(mean_expr$Gene, 
                         levels = sort(target, decreasing = TRUE))
nt_fold <- mean_expr[mean_expr$Library == 2,]$MeanExpr / 
  mean_expr[mean_expr$Library == 1,]$MeanExpr
agm_fold <- mean_expr[mean_expr$Library == 4,]$MeanExpr / 
  mean_expr[mean_expr$Library == 3,]$MeanExpr

fold_change <- data.frame(FoldChange = c(nt_fold, agm_fold),
                          Library = rep(c("NT", "AGM"), each = length(nt_fold)),
                          Gene = factor(target, levels = sort(target, decreasing = TRUE)))
fold_change$Library <- factor(fold_change$Library, levels = c("NT", "AGM"))

ggplot(fold_change, aes(x = Library, y = Gene, fill = FoldChange)) +
  geom_tile() +
  scale_fill_gradientn(colours = c("#C9C9C9", "#1401F5"))

######## Cluster contribution ############

pt <- table(seuratObj.gRNA$libraryID, Idents(seuratObj.gRNA))
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) + theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Cluster") +
  ylab("Proportion") +
  theme(legend.title = element_blank())

######## comaprison of arterial cells between libraries #########

arterial.cell <- subset(seuratObj.gRNA, idents = 4)
Idents(arterial.cell) <- "libraryID"
DimPlot(arterial.cell, split.by = "libraryID", ncol = 2)

arterial.cell.markers <- FindAllMarkers(arterial.cell, only.pos = T)
write.csv(arterial.cell.markers, "arterial_markers_activation.csv")

######## subset hemogenic and blood cells #############

EHT.cells <- subset(seuratObj.gRNA, idents = c(1,3))
DimPlot(EHT.cells, split.by = "libraryID", ncol = 2)

##### Cell cycle ####

# Cell cycle genes
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

# Note: CellCycleScoring does not work with SCT, it needs
# normalisation through NormalizeData/FindVariableFeatures/ScaleData
# See https://satijalab.org/seurat/articles/cell_cycle_vignette.html 
# There are several issues on Seurat's GitHub page about this
# We save the object with a different name to be used only for this
EHT.cells.forCC <- NormalizeData(EHT.cells)
EHT.cells.forCC <- FindVariableFeatures(EHT.cells.forCC, selection.method = "vst")
EHT.cells.forCC <- ScaleData(EHT.cells.forCC, features = rownames(EHT.cells.forCC))

EHT.cells.forCC <- CellCycleScoring(EHT.cells.forCC, s.features = s.genes, 
                                    g2m.features = g2m.genes, set.ident = TRUE)
DimPlot(EHT.cells.forCC, split.by = "libraryID", ncol = 2, pt.size = 1) +
  xlab(expression(UMAP[1])) +
  ylab(expression(UMAP[2])) 

# Whole dataset
cells.for.CC <- NormalizeData(seuratObj.gRNA)
cells.for.CC <- FindVariableFeatures(cells.for.CC, selection.method = "vst")
cells.for.CC <- ScaleData(cells.for.CC, features = rownames(cells.for.CC))
cells.for.CC <- CellCycleScoring(cells.for.CC, s.features = s.genes, 
                                 g2m.features = g2m.genes, set.ident = TRUE)
DimPlot(cells.for.CC, split.by = "libraryID", ncol = 2, pt.size = 0.5) +
  xlab(expression(UMAP[1])) +
  ylab(expression(UMAP[2])) 
# Move cell cycle data back into Seurat object
seuratObj.gRNA$CellCycle <- Idents(cells.for.CC)
seuratObj.gRNA$CellCycle <- factor(seuratObj.gRNA$CellCycle, levels = c("G1", "S", "G2M"))
rm(cells.for.CC)

tb <- table(Idents(EHT.cells.forCC), EHT.cells.forCC$libraryID)
tb <- apply(tb, 1, function(x){x/colSums(tb)*100})

perc.cc <- data.frame(tb) %>% 
  rownames_to_column(var = "LibraryID") %>% 
  pivot_longer(-LibraryID, names_to = "Stage", values_to = "Percentage")

perc.cc$Stage <- factor(perc.cc$Stage, levels = c("G1", "S", "G2M"))

ggplot(subset(perc.cc, LibraryID %in% c(2,4)), 
       aes(x = LibraryID, y = Percentage)) +
  geom_col(aes(fill = Stage), position="dodge") +
  theme(axis.title = element_text(size = 16))

# Save EHT cells to file
saveRDS(EHT.cells, "EHT_cells.rds")

######## MONOCLE 3 ########

##### EHT cells only ####

monocle_object <- as.cell_data_set(EHT.cells)
monocle_object <- cluster_cells(cds = monocle_object, reduction_method = "UMAP")
monocle_object <- learn_graph(monocle_object, use_partition = TRUE)
monocle_object <- order_cells(monocle_object)

pData(monocle_object)$CellLine <- factor(pData(monocle_object)$CellLine, levels = c("NT", "AGM"))
pData(monocle_object)$Treatment <- factor(pData(monocle_object)$Treatment, levels = c("-DOX", "+DOX"))

plot_cells(monocle_object,
           color_cells_by = "pseudotime",
           show_trajectory_graph = TRUE,
           cell_size = 0,
           cell_stroke = 2,
           alpha = 0.5) +
  facet_wrap(~libraryID)


seuratObj.gRNA.DOX <- subset(seuratObj.gRNA, Treatment == "+DOX")
seuratObj.gRNA.DOX$pseudotime <- pseudotime(monocle_object)

FeaturePlot(subset(seuratObj.gRNA.DOX, idents = c(1, 3)), "pseudotime") +
  scale_color_gradientn(colours = c("#15078A", "#D6556C", "#F1F326"), 
                        name = "Pseudotime") +
  ggtitle("") +
  xlab(expression(UMAP[1])) +
  ylab(expression(UMAP[2])) 


pseudotime_df <- data.frame(cell_cycle = seuratObj.gRNA.DOX$CellCycle,
                            pseudotime = seuratObj.gRNA.DOX$pseudotime,
                            library_id = seuratObj.gRNA.DOX$libraryID,
                            RUNX1 = GetAssayData(seuratObj.gRNA.DOX)["RUNX1",]
)

ggplot(pseudotime_df, aes(x = pseudotime, y = cell_cycle)) +
  geom_jitter(aes(col = cell_cycle), size = .5, height = 0.1) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  scale_y_discrete(limits = rev) +
  xlab("Pseudotime") +
  ylab("Cell cycle stage") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))
