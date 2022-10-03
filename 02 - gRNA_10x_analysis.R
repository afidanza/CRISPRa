library(ggplot2)
library(Seurat)

####### IMPORT DATA #######

# Change as appropriate
cellranger_output_dir <- "cellranger_output/"
# Read the matrix from the CellRanger output
# This is returned as a list with 2 elements:
# `Gene Expression`
# `CRISPR Guide`
data10x <- Read10X(cellranger_output_dir)

# Import gene expression data into Seurat. This (RNA) will be the default assay 
seuratObj <- CreateSeuratObject(data10x$`Gene Expression`, 
                                min.cells = 3, min.features = 200, 
                                project = "iSAM_2022")

gRNA <- CreateAssayObject(data10x$`CRISPR Guide`)
# Match the cells in the two assays
gRNA <- subset(gRNA, cells = colnames(seuratObj))
# Add the gRNA assay, to seuratObj
seuratObj[["gRNA"]] <- gRNA
# As a sanity check, we can check which assays are present 
# and what they contain in seuratObj@assays
rm(data10x)
rm(gRNA)

####### QC AND DATA FILTERING #######

# Calculate mitochondrial gene %
seuratObj[["percent.mt"]] <- PercentageFeatureSet(seuratObj, pattern = "^MT-")
# QC metrics
VlnPlot(seuratObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0) 
min_pct_mt <- 1
max_pct_mt <- 15
min_feat <- 1000
max_feat <- 7500

# Unfiltered
FeatureScatter(seuratObj, "nFeature_RNA", "percent.mt", pt.size = 0.01 ) +
  ylim(0,20) +
  geom_hline(yintercept = c(min_pct_mt, max_pct_mt)) +
  geom_vline(xintercept = c(min_feat, max_feat))

seuratObj <- subset(seuratObj, 
                    subset = percent.mt > min_pct_mt & percent.mt < max_pct_mt & 
                      nFeature_RNA > min_feat & nFeature_RNA < max_feat)

# Filtered
FeatureScatter(seuratObj, "nFeature_RNA", "percent.mt", pt.size = 0.01 ) +
  ylim(0,20) +
  geom_hline(yintercept = c(min_pct_mt, max_pct_mt)) +
  geom_vline(xintercept = c(min_feat, max_feat))


####### META DATA #######

# Create metadata for the Seurat object
# Library ID (the 18th character of the cell name is the library)
seuratObj[["libraryID"]] <- substr(colnames(seuratObj), 18, 18)

# Cell line (libraries 1 and 2 are NT, 3 and 4 AGM)
seuratObj[["CellLine"]] <- ifelse(seuratObj[["libraryID"]] == "1"|
                                    seuratObj[["libraryID"]] == "2"
                                  , "NT", "AGM")
# Treatment (libraries 1 and 3 are +Dox, 2 and 4 -Dox)
seuratObj[["Treatment"]] <- ifelse(seuratObj[["libraryID"]] == "1"|
                                     seuratObj[["libraryID"]] == "3" , 
                                   "-DOX", "+DOX")

# Plot the QC for each library
FeatureScatter(subset(seuratObj, subset = libraryID == 1),
               "nFeature_RNA", "percent.mt", pt.size = 0.01 ) +
  ylim(0,20) + ggtitle("NT-DOX")

FeatureScatter(subset(seuratObj, subset = libraryID == 2),
               "nFeature_RNA", "percent.mt", pt.size = 0.01 ) +
  ylim(0,20) + ggtitle("NT+DOX")

FeatureScatter(subset(seuratObj, subset = libraryID == 3),
               "nFeature_RNA", "percent.mt", pt.size = 0.01 ) +
  ylim(0,20) + ggtitle("AGM-DOX")

FeatureScatter(subset(seuratObj, subset = libraryID == 4),
               "nFeature_RNA", "percent.mt", pt.size = 0.01 ) +
  ylim(0,20) + ggtitle("AGM+DOX")

####### DATA NORMALIZATION ############
seuratObj <- SCTransform(seuratObj, vars.to.regress = "percent.mt")

####### MULTIDIMENSION REDUCTION ######

# PCA
seuratObj <- FindVariableFeatures(seuratObj)
seuratObj <- RunPCA(seuratObj)
DimPlot(seuratObj, reduction = "pca")
ElbowPlot(seuratObj, ndims = 50)

# UMAP
seuratObj <- RunUMAP(object = seuratObj, dims = 1:24)
DimPlot(seuratObj, reduction = "umap")

# Clustering
seuratObj <- FindNeighbors(object = seuratObj, dims = 1:24, verbose = FALSE)

seuratObj <- FindClusters(seuratObj, dims= 1:24, resolution = 0.2)
DimPlot(seuratObj, reduction = "umap", split.by = "libraryID", ncol = 2)

seuratObj.markers <- FindAllMarkers(seuratObj, only.pos = T)
write.csv(seuratObj.markers, file = "markers.csv")

# Save to file for loading in further steps
saveRDS(seuratObj, "SeuratObject.rds")