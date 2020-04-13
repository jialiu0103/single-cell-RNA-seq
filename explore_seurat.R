#load seurat data
cells=readRDS("/projectnb/bf528/users/group4/others/GSM2230760_seurat.rda")

###############
#random walk
###############
#how many cells in each clusters
table(cells@active.ident)

#take a look at cluster2
head(subset(as.data.frame(cells@active.ident),cells@active.ident=="2"))

#build a tree
#Constructs a phylogenetic tree relating the 'average' cell from each identity class. 
#Tree is estimated based on a distance matrix constructed in either gene expression space or PCA space
cells_tree<-BuildClusterTree(cells)
Tool(object = cells_tree, slot = 'BuildClusterTree')
PlotClusterTree(cells_tree)

#Calculate the Barcode Distribution Inflection
cells_barcode<-CalculateBarcodeInflections(cells_tree)
SubsetByBarcodeInflections(cells_barcode)

#general look at umap
DimPlot(cells, reduction = "umap",label = TRUE)

###############
#Identify marker genes for each cluster
###############
# Finding differentially expressed features (cluster biomarkers) 

# find all markers of cluster 1
cluster1.markers <- FindMarkers(cells, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

#all clusters
cells.markers <- FindAllMarkers(cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top2markers=cells.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
#sst for delta, ppy for gamma, ghrl for epsilon, cpa1 for acinar, rgs5 for quiescent, pdgfra for activated,
#vwf for endothelial, sds for macrophage, tpsab1 for mast, trac for cytotoxic, sox10 for schwann


##############
#visualization
##############
cells_tsne <- RunTSNE(cells, dims = 1:10)
head(cells_tsne@reductions$tsne@cell.embeddings)
DimPlot(cells_tsne, reduction = "tsne",label = TRUE)

#visualize top2 marker gene
FeaturePlot(cells, features = top2markers$gene[1:4],min.cutoff = 0, max.cutoff = 4)
FeaturePlot(cells, features = top2markers$gene[5:8],min.cutoff = 0, max.cutoff = 4)

 
#Assigning cell type identity to clusters
new.cluster.ids <- c("GCG", "INS", "SST", "PPY", "GHRL", "CPA1", 
                     "KRT19", "RGS5", "PDGFRA", 'VWF', 'SDS', 'TPSAB1', 'TRAC', 'SOX10')
names(new.cluster.ids) <- levels(cells)
cells <- RenameIdents(cells, new.cluster.ids)
DimPlot(cells, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
#plot2<-DimPlot(cells, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()

#bp=DimPlot(cells, reduction = "umap", label = TRUE, pt.size = 0.5)
#bp + legend(name="cell types",labels=c("Alpha", "Beta", "Delta",'Gamma','Epsilon',
#                                  'Acinar','Ductal','Quiescent','Activated','Endothelial',
#                                 'Macrophage','Mast','Cytotoxic','Schwann'))


#heatmap mark genes
library(dplyr)
top1 <- cells.markers %>% group_by(cluster) %>% top_n(n = 1, wt=avg_logFC)
DoHeatmap(cells, features = top1$gene)
