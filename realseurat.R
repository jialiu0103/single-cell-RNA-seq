#load seurat data
cells=readRDS('/projectnb2/bf528/users/group4/project4/code/seurat.rda')

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
DimPlot(cells, reduction = "pca",label = TRUE)

###############
#Identify marker genes for each cluster
###############
# Finding differentially expressed features (cluster biomarkers) 

# find all markers of cluster 1
cluster1.markers <- FindMarkers(cells, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

#all clusters
cells.markers <- FindAllMarkers(cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top1markers=cells.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)
top2markers=cells.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

#sst for delta, ppy for gamma, ghrl for epsilon, cpa1 for acinar, rgs5 for quiescent, pdgfra for activated,
#vwf for endothelial, sds for macrophage, tpsab1 for mast, trac for cytotoxic, sox10 for schwann
write.csv(cells.markers,'realmarker.csv')

##############
#visualization
##############
cells_tsne <- RunTSNE(cells, dims = 1:10)
head(cells_tsne@reductions$tsne@cell.embeddings)
DimPlot(cells_tsne, reduction = "tsne",label = TRUE)

cells_umap <- RunUMAP(cells, dims = 1:10)
head(cells_umap@reductions$umap@cell.embeddings)
DimPlot(cells_umap, reduction = "umap",label = TRUE)

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

realmarker <- read.csv("/projectnb/bf528/users/group4/project4/data/realmarker.csv")
AllMarkers=realmarker
# find novel marker genes
# Find differentially expressed features between feature0 and all other cells, only
# search for positive markers
# Pre-filter features that are detected at <50% frequency in either feature0 and others
feature0.markers <- FindMarkers(cells, ident.1 = "0", ident.2 = NULL, only.pos = TRUE)
feature0=feature0.markers[which(feature0.markers$p_val_adj<0.05),]
#find marker genes already exists
marker0=AllMarkers[which(AllMarkers$cluster==0),]
marker0=marker0$gene
#drop marker genes already exist
f0=feature0[-which(row.names(feature0) %in% marker0),]

feature1.markers <- FindMarkers(cells, ident.1 = "1", ident.2 = NULL, only.pos = TRUE)
feature1=feature1.markers[which(feature1.markers$p_val_adj<0.05),]
marker1=AllMarkers[which(AllMarkers$cluster==1),]
marker1=marker1$gene
f1=feature1[-which(row.names(feature1) %in% marker1),]

feature2.markers <- FindMarkers(cells, ident.1 = "2", ident.2 = NULL, only.pos = TRUE)
feature2=feature2.markers[which(feature2.markers$p_val_adj<0.05),]
marker2=AllMarkers[which(AllMarkers$cluster==2),]
marker2=marker2$gene
f2=feature2[-which(row.names(feature2) %in% marker2),]

feature3.markers <- FindMarkers(cells, ident.1 = "3", ident.2 = NULL, only.pos = TRUE)
feature3=feature3.markers[which(feature3.markers$p_val_adj<0.05),]
marker3=AllMarkers[which(AllMarkers$cluster==3),]
f3=feature3[-which(row.names(feature3) %in% marker3$gene),]

feature4.markers <- FindMarkers(cells, ident.1 = "4", ident.2 = NULL, only.pos = TRUE)
feature4=feature4.markers[which(feature4.markers$p_val_adj<0.05),]
marker4=AllMarkers[which(AllMarkers$cluster==4),]
f4=feature4[-which(row.names(feature4) %in% marker4$gene),]

feature5.markers <- FindMarkers(cells, ident.1 = "5", ident.2 = NULL, only.pos = TRUE)
feature5=feature5.markers[which(feature5.markers$p_val_adj<0.05),]
marker5=AllMarkers[which(AllMarkers$cluster==5),]
f5=feature5[-which(row.names(feature5) %in% marker5$gene),]

feature6.markers <- FindMarkers(cells, ident.1 = "6", ident.2 = NULL, only.pos = TRUE)
feature6=feature6.markers[which(feature6.markers$p_val_adj<0.05),]
marker6=AllMarkers[which(AllMarkers$cluster==6),]
f6=feature6[-which(row.names(feature6) %in% marker6$gene),]

feature.markers <- FindMarkers(cells, ident.1 = "7", ident.2 = NULL, only.pos = TRUE)
feature=feature.markers[which(feature.markers$p_val_adj<0.05),]
marker=AllMarkers[which(AllMarkers$cluster==7),]
f7=feature[-which(row.names(feature) %in% marker$gene),]

feature.markers <- FindMarkers(cells, ident.1 = "8", ident.2 = NULL, only.pos = TRUE)
feature=feature.markers[which(feature.markers$p_val_adj<0.05),]
marker=AllMarkers[which(AllMarkers$cluster==8),]
f8=feature[-which(row.names(feature) %in% marker$gene),]

feature.markers <- FindMarkers(cells, ident.1 = "9", ident.2 = NULL, only.pos = TRUE)
feature=feature.markers[which(feature.markers$p_val_adj<0.05),]
marker=AllMarkers[which(AllMarkers$cluster==9),]
f9=feature[-which(row.names(feature) %in% marker$gene),]

feature.markers <- FindMarkers(cells, ident.1 = "10", ident.2 = NULL, only.pos = TRUE)
feature=feature.markers[which(feature.markers$p_val_adj<0.05),]
marker=AllMarkers[which(AllMarkers$cluster==10),]
f10=feature[-which(row.names(feature) %in% marker$gene),]

feature.markers <- FindMarkers(cells, ident.1 = "11", ident.2 = NULL, only.pos = TRUE)
feature=feature.markers[which(feature.markers$p_val_adj<0.05),]
marker=AllMarkers[which(AllMarkers$cluster==11),]
f11=feature[-which(row.names(feature) %in% marker$gene),]

feature.markers <- FindMarkers(cells, ident.1 = "12", ident.2 = NULL, only.pos = TRUE)
feature=feature.markers[which(feature.markers$p_val_adj<0.05),]
marker=AllMarkers[which(AllMarkers$cluster==12),]
f12=feature[-which(row.names(feature) %in% marker$gene),]

feature.markers <- FindMarkers(cells, ident.1 = "13", ident.2 = NULL, only.pos = TRUE)
feature=feature.markers[which(feature.markers$p_val_adj<0.05),]
marker=AllMarkers[which(AllMarkers$cluster==13),]
f13=feature[-which(row.names(feature) %in% marker$gene),]

feature.markers <- FindMarkers(cells, ident.1 = "14", ident.2 = NULL, only.pos = TRUE)
feature=feature.markers[which(feature.markers$p_val_adj<0.05),]
marker=AllMarkers[which(AllMarkers$cluster==14),]
f14=feature[-which(row.names(feature) %in% marker$gene),]

f0$cluster=rep(0,length(f0$p_val))
f1$cluster=rep(1,length(f1$p_val))
f2$cluster=rep(2,length(f2$p_val))
f3$cluster=rep(3,length(f3$p_val))
f4$cluster=rep(4,length(f4$p_val))
f5$cluster=rep(5,length(f5$p_val))
f6$cluster=rep(6,length(f6$p_val))
f7$cluster=rep(7,length(f7$p_val))
f8$cluster=rep(8,length(f8$p_val))
f9$cluster=rep(9,length(f9$p_val))
f10$cluster=rep(10,length(f10$p_val))
f11$cluster=rep(11,length(f11$p_val))
f12$cluster=rep(12,length(f12$p_val))
f13$cluster=rep(13,length(f13$p_val))
f14$cluster=rep(14,length(f14$p_val))

novelmarker=rbind(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14)
write.csv(novelmarker,'realnovelmarker.csv')
