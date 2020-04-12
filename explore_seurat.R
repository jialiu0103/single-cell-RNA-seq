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

#general look at tsne and umap
rdc_cells <- RunTSNE(cells, dims = 1:10)
DimPlot(rdc_cells, reduction = "umap",label = TRUE)
plot2<-DimPlot(rdc_cells, reduction = "tsne",label = TRUE)
CombinePlots(plots = list(plot1, plot2),legend="bottom")

###############
#Identify marker genes for each cluster
###############
