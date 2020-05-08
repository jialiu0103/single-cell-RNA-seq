###############input data from salmon alevin
files <- "/projectnb/bf528/project_4_scrnaseq/GSM2230760__salmon_quant/alevin/quants_mat.gz"
txi <- tximport(files, type = "alevin")
count=txi$counts
count=as.matrix(count)
##############filter our rowsum colsum == 0
new_count <-count[which(rowSums(count) > 0),]
my_count <-new_count[,which(colSums(new_count) > 0)]
#get interval of sum
my_col<-''
my_col=append(my_col,colSums(my_count))
my_col=my_col[-1]
my_col=as.numeric(my_col)
summary(my_col)

my_row<-''
my_row=append(my_row,rowSums(my_count))
my_row=my_row[-1]
my_row=as.numeric(my_row)
summary(my_row)

row_mid=my_row[my_row>5.5]
row_mid=row_mid[row_mid<377] #length=15784

col_mid=my_col[my_col>32]
col_mid=col_mid[col_mid<179]

##################### use seurat to control quanlity
obj=CreateSeuratObject(my_count)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
#obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
obj <- subset(obj, subset = nFeature_RNA > 5 & nFeature_RNA < 5000)
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)


obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(obj), 10)

# plot variable features with and without labels
plot1=VariableFeaturePlot(obj)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#normalize gene expression in each cells
all.genes <- rownames(obj)
obj <- ScaleData(obj, features = all.genes)

#run pca
obj <- RunPCA(obj, features = VariableFeatures(object = obj))

# Examine and visualize PCA results a few different ways
print(obj[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(obj, dims = 1:2, reduction = "pca")
DimHeatmap(obj, dims = 2, cells = 500, balanced = TRUE)
DimHeatmap(obj, dims = 1:15, cells = 500, balanced = TRUE)

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
obj <- JackStraw(obj, num.replicate = 100)
obj <- ScoreJackStraw(obj, dims = 1:20)
#JackStrawPlot(obj, dims = 1:15)
#ElbowPlot(obj)
