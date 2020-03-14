---
title: "ExplorBatchData"
output: pdf_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```


> num_batch=as.matrix(my_batch)
> num_batch[1:3,1:3]
     A1BG A1BG-AS1        A2M
[1,]    0        0 0.05580404
[2,]    0        0 0.06519350
[3,]    0        0 0.08106312
> maxnumbatch=apply(num_batch, 2, max)
> min(maxnumbatch)
[1] 0.01871642
> max(maxnumbatch)
[1] 0.4308522
> minnumbatch=apply(num_batch,2,min)
> min(minnumbatch)
[1] -0.07792684
> max(minnumbatch)
[1] 0

> for (i in minnumbatch) {if(i==0){a = a+1}}
> a
[1] 5544
> a=0
> for (i in num_batch) {if(i==0){a = a+1}}
> a
[1] 25,327,730
> dim(num_batch)
[1] 21252 14917
#about 1/12.5 values in the datatset are negtive



```{r TransToCount, echo=FALSE}
library("readr", lib.loc="/share/pkg.7/r/3.6.0/install/lib64/R/library")
my_batch=read_csv('/restricted/projectnb/camplab/home/jia_liu/human_batch')
hu_batch=t(my_batch)
colnames(hu_batch)=hu_batch[1,]
hu_batch=hu_batch[-1,]
#hu_mb=as.matrix(hu_batch)
#mode(hu_mb) <- 'numeric'
hu_matrix=data.matrix(hu_batch)
mode(hu_matrix) <- 'numeric'
hu_count=round(10000*hu_matrix)
hu_count[hu_count < 0] <- 0
```

> a=0
> for (i in hu_count) {if(i==0){a = a+1}}
> a
[1] 138616646
#here, more than 1/3 of values in dataset are 0

```{r FilterRows, echo=FALSE}
new_count <-hu_count[which(rowSums(hu_count) > 0),]
new_row<-''
new_row=append(new_row,rowSums(new_count))
new_row=new_row[-1]
new_row=as.numeric(new_row)
```

every rowsums are in interval [254,37480121] 
> min(new_row)
[1] 254
> max(new_row)
[1] 37480121

```{r FilterCols, echo=FALSE}
my_count <-new_count[,which(colSums(new_count) > 0)]
my_row<-''
my_row=append(my_row,colSums(my_count))
my_row=my_row[-1]
my_row=as.numeric(my_row)
```

every colsums are in interval [256134,1420936]
> min(my_row)
[1] 256134
> max(my_row)
[1] 1420936

no 0 in colsum and rowsum in dataset
> length(my_row)
[1] 21252
> length(new_row)
[1] 14917
> dim(hu_count)
[1] 14917 21252

```{r SetModel, echo=FALSE}
library("readr", lib.loc="/share/pkg.7/r/3.6.0/install/lib64/R/library")
library("celda", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.6")
my_batch=read_csv('/restricted/projectnb/camplab/home/jia_liu/human_batch')
hu_batch=t(my_batch)
colnames(hu_batch)=hu_batch[1,]
hu_batch=hu_batch[-1,]
hu_matrix=data.matrix(hu_batch)
mode(hu_matrix) <- 'numeric'
hu_count=round(10000*hu_matrix)
hu_count[hu_count < 0] <- 0
model <- celda_CG(counts = hu_count,K = 20, L = 120, verbose = FALSE)
sttsne <- celdaTsne(counts = hu_count, celdaMod = model)
plotDimReduceCluster(dim1 = sttsne[, 1],
                     dim2 = sttsne[, 2],
                     cluster = clusters(model)$z)
```

```{r SetMarkGene, echo=FALSE}
markg1=c('CXCR4','APLNR')
plotDimReduceFeature(dim1 = sttsne[, 1],
                     dim2 = sttsne[, 2],hu_count,features = markg1)

markg2=c('GJA5','ACKR1')
plotDimReduceFeature(dim1 = sttsne[, 1],
                     dim2 = sttsne[, 2],hu_count,features = markg2)
                     
markg3=c('CA4','PROX1')
plotDimReduceFeature(dim1 = sttsne[, 1],
                     dim2 = sttsne[, 2],hu_count,features = markg3)
```

```{r split NEC|TEC, echo=FALSE}
hu_colname=colnames(hu_count)
splitnt=''
for (i in hu_colname) {splitnt=c(splitnt,substr(i,20,23))}
splitnt

spnt=99
for (i in splitnt) {if (i=='TEC'){spnt=c(spnt,1)} else{spnt=c(spnt,0)}}
spnt=spnt[-1]
spnt_count=rbind(spnt,hu_count[1:14916,1:21252])
sp_matrix=data.matrix(spnt_count)
mode(sp_matrix) <- 'numeric'

plotDimReduceFeature(dim1 = sttsne[, 1],
+                      dim2 = sttsne[, 2],sp_matrix,features = 'spnt')
```
