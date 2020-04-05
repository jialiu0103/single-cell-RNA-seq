---
title: "azcount_RandomWalk"
output: html_document
---

```{input az_data}
az_data <- read.delim("/restricted/projectnb/camplab/home/jia_liu/az/scRNA_logCounts?w=")
az_data <- read.delim("/restricted/projectnb/camplab/home/jia_liu/az/scRNA_logCounts?w=", stringsAsFactors=FALSE)
dim(az_data)
```

#get matrix az
```{build matrix}
az_matrix=as.matrix(az_data)
row.names(az_matrix) <- az_matrix[,1]
az=az_matrix[,-1]
mode(az) <- 'numeric'
```

 
#check rowsum & colsum >0, get interval of rowsum & colsum
#dim of rowsum & colsum shows no 0
#get rowsum interval[5.075463,78573.8]
#get colsum interval[501.3506,1049.117]
```{check az}
az_rowsum <-az[which(rowSums(az) > 0),]
dim(az_rowsum)
new_row<-''
new_row=append(new_row,rowSums(az_rowsum))
new_row=new_row[-1]
new_row=as.numeric(new_row)
min(new_row)
max(new_row)


az_colsum <-az[,which(colSums(az) > 0)]
new_col<-''
new_col=append(new_col,colSums(az_colsum))
new_col=new_col[-1]
new_col=as.numeric(new_col)
min(new_col)
max(new_col)
```

#get APOE row, tag of each cell(cname)
```{select APOE col & filtration}
library(NLP)
apoe=az['APOE',]
apoe_matrix=as.matrix(apoe)
apoe_matrix=t(apoe_matrix)


apoe_colname=colnames(az)
cname_str=''
for (i in apoe_colname){cname_str=append(cname_str,as.String(i))}
cname_str=cname_str[-1]

cname=''
for (i in cname_str){cname=append(cname,substr(i,18,24))}
cname=cname[-1]
table(cname)
```


#APOE in all cells
```{get pid}
metedata <- read.delim("/restricted/projectnb/camplab/home/jia_liu/az/scRNA_metadata?w=", stringsAsFactors=FALSE)
pid=metedata$patient
colnames(apoe_matrix)=pid
table(colnames(apoe_matrix))

ad1sum=sum(apoe_matrix[,which(colnames(apoe_matrix)=='AD1')])
ad1_ave=ad1sum/1554

ad2sum=sum(apoe_matrix[,which(colnames(apoe_matrix)=='AD2')])
ad2_ave=ad2sum/1426

ad3sum=sum(apoe_matrix[,which(colnames(apoe_matrix)=='AD3')])
ad3_ave=ad3sum/1243

ad4sum=sum(apoe_matrix[,which(colnames(apoe_matrix)=='AD4')])
ad4_ave=ad4sum/1342

ad5sum=sum(apoe_matrix[,which(colnames(apoe_matrix)=='AD5')])
ad5_ave=ad5sum/484

ad6sum=sum(apoe_matrix[,which(colnames(apoe_matrix)=='AD6')])
ad6_ave=ad6sum/439

ct1sum=sum(apoe_matrix[,which(colnames(apoe_matrix)=='Ct1')])
ct1_ave=ct1sum/750

ct2sum=sum(apoe_matrix[,which(colnames(apoe_matrix)=='Ct2')])
ct2_ave=ct2sum/305

ct3sum=sum(apoe_matrix[,which(colnames(apoe_matrix)=='Ct3')])
ct3_ave=ct3sum/1964

ct4sum=sum(apoe_matrix[,which(colnames(apoe_matrix)=='Ct4')])
ct4_ave=ct4sum/1006

ct5sum=sum(apoe_matrix[,which(colnames(apoe_matrix)=='Ct5')])
ct5_ave=ct5sum/1191

ct6sum=sum(apoe_matrix[,which(colnames(apoe_matrix)=='Ct6')])
ct6_ave=ct6sum/1066
```


#plot average # between samples 
```{visualization}
avg_all=c(ad1_ave,ad2_ave,ad3_ave,ad4_ave,ad5_ave,ad6_ave,ct1_ave,ct2_ave,ct3_ave,ct4_ave,ct5_ave,ct6_ave)
pname=c('ad1_ave','ad2_ave','ad3_ave','ad4_ave','ad5_ave','ad6_ave','ct1_ave','ct2_ave','ct3_ave','ct4_ave','ct5_ave','ct6_ave')
avg_matrix=matrix(avg_all,nrow = 1,byrow = TRUE)
colnames(avg_matrix)=pname
library(RColorBrewer)
colone=brewer.pal(6, "Blues")
coltwo=brewer.pal(6, "Oranges")
color=c(colone,coltwo)
all_bar=barplot(avg_all,names.arg = pname,main = 'APOE in all celltypes',col = color)
```


#APOE in oligo 
```{split}
oligo_id=metedata[,c(1,3,8)]
oligo_id=oligo_id[which(oligo_id$cellType=='oligo'),]
oligo_id$apoe=az['APOE',oligo_id$sampleID]
oligo_apoe=oligo_id[,c(2,4)]
#aa=oligo_apoe[-which(oligo_apoe$apoe==0),]
aa=unstack(oligo_apoe, apoe~patient)
ad1_ave=sum(aa$AD1)/length(aa$AD1)
ad2_ave=sum(aa$AD2)/length(aa$AD2)
ad3_ave=sum(aa$AD3)/length(aa$AD3)
ad4_ave=sum(aa$AD4)/length(aa$AD4)
ad5_ave=sum(aa$AD5)/length(aa$AD5)
ad6_ave=sum(aa$AD6)/length(aa$AD6)
ct1_ave=sum(aa$Ct1)/length(aa$Ct1)
ct2_ave=sum(aa$Ct2)/length(aa$Ct2)
ct3_ave=sum(aa$Ct3)/length(aa$Ct3)
ct4_ave=sum(aa$Ct4)/length(aa$Ct4)
ct5_ave=sum(aa$Ct5)/length(aa$Ct5)
ct6_ave=sum(aa$Ct6)/length(aa$Ct6)

avg_all=c(ad1_ave,ad2_ave,ad3_ave,ad4_ave,ad5_ave,ad6_ave,ct1_ave,ct2_ave,ct3_ave,ct4_ave,ct5_ave,ct6_ave)
#AD=c(ad1_ave,ad2_ave,ad3_ave,ad4_ave,ad5_ave,ad6_ave)
#CT=c(ct1_ave,ct2_ave,ct3_ave,ct4_ave,ct5_ave,ct6_ave)
pname=c('ad1_ave','ad2_ave','ad3_ave','ad4_ave','ad5_ave','ad6_ave','ct1_ave','ct2_ave','ct3_ave','ct4_ave','ct5_ave','ct6_ave')
avg_matrix=matrix(avg_all,nrow = 1,byrow = TRUE)
colnames(avg_matrix)=pname

oligo_bar=barplot(avg_all,names.arg = pname,main = 'APOE in oligo',col = color)

```


#APOE in astro 
```{split}
astro_id=metedata[,c(1,3,8)]
astro_id=astro_id[which(astro_id$cellType=='astro'),]
astro_id$apoe=az['APOE',astro_id$sampleID]
astro_apoe=astro_id[,c(2,4)]
aa=unstack(astro_apoe, apoe~patient)
ad1_ave=sum(aa$AD1)/length(aa$AD1)
ad2_ave=sum(aa$AD2)/length(aa$AD2)
ad3_ave=sum(aa$AD3)/length(aa$AD3)
ad4_ave=sum(aa$AD4)/length(aa$AD4)
ad5_ave=sum(aa$AD5)/length(aa$AD5)
ad6_ave=sum(aa$AD6)/length(aa$AD6)
ct1_ave=sum(aa$Ct1)/length(aa$Ct1)
ct2_ave=sum(aa$Ct2)/length(aa$Ct2)
ct3_ave=sum(aa$Ct3)/length(aa$Ct3)
ct4_ave=sum(aa$Ct4)/length(aa$Ct4)
ct5_ave=sum(aa$Ct5)/length(aa$Ct5)
ct6_ave=sum(aa$Ct6)/length(aa$Ct6)

avg_all=c(ad1_ave,ad2_ave,ad3_ave,ad4_ave,ad5_ave,ad6_ave,ct1_ave,ct2_ave,ct3_ave,ct4_ave,ct5_ave,ct6_ave)
pname=c('ad1_ave','ad2_ave','ad3_ave','ad4_ave','ad5_ave','ad6_ave','ct1_ave','ct2_ave','ct3_ave','ct4_ave','ct5_ave','ct6_ave')
#avg_matrix=matrix(avg_all,nrow = 1,byrow = TRUE)
#colnames(avg_matrix)=pname

astro_bar=barplot(avg_all,names.arg = pname,main = 'APOE in astro',col = color)
```


#APOE in neuron 
```{split}
neuron_id=metedata[,c(1,3,8)]
neuron_id=neuron_id[which(neuron_id$cellType=='neuron'),]
neuron_id$apoe=az['APOE',neuron_id$sampleID]
neuron_apoe=neuron_id[,c(2,4)]
aa=unstack(neuron_apoe, apoe~patient)
ad1_ave=sum(aa$AD1)/length(aa$AD1)
ad2_ave=sum(aa$AD2)/length(aa$AD2)
ad3_ave=sum(aa$AD3)/length(aa$AD3)
ad4_ave=sum(aa$AD4)/length(aa$AD4)
ad5_ave=sum(aa$AD5)/length(aa$AD5)
ad6_ave=sum(aa$AD6)/length(aa$AD6)
ct1_ave=sum(aa$Ct1)/length(aa$Ct1)
ct2_ave=sum(aa$Ct2)/length(aa$Ct2)
ct3_ave=sum(aa$Ct3)/length(aa$Ct3)
ct4_ave=sum(aa$Ct4)/length(aa$Ct4)
ct5_ave=sum(aa$Ct5)/length(aa$Ct5)
ct6_ave=sum(aa$Ct6)/length(aa$Ct6)

avg_all=c(ad1_ave,ad2_ave,ad3_ave,ad4_ave,ad5_ave,ad6_ave,ct1_ave,ct2_ave,ct3_ave,ct4_ave,ct5_ave,ct6_ave)
pname=c('ad1_ave','ad2_ave','ad3_ave','ad4_ave','ad5_ave','ad6_ave','ct1_ave','ct2_ave','ct3_ave','ct4_ave','ct5_ave','ct6_ave')
#avg_matrix=matrix(avg_all,nrow = 1,byrow = TRUE)
#colnames(avg_matrix)=pname

neuron_bar=barplot(avg_all,names.arg = pname,main = 'APOE in neuron',col = color)
```


#APOE in mg 
```{split}
mg_id=metedata[,c(1,3,8)]
mg_id=mg_id[which(mg_id$cellType=='mg'),]
mg_id$apoe=az['APOE',mg_id$sampleID]
mg_apoe=mg_id[,c(2,4)]
aa=unstack(mg_apoe, apoe~patient)
ad1_ave=sum(aa$AD1)/length(aa$AD1)
ad2_ave=sum(aa$AD2)/length(aa$AD2)
ad3_ave=sum(aa$AD3)/length(aa$AD3)
ad4_ave=sum(aa$AD4)/length(aa$AD4)
ad5_ave=sum(aa$AD5)/length(aa$AD5)
ad6_ave=sum(aa$AD6)/length(aa$AD6)
ct1_ave=sum(aa$Ct1)/length(aa$Ct1)
ct2_ave=sum(aa$Ct2)/length(aa$Ct2)
ct3_ave=sum(aa$Ct3)/length(aa$Ct3)
ct4_ave=sum(aa$Ct4)/length(aa$Ct4)
ct5_ave=sum(aa$Ct5)/length(aa$Ct5)
ct6_ave=sum(aa$Ct6)/length(aa$Ct6)

avg_all=c(ad1_ave,ad2_ave,ad3_ave,ad4_ave,ad5_ave,ad6_ave,ct1_ave,ct2_ave,ct3_ave,ct4_ave,ct5_ave,ct6_ave)
pname=c('ad1_ave','ad2_ave','ad3_ave','ad4_ave','ad5_ave','ad6_ave','ct1_ave','ct2_ave','ct3_ave','ct4_ave','ct5_ave','ct6_ave')
#avg_matrix=matrix(avg_all,nrow = 1,byrow = TRUE)
#colnames(avg_matrix)=pname

mg_bar=barplot(avg_all,names.arg = pname,main = 'APOE in mg',col = color)
```


#APOE in endo 
```{split}
endo_id=metedata[,c(1,3,8)]
endo_id=endo_id[which(endo_id$cellType=='endo'),]
endo_id$apoe=az['APOE',endo_id$sampleID]
endo_apoe=endo_id[,c(2,4)]
aa=unstack(endo_apoe, apoe~patient)
ad1_ave=sum(aa$AD1)/length(aa$AD1)
ad2_ave=sum(aa$AD2)/length(aa$AD2)
ad3_ave=sum(aa$AD3)/length(aa$AD3)
ad4_ave=sum(aa$AD4)/length(aa$AD4)
ad5_ave=sum(aa$AD5)/length(aa$AD5)
ad6_ave=sum(aa$AD6)/length(aa$AD6)
ct1_ave=sum(aa$Ct1)/length(aa$Ct1)
ct2_ave=sum(aa$Ct2)/length(aa$Ct2)
ct3_ave=sum(aa$Ct3)/length(aa$Ct3)
ct4_ave=sum(aa$Ct4)/length(aa$Ct4)
ct5_ave=sum(aa$Ct5)/length(aa$Ct5)
ct6_ave=sum(aa$Ct6)/length(aa$Ct6)

avg_all=c(ad1_ave,ad2_ave,ad3_ave,ad4_ave,ad5_ave,ad6_ave,ct1_ave,ct2_ave,ct3_ave,ct4_ave,ct5_ave,ct6_ave)
pname=c('ad1_ave','ad2_ave','ad3_ave','ad4_ave','ad5_ave','ad6_ave','ct1_ave','ct2_ave','ct3_ave','ct4_ave','ct5_ave','ct6_ave')
#avg_matrix=matrix(avg_all,nrow = 1,byrow = TRUE)
#colnames(avg_matrix)=pname

endo_bar=barplot(avg_all,names.arg = pname,main = 'APOE in endo',col = color)
```


#APOE in OPC 
```{split}
OPC_id=metedata[,c(1,3,8)]
OPC_id=OPC_id[which(OPC_id$cellType=='OPC'),]
OPC_id$apoe=az['APOE',OPC_id$sampleID]
OPC_apoe=OPC_id[,c(2,4)]
aa=unstack(OPC_apoe, apoe~patient)
ad1_ave=sum(aa$AD1)/length(aa$AD1)
ad2_ave=sum(aa$AD2)/length(aa$AD2)
ad3_ave=sum(aa$AD3)/length(aa$AD3)
ad4_ave=sum(aa$AD4)/length(aa$AD4)
ad5_ave=sum(aa$AD5)/length(aa$AD5)
ad6_ave=sum(aa$AD6)/length(aa$AD6)
ct1_ave=sum(aa$Ct1)/length(aa$Ct1)
ct2_ave=sum(aa$Ct2)/length(aa$Ct2)
ct3_ave=sum(aa$Ct3)/length(aa$Ct3)
ct4_ave=sum(aa$Ct4)/length(aa$Ct4)
ct5_ave=sum(aa$Ct5)/length(aa$Ct5)
ct6_ave=sum(aa$Ct6)/length(aa$Ct6)

avg_all=c(ad1_ave,ad2_ave,ad3_ave,ad4_ave,ad5_ave,ad6_ave,ct1_ave,ct2_ave,ct3_ave,ct4_ave,ct5_ave,ct6_ave)
pname=c('ad1_ave','ad2_ave','ad3_ave','ad4_ave','ad5_ave','ad6_ave','ct1_ave','ct2_ave','ct3_ave','ct4_ave','ct5_ave','ct6_ave')
#avg_matrix=matrix(avg_all,nrow = 1,byrow = TRUE)
#colnames(avg_matrix)=pname

OPC_bar=barplot(avg_all,names.arg = pname,main = 'APOE in OPC',col = color)
```


#APOE in unID 
```{split}
unID_id=metedata[,c(1,3,8)]
unID_id=unID_id[which(unID_id$cellType=='unID'),]
unID_id$apoe=az['APOE',unID_id$sampleID]
unID_apoe=unID_id[,c(2,4)]
aa=unstack(unID_apoe, apoe~patient)
ad1_ave=sum(aa$AD1)/length(aa$AD1)
ad2_ave=sum(aa$AD2)/length(aa$AD2)
ad3_ave=sum(aa$AD3)/length(aa$AD3)
ad4_ave=sum(aa$AD4)/length(aa$AD4)
ad5_ave=sum(aa$AD5)/length(aa$AD5)
ad6_ave=sum(aa$AD6)/length(aa$AD6)
ct1_ave=sum(aa$Ct1)/length(aa$Ct1)
ct2_ave=sum(aa$Ct2)/length(aa$Ct2)
ct3_ave=sum(aa$Ct3)/length(aa$Ct3)
ct4_ave=sum(aa$Ct4)/length(aa$Ct4)
ct5_ave=sum(aa$Ct5)/length(aa$Ct5)
ct6_ave=sum(aa$Ct6)/length(aa$Ct6)

avg_all=c(ad1_ave,ad2_ave,ad3_ave,ad4_ave,ad5_ave,ad6_ave,ct1_ave,ct2_ave,ct3_ave,ct4_ave,ct5_ave,ct6_ave)
pname=c('ad1_ave','ad2_ave','ad3_ave','ad4_ave','ad5_ave','ad6_ave','ct1_ave','ct2_ave','ct3_ave','ct4_ave','ct5_ave','ct6_ave')
#avg_matrix=matrix(avg_all,nrow = 1,byrow = TRUE)
#colnames(avg_matrix)=pname

unID_bar=barplot(avg_all,names.arg = pname,main = 'APOE in unIDENTIFICATION',col = color)
```
