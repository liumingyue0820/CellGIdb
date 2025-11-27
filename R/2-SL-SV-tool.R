

####    ---> Tool:  Hallmarker Tools <---  ####

# The AddModuleScore function scores the SL and SV network node gene sets
# SL
sl.node.list=SL.propR03.NodeWeight$node
sl.node.list <- sl.node.list[sl.node.list %in% rownames(sec.Malig)]
sl.node.list <- list(as.data.frame(sl.node.list))
sec.Malig <- AddModuleScore(object = sec.Malig, features = sl.node.list, name = 'SL.Node')
# SV
sv.node.list=SV.propR03.NodeWeight$node
sv.node.list <- sv.node.list[sv.node.list %in% rownames(sec.Malig)]
sv.node.list <- list(as.data.frame(sv.node.list))
sec.Malig <- AddModuleScore(object = sec.Malig, features = sv.node.list, name = 'SV.Node')

#  Calculate hallmarker scores
setwd(path)
hallmarker=readLines('h.all.v2024.1.Hs.symbols.gmt')
hallmarker <- strsplit(hallmarker, "\t")
names(hallmarker) <- vapply(hallmarker, function(y) y[1], character(1))
hallmarker <- lapply(hallmarker, "[", -c(1:2))

term=names(hallmarker)
gene_list=hallmarker[[1]]
gene_list <- gene_list[gene_list %in% rownames(sec.Malig)]
gene_list <- list(as.data.frame(gene_list))
sec.Malig <- AddModuleScore(object = sec.Malig, features = gene_list, name = term[1])
for(i in 2:50){
  gene_list=hallmarker[[i]]
  gene_list <- gene_list[gene_list %in% rownames(sec.Malig)]
  gene_list <- list(as.data.frame(gene_list))
  sec.Malig <- AddModuleScore(object = sec.Malig, features = gene_list, name = term[i])
}


# ---> Correlation between SL score and 50 hallmarker scores
a=sec.Malig@meta.data$SL.Node
index=which(colnames(sec.Malig@meta.data)=='SV.Node1')+1 
b=sec.Malig@meta.data[,index]
ind=intersect(which(a!=0), which(b!=0))
aa=a[ind]
bb=b[ind]
test=cor.test(aa,bb,  method = 'spearman')
spearman.result=data.frame(score1='SL.Node',
                           score2=colnames(sec.Malig@meta.data)[index],
                           p.value=test$p.value,
                           corr=test$estimate)
for(j in (index+1):ncol(sec.Malig@meta.data)){  
  b=sec.Malig@meta.data[,j]
  ind=intersect(which(a!=0), which(b!=0))
  if(length(ind)>2){
    aa=a[ind]
    bb=b[ind]
    test=cor.test(aa,bb,  method = 'spearman')
    spearman.result=rbind(spearman.result, data.frame(score1='SL.Node',
                                                      score2=colnames(sec.Malig@meta.data)[j],
                                                      p.value=test$p.value,
                                                      corr=test$estimate))
  }else{
    spearman.result=rbind(spearman.result, data.frame(score1='SL.Node',
                                                      score2=colnames(sec.Malig@meta.data)[j],
                                                      p.value=NA,
                                                      corr=NA))
    
  }
}


# --->  Correlation between SV score and 50 hallmarker scores
a=sec.Malig@meta.data$SV.Node # SV点打分呢
for(j in index:ncol(sec.Malig@meta.data)){
  b=sec.Malig@meta.data[,j]
  ind=intersect(which(a!=0), which(b!=0))
  if(length(ind)>2){
    aa=a[ind]
    bb=b[ind]
    test=cor.test(aa,bb,  method = 'spearman')
    spearman.result=rbind(spearman.result, data.frame(score1='SV.Node',
                                                      score2=colnames(sec.Malig@meta.data)[j],
                                                      p.value=test$p.value,
                                                      corr=test$estimate))
  }else{
    spearman.result=rbind(spearman.result, data.frame(score1='SV.Node',
                                                      score2=colnames(sec.Malig@meta.data)[j],
                                                      p.value=NA,
                                                      corr=NA))
  }
}
spearman.result[which(spearman.result$p.value==0),3]=2.2e-16
spearman.result[which(spearman.result$p.value<2.2e-16),3]=2.2e-16
spearman.result$log10P=-log10(spearman.result$p.value)
# 整理hallmarker名字
hm=paste(tolower(unlist(strsplit(spearman.result[1,2],'_')))[-1],collapse  = '_')
hm=gsub('1','',hm)
hm=gsub('_',' ',hm)
hm=stringr::str_to_title(hm)
for(i in 2:nrow(spearman.result)){
  temp=paste(tolower(unlist(strsplit(spearman.result[i,2],'_')))[-1],collapse  = '_')
  temp=gsub('1','',temp)
  temp=gsub('_',' ',temp)
  temp=stringr::str_to_title(temp)
  hm=c(hm, temp)
}
spearman.result$hallmarker=hm
spearman.result=spearman.result[,c(1,6,4,5)]




####  --->   Tool: Survival Tools   <-----  #####
# In the TCGA cancer dataset, univariate Cox and log-rank tests identified SL (P<0.05, HR<1) as associated with favorable patient prognosis and SV (P<0.05, HR>1) as associated with poor prognosis.

setwd(path)
expre.tcga.luad=read.table('luad_data_RNA_Seq_v2_expression_median.txt', header=T,sep='\t')
survival.pancancer=read.table('GDC-PANCAN.survival.tsv', header = T, sep = '\t')
survival.pancancer$sample=substr(survival.pancancer[,1], 1,15)
survival.pancancer$sample=gsub('-','\\.',survival.pancancer$sample)
dim(survival.pancancer)  #  [1] 18492     4
length(intersect(survival.pancancer$sample,  colnames(expre.tcga.luad)))#497
survival497=survival.pancancer[pmatch(intersect(survival.pancancer$sample,  colnames(expre.tcga.luad)),
                                      survival.pancancer$sample),]
dim(survival497) # [1] 497   4

network.sl=SL.propR03
network.sv=SV.propR03

AllGene=unique(c(network.sl$symbol1, network.sl$symbol2, network.sv$symbol1, network.sv$symbol2))
AllGene <- AllGene[AllGene %in% expre.tcga.luad$Hugo_Symbol]
expre.tcga.luad.slsv=expre.tcga.luad[pmatch(AllGene, expre.tcga.luad$Hugo_Symbol),]
rownames(expre.tcga.luad.slsv)=expre.tcga.luad.slsv$Hugo_Symbol
expre.tcga.luad.slsv=expre.tcga.luad.slsv[,-c(1:2)]
expre.tcga.luad.slsv=expre.tcga.luad.slsv[,intersect(survival.pancancer$sample,  colnames(expre.tcga.luad.slsv))]
dim(expre.tcga.luad.slsv)  # [1] 5253 497

library(survival)
library(survminer) 
#  ---  SL  ---
dim(network.sl)  # [1] 123  23  
#  ---  SV  ---
dim(network.sv)  # [1] 42 23  
#  *****  The function "tcga.slsv.survival" identifies SLSV gene pairs associated with prognosis.    *****
result.survival=tcga.slsv.survival(network.sl, network.sv, expre.tcga.luad.slsv, survival497)
# SL genes related to poor prognosis
SL.prognosis=colnames(result.survival)[grep('SL', colnames(result.survival))]
# SV genes related to a good prognosis
SV.prognosis=colnames(result.survival)[grep('SV', colnames(result.survival))]



####  --->   Tool: Drug Tools      <-----  ##### 
load('Tool5-Lung-drug-SLSV-list.RData')
dim(Drug.sl.list) # [1] 1080246        10
dim(Drug.sv.list) # [1] 189037         10

#   1. Single-cell data were discretized into gene activation/inactivation profiles
Malig.count=sec.Malig@assays$RNA$counts
allGene=c(Drug.sl.list$symbol1,Drug.sl.list$symbol2,Drug.sv.list$symbol1,Drug.sv.list$symbol2)
allGene=unique(allGene)
allGene=intersect(rownames(Malig.count), allGene)
Malig.count.slsv=Malig.count[allGene,]
# . A gene is defined as inactive (respectively, overactive) if its expression level is less (greater) than the 33rd percentile (67th percentile) across samples for each cancer type (to control for cancer type). Otherwise, it is considered to have a normal activation level. 
# Remove genes that are expressed in less than 10 cells
C1=length(which(Malig.count.slsv[1,]!=0))
for(i in 2:nrow(Malig.count.slsv)){C1=c(C1,length(which(Malig.count.slsv[i,]!=0)))}
Malig.count.slsv=Malig.count.slsv[-which(C1<10),]
Malig_act.inact=matrix(0,nrow(Malig.count.slsv),ncol(Malig.count.slsv))
rownames(Malig_act.inact)=rownames(Malig.count.slsv)
colnames(Malig_act.inact)=colnames(Malig.count.slsv)
for(i in 1:nrow(Malig_act.inact)){
  INDEX=which(Malig.count.slsv[i,]!=0)
  expre=Malig.count.slsv[i,INDEX]
  FW=quantile(as.numeric(expre),c(0.33,0.67)) 
  if(FW[1]==FW[2] & FW[2]==1){
    all.expre=as.numeric(Malig.count.slsv[i,])
    ind1=intersect(which(all.expre>=FW[1]), INDEX)
    ind2=setdiff(1:ncol(Malig.count.slsv), ind1)
    Malig_act.inact[i,ind1]=1
    Malig_act.inact[i,ind2]=(-1)
  }
  if(FW[1]!=FW[2]){
    all.expre=as.numeric(Malig.count.slsv[i,])
    ind1=intersect(which(all.expre>=FW[2]), INDEX)
    ind2=intersect(which(all.expre<=FW[1]), INDEX)
    Malig_act.inact[i,ind1]=1
    Malig_act.inact[i,ind2]=(-1)
  }
}

#  2.  Calculate the efficacy score DDS integrating SL and SV based on the single-cell discrete activation-inactivation spectrum    ~~~~~~~
drugID=intersect(unique(Drug.sl.list$drug_name2), unique(Drug.sv.list$drug_name2))
i=1  
sl.list=Drug.sl.list[which(Drug.sl.list$drug_name2==drugID[i]),]
intG=intersect(union(sl.list$symbol1,sl.list$symbol2),rownames(Malig_act.inact))
tempSL.act.inact=Malig_act.inact[pmatch(intG, rownames(Malig_act.inact)),]
Gene1=setdiff(union(sl.list$symbol1,sl.list$symbol2),intG)

sv.list=Drug.sv.list[which(Drug.sv.list$drug_name2==drugID[i]),]
intG=intersect(union(sv.list$symbol1,sv.list$symbol2),rownames(Malig_act.inact))
tempSV.act.inact=Malig_act.inact[pmatch(intG, rownames(Malig_act.inact)),]
Gene1=setdiff(union(sv.list$symbol1,sv.list$symbol2),intG)

# ********  Call the function  ********
result = DDS.alter.score.tool5(sl.list, tempSL.act.inact, sv.list, tempSV.act.inact)
# Add results for the first drug to malignant epithelial cell Seurat object
sec.Malig@meta.data$Tool5.drug1 = result$score
colnames(sec.Malig@meta.data)[dim(sec.Malig@meta.data)[2]] = drugID[i]

# Loop through remaining drugs
for(i in 2:length(drugID)){
  # Extract SL list and activation/inactivation profiles for the i-th drug
  sl.list = Drug.sl.list[which(Drug.sl.list$drug_name2 == drugID[i]), ]
  print(paste(c('Drug ', i, ' has ', nrow(sl.list), ' SL gene pairs'), collapse = ''))
  # Extract SL-related genes from lung single-cell activation/inactivation profile
  intG = intersect(union(sl.list$symbol1, sl.list$symbol2), rownames(Malig_act.inact))
  tempSL.act.inact = Malig_act.inact[pmatch(intG, rownames(Malig_act.inact)), ]
  Gene1 = setdiff(union(sl.list$symbol1, sl.list$symbol2), intG)

  # Extract SV list and activation/inactivation profiles for the i-th drug
  sv.list = Drug.sv.list[which(Drug.sv.list$drug_name2 == drugID[i]), ]
  print(paste(c('Drug ', i, ' has ', nrow(sv.list), ' SV gene pairs'), collapse = ''))
  intG = intersect(union(sv.list$symbol1, sv.list$symbol2), rownames(Malig_act.inact))
  tempSV.act.inact = Malig_act.inact[pmatch(intG, rownames(Malig_act.inact)), ]
  # Extract expressed SV gene pairs list
  Gene1 = setdiff(union(sv.list$symbol1, sv.list$symbol2), intG)

  # ********  Call function DDS.alter.score.tool5  ********
  result = DDS.alter.score.tool5(sl.list, tempSL.act.inact, sv.list, tempSV.act.inact)
  # Add results to malignant epithelial cell Seurat object
  print(identical(rownames(sec.Malig@meta.data), result$cell)) #[1] TRUE
  sec.Malig@meta.data$Tool5.drug1 = result$score
  index = dim(sec.Malig@meta.data)[2] # [1] 3011   59
  colnames(sec.Malig@meta.data)[index] = drugID[i]
  print(paste(c('Drug ', i, ' DDS scoring completed!!!'), collapse = ''))
}

library(Seurat)
FeaturePlot(object = sec.Malig, features = "GDSC1_CGP-60474")



####  --->   Tool: Immunotherapy Tools      <-----  ##### 
dim(network.sl)  # [1] 12762    22  (&&&&&)
dim(network.sv)  # [1] 1588   23    (&&&&&)

###   dataset1: Ravi2023  ---->
load('Immuno.Ravi2023.NSCLC.Expre.Response.RData')
dim(express.57) # [1] 19227    58 (first column is gene symbol)
length(resp.sample) # 42 samples
length(nonresp.sample) # 15 samples

# Extract all unique genes from SL and SV networks
allGene = unique(c(network.sl$symbol1, network.sl$symbol2,
                   network.sv$symbol1, network.sv$symbol2))
length(allGene) 

# Define gene activation/inactivation states using tertiles
# Fisher's test to check if SL gene pair alterations are enriched in response group
# Using Ravi2023 melanoma dataset to demonstrate SL and SV gene pairs
dim(express.57) # 19227    58
intG = intersect(allGene, express.57[, 1]) 
index = pmatch(allGene, express.57[, 1])
express.Gene = express.57[index[!is.na(index)], ]
rownames(express.Gene) = express.Gene$Name
express.Gene = express.Gene[, -1]
dim(express.Gene) 
express.Gene = as.matrix(express.Gene)

# Create gene activation/inactivation profile across all samples
express.act.inact = matrix(0, nrow(express.Gene), ncol(express.Gene))
rownames(express.act.inact) = rownames(express.Gene)
colnames(express.act.inact) = colnames(express.Gene)
for(i in 1:nrow(express.act.inact)){
  # Calculate tertiles
  FW = quantile(express.Gene[i, ], c(0.33, 0.76))
  ind1 = which(express.Gene[i, ] < FW[1]) # indices below lower tertile
  ind2 = which(express.Gene[i, ] > FW[2]) # indices above upper tertile
  express.act.inact[i, ind1] = (-1) # inactivation
  express.act.inact[i, ind2] = 1    # activation
}

length(resp.sample) # 42 samples
length(nonresp.sample) # 15 samples

# Check input data dimensions
dim(express.act.inact) 
dim(network.sl)  
dim(network.sv)  

#  ---->  Call Fisher.Immunotherapy.SLSV function to identify immunotherapy response-related SL/SV pairs
immunotherapy.result = Fisher.Immunotherapy.SLSV(express.act.inact, network.sl, network.sv)
# Expected output:
# [1] "4 immunotherapy response-related SL gene pairs!"
# [1] "2 immunotherapy response-related SV gene pairs!"





