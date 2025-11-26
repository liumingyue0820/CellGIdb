####  Using malignant epithelial cells from the lung cancer dataset as an example, identify cell type-specific SL and SV networks and compute node weights.
install.packages('Seurat')
library(Seurat)



######    Step1: Identify genes differentially overexpressed in malignant epithelial cells    #######

setwd('F:\\cell type specific gi\\202502-202505\\CellGIdb-GitHub\\data')
sec=readRDS('GSE149655_Seurat.rds')
# table(sec@active.ident)
# Malignant epithelial Epithelial cell      T cell               Mast cell           Macrophage 
# 918                  123                  678                  261                  465 
# Endothelial cell     Fibroblast           Plasma cell          B cell 
# 278                  369                  164                   42 
celltype='Malignant epithelial'
diff_gene <- FindMarkers(sec, min.pct = 0.01, 
                         logfc.threshold = 0.25,
                         ident.1 =celltype)
# Genes with log2FC > 0.5 are significantly up-regulated.
diff.up=rownames(diff_gene)[which(diff_gene$avg_log2FC>0.5)]

# Load SL and SV gene pairs
load('SLSV.RData')
dim(SL) #  [1] 93473    15   
sl.gene=unique(c(SL$symbol1,SL$symbol2)) # [1] 13798
dim(SV) #  [1] 14038    15  
sv.gene=unique(c(SV$symbol1,SV$symbol2)) # [1] 6713

# Upregulate the SL gene associated with genes   
up.sl=intersect(diff.up,sl.gene)
index=data.frame(Sym1.Up=SL$symbol1 %in% up.sl,
                 Sym2.Up=SL$symbol2 %in% up.sl)
index$sum=apply(index,1,sum)
SL.upgene=cbind(SL,index)
SL.upgene=SL.upgene[-which(SL.upgene$sum==0),]
SL.upgene$symbol1_symbol2=paste(SL.upgene$symbol1,SL.upgene$symbol2,sep='_')
# Upregulate the SV gene associated with genes
up.sv=intersect(diff.up,sv.gene)
index=data.frame(Sym1.Up=SV$symbol1 %in% up.sv,
                 Sym2.Up=SV$symbol2 %in% up.sv)
index$sum=apply(index,1,sum)
SV.upgene=cbind(SV,index)
SV.upgene=SV.upgene[-which(SV.upgene$sum==0),]
SV.upgene$symbol1_symbol2=paste(SV.upgene$symbol1,SV.upgene$symbol2,sep='_')



####    Step2:    propR identifies co-expressed SL and SV gene pairs        #####
#install.packages('corpcor')
library(corpcor)
#install.packages('ppcor')
library(ppcor)
#install.packages('propr')
library(propr)
method='propR'

#  2、生成特定细胞类型的Seurat对象
sec.Malig=subset(x = sec, idents = 'Malignant epithelial')
slGene=intersect(unique(c(SL.upgene$symbol1,SL.upgene$symbol2)),rownames(sec.Malig@assays$RNA$counts)) 
svGene=intersect(unique(c(SV.upgene$symbol1,SV.upgene$symbol2)),rownames(sec.Malig@assays$RNA$counts)) 

# Input the expression for the count matrix
count.exp=as.matrix(sec.Malig@assays$RNA$counts)
sl.count.exp=count.exp[slGene,]  # SL Gene Expression Matrix
sv.count.exp=count.exp[svGene,]  # SV Gene Expression Matrix
method='propR'

###    1.   SL    ###
# Calculate Propr matrix using correlation coefficient
SL.propR.matrix <- propr(t(sl.count.exp), metric = "cor", ivar = "clr")
index=data.frame(Symind1=match(SL.upgene[,4],rownames(SL.propR.matrix@matrix)),
                 Symind2=match(SL.upgene[,5],colnames(SL.propR.matrix@matrix)))
esti=SL.propR.matrix@matrix[index[1,1],index[1,2]]
for(i in 2:nrow(index)){esti=c(esti,SL.propR.matrix@matrix[index[i,1],index[i,2]])}
SL.propR=cbind(data.frame(esti=esti),SL.upgene)
#Remove NA
SL.propR=SL.propR[-which(is.na(SL.propR$esti)==T),]
SL.propR03=SL.propR[which(abs(SL.propR$esti)>0.3),]

###    2.   SV    ###
# Calculate Propr matrix using correlation coefficient
SV.propR.matrix <- propr(t(sv.count.exp), metric = "cor", ivar = "clr")
index=data.frame(Symind1=match(SV.upgene[,4],rownames(SV.propR.matrix@matrix)),
                 Symind2=match(SV.upgene[,5],colnames(SV.propR.matrix@matrix)))
esti=SV.propR.matrix@matrix[index[1,1],index[1,2]]
for(i in 2:nrow(index)){esti=c(esti,SV.propR.matrix@matrix[index[i,1],index[i,2]])}
SV.propR=cbind(data.frame(esti=esti),SV.upgene)
#RemoveNA
SV.propR=SV.propR[-which(is.na(SV.propR$esti)==T),]
SV.propR03=SV.propR[which(abs(SV.propR$esti)>0.3),]





####    Step3:    Weight of co-expressed network nodes in propR       #####
celltype # "Malignant epithelial"
dataset='GSE149655'  # GSE149655
propR03=SL.propR03
cancer='Lung cancer'

Node.Weight=function(propR03){   ####     Function-------------->   START
  
  node=unique(c(propR03$symbol1,propR03$symbol2))
  print(paste(c('The number of nodes in the network is:',length(node)),collapse = ''))  
  #Generate a random network
  weight.Random1=sample(propR03$esti,nrow(propR03))
  propR03$weight.Random1=weight.Random1
  for(i in 2:1000){
    propR03[,(i+22)]=sample(propR03$esti,nrow(propR03))
  }
  print(dim(propR03)) 
  colnames(propR03)[24:ncol(propR03)]=paste('weight.Random',2:1000,sep = '')
  
  # Calculate the weight of each point according to the formula.
  # top(S)=[w(i)-μR(i)]/σR(i)  
  # w(i): total strength of their local neighbors, represented as w(celltype)(i)
  # μR(i): The mean weight of gene i in a random network.
  # σR(i): Weight variance of gene i in a random network.
  # Subsequently, a random model was constructed to preserve the underlying network topology while uniformly reshufing the edge weights.  
  node.weight=data.frame(node='',neighbour='',degree='',
                         w='',μR='',σR='',topS='') 
  for(i in 1:length(node)){
    index1=which(propR03$symbol1==node[i])
    index2=which(propR03$symbol2==node[i])
    partner=c(propR03[index1,6],propR03[index2,5])
    degree=length(unique(partner))
    estimate=abs(propR03[c(index1,index2),1]) # neighbors edge weight
    w=sum(estimate) # # w(i): total strength of their local neighbors
    # Gene i's weight in a random network
    index=c(index1,index2)
    if(length(index)==1){
      random.weight=as.numeric(propR03[index[1],24:ncol(propR03)])
    }
    if(length(index)>1){
      random.weight=as.numeric(propR03[index[1],24:ncol(propR03)])
      for(m in 2:length(index)){
        random.weight=c(random.weight,
                        as.numeric(propR03[index[m],24:ncol(propR03)]))
      }
    }
    # μR(i):Mean weight of gene i in a random network
    μR=mean(abs(random.weight))
    # σR(i): Weight variance of gene i in a random network
    σR=var(abs(random.weight))
    # Node Weight
    temp=data.frame(node=node[i],neighbour=paste(partner,collapse = '/'),degree=degree,
                    w=w,μR=μR,σR=σR,topS=(w-μR)/σR) 
    node.weight=rbind(node.weight,temp)
  }
  node.weight=node.weight[-1,]
  node.weight=node.weight[order(as.numeric(node.weight$degree), decreasing = T),]
  propR.node.weight=node.weight
  propR.node.weight$celltype=celltype
  propR.node.weight$dataset=paste(c(dataset,cancer),collapse = '-')
  
  #   ** Node weighting results returned **  
  return(propR.node.weight)
  
}
###  Function-------------->  END


# *** Compute the gene node weights for SL and SV ***
type='SL'
SL.propR.Node.Weight=Node.Weight(SL.propR03)
type='SV'
SV.propR.Node.Weight=Node.Weight(SV.propR03)


