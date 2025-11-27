####  --->   Tool4: Survival Tools   <-----  #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ function ~~~~~START
##  Identification of SLSVs Associated with Patient Survival
##  Input: Cell-type-specific SL and SV within a dataset
##  Output: SL pairs associated with favorable patient survival and SV pairs associated with poor survival
# If the co-expression correlation coefficient > 0, G1 and G2 are activated simultaneously (1)
# If the co-expression correlation coefficient < 0, one of G1 and G2 is activated (1) while the other is deactivated (-1)

library(survival)
library(survminer)

tcga.slsv.survival=function(network.sl, network.sv, tcga.expre, sur.dat){
  
  #  ******  SL   ******
  # --- 1.1 --- Constructing a Gene Activation-Deactivation Profile
  Gene=unique(c(network.sl$symbol1, network.sl$symbol2))
  Gene <- Gene[Gene %in% rownames(tcga.expre)]
  # Extract the gene expression profile
  expre1=tcga.expre[Gene, ]
  expre_act.inact=matrix(0,nrow(expre1),ncol(expre1))
  rownames(expre_act.inact)=rownames(expre1)
  colnames(expre_act.inact)=colnames(expre1)
  for(i in 1:nrow(expre1)){
    expre=expre1[i,which(expre1[i,]!=0)]
    expre.value=as.numeric(expre[1,])
    FW=quantile(expre.value,c(0.33,0.67)) 
    index1=which(expre.value>FW[2])
    index2=which(expre.value<FW[1])
    expre_act.inact[i,index1]=1
    expre_act.inact[i,index2]=(-1)
  }
  
  #  --- 1.2 ---  Gene pair alteration profile, row gene pairs, list samples 
  #  ******   The correlation coefficient of gene co-expression is greater than 0, and G1 and G2 are activated simultaneously (1)     *******
  coact.list=network.sl[which(network.sl$esti>0),]
  coact.alter=matrix(0, nrow(coact.list), ncol(expre_act.inact))
  rownames(coact.alter)=coact.list$symbol1.symbol2 #The line name is the gene pair
  colnames(coact.alter)=colnames(expre_act.inact) #The column name is the cell ID
  for(m in 1:nrow(coact.list)){  
    genepair_act.inact=expre_act.inact[pmatch(coact.list[m,5:6], rownames(expre_act.inact)),]#
    s1=apply(genepair_act.inact,2,sum)
    coact.alter[m,which(s1==2)]=1  # ****  Simultaneously activating the gene pair (plus =2) changes the state to 1  ****
  }
  
  #  ******   The correlation coefficient between genes and co-expression is less than 0, with one being activated (1) and the other inactivated (-1).    *******
  act.inact.list=network.sl[which(network.sl$esti<0),]
  act.inact.list=act.inact.list[which(act.inact.list$sum==1),]
  act.inact.alter=matrix(0, nrow(act.inact.list), ncol(expre_act.inact))
  rownames(act.inact.alter)=act.inact.list$symbol1.symbol2 #The line name is the gene pair
  colnames(act.inact.alter)=colnames(expre_act.inact) #The column name is the cell ID
  for(m in 1:nrow(act.inact.list)){  
    genepair_act.inact=expre_act.inact[pmatch(act.inact.list[m,5:6], rownames(expre_act.inact)),]#
    if(act.inact.list[m,17]==TRUE){ 
      index1=intersect(which(genepair_act.inact[1,]==1),which(genepair_act.inact[2,]==(-1)))
      act.inact.alter[m,index1]=1  # **** The altered state of cells that simultaneously meet the conditions of gene1 up-regulation and gene2 down-regulation is 1    ****
    }
    if(act.inact.list[m,18]==TRUE){  
      index1=intersect(which(genepair_act.inact[1,]==(-1)),which(genepair_act.inact[2,]==1))
      act.inact.alter[m,index1]=1  # ****  The altered state of cells that simultaneously meet the conditions of gene1 up-regulation and gene2 down-regulation is 1   ****
    }
    
  }
  
  # Gene pair alter profile
  alteration=rbind(coact.alter, act.inact.alter)
  # Extract SL gene pairs that have changed in at least three samples
  s1=apply(alteration,1,sum)
  alteration=alteration[which(s1>=3),]

  #  ----  1.3 ----   log-rank test identifies gene pairs related to prognosis
  df=cbind(sur.dat, t(alteration))
  pfilter <- 0.05    
  log.rank.result <- data.frame()  
  for(i in 5:ncol(df)){   
    # The difference in survival between altered and unaltered genes
    surv_diff <- survdiff(Surv(OS.time, OS) ~ df[,i], data = df)
    pvalue <- (1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1))
    # Extract the median survival time
    fit <- survfit(Surv(OS.time, OS) ~ df[,i], data = df)
    med_time=median(fit)
    
    if(pvalue<pfilter){ 
      log.rank.result <- rbind(log.rank.result,
                               cbind(gene_pair=colnames(df)[i],
                                     log.rank.pvalue= pvalue ,
                                     time_0=med_time[1,1],
                                     time_1=med_time[2,1]
                               ))
    }
  }   
  
  # Significant SL gene pairs were extracted, and the survival time of patients in the change group was longer
  ind=which(is.na(log.rank.result$time_1)==T)
  log.rank.result$time_1.new=log.rank.result$time_1
  log.rank.result[ind,5]=as.numeric(log.rank.result[ind,3])+1
  log.rank.result=log.rank.result[which(as.numeric(log.rank.result$time_1)>as.numeric(log.rank.result$time_0)),1:4]
  print(paste(c('SL genes related to poor prognosis:',dim(log.rank.result)[1],'个！！！'),collapse =''))
  if(dim(log.rank.result)[1]>0){
    survival.sl.df=df[,c(1:4,pmatch(log.rank.result$gene_pair, colnames(df)))]
    colnames(survival.sl.df)[5:ncol(survival.sl.df)]=gsub('/',':',colnames(survival.sl.df)[5:ncol(survival.sl.df)])
    colnames(survival.sl.df)[5:ncol(survival.sl.df)]=paste('SL',colnames(survival.sl.df)[5:ncol(survival.sl.df)], sep = '_')
  }else{
    survival.sl.df=df[,1:4]
  }
  
  
  #  ******  SV   ******
  # --- 1.1 --- Constructing a Gene Activation-Deactivation Profile
  Gene=unique(c(network.sv$symbol1, network.sv$symbol2))
  Gene <- Gene[Gene %in% rownames(tcga.expre)]
  # Extract the gene expression profile
  expre1=tcga.expre[Gene, ]
  expre_act.inact=matrix(0,nrow(expre1),ncol(expre1))
  rownames(expre_act.inact)=rownames(expre1)
  colnames(expre_act.inact)=colnames(expre1)
  for(i in 1:nrow(expre1)){
    expre=expre1[i,which(expre1[i,]!=0)]
    expre.value=as.numeric(expre[1,])
    FW=quantile(expre.value,c(0.33,0.67)) 
    index1=which(expre.value>FW[2])
    index2=which(expre.value<FW[1])
    expre_act.inact[i,index1]=1
    expre_act.inact[i,index2]=(-1)
  }
  
  #  --- 1.2 ---   Gene pair alteration profile, row gene pairs, list samples
  #  ******   The correlation coefficient of gene co-expression is greater than 0, and G1 and G2 are activated simultaneously (1)     *******
  coact.list=network.sv[which(network.sv$esti>0),]
  coact.alter=matrix(0, nrow(coact.list), ncol(expre_act.inact))
  rownames(coact.alter)=coact.list$symbol1.symbol2 #The line name is the gene pair
  colnames(coact.alter)=colnames(expre_act.inact) #The column name is the cell ID
  for(m in 1:nrow(coact.list)){   
    genepair_act.inact=expre_act.inact[pmatch(coact.list[m,5:6], rownames(expre_act.inact)),]#
    s1=apply(genepair_act.inact,2,sum)
    coact.alter[m,which(s1==2)]=1  
  }
  
  #  ******   The correlation coefficient between genes and co-expression is less than 0, with one being activated (1) and the other inactivated (-1).    *******
  act.inact.list=network.sv[which(network.sv$esti<0),]
  act.inact.list=act.inact.list[which(act.inact.list$sum==1),]
  act.inact.alter=matrix(0, nrow(act.inact.list), ncol(expre_act.inact))
  rownames(act.inact.alter)=act.inact.list$symbol1.symbol2
  colnames(act.inact.alter)=colnames(expre_act.inact) 
  for(m in 1:nrow(act.inact.list)){   
    genepair_act.inact=expre_act.inact[pmatch(act.inact.list[m,5:6], rownames(expre_act.inact)),]#
    if(act.inact.list[m,17]==TRUE){
      index1=intersect(which(genepair_act.inact[1,]==1),which(genepair_act.inact[2,]==(-1)))
      act.inact.alter[m,index1]=1  # **** The altered state of cells that simultaneously meet the conditions of gene1 up-regulation and gene2 down-regulation is 1    ****
    }
    if(act.inact.list[m,18]==TRUE){ 
      index1=intersect(which(genepair_act.inact[1,]==(-1)),which(genepair_act.inact[2,]==1))
      act.inact.alter[m,index1]=1  # **** The altered state of cells that simultaneously meet the conditions of gene1 up-regulation and gene2 down-regulation is 1    ****
    }
    
  }
  
  # Gene pair alter profile
  alteration=rbind(coact.alter, act.inact.alter)
  #  Extract SL gene pairs that have changed in at least three samples
  s1=apply(alteration,1,sum)
  alteration=alteration[which(s1>=3),]

  #  ----  1.3 ----   log-rank test identifies gene pairs related to prognosis
  df=cbind(sur.dat, t(alteration))
  pfilter <- 0.05    
  log.rank.result <- data.frame() 
  for(i in 5:ncol(df)){   
    # The difference in survival between altered and unaltered genes
    surv_diff <- survdiff(Surv(OS.time, OS) ~ df[,i], data = df)
    pvalue <- (1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1))
    #  Extract the median survival time
    fit <- survfit(Surv(OS.time, OS) ~ df[,i], data = df)
    med_time=median(fit)
    
    if(pvalue<pfilter){ 
      log.rank.result <- rbind(log.rank.result,
                               cbind(gene_pair=colnames(df)[i],
                                     log.rank.pvalue= pvalue ,
                                     time_0=med_time[1,1],
                                     time_1=med_time[2,1]
                               ))
    }
  }   
  
  # Significant SL gene pairs were extracted, and the survival time of patients in the change group was longer
  ind=which(is.na(log.rank.result$time_1)==T)
  log.rank.result$time_1.new=log.rank.result$time_1
  log.rank.result[ind,5]=as.numeric(log.rank.result[ind,3])+1
  log.rank.result=log.rank.result[which(as.numeric(log.rank.result$time_1)<as.numeric(log.rank.result$time_0)),1:4]
  print(paste(c('SV genes related to good prognosis:',dim(log.rank.result)[1],'！！！'),collapse =''))
  if(dim(log.rank.result)[1]>0){
    survival.sv.df=df[,c(1:4,pmatch(log.rank.result$gene_pair, colnames(df)))]
    colnames(survival.sv.df)[5:ncol(survival.sv.df)]=gsub('/',':',colnames(survival.sv.df)[5:ncol(survival.sv.df)])
    colnames(survival.sv.df)[5:ncol(survival.sv.df)]=paste('SV',colnames(survival.sv.df)[5:ncol(survival.sv.df)], sep = '_')
    
  }else{
    survival.sv.df=df[,1:6]
  }
  
  
  #  ----->  Return the altered profile of SL and SV genes related to prognosis
  result.survival=cbind(survival.sl.df, survival.sv.df[,5:ncol(survival.sv.df)])
  return(result.survival)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ function ~~~~~END







####          --->      Tool5: DISCERN Tools      <-----  ##### 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ function ~~~~~START
######        Define function to calculate DDS integrated score: 
######        (SL alterations + reversed SV alterations)       
# Core concept: Calculate the proportion of genetic interactions (GI) that occurred,
# (SL alterations + reversed SV alterations) / N (total number of drug-related SL and SV pairs)
# Input: SL gene pair list, single-cell SL gene activation/inactivation profile,
#        SV gene pair list, single-cell SV gene activation/inactivation profile.
# Output: Score for altered gene pairs in all cells: alter(SL+!SV)/total(SL+SV)


DDS.alter.score.tool5 = function(sl.list, tempSL.act.inact, sv.list, tempSV.act.inact){
  
  #  ******  Process SL gene pair alteration profiles (considering three alteration types)    *******
  
  #   1.  Co-inactivated gene pairs
  sl.coinact.alter = matrix(0, 1, ncol(tempSL.act.inact))
  if(length(which(sl.list$group == 'co.inact')) > 0){
    sl.coinact.list = sl.list[which(sl.list$group == 'co.inact'), ]
    # Initialize co-inactivation alteration matrix for all cells
    sl.coinact.alter = matrix(0, nrow(sl.coinact.list), ncol(tempSL.act.inact))
    rownames(sl.coinact.alter) = apply(sl.coinact.list[, 1:2], 1, paste, collapse = ';') # Row names: co-inact SL genes
    colnames(sl.coinact.alter) = colnames(tempSL.act.inact) # Column names: cell IDs
    for(m in 1:nrow(sl.coinact.list)){   # m: loop through gene pairs
      genepair_act.inact = tempSL.act.inact[pmatch(sl.coinact.list[m, 1:2], rownames(tempSL.act.inact)), ]
      dim(genepair_act.inact) # [1] 2 12842 - alteration profile of the m-th gene pair in all cells
      s1 = apply(genepair_act.inact, 2, sum)
      sl.coinact.alter[m, which(s1 == (-2))] = 1  # **** Simultaneous inactivation of SL pair (sum = -2) -> alteration status = 1 ****
    }
  }
  
  #   2.  G1Act.G2Inact gene pairs
  sl.G1ActG2Inact.alter = matrix(0, 1, ncol(tempSL.act.inact))
  if(length(which(sl.list$group == 'G1Act.G2Inact')) > 0){
    sl.G1ActG2Inact.list = sl.list[which(sl.list$group == 'G1Act.G2Inact'), ]
    # Initialize alteration matrix for all cells
    sl.G1ActG2Inact.alter = matrix(0, nrow(sl.G1ActG2Inact.list), ncol(tempSL.act.inact))
    rownames(sl.G1ActG2Inact.alter) = apply(sl.G1ActG2Inact.list[, 1:2], 1, paste, collapse = ';') # Row names: co-inact SL genes
    colnames(sl.G1ActG2Inact.alter) = colnames(tempSL.act.inact) # Column names: cell IDs
    for(m in 1:nrow(sl.G1ActG2Inact.list)){   # m: loop through gene pairs
      genepair_act.inact = tempSL.act.inact[pmatch(sl.G1ActG2Inact.list[m, 1:2], rownames(tempSL.act.inact)), ]
      dim(genepair_act.inact) # [1] 2 12842 - alteration profile of the m-th gene pair in all cells
      ind1 = which(genepair_act.inact[1, ] == 1) # G1 activated
      ind2 = which(genepair_act.inact[2, ] == (-1)) # G2 inactivated
      int.index = intersect(ind1, ind2)
      sl.G1ActG2Inact.alter[m, int.index] = 1  # **** G1 activation + G2 inactivation -> alteration status = 1 ****
    }
  }
  
  #   3.  G1Inact.G2Act gene pairs
  sl.G1Inact.G2Act.alter = matrix(0, 1, ncol(tempSL.act.inact))
  if(length(which(sl.list$group == 'G1Inact.G2Act')) > 0){
    sl.G1Inact.G2Act.list = sl.list[which(sl.list$group == 'G1Inact.G2Act'), ]
    # Initialize alteration matrix for all cells
    sl.G1Inact.G2Act.alter = matrix(0, nrow(sl.G1Inact.G2Act.list), ncol(tempSL.act.inact))
    rownames(sl.G1Inact.G2Act.alter) = apply(sl.G1Inact.G2Act.list[, 1:2], 1, paste, collapse = ';') # Row names: co-inact SL genes
    colnames(sl.G1Inact.G2Act.alter) = colnames(tempSL.act.inact) # Column names: cell IDs
    for(m in 1:nrow(sl.G1Inact.G2Act.list)){   # m: loop through gene pairs
      genepair_act.inact = tempSL.act.inact[pmatch(sl.G1Inact.G2Act.list[m, 1:2], rownames(tempSL.act.inact)), ]
      ind1 = which(genepair_act.inact[1, ] == (-1)) # G1 inactivated
      ind2 = which(genepair_act.inact[2, ] == 1) # G2 activated
      int.index = intersect(ind1, ind2)
      sl.G1Inact.G2Act.alter[m, int.index] = 1  # **** G1 inactivation + G2 activation -> alteration status = 1 ****
    }
  }
  
  #   4.  Summarize SL gene pair alteration counts
  sl.alter = rbind(sl.coinact.alter, sl.G1ActG2Inact.alter, sl.G1Inact.G2Act.alter)
  # Number of altered SL pairs per cell
  sl.alter.num = apply(sl.alter, 2, sum)
  # Integrated altered SL gene pairs per cell
  sl.alter.gp = rep(NA, length(sl.alter.num))
  for(n in 1:ncol(sl.alter)){sl.alter.gp[n] = paste(rownames(sl.alter)[which(sl.alter[, n] == 1)], collapse = '//')}
  
  
  #  ******  Process SV gene pair alteration profiles    *******
    #   1. Co-inactivated gene pairs (considering three alteration types)
  sv.coinact.alter = matrix(0, 1, ncol(tempSV.act.inact))
  if(length(which(sv.list$group == 'co.inact')) > 0){
    #   1. Co-inactivated gene pairs
    sv.coinact.list = sv.list[which(sv.list$group == 'co.inact'), ]
    # Initialize co-inactivation alteration matrix for all cells
    sv.coinact.alter = matrix(0, nrow(sv.coinact.list), ncol(tempSV.act.inact))
    rownames(sv.coinact.alter) = apply(sv.coinact.list[, 1:2], 1, paste, collapse = ';') # Row names: co-inact SV genes
    colnames(sv.coinact.alter) = colnames(tempSV.act.inact) # Column names: cell IDs
    for(m in 1:nrow(sv.coinact.list)){   # m: loop through gene pairs
      genepair_act.inact = tempSV.act.inact[pmatch(sv.coinact.list[m, 1:2], rownames(tempSV.act.inact)), ]
      dim(genepair_act.inact) # [1] 2 12842 - alteration profile of the m-th gene pair in all cells
      s1 = apply(genepair_act.inact, 2, sum)
      # identical(names(s1), colnames(coinact.alter)) [1] TRUE
      # **** Reverse co-inactivated SV pairs: sum = -2 is SV, sum = -1 (only one gene altered) reverses SV -> alteration status = 1, reference DU, DD SR model. ****
      sv.coinact.alter[m, which(s1 == (-1))] = 1  
    }
  }
  
  #   2. G1Act.G2Inact gene pairs
  sv.G1ActG2Inact.alter = matrix(0, 1, ncol(tempSV.act.inact))
  if(length(which(sv.list$group == 'G1Act.G2Inact')) > 0){
    sv.G1ActG2Inact.list = sv.list[which(sv.list$group == 'G1Act.G2Inact'), ]
    # Initialize alteration matrix for all cells
    sv.G1ActG2Inact.alter = matrix(0, nrow(sv.G1ActG2Inact.list), ncol(tempSV.act.inact))
    rownames(sv.G1ActG2Inact.alter) = apply(sv.G1ActG2Inact.list[, 1:2], 1, paste, collapse = ';') # Row names: co-inact SV genes
    colnames(sv.G1ActG2Inact.alter) = colnames(tempSV.act.inact) # Column names: cell IDs
    for(m in 1:nrow(sv.G1ActG2Inact.list)){   # m: loop through gene pairs
      genepair_act.inact = tempSV.act.inact[pmatch(sv.G1ActG2Inact.list[m, 1:2], rownames(tempSV.act.inact)), ]
      # **** Reverse G1 activation + G2 inactivation SV pairs, reference DU, DD SR model, reversed state: G1 inactivation ****
      ind1 = which(genepair_act.inact[1, ] == (-1)) # G1 inactivated
      ind2 = which(genepair_act.inact[2, ] == (-1)) # G2 inactivated
      int.index = intersect(ind1, ind2)
      #print(paste(c('Gene pair ', m, ' of SV: ', length(int.index)), collapse = ''))
      # identical(names(s1), colnames(coinact.alter)) [1] TRUE
      sv.G1ActG2Inact.alter[m, int.index] = 1  # **** G1 inactivation + G2 inactivation -> alteration status = 1 ****
    }
  }
  
  #   3. G1Inact.G2Act gene pairs
  sv.G1Inact.G2Act.alter = matrix(0, 1, ncol(tempSV.act.inact))
  if(length(which(sv.list$group == 'G1Inact.G2Act')) > 0){
    sv.G1Inact.G2Act.list = sv.list[which(sv.list$group == 'G1Inact.G2Act'), ]
    # Initialize alteration matrix for all cells
    sv.G1Inact.G2Act.alter = matrix(0, nrow(sv.G1Inact.G2Act.list), ncol(tempSV.act.inact))
    rownames(sv.G1Inact.G2Act.alter) = apply(sv.G1Inact.G2Act.list[, 1:2], 1, paste, collapse = ';') # Row names: co-inact SV genes
    colnames(sv.G1Inact.G2Act.alter) = colnames(tempSV.act.inact) # Column names: cell IDs
    for(m in 1:nrow(sv.G1Inact.G2Act.list)){   # m: loop through gene pairs
      genepair_act.inact = tempSV.act.inact[pmatch(sv.G1Inact.G2Act.list[m, 1:2], rownames(tempSV.act.inact)), ]
      #dim(genepair_act.inact) # [1] 2 12842 - alteration profile of the m-th gene pair in all cells
      # **** Reverse G1 inactivation + G2 activation SV pairs, reference DU, DD SR model, reversed state: G2 inactivation ****
      ind1 = which(genepair_act.inact[1, ] == (-1)) # G1 inactivated
      ind2 = which(genepair_act.inact[2, ] == (-1)) # G2 inactivated
      int.index = intersect(ind1, ind2)
      #print(paste(c('Gene pair ', m, ' of SV: ', length(int.index)), collapse = ''))
      # identical(names(s1), colnames(coinact.alter)) [1] TRUE
      sv.G1Inact.G2Act.alter[m, int.index] = 1  # **** G1 inactivation + G2 inactivation -> alteration status = 1 ****
    }
  }
  
  #   4. Summarize SV gene pair alteration counts
  sv.alter = rbind(sv.coinact.alter, sv.G1ActG2Inact.alter, sv.G1Inact.G2Act.alter)
  # Number of altered SV pairs per cell
  sv.alter.num = apply(sv.alter, 2, sum)
  # Integrated altered SV gene pairs per cell
  sv.alter.gp = rep(NA, length(sv.alter.num))
  for(n in 1:ncol(sv.alter)){sv.alter.gp[n] = paste(rownames(sv.alter)[which(sv.alter[, n] == 1)], collapse = '//')}
  
  #  **** Integrate and return results: DDS efficacy score for all cells  ****
  result = data.frame(cell = colnames(sl.coinact.alter), 
                      sl.alter.num = sl.alter.num, 
                      sl.alter = sl.alter.gp,
                      sv.Reverse.alter.num = sv.alter.num, 
                      sv.Reverse.alter = sv.alter.gp,
                      N = sum(nrow(sl.list), nrow(sv.list)))
  result$score = (result$sl.alter.num + result$sv.Reverse.alter.num) / result$N
  
  #  ****** Return integrated SL and SV efficacy score results for all cells of the i-th drug  *****
  return(result)   
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ function ~~~~~END





####          --->      Tool6: Immunotherapy Tools      <-----  ##### 
#####   Immunotherapy Tools, Function1    #####
#####   Fisher's Test to Identify Therapy Response-Related SL/SV Gene Pairs   #####

# Input data dimensions
dim(response.57) # [1] 57 32
dim(express.Gene.57) # [1] 1187   57
dim(network.sl)  # [1] 1249   23
dim(network.sv)  # [1] 206    23


# Define function to perform Fisher's test and return results
perform_fisher_test <- function(data, analysis_name) {
  # Perform Fisher's exact test
  test_result <- fisher.test(data)
  
  # Extract results
  result_df <- data.frame(
    Analysis = analysis_name,
    P_value = test_result$p.value,
    Odds_Ratio = test_result$estimate,
    OR_Lower_95 = test_result$conf.int[1],
    OR_Upper_95 = test_result$conf.int[2]
  )
  
  return(result_df)
}


Fisher.Immunotherapy.SLSV = function(express.act.inact, network.sl, network.sv){
  # Initialize results for SL pairs
  result = data.frame(gene1 = '', gene2 = '', label1 = '', a = 'a', b = 'b', c = 'c', d = 'd',
                      Analysis = '', P_value = '',   
                      Odds_Ratio = '', OR_Lower_95 = '', OR_Upper_95 = '')
  
  # Loop through SL gene pairs
  for(i in 1:nrow(network.sl)){
    index = pmatch(network.sl[i, 5:6], rownames(express.act.inact))
    if(sum(is.na(index)) == 0){
      dat = data.frame(gene1 = express.act.inact[index[1], ],
                       gene2 = express.act.inact[index[2], ])
      colnames(dat)[1:2] = network.sl[i, 5:6]
      # table(dat[, 2:3])
      dat$sum = apply(dat[, 1:2], 1, sum)
      
      # Four alteration patterns for gene1-gene2 pairs
      g1g2.inact = rownames(dat)[which(dat$sum == (-2))] # ---> altered (gene1-gene2 co-inactivation)
      other = setdiff(rownames(dat), g1g2.inact)  # unaltered
      a = length(intersect(g1g2.inact, resp.sample)) # responders with alteration
      b = length(g1g2.inact) - a  # non-responders with alteration = total altered - responders with alteration
      c = length(intersect(resp.sample, other))  # responders without alteration
      d = length(other) - c  # non-responders without alteration = total unaltered - responders without alteration
      sl_data <- matrix(c(a, b,  # Response group: SL altered = 15, SL wildtype = 5
                          c, d), # Non-response group: SL altered = 8, SL wildtype = 12
                        nrow = 2, byrow = F)
      rownames(sl_data) <- c("Response", "NoResponse")
      colnames(sl_data) <- c("SL_Altered", "SL_Wildtype")
      
      # Perform test on SL data
      sl_result <- perform_fisher_test(sl_data, "SL_Alteration_vs_Response")
      # Compile results for this pattern
      temp = cbind(data.frame(gene1 = colnames(dat)[1],
                              gene2 = colnames(dat)[2],
                              label1 = 'g1g2.inact',
                              a = a, b = b, c = c, d = d),
                   sl_result)
      
      g1g2.act = rownames(dat)[which(dat$sum == 2)] # ---> altered (gene1-gene2 co-activation)
      other = setdiff(rownames(dat), g1g2.act)  # ---> unaltered
      a = length(intersect(g1g2.act, resp.sample)) # responders with alteration
      b = length(g1g2.act) - a  # non-responders with alteration = total altered - responders with alteration
      c = length(intersect(resp.sample, other))  # responders without alteration
      d = length(other) - c  # non-responders without alteration = total unaltered - responders without alteration
      sl_data <- matrix(c(a, b,  # Response group: SL altered = 15, SL wildtype = 5
                          c, d), # Non-response group: SL altered = 8, SL wildtype = 12
                        nrow = 2, byrow = F)
      rownames(sl_data) <- c("Response", "NoResponse")
      colnames(sl_data) <- c("SL_Altered", "SL_Wildtype")
      # Perform test on SL data
      sl_result <- perform_fisher_test(sl_data, "SL_Alteration_vs_Response")
      # Compile results for this pattern
      temp.2 = cbind(data.frame(gene1 = colnames(dat)[1],
                                gene2 = colnames(dat)[2],
                                label1 = 'g1g2.act',
                                a = a, b = b, c = c, d = d),
                     sl_result)
      # Compile results for all patterns
      temp = rbind(temp, temp.2)
      
      g1Act.g2Inact = rownames(dat)[intersect(which(dat[, 2] == 1), which(dat[, 3] == (-1)))] # ---> gene1 activation - gene2 inactivation
      other = setdiff(rownames(dat), g1Act.g2Inact) # ---> unaltered
      a = length(intersect(g1Act.g2Inact, resp.sample)) # responders with alteration
      b = length(g1Act.g2Inact) - a  # non-responders with alteration = total altered - responders with alteration
      c = length(intersect(resp.sample, other))  # responders without alteration
      d = length(other) - c  # non-responders without alteration = total unaltered - responders without alteration
      sl_data <- matrix(c(a, b,  # Response group: SL altered = 15, SL wildtype = 5
                          c, d), # Non-response group: SL altered = 8, SL wildtype = 12
                        nrow = 2, byrow = F)
      rownames(sl_data) <- c("Response", "NoResponse")
      colnames(sl_data) <- c("SL_Altered", "SL_Wildtype")
      # Perform test on SL data
      sl_result <- perform_fisher_test(sl_data, "SL_Alteration_vs_Response")
      # Compile results for this pattern
      temp.2 = cbind(data.frame(gene1 = colnames(dat)[1],
                                gene2 = colnames(dat)[2],
                                label1 = 'g1Act.g2Inact',
                                a = a, b = b, c = c, d = d),
                     sl_result)
      # Compile results for all patterns
      temp = rbind(temp, temp.2)
      
      g1Inact.g2Act = rownames(dat)[intersect(which(dat[, 2] == (-1)), which(dat[, 3] == 1))] # ---> gene1 inactivation - gene2 activation
      other = setdiff(rownames(dat), g1Inact.g2Act)  # ---> unaltered
      a = length(intersect(g1Inact.g2Act, resp.sample)) # responders with alteration
      b = length(g1Inact.g2Act) - a  # non-responders with alteration = total altered - responders with alteration
      c = length(intersect(resp.sample, other))  # responders without alteration
      d = length(other) - c  # non-responders without alteration = total unaltered - responders without alteration
      sl_data <- matrix(c(a, b,  # Response group: SL altered = 15, SL wildtype = 5
                          c, d), # Non-response group: SL altered = 8, SL wildtype = 12
                        nrow = 2, byrow = F)
      rownames(sl_data) <- c("Response", "NoResponse")
      colnames(sl_data) <- c("SL_Altered", "SL_Wildtype")
      # Perform test on SL data
      sl_result <- perform_fisher_test(sl_data, "SL_Alteration_vs_Response")
      # Compile results for this pattern
      temp.2 = cbind(data.frame(gene1 = colnames(dat)[1],
                                gene2 = colnames(dat)[2],
                                label1 = 'g1Inact.g2Act',
                                a = a, b = b, c = c, d = d),
                     sl_result)
      # Compile results for all patterns
      temp = rbind(temp, temp.2)
      
      result = rbind(result, temp)
    }
  }
  
  # Extract significant results
  result005 = result[which(result$P_value < 0.05), ]
  result005 = result005[-1, ]
  result005 = result005[which(result005$Odds_Ratio > 1), ]
  result005.SL = result005
  print(paste(c(dim(result005.SL)[1], ' immunotherapy response-related SL gene pairs!'), collapse = ''))
  
  # Initialize results for SV pairs
  result = data.frame(gene1 = '', gene2 = '', label1 = '', a = 'a', b = 'b', c = 'c', d = 'd',
                      Analysis = '', P_value = '',   
                      Odds_Ratio = '', OR_Lower_95 = '', OR_Upper_95 = '')
  
  # Loop through SV gene pairs
  for(i in 1:nrow(network.sv)){
    index = pmatch(network.sv[i, 5:6], rownames(express.act.inact))
    if(sum(is.na(index)) == 0){
      dat = data.frame(gene1 = express.act.inact[index[1], ],
                       gene2 = express.act.inact[index[2], ])
      colnames(dat)[1:2] = network.sv[i, 5:6]
      # table(dat[, 2:3])
      dat$sum = apply(dat[, 1:2], 1, sum)
      
      # Four alteration patterns for gene1-gene2 pairs
      g1g2.inact = rownames(dat)[which(dat$sum == (-2))] # ---> gene1-gene2 co-inactivation
      other = setdiff(rownames(dat), g1g2.inact)  # ---> unaltered
      a = length(intersect(g1g2.inact, nonresp.sample)) # non-responders with alteration
      b = length(g1g2.inact) - a  # responders with alteration = total altered - non-responders with alteration
      c = length(intersect(nonresp.sample, other))  # non-responders without alteration
      d = length(other) - c  # responders without alteration = total unaltered - non-responders without alteration
      sv_data <- matrix(c(a, b,  # Response group: SL altered = 15, SL wildtype = 5
                          c, d), # Non-response group: SL altered = 8, SL wildtype = 12
                        nrow = 2, byrow = F)
      rownames(sv_data) <- c("NoResponse", "Response")
      colnames(sv_data) <- c("SV_Altered", "SV_Wildtype")
      
      # Perform test on SV data
      sv_result <- perform_fisher_test(sv_data, "SV_Alteration_vs_Response")
      # Compile results for this pattern
      temp = cbind(data.frame(gene1 = colnames(dat)[1],
                              gene2 = colnames(dat)[2],
                              label1 = 'g1g2.inact',
                              a = a, b = b, c = c, d = d),
                   sv_result)
      
      g1g2.act = rownames(dat)[which(dat$sum == 2)] # ---> altered: gene1-gene2 co-activation
      other = setdiff(rownames(dat), g1g2.act)  # ---> unaltered
      a = length(intersect(g1g2.act, nonresp.sample)) # non-responders with alteration
      b = length(g1g2.act) - a  # responders with alteration = total altered - non-responders with alteration
      c = length(intersect(nonresp.sample, other))  # non-responders without alteration
      d = length(other) - c  # responders without alteration = total unaltered - non-responders without alteration
      sv_data <- matrix(c(a, b,  # Response group: SL altered = 15, SL wildtype = 5
                          c, d), # Non-response group: SL altered = 8, SL wildtype = 12
                        nrow = 2, byrow = F)
      rownames(sv_data) <- c("NoResponse", "Response")
      colnames(sv_data) <- c("SV_Altered", "SV_Wildtype")
      
      # Perform test on SV data
      sv_result <- perform_fisher_test(sv_data, "SV_Alteration_vs_Response")
      # Compile results for this pattern
      temp.2 = cbind(data.frame(gene1 = colnames(dat)[1],
                                gene2 = colnames(dat)[2],
                                label1 = 'g1g2.inact',
                                a = a, b = b, c = c, d = d),
                     sv_result)
      # Compile results for all patterns
      temp = rbind(temp, temp.2)
      
      g1Act.g2Inact = rownames(dat)[intersect(which(dat[, 2] == 1), which(dat[, 3] == (-1)))] # ---> gene1 activation - gene2 inactivation
      other = setdiff(rownames(dat), g1Act.g2Inact)
      a = length(intersect(g1Act.g2Inact, nonresp.sample)) # non-responders with alteration
      b = length(g1Act.g2Inact) - a  # responders with alteration = total altered - non-responders with alteration
      c = length(intersect(nonresp.sample, other))  # non-responders without alteration
      d = length(other) - c  # responders without alteration = total unaltered - non-responders without alteration
      sv_data <- matrix(c(a, b,  # Response group: SL altered = 15, SL wildtype = 5
                          c, d), # Non-response group: SL altered = 8, SL wildtype = 12
                        nrow = 2, byrow = F)
      rownames(sv_data) <- c("NoResponse", "Response")
      colnames(sv_data) <- c("SV_Altered", "SV_Wildtype")
      
      # Perform test on SV data
      sv_result <- perform_fisher_test(sv_data, "SV_Alteration_vs_Response")
      # Compile results for this pattern
      temp.2 = cbind(data.frame(gene1 = colnames(dat)[1],
                                gene2 = colnames(dat)[2],
                                label1 = 'g1g2.inact',
                                a = a, b = b, c = c, d = d),
                     sv_result)
      # Compile results for all patterns
      temp = rbind(temp, temp.2)
      
      g1Inact.g2Act = rownames(dat)[intersect(which(dat[, 2] == (-1)), which(dat[, 3] == 1))] # ---> gene1 inactivation - gene2 activation
      other = setdiff(rownames(dat), g1Inact.g2Act)
      a = length(intersect(g1Inact.g2Act, nonresp.sample)) # non-responders with alteration
      b = length(g1Inact.g2Act) - a  # responders with alteration = total altered - non-responders with alteration
      c = length(intersect(nonresp.sample, other))  # non-responders without alteration
      d = length(other) - c  # responders without alteration = total unaltered - non-responders without alteration
      sv_data <- matrix(c(a, b,  # Response group: SL altered = 15, SL wildtype = 5
                          c, d), # Non-response group: SL altered = 8, SL wildtype = 12
                        nrow = 2, byrow = F)
      rownames(sv_data) <- c("NoResponse", "Response")
      colnames(sv_data) <- c("SV_Altered", "SV_Wildtype")
      
      # Perform test on SV data
      sv_result <- perform_fisher_test(sv_data, "SV_Alteration_vs_Response")
      # Compile results for this pattern
      temp.2 = cbind(data.frame(gene1 = colnames(dat)[1],
                                gene2 = colnames(dat)[2],
                                label1 = 'g1g2.inact',
                                a = a, b = b, c = c, d = d),
                     sv_result)
      # Compile results for all patterns
      temp = rbind(temp, temp.2)
      
      result = rbind(result, temp)
    }
  }
  
  # Extract significant results
  result005 = result[which(result$P_value < 0.05), ]
  result005 = result005[-1, ]
  result005 = result005[which(result005$Odds_Ratio > 1), ]
  result005.SV = result005
  print(paste(c(dim(result005.SV)[1], ' immunotherapy response-related SV gene pairs!'), collapse = ''))
  
  # Return results
  immunotherapy.result = list(result005.SL, result005.SV)
  return(immunotherapy.result)
}

