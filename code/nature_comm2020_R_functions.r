#===========================================================
# This function retrieve official gene symbols, entrezID and ENSG for Darmanis data
#===========================================================
updata.darmanis.geneIDs20200416.1641=function(dar){
  # Get Gene symbols' Entrez and ENSG ids
  sym=dar$V1
  # Ensembl DB
  ens75 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org",
                  path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
  ens38=useEnsembl('ensembl',dataset = "hsapiens_gene_ensembl")
  sym1=getBM(attributes = c('hgnc_symbol','ensembl_gene_id','entrezgene'),
             filters = 'hgnc_symbol',values = sym,mart=ens75)
  sym2=getBM(attributes = c('hgnc_symbol','ensembl_gene_id','entrezgene'),
             filters = 'hgnc_symbol',values = sym,mart=ens38)
  sym3=as.data.table(rbind(sym1,sym2))
  sym3=sym3[,lapply(.SD,function(x){paste(collapse = ",",unique(x))}),by='hgnc_symbol']
  sym3$entrezgene=gsub(sym3$entrezgene,pattern = 'NA',replacement = "")
  sym3$entrezgene=unlist(llply(sym3$entrezgene,.progress = 'text',
                               function(x){paste(collapse = ",",unlist(strsplit(x,",")))}))
  sym3=cbind(sym3$hgnc_symbol,sym3);colnames(sym3)[1:2]=c('symbol','official_symbol')
  
  # NCBI seraching
  sym=setdiff(sym,sym3$hgnc_symbol)
  sym4=ldply(sym,.progress = 'text',
             function(x){
               tryCatch(Update_EntrezID_Officialsymbol_From_Symbol(symbol = x),error=function(e){})
             })
  colnames(sym4)=c('symbol','official_symbol','entrezgene')
  # Get ENSG id from entrez
  sym4a=getBM(attributes = c('entrezgene','ensembl_gene_id'),filters = 'entrezgene',
              values=sym4$entrezgene,mart=ens38)
  sym4b=getBM(attributes = c('entrezgene','ensembl_gene_id'),filters = 'entrezgene',
              values=sym4$entrezgene,mart=ens75)
  sym4c=as.data.table(rbind(sym4a,sym4b))
  sym4c=sym4c[,lapply(.SD,function(x){paste(collapse = ",",unique(x))}),by='entrezgene']
  # Combine the gene matrix
  sym5=merge(sym4,sym4c,by='entrezgene',all=T)
  head(sym5)
  sym5a=sym5[,c(2,3,4,1)]
  head(sym5a)
  sym6=rbind(sym5a,sym3)
  
  # Manual searching  
  sym=setdiff(sym, c(sym6$symbol, sym6$official_symbol))
  #BDAG1->CCDC180,100499483,ENSG00000197816
  #C2orf28->ATRAID,51374,ENSG00000138085
  #CCL14-CCL15->CCL15-CCL14,348249,ENSG00000275688
  #LOC100130000->LOC100132057,100132057,ENSG00000231551
  #LOC100286793->LINC01138,388685,ENSG00000274020
  #MGC34034->LINC01312,154089,ENSG00000223586
  #PP14571->LOC100130449,100130449,ENSG00000218416
  #RN45S->RNA45S,?,?
  #RP1-177G6.2->LINC00632,286411,ENSG00000203930
  #RP11-165H20.1->CHIAP2,149620,ENSG00000203878
  sym6a=rbind(c('BDAG1','CCDC180','100499483','ENSG00000197816'),
              c('C2orf28','ATRAID','51374','ENSG00000138085'),
              c('CCL14-CCL15','CCL15-CCL14','348249','ENSG00000275688'),
              c('LOC100130000','LOC100132057','100132057','ENSG00000231551'),
              c('LOC100286793','LINC01138','388685','ENSG00000274020'),
              c('MGC34034','LINC01312','154089','ENSG00000223586'),
              c('PP14571','LOC100130449','100130449','ENSG00000218416'),
              c('RN45S','RNA45S','',''),
              c('RP1-177G6.2','LINC00632','286411','ENSG00000203930'),
              c('RP11-165H20.1','CHIAP2','149620','ENSG00000203878')
  )
  colnames(sym6a)=c('symbol','official_symbol','entrezgene','ensembl_gene_id')
  sym6a=sym6a[,c(1:2,4,3)]
  # Combine the all results
  sym7=rbind(sym6,sym6a)
  
  #NCBI searching for genes having ensembl_gene_id but no entrezgene
  sym_ensg=unique(unlist(strsplit(sym7$ensembl_gene_id[is.na(sym7$entrezgene)|sym7$entrezgene==""],",")))
  sym7a=ldply(sym_ensg,.progress='text',
              function(x){tryCatch(Update_EntrezID_Officialsymbol_From_Symbol(x),error=function(e){
              })})
  sym_ensg=setdiff(sym_ensg,sym7a$V1)
  sym7b=cbind(sym7a[,c(2,3,1)]);colnames(sym7b)=c('official_symbol','entrezgene','ensembl_gene_id')
  # put information to the matrix
  put_last_info_to_sym7=function(x=as.data.frame(sym7,stringsAsFactors=F),
                                 y=as.data.frame(sym7b,stringsAsFactors=F)){
    for(i in 1:nrow(y)){
      wh=grep(x$ensembl_gene_id,pattern = y$ensembl_gene_id[i])
      tmp=paste(x$entrezgene[wh],y$entrezgene[i],sep = ",")
      x$entrezgene[wh]=tmp
    }
    x$entrezgene=unlist(llply(x$entrezgene,.progress = 'text',function(k){
      tmp=unlist(strsplit(k,","))
      tmp=paste(unique(tmp[tmp!=""]),collapse=",")
      return(tmp)
    }))
    return(x)
  }
  sym8=put_last_info_to_sym7()
  # Check whether all symbols are matched
  setdiff(dar$V1,sym8$symbol)
  # Merge the table
  dar2=merge(sym8,dar,by.y='V1',by.x='symbol',all.y=T)
  # Remove LRG~~~id in dar2
  dar2$ensembl_gene_id=unlist(llply(dar2$ensembl_gene_id,.progress='text',
                                    function(x){
                                      y=unlist(strsplit(x,","))
                                      y=y[grep(y,pattern = 'ENSG')]
                                      y=paste(y,collapse = ",")
                                      return(y)
                                    }))
  return(dar2)
}


#===========================================================
# This function retrieve official gene symbols using retnrez package
#===========================================================
Update_EntrezID_Officialsymbol_From_Symbol=function (symbol) 
{
  pkgs = c("curl", "haven", "vctrs", "openxlsx", "scales", 
           "prodlim", "RSQLite", "reticulate", "httr", "curl", "rentrez")
  uninstalled.pkgs = setdiff(pkgs, rownames(installed.packages()))
  if (length(uninstalled.pkgs) > 0) {
    install.packages(pkgs)
  }
  library(httr)
  library(curl)
  library(rentrez)
  query1 <- paste("9606", "[TID] AND ", symbol, "[GENE]", sep = "")
  res1 <- rentrez::entrez_search(db = "gene", term = query1, 
                                 retmax = 1, version = "2.0")
  cv1_1 <- tryCatch(rentrez::entrez_summary(db = "gene", id = res1$ids), 
                    error = function(e) {
                      return(NULL)
                    })
  cv1_2 <- as.vector(unlist(tryCatch(rentrez::extract_from_esummary(cv1_1, 
                                                                    c("name", "uid")), error = function(e) {
                                                                      return(NULL)
                                                                    })))
  if (!is.null(cv1_2) & length(cv1_2) == 2) {
    mt1 <- c(symbol, cv1_2)
    return(mt1)
  }
  if (is.null(cv1_2)) {
    query2 <- paste("9606", "[TID] AND ", symbol, "[ALL]", 
                    sep = "")
    res2 <- rentrez::entrez_search(db = "gene", term = query2, 
                                   retmax = 1)
    cv2_1 <- tryCatch(rentrez::entrez_summary(db = "gene", 
                                              id = res2$ids), error = function(e) {
                                                return(NULL)
                                              })
    cv2_2 <- as.vector(unlist(tryCatch(rentrez::extract_from_esummary(cv2_1, 
                                                                      c("name", "uid")), error = function(e) {
                                                                        return(NULL)
                                                                      })))
    if (!is.null(cv2_2) & length(cv2_2) == 2) {
      mt2 <- c(symbol, cv2_2)
      return(mt2)
    }
    if (is.null(cv2_2)) {
      unigene.query <- paste("9606", "[TXID] AND ", symbol, 
                             "[ALL]", sep = "")
      unigene.res2 <- rentrez::entrez_search(db = "unigene", 
                                             term = unigene.query, retmax = 1)
      unigene.1 <- rentrez::entrez_summary(db = "unigene", 
                                           id = unigene.res2$ids)
      unigene.2 <- as.vector(unlist(rentrez::extract_from_esummary(unigene.1, 
                                                                   c("gene"))))
      query3 <- paste("9606", "[TID] AND ", unigene.2, 
                      "[GENE]", sep = "")
      res3 <- rentrez::entrez_search(db = "gene", term = query3, 
                                     retmax = 1)
      cv3_1 <- rentrez::entrez_summary(db = "gene", id = res3$ids)
      cv3_2 <- as.vector(unlist(rentrez::extract_from_esummary(cv3_1, 
                                                               c("name", "uid"))))
      if (!is.null(cv3_2) & length(cv3_2) == 2) {
        mt3 <- c(symbol, cv3_2)
        return(mt3)
      }
    }
  }
}

#===========================================================
# This function downloads gene expression data of TCGA
#===========================================================
Download.tcga.expr.data20200416.1405=function(){
  check.n.install.lib(lib.name = c('TCGAbiolinks','SummarizedExperiment'),
                      lib.type = rep('bioc',2))
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  # download gene expression data
  query.exp <- GDCquery(project = "TCGA-GBM",
                        data.category = "Transcriptome profiling",
                        data.type = "Gene expression quantification",
                        #platform = "Illumina HiSeq",
                        workflow.type = "HTSeq - FPKM",
                        experimental.strategy='RNA-seq',
                        legacy = F)
  
  GDCdownload(query.exp)
  gbm_exp=GDCprepare(query.exp)
  
  #Expression matirx
  exp=assay(gbm_exp)
  # Check whether there is normal
  # The sample annotation is 4th id, Normal is 10-19 and control sample is 20-29
  # tumor is 01-09
  # visit https://wiki.nci.nih.gov/display/TCGA/TCGA+barcode
  tums=paste0(0,1:9)
  wn=c();for(i in 1:ncol(exp)){
    tmp=unlist(strsplit(colnames(exp)[i],'-'))[4]
    tmp=gsub(tmp,pattern = '[A-Z]',replacement = "")
    if(tmp %in% tums==F){
      wn=c(wn,i)
    }
  }
  ex2=exp[,-wn] # Exclude normal and control samples
  # Save the data
  return(ex2)
}

#===========================================================
# This function retrieves official genesymbols and entrez symbols
#===========================================================
attach.gene.symbol2TCGA.exp20190117.1525=function(x){
  check.n.install.lib(lib.name=c('org.Hs.eg.db','EnsDb.Hsapiens.v75','biomaRt','plyr'),
                      lib.type=c(rep('bioc',3),'cran'))
  # Combine gbm expression and lgg expression data
  x1=x
  # Get gene symbols
  # Org.hs.eg.db
  library(org.Hs.eg.db)
  org.sym=select(org.Hs.eg.db,keys = rownames(x1),columns = c("ENSEMBL","ENTREZID","SYMBOL"),keytype = 'ENSEMBL')
  org.sym=na.omit(org.sym)
  fail.1=setdiff(rownames(x1),org.sym$ENSEMBL)
  # ENSDB 75
  library(EnsDb.Hsapiens.v75)
  ens.sym=select(EnsDb.Hsapiens.v75,keys=fail.1,columns = c("ENTREZID","SYMBOL"),keytype = "GENEID")
  ens.sym=na.omit(ens.sym)
  fail.2=setdiff(fail.1,ens.sym$GENEID)
  # biomaRt
  library(biomaRt)
  if(!exists('ens38')){
    ens38 <<- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
    ens37 <<- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",GRCh = 37)
  }
  ens38.sym=getBM(mart=ens38,values = fail.2,attributes =c('ensembl_gene_id','entrezgene','hgnc_symbol'),filters ='ensembl_gene_id')
  ens38.sym$entrezgene[which(is.na(ens38.sym$entrezgene))]=''
  ens38.sym=ens38.sym[which(ens38.sym$hgnc_symbol!='' | ens38.sym$entrezgene!=''),]
  fail.3=setdiff(fail.2,ens38.sym$ensembl_gene_id)
  ens37.sym=getBM(mart=ens37,values = fail.3,attributes =c('ensembl_gene_id','entrezgene','hgnc_symbol'),filters ='ensembl_gene_id')
  ens37.sym$entrezgene[which(is.na(ens37.sym$entrezgene))]=''
  ens37.sym=ens37.sym[which(ens37.sym$hgnc_symbol!=''|ens37.sym$entrezgene!=''),]
  fail.4=setdiff(fail.3,ens37.sym$ensembl_gene_id)
  # Search using entrezR
  library(plyr)
  entre.sym=ldply(1:length(fail.4),.progress='text',function(k){
    tryCatch({
      tmp=Update_EntrezID_Officialsymbol_From_Symbol(fail.4[k])
      return(tmp)
    },error=function(e){})
  })
  fail.5=setdiff(fail.4,entre.sym[,1])
  colnames(entre.sym)=c('ENSG','symbol','Entrez')
  entre.sym=entre.sym[,c('ENSG','Entrez','symbol')]
  head(entre.sym)
  # Combine the symbol results
  colnames(org.sym)=colnames(ens.sym)=colnames(ens38.sym)=colnames(ens37.sym)=colnames(entre.sym)
  genes=rbind(org.sym,ens.sym,ens38.sym,ens37.sym,entre.sym)
  # Shirink the table to reduce redundancy
  genes=dup.matrix(genes,key_column = 'ENSG',my.Function = function(k){paste(collapse=';',unique(k))})
  # attach the symbols to table
  x1=as.data.frame(x1,stringsAsFactors=F)
  x1$ENSG=rownames(x1)
  x1=log2(x1+1)
  x2=mat.merge(x=genes,y=x1,by='ENSG',all=T)
  return(list(merged_tab=x2,symbols=genes,failed.ENSG=fail.5))
}

#===========================================================
# This function loads an excel file into R environment.
# x : File path.
# sheet : Sheet name of sheet number.
#===========================================================
read.xl=function(x='Supplementary Data 1.xlsx',sheet = 2){
  check.n.install.lib('readxl')
  library(readxl)
  x1=readxl::read_excel(path = x,sheet = sheet)
  x1=as.data.frame(x1,stringsAsFactors=F)
  return(x1)
}

#===========================================================
# This function expands the number of rows of provided matrix or data.frame according to the designated key column.
# @param key_colnames (default=NA) : Column name of matrix containing key values to expand table.
# @param mat (default=NA) : Matrix to be expanded.
# @param sep (defaule=';') : Character separater to separate text in the field.
# (ex. TP53;EGFR -> [1]TP53 [2]EGFR)
#===========================================================
table_expand=function(key_colnames=NA,mat=c(),sep=';'){
  check.n.install.lib('tidyr')
  library('tidyr')
  if(is.na(key_colnames)){
    cat('key_colnames : column name of matrix to expand the table\n')
    cat('mat : matrix\n')
    cat('sep : separator character will separate the values\n')
    cat('ex. ERBB1;EGFR->ERBB1, EGFR-> 2 rows matrix')
    break()
  }else{
    # Select rows which have designated separator value
    mat=as.data.frame(mat,stringsAsFactors=F)
    wh=which(colnames(mat)==key_colnames)
    gr=grep(mat[,wh],pattern=sep,fixed=T)
    if(length(gr)>0){
      # Split the line
      mat.1=tidyr::separate_rows(mat[gr,],key_colnames,sep=sep,convert = T)
      # Combine matrix
      mat.2=mat[-gr,]
      mat.3=rbind(mat.1,mat.2)
    }else{mat.3=mat}
    return(mat.3)
  }
}

#===========================================================
# This function handles tables containing duplicates.
# @param mm (default=NA) : Matrix to be processed.
# @param key_column (default=NA) : Name of column with duplicates.
# @param my.Function (default=NA) : Function executed for duplicates.
#===========================================================
dup.matrix=function(mm=NA,key_column=NA,my.Function=NA){
  check.n.install.lib(c('data.table','plyr'))
  library(data.table);library(plyr)
  
  # Check duplicates in key_column
  print('Check duplciates')
  dups=mm[,key_column]
  dups=unique(dups[which(duplicated(dups))])
  if(length(dups)>0){
    # Separate duplicates into new matrix
    print('Check duplicate matrix')
    wh=which(mm[,key_column] %in% dups)
    if(nrow(mm)-length(wh)!=1){
      mm.dup=mm[wh,]
      mm.uniq=mm[-wh,]
    }else{mm.dup=mm;mm.uniq=c()}
    # Apply function to the duplicate matrix
    print('Apply function to the duplicate matrix')
    mm.dup=data.table::as.data.table(mm.dup)
    mm.dup2=mm.dup[,lapply(.SD,{my.Function}),by=key_column]
    mm.dup2=as.data.frame(mm.dup2,stringsAsFactors=F)
    # Combine unique and processed duplicate matrix
    print('Combine unique and processed duplicate matrix')
    mm2=rbind(mm.uniq,mm.dup2)
    # Export the results
    gc()
    return(mm2)
  }else{
    # If there are no duplicates
    print('No duplicates')
    return(mm)
  }
}

#===========================================================
# This function evaluates adequacy of class labeling.
# Suppose, there is a gene expression matrix with 100 samples classified into two subtypes, A and B.
# If the subtyping is appropriate,
# distances samples of different subtypes (inter) will be distant while that between samples classified to the same subtype (intra) will be close.
# This function evaluates intra- and inter-subtype distances.
# Pseudocode is as following:
# function{
#   a <- Calculate distances of all samples
#   a1 <- Select intra-connectivity (i.e Distances between samples of same subtypes)
#   a2 <- Separate intra-connectivity into lists according to the number of subtypes.
#   a3 <- Get average intra-connectivity (i.e average distance from samples of identical subtype)
#   b1 <- Select inter-connectivity (i.e Distances between samples in different subtypes)
#   b2 <- Separate inter-connectivity into lists according to the number of subtypes.
#   b3 <- Get average inter-connectivity (i.e average distance from samples of different subtypes)
#   return(sum(a3)-sum(b3))
# }
# @param mm (default=NA) : matrix to be processed.
# @param class.list : List object that has sample name and classes.
# (ex. list('subtypeA'=c('A1','A2','A3'),'subtypeB'=c('b1','b2','b3')))
# @param distance.method (default='euclidean') : Calculation method for distance.
#												 Allowed values are 'euclidean' and 'correlation'.
#===========================================================
calc_cluster.index=function(mm,class.list,distance.method='euclidean'){
  check.n.install.lib(c('wordspace','plyr'))
  x=mm[,unlist(class.list)]
  # Calculate distance
  if(dim(x)[1]!=1){
    x1=switch(distance.method,
              euclidean=as.matrix(wordspace::dist.matrix(t(x),as.dist=F)),
              correlation=as.matrix(1-cor(x))
    )
  }else{
    x1=dist(t(x),diag = T,upper = T)
  }
  # Calculate intra connectivity
  intra.con=plyr::llply(class.list,function(k){
    k1=x1[k,k]
    k1=k1[upper.tri(k1)]
    return(mean(k1))
  })
  # Calculate inter connectivity
  inter.con=plyr::llply(1:(length(class.list)-1),function(a){
    inter.con=c()
    for(b in (a+1):length(class.list)){
      a1=class.list[[a]]
      b1=class.list[[b]]
      k1=x1[a1,b1]
      inter.con=c(inter.con,mean(k1))
    }
    return(inter.con)
  })
  # Calculate cluster index
  cluster.index=sum(unlist(intra.con))-sum(unlist(inter.con))
  # Export result
  return(cluster.index)
}

#===========================================================
# This function evaluates adequacy of class labeling by permutation test.
# Pseudocode is written below.
# function{
#   real_cluster_index_value <- calc_cluster.index(mm=matrix,class.list=subtype.list)
#	random_cluster_index_values<-unlist(llply(i in 1:10000){
#				permutated_matrix <- Random shuffling of sample labels.
#				random_cluster_index_value <-calc_cluster.index(mm=permutated_matrix,class.list=subtype.list)
#				return(random_cluster_index_value)
#			})
#	p.value <- sum(random_cluster_index_values < real_cluster_index_value)/10000
#	return(p.value)
# }
# @param mm (default=NA) : matrix to be processed.
# @param class.list : List object that has sample names and classes.
# (ex. list('subtypeA'=c('A1','A2','A3'),'subtypeB'=c('b1','b2','b3')))
# @param distance.method (default='euclidean') : Calculation method for distance.
#												 Allowed values are 'euclidean' and 'correlation'.
#===========================================================
cluster_index_P=function(mm,class.list,distance.method='euclidean',perm=10000,plot=F){
  check.n.install.lib('plyr')
  library(plyr)
  # Get real value
  real.value=calc_cluster.index(mm=mm,class.list = class.list,distance.method = distance.method)
  # Calculate permutated values
  perm.value=unlist(plyr::llply(1:perm,.progress = 'text',function(k){
    # Permutation labels
    perm.label=colnames(mm)[sample(1:ncol(mm),ncol(mm))]
    colnames(mm)=perm.label
    x1=calc_cluster.index(mm=mm,class.list = class.list,distance.method = distance.method)
    return(x1)
  }))
  # Calculate p-value
  p.val=sum(unlist(perm.value)<real.value)/perm
  # Draw plot
  if(plot){
    plot(density(perm.value),main='Permutation result')
    abline(v=real.value,col='red')
    legend(lty=1,legend = 'Real cluster index',col='red','topleft')
  }
  # Return p-value
  return(p.val)
}


#===========================================================
# This function installs R packages
# @param lib.name : Name of package. You can provide them as character vector.
# @param lib.type : Source of package. Allowed options are 'cran' and 'bioc' (bioconductor).
# The number of items in lib.type should be identical with the number of items in lib.name.
#===========================================================
check.n.install.lib=function(lib.name='',lib.type='cran'){
  x1=lib.name %in% rownames(installed.packages())
  x2=lib.name[which(!x1)]
  y2=lib.type[which(!x1)]
  if(length(x2)>0){
    for(i in 1:length(x2)){
      if(y2[i]=='cran'){
        install.packages(x2[i])
      }else{
        BiocManager::install(x2[i])
      }
    }
  }
}

#===========================================================
# This function calculates PAC score.
# @param ccl.stat : Object provided by ConsensusClusterPlus function of ConsensusClusterPlus package.
# @param kvec : Integer vector containing the number of clusters (K) to be evaluated.
# @param threshold (default=c(0.1,0.9)) : Threshold range defining the intermediate sub-interval.
# @param barplot (default=T) : Whether PAC score bar plot should be plotted.
#===========================================================
pac_score=function(ccl.stat=NA,kvec=NA,threshold=c(0.1,0.9),barplot=T){
  if(all(is.na(ccl.stat))){
    stop('Provide ccl.stat!')
  }
  # Set K-vectors
  if(is.na(kvec)){
    kvec=2:length(ccl.stat)
  }
  # Calculate PAC score
  message('Calculate PAC score')
  x1 = threshold[1]; x2 = threshold[2]# threshold defining the intermediate sub-interval
  PAC = rep(NA,length(kvec))
  names(PAC) = paste("K=",kvec,sep="") # from 2 to maxK
  for(i in kvec){
    M = ccl.stat[[i]]$consensusMatrix
    Fn = ecdf(M[lower.tri(M)])
    PAC[i-1] = Fn(x2) - Fn(x1)
  }
  # Return PAC scores
  # Drawing barplot
  if(barplot){
    tmp=barplot(PAC,main='PAC score',ylab='PAC')
    points(tmp,PAC,'b',pch=19)
  }
  # Report optimal K-vector
  optk = kvec[which.min(PAC)]
  message(paste0('Optimal K is k=',optk))
  # Return the values
  return(list(pac.score=PAC))
}

#===========================================================
# This function is a modification of gsub function in R.
# @param input : Character vector to be processed.
# @param pattern : Pattern of characters to be modified.
# @param replacement : Replacement of the pattern. The number of pattern and replacement have to be identical.
# @param exact : Whether only values which exactly match given patterns are to be changed. If true, this function does not use gsub.
# (ex. multi_sub(input = c("ham", "hamburger"), pattern = "ham", replacement = "HELLO", exact = T) returns c("HELLO", "hamburger")
# @param partial : Whether character having patterns are partially modified.
# (ex. multi_sub(input = c("ham", "hamburger"), pattern = "ham", replacement = "HAM", partial = T) returns c("HAM", "HAMburger")
# @fixed : Same as gsub function. See ?gsub
# @ignore.case : Same as gsub function. See ?gsub 
#===========================================================
multi_sub=function(input=NULL,
                   pattern=NULL,
                   replacement=NULL,
                   exact=F,
                   partial=F,
                   fixed=T,
                   ignore.case=T){
  input2=input
  # Partial
  if(partial){
    message('partial=T. Input vector will be partially changed using gsub function.')
    for(i in 1:length(pattern)){
      input2=gsub(input2,pattern = pattern[i],replacement = replacement[i],ignore.case = ignore.case,
                  fixed = fixed)
    }
  }else{
    # Exact matching
    message('partial=F, Input vector will be completely changed using which/grepl function.')
    if(exact){message('exact=T. Which function is used to find pattern.')
      # Which
      for(i in 1:length(pattern)){
        wh=which(input2==pattern[i])
        if(length(wh)!=0){input2[wh]=replacement[i]}
      }
    }else{
      # grepl
      message('exact=F. grepl function is used to find pattern.')
      for(i in 1:length(pattern)){
        wh=which(grepl(input2,pattern = pattern[i]))
        if(length(wh)!=0){input2[wh]=replacement[i]}
      }
    }
  }
  return(input2)
}

#===========================================================
# This function is a modification of grep function in R.
#===========================================================
multi_grep=function(input=c(),pattern=c(),ignore.case=T,fixed=F){
  if(is.null(input)){
    print('input : input character vector')
    print('pattern : pattern list character vector')
    print('ignore.case (default T) : same with ignore.case')
    print('Output will be unique list of sites')
  }else{
    wh=c();for(j in 1:length(pattern)){
      wh=c(wh,grep(input,pattern = pattern[j],ignore.case = ignore.case,fixed = fixed))
    }
    wh=unique(wh)
    return(wh)
  }
}

#===========================================================
# This function retrieves domain ids of genes using biomaRt.
#===========================================================
get_domain_info=function(symbol='EGFR'){
  check.n.install.lib(c('plyr','rjson','RCurl'))
  # Start Ensembl
  if(!exists('ensemblDB')){
    message('ensemblDB creating')
    library(biomaRt)
    ensemblDB <<- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37, version=75)
  }
  # Get ENSP id
  message(paste0('Get ENSP id and pfam information for ',symbol))
  enspID=getBM(attributes = c('ensembl_peptide_id','hgnc_symbol','peptide'),values = symbol,
               filters = 'hgnc_symbol',mart = ensemblDB)
  enspID=enspID[which(enspID$ensembl_peptide_id!=''),]
  enspID$length=nchar(enspID$peptide)
  # Get pfam id for the selected ENSP id
  pfamID=getBM(attributes = c('ensembl_peptide_id','hgnc_symbol','pfam','pfam_start','pfam_end'),
               values=unique(enspID$ensembl_peptide_id),mart = ensemblDB,filters = 'ensembl_peptide_id')
  # Remove ENSP id without pfam id
  pfamID2=na.omit(pfamID)
  # Take largest ENSP protein
  enspID2=mat.merge(pfamID2,enspID,by="ensembl_peptide_id")
  large.ensp=enspID2$ensembl_peptide_id[which.max(enspID2$length)]
  pfamID3=pfamID2[which(pfamID2$ensembl_peptide_id==large.ensp),]
  # Get human readable domain name
  pfam.name=domain_annot(ids = unique(pfamID3$pfam),idtype='pfam')
  pfam.name=data.frame('pfam'=names(pfam.name),'domain.name'=pfam.name,stringsAsFactors = F)
  # Merge tables
  pfam.table=mat.merge(pfamID3,pfam.name,all = F,by='pfam')
  pfam.table$length=enspID2$length[which.max(enspID2$length)]
  # Return the table
  return(pfam.table)
}

#===========================================================
# This function combines two matrices.
# Parameters are identical to merge function in R
#===========================================================
mat.merge=function(x,y,all=F,all.x=F,all.y=F,by=NULL,by.x=NULL,by.y=NULL){
  check.n.install.lib('dplyr')
  library(dplyr)
  # Convert matrices
  if(!is.data.frame(x) | length(class(x))!=1){x=as.data.frame(x,stringsAsfactors=F)}
  if(!is.data.frame(y) | length(class(y))!=1){y=as.data.frame(y,stringsAsfactors=F)}
  # Which column is key for merging two tables?
  if(!is.null(by.x) & !is.null(by.y)){
    colnames(y)=suppressMessages(multi_sub(colnames(y),pattern = by.y,replacement = by.x,exact = T))
    by=by.x
  }
  # Change key features into character value
  y[,by]=as.character(y[,by])
  x[,by]=as.character(x[,by])
  # All data should be kept
  if(all){
    all.x=all.y=T
    x1=dplyr::full_join(x,y,by=by)
  }
  # Only intersect should be kept
  if(!all.x & !all.y){
    x1=dplyr::inner_join(x,y,by=by)
  }
  # Only matrix X should be kept
  if(all.x & !all.y){
    x1=dplyr::left_join(x,y,by=by)
  }
  # Only matrix Y should be kept
  if(!all.x & all.y){
    x1=dplyr::right_join(x,y,by=by)
  }
  # Return matrix
  x1=as.data.frame(x1,stringsAsFactors=F)
  # Sort the table. 1st column should be the key column
  wh=which(colnames(x1)==by)
  x1=x1[,c(wh,setdiff(1:ncol(x1),wh))]
  return(x1)
}

#===========================================================
# This function retrieves information of protein domains.
# It requires protein domain ids from pfam or uniprot
#===========================================================
domain_annot=function(ids='PF13900',idtype='pfam'){
  check.n.install.lib(lib.name = 'trackViewer',lib.type = 'bioc')
  if(is.null(ids)){
    print('ids : input character vector for pfam_id or uniprot_id')
    print('idtype : pfam(default) or uniprot')
    print('use biomart to get pfam id and uniprot id')
    print('try bellow')
    cat("ensembl=useMart(biomart='ensembl',dataset = 'hsapiens_gene_ensembl')\n")
    cat("tmp=getBM(attributes = c('hgnc_symbol','ensembl_peptide_id','pfam','pfam_start','pfam_end'),
          filters= 'ensembl_peptide_id',values='ENSP00000428056',mart=ensembl)")
  }else{
    library('rjson');library('RCurl');library(plyr)
    ids=ids
    if(length(idtype)>1){warning('Provide unique id (pfam or uniprot)');break}
    if(idtype=='pfam'){
      print('pfam initiating')
      id2=unlist(llply(ids,.progress='text',function(x){
        #x=ids
        cat('\r',x)
        url=getURL(paste0('http://pfam.xfam.org/family/',x))
        info=unlist(strsplit(url,"\""))
        info=info[grep(info,pattern = 'Summary:')]
        if(length(info)!=0){
          info=unlist(strsplit(info,": ",fixed = T))[2]
          info=unlist(strsplit(info,"<",fixed = T))[1]}else{info='Not_available'}
      }))
      print('Conversion done')
      if(length(id2)==length(ids)){
        names(id2)=ids
      }
    }
    
    if(idtype=='uniprot'){
      print('Uniprot domain information consists of position and names')
      id2=ldply(ids,.progress='text',function(x){
        cat('\r',x)
        url=readLines(paste0('http://www.uniprot.org/uniprot/',x))
        info=url[grep(url,pattern = 'key=Domain',fixed = F)]
        info=unlist(strsplit(info,split="=",perl = T))
        # Location
        inf2=info[grep(info,pattern = x)]
        inf2=unlist(lapply(inf2,function(k){unlist(strsplit(k,"&"))[1]})) 
        inf2=gsub(inf2,pattern = x,replacement = "")
        inf2=gsub(inf2,pattern = ']',replacement = "")
        inf2=gsub(inf2,pattern = "[",replacement = "",fixed = T)
        # names
        inf3=info[grep(info,pattern = "</span><span class")]
        inf3=unlist(lapply(inf3,function(k){unlist(strsplit(k,">"))[2]}))
        inf3=unlist(lapply(inf3,function(k){unlist(strsplit(k,"<"))[1]}))
        
        res=cbind(rep(x,length(inf2)),inf2,inf3)
        colnames(res)=c('uniprotID','Position_in_protein','Domain_names')
        mat=res
      })
    }
  }
  return(id2)
}

#===========================================================
# Function to identify distributions of p-values and prognosis.
# @param x : Expression matrix. Rows and columns are genes and samples, respectively.
# @param y : Survival information matrix. Column names should be 'sample', 'surv', and 'status'.
#===========================================================
sliding_log_rank_test<-function(x=NULL,y=NULL,bias=NULL){
  check.n.install.lib(lib.name = c('plyr','matrixStats','survminer'),
                      lib.type = c('cran','cran','cran'))
  library(plyr);library(matrixStats);library(survival)
  # Remove 0 exp genes
  zero.gene=matrixStats::rowMeans2(as.matrix(x))
  x1=x[which(zero.gene!=0),]
  # Remove samples without survival data
  x1=x1[,which(colnames(x1) %in% y$sample)]
  y1=y[which(y$sample %in% colnames(x1)),]
  # sorting the tables
  x1=x1[,order(colnames(x1))]
  y1=y1[order(y$sample),]
  # Conduct log-rank test. Assign at least 5 samples per group #k=1
  x2=plyr::ldply(1:nrow(x),.progress='text',function(k){
    df=data.frame(exp=unlist(x1[k,]),os=y1$surv,stat=y1$status,stringsAsFactors=F)
    df$exp=as.numeric(as.character(df$exp))
    df=df[order(df$exp,decreasing=F),]
    df$group='High'
    k1=plyr::ldply(5:(nrow(df)-5),function(j){
      df$group[1:j]='Low'
      j1=survival::Surv(df$os,df$stat)~df$group
      j2=survival::survfit(j1)
      m.surv=survminer::surv_median(j2)
      os=survival::survdiff(j1)
      pv=pchisq(os$chisq,length(os$n)-1,lower.tail = F)
      harzard.ratio=survival::coxph(j1)
      harzard.ratio2=exp(harzard.ratio$coefficients)
      if(grepl(names(harzard.ratio2),pattern = 'Low')){
        res=ifelse(harzard.ratio2>1,'Favorable','Unfavorable')
        res2=1/harzard.ratio2
      }else{
        res=ifelse(harzard.ratio2<1,'Favorable','Unfavorable')
        res2=harzard.ratio2
      }
      res3=c(pv,res,df$exp[j],res2)
      return(res3)
    })
    k1=as.data.frame(k1,stringsAsFactors=F)
    k1$V1=as.numeric(k1$V1)
    k1$V2=as.numeric(k1$V2)
    if(is.null(bias)){
      k2=k1
    }else{
      if(bias=='good'){
        k2=k1[which(k1[,2]=='Favorable'),]
      }else{
        k2=k1[which(k1[,2]=="Unfavorable"),]
      }
    }
    if(nrow(k2)==0){
      k2=data.frame(matrix(rep(NA, 4), nrow = 1, ncol = 4), stringsAsFactors = F)
    }
    colnames(k2)=c('p.value','prognosis','cut_off','Hazard.ratio')
    k2$gene<-rownames(x)[k]
    return(k2)
  })
  x2$cut_off=as.numeric(x2$cut_off)
  return(x2)
}

#===========================================================
# Survival analysis with Kaplan-Meier plots.
#===========================================================
  survival_analysis=function(class=NULL,survival_time=NULL,status=NULL,plot=T,colors=NULL,
                             xlab=NULL,title=NULL,export=T,scientific=F,nlarge=5){
    frms=survival::Surv(survival_time,status)~class
    os=survival::survfit(frms)
    # Calculate Log-rank p-value
    pv=survival::survdiff(frms)
    pv=pchisq(pv$chisq,length(pv$n)-1,lower.tail = F)
    # Plotting?
    if(plot){
      plot(os,ylab='Survival probability',xlab=xlab,main=title,mark.time=T,cex=5,lwd=5,col=colors)
      # p-value formating
      if(scientific){
        if(pv<10^(-nlarge)){
          pv2=toupper(format(pv,scientific=T,nsmall=3,digits=3))
          pv2=gsub(pv2,pattern = '-0',replacement = '-')
        }else{
          pv2=format(pv,scientific=F,nsmall=3,digits=3)
        }
      }else{
        pv2=format(pv,scientific=F,nsmall=3,digits=3)
      }
      mtext(paste0('Log-rank p=',pv2),adj=1,side=3)
      # Legend
      legend.txt=gsub(paste0(names(os$strata),'(N=',os$n,')'),pattern='class=',replacement = "",fixed=T)
      legend('topright',legend = legend.txt,lwd=5,col=colors)
    }
    # Calculate hazard ratio
    # Get a subtype showing the best prognosis (longest median survival)
    median.os=quantile(os)$quantile[,2]
    median.os=sort(median.os,decreasing = T) # Sort the classes according to median survival period
    sorted.su=gsub(names(median.os),pattern = 'class=',replacement = '',fixed = T)
    # Change class name
    class2=as.factor(multi_sub(as.character(class),pattern = sorted.su,
                               replacement = 1:length(sorted.su),exact=T))
    frms2=survival::Surv(survival_time,status)~class2
    hazard.ratio=coxph(frms2)
    hazard.ratio2=exp(hazard.ratio$coefficients)
    # Return the hazard ratio
    hazard.ratio3=c(1,hazard.ratio2)
    names(hazard.ratio3)=sorted.su
    
    # Export result?
    if(export){return(list('surv_formula'=frms,'survfit'=os,'p.value'=pv,'hazard.ratio'=hazard.ratio3))}
  }

#==================================================
# Draw box plot with statistics.
#==================================================
sig_boxplot=function(input_list=NULL,test_method='t.test',one.sided=F,
                     ns_visualize=F,
                     nsmall=3,digits=3,yaxis_lab=NULL,colors=NULL,
                     title=NULL,jitter=T,n_of_numbers=5,violin=F,drawRect=F,
                     jitter.only=F){
  if(is.null(input_list)){stop('Input data is empty. Terminate process')}
  check.n.install.lib('vioplot')
  # Test the categories pairwise
  test_m=match.arg(test_method,choices = c('ks.test','t.test','wilcox.test'))
  test_n=length(input_list)-1
  test_res=c()
  for(i in 1:test_n){for(j in (i+1):length(input_list)){
    x=input_list[[i]];y=input_list[[j]]
    if(length(x)>1 & length(y)>1){
      if(!one.sided){ # Two-sided test
        res=switch(test_m,
                   t.test=t.test(x,y)$p.value,
                   wilcox.test=wilcox.test(x,y)$p.value,
                   ks.test=ks.test(x,y)$p.value)
      }else{
        res=one.sided_test(x=list(x,y),test_method = test_method)
      }
      tmp=format(res,nsmall=nsmall,digits=1,scientific=F)
      if(nchar(tmp)>=n_of_numbers){tmp=format(res,nsmall=nsmall,digits=digits,scientific=T)}
    }else{tmp=1}
    test_res=c(test_res,tmp)
    names(test_res)[length(test_res)]=paste0(names(input_list)[i],'_&&_',names(input_list)[j])
  }}
  # Remove non-sig pairs?
  if(!ns_visualize){test_res=test_res[as.numeric(test_res)<0.05]}else{
    test_res[as.numeric(test_res)>0.05]=paste0('NS(',test_res[as.numeric(test_res)>0.05],')')}  
  # e-0 to E-
  test_res=gsub(test_res,pattern = 'e-0',replacement = "E-")
  # Draw box plot
  y_range=range(unlist(input_list));y_top=y_range[2]*(1+0.2*length(test_res))
  if(!jitter.only){
    if(violin){ # Violin box plot?
      plot(NULL,xlim=c(0.5,length(input_list)+0.5),
           ylim=c(y_range[1],y_top),xaxt='n',xlab='',
           ylab=yaxis_lab,main=title)
      lapply(1:length(input_list),function(k){vioplot::vioplot(input_list[[k]],at=k,col = colors[k],add=T
                                                               ,drawRect=T)})
    }else{
      boxplot(input_list,ylab=yaxis_lab,col=colors,
              ylim=c(y_range[1],y_top),
              outline=(!jitter),ann=T,main=title,xaxt='n')
    }
  }else{
    jitter=F
    plot(NULL,xlim=c(0.5,length(input_list)+0.5),
         ylim=c(y_range[1],y_top),xaxt='n',
         xlab='',ylab=yaxis_lab,main=title)
    colors=rep(colors,length(input_list)) # To prevent jitter plot without color, increase color vector
    tmp=boxplot(input_list,plot=F)
    tmp=tmp$stats[c(1,5,3,1,5),]
    tryCatch({bxp(list(stats=tmp,n=seq_len(nrow(tmp))),add=T,lty=1,boxlty=0)},error=function(e){})
    l_ply(1:length(input_list),function(x){
      stripchart(add=T,vertical = T,method = 'jitter',
                 pch=16,at=x,input_list[[x]],bg = 'black',col=colors[x])
    })
  }
  axis(side=1,at=1:length(input_list),labels = F)
  # Draw jitter plot
  if(jitter){
    l_ply(1:length(input_list),function(x){
      stripchart(add=T,vertical = T,method = 'jitter',
                 pch=16,at=x,input_list[[x]],bg = 'black')
    })}
  # Draw sigbar
  # Draw clamps for pairs
  label_y_gaps=seq(y_range[2]*1.04,y_top*0.96,length.out = length(test_res))
  test_label=unlist(strsplit(test_m,'.',fixed = T))
  test_label=paste(sep = "-",toupper(test_label[1]),test_label[2])
  for(i in 1:length(test_res)){
    fr_p=unlist(strsplit(names(test_res)[i],"_&&_"))
    to_p=which(names(input_list)==fr_p[2]);fr_p=which(names(input_list)==fr_p[1])
    if(length(label_y_gaps)!=0 & length(fr_p)!=0 & length(to_p)!=0){
      segments(x0=fr_p,x1=to_p,y0=label_y_gaps[i],lwd = 3)
      if(test_res[i]!='NS'){
        text(x=(fr_p+to_p)/2,y=label_y_gaps[i],labels =paste0(test_label,' p=',test_res[i]),pos=3)
      }else{
        text(x=(fr_p+to_p)/2,y=label_y_gaps[i],labels = 'NS',pos=3)
      }
    }
  }
  # Add label
  y_min=par('usr')[3]
  text(names(input_list),srt=-45,xpd=T,x=1:length(input_list),y=y_min,adj=0)
}

#=============================
# Conduct one-sided test.
#=============================
one.sided_test=function(x=NULL,test_method='t.test',ignore.NA=T){
  # Select test-method
  test=match.arg(test_method,c('ks.test','t.test','wilcox.test'))
  message('Selected test method : ',test)
  # Do test
  x1=get(test)(x[[1]],x[[2]],na.rm=ignore.NA,alternative='less')$p.value
  x2=get(test)(x[[1]],x[[2]],na.rm=ignore.NA,alternative='greater')$p.value
  # Get minimum p-value
  x3=min(c(x1,x2))
  dir=which.min(c(x1,x2))
  dir=ifelse(dir==2,'greater','less')
  names(x3)=paste0(names(x)[1],'_',dir)
  # Return the result
  return(x3)
}

#===========================================================
# Modification of fread function in data.table package.
#===========================================================
fread2=function(file = '',stringsAsFactors = F,data.table = F,skip=0){
  check.n.install.lib(lib.name='data.table')
  library(data.table)
  x1=data.table::fread(file = file,
                       stringsAsFactors = stringsAsFactors,
                       data.table = data.table,
                       skip=skip)
  return(x1)
}

#===========================================================
# This function draws a scatter plot with correlation coefficient and p-value.
# @param x : Values in x-axis.
# @param y : Values in y-axis.
# @param xlabel : Label of x-axis.
# @param ylabel : Label of y-axis.
# @param title : Title.
# @param xlim : Range of x-axis.
# @param ylim : Range of y-axis.
# @param pch : Point style.
# @param cex : Point size.
# @param col : Point color.
# @param cex.axis : Size of label sizes.
# @param cex.main : Size of title.
# @param cor.method : Correlation test method. Allowed values are 'spearman' and 'pearson'.
# @param legend.param : List object of legend parameters.
# Default format is list(position='topright',fill=NULL,pch=NULL,legend=NULL,ncol=1).
#===========================================================
scatterPlot=function(x=NULL,y=NULL,xlabel=NULL,ylabel=NULL,title=NULL,
                     cex.lab=1,cex.axis=1,cex.main=1,ylim=NULL,xlim=NULL,
                     pch=19,cex=1,col='black',cor.method='spear',
                     legend.param=list(poistion='topright',fill=NULL,pch=NULL,legend = NULL,ncol=1)){
  # Drawing part
  plot(x,y,xlab=xlabel,ylab=ylabel,ylim=ylim,xlim,xlim,
       main=title,pch=pch,cex=cex,cex.lab=cex.lab,
       cex.axis=cex.axis,cex.main=cex.main,col=col)
  abline(lm(y~x),col='red',lty=2,lwd=5)
  # Calculate correlation coefficient and p-value
  pv=cor.test(x,y,method=cor.method)
  rho=format(c(pv$estimate,pv$p.value),nsmall=3,digits=1)
  cor.m=toupper(cor.method)
  mtext(side=3,adj=1,paste0(cor.m,' cor=',rho[1],' ',cor.m,' p.value=',rho[2]))
  # Add legend
  if(!is.null(legend.param$legend)){
    legend(x=legend.param$position,col = legend.param$fill,legend = legend.param$legend,ncol = legend.param$ncol,
           pch = legend.param$pch)
  }
  return(list('Cor'=rho[1],'P.value'=rho[2],'Test.method'=cor.method))
}

#===========================================================
# This function applies as.numeric function to the provided matrix.
# @param mm (default=NA) : Matrix or data.frame to be processed.
#===========================================================
as.num.mat=function(mm=NA){
  # Get rownames
  x=rownames(mm)
  mm=as.data.frame(mm,stringsAsFactors=F)
  # Change values
  x1=apply(mm,2,function(k){as.numeric(as.character(k))})
  # Give rownames
  rownames(x1)=x
  # Return matrix
  return(x1)
}

#===========================================================
# Align peptides on the PKM sequence.
#===========================================================
align_peptide.on.pkm_isoforms20181010.1812=function(x=unique(pkm.pep$SEQUENCE.X)){
  # matching the isoform from Uniprot SITE, pkm1_seq=P14618-2, pkm2_seq=P14618-1
  pkm1_seq='MSKPHSEAGTAFIQTQQLHAAMADTFLEHMCRLDIDSPPITARNTGIICTIGPASRSVETLKEMIKSGMNVARLNFSHGTHEYHAETIKNVRTATESFASDPILYRPVAVALDTKGPEIRTGLIKGSGTAEVELKKGATLKITLDNAYMEKCDENILWLDYKNICKVVEVGSKIYVDDGLISLQVKQKGADFLVTEVENGGSLGSKKGVNLPGAAVDLPAVSEKDIQDLKFGVEQDVDMVFASFIRKASDVHEVRKVLGEKGKNIKIISKIENHEGVRRFDEILEASDGIMVARGDLGIEIPAEKVFLAQKMMIGRCNRAGKPVICATQMLESMIKKPRPTRAEGSDVANAVLDGADCIMLSGETAKGDYPLEAVRMQHLIAREAEAAMFHRKLFEELVRASSHSTDLMEAMAMGSVEASYKCLAAALIVLTESGRSAHQVARYRPRAPIIAVTRNPQTARQAHLYRGIFPVLCKDPVQEAWAEDVDLRVNFAMNVGKARGFFKKGDVVIVLTGWRPGSGFTNTMRVVPVP'
  pkm2_seq='MSKPHSEAGTAFIQTQQLHAAMADTFLEHMCRLDIDSPPITARNTGIICTIGPASRSVETLKEMIKSGMNVARLNFSHGTHEYHAETIKNVRTATESFASDPILYRPVAVALDTKGPEIRTGLIKGSGTAEVELKKGATLKITLDNAYMEKCDENILWLDYKNICKVVEVGSKIYVDDGLISLQVKQKGADFLVTEVENGGSLGSKKGVNLPGAAVDLPAVSEKDIQDLKFGVEQDVDMVFASFIRKASDVHEVRKVLGEKGKNIKIISKIENHEGVRRFDEILEASDGIMVARGDLGIEIPAEKVFLAQKMMIGRCNRAGKPVICATQMLESMIKKPRPTRAEGSDVANAVLDGADCIMLSGETAKGDYPLEAVRMQHLIAREAEAAIYHLQLFEELRRLAPITSDPTEATAVGAVEASFKCCSGAIIVLTKSGRSAHQVARYRPRAPIIAVTRNPQTARQAHLYRGIFPVLCKDPVQEAWAEDVDLRVNFAMNVGKARGFFKKGDVVIVLTGWRPGSGFTNTMRVVPVP'
  # Find location
  check.n.install(lib.name='stringr')
  library(stringr)
  # Location on PKM1 isoform
  pkm1.location=ldply(1:length(x),function(j){
    x1=as.data.frame.matrix(stringr::str_locate_all(pattern = x[j],pkm1_seq)[[1]])
    if(nrow(x1)!=0){x1$peptide=x[j]}
    return(x1)
  })
  pkm1.location$isoform=paste0('PKM1','(',pkm1.location$start,'_',pkm1.location$end,')')
  # Location on PKM2 isoform
  pkm2.location=ldply(1:length(x),function(j){
    x1=as.data.frame.matrix(stringr::str_locate_all(pattern = x[j],pkm2_seq)[[1]])
    if(nrow(x1)!=0){x1$peptide=x[j]}
    return(x1)
  })
  pkm2.location$isoform=paste0('PKM2','(',pkm2.location$start,'_',pkm2.location$end,')')
  # Merge the tables
  pkm.tab=merge(pkm1.location,pkm2.location,by='peptide',all=T)
  # Combine the isoform tables
  pkm.tab$isoform=paste0(pkm.tab$isoform.x,';',pkm.tab$isoform.y)
  pkm.tab$isoform=gsub(pkm.tab$isoform,pattern = 'NA;',replacement = "")
  pkm.tab$isoform=gsub(pkm.tab$isoform,pattern = ';NA',replacement = "")
  # Select columns
  pkm.tab=pkm.tab[,c('peptide','isoform')]
  # Return the table
  return(pkm.tab)
}

#===========================================================
# Data processing of PKM data.
#===========================================================
make.pkm.matrix_into_proper.form.to.anaylze.pkm.intensity.by.pos20181010.1917=function(
  x=pkm.pep2,
  samples=colnames(pkm.pep2)[12:61],
  pkm.isoform='pkm1'){
  pkm1_seq='MSKPHSEAGTAFIQTQQLHAAMADTFLEHMCRLDIDSPPITARNTGIICTIGPASRSVETLKEMIKSGMNVARLNFSHGTHEYHAETIKNVRTATESFASDPILYRPVAVALDTKGPEIRTGLIKGSGTAEVELKKGATLKITLDNAYMEKCDENILWLDYKNICKVVEVGSKIYVDDGLISLQVKQKGADFLVTEVENGGSLGSKKGVNLPGAAVDLPAVSEKDIQDLKFGVEQDVDMVFASFIRKASDVHEVRKVLGEKGKNIKIISKIENHEGVRRFDEILEASDGIMVARGDLGIEIPAEKVFLAQKMMIGRCNRAGKPVICATQMLESMIKKPRPTRAEGSDVANAVLDGADCIMLSGETAKGDYPLEAVRMQHLIAREAEAAMFHRKLFEELVRASSHSTDLMEAMAMGSVEASYKCLAAALIVLTESGRSAHQVARYRPRAPIIAVTRNPQTARQAHLYRGIFPVLCKDPVQEAWAEDVDLRVNFAMNVGKARGFFKKGDVVIVLTGWRPGSGFTNTMRVVPVP'
  pkm2_seq='MSKPHSEAGTAFIQTQQLHAAMADTFLEHMCRLDIDSPPITARNTGIICTIGPASRSVETLKEMIKSGMNVARLNFSHGTHEYHAETIKNVRTATESFASDPILYRPVAVALDTKGPEIRTGLIKGSGTAEVELKKGATLKITLDNAYMEKCDENILWLDYKNICKVVEVGSKIYVDDGLISLQVKQKGADFLVTEVENGGSLGSKKGVNLPGAAVDLPAVSEKDIQDLKFGVEQDVDMVFASFIRKASDVHEVRKVLGEKGKNIKIISKIENHEGVRRFDEILEASDGIMVARGDLGIEIPAEKVFLAQKMMIGRCNRAGKPVICATQMLESMIKKPRPTRAEGSDVANAVLDGADCIMLSGETAKGDYPLEAVRMQHLIAREAEAAIYHLQLFEELRRLAPITSDPTEATAVGAVEASFKCCSGAIIVLTKSGRSAHQVARYRPRAPIIAVTRNPQTARQAHLYRGIFPVLCKDPVQEAWAEDVDLRVNFAMNVGKARGFFKKGDVVIVLTGWRPGSGFTNTMRVVPVP'
  # Get length of pkm isoforms
  pkm1.length=nchar(pkm1_seq)
  pkm2.length=nchar(pkm2_seq)
  # Make expression matrix into list according to samples
  x1=llply(1:length(samples),function(j){
    wh=which(colnames(x) %in% c('isoform',samples[j]))
    x1=x[,c(wh)]
    # Make identifier and reduce redundancy
    x1$dup=paste0(x1[,1],x1[,2])
    x1=x1[which(!duplicated(x1$dup)),]
    x1=x1[,c(1,2)]
    return(x1)
  })
  names(x1)=samples
  # Create matrix (number of peptides x amino acid length) for pkm1 and pkm2
  if(pkm.isoform=='pkm1'){pkm_seq=pkm1_seq}else{pkm_seq=pkm2_seq}
  x2=llply(1:length(samples),function(j){
    tmp=x1[[j]]
    # Select PKM isoform
    tmp=tmp[grep(tmp$isoform,pattern = pkm.isoform,ignore.case = T),]
    # Remove other pkm isoform information from the table
    for(k in 1:nrow(tmp)){
      tmp2=unlist(strsplit(tmp$isoform[k],";"))
      tmp2=tmp2[grep(tmp2,pattern = pkm.isoform,ignore.case = T)]
      tmp$isoform[k]=tmp2
    }
    # Add 'start' and 'end' columns
    tmp$start=NA;tmp$end=NA
    for(k in 1:nrow(tmp)){
      tmp2=unlist(strsplit(tmp$isoform[k],'(',fixed = T))[2]
      tmp2=unlist(strsplit(tmp2,'_',fixed = T))
      tmp2=gsub(tmp2,pattern = ')',replacement = '',fixed = T)
      tmp2=as.numeric(tmp2)
      tmp[k,3:4]=tmp2
    }
    # Create matrix form for pkm analysis
    mat=matrix(NA,ncol=nchar(pkm_seq),nrow=nrow(tmp))
    for(k in 1:nrow(tmp)){
      st=tmp$start[k];end=tmp$end[k]
      mat[k,st:end]=tmp[k,2]
    }
    rownames(mat)=tmp$isoform
    colnames(mat)=unlist(strsplit(pkm_seq,''))
    return(mat)
  })
  names(x2)=names(x1)
  # Export the tables
  return(list(listed.table=x1,listed.parsed.table.for.pkm.analysis=x2,type=pkm.isoform))
}

#===========================================================
# This function retrieves peptide sequence of PKM using biomaRT package.
#===========================================================
getting_ensg.id_and_pkm.sequence20180816.1410=function(genes=c('PKM')){
  check.n.install.lib('biomaRt',lib.type = 'bioc')
  library(biomaRt)
  if(!exists('ensembl')){
    ensembl <<- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37, version=75)
  }
  
  human.prot = getBM(attributes = c("peptide",'ensembl_gene_id','hgnc_symbol','ensembl_transcript_id'),values = genes,mart=ensembl,filters = 'hgnc_symbol')
  # Remove ENSG ids and ENST ids without peptide sequences
  human.prot=human.prot[human.prot$peptide!="Sequence unavailable",]
  # Remove asterisks
  human.prot$peptide=gsub(human.prot$peptide,pattern = '*',fixed = T,replacement = "")
  human.prot$peptide_length=unlist(llply(1:nrow(human.prot),function(x){nchar(human.prot$peptide[x])}))
  # match the isoform from Uniprot SITE, pkm1_seq=P14618-2, pkm2_seq=P14618-1
  pkm1_seq='MSKPHSEAGTAFIQTQQLHAAMADTFLEHMCRLDIDSPPITARNTGIICTIGPASRSVETLKEMIKSGMNVARLNFSHGTHEYHAETIKNVRTATESFASDPILYRPVAVALDTKGPEIRTGLIKGSGTAEVELKKGATLKITLDNAYMEKCDENILWLDYKNICKVVEVGSKIYVDDGLISLQVKQKGADFLVTEVENGGSLGSKKGVNLPGAAVDLPAVSEKDIQDLKFGVEQDVDMVFASFIRKASDVHEVRKVLGEKGKNIKIISKIENHEGVRRFDEILEASDGIMVARGDLGIEIPAEKVFLAQKMMIGRCNRAGKPVICATQMLESMIKKPRPTRAEGSDVANAVLDGADCIMLSGETAKGDYPLEAVRMQHLIAREAEAAMFHRKLFEELVRASSHSTDLMEAMAMGSVEASYKCLAAALIVLTESGRSAHQVARYRPRAPIIAVTRNPQTARQAHLYRGIFPVLCKDPVQEAWAEDVDLRVNFAMNVGKARGFFKKGDVVIVLTGWRPGSGFTNTMRVVPVP'
  pkm2_seq='MSKPHSEAGTAFIQTQQLHAAMADTFLEHMCRLDIDSPPITARNTGIICTIGPASRSVETLKEMIKSGMNVARLNFSHGTHEYHAETIKNVRTATESFASDPILYRPVAVALDTKGPEIRTGLIKGSGTAEVELKKGATLKITLDNAYMEKCDENILWLDYKNICKVVEVGSKIYVDDGLISLQVKQKGADFLVTEVENGGSLGSKKGVNLPGAAVDLPAVSEKDIQDLKFGVEQDVDMVFASFIRKASDVHEVRKVLGEKGKNIKIISKIENHEGVRRFDEILEASDGIMVARGDLGIEIPAEKVFLAQKMMIGRCNRAGKPVICATQMLESMIKKPRPTRAEGSDVANAVLDGADCIMLSGETAKGDYPLEAVRMQHLIAREAEAAIYHLQLFEELRRLAPITSDPTEATAVGAVEASFKCCSGAIIVLTKSGRSAHQVARYRPRAPIIAVTRNPQTARQAHLYRGIFPVLCKDPVQEAWAEDVDLRVNFAMNVGKARGFFKKGDVVIVLTGWRPGSGFTNTMRVVPVP'
  human.prot$pkm_isoform=NA
  human.prot$pkm_isoform[which(human.prot$peptide==pkm1_seq)]='PKM1'
  human.prot$pkm_isoform[which(human.prot$peptide==pkm2_seq)]='PKM2'
  # Export result
  return(human.prot)
}

#===========================================================
# This function determines the optimal cut-off point which gives maximum survival difference between expression low and high groups
# The test method is log-rank test.
# @param x (default=NULL) : Gene expression table, rows are genes, columns are samples
# @param y (default=NULL) : Survival matrix. Column names should be 'sample', 'surv', and 'status'.
# @param bias (default=NULL) : Whether this function should locate survival point based on user provided direction
# options are 'good', 'bad'. Default (NULL) provides unbiased result.
#===========================================================
get_optExp_survPoint<-function(x=NULL,y=NULL,bias=NULL){
  check.n.install.lib(lib.name = c('plyr','matrixStats','survminer'))
  library(plyr);library(matrixStats)
  # Remove 0 exp genes
  zero.gene=matrixStats::rowMeans2(as.matrix(x))
  x1=x[which(zero.gene!=0),]
  # Remove samples without survival data
  x1=x1[,which(colnames(x1) %in% y$sample)]
  y1=y[which(y$sample %in% colnames(x1)),]
  # sort the tables
  x1=x1[,order(colnames(x1))]
  y1=y1[order(y$sample),]
  # Conduct log-rank test. Assign at least 5 samples per group #k=1
  x2=plyr::ldply(1:nrow(x),.progress='text',function(k){
    df=data.frame(exp=unlist(x1[k,]),os=y1$surv,stat=y1$status,stringsAsFactors=F)
    df$exp=as.numeric(as.character(df$exp))
    df=df[order(df$exp,decreasing=F),]
    df$group='High'
    k1=plyr::ldply(5:(nrow(df)-5),function(j){
      df$group[1:j]='Low'
      j1=survival::Surv(df$os,df$stat)~df$group
      j2=survival::survfit(j1)
      m.surv=survminer::surv_median(j2)
      os=survival::survdiff(j1)
      pv=pchisq(os$chisq,length(os$n)-1,lower.tail = F)
      harzard.ratio=survival::coxph(j1)
      harzard.ratio2=exp(harzard.ratio$coefficients)
      if(grepl(names(harzard.ratio2),pattern = 'Low')){
        res=ifelse(harzard.ratio2>1,'Favorable','Unfavorable')
        res2=1/harzard.ratio2
      }else{
        res=ifelse(harzard.ratio2<1,'Favorable','Unfavorable')
        res2=harzard.ratio2
      }
      res3=c(pv,res,df$exp[j],res2)
      return(res3)
    })
    k1=as.data.frame(k1,stringsAsFactors=F)
    k1$V1=as.numeric(k1$V1)
    k1$V2=as.numeric(k1$V2)
    # Whether survival should consider user-wanted direction
    if(is.null(bias)){
      k2=k1[which.min(k1[,1]),]
    }else{
      if(bias=='good'){
        k2=k1[which(k1[,2]=='Favorable'),]
        k2=k2[which.min(k2[,1]),]
      }else{
        k2=k1[which(k1[,2]=="Unfavorable"),]
        k2=k2[which.min(k2[,1]),]
      }
    }
    if(nrow(k2)==0){
      k2=data.frame(matrix(rep(NA, 4), nrow = 1, ncol = 4), stringsAsFactors = F)
    }
    colnames(k2)=c('p.value','prognosis','cut_off','Hazard.ratio')
    k2$gene<-rownames(x)[k]
    return(k2)
  })
  x2$cut_off=as.numeric(x2$cut_off)
  # Conduct cox regression test # head(x[,1:3]);head(y)
  a1=t(x)
  a2=plyr::ldply(1:ncol(a1),function(k){ #k=1
    k1=survival::Surv(y$surv,y$status) ~ a1[,k]
    k2=survival::coxph(k1)
    k3=summary(k2)
    k4=coef(k3)
    k4.1=k4[2] # HR
    k4.2=k3$waldtest[3]
    k4.3=k3$logtest[3]
    k4.4=k3$sctest[3]
    c(k4.1,k4.2,k4.3,k4.4)
  })
  colnames(a2)=c('Hazard.ratio_coxReg','Cox_Wald.test','Cox_Loglikelihood.test','Cox_Score.test')
  a2$Cox_prognosis=ifelse(a2$Hazard.ratio_coxReg>1,'Unfavorable','Favorable')
  a2$gene=colnames(a1)
  
  # Combine matrix
  x3=mat.merge(x2,a2,by='gene')
  # Return the result
  return(x3)
}

#===========================================================
# This function plots a density plot with colors.
# @param input_list (default=NULL) : Input list object with names.
# @param colors (default=NULL) : Density plot color. Provide it as vector.
# @param xaxis_lab (default=NULL) : Label of x-axis.
# @param alpha (default=NULL) : See ?adjustcolor.
# @param title (default=NULL) : Title of the plot.
# @param legend.pos (default='topleft') : Decides where legend box should be placed.
#===========================================================
density_plot=function(input_list=NULL,xaxis_lab=NULL,alpha=0.5,
                      colors=NULL,title=NULL,legend.pos='topleft'){
  if(is.null(input_list)){stop('Input data is empty. Terminate process')}
  
  # Calculate densities to get plot size
  x=input_list
  x1=llply(x,function(k){
    k1=density(k)
    k1=k1[c('x','y')]
    return(k1)
  })
  # Get x axis range
  x.range=llply(x1,function(k){
    k1=k$x
    k2=range(k1)
    return(k2)
  })
  x.range2=range(unlist(x.range))
  # Get y axis range
  y.range=llply(x1,function(k){
    k1=k$y
    k2=range(k1)
    return(k2)
  })
  y.range2=range(unlist(y.range))
  # Generate plot
  plot(NULL,xlim=x.range2,ylim=y.range2,xlab=xaxis_lab,ylab='Density',main=title)
  # Modify color
  check.n.install(lib.name='grDevices')
  library(grDevices)
  colors2=grDevices::adjustcolor(colors,alpha.f = alpha)
  # Draw the plots with colors #i=1
  for(i in 1:length(x)){
    i1=x1[[i]]
    graphics::polygon(i1$x,i1$y,col = colors2[i])
  }
  # Add legend
  legend(legend.pos,fill=colors2,legend=names(x))
}

#===========================================================
# This function draws a bar plot with error bars. 
# @param input_list (default=NULL) : Input list object with names.
# @param yaxis_lab (default=NULL) : Label of y-axis.
# @param title (default=NULL) : Title of the plot.
#===========================================================
sig_barplot=function(input_list=NULL,yaxis_lab=NULL,xaxis_lab=NULL,
                     title=NULL){
  if(is.null(input_list)){stop('Input data is empty. Terminate process')}
  #================================
  # Draw bar plot using GGPlot2
  #================================
  check.n.install(lib.name=c('ggplot2','ggsignif'),lib.type=c('cran','cran'))
  library(ggplot2)
  library(ggsignif)
  # Make list into data.
  x=input_list
  x1=ldply(1:length(x),function(k){ #k=1
    k1=x[[k]]
    k2=data.frame(Group=names(x)[k],Value=mean(k1,na.rm=T),Std=sd(k1,na.rm = T))
    return(k2)
  })
  x1[,-1]=as.num.mat(x1[,-1])
  # Bar plot main
  # Y axis range
  x2=ggplot(x1,aes(Group,Value))+
    geom_bar(aes(fill = Group),stat='identity', position = 'dodge',width = .5)+
    geom_errorbar(aes(ymin=Value-Std,ymax=Value+Std),width=0.2,
                  position = position_dodge(0.9))+
    theme_classic()
  
  # Y axis label
  x3=x2 +labs(title = title,x=xaxis_lab,y=yaxis_lab)
  return(x3)
}

#===========================================================
# Figure 4h data parsing.
#===========================================================
fig4h.data.parsing20200408.1729=function(x=fig4h.data){
  # Split cells
  x1=x[-1,-1]
  rownames(x1)=1:nrow(x1)
  colnames(x1)=paste0('V',1:ncol(x1))
  hs683=x1[1:16,]
  snu1105=x1[23:38,]
  rownames(hs683)=rownames(snu1105)=1:nrow(hs683)
  
  # Split data
  # HS683
  hs683.activity=hs683[1:2,3:8]
  hs683.activity2=hs683.activity[-1,]
  colnames(hs683.activity2)=hs683.activity[1,]
  hs683.distance=hs683[4:nrow(hs683),3:8]
  hs683.distance2=hs683.distance[-1,]
  colnames(hs683.distance2)=hs683.distance[1,]
  # SNU1105
  snu1105.activity=snu1105[1:2,3:8]
  snu1105.activity2=snu1105.activity[-1,]
  colnames(snu1105.activity2)=snu1105.activity[1,]
  snu1105.distance=snu1105[4:nrow(snu1105),3:8]
  snu1105.distance2=snu1105.distance[-1,]
  colnames(snu1105.distance2)=snu1105.distance[1,]
  # return the result
  hs683.datas=list('activity'=hs683.activity2,'invasion'=hs683.distance2)
  snu1105.datas=list('activity'=snu1105.activity2,'invasion'=snu1105.distance2)
  datas=list('HS683'=hs683.datas,'SNU1105'=snu1105.datas)
  return(datas)
}

#===========================================================
# Extended Figure 4d data parsing.
#===========================================================
Extended_Figure4d.data.parsing20200408.1615=function(x=exfig4d.data){
  # Split cells
    x1.snu201=x[1:14,]
    x1.kns81=x[21:34,]
  # SNU201 data parsing
    x2.snu201=x1.snu201[-c(1),]
    colnames(x2.snu201)=x1.snu201[1,]
  # KNS81 data parsing
    x2.kns81=x1.kns81[-c(1:2),]
    rownames(x2.kns81)=1:nrow(x2.kns81)
    colnames(x2.kns81)=x1.kns81[2,]
  # Export tables
  return(list('SNU201'=x2.snu201,'KNS81'=x2.kns81))
}

#===========================================================
# This function draws a bar plot with annotation.
# @param mm (default=NA) : Input matrix. Rows are categories and columns are samples
# @param color (default=NULL) : Color for categories
# @param color.sample (default=NULL) : Color vector for samples
# @param color.sample.legend : Color legend for samples
# @param legend.sample (default=NULL) : Legends for sample annotation
# @param beside (default=F) : Whether barplot should be stacked (F) or not(T)
# @param bar.legend.ncol (default=2) : The number of legend column
# @param ylab  (default=NULL) : Y-axis label
# @param xlab  (default=NULL) : X-axis label
# @param title  (default=NULL) : Title
# @param xaxt (default=NULL) : Whether sample label should be written? ('n' for no)
# @param cex.lab (default=1.5) : Size of label
# @param cex.axis (default=1.5) : Size of axis label
# @param las (default=2) : Vertical (1) or horizontal(2)?
# @keywords barplot
# @export
# @examples
# mm=as.data.frame(matrix(sample(1:100,100),5,20),stringsAsFactors=F)
# rownames(mm)=LETTERS[1:5]
# mm=apply(mm,2,function(k){k/sum(k)})
# title=xlab=ylab='example';beside=F
# color=rainbow(5)
# color.sample=sample(rainbow(5),20,replace=T)
# names(color.sample)=LETTERS[1:length(color.sample)]
# annotated_barplot(mm=mm,color=color,color.sample=color.sample,ylab=ylab,xlab=xlab,title=title,
# beside=F,bar.legend.ncol=2)
#===========================================================
annotated_barplot=function(mm=NA,color,color.sample,color.sample.legend,xlab,ylab,title=NULL,beside=F,
                           bar.legend.ncol=2,xaxt='n',las=2,cex.lab=1.5,cex.axis=1.5){
  # Get categories
  x.coords=barplot(mm,plot = F)
  # Put annotation on the barplot
  ytop=max(colSums(mm))
  barplot(mm,plot = T,col = color,ylim=c(0,ytop*1.4),ylab=ylab,main=title,xlab=xlab,
          beside = beside,
          xaxt=xaxt,las=las,cex.lab=cex.lab,cex.axis=cex.axis)
  legend('bottomright',legend = rownames(mm),fill = color,ncol = bar.legend.ncol,title = 'Barplot legend')
  points(x.coords,rep(ytop*1.05,length(x.coords)),bg=color.sample,pch=22)
  # Put legend
  legend('top',legend=names(color.sample.legend),fill=color.sample.legend,
         ncol=length(color.sample.legend),
         title = 'Sample legend')
}

#===========================================================
# single sample Gene Set Enrichment Analysis (ssGSEA)
# Source code from https://gist.github.com/gaoce/39e0907146c752c127728ad74e123b33
# @param mm matrix. Rows are genes. Columns are samples. Row names are symbols.
# @param gene_sets list. Each element is a string vector with gene symbols.
# @param alpha numeric. Parameter for ssGSEA, the default is 0.25
# @param scale logical. If True, normalize the scores by number of genes in the gene sets.
# @param norm logical. If True, normalize the scores by the absolute difference between max and min values.
# @param single logical. If True, use ssGSEA algorithm, otherwise use GSEA.
#===========================================================
ssgsea = function(mm, gene_sets, alpha = 0.25, scale = T, norm = F, single = T) {
  row_names = rownames(mm)
  num_genes = nrow(mm)
  gene_sets = lapply(gene_sets, function(genes) {which(row_names %in% genes)})
  
  # Ranks for genes
  R = matrixStats::colRanks(mm, preserveShape = T, ties.method = 'average')
  
  # Calculate enrichment score (es) for each sample (column)
  es = apply(R, 2, function(R_col) {
    gene_ranks = order(R_col, decreasing = TRUE)
    
    # Calc es for each gene set
    es_sample = sapply(gene_sets, function(gene_set_idx) {
      # pos: match (within the gene set)
      # neg: non-match (outside the gene set)
      indicator_pos = gene_ranks %in% gene_set_idx
      indicator_neg = !indicator_pos
      
      rank_alpha  = (R_col[gene_ranks] * indicator_pos) ^ alpha
      
      step_cdf_pos = cumsum(rank_alpha)    / sum(rank_alpha)
      step_cdf_neg = cumsum(indicator_neg) / sum(indicator_neg)
      
      step_cdf_diff = step_cdf_pos - step_cdf_neg
      
      # Normalize by gene number
      if (scale) step_cdf_diff = step_cdf_diff / num_genes
      
      # Use ssGSEA or not
      if (single) {
        sum(step_cdf_diff)
      } else {
        step_cdf_diff[which.max(abs(step_cdf_diff))]
      }
    })
    unlist(es_sample)
  })
  
  if (length(gene_sets) == 1) es = matrix(es, nrow = 1)
  
  # Normalize by absolute diff between max and min
  if (norm) es = es / diff(range(es))
  
  # Prepare output
  rownames(es) = names(gene_sets)
  colnames(es) = colnames(mm)
  return(es)
}

#===========================================================
# VAF filter function for Fig 1a.
#===========================================================
vaf_filter20180918.1050=function(x=smc1.snv,vaf=0.05){
  #message(paste0('Export VCF above VAF=',vaf))
  x1=llply(x,.progress = 'text',function(j){
    j$vaf=as.numeric(j$vaf)
    j=j[j$vaf>vaf,]
    return(j)
  })
  return(x1)
}

#===========================================
# This function calculates sample-pairwise jaccard coefficient.
#===========================================
calc_jaccard.coef20180918.1113=function(x=x1){
  x1=llply(x,.progress='text',function(j){
    j$feature=paste0(j[,1],'@',j[,2],'@',j[,4],'@',j[,5])
    return(j$feature)
  })
  # Create jaccard coef matrix
  n_sample=length(x1)
  x2=matrix(0,dimnames = list(names(x1),names(x1)),ncol=n_sample,nrow=n_sample)
  # Calculate and add Jaccard coefficient
  message('Calc jc.coef')
  for(i in 1:nrow(x2)){for(j in 1:ncol(x2)){
    jc=length(intersect(x1[[i]],x1[[j]]))/length(union(x1[[i]],x1[[j]]))
    x2[i,j]=jc
  }}
  # 1-jc.mat=distance
  x3=as.dist(1-x2)
  return(list(jc.coef=x2,dist.mat=x3))
}

#===========================================================
# This function conducts Gene Enrichment Analysis similarly to the DAVID site.
# P-value will be calculated by one-sided hypergeometric test.
# phyper(N_gene-1,N_all_included,N_backgroud-N_all_included,N_gene,lower.tail=T)
# @param input.gene : Gene vector to be examined.
# @param backgroud.gene : Background gene vector to be examined.
# @param geneset.list : Geneset list. Its format should be list.
# example.list=list(Aset=c('GENEA','GENEB','GENEC','GENED'),
# proteasome=c('GENEA','GENEB','GENEF'))
#===========================================================
Gene.Enrichment.Analysis=function(input.gene,background.gene,geneSet.list){
  message('Conducting Gene.Enrichment.Analysis')
  results=plyr::ldply(1:length(geneSet.list),function(k){
    # Add numbers to the contingency table
    # input.gene
    hit.genes=paste(collapse=',',intersect(input.gene,geneSet.list[[k]]))
    in.input=length(intersect(input.gene,geneSet.list[[k]]))
    out.input=length(setdiff(input.gene,geneSet.list[[k]]))
    # Backgroud.Gene
    in.back=length(intersect(background.gene,geneSet.list[[k]]))
    out.back=length(setdiff(background.gene,geneSet.list[[k]]))
    if(T){
      #P-value will be calculated by one-sided hypergeometric test.
      #phyper(N_gene-1,N_all_included,N_backgroud-N_all_included,N_gene,lower.tail=T)
      p.val=phyper(in.input-1,length(input.gene),length(background.gene),in.input+in.back,lower.tail = F)
    }
    # Return the results
    results=c(p.val,in.input,out.input,in.back,out.back,hit.genes)
    return(results)
  })
  names(results)=c('Pvalue','In','Out','In.Background','Out.Backgroud','Hits')
  # Get Enriched cell only
  results=as.data.frame.matrix(results,stringsAsFactors=F)
  results$Geneset.Size=unlist(llply(geneSet.list,function(k){{length(k)}}))
  results$FDR=p.adjust(as.numeric(results$Pvalue),'fdr')
  results$GeneSet=names(geneSet.list)
  # Sort the table
  results=results[,c('GeneSet','Geneset.Size','Pvalue','FDR','In','Out',
                     'In.Background','Out.Backgroud','Hits')]
  # Export the table
  return(results)
}

#=========================================================== 
# This function performs statistical tests (ks.test, t.test, wilcox.test) on two groups in the provided matrix.
# @param mm (default=NULL) : Input matrix (gene x sample). 
# It is not necssary for the number of columns to be identical with the number of samples.
# @param groups (default=list) : List of two groups. Do not provide more than two groups.
# ex) list('GPC1'=c('sampleA','sampleB','sampleC'),'GPC2'=c('sampleX','sampleY','sampleZ'))
# @param feature (default=NULL) : Name of feature in the matrix. (such as 'symbol','gene') 
# If it is not provided, rownames will be used as name of feature.
# @param test_method (default='t.test') : Test-method. Allowed values are 'ks.test', 't.test', 'wilcox.test'.
# @param one.sided (default=F) : Whether test should be conducted as one-sided.
# @param ignore.NA (default=F) : If true, average of group will be calculated regardless of NA presence.
#===========================================================
matrix_test=function(mm=NULL,groups=NULL,feature=NULL,test_method='t.test',
                     ignore.NA=F,one.sided=F){
  if(is.null(mm)){warnings('Provide matrix');break}
  if(is.null(groups)){warnings('Provide group list');break}
  # Convert matrix into data.frame.matrix
  mm=as.data.frame.matrix(mm)
  # Test
  # Select test method
  test=match.arg(test_method,c('ks.test','t.test','wilcox.test'))
  message('Selected test method : ',test)
  # Test
  wh.x=colnames(mm) %in% groups[[1]]
  wh.y=colnames(mm) %in% groups[[2]]
  if(one.sided){message(paste0('One-sided ',test_method))}
  res=ldply(1:nrow(mm),.progress = 'text',function(i){
    x1=as.numeric(mm[i,wh.x]);x2=as.numeric(mm[i,wh.y])
    tryCatch({
      if(one.sided){
        pv=suppressMessages(one.sided_test(list(x1,x2),test_method = test,ignore.NA = ignore.NA)) 
      }else{
        pv=get(test)(x1,x2,na.rm=ignore.NA,)$p.value
      }
      return(c(pv,mean(x1,na.rm=ignore.NA),mean(x2,na.rm=ignore.NA),mean(c(x1,x2),na.rm=ignore.NA)))
    },error=function(e){return(rep(NA,4))})
  })
  colnames(res)=c(paste0(test,'_pvalue'),paste0(names(groups),'_mean'),'All_mean')
  if(is.null(feature)){res$feature=rownames(mm)}else{res$feature=mm[,feature]}
  res$fdr=p.adjust(res[,1],'fdr')
  # Order table
  res=res[,c('feature',paste0(test,'_pvalue'),'fdr',paste0(names(groups),'_mean'),'All_mean')]
  colnames(res)[3]=paste0(test,'_FDR')
  return(res)
}

# Writter : Sejin Oh
# English edditing : Hwanho Lee, Myung Joon Oh