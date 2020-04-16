#===========================
# Load data
#===========================
load.data20200226.0011=function(wd=wd){
  setwd(wd)
  smc1.info=read.xl('Supplementary Data 1.xlsx',sheet = 2)
  smc1.info<<-smc1.info
  smc1.prot=read.xl('Supplementary Data 2.xlsx',sheet = 2)
  colnames(smc1.prot)=gsub(colnames(smc1.prot),pattern = ' ',replacement = '')
  smc1.prot<<-smc1.prot
  # Load SMC1 data
  load('SMC1-RNAseq-data-for-Fig1d-3b-4c-6b-6d-FigS1b-S3b-S6a.Rdata',envir = .GlobalEnv)
  # Phospho datasets
    smc1.phos=read.xl('Supplementary Data 2.xlsx',sheet = 3)
    smc1.phos<<-smc1.phos
  # Drug response
    smc1.auc=read.xl('Supplementary Data 6.xlsx',sheet = 3)
    smc1.ed=read.xl('Supplementary Data 6.xlsx',sheet = 4)
  for(i in 2:ncol(smc1.auc)){
    smc1.auc[,i]=gsub(smc1.auc[,i],pattern = "'",replacement = '')
    smc1.ed[,i]=gsub(smc1.ed[,i],pattern = "'",replacement = '')
  }
  smc1.auc[,-1]=as.num.mat(smc1.auc[,-1])
  smc1.auc<<-smc1.auc
  smc1.ed[,-1]=as.num.mat(smc1.ed[,-1])
  smc1.ed<<-smc1.ed
  drug.info=read.xl('Supplementary Data 6.xlsx',sheet = 2)
  drug.info<<-drug.info
  # Load TCGA data
  load('TCGA-cohort-data-RNAseq-for-Fig4c-4e.datas.Rdata',envir = .GlobalEnv) # GBM clinical and expression data
  # download msigDB version 6.2 and corum v3.0
  # We made genesets into list object
   load('CORUM.v3.0_geneset-for-Fig3a.Rdata',envir = .GlobalEnv)
   corum=corum.v3.0.proteinComplex
   corum<<-corum
   load('msigDB.v6.2-for-Fig3a-3b-6d-FigS1b-S1e-S3a.Rdata')
   msig=updated.symbol.msigDB.v6.2_20190225
   msig<<-msig
 # download CCLE transcriptome and metabolome data
   #https://portals.broadinstitute.org/ccle/data
   # CCLE_metabolomics_20190502.csv
   # CCLE_RNAseq_genes_rpkm_20180929.gct.gz
   # CCLE_DepMap_18Q1_maf_20180207.txt
    load('CCLE-data-for-Fig3e.Rdata',envir = .GlobalEnv)
 # download Darmanis data
   # I made the datas into Rdata format to facilitate loading process.
   load('Darmanis_meta_info-for-Fig5.Rdata',envir = .GlobalEnv)
   # Download and load Darmanis single cell RNAseq data
   message('Download load Single-RNA-seq data of Daramanis into R')
  # BRCANESS formula PMC2917311
   brcaness.formula=read.xl('BRCAness_formula_PMC2917311.xlsx',sheet = 1)
   brcaness.formula<<- brcaness.formula
# Microarray datas
   load('ANOCEF-cohort-data-microarray-for-Fig4c-FigS5a.Rdata',envir = .GlobalEnv)
   load('Yonsei-cohort-data-microarray-for-Fig4c-FigS5a.Rdata',envir = .GlobalEnv)
 # Peptide intensities
   set1.peptide=fread2('SMC1-Set1-Globalpeptide-quantification-data-for-Fig3d.txt')
   set2.peptide=fread2('SMC1-Set2-Globalpeptide-quantification-data-for-Fig3d.txt')
   set1.peptide<<-set1.peptide
   set2.peptide<<-set2.peptide
  
# SMC-TMA info
   smc.tma=read.xl('Supplementary Data 1.xlsx',sheet = 4)
   smc.tma<<-smc.tma
# Load SMC2data
   load('SMC2-RNAseq-data-for-Fig4c.Rdata',envir=.GlobalEnv)
# Load SMC-TMA DATA
   load('SMC-TMA-cohort-data-for-Fig5g-FigS5b.Rdata',envir = .GlobalEnv)
# Load SMC1.cnv data
   load('SMC1-CNV-data-for-FigS1a-FigS6a.Rdata',envir = .GlobalEnv)
# Figure 4h enzyme and invasion distance
   fig4h.data=as.data.frame(readxl::read_xlsx('Figure4h_source_data.xlsx',sheet = 1),stringsAsFactors=F)
   fig4h.data2=fig4h.data.parsing20200408.1729(x=fig4h.data)
   fig4h.data<<-fig4h.data2
# Extended data Figure 4c enzyme activity data
   exfig4c.data=read.xl('Extended_data_figure_4c.xls',sheet = 1)
   exfig4c.data<<-exfig4c.data
# Extended data figure 4d
   exfig4d.data=read.xl('Extended_data_figure_4d.xlsx',sheet = 1)
   exfig4d.data2=Extended_Figure4d.data.parsing20200408.1615(x=exfig4d.data)
   exfig4d.data<<-exfig4d.data2
# SMC1-SNV
   load('SMC1-SNV-data-for-Fig1a-1c.Rdata',envir = .GlobalEnv)
# SMC1-SNV.by.mass
   load('SMC1-SNV.by.mass-data-for-Fig1c.Rdata',envir = .GlobalEnv)
}

#===========================
# Figure 1a panel
#===========================
#x=smc1.snv;y=smc1.info
Figure1a20200409.1046=function(x=smc1.snv,y=smc1.info){
  check.n.install.lib(c('plyr','pheatmap'),lib.type = c('cran','cran'))
  library(plyr);library(pheatmap)
  # Get VAF>5% SNVs
    x1=vaf_filter20180918.1050(x=x,vaf=0.05)
  # Calculate Jaccard coefficient
    x2=calc_jaccard.coef20180918.1113(x=x1)
  # Make annotation tables
    sort(colnames(y))
    y$Multi_samples=gsub(y$Sample_name,pattern = "'",replacement = '')
    y$Multi_samples=unlist(llply(y$Multi_samples,function(k){
      unlist(strsplit(k,'-'))[1]
    }))
    annot=y[,rev(c("NF1_mut","IDH1_mut","EGFR_mut","EGFRvIII",
               "TP53_mut","RB1_mut","PTEN_mut","PIK3CA_mut",
               "PIK3R1_mut","1p/19q codel","RNA_subtype4",
               "Grade","Initial_or_Recurred_sample",'Multi_samples',
               "Features"))]
    rownames(annot)=gsub(y$Sample_name,pattern = "'",replacement = '')
    colnames(annot)=multi_sub(colnames(annot),pattern = c("Initial_or_Recurred_sample"),
                              replacement = c('Primary/Relapse'),
                              exact = T)
    # Parse the annotation data
    # Initial recur
    annot$`Primary/Relapse`[annot$`Primary/Relapse`=='Secondary']='R'
    # Mutation hotspot
    # NF1
    annot$NF1_mut=ifelse(is.na(annot$NF1_mut),'WT','MUT')
    # PTEN
    annot$PTEN_mut=ifelse(is.na(annot$PTEN_mut),'WT','MUT')
    annot$PTEN_mut[y$PTEN_mut_hotspot!='']='HOT'
    # EGFR
    annot$EGFR_mut=ifelse(is.na(annot$EGFR_mut),'WT','MUT')
    annot$EGFR_mut[y$EGFR_mut_hotspot!=""]='HOT'
    # TP53
    annot$TP53_mut=ifelse(is.na(annot$TP53_mut),'WT','MUT')
    annot$TP53_mut[y$TP53_mut_hotspot!='']='HOT'
    # RB1
    annot$RB1_mut=ifelse(is.na(annot$RB1_mut),'WT','MUT')
    # PIK3R1
    annot$PIK3R1_mut=ifelse(is.na(annot$PIK3R1_mut),'WT','MUT')
    # PIK3CA
    annot$PIK3CA_mut=ifelse(is.na(annot$PIK3CA_mut),'WT','MUT')
    annot$PIK3CA_mut[y$PIK3CA_mut_hotspot!='']='HOT'
    # IDH1
    annot$IDH1_mut=ifelse(is.na(annot$IDH1_mut),'WT','MUT')
    annot$IDH1_mut[y$IDH1_mut_hotspot!='']='HOT'
    # EGFR vIII
    annot$EGFRvIII[which(is.na(annot$EGFRvIII))]='WT'
    annot$`1p/19q codel`[which(is.na(annot$`1p/19q codel`))]='WT'
    # Duplicated samples
    dup=na.omit(unique(annot$Multi_samples[duplicated(annot$Multi_samples)]))
    annot$Multi_samples[!annot$Multi_samples %in% dup]='Unique'
    # Feature
    annot$Features[which(is.na(annot$Features))]='Nothing'
    annot$Features[annot$Features %in% c("locally adjacent")]='Adjacent'
    
    annot$Features[annot$Features=="5-ALA(++)"]='5-ALA(+)'
    annot$Features[annot$Features=='main']='Main'
    annot$Features[annot$Features=='margin']='Margin'
    # Color setting
    mut_col=c('WT'='grey80','MUT'='black','HOT'='tomato')
    relapse_col=c('I'='grey80','R'='yellow')
    patient_col=c("2"='orange','6'='red','4'='darkgreen','1'='blue',
                  '3'='cyan','5'='brown','8'='skyblue','7'='purple',
                  'Unique'='grey80')
    feature_col=c('5-ALA(+)'='red','5-ALA(-)'='blue','Nothing'='grey80',
                  'Main'='darkgreen','Margin'='green','Adjacent'='Orange')
    grade_col=c('II'='darkgreen','III'='purple','IV'='grey80')
    rna.Subtype_col=c("MESENCHYMAL"='firebrick','CLASSICAL'='blue','NEURAL'='darkgreen','PRONEURAL'='purple')
    codel.col=c('1p/19 codel'='black','WT'='grey80')
    EGFRviii.col=c('EGFRvIII'='black','WT'='grey80')
    colnames(annot)=multi_sub(colnames(annot),pattern = c('RNA_subtype4'),
                              replacement = c('RNA_subtype'),exact=T)
    colnames(annot)=gsub(colnames(annot),pattern = "_mut",replacement = "")
    pheatmap(x2$jc.coef[1:2,],clustering_distance_cols = x2$dist.mat,clustering_distance_rows = x$dist.mat,
             annotation_col = annot,show_rownames = F,border_color = 'black',
             annotation_colors = list(NF1=mut_col,PTEN=mut_col,EGFR=mut_col,
                                      TP53=mut_col,RB1=mut_col,PIK3R1=mut_col,
                                      PIK3CA=mut_col,IDH1=mut_col,'1p/19q codel'=codel.col,
                                      'Primary/Relapse'=relapse_col,RNA_subtype=rna.Subtype_col,
                                      Multi_samples=patient_col,Features=feature_col,
                                      Grade=grade_col,'EGFRvIII'=EGFRviii.col),cluster_rows = F,
             title='Figure 1a')
}
#===========================
# Figure 1c panel
#===========================
#x=smc1.snv;y=smc1.snv.by.mass
Figure1c20200409.1354=function(x=smc1.snv,y=smc1.snv.by.mass){
  # Get Somatic 5% VAFs
    x1=vaf_filter20180918.1050(x=x,vaf=0.05)
  # Make coordinates k=x1[[1]]
    x2=ldply(x1,function(k){
      k1=paste0(k[,1],'#',k[,2])
      k2=paste(collapse = ';',k1)
      return(k2)
    })
    colnames(x2)=colnames(y)=c('sample','variant')
  # Calulcate jaccard coefficients
    x2=x2[order(x2$sample),]
    y=y[order(y$sample),]
  # Select samples having variants identified by mass
    x3=x2[which(x2$sample %in% y$sample),]
  # Jaccard matrix
    jc.mat=matrix(0,nrow=nrow(x3),ncol=nrow(x3),
                  dimnames = list(x3$sample,x3$sample))
    for(i in 1:nrow(jc.mat)){ #i=j=1
      i1=unlist(strsplit(x3[i,2],';'))
      for(j in 1:ncol(jc.mat)){
        j1=unlist(strsplit(y[j,2],';'))
        intersect.v=length(intersect(i1,j1))
        union.v=length(union(i1,j1))
        jc.mat[i,j]=intersect.v/union.v
      }
    }
    library(pheatmap)
    pheatmap(jc.mat,cluster_rows = F,cluster_cols = F,
             color = colorRampPalette(c('white','red'))(50),
             main = 'Figure 1c')
}
#===========================
# Figure 1d panel
#===========================
#x=smc1.prot;y=smc1.rna
Figure1d.20200226.0016=function(x=smc1.prot,y=smc1.rna){
  library(plyr)
  # Get symbol expression for protein
    x1=table_expand(x,key_colnames = 'Symbol',sep=';')
    x1$Symbol=gsub(x1$Symbol,pattern = "'",replacement = '')
    x2=dup.matrix(x1,key_column = 'Symbol',max)
  # Get symbol expression for RNA
    y=y$ENSG_id
    colnames(y)=multi_sub(colnames(y),pattern = 'symbol',replacement = 'Symbol',exact=T)
    y1=dup.matrix(y,key_column = 'Symbol',max)
    y1$Symbol=gsub(y1$Symbol,pattern = "'",replacement = '')
  # Select intersect genes
    int.gene=intersect(x1$Symbol,y1$Symbol)
    x3=x2[which(x2$Symbol %in% int.gene),]
    y2=y1[which(y1$Symbol %in% int.gene),]
  # Sort the tables
    x3=x3[order(x3$Symbol),]
    y3=y2[order(y2$Symbol),]    
  # Make column names identical
    colnames(x3)=gsub(colnames(x3),pattern = ' ',replacement = '')
    colnames(x3)=gsub(colnames(x3),pattern = "'",replacement = '')
    int.sample=setdiff(intersect(colnames(x3),colnames(y3)),'Symbol')
    genes=x3$Symbol
    x4=x3[,int.sample]
    y4=y3[,int.sample]
  # Coduct correlation test #k=1
    rho=ldply(1:nrow(x4),.progress='text',function(k){ 
      k1=as.numeric(x4[k,])
      j1=as.numeric(y4[k,])
      k2=cor.test(k1,j1,method = 'spearman')
      return(c('rho'=k2$estimate,'p'=k2$p.value))
    })
    colnames(rho)=c('rho','p')
    x5=data.frame('symbol'=genes,'rho'=rho$rho,'p'=rho$p)
    x5=na.omit(x5)
    x5$fdr=p.adjust(x5$p,'fdr')
  # Get FDR < 5% positive cor
    strong.cor=min(x5$rho[which(x5$rho>0 & x5$fdr<0.05)])
  # Drawing
    tmp=hist(x5$rho,breaks = 100,xlab='Spearman correlation coefficient',main='Figure 1d',freq = F,plot = F)
    # Set colors
    colors=ifelse(tmp$breaks<0,'darkgreen','orange')
    colors=ifelse(tmp$breaks>strong.cor,'red',colors)
    hist(x5$rho,breaks = 100,xlab='Spearman correlation coefficient',main='Figure 1d',freq = F,plot = T,
         col = colors,border = colors)
    lines(density(x5$rho),lwd=2)
}

#===========================
# Figure 1e panel
#===========================
#x=smc1.prot;y=smc1.info
Figure1e.20200226.0101=function(x=smc1.prot,y=smc1.info){
  # Get tumors
    y1=y[which(y$sample.type=='tumor'),]
  # Calculate protein distance
    colnames(x)=gsub(colnames(x),pattern = "'",replacement = '')
    x1=x[,intersect(colnames(x),y1$Sample_name)]
    x2=1-cor(x1)
    x3=hclust(as.dist(x2),method = 'complete')
    plot(x3,hang=-1,main='Figure 3e dendrogram')
    print('Dendrogram plotted')
  # Make table for annotation
    sample.order=x3$labels[x3$order]
    rownames(y1)=y1$Sample_name
  # Annotate hotspot mutations
    genes=c('IDH1','EGFR','TP53','PTEN','PIK3CA')
    for(i in genes){ #i='IDH1'
      wh=paste0(i,'_mut')
      wh2=paste0(wh,'_hotspot')
      wh3=which(!is.na(y1[,wh2]))
      y1[wh3,wh]='Hotspot'
    }
    y1$Multi_sample=unlist(llply(y1$Sample_name,function(k){
      unlist(strsplit(k,'-'))[1]
    }))
    colnames(y1)=gsub(colnames(y1),pattern = '...40',replacement = '')
    dup=y1$Multi_sample[which(duplicated(y1$Multi_sample))]
    y1$Multi_sample[which(!y1$Multi_sample %in% dup)]=NA
    y2=y1[sample.order,rev(c('NF1_mut','IDH1_mut','EGFR_mut','EGFRvIII',
                         'TP53_mut','RB1_mut','PTEN_mut','PIK3CA_mut',
                         'PIK3R1_mut',
                         '1p/19q codel','RNA_subtype4','Grade',
                         "Initial_or_Recurred_sample","Multi_sample",
                         "Features"))]
    for(i in c('NF1_mut','IDH1_mut','EGFR_mut','EGFRvIII',
        'TP53_mut','RB1_mut','PTEN_mut','PIK3CA_mut',
        'PIK3R1_mut')){
          wh=which(!is.na(y2[,i]) & y2[,i]!='Hotspot')
          y2[wh,i]='Mutation'
          wh2=which(is.na(y2[,i]))
          y2[wh2,i]='Wildtype'
    }
    y2$Initial_or_Recurred_sample=gsub(y2$Initial_or_Recurred_sample,pattern = 'Secondary',replacement = 'R')
    unique(y2$RNA_subtype4)
    y2$Features=gsub(y2$Features,pattern = "5-ALA(++)",replacement = "5-ALA(+)",fixed = T)
  # Select necessary columns
    mut.col=c('Mutation'='black','Hotspot'='tomato','Wildtype'='grey')
    p1q19.codel=c('1p/19 codel'='black')
    grade.col=c('IV'='grey','III'='purple','II'='darkgreen')
    rnasubtype.col=c('CLASSICAL'='blue','MESENCHYMAL'='red','NEURAL'='darkgreen','PRONEURAL'='purple')
    primary.relapse.col=c('I'='grey','R'='yellow')
    feature.col=c('5-ALA(+)'='red','5-ALA(-)'='blue','main'='darkgreen','margin'='green','locally adjacent'='orange')
    unique(y2$Features)
    library(pheatmap)
    pheatmap(x2[1:2,sample.order],cluster_rows = F,cluster_cols = F,annotation = y2,
             annotation_colors = list('NF1_mut'=mut.col,
                                        'IDH1_mut'=mut.col,
                                        'EGFR_mut'=mut.col,
                                        'EGFRvIII'=mut.col,
                                        'TP53_mut'=mut.col,
                                        'RB1_mut'=mut.col,
                                        'PTEN_mut'=mut.col,
                                        'PIK3CA_mut'=mut.col,
                                        'PIK3R1_mut'=mut.col,
                                        '1p/19q codel'=p1q19.codel,
                                        'RNA_subtype4'=rnasubtype.col,
                                        'Grade'=grade.col,
                                        "Initial_or_Recurred_sample"=primary.relapse.col,
                                        'Features'=feature.col
                                        ))
}

#===========================
# Permutation test for clinical features
#===========================
#x=smc1.prot;y=smc1.info
perm.test.4.features20200226.2044=function(x=smc1.prot,y=smc1.info){
  check.n.install.lib('wordspace')
  library(wordspace)
  y1=y[which(y$sample.type=='tumor'),]
  y1$Sample_name2=paste0("'",y1$Sample_name)
  x1=x[,y1$Sample_name2]
  # RNA subtype
    rna.subtype=as.list(unique(y1$RNA_subtype4))
    names(rna.subtype)=unique(y1$RNA_subtype4)
    rna.subtype2=llply(rna.subtype,function(k){
      y1$Sample_name2[which(y1$RNA_subtype4==k)]
    })
    rna.subtype.p=cluster_index_P(mm = x1,class.list = rna.subtype2)
  # 5-ALA positivity
    ala.positivity=list('5-ALA(+)','5-ALA(-)')
    names(ala.positivity)=unlist(ala.positivity)
    y1$Features=gsub(y1$Features,pattern = '5-ALA(++)',replacement = '5-ALA(+)',fixed = T)
    ala.positivity2=llply(ala.positivity,function(k){
      y1$Sample_name2[which(y1$Features==k)]
    })
    ala.pos.p=cluster_index_P(mm = x1,class.list = ala.positivity2)
  # Primary relapse
    primary.relapse=list('R','I')
    names(primary.relapse)=unlist(primary.relapse)
    y1$Initial_or_Recurred_sample=gsub(y1$Initial_or_Recurred_sample,
                                       pattern = 'Secondary',replacement = 'R',fixed = T)
    primary.relapse2=llply(primary.relapse,function(k){
      y1$Sample_name2[which(y1$Initial_or_Recurred_sample==k)]
    })
    primary.p=cluster_index_P(mm = x1,class.list = primary.relapse2)
  # Return the result
    p.result=list('RNA.subtype'=rna.subtype.p,'5-ALA'=ala.pos.p,'Primary/Relapse'=primary.p)
    return(p.result)
}

#===========================
# Extended figure 1a
#===========================
# x=smc1.cnv
prepare_Extended_Figure1a20200412.1717=function(x=smc1.cnv){
  message('This is running CBS algorithm not plotting process')
  x$maploc=paste0(x$Chr,'@',as.integer(rowMeans(as.num.mat(x[,2:3]))))
  # Remove duplicated regions
    x1=x[which(!duplicated(x$maploc)),]
    x1$maploc=unlist(llply(x1$maploc,function(k){
      unlist(strsplit(k,'@'))[2]
    }))
    x1$maploc=as.numeric(x1$maploc)
  # Select sample CNV values
    remove.col=c("Chr","Start","End","Symbol","ENSG_id","Entrez_id","maploc")
    x2=x1[,setdiff(colnames(x1),remove.col)]
    check.n.install.lib('DNAcopy',lib.type = 'bioc')
    library(DNAcopy)
  # Run CBS algorithm
    x3=CNA(genomdat = x2,chrom = x1$Chr,maploc = x1$maploc,sampleid = colnames(x2))
    x4=smooth.CNA(x3)
    x5=segment(x4)
    dnacopy.out=x5$output
  # Make maker position
    map=data.frame('Marker_name'=paste0('A',1:nrow(x1)),
                   'Chromosome'=x1$Chr,
                   'Marker_position'=x1$maploc,stringsAsFactors=F)
  # Export tables
    check.n.install.lib('data.table')
    library(data.table)
    message('Expected_Extended_Fig1a_CBS_output.txt was made.')
    fwrite(dnacopy.out,file='Expected_Extended_Fig1a_CBS_output.txt',sep='\t')
    message('Expected_Extended_Fig1a_Marker_position_output.txt was made.')
    fwrite(map,file='Expected_Extended_Fig1a_Marker_position_output.txt',sep='\t')
    message('Go genepattern site and run GISTIC2 with default options')
    print('Parameters are following')
    print('refgene file : Human_Hg19.mat')
    print('seg file : Expected_Extended_Fig1a_CBS_output.txt')
    print('markers file : Expected_Extended_Fig1a_Marker_position_output.txt')
    print('After GISTIC2 running, download amp_qplot.pdf and del_qplot.pdf from Genepattern site.')
}
#===========================
# Extended figure 1b
#===========================
#x=smc1.rna;y=msig
Extended_Figure1b20200408.2057=function(x=smc1.rna,y=msig){
  check.n.install.lib('matrixStats')
  library(matrixStats)
  x=x$ENSG_id
  # Set Genesets #head(x)
  y1=y$updated.GeneSet
  y2=y1[multi_grep(names(y1),pattern = 'VERHAAK_GLIOBLASTOMA')]
  
  # Take max values for redundant symbols
  colnames(x)=multi_sub(colnames(x),pattern = 'symbol',replacement = 'Symbol')
  x1=dup.matrix(x,key_column = 'Symbol',my.Function = max)
  rownames(x1)=gsub(x1$Symbol,pattern = "'",replacement = '')
  
  # change expression values into numeric values
  x1[,-c(1:4)]=as.num.mat(x1[,-c(1:4)])
  x2=t(scale(t(log2(x1[,-c(1:4)]+1)),center=T,scale=T))
  colnames(x2)=gsub(colnames(x2),pattern = "'",replacement = '')
  # Remove all zero genes to increase speed
  zero.genes=rowMeans(x2)
  x3=x2[which(zero.genes!=0),]
  x4=ssgsea(mm = as.num.mat(x3),gene_sets = y2)
  # Sort table
  rownames(x4)=gsub(rownames(x4),pattern = "VERHAAK_GLIOBLASTOMA_",replacement = '')
  rownames(x4)=multi_sub(rownames(x4),pattern = c("PRONEURAL","NEURAL","CLASSICAL","MESENCHYMAL"),
                         replacement = c('PNE','NEU','CLA','MES'),exact = T)
  x5=as.data.frame(t(x4),stringsAsFactors=F)
  # Get class
  subtypes=unlist(llply(1:nrow(x5),function(j){ #j=1
      names(which.max(x5[j,]))
    }))
  x5$subtype=subtypes
  x5$sample=rownames(x5)
  x6=ldply(c('CLA','MES','NEU','PNE'),function(k){ #k='CLA'
    x5[which(x5$subtype==k),]
  })
  rownames(x6)=x6$sample
  pheatmap(x6[,rev(c('CLA','MES','NEU','PNE'))],
           cluster_rows = F,cluster_cols = F,scale = 'row',
           color = colorRampPalette(c('blue','white','red'))(51),
           breaks = seq(-1.5,1.5,length.out = 51),main='Extended Figure 1b')
}

#===========================
# Figure 2A
#===========================
#x=smc1.prot;y=smc1.info
Figure2a20200228.0104=function(x=smc1.prot,y=smc1.info){
  library(ConsensusClusterPlus)
  # Select IDH-WT GBMs only
    y1=y[which(y$Grade=='IV' & is.na(y$IDH1_mut)),]
  # Select IDH-WT GBM expression only
    y1$Sample_name2=paste0("'",y1$Sample_name)
    x1=x[,y1$Sample_name2]
  # Conduct consensus clustering
    x2=ConsensusClusterPlus(d=as.matrix(x1),maxK = 5,plot = F,reps = 1000,pItem = 0.7)
  
  # Draw heatmap
    x3=x2[[2]]
    x4=x3$consensusMatrix
    dimnames(x4)=list(x3$consensusClass,names(x3$consensusClass))
    # Annotation heatmap
    annot=data.frame('GPC'=paste0('GPC',y1[,c("GPC_subtype")]))
    rownames(annot)=y1$Sample_name2
    gpc.col=c('GPC1'='skyblue','GPC2'='tomato')
    
  # Pac score barplot
    pac_score(ccl.stat = x2,kvec = 2:5,threshold=c(0.1,0.9),barplot=T)
    library(grid)
    p<-pheatmap(x4,color = colorRampPalette(c('white','blue'))(50),annotation_col = annot,
                annotation_colors = list('GPC'=gpc.col),silent = T,main='Figure2a')
    return(p)
}
#===========================
# Figure 2B
#===========================
#x=smc1.prot;y=smc1.info
Figure2b20200228.0159=function(x=smc1.prot,y=smc1.info){
  # Select IDH-WT GBMs only
  y1=y[which(y$Grade=='IV' & is.na(y$IDH1_mut)|y$sample.type=='normal'),]
  # Select IDH-WT GBM expression only
  y1$Sample_name2=paste0("'",y1$Sample_name)
  y1=y1[order(y1$GPC_subtype),]
  x1=x[,y1$Sample_name2]
  # Make annotation table
  annot=y1[,c('GPC_subtype','RNA_subtype4')]
  annot$GPC_subtype[which(!is.na(annot$GPC_subtype))]=paste0('GPC',annot$GPC_subtype[which(!is.na(annot$GPC_subtype))])
  rownames(annot)=y1$Sample_name2
  rna.subtype.col=c('MESENCHYMAL'='red','CLASSICAL'='blue','NEURAL'='darkgreen','PRONEURAL'='purple')
  gpc.subtype.col=c('GPC1'='skyblue','GPC2'='tomato')
  pheatmap(x1,color = colorRampPalette(c('blue','white','red'))(51),scale = 'row',
           cluster_rows = T,cluster_cols = F,show_rownames = F,annotation_col = annot,
           annotation_colors = list('RNA_subtype4'=rna.subtype.col,'GPC_subtype'=gpc.subtype.col),
           main = 'Figure 2b')
}
#===========================
# Figure 2c
#===========================
Figure2c20200228.0715=function(x=smc1.info){
  check.n.install.lib('riverplot')
  library(riverplot)
  x=x[which(is.na(x$IDH1_mut) & x$Grade=='IV'),]
  # Making subtype column
  x1=as.data.frame(cbind(x$Sample_name,x$RNA_subtype4,x$GPC_subtype),stringsAsFactors = F)
  colnames(x1)=c('name','rna','pro')
  x2=x1[which(!is.na(x1$rna) & !is.na(x1$pro)),]
  x2$pro=paste0('GPC',x2$pro)
  # Count RNAsubtype and Protein subtype
  rna.count=table(x2$rna)
  pro.count=table(x2$pro)
  
  for(i in names(rna.count)){
    x2$rna[which(x2$rna==i)]=paste0(i,' (',rna.count[i],')')
  }
  for(i in names(pro.count)){
    x2$pro[which(x2$pro==i)]=paste0(i,' (',pro.count[i],')')
  }
  
  # River plot for tracking the subtype transcriptome vs GBC subtype
  # Draw the plot
  node=data.frame(ID=c(unique(x2$rna),sort(unique(x2$pro))),y=c(1:4,c(1.3,3.85)))
  node$x=c(rep(1,4),rep(2,2));rownames(node)=node$ID
  pair=paste0(x2$rna,";",x2$pro)
  
  edge=data.frame(N1=unique(pair),N2=NA,Value=NA,stringsAsFactors = F)
  for(i in 1:nrow(edge)){
    tmp=unlist(strsplit(edge$N1[i],";"))
    edge[i,1:2]=tmp
    edge$Value[i]=length(which(pair %in% paste(unlist(edge[i,1:2]),collapse=";")))
  }
  # Atach numbers
  palette = adjustcolor(c('darkgreen','purple',"red",'blue','skyblue','tomato'),alpha.f = 1)
  styles = lapply(1:length(node$y), function(n) {list(col = palette[n], lty = 0, textcol = "black",srt=0)})
  names(styles) = as.character(node$ID)
  rp=list(nodes=node,edges=edge,styles=styles)
  class(rp)=c(class(rp),'riverplot')
  plot(rp,plot_area=1,yscale=0.04,xscale=1)
  title('Figure 2c')
}

#===========================
# Figure 2d
#===========================
#x=smc1.info
Figure2d20200228.0741=function(x=smc1.info){
  check.n.install.lib('digest')
  library(digest)
  x=x[order(x$GPC_subtype),]
  x1=x[which(is.na(x$IDH1_mut) & x$Grade=='IV'),]
  # Take highest impact mutation
  # I decided to stop > frame > missense > splice donor
  imp_ords=c('Hotspot','stop_gained','frameshift','missense','inframe','splice_donor')
  genes=c('TP53','EGFR','PIK3CA','NF1','RB1','PTEN')
  # Annotate hotspot mutation
  for(i in 1:length(genes)){ #i=1
    wh=setdiff(grep(colnames(x1),pattern = genes[i]),which(colnames(x1)=='EGFRvIII'))
    if(length(wh)>1){
      # Check Hotspot
      for(j in 1:nrow(x1)){
        tmp=unlist(strsplit(x1[j,wh[1]],';'))
        tmp2=unlist(strsplit(x1[j,wh[2]],';'))
        if(any(tmp %in% tmp2) & length(tmp2)>0){
          tmp=multi_sub(tmp,pattern = tmp2,replacement = paste0(tmp2,'(Hotspot)'),exact = T)
          x1[j,wh[1]]=paste(tmp,collapse = ';')
        }
      }
    }
  }
  # Select mutation according to priority
  genes2=paste0(genes,'_mut')
  for(i in 1:length(genes)){
    wh=grep(colnames(x1),pattern = genes2[i])
    for(j in 1:nrow(x1)){
      tmp=unlist(strsplit(x1[j,wh[1]],';'))
      if(length(tmp)>0){
        for(k in 1:length(imp_ords)){
          if(any(grepl(tmp,pattern = imp_ords[k]))){
            tmp=tmp[which(grepl(tmp,pattern = imp_ords[k]))]
          }
        }
        x1[j,wh[1]]=tmp[1]
      }
    }
  }
  # Select necessary columns
  x2=x1[,multi_grep(colnames(x1),pattern = genes)]
  x2=x2[,-multi_grep(colnames(x2),pattern = 'hotspot')]
  rownames(x2)=x1$Sample_name
  # Change names
  for(i in 1:ncol(x2)){
    x2[,i]=multi_sub(x2[,i],pattern = imp_ords,replacement = imp_ords)
  }
  cols=c('tomato','purple','chartreuse4','cornflowerblue','gold')
  x2$EGFRvIII=gsub(x2$EGFRvIII,pattern = 'EGFRvIII',replacement = 'inframe')
  names(cols)=imp_ords[-length(imp_ords)]
  cols=c(cols,'WT'='grey')
  # Calculate p-value
  gpc_list=list(GPC1=x1$Sample_name[which(x1$GPC_subtype==1)],
                GPC2=x1$Sample_name[which(x1$GPC_subtype==2)])
  x3=x2
  x3[is.na(x3)|x3=='NA']='WT'
  
  # Split heatmap into two
  library(ComplexHeatmap)
  x4=as.matrix(x3[which(rownames(x3) %in% gpc_list$GPC1),])
  x5=as.matrix(x3[which(rownames(x3) %in% gpc_list$GPC2),])
  x4=ComplexHeatmap::Heatmap(t(x4),col=cols,rect_gp = gpar(col = "white", lwd = 2),column_title = 'GPC1')
  x5=ComplexHeatmap::Heatmap(t(x5),col=cols,rect_gp = gpar(col = "white", lwd = 2),column_title = 'GPC2')
  draw(x4+x5)
}

#===========================
# Figure 2d lolliplot
#===========================
#x=smc1.info;y='EGFR';y='PIK3CA'
Figure2d.lolliplot20200228.0830=function(x=smc1.info,y='EGFR'){
  check.n.install.lib('biomaRt',lib.type = 'bioc')
  library(biomaRt)
    # Get domain info for genes
      domain.info=get_domain_info(symbol = y)
    # Get GPC1 vs GPC2 mutation info
      x1=x[,multi_grep(colnames(x),pattern = c('GPC_subtype',y))]
    # Make mutation table for lolliplot
      col.names=colnames(x1)[grep(colnames(x1),pattern = y)]
      x2=tidyr::unite(x1,col='variant',col.names,sep=';',remove=F)
      wh=grep(colnames(x2),pattern = 'hot')
      for(i in 1:nrow(x2)){
        tmp=unlist(strsplit(x2[i,'variant'],';'))
        tmp2=unlist(strsplit(x2[i,wh],';'))
        if(length(tmp2)>0){
          tmp=unique(multi_sub(tmp,pattern = tmp2,replacement = paste0(tmp2,'(hot)'),exact=T))
        }
        tmp=unlist(strsplit(tmp,';'))
        tmp=tmp[which(tmp!='')]
        x2[i,'variant']=paste(tmp,collapse=';')
      }
      x2$variant=gsub(x2$variant,pattern = 'EGFRvIII',replacement = 'EGFRvIII(267)')
    # Make mutation table for lolliplot-2
      x3=table_expand(x2,key_colnames = 'variant',sep=';')
    # Remove 0 mutation samples
      x3=x3[which(x3$variant!='NA'),]
      x4=as.data.frame(table(x3$variant,paste0('GPC',x3$GPC_subtype)),stringsAsFactors=F)
      colnames(x4)=c('variant','subtype','number')
      x4=x4[which(x4$number!=0),]
    # make other columns
    # Symbol
    x4$symbol=y
    # Position
    x4$pos=gsub(x4$variant,pattern = '[a-z]',replacement = '')
    x4$pos=gsub(x4$pos,pattern = '[[:punct:]]',replacement = '')
    x4$pos=as.numeric(gsub(x4$pos,pattern = '[A-Z]',replacement = ''))
    # conseqeunce
    tmp=unlist(strsplit(x4$variant,')',fixed = T))
    x4$consq=gsub(x4$variant,pattern = '(hot)',replacement = '',fixed = T)
    x4$consq=multi_sub(x4$consq,pattern = c('missense(','frameshift'),
                       replacement = c('','frs'),partial = T,)
    # hotspot
    x4$hotspot='black'
    x4$hotspot[grep(x4$variant,pattern = 'hot')]='yellow'
    # Select columns
    x4=x4[,c('symbol','pos','consq','subtype','number','hotspot')]
    # Subtype
    x4$subtype=ifelse(x4$subtype=='GPC1','skyblue','tomato')
    # Make Grange object
    library(trackViewer)
    snp_gr=GRanges(seqnames = x4$symbol,
                   IRanges(x4$pos,width = 1,
                           names=x4$consq),
                   color=x4$subtype,
                   border=x4$hotspot)
    snp_gr$label.parameter.rot=45
    snp_gr$label.parameter.gp=gpar(col='black')
    snp_gr$label=NA;snp_gr$label[which(snp_gr$border=='yellow')]='H'
    snp_gr$label.col='white'
    snp_gr$score=as.numeric(x4$number)
    # domain Granges
    if(nrow(domain.info)>0){
      domain.info$color=as.factor(domain.info$domain.name)
      dom_gr=GRanges(seqnames = rep(y,nrow(domain.info)+2),
                     IRanges(start=c(1,domain.info$pfam_start,
                                     as.numeric(unique(domain.info$length))),
                             width=c(1,domain.info$pfam_end-domain.info$pfam_start,1),
                             names=c("",domain.info$domain.name,"")),
                     fill=c('black',domain.info$color,'black'),
                     height=rep(0.05,nrow(domain.info)+2)
      )
    }else{ # In case of pfam id doesn't exists
      dom_gr=GRanges(seqnames = rep(syms[i],2),
                     IRanges(start = c(1,as.numeric(unique(dmat$length[wh]))),
                             width = c(1,1),
                             names=c("","")),
                     fill=c('black','black'),
                     height=rep(0.05,2))
    }
    trackViewer::lolliplot(snp_gr,dom_gr,ylab=y)
    }

#===========================
# Extended Figure2A
#===========================
#x=smc1.prot;y=smc1.info
Extended_Figure2a20200406.1932=function(x=smc1.prot,y=smc1.info){
  check.n.install.lib('ConsensusClusterPlus')
    library(ConsensusClusterPlus)
    # Select IDH-WT GBMs only
    y1=y[which(y$sample.type=='tumor'),]
    # Select IDH-WT GBM expression only
    y1$Sample_name2=paste0("'",y1$Sample_name)
    x1=x[,y1$Sample_name2]
    # Conduct consensus clustering
    x2=ConsensusClusterPlus(d=as.matrix(x1),maxK = 5,plot = F,reps = 1000,pItem = 0.7)
    
    # Draw heatmap
    x3=x2[[2]]
    x4=x3$consensusMatrix
    dimnames(x4)=list(x3$consensusClass,names(x3$consensusClass))
    # Annotation heatmap
    annot=data.frame('GPC'=paste0('A',y1[,c("GPC_subtype")]))
    rownames(annot)=y1$Sample_name2
    gpc.col=c('A1'='blue','A2'='magenta')
    p<-pheatmap(x4,color = colorRampPalette(c('white','blue'))(50),annotation_col = annot,
             annotation_colors = list('GPC'=gpc.col),silent = T)
    # Pac score barplot
    pac_score(ccl.stat = x2,kvec = 2:5,threshold=c(0.1,0.9),barplot=T)
    #return results
    return(list(cls=x3$consensusClass,heatmap=p))
}

#===========================
# Extended Figure2B
#===========================
#x=smc1.prot;y=smc1.info
Extended_Figure2b20200406.1935=function(x=smc1.prot,y=smc1.info){
  check.n.install.lib('ConsensusClusterPlus')
  library(ConsensusClusterPlus)
  # Select non redundant samples
  y1=y[which(y$sample.type=='tumor'),]
  y1$patient=unlist(llply(y1$Sample_name,function(k){
    unlist(strsplit(k,'-'))[1]
  }))
  dup.sample=y1$patient[which(duplicated(y1$patient))]
  y2=y1[which(!y1$patient %in% dup.sample),]
  # Select IDH-WT GBM expression only
  y2$Sample_name2=paste0("'",y2$Sample_name)
  x1=x[,y2$Sample_name2]
  # Conduct consensus clustering
  x2=ConsensusClusterPlus(d=as.matrix(x1),
                          maxK = 5,plot = F,reps = 1000,pItem = 0.7)
  
  # Draw heatmap
  x3=x2[[2]]
  x4=x3$consensusMatrix
  dimnames(x4)=list(x3$consensusClass,names(x3$consensusClass))
  # Annotation heatmap
  annot=data.frame('GPC'=paste0('B',x3$consensusClass))
  rownames(annot)=y2$Sample_name2
  gpc.col=c('B1'='dodgerblue','B2'='brown1')
  p=pheatmap(x4,color = colorRampPalette(c('white','blue'))(50),annotation_col = annot,
           annotation_colors = list('GPC'=gpc.col),silent = T)
  # Pac score barplot
  pac_score(ccl.stat = x2,kvec = 2:5,threshold=c(0.1,0.9),barplot=T)
  return(list(cls=x3$consensusClass,heatmap=p))
}

#===========================
# Extended Figure2c
#===========================
#x=smc1.prot;y=smc1.info
Extended_Figure2c20200406.1947=function(x=smc1.prot,y=smc1.info){
  check.n.install.lib('ConsensusClusterPlus')
  library(ConsensusClusterPlus)
  # Select GBM samples
  y1=y[which(y$sample.type=='tumor' & y$Grade=='IV'),]
  y2=y1
  # Select IDH-WT GBM expression only
  y2$Sample_name2=paste0("'",y2$Sample_name)
  x1=x[,y2$Sample_name2]
  # Conduct consensus clustering
  x2=ConsensusClusterPlus(d=as.matrix(x1),
                          maxK = 5,plot = F,reps = 1000,pItem = 0.6)
  
  # Draw heatmap
  x3=x2[[2]]
  x4=x3$consensusMatrix
  dimnames(x4)=list(x3$consensusClass,names(x3$consensusClass))
  # Annotation heatmap
  annot=data.frame('GPC'=paste0('C',x3$consensusClass))
  rownames(annot)=y2$Sample_name2
  gpc.col=c('C1'='deepskyblue4','C2'='darkred')
  p=pheatmap(x4,color = colorRampPalette(c('white','blue'))(50),annotation_col = annot,
           annotation_colors = list('GPC'=gpc.col),silent = T)
  # Pac score barplot
  pac_score(ccl.stat = x2,kvec = 2:5,threshold=c(0.1,0.9),barplot=T)
  return(list(cls=x3$consensusClass,heatmap=p))
}

#===========================
# Extended Figure2d
#===========================
#x=smc1.info;a=ex2a;b=ex2b;d=ex2c
Extended_Figure2d20200406.1950=function(x=smc1.info,
                                        a=ex2a,b=ex2b,d=ex2c){
  check.n.install.lib('dplyr')
  library(dplyr)
  # Merge matrix
    a1=data.frame(sample=names(a$cls),GPC=a$cls)
    b1=data.frame(sample=names(b$cls),GPC=b$cls)
    d1=data.frame(sample=names(d$cls),GPC=d$cls)
    x1=x[which(x$Grade=='IV' & is.na(x$IDH1_mut)),]
    x2=data.frame(sample=paste0("'",x1$Sample_name),GPC=x1$GPC_subtype)
    x3=mat.merge(a1,b1,by='sample',all = T)
    x4=mat.merge(x3,d1,by='sample',all = T)
    x5=mat.merge(x4,x2,by='sample',all = T)
  # Draw heatmap
    x6=t(x5[,-1])
    colnames(x6)=x5[,1]
    rownames(x6)=c('All','Unique','GBM','IDH-WT GBM')
    x6=x6[,order(x6[1,])]
    pheatmap(x6,cluster_rows = F,cluster_cols = F,main = 'Ex Figure4d')
}

#===========================
# Figure 3a PCA plot
#===========================
#x=smc1.prot;y=smc1.info
Figure3a.PCAplot20200309.2120=function(x=smc1.prot,y=smc1.info){
  # Select IDH-WT GBMs
    y1=y[which(y$Grade=='IV' & is.na(y$IDH1_mut)),]
    rownames(x)=paste0('A',1:nrow(x),'@',x$Symbol)
    rownames(x)=gsub(rownames(x),pattern = "'",replacement = '')
    colnames(x)=gsub(colnames(x),pattern = "'",replacement = '')
    x1=x[,y1$Sample_name]
  # Conduct PCA
    x2=prcomp(t(x1),center = T,scale. = T)
  # Get coordinates
    x3=x2$x
  # Set colors
    y1$color=ifelse(y1$GPC_subtype==1,'skyblue','tomato')
    plot(x3[,1],-x3[,2],xlab='PC1',ylab='PC2',col=y1$color,pch=16,cex=1.5,
         main='Figure 3a left')
  # Get enrichment
    pc1.loadings=x2$rotation[,1]
    pc1.abs.loadings=sort(abs(pc1.loadings),decreasing = T)
    top.pc1.loadings=pc1.abs.loadings[1:(length(pc1.abs.loadings)*0.1)]
    backgroud=pc1.abs.loadings
  # Return the result
    return(list(input=top.pc1.loadings,back=backgroud))
}

#===========================
# Figure 3a right barplot 
#===========================
#x=pc1.top.ten;y=corum.v3.0.proteinComplex;z=updated.symbol.msigDB.v6.2_20190225
Figure3a.barPlot20200309.2120=function(x=pc1.top.ten,
                                       y=corum.v3.0.proteinComplex,
                                       z=updated.symbol.msigDB.v6.2_20190225){
  # Select optimal genesets
    y1=y$updated.GeneSet
    names(y1)=paste0('CORUM_',names(y1))
    z1=z$updated.GeneSet
    names(z1)=paste0('mSIG_',names(z1))
    y2=c(y1,z1)
  # Filtering the genesets
    y2.size=unlist(llply(y2,function(k){
      length(k)
    }))
    y3=y2[which(y2.size<151 & y2.size>4)]
  # Get genes
    x1=llply(x,function(k){ #k=x[[1]]
      names(k)
    })
    x2=llply(x1,function(k){ #k=x1[[2]]
      k1=unlist(llply(k,function(j){ #j=k[110]
        unlist(strsplit(j,'@'))[2]
      }))
      return(k1)
    })
    x3=llply(x2,function(k){
      unlist(strsplit(k,';'))
    })
    input=x3$input
    backgroud=setdiff(x3$back,x2$input)
  # Conduct GSA
    x4=Gene.Enrichment.Analysis(input.gene = input,background.gene = backgroud,geneSet.list = y3)
  # Select genesets at least four hits
    x4$In=as.numeric(x4$In)
    x4$Pvalue=as.numeric(x4$Pvalue)
    x5=x4[which(x4$In>=4),]
    x5$FDR=p.adjust(x5$Pvalue,'fdr')
  # Sort table
    x5=x5[order(x5$Pvalue),]
  # Select top 5 genesets
    x6=x5[1:5,]
    x6$logp=-log10(x6$Pvalue)
    barplot(x6$logp,col = 'pink',names.arg = x6$GeneSet,las=2,ylim=c(0,35),ylab='-Log (P)',
            main='Figure 3a right')
  # Return result
    return(x5)
}

#===========================
# Figure 3b heatmap
#===========================
#x=smc1.prot;y=smc1.info;a=smc1.rna;z=msig
Figure3b.heatmap20200309.2235=function(x=smc1.prot,y=smc1.info,
                                       z=msig,
                                       a=smc1.rna){
  a=a$ENSG_id
  # Select IDH-WT GBM and normals
    y1=y[which(is.na(y$IDH1_mut) & y$Grade=='IV' | y$sample.type=='normal'),]
    colnames(x)=gsub(colnames(x),pattern = "'",replacement = '')
    colnames(a)=gsub(colnames(a),pattern = "'",replacement = '')
    colnames(a)=multi_sub(colnames(a),pattern = 'symbol',replacement = 'Symbol',exact = T)
    x1=x[,c('Symbol',y1$Sample_name)]
    a1=a[,c("Symbol",intersect(colnames(a),y1$Sample_name))]
  # Select geneset
    z1=z$updated.GeneSet["KEGG_OXIDATIVE_PHOSPHORYLATION"]
  # Sort table
    y1=y1[order(y1$GPC_subtype),]
    x2=x1[multi_grep(x1$Symbol,pattern = z1$KEGG_OXIDATIVE_PHOSPHORYLATION),]
    a2=a1[multi_grep(a1$Symbol,pattern = z1$KEGG_OXIDATIVE_PHOSPHORYLATION),]
    x3=table_expand(x2,key_colnames = 'Symbol')
    a3=a2
    x3$Symbol=gsub(x3$Symbol,pattern = "'",replacement = '')
    a3$Symbol=gsub(a3$Symbol,pattern = "'",replacement = '')
  # Select genes
    z2=unlist(z1)
    z3=intersect(z2,x3$Symbol)
    z4=intersect(z3,a3$Symbol)
    x4=x3[which(x3$Symbol %in% z4),y1$Sample_name]
    a4=a3[which(a3$Symbol %in% z4),]
    a4=dup.matrix(a4,key_column = 'Symbol',max)
    a5=log2(a4[,setdiff(y1$Sample_name,c("586N","655N","753N","608N"))]+1)
  # Make annotation table
    annot=data.frame('Subtype'=paste0('GPC',y1$GPC_subtype),stringsAsFactors=F)
    annot$Subtype[which(annot$Subtype=='GPCNA')]='Normal_tissue'
    rownames(annot)=y1$Sample_name
    subtype.col=c('GPC1'='skyblue','GPC2'='tomato','Normal_tissue'='grey')
    pheatmap(x4,color = colorRampPalette(c('blue','white','red'))(31),scale = 'row',cluster_cols = F,
             annotation = annot,annotation_colors = list('Subtype'=subtype.col),show_rownames = F,
             main = 'OXPHOS protein (Z-score)')
    
    pheatmap(a5,color = colorRampPalette(c('blue','white','red'))(31),scale = 'row',cluster_cols = F,
             annotation = annot,annotation_colors = list('Subtype'=subtype.col),show_rownames = F,
             main = 'OXPHOS mRNA (Z-score)')
}

#===========================
# Figure 3b ssGSEA boxplot
#===========================
#x=smc1.prot;y=smc1.info;a=smc1.rna
#z=msig
Figure3b.boxplot20200309.2352=function(x=smc1.prot,y=smc1.info,
                                       z=updated.symbol.msigDB.v6.2_20190225,
                                       a=smc1.rna){
  a=a$ENSG_id
  # Select geneset
    z1=z$updated.GeneSet["KEGG_OXIDATIVE_PHOSPHORYLATION"]
    y1=y[which(is.na(y$IDH1_mut) & y$Grade=='IV' | y$sample.type=='normal'),]
    colnames(x)=gsub(colnames(x),pattern = "'",replacement = '')
    colnames(a)=gsub(colnames(a),pattern = "'",replacement = '')
    x1=x[,c('Symbol',y1$Sample_name)]
    colnames(a)=multi_sub(colnames(a),pattern = 'symbol',replacement = 'Symbol',exact = T)
    a1=a[,c("Symbol",intersect(colnames(a),y1$Sample_name))]
  # Formatting matrix
    x1$Symbol=gsub(x1$Symbol,pattern = "'",replacement = '')
    a1$Symbol=gsub(a1$Symbol,pattern = "'",replacement = '')
    x2=table_expand(x1,key_colnames = 'Symbol')
    a2=table_expand(a1,key_colnames = 'Symbol')
    a2[,-1]=log2(a2[,-1]+1)
  # Take max values
    x3=dup.matrix(x2,key_column = 'Symbol',my.Function = max)
    a3=dup.matrix(a2,key_column = 'Symbol',my.Function = max)
    x4=x3[,-1]
    rownames(x4)=x3$Symbol
    a4=t(scale(t(a3[,-1]),center=T,scale=F))
    rownames(a4)=a3$Symbol
  # Conduct ssGSEA
    x5=ssgsea(mm = as.matrix(x4),gene_sets = z1)
    a5=ssgsea(mm = as.matrix(a4),gene_sets = z1)
  # Conduct wilcoxon rank sum test
    wh.GPC1=which(y1$GPC_subtype==1)
    wh.GPC2=which(y1$GPC_subtype==2)
    wh.norm=which(is.na(y1$GPC_subtype))
    gpc.list=list('GPC1'=y1$Sample_name[wh.GPC1],
                  'GPC2'=y1$Sample_name[wh.GPC2],
                  'Normal'=y1$Sample_name[wh.norm])
    par(mfrow=c(1,2))
  # Protein
    x6=llply(gpc.list,function(k){
      as.numeric(x5[,k])
    })
    sig_boxplot(x6,test_method = 'wilcox',colors = c('skyblue','tomato','grey'),
                yaxis_lab = 'OXPHOS protein (ssGSEA score)',
                title='Figure 3b (left)')
  # Protein
    x6.mrna=llply(gpc.list[1:2],function(k){
      as.numeric(a5[,k])
    })
    sig_boxplot(x6.mrna,test_method = 'wilcox',colors = c('skyblue','tomato','grey'),
                yaxis_lab = 'OXPHOS mRNA (ssGSEA score)',ns_visualize = T,
                title='Figure 3b (left)')
    par(mfrow=c(1,1))
}

#===========================
# Extended Figure3b
#===========================
#x=smc1.rna;y=smc1.info
Extended_Figure3b20200404.1544=function(x=smc1.rna,y=smc1.info){
  x=x$ENST_id
  par(mfrow=c(1,2))
  # Message Get PKM isoform informations
    pkm.isoform=getting_ensg.id_and_pkm.sequence20180816.1410(genes=c('PKM'))
    pkm.isoform2=na.omit(pkm.isoform)
  # Get expression values for PKM isoforms
    x1=x[which(x$ENST_id %in% pkm.isoform2$ensembl_transcript_id),]
  # Select IDH-WT GBMs
    y1=y[which(is.na(y$IDH1_mut) & y$Grade=='IV'),]
    y1$Sample_name2=paste0("'",y1$Sample_name)
    rownames(x1)=x1$ENST_id
    x2=x1[,y1$Sample_name2]
  # Get average PKM isoform expression values
    wh=which(y1$GPC_subtype==1)
    gpc.list=list('GPC1'=y1$Sample_name2[wh],'GPC2'=y1$Sample_name2[-wh])
    for(i in c('PKM1','PKM2')){ #i='PKM1'
      wh=which(pkm.isoform2$pkm_isoform==i)
      wh2=pkm.isoform2$ensembl_transcript_id[wh]
      gpc.list2=llply(gpc.list,function(k){
        colMeans(x2[wh2,k])
      })
      ylabel=paste0(i,' (Log2 (FPKM+1))')
      sig_boxplot(gpc.list2,yaxis_lab = ylabel,colors = c('skyblue','tomato'),
                  title='Ex Figure3b',ns_visualize = T)
    }
    par(mfrow=c(1,1))
}

#===========================
# Figure 3e metabolomics data
#===========================
#x=ccle.datas;p_cut=0.05
Figure3e20200407.1627=function(x=ccle.datas,p_cut=0.05){
  # Get sGPC infos
    a1=x$exp
    b1=x$sGPC.info$result
    c1=x$metabolome
    d1=x$sGPC.info$classifiers
  # set sGPC lists
    b1$minP=ifelse(b1$GPC1.p<b1$GPC2.p,b1$GPC1.p,b1$GPC2.p)
    b2=b1[which(b1$minP<p_cut),]
    wh=which(b2$sGPC=='sGPC1')
    sGPC.list.ccle=list('sGPC1'=b2$sample[wh],'sGPC2'=b2$sample[-wh])
  # Draw heatmap
    rownames(a1)=unlist(llply(as.character(a1$Name),function(k){
      unlist(strsplit(k,'.',fixed = T))[1]
    }))
    b1$sGPC=ifelse(b1$minP<p_cut,b1$sGPC,'Unknown')
    a2=a1[unlist(d1),as.character(b1$sample)]
    annot.col=data.frame('Predicted.subtype'=b1$sGPC)
    rownames(annot.col)=b1$sample
    annot.row=data.frame('Classifier.genes'=rep(c('GPC1-high','GPC2-high'),each=100))
    rownames(annot.row)=unlist(d1)
    sGPC.col=c('sGPC1'='skyblue','sGPC2'='tomato','Unknown'='grey')
    gpc.high.col=c('GPC1-high'='dodgerblue4','GPC2-high'='red')
    figure3e.left.heatmap<-
      pheatmap(as.num.mat(a2),cluster_rows = F,cluster_cols = F,scale='row',
             show_colnames = F,show_rownames = F,
             color = colorRampPalette(c('blue','white','red'))(51),
             breaks = seq(-3,3,length.out = 51),
             annotation_col = annot.col,annotation_row = annot.row,
             annotation_colors = list('Predicted.subtype'=sGPC.col,
                                      'Classifier.genes'=gpc.high.col),
             gaps_col = c(13,35),main='Figure 3e left')
    
  # Select Latate level and compare
    sGPC.list.ccle2=llply(sGPC.list.ccle,function(k){
      wh=which(c1$CCLE_ID %in% k)
      return(c1$lactate[wh])
    })
  # Draw sigboxplot
    sig_boxplot(sGPC.list.ccle2,colors = c('skyblue','tomato'),
                yaxis_lab = 'Lactate level',title = 'Figure 3e right')
    axis(side=2,at=c(5.25,5.75,6.25),las=2)
    return(figure3e.left.heatmap)
}

#===========================
# Extended Figure3c
#===========================
#x=smc1.prot;y=smc1.info
Extended_Figure3c20200407.1326=function(x=smc1.prot,y=smc1.info){
  # Select samples IDH-WT GBMs and normals
    y1=y[which((is.na(y$IDH1_mut) & y$Grade=='IV')|y$sample.type=='normal'),]
  # Set gpc subtypes
    colnames(x)=gsub(colnames(x),pattern = "'",replacement = '')
    wh=which(y1$sample.type=='normal')
    wh.gpc1=which(y1$GPC_subtype==1)
    wh.gpc2=which(y1$GPC_subtype==2)
    gpc.list=list('Normal'=y1$Sample_name[wh],
                  'GPC1'=y1$Sample_name[wh.gpc1],
                  'GPC2'=y1$Sample_name[wh.gpc2])
  # Get IDH expression
    x1=x[grep(x$Symbol,pattern = 'IDH1'),]
    gpc.list2=llply(gpc.list,function(k){
      as.numeric(x1[1,k])
    })
  # Draw sig boxplot
    sig_boxplot(gpc.list2,colors = c('white','skyblue','tomato'),
                yaxis_lab = 'IDH1 (Log2)',digits = 2,title='Ex Fig3c')
}

#================================
# Select necessary column and convert it into log2-ratio
#================================
data.parse_and_log2.protein.exp.ratio20200310.0045=function(x=pep.set1,y=pep.set2,z=sam2){
  # Remove reverse peptide sequences
  x=x[x$Reverse!='+',];y=y[y$Reverse!='+',]
  # Remove potential contaminant peptides
  x=x[x$`Potential contaminant`!='+',];y=y[y$`Potential contaminant`!='+',]
  # Select necessary columns (peptide sequence, modifications, proteins, gene name,Oxidation (M) site IDs, intensities)
  select.columns.x=c('Sequence','Modifications',"Proteins","Gene Names","Oxidation (M) site IDs",
                     colnames(x)[grep(colnames(x),pattern = "Reporter intensity corrected")])
  select.columns.y=c('Sequence','Modifications',"Proteins","Gene Names","Oxidation (M) site IDs",
                     colnames(y)[grep(colnames(y),pattern = "Reporter intensity corrected")])
  x=x[,select.columns.x];y=y[,select.columns.y]
  # Reporter inensity corrected -> TMT
  colnames(x)=gsub(colnames(x),pattern = "Reporter intensity corrected ",replacement = "TMT")
  colnames(y)=gsub(colnames(y),pattern = "Reporter intensity corrected ",replacement = "TMT")
  # Remove TMT only column
  exclude.columns=paste0('TMT',0:5)
  x=x[,which(!colnames(x) %in% exclude.columns)];y=y[,which(!colnames(y) %in% exclude.columns)]
  # Exclude peptide having no protein groups and protein ids
  x=x[which(!is.na(x$Proteins)),];y=y[which(!is.na(y$Proteins)),]
  # Merge the two tables
  # Make identifiers (sequence, modifications, oxidation)
  x$identifier=paste0(x$Sequence,'@',x$Modifications,'@',x$`Oxidation (M) site IDs`)
  y$identifier=paste0(y$Sequence,'@',y$Modifications,'@',y$`Oxidation (M) site IDs`)
  # Merge the tables according to their identifiers  
  x1=merge(x,y,by='identifier',all=T)
  # Sort the table
  x1=x1[,c(c(1:6,43:47),setdiff(c(1:ncol(x1)),c(1:6,43:47)))]
  # Remove peptides that were identified but not quantified in GIS.
  parsed.data=x1
  quantified.peptide.in.gis=which(x1$`TMT0 Set1`!=0 & x1$`TMT0 Set2`!=0 & x1$`TMT0 Set3`!=0 &
                                    x1$`TMT0 Set4`!=0 & x1$`TMT0 Set5`!=0 & x1$`TMT0 Set6`!=0 &
                                    x1$`TMT5 set7`!=0 & x1$`TMT0 set8`!=0 & x1$`TMT0 set9`!=0 &
                                    x1$`TMT0 set10`!=0 & x1$`TMT0 set11`!=0)
  x2=x1[quantified.peptide.in.gis,]
  # Calculate log2-ratio
  # Batch1 log2-ratio
  bat1s=expand.grid(TMT=paste0('TMT',1:5),SET=paste0('Set',1:6),GIS='TMT0',stringsAsFactors = F)
  sample1s=c(paste0(bat1s$TMT,' ',bat1s$SET),paste0('TMT',1:4, ' set7'))
  gis1s=c(paste0(bat1s$GIS,' ',bat1s$SET),rep('TMT5 set7',4))
  for(i in 1:length(gis1s)){
    sample.wh=which(colnames(x2)==sample1s[i])
    gis.wh=which(colnames(x2)==gis1s[i])
    x2[,sample.wh]=log2(x2[,sample.wh]/x2[,gis.wh])
  }
  # Batch2 log2-ratio
  bat2s=expand.grid(TMT=paste0('TMT',1:5),SET=paste0('set',8:11),GIS='TMT0',stringsAsFactors = F)
  sample2s=c(paste0(bat2s$TMT,' ',bat2s$SET))
  gis2s=c(paste0(bat2s$GIS,' ',bat2s$SET))
  for(i in 1:length(gis2s)){
    sample.wh=which(colnames(x2)==sample2s[i])
    gis.wh=which(colnames(x2)==gis2s[i])
    x2[,sample.wh]=log2(x2[,sample.wh]/x2[,gis.wh])
  }
  # Select samples exp
  z$TMT=z$TMT-126
  z$tmt_sets=paste0('TMT',z$TMT,' SET',z$Set)
  colnames(x2)=toupper(colnames(x2))
  x2=x2[,c(1:11,which(colnames(x2) %in% z$tmt_sets))]
  # Change column name
  colnames(x2)=multi_sub(colnames(x2),pattern = z$tmt_sets,replacement = z$Sample_name,exact = T)
  # Return the tables
  return(list(parsed.data=parsed.data,log2.table=x2))
}

#===========================
# Get PKM expression and parsed the expression data
#===========================
select_pkm.peptide.exp_and_table.expand_according.to.pID20181010.1701=function(x=pep.exp$log2.table,gene='PKM'){
  # Select peptide for target gene
  wh=unique(c(multi_grep(x$`GENE NAMES.X`,pattern = gene),multi_grep(x$`GENE NAMES.Y`,pattern = gene)))
  x1=x[wh,]
  # Do table expand
  x2=table_expand(key_colnames = "PROTEINS.X",mat = x1,sep = ';')
  x2=table_expand(key_colnames = "GENE NAMES.X",mat = x2,sep = ';')
  x2=x2[which(x2$`GENE NAMES.X` %in% gene),]
  # Return the table
  return(x2)
}




#========================================
# PKM jitter plot
#========================================
#x=pkm.pep2;y=smc1.info;title='IDH-WT GBM';jitter=F
pkm_jitterplot_specific.peptide_for_pkm.isoform20181015.1223=function(x=pkm.pep2,y=smc1.info,title='IDH-WT GBM',
                                                                      jitter=F){
  x=x[-grep(x$isoform,pattern = ';'),]
  # Split table into two tables for PKM1 and PKM2
  pkm1i=x[grep(x$isoform,pattern = 'PKM1'),]
  pkm2i=x[grep(x$isoform,pattern = 'PKM2'),]
  # Remove duplicates
  pkm1i=pkm1i[!duplicated(pkm1i$IDENTIFIER),]
  pkm2i=pkm2i[!duplicated(pkm2i$IDENTIFIER),]
  # Calculate mean
  pkm1i.exp=colMeans(pkm1i[,y$Sample_name],na.rm=T)
  pkm2i.exp=colMeans(pkm2i[,y$Sample_name],na.rm=T)
  # Draw sigboxplot
  gpc.list=list(GPC1=y$Sample_name[which(y$GPC_subtype==1 & is.na(y$IDH1_mut) & y$Grade=='IV')],
                GPC2=y$Sample_name[which(y$GPC_subtype==2 & is.na(y$IDH1_mut) & y$Grade=='IV')])
  # PKM1 isoform
  gpc.list.pkm1.exp=list(GPC1=pkm1i.exp[gpc.list[[1]]],GPC2=pkm1i.exp[gpc.list[[2]]])
  sig_boxplot(gpc.list.pkm1.exp,colors = c('skyblue','tomato'),yaxis_lab = 'Peptide expression(Log2)',
              title = paste0('PKM1 specific\npeptide expression\n(Exon9) ',title),jitter.only = jitter)
  
  gpc.list.pkm2.exp=list(GPC1=pkm2i.exp[gpc.list[[1]]],GPC2=pkm2i.exp[gpc.list[[2]]])
  sig_boxplot(gpc.list.pkm2.exp,colors = c('skyblue','tomato'),yaxis_lab = 'Peptide expression(Log2)',
              title = paste0('PKM2 specific\npeptide expression\n(Exon10) ',title),jitter.only = jitter)
}

#===========================
# Figure 3d PKM boxplot
#===========================
#x=set1.peptide;y=set2.peptide;z=smc1.info
Figure3d.PKM.boxplot20200310.0031=function(x=NULL,y=NULL,z=smc1.info){
  if(is.null(x)|is.null(y)){
    message('Put peptide quantification data from maxquant')
    stop()
  }else{
  pep.exp=data.parse_and_log2.protein.exp.ratio20200310.0045(x=x,y=y,z=smc1.info)
  
  #===========================
  # Get PKM expression and table.expand according to protein ID
  #===========================
  pkm.pep=select_pkm.peptide.exp_and_table.expand_according.to.pID20181010.1701(x=pep.exp$log2.table,gene='PKM')
  #===========================
  # Align the peptide sequence on pkm isoform
  # PKM isoform sequences were derived from Uniprot
  #===========================
  pkm.loc=align_peptide.on.pkm_isoforms20181010.1812(x=unique(pkm.pep$SEQUENCE.X))
  # Merge the location and peptide info
  pkm.pep2=merge(pkm.pep,pkm.loc,by.x='SEQUENCE.X',by.y='peptide')
  
  #===========================
  # Draw boxplot
  #===========================
  z1=z[which(z$Grade=='IV' & is.na(z$IDH1_mut)),]
  par(mfrow=c(1,2))
  pkm_jitterplot_specific.peptide_for_pkm.isoform20181015.1223(x=pkm.pep2,y=z1)
  }
  par(mfrow=c(1,1))
}

#===========================
# Figure 4c SMC1
#===========================
#x=smc1.rna;y=smc1.info
Figure4c.SMC1.20200404.1718=function(x=smc1.rna,y=smc1.info){
  x=x$ENSG_id
  colnames(x)=multi_sub(colnames(x),pattern = 'symbol',replacement = 'Symbol',exact = T)
  x$Symbol=gsub(x$Symbol,pattern = "'",replacement = '')
  colnames(x)=gsub(colnames(x),pattern = "'",replacement = '')
  x1=x[which(x$Symbol %in% c('CD274','PDCD1LG2')),]
  # Get IDH-WT GBMs
    y1=y[which(is.na(y$IDH1_mut) & y$Grade=='IV'),]
    #y1=y[which(y$sample.type=='tumor'),]
    wh=which(y1$GPC_subtype==1)
    gpc.list=list('GPC1'=y1$Sample_name[wh],'GPC2'=y1$Sample_name[-wh])
  # Sig boxplot
    par(mfrow=c(1,2))
    for(i in 1:nrow(x1)){
    gpc.list2=llply(gpc.list,function(k){
      (as.numeric(x1[i,k]))
    })
    ylabel=paste0(x1$Symbol[i],' (Log2 (FPKM+1))')
    sig_boxplot(gpc.list2,colors = c('skyblue','tomato'),yaxis_lab = ylabel,
                title='Figure 4c-SMC1',digits = 2)
    }
    par(mfrow=c(1,1))
}

#===========================
# Figure 4c SMC2
#===========================
#x=smc2.datas
Figure4c.SMC2.20200406.1352=function(x=smc2.datas){
  # RNAs
    x1=x$exp
    y1=x$info
  # Set GPC lists
    wh=which(y1$sGPC=='sGPC1')
    smc2.sgpc.list=list('sGPC1'=y1$sample2[wh],'sGPC2'=y1$sample2[-wh])
  # Select genes
  x2=x1[which(x1$SYMBOL %in% c('CD274','PDCD1LG2')),]
  rownames(x2)=x2$SYMBOL
  par(mfrow=c(1,2))
  for(i in c('CD274','PDCD1LG2')){ #i='CD274'
    smc2.sgpc.list2=llply(smc2.sgpc.list,function(k){
      (as.numeric(x2[i,k]))
    })
    ylabel=paste0(i,' (Log2 (FPKM+1))')
    sig_boxplot(smc2.sgpc.list2,colors = c('skyblue','tomato'),yaxis_lab = ylabel,
                title='Figure 4c-SMC2',digits = 2,ns_visualize = T)
  }
  par(mfrow=c(1,1))
}

#===========================
# Figure 4c TCGA
#===========================
#x=tcga.datas
Figure4c.TCGA.20200406.1416=function(x=tcga.datas){
    x1=x$exp
    y1=x$sGPC.info
    y1$sGPC=ifelse(y1$sGPC1.p>y1$sGPC2.p,
                   'sGPC1','sGPC2')
  # Set sGPC subtypes
    wh=which(y1$sGPC=='sGPC1')
    tcga.sgpc.list=list('sGPC1'=y1$sample[wh],'sGPC2'=y1$sample[-wh])
    x2=x1[which(x1$symbol %in% c('CD274','PDCD1LG2')),]
    rownames(x2)=x2$symbol
  # Sig boxplot
    par(mfrow=c(1,2))
    for(i in c('CD274','PDCD1LG2')){
      tcga.sgpc.list2=llply(tcga.sgpc.list,function(k){
        as.numeric(x2[i,k])
      })
      ylabel=paste0(i,' (Log2 (FPKM+1))')
      sig_boxplot(tcga.sgpc.list2,colors = c('skyblue','tomato'),
                  yaxis_lab = ylabel,
                  title='Figure 4c-TCGA',digits = 2)
    }
    par(mfrow=c(1,1))
}

#===========================
# Figure 4c Yonsei
#===========================
#x=yonsei.datas
Figure4c.YONSEI.20200406.1428=function(x=yonsei.datas){
  x1=x$exp
  y1=x$info
  y2=y1[which(y1$min_P<0.05 & y1$Pathology_IDH=="IDH-wildtype"),]
  # Set sGPC subtypes
  wh=which(y2$subtype=='sGPC1')
  y.sgpc.list=list('sGPC1'=y2$sample[wh],'sGPC2'=y2$sample[-wh])
  x2=x1[which(x1$SYMBOL %in% c('CD274','PDCD1LG2')),]
  rownames(x2)=x2$SYMBOL
  # Sig boxplot
  par(mfrow=c(1,2))
  for(i in c('CD274','PDCD1LG2')){
    y.sgpc.list2=llply(y.sgpc.list,function(k){
      as.numeric(x2[i,k])
    })
    ylabel=paste0(i,' (Log2)')
    sig_boxplot(y.sgpc.list2,colors = c('skyblue','tomato'),
                yaxis_lab = ylabel,
                title='Figure 4c-Yonsei',digits = 2,ns_visualize = T)
  }
  par(mfrow=c(1,1))
}

#===========================
# Figure 4c ANOCEF
#===========================
#x=anocef.datas
Figure4c.ANOCEF.20200406.1447=function(x=yonsei.datas){
  x1=x$exp
  y1=x$info
  y2=y1[which(y1$min_P<0.05),]
  colnames(x1)=gsub(colnames(x1),pattern = '_',replacement = '')
  # Set sGPC subtypes
  wh=which(y2$sGPC=='sGPC1')
  a.sgpc.list=list('sGPC1'=y2$sample[wh],'sGPC2'=y2$sample[-wh])
  x2=x1[which(x1$NAME %in% c('CD274','PDCD1LG2')),]
  rownames(x2)=x2$NAME
  # Sig boxplot
  par(mfrow=c(1,2))
  for(i in c('CD274','PDCD1LG2')){
    a.sgpc.list2=llply(a.sgpc.list,function(k){
      as.numeric(x2[i,k])
    })
    ylabel=paste0(i,' (Log2)')
    sig_boxplot(a.sgpc.list2,colors = c('skyblue','tomato'),
                yaxis_lab = ylabel,
                title='Figure 4c-ANOCEF',digits = 2,ns_visualize = T)
  }
  par(mfrow=c(1,1))
}

#===========================
# Figure 4d
#===========================
#x=smc1.info;y=smc1.prot
Figure4d20200310.1447=function(x=smc1.info,y=smc1.prot){
  # Select IDH-WT GBMs
    x1=x[which(is.na(x$IDH1_mut) & x$Grade=='IV'),]
    x1$Sample_name2=paste0("'",x1$Sample_name)
  # Take max values for gene symbols
    y1=table_expand(y,key_colnames = 'Symbol')
    y2=dup.matrix(y1,key_column = 'Symbol',max)
  # Select PHGDH, RFTN2, FKBP9
    y2$Symbol=gsub(y2$Symbol,pattern = "'",replacement = '')
    y3=y2[which(y2$Symbol %in% c('PHGDH','RFTN2','FKBP9')),]
  # Take mean values for multi sector samples
    x2=x1[,c('Sample_name2','OS_from_initial_to_death_or_last_fu','Expired(dead=1)')]
    colnames(x2)=c('sample','surv','status')
    x2$patient=unlist(llply(x2$sample,function(k){
      unlist(strsplit(k,'-'))[1]
    }))
    dup.pat=x2$patient[which(duplicated(x2$patient))]
    dup.pat=unique(dup.pat)
    y4=y3
    for(i in dup.pat){ #i=dup.pat[1]
      wh=x2$sample[which(x2$patient==i)]
      wh2=which(colnames(y4) %in% wh)
      y4=cbind(y4[,-wh2],rowMeans(y4[,wh2]))
      colnames(y4)[ncol(y4)]=i
    }
    rownames(y4)=y4$Symbol
  # Conduct analysis
    x3=dup.matrix(x2,key_column = 'patient',max)
    x4=x3[,c('patient','surv','status')]
    colnames(x4)[1]='sample'
    y5=y4[,x4$sample]
    x5=sliding_log_rank_test(x = y5,y=x4)
  # Select minimum p-value point
    genes=unique(x5$gene)
    x6=ldply(genes,function(k){ #k=genes[1]
      k1=x5[which(x5$gene==k),]
      k2=k1[which.min(k1$p.value),]
      return(k2)
    })
  # Merge table
    y6=as.data.frame(t(y5),stringsAsFactors=F)
    y6$sample=rownames(y6)
    y7=mat.merge(y6,x4,by='sample')
  # Sort table
    y7=y7[,c('sample','PHGDH','RFTN2','FKBP9','surv','status')]
  # Draw survival plot with best cutoff
    par(mfrow=c(1,3))
    genes=c('PHGDH','RFTN2','FKBP9')
    l_ply(genes,function(k){ #k=genes[1]
      k1=x6[which(x6$gene==k),]
      k2=k1$cut_off
      # Set groups
       class=ifelse(y7[,k]>k2,'High','Low')
       survival_analysis(class=class,survival_time = y7$surv/30,status = y7$status,
                         colors = c('tomato','black'),title = k,xlab = 'Month')
    })
    par(mfrow=c(1,1))
}

#===========================
# Figure 4E
#===========================
#x=tcga.datas;y=tcga.exp.cli
Figure4e20200406.1448=function(x=tcga.datas){
  x1=x$exp;head(x1)[,1:5]
  y1=x$info
  # Select IDH-WT GBMs
  colnames(y1)=multi_sub(colnames(y1),pattern = 'patient',
                         replacement = 'Case',exact = T)
  y1$IDH.status=as.character(y1$IDH.status)
  y1$Grade=as.character(y1$Grade)
  y2=y1[which(y1$IDH.status=='WT' & y1$Grade=='G4' & !is.na(y1$Survival..months.)),]
  x2=x1[which(x1$symbol %in% c('PHGDH','RFTN2','FKBP9')),
        c(3,multi_grep(colnames(x1),pattern=y2$Case))]
  #x2=x1[which(x1$symbol %in% c('PHGDH','RFTN2','FKBP9')),]
  # Remove normal or control samples
  # Get average expression for the duplicates
  pat=unlist(llply(colnames(x2)[2:ncol(x2)],function(k){
    k1=unlist(strsplit(k,'-'))[1:3]
    k2=paste(k1,collapse = '-')
    return(k2)
  }))
  dup.pat=unique(pat[which(duplicated(pat))])
  for(i in dup.pat){ #i=dup.pat[1]
    wh=grep(colnames(x2),pattern = i)
    x2=cbind(x2[,-wh],rowMeans(x2[,wh]))
    colnames(x2)[ncol(x2)]=i
  }
  colnames(x2)[2:ncol(x2)]=unlist(llply(colnames(x2)[2:ncol(x2)],function(k){
    paste(unlist(strsplit(k,'-'))[1:3],collapse = '-')
  }))
  h1=x$Survival_cutoff_PMID28818916
  h2=h1[which(h1$Symbols %in% c('PHGDH','RFTN2','FKBP9')),]
  cut_offs=h2$`Expression Cutoffs`
  names(cut_offs)=h2$Symbols # PMID28818916, Table S6 glioma
  # Matrix merge
  rownames(x2)=x2$symbol
  x3=as.data.frame(t(x2[,-1]),stringsAsFactors=F)
  x3$Case=rownames(x3);head(x3)
  x4=mat.merge(y2,x3,by='Case');head(x4)
  x4$PHGDH=ifelse(2^x4$PHGDH-1>=cut_offs['PHGDH'],'High','Low')
  x4$RFTN2=ifelse(2^x4$RFTN2-1>=cut_offs['RFTN2'],'High','Low')
  x4$FKBP9=ifelse(2^x4$FKBP9-1>=cut_offs['FKBP9'],'High','Low');head(x4)
  # Drawing
  par(mfrow=c(1,3))
  for(i in c('PHGDH','RFTN2','FKBP9')){
    survival_analysis(class=x4[,i],
                      survival_time = as.numeric(x4$Survival..months.),
                      status=as.numeric(x4$Vital.status..1.dead.),
                      xlab = 'Month',title = paste0('Figure4e-',i),
                      colors = c('black','tomato'))
  }
  par(mfrow=c(1,1))
}

#===========================
# Figure 4f
#===========================
#x=smc1.prot;y=smc1.info
Figure4f20200314.1559=function(x=smc1.prot,y=smc1.info){
  # Draw boxplots for PHGDH and RFTN2 and FKBP9
      x1=table_expand(x,key_colnames = 'Symbol')
    # Get max value for protein
      x2=dup.matrix(x1,key_column = 'Symbol',
                    my.Function = function(k){max(as.numeric(k))})
      x2$Symbol=gsub(x2$Symbol,pattern = "'",replacement = '')
    # Select genes
      x3=x2[which(x2$Symbol %in% c('PHGDH','RFTN2','FKBP9')),]
    # Select samples
      y1=y[which(is.na(y$IDH1_mut) & y$Grade=='IV'),]
    # Gpc.list
      wh=which(y1$GPC_subtype==1)
      gpc.list=list('GPC1'=y1$Sample_name[wh],'GPC2'=y1$Sample_name[-wh])
      colnames(x3)=gsub(colnames(x3),pattern = "'",replacement = '')
    # Sort table
      rownames(x3)=x3$Symbol
      x3=x3[c('PHGDH','RFTN2','FKBP9'),]
    # Drawing boxplot
      par(mfrow=c(1,3))
      for(i in 1:nrow(x3)){ #i=1
        gpc.list2=llply(gpc.list,function(k){
          as.numeric(x3[i,k])
        })
        print(x3$Symbol[i])
        sig_boxplot(gpc.list2,colors = c('skyblue','tomato'),one.sided = T,
                    title = x3$Symbol[i],yaxis_lab = 'Log2')
      }
      par(mfrow=c(1,1))
}

#===========================
# Figure 4g
#===========================
#x=smc.tma;y=smc.tma.cell.fraction.data
Figure4g20200406.1707=function(x=smc.tma,y=smc.tma.cell.fraction.data){
  # Remove normals
    x1=x[-grep(x$Patient_ID,pattern = '_N'),]
    x2=x1[which(x1$IDH_status=='WT' & x1$Grade=='4'),]
  
  # Get max value for redudant samples
    x2$Case=unlist(llply(x2$Patient_ID,function(k){
      unlist(strsplit(k,'_'))[1]
    }))
    x2$Case=gsub(x2$Case,pattern = '[A-Z]',replacement = '')
    x2$Case=paste0('P',x2$Case)
  # Get total PHGDH value
    y1=ldply('PHGDH+',function(k){
      k1=y[grep(rownames(y),pattern = k,fixed = T),]
      k2=colSums(k1)
      return(k2)
    })
    rownames(y1)='PHGDH+'
    y2=data.frame('sample'=colnames(y1),'PHGDH+'=as.numeric(y1[1,]))
  # Combine tables
    x3=mat.merge(x2,y2,by.x="Sample_ID",by.y='sample')
  # Get Max value for redundant samples
    x4=dup.matrix(x3,key_column = 'Case',my.Function = max)
  # Drawing PHGDH plot
    x4.1=x4[,c("PHGDH+/Nestin-",'PHGDH.')]
    x4.1$`PHGDH+/Nestin-`=gsub(x4.1$`PHGDH+/Nestin-`,pattern = "'",replacement = '')
    x4.1$`PHGDH+/Nestin-`=as.numeric(x4.1$`PHGDH+/Nestin-`)
    rownames(x4.1)=x4$Case
    x4.1=t(x4.1)
    x4.2=x4[,c('Case',"OS_from_initial_to_death_or_last_fu","Deceased")]
    colnames(x4.2)=c('sample','surv','status')
    x4.2$surv=as.numeric(x4.2$surv)/30
    x4.2$status=as.numeric(x4.2$status)
    x5=get_optExp_survPoint(x4.1,x4.2)
    x5=x5[which(x5$gene=='PHGDH.'),]
    x4.2$class=ifelse(x4.1['PHGDH.',]>x5$cut_off,'High','Low')
    survival_analysis(x4.2$class,survival_time = x4.2$surv,status = x4.2$status,
                      colors = c('black','tomato'),title = 'Figure4g',xlab = 'Month')
}

#==========================
# Figure4h
#==========================
#x=fig4h.data
Figure4h20200408.1742=function(x=fig4h.data){
  a1=x$HS683
  b1=x$SNU1105
  # PHGDH activity
    # HS683
    a2.1=a1$activity
    groups=c('Cont','NCT-502 35uM')
    groups2=llply(groups,function(k){
      as.numeric(a2.1[1,grep(colnames(a2.1),pattern = k)])
    })
    names(groups2)=groups
    a3=sig_barplot(groups2,
                title = 'HS683 PHGDH activity',yaxis_lab = 'PHGDH activity (mU/mg)')
    hs683.activity.p.value=t.test(groups2$Cont,groups2$`NCT-502 35uM`)$p.value
    # SNU1105
    b2.1=b1$activity
    groups=c('Cont','NCT-502 3.5uM')
    groups2=llply(groups,function(k){
      as.numeric(b2.1[1,grep(colnames(b2.1),pattern = k)])
    })
    names(groups2)=groups
    snu1105.activity.p.value=t.test(groups2$Cont,groups2$`NCT-502 3.5uM`)$p.value
    b3=sig_barplot(groups2,
                title = 'SNU1105 PHGDH activity',yaxis_lab = 'PHGDH activity (mU/mg)')
    library(ggplot2)
    check.n.install.lib('gridExtra')
    library(gridExtra)
    
  # invasion activity  
    a2.1=as.num.mat(a1$invasion)
    groups=c('Cont','NCT-502 35uM')
    groups2=llply(groups,function(k){
      as.numeric(colMedians(a2.1[,grep(colnames(a2.1),pattern = k)]))
    })
    names(groups2)=groups
    a3.2=sig_barplot(groups2,
                   title = 'HS683 invasion',yaxis_lab = 'Relative invasion length')
    # Calculate ANOVA p for HS683
    library(tidyr)
    a3.1=gather(as.data.frame(a2.1,stringsAsFactors=F),wells,values)
    a3.1$types=multi_sub(a3.1$wells,pattern = paste0('_',1:3),
                         replacement = rep('',3),partial = T)
    hs683.invasion.anovap=aov(values~types*wells,data=a3.1)
    hs683.invasion.anovap2=summary(hs683.invasion.anovap)
    
    # SNU1105 invasion
    b2.1=as.num.mat(b1$invasion)
    groups=c('Cont','NCT-502 3.5uM')
    groups2=llply(groups,function(k){
      as.numeric(colMedians(b2.1[,grep(colnames(b2.1),pattern = k)]))
    })
    names(groups2)=groups
    b3.2=sig_barplot(groups2,
                   title = 'SNU1105 invasion',yaxis_lab = 'Relative invasion length')
    b3.1=gather(as.data.frame(b2.1,stringsAsFactors=F),wells,values)
    b3.1$types=multi_sub(b3.1$wells,pattern = paste0('_',1:3),
                         replacement = rep('',3),partial = T)
    snu1105.invasion.anovap=aov(values~types*wells,data=b3.1)
    snu1105.invasion.anovap2=summary(snu1105.invasion.anovap)
  # Return the results
    plots=grid.arrange(a3,a3.2,b3,b3.2,nrow=2)
    return(list('Plot'=plots,
                'HS683_PHGDH_activity_t.test'=hs683.activity.p.value,
                'HS683_invasion_ANOVAp'=hs683.invasion.anovap2,
                'SNU1105_PHGDH_activity_t.test'=snu1105.activity.p.value,
                'SNU1105_invasion_ANOVAp'=snu1105.invasion.anovap2)
           )
}

#===========================
# Extended Figure 4a
#===========================
#x=yonsei.datas;y=anocef.datas
Extended_Figure4a20200404.1433=function(x=yonsei.datas,y=anocef.datas){
  message('There will be 3 yonsei plots and 3 anocef plots')
  par(mfrow=c(1,3))
  # Get gene tables
    a1=x$exp;a.1=x$info
    b1=y$exp;b.1=y$info
  # Select genes
    genes=c('PHGDH','RFTN2','FKBP9')
    a2=a1[which(a1$SYMBOL %in% genes),]
    a3=dup.matrix(a2,key_column ='SYMBOL',max)
    rownames(a3)=a3$SYMBOL
    b2=b1[which(b1$NAME %in% genes),]
    rownames(b2)=b2$NAME
    colnames(b2)=gsub(colnames(b2),pattern = '_',replacement = '')
  # Select IDH-WT GBMsamples in yonsei cohort
    a.2=a.1[which(a.1$Pathology_IDH=="IDH-wildtype" & a.1$Dx_final=='GBM'),]
  # Select samples with survival info
    a4=a3[,a.2$sample]
  # Take average values for redundant patient sample
    dup.patients=a.2$Patient_ID[which(duplicated(a.2$Patient_ID))]
    for(i in dup.patients){ #i=dup.patients[1]
      wh=which(a.2$Patient_ID==i)
      wh2=a.2$sample[wh]
      wh3=which(colnames(a4) %in% wh2)
      a4=cbind(a4[,-wh3],rowMeans(a4[,wh3]))
      colnames(a4)[ncol(a4)]=wh2[1]
    }
    a.3=a.2[which(a.2$sample %in% colnames(a4)),]
    b3=b2[,which(colnames(b2) %in% b.1$sample)]
    b.2=b.1[which(b.1$sample %in% colnames(b3)),]
  # Sort tables
    a4=a4[,order(colnames(a4))]
    a.4=a.3[order(a.3$sample),c("sample","duration_OS","event_OS")]
    b3=b3[,order(colnames(b3))]
    b.3=b.2[order(b.2$sample),c("sample","surv.m","status")]
    colnames(a.4)=colnames(b.3)=c('sample','surv','status')
    b.3[,-1]=as.num.mat(b.3[,-1])
    a5=sliding_log_rank_test(x=a4,y=a.4)
    b4=sliding_log_rank_test(x=b3,y=b.3)
  # Get optimal points based on P-values
    b5=ldply(genes,function(k){ #k='PHGDH'
      k1=b4[which(b4$gene==k),]
      k2=k1[which.min(k1$p.value),]
      return(k2)
    })
    a6=ldply(genes,function(k){ #k='PHGDH'
      k1=a5[which(a5$gene==k),]
      k2=k1[which.min(k1$p.value),]
      return(k2)
    })
  # Draw survival plot
    # Yonesi cohort #i=genes[1]
      for(i in genes){
        cut_off=a6$cut_off[which(a6$gene==i)]
        # Set groups
          groups=ifelse(as.numeric(a4[i,])>cut_off,'High','Low')
        # Set title
          title=paste0('Ex Fig4a-Yonsei (',i,')')
          survival_analysis(class=groups,survival_time = a.4$surv,status = a.4$status,
                            colors = c('black','tomato'),xlab = 'Month',title = title)
      }
    # ANOCEF cohort #i=genes[1]
      for(i in genes){
        cut_off=b5$cut_off[which(b5$gene==i)]
        # Set groups
        groups=ifelse(as.numeric(b3[i,])>cut_off,'High','Low')
        # Set title
        title=paste0('Ex Fig4a-ANOCEF (',i,')')
        survival_analysis(class=groups,survival_time = b.3$surv,status = b.3$status,
                          colors = c('black','tomato'),xlab = 'Month',title = title)
      }
    par(mfrow=c(1,1))
}

#===========================
# Extended Figure 4c
#===========================
#x=exfig4c.data
Extended_Figure4c20200408.1934=function(x=exfig4c.data){
  par(mfrow=c(1,1))
  colnames(x)[1]='Sample'
  x1=rowMeans(x[,-1])
  names(x1)=x$Sample
  group.list=list('PHGDH-low'=c('SNU466','KNS81','SNU201','SNU626','A172'),
                  'PHGDH-high'=c('HS683','SNU1105','T98G','U87MG'))
  group.list2=llply(group.list,function(k){
    x1[k]
  })
  sig_boxplot(group.list2,colors = c('deepskyblue4','firebrick1'),
              yaxis_lab = 'PHGDH activity (mU/mg)',
              title = 'Extended Figure 4c')
}

#===========================
# Extended Figure 4d
#===========================
#x=exfig4d.data
Extended_Figure4d20200411.1813=function(x=exfig4d.data){
  check.n.install.lib('matrixStats')
  # Draw barplot
    # SNU201
      a=x$SNU201
      groups=c('Vector','PHGDH')
      library(matrixStats)
      groups2=llply(groups,function(k){ #k=groups[1]
        colMedians(as.num.mat(a[,grep(colnames(a),pattern = k)]))
      })
      names(groups2)=groups
      groups2[[2]]=groups2[[2]]/mean(groups2[[1]])
      groups2[[1]]=groups2[[1]]/mean(groups2[[1]])
    # SNU201 bar plot
      snu201.plot=sig_barplot(groups2,yaxis_lab = 'Relative invasion length',
                  title = 'SNU201')
    # KNS81
      b=x$KNS81
      groups=c('Vector','PHGDH')
      library(matrixStats)
      groups2=llply(groups,function(k){ #k=groups[1]
        colMedians(as.num.mat(b[,grep(colnames(b),pattern = k)]))
      })
      names(groups2)=groups
      groups2[[2]]=groups2[[2]]/mean(groups2[[1]])
      groups2[[1]]=groups2[[1]]/mean(groups2[[1]])
      # KNS81 bar plot
      kns81_plot=sig_barplot(groups2,yaxis_lab = 'Relative invasion length',
                  title = 'KNS81')
  # Calculate ANOVA p
    a1=gather(a,wells,values)
    a1$types=multi_sub(a1$wells,pattern = paste0('_',1:3),replacement = rep('',3),
                       partial = T)
    a2=aov(values~types*wells,data=a1)
    a3=summary(a2)
    b1=gather(b,wells,values)
    b1$types=multi_sub(b1$wells,pattern = paste0('_',1:3),replacement = rep('',3),
                       partial = T)
    b2=aov(values~types*wells,data=b1)
    b3=summary(b2)
  # Return result
    library(gridExtra)
    plots=grid.arrange(snu201.plot,kns81_plot,nrow=2)
    x1=list('Plot'=plots,'SNU201_ANOVA_test'=a3,'KNS81_ANOVA_test'=b3)
    return(x1)
}
#===========================
# Extended Figure 4e
#===========================
#x=ccle.exp;y=ccle.met;z=ccle.mut
#x=ccle.datas
Extended_Figure4e20200322.1448=function(x=ccle.datas){
  y=x$metabolome
  z=x$mutation
  x=x$exp
  # Select IDH-WT cell lines
    check.n.install.lib('tidyr')
    z1=z[,c("Hugo_Symbol","Variant_Classification","Tumor_Sample_Barcode")]
    z1=z1[which(!z1$Variant_Classification %in% c('Silent',"Splice_Site")),]
    z1$dup=paste0(z1$Hugo_Symbol,'@',z1$Tumor_Sample_Barcode)
    z2=dup.matrix(z1,key_column = 'dup',my.Function = function(k){paste(unique(k),collapse = ';')})
    z3=z2[,c("Hugo_Symbol","Variant_Classification","Tumor_Sample_Barcode")]
    z4=spread(z3,key="Hugo_Symbol",value = "Variant_Classification")
    z5=z4[which(is.na(z4$IDH1) & is.na(z4$IDH2)),c("Tumor_Sample_Barcode")]
    colnames(z4)[1]
  # GET PHGDH
    x1=x[which(x$Description=='PHGDH'),which(colnames(x) %in% z5)]
  
  # Get 2-HG level
    y1=as.data.frame(y[,c('CCLE_ID','2-hydroxyglutarate')],stringsAsFactors=F)
    rownames(y1)=y1$CCLE_ID
  # Get intersects
    intersect.cells=intersect(colnames(x1),y1$CCLE_ID)
  # Sort table
    x2=t(x1[,intersect.cells])
    y2=y1[intersect.cells,]
    identical(rownames(x2),rownames(y2))  
  # Draw scatter plot
    scatterPlot(x=log2(x2[,1]+1),y=y2[,2],xlabel = 'PHGDH (Log2(RPKM+1))',
                ylabel = '2-HG level',title = 'Extended Figure 4e')
}

#===========================
# Figure 4a
#===========================
#x=smc1.prot;y=smc1.info
Figure4a20200322.1614=function(x=smc1.prot,y=smc1.info){
  # Get maxvalue for gene symbols
    x$Symbol=gsub(x$Symbol,pattern = "'",replacement = '')
    x1=x[which(x$Symbol %in% c('NES','VIM','CD44',
                               'CLDN11','MOG',
                               'SLC1A2','SLC1A3','HEPACAM','ALDH1A1','S100B','GLUL')),]
    x2=dup.matrix(x1,key_column = 'Symbol',my.Function = max)
    rownames(x2)=x2$Symbol
  # Sort table
    y=y[order(y$GPC_subtype),]
    y1=y[which(is.na(y$IDH1_mut) & y$Grade=='IV'),]
    colnames(x2)=gsub(colnames(x2),pattern = "'",replacement = '')
    x3=x2[c('NES','VIM','CD44',
            'CLDN11','MOG',
            'SLC1A2','SLC1A3','HEPACAM','ALDH1A1','S100B','GLUL'),y1$Sample_name]
  # Make annotation table
    row.annot=data.frame(Marker=c(rep('Stem-EMT',3),
                                  rep('Oligo',2),
                                  rep('Astrocyte',6)
                         ))
    rownames(row.annot)=c('NES','VIM','CD44',
                          'CLDN11','MOG',
                          'SLC1A2','SLC1A3','HEPACAM','ALDH1A1','S100B','GLUL')
    col.annot=data.frame('GPC'=paste0('GPC',y1$GPC_subtype))
    rownames(col.annot)=y1$Sample_name
  # Draw plot
    library(pheatmap)
    gpc.col=c('GPC1'='skyblue','GPC2'='tomato')
    marker.col=c('Stem-EMT'='orange','Oligo'='lightblue',
                 'Astrocyte'='violet')
    pheatmap(x3,color = colorRampPalette(c('blue','white','red'))(31),
             scale = 'row',
             annotation_row = row.annot,
             annotation_col = col.annot,
             annotation_colors = list('GPC'=gpc.col,'Marker'=marker.col),
             cluster_cols = F,cluster_rows = F,
             main='Figure 4a')
}

#===========================
# Conduct t-test for bulk RNAseq data
#===========================
#x=smc1.rna;y=smc1.info
get.DEGs.btw.GPC.subtypes20200403.1247=function(x=smc1.rna,y=smc1.info){
  x=x$ENENSG_id
  # GPC lits
    y1=y[which(y$Grade=='IV' & is.na(y$IDH1_mut)),]
    wh.gpc1=which(y1$GPC_subtype==1)
    wh.gpc2=which(y1$GPC_subtype==2)
    gpc.list=list('GPC1'=y1$Sample_name[wh.gpc1],'GPC2'=y1$Sample_name[wh.gpc2])
    colnames(x)=gsub(colnames(x),pattern = "'",replacement = '')
  # Get p-value
    x=dup.matrix(x,key_colnames = 'ENSG_id',max)
    x1=matrix_test(mm=x,feature = 'ENSG_id',groups = gpc.list)
  # Return
    return(x1)
}

#==============================
# Draw figure5a
#==============================
#x=darmanis_sGPC_info
Figure5a20200403.1240=function(x=darmanis_sGPC_info){
  # Color sacle
  check.n.install.lib('circlize')
  library(circlize)
    sGPC1_col=circlize::colorRamp2(breaks = c(0,-log10(1/(1000*10))),colors = c('grey','skyblue'))
    sGPC2_col=circlize::colorRamp2(breaks = c(0,-log10(1/(1000*10))),colors = c('grey','tomato'))
    wh=which(x$del_z>0)
  x$col=NA
  x$col[wh]=sGPC1_col(-log10(x$p.val[wh]))
  x$col[-wh]=sGPC2_col(-log10(x$p.val[-wh]))
  # Drawing
    xrange=max(abs(range(x$del_z)))
    x$p.val[which(x$p.val==0)]=1E-4
    plot(x$del_z,-log10(x$p.val),pch=19,
         ylab='-log10(p-value)',xlab='GPC_delta_Zscore',
         main='Figure 5a',col=x$col,xlim=c(-xrange,xrange),
         cex.axis=1.5,cex.lab=1.5,cex.main=1.5,)
    abline(h=-log10(0.05),col='red',lwd=5,lty=3)
}

#==============================
# Draw figure5b
#==============================
#x=darmanis_sGPC_info;p_cut=0.05
Figure5b20200403.1450=function(x=sGPC.result,p_cut=0.05){
  # Set colors
    x$sGPC.col='grey'
    x$sGPC=ifelse(x$del_z>0,'sGPC1','sGPC2')
    x$sGPC[which(x$p.val>p_cut)]='Unknown'
    wh.gpc1=which(x$sGPC=='sGPC1')
    wh.gpc2=which(x$sGPC=='sGPC2')
    x$sGPC.col[wh.gpc1]='skyblue'
    x$sGPC.col[wh.gpc2]='tomato'
  # Overlay the cell on tSNE map
    plot(NULL,xlim=c(-70,60),ylim=c(-80,60),ylab='tSNE2',xlab='tSNE1',
         main='Figure5b')
    subtype=c('Unknown','sGPC1','sGPC2')
    for(i in subtype){
      wh=which(x$sGPC==i)
      points(x$Dim1[wh],x$Dim2[wh],col=x$sGPC.col[wh],pch=16,cex=0.5)
    }
}

#==============================
# Draw figure5c
#==============================
#x=darmanis_sGPC_info;y=single_cell_darmanis_symbol_updated_CPM
#gene=c('CD44','VIM','OLIG2','MOG','CLDN11','SLC1A2','S100B')
Figure5c20200403.1455=function(x=sGPC.result,y=single_cell_darmanis_symbol_updated_CPM,
                               p_cut=0.05,
                               gene=c('CD44','VIM','OLIG2','MOG','CLDN11','SLC1A2','S100B')){
  # Select genes
    y1=y[which(y$symbol %in% gene),] #head(y1[,1:10])
    rownames(y1)=y1$symbol
    y1=y1[gene,] #head(y1[,1:10])
  # Set gpc groups
    x$sGPC=ifelse(x$del_z>0,'sGPC1','sGPC2')
    x$sGPC[which(x$p.val>=p_cut)]='Unknown'
    print(table(x$sGPC))
    wh.gpc1=which(x$sGPC=='sGPC1')
    wh.gpc2=which(x$sGPC=='sGPC2')
    sgpc.list=list('sGPC1'=x$sample[wh.gpc1],'sGPC2'=x$sample[wh.gpc2])
  # Draw boxplot #i=1 
    par(mfrow=c(1,2))
    for(i in 1:nrow(y1)){
      sgpc.list2=llply(sgpc.list,function(k){
        as.numeric(y1[i,k])
      })
      ylabel=paste0(gene[i],' (Log2 (CPM))')
      title=paste0('Figure5c-',gene[i])
      sig_boxplot(input_list=sgpc.list2,test_method='wilcox.test',one.sided = F,
                 nsmall=3,digits=3,yaxis_lab=ylabel,colors=c('skyblue','tomato'),
                 title=title,jitter=F,violin=T,drawRect=F,
                 jitter.only=F)
    }
}

#==============================
# Draw figure5d
#==============================
#x=darmanis_sGPC_info;p_cut=0.05
Figure5d20200403.1747=function(x=darmanis_sGPC_info,p_cut=0.05){
  # Set sGPC subtype
    x$sGPC=ifelse(x$del_z>0,'sGPC1','sGPC2')
    x$sGPC[which(x$p.val>p_cut)]='Unknown'
  # Remove unknown
    x1=x[which(x$sGPC!='Unknown'),]
  # Tumor core
    x2=x1[which(x1$Location=='Tumor'),]
    x3=as.matrix(table(x2$sGPC,x2$Sample.name))
    par(mfrow=c(1,2))
  # Normalize and get percentage
    for(i in 1:ncol(x3)){
      x3[,i]=x3[,i]/sum(x3[,i])*100
    }
    barplot(x3,ylab = 'Percent (%)',col = c('skyblue','tomato'),las=2,
            legend.text = c('sGPC1','sGPC2'),main = 'Figure 6d-Tumor core')
  # Select neoplatic
    x2.neo=x1[which(x1$Cell_type=='Neoplastic'),]
    x3.neo=as.matrix(table(x2.neo$sGPC,x2.neo$Sample.name))
    for(i in 1:ncol(x3.neo)){
      x3.neo[,i]=x3.neo[,i]/sum(x3.neo[,i])*100
    }

    barplot(x3.neo,ylab = 'Percent (%)',col = c('skyblue','tomato'),las=2,
            legend.text = c('sGPC1','sGPC2'),main = 'Figure 6d-Neoplastic')
    par(mfrow=c(1,1))
}

#==============================
# Draw figure5e
#==============================
#x=darmanis_sGPC_info;p_cut=0.05
Figure5e20200403.1813=function(x=darmanis_sGPC_info,p_cut=0.05){
  # Set sGPC subtype
    x$sGPC=ifelse(x$del_z>0,'sGPC1','sGPC2')
    x$sGPC[which(x$p.val>p_cut)]='Unknown'
  # Remove unknown
    x1=x[which(x$sGPC!='Unknown'),]
  # Set distant as periphery
    x1$Location=gsub(x1$Location,pattern = 'Distant',replacement = 'Periphery')
  # Non neoplastic
    unique(x1$Cell_type)
    x1$Cell_type=as.character(x1$Cell_type)
    x2=x1[which(x1$Cell_type!='Neoplastic'),]
    x3=as.matrix(table(x2$sGPC,x2$Location))
    for(i in 1:ncol(x3)){
      x3[,i]=x3[,i]/sum(x3[,i])*100
    }
  # Draw pieplot
    check.n.install.lib('graphics',lib.type = 'cran')
    library(graphics)
    periphery=x3[,1]
    periphery.label=paste0(names(periphery),' (',round(periphery,0),'%)')
    
    core=x3[,2]
    core.label=paste0(names(core),' (',round(core,0),'%)')
    graphics::pie(periphery,radius = 1,col = c('darkslategray1','darkorange'),
        labels = periphery.label,init.angle = 90,main = 'Figure 6e-Non-neoplastic cells')
    par(new=T)
    graphics::pie(core,radius = 0.7,col = c('deepskyblue2','tomato'),
        labels = core.label,init.angle = 90)
}

#==============================
# Draw figure5f
#==============================
#x=darmanis_sGPC_info;y=single_cell_darmanis_symbol_updated_CPM
#gene=c('CD274')
Figure5f20200403.1830=function(x=darmanis_sGPC_info,
                               y=single_cell_darmanis_symbol_updated_CPM,
                               gene=c('CD274'),p_cut=0.05){
  # Set sGPC subtype
    x$sGPC=ifelse(x$del_z>0,'sGPC1','sGPC2')
    x$sGPC[which(x$p.val>p_cut)]='Unknown'
  # Select neoplastic cells
    x$Cell_type=as.character(x$Cell_type)
    x1=x[which(x$Cell_type=='Neoplastic'),]
  # Set sGPC subtypes
    wh.gpc1=which(x1$sGPC=='sGPC1')
    wh.gpc2=which(x1$sGPC=='sGPC2')
    sgpc.list.neo=list('sGPC1'=x1$sample[wh.gpc1],'sGPC2'=x1$sample[wh.gpc2])
  # Draw boxplot
    y1=y[which(y$symbol %in% gene),]
    for(i in 1:nrow(y1)){
      sgpc.list.neo2=llply(sgpc.list.neo,function(k){
        as.numeric(y1[which(y1$symbol==gene[i]),k])
      })
      ylabel=paste0(gene[i],' (Log2 (CPM))')
      title=paste0('Figure5f-',gene[i],'-Neoplastic')
      sig_boxplot(input_list=sgpc.list.neo2,test_method='wilcox.test',one.sided = F,
                  nsmall=3,digits=3,yaxis_lab=ylabel,colors=c('skyblue','tomato'),
                  title=title,jitter=F,violin=T,drawRect=F,
                  jitter.only=F)
    }
}

#===========================
# Extended Figure 5a
#===========================
#x=single_cell_darmanis_symbol_updated_CPM;y=darmanis_sGPC_info;z=dar.gpc.classfier
Extended_Figure5a20200407.1339=function(x=single_cell_darmanis_symbol_updated_CPM,
                                        y=darmanis_sGPC_info,
                                        z=dar.gpc.classfier,p_cut=0.05){
  x1=x[which(x$ensembl_gene_id %in% unlist(z)),]
  rownames(x1)=x1$ensembl_gene_id
  x1=x1[unlist(z),]
  # Sort table
  y=y[order(y$p.val,decreasing = F),]
  y$subtype=ifelse(y$del_z>0,'sGPC1','sGPC2')
  y$subtype[which(y$p.val>p_cut)]='Unknown'
  y1=ldply(c('sGPC1','Unknown','sGPC2'),function(k){
    k1=y[which(y$subtype==k),]
    return(k1)
  })
  x2=x1[,y1$sample]
  identical(y1$sample,colnames(x2))
  # Draw heatmap
  annot.col=data.frame('Predicted subtype'=y1$subtype)
  rownames(annot.col)=y1$sample
  annot.row=data.frame('Classifier genes'=c(rep('GPC1 high',12),rep('GPC2 high',12)))
  rownames(annot.row)=unlist(z)
  sGPC.col=c('sGPC1'='skyblue','sGPC2'='tomato','Unknown'='grey')
  gpc.high.col=c('GPC1 high'='dodgerblue4','GPC2 high'='red')
  wh=max(which(y1$subtype=='sGPC1'))
  wh2=min(which(y1$subtype=='sGPC2'))-1
  pheatmap(x2,cluster_rows = F,cluster_cols = F,annotation_col = annot.col,
           color = colorRampPalette(c('blue','white','red'))(51),
           annotation_row = annot.row,show_rownames = F,show_colnames = F,
           scale = 'row',gaps_col = c(wh,wh2),
           annotation_colors = list('Predicted.subtype'=sGPC.col,
                                    'Classifier.genes'=gpc.high.col),
           main='Extended Figure 5a')
}

#===========================
# Extended Figure 5b
#===========================
#x=smc.tma
Extended_Figure5b20200407.1402=function(x=smc.tma){
  x$`Nestin+/PHGDH-`=gsub(x$`Nestin+/PHGDH-`,pattern = "'",replacement = '')
  x$`PHGDH+/Nestin-`=gsub(x$`PHGDH+/Nestin-`,pattern = "'",replacement = '')
  x$`PHGDH+/Nestin-`=as.numeric(x$`PHGDH+/Nestin-`)
  x$`Nestin+/PHGDH-`=as.numeric(x$`Nestin+/PHGDH-`)
  # Draw annotation barplot
    # Sort according to Nestin+/PHGDH-
      x1=x[order(x$`Nestin+/PHGDH-`),]
      x1=x1[which(x1$Grade=='4' & x1$IDH_status=='WT'),]
      x1$sGPC_subtype.col='white'
      x1$sGPC_subtype.col[which(x1$sGPC_subtype=='sGPC1')]='skyblue'
      x1$sGPC_subtype.col[which(x1$sGPC_subtype=='sGPC2')]='tomato'
      tmp=t(data.frame('Nestin+/PHGDH-'=x1$`Nestin+/PHGDH-`,stringsAsFactors=F))
      annotated_barplot(mm = tmp,color.sample = x1$sGPC_subtype.col,color = 'grey',
                        xlab='IDH-WT GBM samples',ylab = 'Cell fraction (Nestin+, PHGDH-)',
                        color.sample.legend = c('sGPC1'='skyblue','sGPC2'='tomato'),
                        title = 'Extended Figure 5b top')
    # Sort according to Nestin+/PHGDH-
      x1=x[order(x$`PHGDH+/Nestin-`,decreasing = T),]
      x1=x1[which(x1$Grade=='4' & x1$IDH_status=='WT'),]
      x1$sGPC_subtype.col='white'
      x1$sGPC_subtype.col[which(x1$sGPC_subtype=='sGPC1')]='skyblue'
      x1$sGPC_subtype.col[which(x1$sGPC_subtype=='sGPC2')]='tomato'
      tmp=t(data.frame('PHGDH+/Nestin-'=x1$`PHGDH+/Nestin-`,stringsAsFactors=F))
      annotated_barplot(mm = tmp,color.sample = x1$sGPC_subtype.col,color = 'grey',
                        xlab='IDH-WT GBM samples',ylab = 'Cell fraction (Nestin-, PHGDH+)',
                        color.sample.legend = c('sGPC1'='skyblue','sGPC2'='tomato'),
                        title = 'Extended Figure 5b bottom')
    # Calculate sig boxplot
      wh.gpc1=which(x1$sGPC_subtype=='sGPC1')
      wh.gpc2=which(x1$sGPC_subtype=='sGPC2')
      group.list=list('sGPC1'=wh.gpc1,'sGPC2'=wh.gpc2)
      # Sort according to Nestin+/PHGDH-
      par(mfrow=c(1,2))
      group.list2=llply(group.list,function(k){
        x1$`Nestin+/PHGDH-`[k]
      })
      sig_boxplot(group.list2,test_method = 'wilcox.test',ns_visualize = F,
                  yaxis_lab = 'Cell fraction (Nestin+, PHGDH-)',
                  colors = c('skyblue','tomato'),one.sided = T,
                  title='Extended Figure 5b top')
      
      group.list2=llply(group.list,function(k){
        x1$`PHGDH+/Nestin-`[k]
      })
      sig_boxplot(group.list2,test_method = 'wilcox.test',ns_visualize = T,
                  yaxis_lab = 'Cell fraction (Nestin+, PHGDH-)',
                  colors = c('skyblue','tomato'),one.sided = T,
                  title='Extended Figure 5b bottom',digits = 3,nsmall = 5)
      par(mfrow=c(1,1))
}

#===========================
# Figure 6a prepare
#===========================
#x=smc1.prot;y=smc1.auc;z=smc1.info
cor.test_Drug.x.global20200326.1750=function(x=smc1.prot,y=smc1.auc,z=smc1.info){
  # Get samples drug response data
    z1=z[which(z$sample.type=='tumor'),]
  # Change colnames into proteome label]
    rownames(y)=y$Compound
    y1=as.num.mat(y[,-1])
    #head(y);head(y1)
  # Sort table
    y2=y1[,paste0("'",z1$Sample_name)]
    x1=x[,c('Uniprot',colnames(y2))]
  # Proteome data
  # Redunce redundancy
    rownames(x1)=x1$Uniprot
    x2=as.num.mat(x1[,-1])
    colnames(x2)==colnames(y2)
  # Conduct correlation test
  message('Conduct correlation test')
  cor.res=ldply(1:nrow(x2),.progress='text',function(k){ #k=j=1
    d.x=round(as.numeric(x2[k,]),5)
    cor.res=ldply(1:nrow(y2),function(j){
      d.y=round(as.numeric(y2[j,]),5)
      tmp=cor.test(d.x,d.y,method = 'spearman')
      return(c(tmp$estimate,tmp$p.value))
    })
    cor.res$drug=rownames(y1)
    cor.res$Proteins=rownames(x1)[k]
    cor.res$FDR=p.adjust(as.numeric(as.character(cor.res[,2])),'fdr')
    return(cor.res)
  })
  colnames(cor.res)[1:2]=c('Rho','Spear_pvalue')
  # Return the correlation test results
  cor.res=table_expand(cor.res,key_colnames = 'Proteins',sep = ';')
  return(cor.res)
}

#x=smc1.prot;y=smc1.auc;z=smc1.ed
Prepare.Figure6a20200322.1813=function(x=smc1.prot,y=smc1.auc,z=smc1.ed){
  # Calculate correlation
    message('Calculate correlation')
  tmp.auc=cor.test_Drug.x.global20200326.1750(x=smc1.prot,y=smc1.auc,z=smc1.info)
  tmp.ed=cor.test_Drug.x.global20200326.1750(x=smc1.prot,y=smc1.ed,z=smc1.info)
  # Return the result
    return(list('AUC.cor'=tmp.auc,'ED.cor'=tmp.ed))
}

#===========================
# Figure 6a left
#===========================
#x=drug.prot.cor
Figure6a.left.20200322.1842=function(x=drug.prot.cor){
  x1=x$AUC.cor
  y1=x$ED.cor
  z1=rbind(x1,y1)
  library(tidyr)
  #z2=table_expand(z1,key_colnames = 'protein',sep = ';')
  z1$dup=paste0(z1$drug,'@',z1$Spear_pvalue,'@',z1$FDR)
  z2=z1[which(!duplicated(z1$dup)),]
  z2$Spear_pvalue=as.numeric(z2$Spear_pvalue)
  hist(z2$Spear_pvalue,xlab='P',ylab='Count (X1,000)',breaks = 50,ylim = c(0,18*1000),yaxt='n',
       main='Figure6a left')
  axis(side=2,at=seq(0,17.5,by = 2.5)*1000,las=2)
}

#===========================
# Figure 6a right
#===========================
#x=drug.prot.cor
Figure6a.right.20200325.2142=function(x=drug.prot.cor){
  x1=x$AUC.cor
  y1=x$ED.cor
  x2=rbind(x1,y1)
  # Remove redunant samples
  x2$dup=paste0(x2$drug,'@',x2$Spear_pvalue,'@',x2$FDR)
  library(tidyr)
  x3=x2[which(!duplicated(x2$dup)),]
  #z2=table_expand(z1,key_colnames = 'protein',sep = ';')
  # Set colors
    drugs=sort(unique(x3$drug))
    drug.cols=rainbow(length(drugs))
    drug.cols2=data.frame('drug'=drugs,'colors'=drug.cols,stringsAsFactors=F)
    x4=mat.merge(x3,drug.cols2,by='drug')
  plot(x4$Rho, -log10(x4$FDR),xlab='Rho',ylab='-Log (FDR)',pch=16,col=x4$colors,
       main='Figure 6a right')
}

#===========================
# Figure 6b left protein
#===========================
#x=smc1.prot;y=smc1.auc
Figure6b.protein.20200326.2010=function(x=smc1.prot,y=smc1.auc){
  # Bortezomib and panobinostat
    y1=y[multi_grep(y$Compound,pattern = c('bortezomib','panobinostat')),]
  # Sort tables
    int.samples=intersect(colnames(y1),colnames(x))
    y2=y1[,int.samples]
    rownames(y2)=y1$Compound    
  # Take max values
    x1=x[,c('Symbol',int.samples)]
    x2=table_expand(x1,key_colnames = 'Symbol',sep=';')
    x3=dup.matrix(x2,key_column = 'Symbol',my.Function = function(k){max(as.numeric(k))})
    rownames(x3)=gsub(x3$Symbol,pattern = "'",replacement = '')
    x4=x3[,-1]
  # Calculate average protein expression for 20S proteasome complex
    target=c(paste0('PSMA',1:7),paste0('PSMB',1:7))
    wh=which(rownames(x4) %in% target)
    avg.proteasome=colMeans(x4[wh,])
    hdacs=c('HDAC1','HDAC2')
    wh=which(rownames(x4) %in% hdacs)
    avg.hdacs=colMeans(x4[wh,])
  # Draw scatterplot
    par(mfrow=c(2,1))
    scatterPlot(avg.proteasome,y=as.numeric(y2['Bortezomib',]),
                xlabel = 'Average protein expression (Log2) of 20S proteasomal subunits (N=14)',
                ylabel = 'Bortezomib (AUC)',
                title='Figure6b left top')
    
    scatterPlot(avg.hdacs,y=as.numeric(y2['Panobinostat',]),
                xlabel = 'Average protein expression (Log2) of HDAC1 and HDAC2',
                ylabel = 'Panobinostat (AUC)',
                title='Figure6b left bottom')
}

#===========================
# Figure 6b right mRNA
#===========================
#x=smc1.rna;y=smc1.auc
Figure6b.mRNA.20200326.2038=function(x=smc1.rna,y=smc1.auc){
  x=x$ENSG_id
  # Bortezomib and panobinostat
  y1=y[multi_grep(y$Compound,pattern = c('bortezomib','panobinostat')),]
  colnames(x)=multi_sub(colnames(x),pattern = 'symbol',replacement = 'Symbol',exact = T)
  colnames(y1)=gsub(colnames(y1),pattern = "'",replacement = "")
  # Sort tables
  int.samples=intersect(colnames(y1),colnames(x))
  y2=y1[,int.samples]
  rownames(y2)=y1$Compound    
  # Take max values
  x1=x[,c('Symbol',int.samples)]
  x2=dup.matrix(x1,key_column = 'Symbol',my.Function = function(k){max(as.numeric(k))})
  rownames(x2)=gsub(x2$Symbol,pattern = "'",replacement = '')
  x3=x2[,-1]
  # Calculate average protein expression for 20S proteasome complex
  target=c(paste0('PSMA',1:7),paste0('PSMB',1:7))
  wh=which(rownames(x3) %in% target)
  avg.proteasome=colMeans(x3[wh,])
  hdacs=c('HDAC1','HDAC2')
  wh=which(rownames(x3) %in% hdacs)
  avg.hdacs=colMeans(x3[wh,])
  # Draw scatterplot
  par(mfrow=c(2,1))
  scatterPlot(avg.proteasome,y=as.numeric(y2['Bortezomib',]),
              xlabel = 'Average mRNA expression (Log2) of 20S proteasomal subunits (N=14)',
              ylabel = 'Bortezomib (AUC)',
              title='Figure6b right top')
  
  scatterPlot(avg.hdacs,y=as.numeric(y2['Panobinostat',]),
              xlabel = 'Average mRNA expression (Log2) of HDAC1 and HDAC2',
              ylabel = 'Panobinostat (AUC)',
              title='Figure6b right bottom')
}

#===========================
# Figure 6c left delta AUC
#===========================
#x=smc1.auc;y=smc1.info
Figure6c.left20200404.1228=function(x=smc1.auc,y=smc1.info){
  # Select IDH-WT GBMs
    y1=y[which(is.na(y$IDH1_mut) & y$Grade=='IV'),]
    colnames(x)=gsub(colnames(x),pattern = "'",replacement = '')
    wh=which(y1$GPC_subtype==1)
    gpc.list=list('GPC1'=y1$Sample_name[wh],'GPC2'=y1$Sample_name[-wh])
  # Conduct ks-test
    x2=matrix_test(x,groups = gpc.list,test_method = 'ks.test',ignore.NA = T,feature = 'Compound')
  # Calculate delta AUC
    x2$delta=x2$GPC1_mean-x2$GPC2_mean
  # Make colors
    x2$col='black'
    x2$col[which(x2$delta<0 & x2$ks.test_pvalue<0.05)]='skyblue'
    x2$col[which(x2$delta>0 & x2$ks.test_pvalue<0.05)]='tomato'
  # Draw
    plot(x2$delta,-log10(x2$ks.test_pvalue),col=x2$col,pch=19,xlab='Delta_AUC',ylab='-Log (P)',
         main='Figure6c left')
    abline(v=0,h=-log10(0.05),lty=2,col='red')
  # Place texts
    wh=which(x2$ks.test_pvalue<0.05)
    text(x2$delta[wh]+7,-log10(x2$ks.test_pvalue[wh]),labels = x2$feature[wh])
}

#===========================
# Figure 6c right delta ED50
#===========================
#x=smc1.ed;y=smc1.info
Figure6c.right20200404.1237=function(x=smc1.ed,y=smc1.info){
  # Select IDH-WT GBMs
  y1=y[which(is.na(y$IDH1_mut) & y$Grade=='IV'),]
  colnames(x)=gsub(colnames(x),pattern = "'",replacement = '')
  wh=which(y1$GPC_subtype==1)
  gpc.list=list('GPC1'=y1$Sample_name[wh],'GPC2'=y1$Sample_name[-wh])
  # Conduct ks-test
  x2=matrix_test(x,groups = gpc.list,test_method = 'ks.test',ignore.NA = T,feature = 'Compound')
  # Calculate delta AUC
  x2$delta=x2$GPC1_mean-x2$GPC2_mean
  # Make colors
  x2$col='black'
  x2$col[which(x2$delta<0 & x2$ks.test_pvalue<0.05)]='skyblue'
  x2$col[which(x2$delta>0 & x2$ks.test_pvalue<0.05)]='tomato'
  # Draw
  plot(x2$delta,-log10(x2$ks.test_pvalue),col=x2$col,pch=19,xlab='Delta_ED50',ylab='-Log (P)',
       main='Figure6c right')
  abline(v=0,h=-log10(0.05),lty=2,col='red')
  # Place texts
  wh=which(x2$ks.test_pvalue<0.05 & x2$feature %in% c('AZD2014','LDE225 (NVP-LDE225, Erismodegib)'))
    # AZD2014 and erismodegib
    text(x2$delta[wh],-log10(x2$ks.test_pvalue[wh])-0.1,labels = x2$feature[wh])
    # Canertinib
    wh2=which(x2$ks.test_pvalue<0.05 & x2$feature %in% c('CI-1033 (Canertinib)'))
    text(x2$delta[wh2]-1.5,-log10(x2$ks.test_pvalue[wh2]),labels = x2$feature[wh2])
}

#===========================
# Figure 6d
#===========================
#x=smc1.prot;y=smc1.rna;z=smc1.info;a=msig;b=brcaness.formula
Figure6d20200404.1243=function(x=smc1.prot,y=smc1.rna,z=smc1.info,a=msig,b=brcaness.formula){
  y=y$ENSG_id
  # Select IDH-WT GBMs
    z1=z[which(z$Grade=='IV' & is.na(z$IDH1_mut)),]
    z1$Sample_name=paste0("'",z1$Sample_name)
    colnames(y)=multi_sub(colnames(y),pattern = 'symbol',replacement = 'Symbol',exact = T)
  # Get max symbol value
    x1=table_expand(x,key_colnames = 'Symbol',sep = ';')
    x2=dup.matrix(x1,key_column = 'Symbol',max)
    rownames(x2)=gsub(x2$Symbol,pattern = "'",replacement = "")
    x3=as.num.mat(x2[,z1$Sample_name])
  # Calculate ssGSEA for selected genesets
    a1=a$updated.GeneSet
    genesets=c("GO_PLATELET_DERIVED_GROWTH_FACTOR_RECEPTOR_BINDING",
               "SEIDEN_ONCOGENESIS_BY_MET",
               "GO_TRANSLATIONAL_INITIATION",
               "KEGG_ERBB_SIGNALING_PATHWAY",
               "PID_HEDGEHOG_GLI_PATHWAY")
    a2=a1[genesets]
    x4=ssgsea(mm=as.matrix(x3),gene_sets = a2)
  # Calculate BRCAness score
    y$Symbol=gsub(y$Symbol,pattern = "'",replacement = "")
    y1=y[which(y$Symbol %in% b$`Gene Symbol`),]
    y2=dup.matrix(y1,key_column = 'Symbol',max)
    b1=b[which(b$`Gene Symbol` %in% y2$Symbol),]
    rownames(b1)=b1$`Gene Symbol`
    rownames(y2)=y2$Symbol
    # Sort table
      y2=y2[rownames(b1),]
      identical(rownames(b1),rownames(y2))
      colnames(y2)[5:ncol(y2)]=paste0("'",colnames(y2)[5:ncol(y2)])
    # Calculate brcaness
      y3=y2
      for(i in 5:ncol(y3)){
        y3[,i]=y2[,i]*b1$Weight
      }
      y4=colSums(y3[,z1$Sample_name])
  # Make GPC list
    wh=which(z1$GPC_subtype==1)
    gpc.list=list('GPC1'=z1$Sample_name[wh],'GPC2'=z1$Sample_name[-wh])
  # Combine scores
    identical(colnames(x4),names(y4))
    x5=rbind(x4,y4)
    rownames(x5)[nrow(x5)]='BRCAness'
  # Sort table
    genesets2=c("GO_PLATELET_DERIVED_GROWTH_FACTOR_RECEPTOR_BINDING",
               "SEIDEN_ONCOGENESIS_BY_MET",
               'BRCAness',
               "GO_TRANSLATIONAL_INITIATION",
               "KEGG_ERBB_SIGNALING_PATHWAY",
               "PID_HEDGEHOG_GLI_PATHWAY")
    x6=x5[genesets2,]
  # Draw box-jitter plots
    ylabs=c('GO PDGFR binding (ssGSEA)',
            'Seiden oncogenesis by MET (ssGSEA)',
            'BRCAness score',
            'GO translational initiation (ssGSEA)',
            'KEGG ERBB signaling (ssGSEA)',
            'PID Hedgehog-GLI pathway (ssGSEA)')
    par(mfrow=c(1,3))
    for(i in 1:nrow(x6)){ #i=1
      gpc.list2=llply(gpc.list,function(k){
        as.numeric(x6[i,k])
      })
      sig_boxplot(gpc.list2,test_method = 'wilcox',yaxis_lab = ylabs[i],
                  colors = c('skyblue','tomato'),nsmall = 1,digits = 1,
                  title='Figure6d')
    }
    par(mfrow=c(1,1))
}

#===========================
# Figure 6e
#===========================
#x=smc1.prot;y=smc1.ed;z=smc1.info
Figure6e20200404.1310=function(x=smc1.prot,y=smc1.ed,z=smc1.info){
  # Select genes
    x1=x[multi_grep(x$Symbol,pattern = c('PHGDH','FKBP9')),]
    x1$Symbol=gsub(x1$Symbol,pattern = "'",replacement = '')
    x2=table_expand(x1,key_colnames = 'Symbol')
    x3=x2[which(x2$Symbol %in% c('PHGDH','FKBP9')),]
    rownames(x3)=x3$Symbol
  # Select drugs
    y1=y[which(y$Compound=='AZD2014'),]
  # Select IDH-WT GBMs
    z1=z[which(is.na(z$IDH1_mut) & z$Grade=='IV'),]
    z1$Sample_name=paste0("'",z1$Sample_name)
  # Sort table
    par(mfrow=c(1,2))
    x4=x3[,z1$Sample_name]
    y2=y1[,z1$Sample_name]
    scatterPlot(as.numeric(x4[2,]),log2(as.numeric(y2)),title='Figure6 left',
                xlabel ='PHGDH (Log2)',ylabel = 'AZD2014 (Log2 (uM))')
    scatterPlot(as.numeric(x4[1,]),log2(as.numeric(y2)),title='Figure6 right',
                xlabel ='FKBP9 (Log2)',ylabel = 'AZD2014 (Log2 (uM))')
    par(mfrow=c(1,1))
}

# Extended data figure 6a
#x=smc1.auc;a=smc1.prot;b=smc1.rna;d=smc1.cnv;e=drug.info
Prepare_Extended_Figure6a20200407.1739=function(x=smc1.auc,
                                        a=smc1.prot,
                                        b=smc1.rna,
                                        d=smc1.cnv,
                                        e=drug.info){
  b=b$ENSG_id
  colnames(b)=multi_sub(colnames(b),pattern = 'symbol',replacement = 'Symbol',exact=T)
  a$Symbol=gsub(a$Symbol,pattern = "'",replacement = '')
  b$Symbol=gsub(b$Symbol,pattern = "'",replacement = '')
  colnames(a)=gsub(colnames(a),pattern = "'",replacement = '')
  colnames(b)=gsub(colnames(b),pattern = "'",replacement = '')
  # Select intersecting genes
    a1=unlist(strsplit(a$Symbol,';'))
    b1=unlist(strsplit(b$Symbol,';'))
    d1=unlist(strsplit(d$Symbol,';'))
    int.genes=intersect(a1,b1)  
    int.genes2=intersect(int.genes,d1)
    #int.genes=intersect(a$Symbol,b1$Symbol)  
    #int.genes2=intersect(int.genes,d1$Symbol)
  # Select genes
    a1=table_expand(a,key_colnames = 'Symbol')
    a2=a1[which(a1$Symbol %in% int.genes2),]
    a3=dup.matrix(a2,key_column = 'Symbol',my.Function = function(k){max(as.numeric(as.character(k)))})
    b1=b[which(b$Symbol %in% int.genes2),]
    b2=dup.matrix(b1,key_column = 'Symbol',my.Function = function(k){max(as.numeric(as.character(k)))})
    d1=d[which(d$Symbol %in% int.genes2),]
    d1=dup.matrix(d1,key_column = 'Symbol',my.Function = function(k){max(as.numeric(as.character(k)))})
  # Get correlation coefficients
    # Sort table
      colnames(x)=gsub(colnames(x),pattern = "'",replacement = '')
      rownames(a3)=a3$Symbol
      rownames(b2)=b2$Symbol
      rownames(d1)=d1$Symbol
      a4=a3[int.genes2,colnames(x)[-1]]
      b3=b2[int.genes2,colnames(x)[-1]]
      d2=d1[int.genes2,colnames(x)[-1]]
      x1=x[,colnames(x)[-1]]
      #identical(colnames(d2),colnames(x1))
    # Conduct analysis #k=j=1
      message('Calculate correlation coefficients')
      cor.res=ldply(1:nrow(a4),.progress='text',function(k){
          p.gene=as.numeric(a4[k,])
          r.gene=as.numeric(b3[k,])
          c.gene=as.numeric(d2[k,])
        k1=ldply(1:nrow(x),function(j){
          drug=round(as.numeric(x1[j,]),5)
          wh=which(is.na(drug))
          if(length(wh)>0){
            jp=cor.test(p.gene[-wh],drug[-wh],method = 'spearman')$estimate
            jr=cor.test(r.gene[-wh],drug[-wh],method = 'spearman')$estimate
            jc=cor.test(c.gene[-wh],drug[-wh],method = 'spearman')$estimate
          }else{
            jp=cor.test(p.gene,drug,method = 'spearman')$estimate
            jr=cor.test(r.gene,drug,method = 'spearman')$estimate
            jc=cor.test(c.gene,drug,method = 'spearman')$estimate
          }
          c(jp,jr,jc)
        })
        k1$compound=x$Compound
        k1$gene=rownames(a4)[k]
        colnames(k1)[1:3]=c('prot','rna','cnv')
        return(k1)
      })
    # Attach targets
      e$Drug2=e$Drug
      if(F){
      e$Drug2=multi_sub(e$Drug,
                        pattern = c("ABT-199 (GDC-0199 )","AEE788 (NVP-AEE788)","Afatinib (BIBW2992)",
                                    "AUY922 (NVP-AUY922)","AZD6244 (Selumetinib)","BGJ398 (NVP-BGJ398)",
                                    "BKM120 (NVP-BKM120)","BMS-599626 (AC480)","Bortezomib (Velcade)",
                                    "Cabozantinib (XL184)","CI-1033 (Canertinib)","Cediranib (AZD2171)",
                                    "CO-1686 (Rociletinib)","Crizotinib (PF-02341066)","Dacomitinib (PF299804,PF-00299804)",
                                    "Dasatinib (BMS-354825)","Dovitinib (TKI-258)","LDE225 (NVP-LDE225, Erismodegib)",
                                    "Erlotinib HCl","Everolimus (RAD001)","Foretinib (XL880)",
                                    "Gefitinib (Iressa)","Imatinib (Gleevec)","Neratinib (HKI-272)",
                                    "Nilotinib (AMN-107)","Olaparib (AZD2281)","Panobinostat (LBH589)",
                                    "Pazopanib HCl","PD 0332991 (Palbociclib) HCl","PF-05212384 (PKI-587)",
                                    "Sunitinib Malate (Sutent)","Tandutinib (MLN518)","Tivozanib (AV-951)"),
                        
                        replacement = c("GDC-0199","AEE788","Afatinib",
                                        "AUY922","Selumetinib","BGJ398",
                                        "BKM120","BMS-599626","Bortezomib",
                                        "Cabozantinib","Canertinib","Cediranib",
                                        "CO-1686","Crizotinib","Dacomitinib",
                                        "Dasatinib","Dovitinib","Erismodegib",
                                        "Erlotinib","Everolimus","Foretinib",
                                        "Gefitinib","Imatinib","Neratinib",
                                        "Nilotinib","Olaparib","Panobinostat",
                                        "Pazopanib","Palbociclib","PKI-587",
                                        "Sunitinib","Tandutinib","Tivozanib"),
                        exact=T)
      }
    # Get Target and durg pairs
      cor.res$target='X'
      for(i in 1:nrow(e)){ #i=1
        wh=which(cor.res$compound==e$Drug2[i])
        targets=unlist(strsplit(e$Drug_target[i],','))
        wh2=which(cor.res$gene %in% targets)
        wh3=intersect(wh,wh2)
        if(length(wh3)>0){
          cor.res$target[wh3]='O'
        }
      }
      table(cor.res$target)
      #head(x$pro$spearman.rho)
      cor.res$prot=round(cor.res$prot,digits = 5)
      cor.res$cnv=round(cor.res$cnv,digits = 5)
      cor.res$rna=round(cor.res$rna,digits = 5)
    # Return the result
      return(cor.res)
}

# Extended data figure 6a
#x=ex.fig6a.cor.result
Extended_Figure6a20200407.1835=function(x=ex.fig6a.cor.result){
  par(mfrow=c(1,3))
  # Draw density plot for CNV drug target
    wh=which(x$target=='O')
    cnv.list=list('Others'=x$cnv[-wh],'Target'=x$cnv[wh])
    #sort(cnv.list$Target)==sort(on.target)
    #cnv.list=list('Target'=tmp.cnv$spearman.rho[wh],'Others'=tmp.cnv$spearman.rho[-wh])
    density_plot(cnv.list,xaxis_lab = 'Spearman Rho',colors = c('grey','tomato'),
                 title = 'Ex Figure 6a (DNA copy number)')
    c.pv=wilcox.test(cnv.list$Target,cnv.list$Others)$p.value
    mtext(paste0('Wilcox-test p = ',format(c.pv,scientific=T,nsmall=1,digits=1)))
  # Draw density plot for RNA drug target
    wh=which(x$target=='O')
    rna.list=list('Others'=na.omit(x$rna[-wh]),'Target'=na.omit(x$rna[wh]))
    density_plot(rna.list,xaxis_lab = 'Spearman Rho',colors = c('grey','tomato'),
                 title = 'Ex Figure 6a (mRNA level)')
    r.pv=wilcox.test(rna.list$Target,rna.list$Others)$p.value
    mtext(paste0('Wilcox-test p = ',format(r.pv,scientific=T,nsmall=1,digits=1)))
  # Draw density plot for RNA drug target
    wh=which(x$target=='O')
    prot.list=list('Others'=na.omit(x$prot[-wh]),'Target'=na.omit(x$prot[wh]))
    density_plot(prot.list,xaxis_lab = 'Spearman Rho',colors = c('grey','tomato'),
                 title = 'Ex Figure 6a (Protein level)')
    p.pv=wilcox.test(prot.list$Target,prot.list$Others)$p.value
    mtext(paste0('Wilcox-test p = ',format(p.pv,scientific=T,nsmall=1,digits=1)))  
    par(mfrow=c(1,1))
}

#===========================
# Extended figure 6b
#===========================
#x=smc1.phos;y=smc1.auc;z=smc1.ed
Prepare_Extended_Figure6b20200410.1428=function(x=smc1.phos,y=smc1.auc,z=smc1.ed){
  # Calculate phospho correlation
    # Sort tables
      int.samples=intersect(colnames(x),colnames(y))
      rownames(x)=x$repPsite
      x1=x[,int.samples]
      rownames(y)=y$Compound
      y1=y[,int.samples]
      rownames(z)=z$Compound
      z1=z[rownames(y1),int.samples]
      
    # Get values
      for(i in 1:ncol(y1)){
        y1[,i]=as.numeric(gsub(y1[,i],pattern = "'",replacement = ''))
        z1[,i]=as.numeric(gsub(z1[,i],pattern = "'",replacement = ''))
        x1[,i]=as.numeric(x1[,i])
      }
    # Conduct correlation test #k=j=1
      x2=ldply(1:nrow(x1),.progress='text',function(k){
        exprs=as.numeric(x1[k,])
        wh.exprs=which(!is.na(exprs))
        j1=ldply(1:nrow(z1),function(j){
          aucs=as.numeric(y1[j,])
          wh.aucs=which(!is.na(aucs))
          wh=intersect(wh.exprs,wh.aucs)
          auc.cor=cor.test(exprs[wh],aucs[wh],method='spearman')
          auc.rho=auc.cor$estimate
          auc.p=auc.cor$p.value
          
          eds=as.numeric(z1[j,])
          wh.eds=which(!is.na(eds))
          wh2=intersect(wh.exprs,wh.eds)
          ed.cor=cor.test(exprs[wh2],eds[wh2],method='spearman')
          ed.rho=ed.cor$estimate
          ed.p=ed.cor$p.value
          c(auc.rho,auc.p,ed.rho,ed.p)
        })
        colnames(j1)=c('auc.r','auc.p','ed.r','ed.p')
        j1$auc.fdr=p.adjust(j1$auc.p,'fdr')
        j1$ed.fdr=p.adjust(j1$ed.p,'fdr')
        j1$drugs=rownames(z1)
        j1$site=rownames(x1)[k]
        return(j1)
      })
      return(x2)
}

#===========================
# Extended figure 6b
#===========================
#x=phos.cor;y=smc1.phos
Extended_Figure6b20200410.1522=function(x=phos.cor){
  pv=c(x$auc.p,x$ed.p)
  hist(pv,breaks=50,ylim=c(0,14)*1000,yaxt='n',
       ylab='Count (X,1000)',xlab='P',main='Extended Figure 6b (left)')
  axis(side=2,at=seq(0,14,by=2)*1000,las=2,labels = seq(0,14,by=2))
  # Volcano plot
    # Get FDRs
      fdrs=c(x$auc.fdr,x$ed.fdr)
      rhos=c(x$auc.r,x$ed.r)
    # Set colors
      drug.col=rainbow(length(unique(x$drugs)))
      drug.cols=data.frame('drugs'=unique(x$drugs),'drug.col'=drug.col,stringsAsFactors=F)
    x1=mat.merge(x,drug.cols,by='drugs')
  plot(rhos,-log10(fdrs),col=x1$drug.col,pch=16,xlab='Rho',ylab='-Log (FDR)',
       main='Extended Figure 6b (right)')
}

#==========================
# Extended Figure 6c
#==========================
#x=smc1.phos;y=smc1.ed
Extended_Figure6c20200410.1552=function(x=smc1.phos,y=smc1.ed){
  # Select markers
  x1=x[multi_grep(x$SYMBOL,pattern = c('EGFR','SRC')),]
  x2=x1[multi_grep(x1$Feature,pattern = c('Y1197','S17','Y419')),]
  # Table expand
  x2$Feature=gsub(x2$Feature,pattern = "|",replacement = ';',fixed = T)
  # Select genes
  x3=table_expand(x2,key_colnames = 'Feature',sep = ';')
  features=c('P00533(Y1197)','P12931(S17)','P12931(Y419)')
  x4=x3[which(x3$Feature %in% features),]
  rownames(x4)=x4$Feature
  # Select drugs
  drugs=c('Afatinib','Bosutinib','Bosutinib')
  y1=y[multi_grep(y$Compound,pattern = drugs),]
  for(i in 1:ncol(y1)){
    y1[,i]=gsub(y1[,i],pattern = "'",replacement = '')
  }
  y1[,-1]=as.num.mat(y1[,-1])
  rownames(y1)=y1$Compound
  # Sort tables
  int.samples=intersect(colnames(x4),colnames(y1))
  x5=x4[,int.samples]
  y2=y1[,int.samples]
  # Get corrleations
  xlabels=paste0(c('EGFR-pY1197','SRC-pS17','SRC-pY419'), ' (Log2)')
  ylabels=paste0(drugs,' (ED50, uM)')
  par(mfrow=c(1,3))
  for(i in 1:length(drugs)){ #i=1
    x.value=as.numeric(x5[grep(rownames(x5),pattern = features[i],fixed = T),])
    y.values=as.numeric(y2[grep(rownames(y2),pattern = drugs[i]),])
    scatterPlot(x.value,y.values,xlabel = xlabels[i],ylabel = ylabels[i],
                title = 'Extended Figure 6c')
  }
}
#==========================
# Extended Figure 6d
#==========================
#x=smc1.phos;y=smc1.ed
Extended_Figure6d20200410.1653=function(x=smc1.phos,y=smc1.ed){
  # Select markers
  x1=x[multi_grep(x$SYMBOL,pattern = c('EGFR','RAB4B','SNAP91','ANK2')),]
  x2=x1[multi_grep(x1$Feature,pattern = c('T693','S193','T309','S2516')),]
  # Table expand #head(x2[,1:5])
  x2$Feature=gsub(x2$Feature,pattern = "|",replacement = ';',fixed = T)
  # Select genes
  x3=table_expand(x2,key_colnames = 'Feature',sep = ';')
  features=c('P00533(T693)','P61018(S193)','O60641(T309)','Q01484(S2516)')
  x4=x3[which(x3$Feature %in% features),]
  rownames(x4)=x4$Feature
  # Select drugs
  drugs=rep('Lapatinib',4)
  y1=y[multi_grep(y$Compound,pattern = drugs),]
  for(i in 1:ncol(y1)){
    y1[,i]=gsub(y1[,i],pattern = "'",replacement = '')
  }
  rownames(y1)=y1$Compound
  # Sort tables
  int.samples=intersect(colnames(x4),colnames(y1))
  x5=x4[,int.samples]
  y2=y1[,int.samples]
  # Get corrleations
  xlabels=paste0(c('EGFR-pT693','RAB4B-pS193',
                   'SNAP91-pT309','ANK2-pS2516'), ' (Log2)')
  ylabels=paste0(drugs,' (ED50, uM)')
  par(mfrow=c(2,2))
  for(i in 1:length(drugs)){ #i=2
    x.value=as.numeric(x5[grep(rownames(x5),pattern = features[i],fixed = T),])
    y.values=as.numeric(y2[grep(rownames(y2),pattern = drugs[i]),])
    scatterPlot(x.value,y.values,xlabel = xlabels[i],ylabel = ylabels[i],
                title = 'Extended Figure 6d',cor.method = 'spearman',)
  }
}

#==========================
# Extended Figure 6e
#==========================
#x=smc1.phos;y=smc1.ed;z=smc1.auc
Extended_Figure6e20200410.1701=function(x=smc1.phos,y=smc1.ed,z=smc1.auc){
  # Select markers
  x1=x[multi_grep(x$SYMBOL,pattern = c('TSC2','TOMM22')),]
  x2=x1[multi_grep(x1$Feature,pattern = c('S1420','S45')),]
  # Table expand #head(x2[,1:5])
  x2$Feature=gsub(x2$Feature,pattern = "|",replacement = ';',fixed = T)
  # Select genes
  x3=table_expand(x2,key_colnames = 'Feature',sep = ';')
  features=c('P49815(S1420)','Q9NS69(S45)')
  x4=x3[which(x3$Feature %in% features),]
  rownames(x4)=x4$Feature
  # Select drugs
  y1=y[multi_grep(y$Compound,pattern = 'Bosutinib'),]
  for(i in 1:ncol(y1)){
    y1[,i]=gsub(y1[,i],pattern = "'",replacement = '')
  }
  rownames(y1)=y1$Compound
  z1=z[multi_grep(z$Compound,pattern = 'AZD4547'),]
  for(i in 1:ncol(z1)){
    z1[,i]=gsub(z1[,i],pattern = "'",replacement = '')
  }
  rownames(z1)=z1$Compound
  # Sort tables
  int.samples=intersect(colnames(x4),colnames(y1))
  x5=x4[,int.samples]
  y2=y1[,int.samples]
  z2=z1[,int.samples]
  # Get corrleations
  par(mfrow=c(1,2))
  scatterPlot(as.numeric(x5['P49815(S1420)',]),as.numeric(y2['Bosutinib',]),
              title='Extended Figure 6e',
              xlabel = 'TSC2-pS1420 (Log2)',
              ylabel = 'Bosutinib (ED50, uM)')
  scatterPlot(as.numeric(x5['Q9NS69(S45)',]),as.numeric(z2['AZD4547',]),
              title='Extended Figure 6e',
              xlabel = 'TOMM22-pS45 (Log2)',
              ylabel = 'AZD4547 (AUC)')
}
