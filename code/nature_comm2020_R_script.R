# Set a working directory
wd='C:/Users/osj118/Desktop/2015gbm_revision_code/Source_files/3rd' 
setwd(wd)

# Load R functions
source('nature_comm2020_R_functions.r')
source('nature_comm2020_R_sub.R')

# Download TCGA data
tcga.exp=Download.tcga.expr.data20200416.1405()
tcga.exp2=attach.gene.symbol2TCGA.exp20190117.1525(x=tcga.exp)

# Preprocessing darmanis single cell data
updata.darmanis.geneIDs20200416.1641(dar='Darmanis_RNAseq_data')->single_cell_darmanis_symbol_updated_CPM

# Load datas
load.data20200226.0011(wd=wd)
# Prepare datas
ccle.datas=list('exp'=Here_ccle_expression_data,
                'metabolome'=Here_metabolome_data,
                'mutation'=Here_mutation_data,
                'sGPC.info'=ccle.datas$sGPC.info)
tcga.datas=list('exp'=tcga.exp2,
                'info'=tcga.datas$info,
                'sGPC.info'=tcga.datas$sGPC.info,
                'Survival_cutoff_PMID28818916'=tcga.datas$Survival_cutoff_PMID28818916)

#=======================
# Figure 1
#=======================
# Figure 1a
pdf('Fig1a_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
Figure1a20200409.1046(x=smc1.snv,y=smc1.info)
dev.off()

# Figure 1c
pdf('Fig1c_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
Figure1c20200409.1354(x=smc1.snv,y=smc1.snv.by.mass)
dev.off()

# Figure 1d
pdf('Fig1d_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
Figure1d.20200226.0016(x=smc1.prot,y=smc1.rna)
dev.off()

# Figure 1e # two plots : dendrogram and annotation heatmap
pdf('Fig1e_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
Figure1e.20200226.0101(x=smc1.prot,y=smc1.info) 
dev.off()

# Permutation test to evaluate clustering for selected features
# RNA subtypes, clinical phenotypes, 3~5 minutes
permutation.result=perm.test.4.features20200226.2044(x=smc1.prot,y=smc1.info)

# Extended Figure 1a # 3~5minutes
prepare_Extended_Figure1a20200412.1717(x=smc1.cnv)

# Extended Figure 1b # 1 minute
pdf('Extended_Fig1b_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
Extended_Figure1b20200408.2057(x=smc1.rna,y=msig)
dev.off()

#===============================
# Figure 2
#===============================
# Figure 2a
pdf('Fig2a_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
  # PAC score barplot
    Figure2a20200228.0104(x=smc1.prot,y=smc1.info)->F2.heatmap 
  # Heatmap
    plot(F2.heatmap$gtable)
dev.off()

# Figure 2b
pdf('Fig2b_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
Figure2b20200228.0159(x=smc1.prot,y=smc1.info)
dev.off()

# Figure 2c
pdf('Fig2c_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
Figure2c20200228.0715(x=smc1.info)
dev.off()

# Figure 2d # This may require 2 ~ 3 minutes
pdf('Fig2d_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
Figure2d20200228.0741(x=smc1.info)
Figure2d.lolliplot20200228.0830(x=smc1.info,y='PIK3CA')
Figure2d.lolliplot20200228.0830(x=smc1.info,y='EGFR')
dev.off()

# Extended data Figure 2a
Extended_Figure2a20200406.1932(x=smc1.prot,y=smc1.info)->ex2a
pdf('Extended_Fig2a_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
plot(ex2a$heatmap$gtable)
dev.off()

# Extended data Figure 2b
Extended_Figure2b20200406.1935(x=smc1.prot,y=smc1.info)->ex2b
pdf('Extended_Fig2b_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
plot(ex2b$heatmap$gtable)
dev.off()

# Extended data Figure 2c
Extended_Figure2c20200406.1947(x=smc1.prot,y=smc1.info)->ex2c
pdf('Extended_Fig2c_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
plot(ex2c$heatmap$gtable)
dev.off()

# Extended data Figure 2d
pdf('Extended_Fig2d_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
Extended_Figure2d20200406.1950(x=smc1.info,a=ex2a,b=ex2b,d=ex2c)
dev.off()

#===============================
# Figure 3
#===============================
# Figure 3a
  # Figure 3a left
    pdf('Fig3a_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
    Figure3a.PCAplot20200309.2120(x=smc1.prot,y=smc1.info)->pc1.top.ten
  # Figure 3a right 1~2minutes
    Figure3a.barPlot20200309.2120(x=pc1.top.ten,
                                 y=corum.v3.0.proteinComplex,
                                 z=msig)->pc1.gea.result
    dev.off()

# Figure 3b
  # Figure 3b left
    pdf('Fig3b_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
        Figure3b.heatmap20200309.2235(x=smc1.prot,y=smc1.info,
                                       z=msig,a=smc1.rna)
  # Figure 3b right
        Figure3b.boxplot20200309.2352(x=smc1.prot,y=smc1.info,
                                       z=msig,
                                       a=smc1.rna)
        dev.off()
# Figure 3d 1~2 minutes
  pdf('Fig3d_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
  Figure3d.PKM.boxplot20200310.0031(x=set1.peptide,y=set2.peptide,z=smc1.info)
  dev.off()
  
# Figure 3e
  pdf('Fig3e_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
  Figure3e20200407.1627(x=ccle.datas,p_cut=0.05)->heatmap.fig3e
  dev.off()
  
# Extended Figure 3b
  pdf('Extended_Fig3b_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
  Extended_Figure3b20200404.1544(x=smc1.rna,y=smc1.info)
  dev.off()
# Extended Figure 3c
  pdf('Extended_Fig3c_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
  Extended_Figure3c20200407.1326(x=smc1.prot,y=smc1.info)
  dev.off()
#===============================
# Figure 4
#===============================
# Figure 4a
  pdf('Fig4a_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
  Figure4a20200322.1614(x=smc1.prot,y=smc1.info)
  dev.off()
  
# Figure 4c
  pdf('Fig4c_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
  Figure4c.SMC1.20200404.1718(x=smc1.rna,y=smc1.info)
  Figure4c.SMC2.20200406.1352(x=smc2.datas)
  Figure4c.TCGA.20200406.1416(x=tcga.datas)
  Figure4c.YONSEI.20200406.1428(x=yonsei.datas)
  Figure4c.ANOCEF.20200406.1447(x=anocef.datas)
  dev.off()
  
# Figure 4d # 1 minute
  pdf('Fig4d_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
  Figure4d20200310.1447(x=smc1.info,y=smc1.prot)
  dev.off()
  
# Figure 4e
  pdf('Fig4e_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
  Figure4e20200406.1448(x=tcga.datas)
  dev.off()
  
# Figure 4f # 1 ~ 2minutes
  pdf('Fig4f_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
  Figure4f20200314.1559(x=smc1.prot,y=smc1.info)
  dev.off()
  
# Figure 4g
  pdf('Fig4g_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
  Figure4g20200406.1707(x=smc.tma,y=smc.tma.cell.fraction.data)
  dev.off()
# Figure 4h
  pdf('Fig4h_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
  Figure4h20200408.1742(x=fig4h.data)->Figure4h_plot_n_stats
  dev.off()
  # HS683 PHGDH activity t.test result
    Figure4h_plot_n_stats$HS683_PHGDH_activity_t.test
  # HS683 invasion ANOVA test results
    Figure4h_plot_n_stats$HS683_invasion_ANOVAp
  # SNU1105 PHGDH activity t.test result
    Figure4h_plot_n_stats$SNU1105_PHGDH_activity_t.test
  # SNU1105 invasion ANOVA test results
    Figure4h_plot_n_stats$SNU1105_invasion_ANOVAp

# Extended data figure 4a
  pdf('Extended_Fig4a_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
  Extended_Figure4a20200404.1433(x=yonsei.datas,y=anocef.datas)
  dev.off()
  
# Extended data figure 4c
  pdf('Extended_Fig4c_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
  Extended_Figure4c20200408.1934(x=exfig4c.data)
  dev.off()
  
# Extended data figure 4d
  pdf('Extended_Fig4d_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
  Extended_Figure4d20200411.1813(x=exfig4d.data)
  dev.off()
  
# Extended data figure 4e # 1 ~ 3 minutes
  pdf('Extended_Fig4e_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
  Extended_Figure4e20200322.1448(x=ccle.datas)
  dev.off()

#===========================
# Figure 5 
#===========================
# Figure 5a
pdf('Fig5a_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
Figure5a20200403.1240(x=darmanis_sGPC_info)
dev.off()

# Figure 5b
pdf('Fig5b_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
Figure5b20200403.1450(x=darmanis_sGPC_info,p_cut=0.05)
dev.off()

# Figure 5C
pdf('Fig5c_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
Figure5c20200403.1455(x=darmanis_sGPC_info,
                      y=single_cell_darmanis_symbol_updated_CPM,
                      p_cut=0.05,
                      gene=c('CD44','VIM','OLIG2','MOG','CLDN11','SLC1A2','S100B'))
dev.off()

# Figure 5d
pdf('Fig5d_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
Figure5d20200403.1747(x=darmanis_sGPC_info,p_cut=0.05)
dev.off()

# Figure 5e
pdf('Fig5e_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
Figure5e20200403.1813(x=darmanis_sGPC_info,p_cut=0.05)
dev.off()

# Figure 5f
pdf('Fig5f_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
Figure5f20200403.1830(x=darmanis_sGPC_info,
                    y=single_cell_darmanis_symbol_updated_CPM,
                    gene=c('CD274'))
dev.off()

# Extended Figure 5a
pdf('Extended_Fig5a_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
Extended_Figure5a20200407.1339(x=single_cell_darmanis_symbol_updated_CPM,
                              y=darmanis_sGPC_info,
                              z=dar.gpc.classfier,p_cut=0.05)
dev.off()

# Extended Figure 5b 
pdf('Extended_Fig5b_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
Extended_Figure5b20200407.1402(x=smc.tma)
dev.off()

#===========================
# Figure 6
#===========================
# Figure 6a
pdf('Fig6a_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
drug.prot.cor=Prepare.Figure6a20200322.1813(x=smc1.prot,y=smc1.auc,z=smc1.ed) # 3 ~ 5 minutes
Figure6a.left.20200322.1842(x=drug.prot.cor) # 30 seconds
Figure6a.right.20200325.2142(x=drug.prot.cor) # 2 minutes
dev.off()

# Figure 6b
pdf('Fig6b_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
Figure6b.protein.20200326.2010(x=smc1.prot,y=smc1.auc)
Figure6b.mRNA.20200326.2038(x=smc1.rna,y=smc1.auc)
dev.off()

# Figure 6c
pdf('Fig6c_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
Figure6c.left20200404.1228(x=smc1.auc,y=smc1.info)
Figure6c.right20200404.1237(x=smc1.ed,y=smc1.info)
dev.off()

# Figure 6d # Sixplots will be made
pdf('Fig6d_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
Figure6d20200404.1243(x=smc1.prot,y=smc1.rna,z=smc1.info,a=msig,b=brcaness.formula)
dev.off()

# Figure 6e
pdf('Fig6e_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
Figure6e20200404.1310(x=smc1.prot,y=smc1.ed,z=smc1.info)
dev.off()

# Extended Figure 6a # 3~5minutes
pdf('Extended_Fig6a_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
ex.fig6a.cor.result=Prepare_Extended_Figure6a20200407.1739(x=smc1.auc,
                                      a=smc1.prot,b=smc1.rna,
                                      d=smc1.cnv,e=drug.info)
Extended_Figure6a20200407.1835(x=ex.fig6a.cor.result)
dev.off()

# Extended Figure 6b # 3 ~ 5 minutes
pdf('Extended_Fig6b_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
phos.cor=Prepare_Extended_Figure6b20200410.1428(x=smc1.phos,y=smc1.auc,z=smc1.ed)
Extended_Figure6b20200410.1522(x=phos.cor)
dev.off()

# Extended Figure 6c
pdf('Extended_Fig6c_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
Extended_Figure6c20200410.1552(x=smc1.phos,y=smc1.ed)
dev.off()  

# Extended Figure 6d
pdf('Extended_Fig6d_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
Extended_Figure6d20200410.1653(x=smc1.phos,y=smc1.ed)
dev.off()

# Extended Figure 6e
pdf('Extended_Fig6e_expected_result.pdf',colormodel = 'rgb',paper = 'letter')
Extended_Figure6e20200410.1701(x=smc1.phos,y=smc1.ed,z=smc1.auc)
dev.off()
