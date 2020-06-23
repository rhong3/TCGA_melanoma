# Data preparation
library(readxl)
library(biomaRt)
library(org.Hs.eg.db)
library(stringr)
library(readr)
library(plyr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)


SAAV <- read_excel("Data/Melanoma proteome/MM List of peptides with SAAVs.xlsx")
SAAV$UNIPROT = SAAV$Canonical.Accession
annotation <- select(org.Hs.eg.db, keys=unique(SAAV$UNIPROT), columns=c('UNIPROT', 'SYMBOL', 'ENTREZID'), keytype="UNIPROT")
SAAV.out = merge(SAAV, annotation, by="UNIPROT")

write.csv(SAAV.out, 'Data/Melanoma proteome/SAAV.csv', row.names = FALSE)



SAAV <- read.csv("~/Documents/TCGA_melanoma/Data/Melanoma proteome/SAAV.csv")
rnaseq <- read.delim("~/Documents/TCGA_melanoma/Data/skcm_tcga_pan_can_atlas_2018/data_RNA_Seq_v2_expression_median.txt")
rnaseq[, 3:445] = log2(as.matrix(rnaseq[, 3:445])+1)
library(matrixStats)
library(dplyr)
rnaseq[, 3:445]  = mutate_all(rnaseq[, 3:445] , function(x) as.numeric(as.character(x)))
write.csv(rnaseq, 'Data/log2_RNAseq.csv', row.names = FALSE)
rnaseq.saav = rnaseq[(rnaseq$Hugo_Symbol %in% unique(SAAV$SYMBOL)),]

clinical <- read.delim("~/Documents/TCGA_melanoma/Data/skcm_tcga_pan_can_atlas_2018/data_clinical_sample.txt", comment.char="#")
clinical$SAMPLE_ID=gsub('-','.', clinical$SAMPLE_ID)
clinical = clinical[ , c('PATIENT_ID', 'SAMPLE_ID', 'CANCER_TYPE_DETAILED', 'ANEUPLOIDY_SCORE', 'SAMPLE_TYPE')]
patient = read.delim("~/Documents/TCGA_melanoma/Data/skcm_tcga_pan_can_atlas_2018/data_clinical_patient.txt", comment.char="#")
patient = patient[, c('PATIENT_ID', 'AGE', 'SEX', 'AJCC_PATHOLOGIC_TUMOR_STAGE', 'WEIGHT', 'DSS_MONTHS')]

patient.clinical = merge(patient, clinical, by="PATIENT_ID")
row.names(patient.clinical) = patient.clinical$SAMPLE_ID
patient.clinical = patient.clinical[,-c(1, 7)]
patient.clinical.trans = data.frame(t(patient.clinical))

row.names(rnaseq.saav) = rnaseq.saav$Hugo_Symbol
rnaseq.saav = rnaseq.saav[, -c(1,2)]

patient.clinical.trans = patient.clinical.trans[, colnames(patient.clinical.trans) %in% colnames(rnaseq.saav)]

combined = rbind(rnaseq.saav, patient.clinical.trans)
combined = rbind(combined[748:755,], combined[1:747,])
write.csv(combined, 'Data/SAAV_RNAseq.csv', row.names = TRUE)

# col_fun1 = colorRamp2(c(0, 0.05), c("white", "green"))
# col_fun2 = colorRamp2(c(0, 10000), c("white","blue"))
# col_fun3 = colorRamp2(c(1,4), c("white", "orange"))
# col_fun4 = colorRamp2(c(0, 1), c("white", "purple"))
gn = rowAnnotation(gene.name = anno_text(row.names(combined)[9:755],
                                         location = 0.5, just = "center"))


anno = HeatmapAnnotation(age = as.numeric(as.matrix(patient.clinical.trans[1,])), 
                         gender = as.character(as.matrix(patient.clinical.trans[2,])), 
                         AJCC_stage = as.character(as.matrix(patient.clinical.trans[3,])),
                         weight = as.numeric(as.matrix(patient.clinical.trans[4,])), 
                         survival_months = as.numeric(as.matrix(patient.clinical.trans[5,])), 
                         cancer_type = as.character(as.matrix(patient.clinical.trans[6,])),
                         aneuploidy_score = as.numeric(as.matrix(patient.clinical.trans[7,])), 
                         sample_type = as.character(as.matrix(patient.clinical.trans[8,])) )
                         # col = list(days = col_fun2, stage = col_fun3, NRAS = col_fun4))

breaksList = seq(min(rnaseq.saav), max(rnaseq.saav), by=2)
col = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList))
pdf("Data/HM_SAAV_RNAseq.pdf", height = 20, width = 30)
hp = Heatmap(as.matrix(rnaseq.saav), col = col, column_title = paste("SAAV-related RNAseq log2-transformed Median Expression"), top_annotation = anno,
             cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = FALSE, show_column_names = FALSE, name = "value", heatmap_legend_param = list(direction = "vertical"))
draw(hp, heatmap_legend_side = "bottom", 
     annotation_legend_side = "right", merge_legend = TRUE,)
dev.off()



non.rnaseq.saav = rnaseq[-(rnaseq$Hugo_Symbol %in% unique(SAAV$SYMBOL)), 3:445]
MW = data.frame(RNAseq=c(as.matrix(rnaseq.saav)), group='SAAV')
MWn = data.frame(RNAseq=c(as.matrix(non.rnaseq.saav)), group='nonSAAV')
MW = rbind(MW, MWn)

# Mann-Whitney test
library(ggplot2)
library(ggpubr)
library(gridExtra)
pp = ggboxplot(MW, x = "group", y = "RNAseq",
               color = "black", fill = "group", palette = "jco")+ 
  stat_compare_means(method.args = list(alternative = "two.sided"), aes(group = group), label = "p.signif", label.y = 22) + 
  stat_compare_means(method.args = list(alternative = "two.sided"), aes(group = group), label = "p.format", label.y = 23)

pdf(file="Data/SAAV_nonSAAV.pdf", width=8,height=8)
grid.arrange(pp,nrow=1)
dev.off()

# 
# library(ggplot2)
# rnaseq.saav$mean = rowMeans(rnaseq.saav)
# rnaseq.saav$gene = row.names(rnaseq.saav)
# 
# ggplot(rnaseq.saav, aes(x=gene, y=mean)) +
#   geom_bar(stat="identity", fill='#A4A4A4', color="black")+
#   theme_classic()

# Heatmap for all RNA
RNAseq <- read.csv("Data/log2_RNAseq.csv")
RNAseq$CORUM = 0
library(MoonFinder)
data("corumModuleList")
corum_gene=c()
for (corum in corumModuleList){
  corum_gene=c(corum_gene, corum)
}

for (row in 1:nrow(RNAseq)){
  if (RNAseq[row, "Hugo_Symbol"] %in% corum_gene){
    RNAseq[row, "CORUM"] = 1
  }
}

clinical <- read.delim("~/Documents/TCGA_melanoma/Data/skcm_tcga_pan_can_atlas_2018/data_clinical_sample.txt", comment.char="#")
clinical$SAMPLE_ID=gsub('-','.', clinical$SAMPLE_ID)
clinical = clinical[ , c('PATIENT_ID', 'SAMPLE_ID', 'CANCER_TYPE_DETAILED', 'ANEUPLOIDY_SCORE', 'SAMPLE_TYPE')]
patient = read.delim("~/Documents/TCGA_melanoma/Data/skcm_tcga_pan_can_atlas_2018/data_clinical_patient.txt", comment.char="#")
patient = patient[, c('PATIENT_ID', 'AGE', 'SEX', 'AJCC_PATHOLOGIC_TUMOR_STAGE', 'WEIGHT', 'DSS_MONTHS')]

patient.clinical = merge(patient, clinical, by="PATIENT_ID")
row.names(patient.clinical) = patient.clinical$SAMPLE_ID
patient.clinical = patient.clinical[,-c(1, 7)]
patient.clinical.trans = data.frame(t(patient.clinical))

row.names(RNAseq) = RNAseq$Entrez_Gene_Id
RNAseq = RNAseq[, -c(1,2)]

patient.clinical.trans = patient.clinical.trans[, colnames(patient.clinical.trans) %in% colnames(RNAseq)]
patient.clinical.trans$CORUM = NA
combined = rbind(RNAseq, patient.clinical.trans)


combined = rbind(combined[20524:20531,], combined[1:20523,])
write.csv(combined, 'Data/clinical_RNAseq.csv', row.names = TRUE)

gn = rowAnnotation(CORUM = combined$CORUM)


anno = HeatmapAnnotation(age = as.numeric(as.matrix(patient.clinical.trans[1,1:443])), 
                         gender = as.character(as.matrix(patient.clinical.trans[2,1:443])), 
                         AJCC_stage = as.character(as.matrix(patient.clinical.trans[3,1:443])),
                         weight = as.numeric(as.matrix(patient.clinical.trans[4,1:443])), 
                         survival_months = as.numeric(as.matrix(patient.clinical.trans[5,1:443])), 
                         cancer_type = as.character(as.matrix(patient.clinical.trans[6,1:443])),
                         aneuploidy_score = as.numeric(as.matrix(patient.clinical.trans[7,1:443])), 
                         sample_type = as.character(as.matrix(patient.clinical.trans[8,1:443])) )
# col = list(days = col_fun2, stage = col_fun3, NRAS = col_fun4))

breaksList = seq(min(RNAseq), max(RNAseq), by=2)
col = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList))
pdf("Data/HM_RNAseq.pdf", height = 30, width = 30)
hp = Heatmap(as.matrix(RNAseq[,1:443]), col = col, column_title = paste("RNAseq log2-transformed Median Expression"), top_annotation = anno, right_annotation=gn,
             cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = FALSE, show_column_names = FALSE, name = "value", heatmap_legend_param = list(direction = "vertical"))
draw(hp, heatmap_legend_side = "bottom", 
     annotation_legend_side = "right", merge_legend = TRUE,)
dev.off()

