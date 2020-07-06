# Mitochondria validation
RNAseq <- read.csv("Data/log2_RNAseq.csv")
clinical <- read.delim("~/Documents/TCGA_melanoma/Data/skcm_tcga_pan_can_atlas_2018/data_clinical_sample.txt", comment.char="#")
clinical$SAMPLE_ID=gsub('-','.', clinical$SAMPLE_ID)
clinical = clinical[ , c('PATIENT_ID', 'SAMPLE_ID', 'TUMOR_TISSUE_SITE', 'ANEUPLOIDY_SCORE', 'SAMPLE_TYPE')]
patient = read.delim("~/Documents/TCGA_melanoma/Data/skcm_tcga_pan_can_atlas_2018/data_clinical_patient.txt", comment.char="#")
patient = patient[, c('PATIENT_ID', 'AGE', 'SEX', 'AJCC_PATHOLOGIC_TUMOR_STAGE', 'WEIGHT', 'DSS_MONTHS', 
                      'HISTORY_NEOADJUVANT_TRTYN', 'NEW_TUMOR_EVENT_AFTER_INITIAL_TREATMENT')]
patient.clinical = merge(patient, clinical, by="PATIENT_ID")
row.names(patient.clinical) = patient.clinical$SAMPLE_ID
patient.clinical$TUMOR_TISSUE_SITE = gsub('Trunk', '', patient.clinical$TUMOR_TISSUE_SITE)
patient.clinical$TUMOR_TISSUE_SITE = gsub('Other', '', patient.clinical$TUMOR_TISSUE_SITE)
patient.clinical$TUMOR_TISSUE_SITE = gsub('Extremities', '', patient.clinical$TUMOR_TISSUE_SITE)
patient.clinical$TUMOR_TISSUE_SITE = gsub('Head and Neck', '', patient.clinical$TUMOR_TISSUE_SITE)
patient.clinical$TUMOR_TISSUE_SITE = gsub("NA", '', patient.clinical$TUMOR_TISSUE_SITE)
patient.clinical$TUMOR_TISSUE_SITE = gsub("\\|", '', patient.clinical$TUMOR_TISSUE_SITE)
patient.clinical[patient.clinical==""]<-NA

patient.clinical = patient.clinical[,-c(1, 9)]
patient.clinical.trans = data.frame(t(patient.clinical))

row.names(RNAseq) = RNAseq$Entrez_Gene_Id
RNAseq = RNAseq[, -c(1,2)]

patient.clinical.trans = patient.clinical.trans[, colnames(patient.clinical.trans) %in% colnames(RNAseq)]
combined = rbind(RNAseq, patient.clinical.trans)
combined = rbind(combined[20532:20541,], combined[1:20531,])
write.csv(combined, 'Data/big_clinical_RNAseq.csv', row.names = TRUE)


# MCM clustering
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
RNAseq <- read.csv("Data/log2_RNAseq.csv")
patient.clinical = read.csv('Data/big_clinical_RNAseq.csv')
row.names(patient.clinical) = patient.clinical$X
patient.clinical = patient.clinical[1:10, -1]
MCM = RNAseq[RNAseq$Hugo_Symbol %in% c('MCM2', 'MCM3', 'MCM4', 'MCM5', 'MCM6', 'MCM7'), ]
row.names(MCM) = MCM$Hugo_Symbol
MCM = MCM[, -c(1,2)]
MCM.numeric <- as.matrix(sapply(MCM, as.numeric))
row.names(MCM.numeric) = row.names(MCM)
MCM.clinical = rbind(patient.clinical, MCM.numeric)
MCM.numeric = MCM.clinical[11:16, ]
MCM.numeric <- as.matrix(sapply(MCM, as.numeric))
row.names(MCM.numeric) = row.names(MCM)

anno = HeatmapAnnotation(age = as.numeric(as.matrix(MCM.clinical[1,1:443])), 
                         gender = as.character(as.matrix(MCM.clinical[2,1:443])), 
                         tumor_stage = as.character(as.matrix(MCM.clinical[3,1:443])),
                         weight = as.numeric(as.matrix(MCM.clinical[4,1:443])), 
                         survival_months = as.numeric(as.matrix(MCM.clinical[5,1:443])), 
                         treatment = as.character(as.matrix(MCM.clinical[6,1:443])),
                         new_tumor_after_treatment = as.character(as.matrix(MCM.clinical[7,1:443])),
                         sample_origin = as.character(as.matrix(MCM.clinical[8,1:443])), 
                         aneuploidy_score = as.numeric(as.matrix(MCM.clinical[9,1:443])), 
                         sample_type = as.character(as.matrix(MCM.clinical[10,1:443])) )

breaksList = seq(min(MCM.numeric), max(MCM.numeric), by=2)
col = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList))

hp = Heatmap(MCM.numeric, col = col, column_title = paste("MCM RNAseq log2-transformed Median Expression"), top_annotation = anno, column_km = 4, column_km_repeats = 100,
             cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = TRUE, show_column_names = FALSE, name = "value", heatmap_legend_param = list(direction = "vertical"))

MCM.clinical.cluster = data.frame(matrix(nrow = 1, ncol = ncol(MCM.clinical)))
colnames(MCM.clinical.cluster) = colnames(MCM.clinical)
row.names(MCM.clinical.cluster) = c('Proliferation.status')
for (i in column_order(hp)$`1`){
  MCM.clinical.cluster[1, i] = 'Low'
}
for (i in column_order(hp)$`2`){
  MCM.clinical.cluster[1, i] = 'Medium_Low'
}
for (i in column_order(hp)$`3`){
  MCM.clinical.cluster[1, i] = 'Medium_High'
}
for (i in column_order(hp)$`4`){
  MCM.clinical.cluster[1, i] = 'High'
}
MCM.clinical = rbind(MCM.clinical.cluster, MCM.clinical)
write.csv(MCM.clinical, 'Data/MCM_RNAseq.csv')

MCM.numeric = MCM.clinical[12:17, ]
MCM.numeric <- as.matrix(sapply(MCM, as.numeric))
row.names(MCM.numeric) = row.names(MCM)

anno = HeatmapAnnotation(age = as.numeric(as.matrix(MCM.clinical[2,1:443])), 
                         gender = as.character(as.matrix(MCM.clinical[3,1:443])), 
                         tumor_stage = as.character(as.matrix(MCM.clinical[4,1:443])),
                         weight = as.numeric(as.matrix(MCM.clinical[5,1:443])), 
                         survival_months = as.numeric(as.matrix(MCM.clinical[6,1:443])), 
                         treatment = as.character(as.matrix(MCM.clinical[7,1:443])),
                         new_tumor_after_treatment = as.character(as.matrix(MCM.clinical[8,1:443])),
                         sample_origin = as.character(as.matrix(MCM.clinical[9,1:443])), 
                         aneuploidy_score = as.numeric(as.matrix(MCM.clinical[10,1:443])), 
                         sample_type = as.character(as.matrix(MCM.clinical[11,1:443])))
banno = HeatmapAnnotation(proliferation_status = as.character(as.matrix(MCM.clinical[1,1:443])))

breaksList = seq(min(MCM.numeric), max(MCM.numeric), by=2)
col = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList))


pdf("Results/Mitochondria/MCM_HM_RNAseq.pdf", height = 8, width = 30)
hp = Heatmap(MCM.numeric, col = col, column_title = paste("MCM RNAseq log2-transformed Median Expression"), top_annotation = anno, column_km = 4, column_km_repeats = 100, bottom_annotation = banno,
             cluster_rows = FALSE, cluster_columns = TRUE, show_row_names = TRUE, show_column_names = FALSE, name = "value", heatmap_legend_param = list(direction = "vertical"))
draw(hp, heatmap_legend_side = "bottom", 
     annotation_legend_side = "right", merge_legend = TRUE,)
dev.off()


# Heatmap Clustering
library(readxl)
library(biomaRt)
library(org.Hs.eg.db)
library(stringr)
library(readr)
library(plyr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(robustHD)
RNAseq <- read.csv("Data/log2_RNAseq.csv")
patient.clinical = read.csv('Data/MCM_RNAseq.csv')
patient.clinical = patient.clinical[1:11,]
row.names(patient.clinical) = patient.clinical$X
patient.clinical = patient.clinical[1:11, -1]

Prospective <- read_excel("Data/Melanoma proteome/MitochondrialProteinListProspectivePostmortem.xlsx", sheet = "Prospective")
Postmortem <- read_excel("Data/Melanoma proteome/MitochondrialProteinListProspectivePostmortem.xlsx", sheet = "Postmortem")
Prospective = Prospective[-c(1:11), -c(78:84, 86:91)]
Postmortem = Postmortem[-c(1:13), -c(75:81, 83:87)]
Prospective[,1:77] <- lapply(Prospective[,1:77], function(x) as.numeric(gsub("NaN", NA, x)))
Postmortem[,1:74] <- lapply(Postmortem[,1:74], function(x) as.numeric(gsub("NaN", NA, x)))
temp = strsplit(Prospective$PG.ProteinAccession, ';')
for (i in 1:length(temp)){Prospective[i,78] = temp[[i]][1]}
temp = strsplit(Postmortem$`0PG.ProteinAccessions`, ';')
for (i in 1:length(temp)){Postmortem[i,75] = temp[[i]][1]}

annotation <- select(org.Hs.eg.db, keys=unique(Prospective$PG.ProteinAccession), columns=c('UNIPROT', 'SYMBOL'), keytype="UNIPROT")
colnames(Prospective) = gsub('PG.ProteinAccession', "UNIPROT", colnames(Prospective))
Prospective = merge(Prospective, annotation, by="UNIPROT")

annotation <- select(org.Hs.eg.db, keys=unique(Postmortem$`0PG.ProteinAccessions`), columns=c('UNIPROT', 'SYMBOL'), keytype="UNIPROT")
colnames(Postmortem) = gsub('0PG.ProteinAccessions', "UNIPROT", colnames(Postmortem))
Postmortem = merge(Postmortem, annotation, by="UNIPROT")

RNA_prospective = RNAseq[RNAseq$Hugo_Symbol %in% unique(Prospective$SYMBOL), ]
RNA_Postmortem = RNAseq[RNAseq$Hugo_Symbol %in% unique(Postmortem$SYMBOL), ]
row.names(RNA_prospective) = RNA_prospective$Hugo_Symbol
row.names(RNA_Postmortem) = RNA_Postmortem$Hugo_Symbol
RNA_prospective= RNA_prospective[,-c(1,2)]
RNA_Postmortem = RNA_Postmortem[, -c(1,2)]
RNA_prospective_clinical = rbind(as.matrix(patient.clinical), as.matrix(RNA_prospective))
write.csv(RNA_prospective_clinical, 'Data/RNA_prospective_clinical.csv')
RNA_Postmortem_clinical = rbind(as.matrix(patient.clinical), as.matrix(RNA_Postmortem))
write.csv(RNA_prospective_clinical, 'Data/RNA_Postmortem_clinical.csv')

RNA_prospective.numeric = RNA_prospective_clinical[12:nrow(RNA_prospective_clinical), ]
RNA_prospective.numeric <- as.matrix(sapply(RNA_prospective, as.numeric))
row.names(RNA_prospective.numeric) = row.names(RNA_prospective_clinical)[12:nrow(RNA_prospective_clinical)]

RNA_Postmortem.numeric = RNA_Postmortem_clinical[12:nrow(RNA_Postmortem_clinical), ]
RNA_Postmortem.numeric <- as.matrix(sapply(RNA_Postmortem, as.numeric))
row.names(RNA_Postmortem.numeric) = row.names(RNA_Postmortem_clinical)[12:nrow(RNA_Postmortem_clinical)]

anno = HeatmapAnnotation(proliferation_status = as.character(as.matrix(patient.clinical[1,1:443])),
                         age = as.numeric(as.matrix(patient.clinical[2,1:443])), 
                         gender = as.character(as.matrix(patient.clinical[3,1:443])), 
                         tumor_stage = as.character(as.matrix(patient.clinical[4,1:443])),
                         weight = as.numeric(as.matrix(patient.clinical[5,1:443])), 
                         survival_months = as.numeric(as.matrix(patient.clinical[6,1:443])), 
                         treatment = as.character(as.matrix(patient.clinical[7,1:443])),
                         new_tumor_after_treatment = as.character(as.matrix(patient.clinical[8,1:443])),
                         sample_origin = as.character(as.matrix(patient.clinical[9,1:443])), 
                         aneuploidy_score = as.numeric(as.matrix(patient.clinical[10,1:443])), 
                         sample_type = as.character(as.matrix(patient.clinical[11,1:443])))

breaksList = seq(min(RNA_prospective.numeric), max(RNA_prospective.numeric), by=2)
col = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList))

pdf("Results/Mitochondria/Prospective_HM_RNAseq.pdf", height = 30, width = 30)
hp = Heatmap(RNA_prospective.numeric, col = col, column_title = paste("Prospective RNAseq log2-transformed Median Expression"), top_annotation = anno, column_km = 3, column_km_repeats = 100, row_km = 3, row_km_repeats = 100, 
             cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = FALSE, show_column_names = FALSE, name = "value", heatmap_legend_param = list(direction = "vertical"))
draw(hp, heatmap_legend_side = "bottom", 
     annotation_legend_side = "right", merge_legend = TRUE,)
dev.off()


breaksList = seq(min(RNA_Postmortem.numeric), max(RNA_Postmortem.numeric), by=2)
col = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList))

pdf("Results/Mitochondria/Postmortem_HM_RNAseq.pdf", height = 30, width = 30)
hp = Heatmap(RNA_Postmortem.numeric, col = col, column_title = paste("Postmortem RNAseq log2-transformed Median Expression"), top_annotation = anno, column_km = 4, column_km_repeats = 100, row_km = 3, row_km_repeats = 100, 
             cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = FALSE, show_column_names = FALSE, name = "value", heatmap_legend_param = list(direction = "vertical"))
draw(hp, heatmap_legend_side = "bottom", 
     annotation_legend_side = "right", merge_legend = TRUE,)
dev.off()


# Quantify ADP/ATP translocases/Mitochondrial fission/Mitochondrial fusion/Mitochondrial surtuins
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(dplyr)

RNAseq <- read.csv("Data/log2_RNAseq.csv")
patient.clinical = read.csv('Data/MCM_RNAseq.csv')
patient.clinical = patient.clinical[1:11,]
row.names(patient.clinical) = patient.clinical$X
patient.clinical = patient.clinical[1:11, -1]

ADPATP = c('SLC25A4', 'SLC25A5', 'SLC25A6')
fission = c('DNM1L', 'MFF', 'FIS1', 'MID51')
fusion = c('MFN1', 'MFN2', 'OPA1', 'OPA3')
sirtuins = c('SIRT3', 'SIRT5')

RNA_ADPATP = RNAseq[RNAseq$Hugo_Symbol %in% ADPATP, ]
row.names(RNA_ADPATP) = RNA_ADPATP$Hugo_Symbol
RNA_fission = RNAseq[RNAseq$Hugo_Symbol %in% fission, ]
row.names(RNA_fission) = RNA_fission$Hugo_Symbol
RNA_fusion = RNAseq[RNAseq$Hugo_Symbol %in% fusion, ]
row.names(RNA_fusion) = RNA_fusion$Hugo_Symbol
RNA_sirtuins = RNAseq[RNAseq$Hugo_Symbol %in% sirtuins, ]
row.names(RNA_sirtuins) = RNA_sirtuins$Hugo_Symbol


RNA_ADPATP_clinical = t(rbind(as.matrix(patient.clinical), as.matrix(RNA_ADPATP[, 3:445])))
RNA_ADPATP_plot = rbind(RNA_ADPATP_clinical[,c(9,12)], RNA_ADPATP_clinical[,c(9,13)], RNA_ADPATP_clinical[,c(9,14)])
colnames(RNA_ADPATP_plot)[2] = 'Expression'
RNA_ADPATP_plot = as.data.frame(RNA_ADPATP_plot)
RNA_ADPATP_plot$Gene = c(rep('SLC25A4', 443), rep('SLC25A5', 443), rep('SLC25A6', 443))
RNA_ADPATP_plot[, 2] <- sapply(RNA_ADPATP_plot[, 2], as.character)
RNA_ADPATP_plot[, 2] <- sapply(RNA_ADPATP_plot[, 2], as.numeric)
RNA_ADPATP_plot = na.omit(RNA_ADPATP_plot)

RNA_fission_clinical = t(rbind(as.matrix(patient.clinical), as.matrix(RNA_fission[, 3:445])))
RNA_fission_plot = rbind(RNA_fission_clinical[,c(9,12)], RNA_fission_clinical[,c(9,13)], RNA_fission_clinical[,c(9,14)])
colnames(RNA_fission_plot)[2] = 'Expression'
RNA_fission_plot = as.data.frame(RNA_fission_plot)
RNA_fission_plot$Gene = c(rep('DNM1L', 443), rep('FIS1', 443), rep('MFF', 443))
RNA_fission_plot[, 2] <- sapply(RNA_fission_plot[, 2], as.character)
RNA_fission_plot[, 2] <- sapply(RNA_fission_plot[, 2], as.numeric)
RNA_fission_plot = na.omit(RNA_fission_plot)

RNA_fusion_clinical = t(rbind(as.matrix(patient.clinical), as.matrix(RNA_fusion[, 3:445])))
RNA_fusion_plot = rbind(RNA_fusion_clinical[,c(9,12)], RNA_fusion_clinical[,c(9,13)], RNA_fusion_clinical[,c(9,14)], RNA_fusion_clinical[,c(9,15)])
colnames(RNA_fusion_plot)[2] = 'Expression'
RNA_fusion_plot = as.data.frame(RNA_fusion_plot)
RNA_fusion_plot$Gene = c(rep('MFN1', 443), rep('MFN2', 443), rep('OPA1', 443), rep('OPA3', 443))
RNA_fusion_plot[, 2] <- sapply(RNA_fusion_plot[, 2], as.character)
RNA_fusion_plot[, 2] <- sapply(RNA_fusion_plot[, 2], as.numeric)
RNA_fusion_plot = na.omit(RNA_fusion_plot)

RNA_sirtuins_clinical = t(rbind(as.matrix(patient.clinical), as.matrix(RNA_sirtuins[, 3:445])))
RNA_sirtuins_plot = rbind(RNA_sirtuins_clinical[,c(9,12)], RNA_sirtuins_clinical[,c(9,13)])
colnames(RNA_sirtuins_plot)[2] = 'Expression'
RNA_sirtuins_plot = as.data.frame(RNA_sirtuins_plot)
RNA_sirtuins_plot$Gene = c(rep('SIRT3',  443), rep('SIRT5', 443))
RNA_sirtuins_plot[, 2] <- sapply(RNA_sirtuins_plot[, 2], as.character)
RNA_sirtuins_plot[, 2] <- sapply(RNA_sirtuins_plot[, 2], as.numeric)
RNA_sirtuins_plot = na.omit(RNA_sirtuins_plot)


RNA_ADPATP_summary = RNA_ADPATP_plot %>%
  group_by(TUMOR_TISSUE_SITE, Gene) %>%
  summarise(expression = mean(Expression), sd = sd(Expression))

pa<- ggplot(RNA_ADPATP_summary, aes(x=Gene, y=expression, fill=TUMOR_TISSUE_SITE)) +
  xlab("") + ggtitle('ADP/ATP translocases') + 
  geom_bar(stat="identity", color="black",
           position=position_dodge())+ coord_cartesian(ylim = c(0,16)) + 
  scale_fill_manual(values=c("#e41a1c", "#377eb8", "#984ea3", "#ff7f00")) +
  geom_errorbar(aes(ymin=expression-sd, ymax=expression+sd), width=.2,
                position=position_dodge(.9)) +
  theme_bw()+ 
  theme(axis.text.x = element_text(size = 15, angle = 0, hjust = 1), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), legend.position = "none")


RNA_fission_summary = RNA_fission_plot %>%
  group_by(TUMOR_TISSUE_SITE, Gene) %>%
  summarise(expression = mean(Expression), sd = sd(Expression))

pb<- ggplot(RNA_fission_summary, aes(x=Gene, y=expression, fill=TUMOR_TISSUE_SITE)) +
  xlab("") + ggtitle('Mitochondrial fission') + 
  geom_bar(stat="identity", color="black",
           position=position_dodge())+ coord_cartesian(ylim = c(0,16)) + 
  scale_fill_manual(values=c("#e41a1c", "#377eb8", "#984ea3", "#ff7f00")) +
  geom_errorbar(aes(ymin=expression-sd, ymax=expression+sd), width=.2,
                position=position_dodge(.9)) +
  theme_bw()+ 
  theme(axis.text.x = element_text(size = 15, angle = 0, hjust = 1), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), legend.position = "none")


RNA_fusion_summary = RNA_fusion_plot %>%
  group_by(TUMOR_TISSUE_SITE, Gene) %>%
  summarise(expression = mean(Expression), sd = sd(Expression))

pc<- ggplot(RNA_fusion_summary, aes(x=Gene, y=expression, fill=TUMOR_TISSUE_SITE)) +
  xlab("") + ggtitle('Mitochondrial fusion') + 
  geom_bar(stat="identity", color="black",
           position=position_dodge())+ coord_cartesian(ylim = c(0,16)) + 
  scale_fill_manual(values=c("#e41a1c", "#377eb8", "#984ea3", "#ff7f00")) +
  geom_errorbar(aes(ymin=expression-sd, ymax=expression+sd), width=.2,
                position=position_dodge(.9)) +
  theme_bw()+ 
  theme(axis.text.x = element_text(size = 15, angle = 0, hjust = 1), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), legend.position = "none")


RNA_sirtuins_summary = RNA_sirtuins_plot %>%
  group_by(TUMOR_TISSUE_SITE, Gene) %>%
  summarise(expression = mean(Expression), sd = sd(Expression))

pd<- ggplot(RNA_sirtuins_summary, aes(x=Gene, y=expression, fill=TUMOR_TISSUE_SITE)) +
  xlab("") + ggtitle('Mitochondrial sirtuins') + 
  geom_bar(stat="identity", color="black",
           position=position_dodge())+ coord_cartesian(ylim = c(0,16)) + 
  scale_fill_manual(values=c("#e41a1c", "#377eb8", "#984ea3", "#ff7f00")) +
  geom_errorbar(aes(ymin=expression-sd, ymax=expression+sd), width=.2,
                position=position_dodge(.9)) +
  theme_bw()+ 
  theme(axis.text.x = element_text(size = 15, angle = 0, hjust = 1), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), legend.position = "right", legend.text = element_text(size = 6), legend.title = element_blank())


pdf(file='Results/Mitochondria/quantified_features.pdf',
    width=10,height=10)
grid.arrange(pa, pb, pc, pd, nrow=2, ncol=2)
dev.off()










