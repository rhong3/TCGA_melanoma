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
pdf("Results/Mitochondria/MCM_HM_RNAseq.pdf", height = 8, width = 30)
hp = Heatmap(as.matrix(MCM.clinical[11:16,]), col = col, column_title = paste("MCM RNAseq log2-transformed Median Expression"), top_annotation = anno,
             cluster_rows = FALSE, cluster_columns = TRUE, show_row_names = TRUE, show_column_names = FALSE, name = "value", heatmap_legend_param = list(direction = "vertical"))
draw(hp, heatmap_legend_side = "bottom", 
     annotation_legend_side = "right", merge_legend = TRUE,)
dev.off()

