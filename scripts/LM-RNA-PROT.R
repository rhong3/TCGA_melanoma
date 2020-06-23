# LM RNA-PROT
library(robustHD)
# Combine proteomics data
Erika <- read_excel("Data/Melanoma proteome/Tumor samples/DDA/Erika Pilot raw and normalized-no normal tissue.xlsx", 
                    sheet = "normalized")
Erika = Erika[, -2]

Lotta <- read_excel("Data/Melanoma proteome/Tumor samples/DDA/Lotta Fractionation SCX raw and normalized.xlsx", 
                    sheet = "normalized")
Lotta = Lotta[,-2]
proteomics = merge(Erika, Lotta, by=c('Accession', 'Gene name'),  all=TRUE, suffixes=c('Erika', 'Lotta'))

Masha <- read_excel("Data/Melanoma proteome/Tumor samples/DDA/Masha-Eupa-fractionation_raw and normalized.xlsx", 
                    sheet = "normalized")
colnames(Masha) = gsub('Gene Symbol', 'Gene name', colnames(Masha))
Masha = Masha[,-2]
proteomics = merge(proteomics, Masha, by=c('Accession', 'Gene name'),  all=TRUE, suffixes=c('x', 'Masha'))

MEDAFASP <- read_excel("Data/Melanoma proteome/Tumor samples/DDA/MEDAFASP-MM_pilot_raw and normalized.xlsx", 
                       sheet = "normalized")
colnames(MEDAFASP) = gsub('Gene Symbol', 'Gene name', colnames(MEDAFASP))
MEDAFASP = MEDAFASP[,-2]
proteomics = merge(proteomics, MEDAFASP, by=c('Accession', 'Gene name'),  all=TRUE, suffixes=c('x', 'MEDAFASP'))

Pilot <- read_excel("Data/Melanoma proteome/Tumor samples/DDA/Pilot Jeo DDA raw and normalized.xlsx", 
                    sheet = "normalized")
colnames(Pilot) = gsub('Gene Symbol', 'Gene name', colnames(Pilot))
Pilot = Pilot[, -22]
proteomics = merge(proteomics, Pilot, by=c('Accession', 'Gene name'),  all=TRUE, suffixes=c('x', 'Pilot'))

Primero <- read_excel("Data/Melanoma proteome/Tumor samples/DDA/Primero_Proteins_raw_normalized.xlsx", 
                      sheet = "normalized")[,2:115]
Primero = Primero[, -2]
proteomics = merge(proteomics, Primero, by=c('Accession', 'Gene name'),  all=TRUE, suffixes=c('x', 'Primero'))

Segundo <- read_csv("Data/Melanoma proteome/Tumor samples/DDA/Segundo_proteomics.csv")
Segundo = Segundo[,-112]
proteomics = merge(proteomics, Segundo, by=c('Accession', 'Gene name'),  all=TRUE, suffixes=c('x', 'Segundo'))

Cuarto <- read_excel("Data/Melanoma proteome/Tumor samples/DIA/Cuarto(IV) raw and normalized-no normal tissue.xlsx", 
                     sheet = "Sheet1")
colnames(Cuarto) = gsub('PG.ProteinAccessions', 'Accession', colnames(Cuarto))
Cuarto = Cuarto[,-c(2, 70)]
proteomics = merge(proteomics, Cuarto, by=c('Accession', 'Gene name'),  all=TRUE, suffixes=c('x', 'Cuarto'))

Tercero <- read_excel("Data/Melanoma proteome/Tumor samples/DIA/Tercero (III)_raw and normalized.xlsx", 
                      sheet = "normalized")
Tercero[,3:76] = mutate_all(Tercero[,3:76], function(x) as.numeric(as.character(x)))
Tercero[,3:76] = (Tercero[,3:76] - mean(as.matrix(Tercero[,3:76]), na.rm = TRUE))/sd(as.matrix(Tercero[,3:76]), na.rm = TRUE)
proteomics = merge(proteomics, Tercero, by=c('Accession', 'Gene name'),  all=TRUE, suffixes=c('x', 'Tercero'))

cells <- read_excel("Data/Melanoma proteome/Cells/DDA/cells Pilot Jeo DDA raw and normalized.xlsx", 
                    sheet = "normalized")
colnames(cells) = gsub('Gene Symbol', 'Gene name', colnames(cells))
cells = cells[, -2]
proteomics = merge(proteomics, cells, by=c('Accession', 'Gene name'),  all=TRUE, suffixes=c('x', 'cells'))

RPMI <- read_excel("Data/Melanoma proteome/Cells/DDA/pilot_RPMI cell_raw and normalized.xlsx", 
                   sheet = "normalized")
colnames(RPMI) = gsub('Gene Symbol', 'Gene name', colnames(RPMI))
RPMI = RPMI[, -2]
proteomics = merge(proteomics, RPMI, by=c('Accession', 'Gene name'),  all=TRUE, suffixes=c('x', 'RPMI'))

CellsAcetyl <- read_excel("Data/Melanoma proteome/Cells/DIA/CellsAcetyl_DIA_ArgC_raw and normalized.xlsx", 
                          sheet = "normalized")
colnames(CellsAcetyl) = gsub('PG.ProteinAccessions', 'Accession', colnames(CellsAcetyl))
CellsAcetyl = CellsAcetyl[,-2]
proteomics = merge(proteomics, CellsAcetyl, by=c('Accession', 'Gene name'),  all=TRUE, suffixes=c('x', 'CellsAcetyl'))

CellsTrypKR <- read_excel("Data/Melanoma proteome/Cells/DIA/CellsTrypKR_DIA_raw and normalized.xlsx", 
                          sheet = "normalized")
colnames(CellsTrypKR) = gsub('PG.ProteinAccessions', 'Accession', colnames(CellsTrypKR))
colnames(CellsTrypKR) = gsub('PG.Genes', 'Gene name', colnames(CellsTrypKR))
CellsTrypKR = CellsTrypKR[,-3]
proteomics = merge(proteomics, CellsTrypKR, by=c('Accession', 'Gene name'),  all=TRUE, suffixes=c('x', 'CellsTrypKR'))

write_csv(proteomics, "Data/proteomics_all.csv")

#-------------------------------------------
# CORUM protein complex
library(MoonFinder)
data("corumModuleList")
corum_gene=c()
for (corum in corumModuleList){
  corum_gene=c(corum_gene, corum)
}

RNAseq <- read_csv("Data/log2_RNAseq.csv")
proteomics <- read_csv("Data/proteomics_all.csv")

# Venn plot of RNA, PROT
library(VennDiagram)
PROT_unq = unique(na.omit(proteomics$`Gene name`))
RNA_unq = unique(na.omit(RNAseq$Hugo_Symbol))
venn.diagram(
  x = list(PROT_unq, RNA_unq, corum_gene),
  category.names = c("Genes in proteomics (12770)" , "Genes in transcriptomics (20501)", "Genes in complex (1670)"),
  filename = 'Results/PROT_RNA_CORUM_venn.png',
  output = TRUE ,
  imagetype="png" ,
  height = 700 , 
  width = 700 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#440154ff", '#21908dff', '#e6550d'),
  fill = c(alpha("#440154ff",0.5), alpha('#21908dff',0.5), alpha('#e6550d',0.5)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#440154ff", '#21908dff', '#e6550d')
)

# Take mean
library(matrixStats)
library(dplyr)
RNAseq$RNA_Mean = rowMeans(RNAseq[,3:445], na.rm = TRUE)
RNAseq$RNA_Median = rowMedians(as.matrix(RNAseq[,3:445]), na.rm = TRUE)
RNAseq$RNA_STD = rowSds(as.matrix(RNAseq[,3:445]), na.rm = TRUE)
RNAseq$RNA_absVMR = abs(RNAseq$RNA_STD*RNAseq$RNA_STD/RNAseq$RNA_Mean)
RNAseq = RNAseq[,-c(2:445)]
RNAseq = na.omit(RNAseq)
colnames(RNAseq) = gsub('Hugo_Symbol', 'Gene name', colnames(RNAseq))

proteomics$PROT_Mean = rowMeans(proteomics[,3:448], na.rm = TRUE)
proteomics$PROT_Median = rowMedians(as.matrix(proteomics[,3:448]), na.rm = TRUE)
proteomics$PROT_STD = rowSds(as.matrix(proteomics[,3:448]), na.rm = TRUE)
proteomics$PROT_absVMR = abs(proteomics$PROT_STD*proteomics$PROT_STD/proteomics$PROT_Mean)
proteomics = proteomics[,-c(1,3:448)]
proteomics = na.omit(proteomics)

RNA_PROT = merge(RNAseq, proteomics, by = 'Gene name')

RNA_PROT_sum = RNA_PROT  %>%
  group_by(`Gene name`) %>%
  summarize(
    RNA_Mean = mean(RNA_Mean)
    , RNA_Median = mean(RNA_Median)
    , RNA_STD = mean(RNA_STD)
    , RNA_absVMR = mean(RNA_absVMR)
    , PROT_Mean = mean(PROT_Mean)
    , PROT_Median = mean(PROT_Median)
    , PROT_STD = mean(PROT_STD)
    , PROT_absVMR = mean(PROT_absVMR)
  )

write.csv(RNA_PROT_sum, "Data/RNA_PROT_sum.csv", row.names=FALSE)

# Correlation plots
library(ggplot2)
library(ggpubr)
RNA_PROT_sum = read_csv("Data/RNA_PROT_sum.csv")
RNA_PROT_sum$Mean_outlier = "NA"
RNA_PROT_sum$Median_outlier = "NA"
RNA_PROT_sum$CORUM = "No"
for (row in 1:nrow(RNA_PROT_sum)){
  if (RNA_PROT_sum[row, "PROT_Mean"] < -10 | RNA_PROT_sum[row, "PROT_Mean"] > 10){
    RNA_PROT_sum[row, "Mean_outlier"] = RNA_PROT_sum[row, "Gene name"]
  }
  if (RNA_PROT_sum[row, "PROT_Median"] < -10 | RNA_PROT_sum[row, "PROT_Mean"] > 10){
    RNA_PROT_sum[row, "Median_outlier"] = RNA_PROT_sum[row, "Gene name"]
  }
  if (RNA_PROT_sum[row, "Gene name"] %in% corum_gene){
    RNA_PROT_sum[row, "CORUM"] = "Yes"
  }
}
RNA_PROT_sum[RNA_PROT_sum=="NA"]=NA
mean_cor = cor(RNA_PROT_sum$RNA_Mean, RNA_PROT_sum$PROT_Mean, method='pearson')
median_cor = cor(RNA_PROT_sum$RNA_Median, RNA_PROT_sum$PROT_Median, method='pearson')

pdf(file='Results/Mean_cor_all.pdf',
    width=15,height=10,useDingbats=FALSE)
ggscatter(RNA_PROT_sum, x = "RNA_Mean", y = "PROT_Mean", 
          add = "reg.line", conf.int = TRUE, 
          xlab = "log2 RNA Mean", ylab = "Protein Mean", size=0.05, title="Correlation of RNA and Protein Mean Expression")+
          stat_cor(method = "pearson", label.x = 3, label.y = 15) + 
          geom_point(color='#21908dff', size=0.05)+ 
          geom_text(label=RNA_PROT_sum$Mean_outlier, nudge_x = 0.35, size = 3)+
          theme(plot.title = element_text(hjust=0.5))
dev.off()

pdf(file='Results/Median_cor_all.pdf',
    width=15,height=10,useDingbats=FALSE)
ggscatter(RNA_PROT_sum, x = "RNA_Median", y = "PROT_Median", 
          add = "reg.line", conf.int = TRUE, 
          xlab = "log2 RNA Median", ylab = "Protein Median", size=0.05, title="Correlation of RNA and Protein Median Expression")+
  stat_cor(method = "pearson", label.x = 3, label.y = 15) + 
  geom_point(color='#21908dff', size=0.05)+ 
  geom_text(label=RNA_PROT_sum$Median_outlier, nudge_x = 0.35, size = 3)+
  theme(plot.title = element_text(hjust=0.5))
dev.off()


RNA_PROT_sum_corum = RNA_PROT_sum[RNA_PROT_sum$CORUM == "Yes",]
RNA_PROT_sum_non_corum= RNA_PROT_sum[RNA_PROT_sum$CORUM == "No",]

mean_cor = cor(RNA_PROT_sum_corum$RNA_Mean, RNA_PROT_sum_corum$PROT_Mean, method='pearson')
median_cor = cor(RNA_PROT_sum_corum$RNA_Median, RNA_PROT_sum_corum$PROT_Median, method='pearson')

pdf(file='Results/Mean_cor_CORUM.pdf',
    width=15,height=10,useDingbats=FALSE)
ggscatter(RNA_PROT_sum_corum, x = "RNA_Mean", y = "PROT_Mean", 
          add = "reg.line", conf.int = TRUE, 
          xlab = "log2 RNA Mean", ylab = "Protein Mean", size=0.05, title="Correlation of RNA and Protein Mean Expression in Complex")+
  stat_cor(method = "pearson", label.x = 3, label.y = 10) + 
  geom_point(color='#21908dff', size=0.05)+ 
  geom_text(label=RNA_PROT_sum_corum$Mean_outlier, nudge_x = 0.35, size = 3)+
  theme(plot.title = element_text(hjust=0.5))
dev.off()

pdf(file='Results/Median_cor_CORUM.pdf',
    width=15,height=10,useDingbats=FALSE)
ggscatter(RNA_PROT_sum_corum, x = "RNA_Median", y = "PROT_Median", 
          add = "reg.line", conf.int = TRUE, 
          xlab = "log2 RNA Median", ylab = "Protein Median", size=0.05, title="Correlation of RNA and Protein Median Expressionin in Complex")+
  stat_cor(method = "pearson", label.x = 3, label.y = 10) + 
  geom_point(color='#21908dff', size=0.05)+ 
  geom_text(label=RNA_PROT_sum_corum$Median_outlier, nudge_x = 0.35, size = 3)+
  theme(plot.title = element_text(hjust=0.5))
dev.off()


mean_cor = cor(RNA_PROT_sum_non_corum$RNA_Mean, RNA_PROT_sum_non_corum$PROT_Mean, method='pearson')
median_cor = cor(RNA_PROT_sum_non_corum$RNA_Median, RNA_PROT_sum_non_corum$PROT_Median, method='pearson')

pdf(file='Results/Mean_cor_nonCORUM.pdf',
    width=15,height=10,useDingbats=FALSE)
ggscatter(RNA_PROT_sum_non_corum, x = "RNA_Mean", y = "PROT_Mean", 
          add = "reg.line", conf.int = TRUE, 
          xlab = "log2 RNA Mean", ylab = "Protein Mean", size=0.05, title="Correlation of RNA and Protein Mean Expression not in Complex")+
  stat_cor(method = "pearson", label.x = 3, label.y = 15) + 
  geom_point(color='#21908dff', size=0.05)+ 
  geom_text(label=RNA_PROT_sum_non_corum$Mean_outlier, nudge_x = 0.35, size = 3)+
  theme(plot.title = element_text(hjust=0.5))
dev.off()

pdf(file='Results/Median_cor_nonCORUM.pdf',
    width=15,height=10,useDingbats=FALSE)
ggscatter(RNA_PROT_sum_non_corum, x = "RNA_Median", y = "PROT_Median", 
          add = "reg.line", conf.int = TRUE, 
          xlab = "log2 RNA Median", ylab = "Protein Median", size=0.05, title="Correlation of RNA and Protein Median Expressionin not in Complex")+
  stat_cor(method = "pearson", label.x = 3, label.y = 15) + 
  geom_point(color='#21908dff', size=0.05)+ 
  geom_text(label=RNA_PROT_sum_non_corum$Median_outlier, nudge_x = 0.35, size = 3)+
  theme(plot.title = element_text(hjust=0.5))
dev.off()




theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

