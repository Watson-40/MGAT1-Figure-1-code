#3/9/2025
#Josh Watson
#MGAT1-Figure-1-code
#MGAT1-Mediated Glycosylation Orchestrates Immune Checkpoints and Antitumor Immunity

#MGAT Figure 1A
#Heatmap of Immune Genes in BLBC with MGAT1 TPM as top annotation
#
library(readxl)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)

#will need to acquire rawTPM data from TCGA GDC for BRCA
#https://portal.gdc.cancer.gov/projects/TCGA-BRCA
TCGA <- read_excel("~/Downloads/TCGA_BRCA_FPKMandTPM.xlsx", sheet = "rawTPM")

BLBC_TCGA <- TCGA[,c(1,which(TCGA[1,] == "Basal"))]

#For inital clustering- run script with Immunegenes_Em.txt
#For Figure 1A- run script with Emory_70immunes.txt

#immune <- read.table("Downloads/Immunegenes_Em.txt",header = T)
immune <- read.table("~/Downloads/Emory_70immunes.txt",header = T)
meta <- read_excel(path = "~/Downloads/11.8.2023 function category MGAT1 gene list.xlsx", sheet = "Sheet1")
#change BLBC_TCGA and TCGA
immune_BLBC <- subset(TCGA, TCGA %in% immune$ImmuneGenes)
immune_BLBC <- data.frame(immune_BLBC, row.names = 1)
matimmune  <- as.matrix(immune_BLBC)
matimmune <- apply(matimmune, 2, as.numeric)
matimmune <- t(scale(t(matimmune)))
rownames(matimmune) <- rownames(immune_BLBC)
#change BLBC_TCGA and TCGA
immune_BLBCmeta <- data.frame(t(data.frame(subset(TCGA, TCGA %in% c("PAM50", "MGAT1")),row.names = 1)))
immune_BLBCmeta$MGAT1 <- as.numeric(immune_BLBCmeta$MGAT1)

MGAT1_color_breaks <- c(16, 130, 242)
MGAT1_color_panel <- colorRamps::blue2yellow(length(MGAT1_color_breaks))
MGAT1_colors <- colorRamp2(rev(MGAT1_color_breaks), rev(MGAT1_color_panel));

PAM50 <- c("blue","skyblue","red","hotpink","green")
names(PAM50) <- unique(immune_BLBCmeta$PAM50)
colors <- list(MGAT1 = MGAT1_colors,
               PAM50 = PAM50)
immune_BLBCmeta$PAM50 <- factor(immune_BLBCmeta$PAM50, levels = c("Her2","Normal","LumB","Basal","LumA"
))
a0 <- HeatmapAnnotation(#MGAT1= immune_BLBCmeta$MGAT1,
  PAM50 = immune_BLBCmeta$PAM50,
  col =list(
    #MGAT1 = colors$MGAT1,
    PAM50 = colors$PAM50
  ),
  show_annotation_name = TRUE,
  annotation_name_side = "right")
colors2 <- c("#648ace",
             "#ab973d",
             "red",
             "#5ba962",
             "yellow",
             "navy")
names(colors2) <- unique(meta$`function`)
meta <- subset(meta, Gene %in% rownames(matimmune))
meta$`function` <- factor(meta$`function`, levels = unique(meta$`function`))
a2 <- rowAnnotation(#MGAT1= immune_BLBCmeta$MGAT1,
  Function = meta$`function`,
  col =list(
    Function = colors2
  ),
  show_annotation_name = TRUE)
matimmune <- matimmune[match(meta$Gene, rownames(matimmune)),]
matimmune[which(matimmune > 3)] <-3

col_fun = colorRamp2(c(-2, 0, 3), c("blue", "white", "red"))
col_fun(seq(-2, 3))

#initial clustering block
#only run following for initial clustering#
library(dendextend)
hclust_BLBC_all129 <- hclust(dist(t(matimmune)))
treeee <- cutree(hclust_BLBC_all129, k = 4)
table(treeee)
dendro <- as.dendrogram(hclust_BLBC_all129) 
dendrocol <- color_branches(dendro, k=4)
#initial clustering block

a <- Heatmap(matrix = matimmune, 
             top_annotation = a0,
             left_annotation = a2,
             cluster_rows = FALSE,
             #row_order = match(meta$Gene, rownames(matimmune)),
             show_row_dend = FALSE,
             col = col_fun,
             #uncomment following for initial clustering 
             #cluster_columns = dendrocol,
             column_split = immune_BLBCmeta$PAM50,
             #column_order = order(immune_BLBCmeta$MGAT1),
             show_row_names = TRUE,show_column_names = FALSE,
             show_heatmap_legend = TRUE,
             row_names_gp = gpar(fontsize = 9),
             heatmap_legend_param = list(title = "Z-score")
) 
pdf("Downloads/1Araw.pdf",width = 14, height = 10)
draw(a)
dev.off()


#MGAT1 Figure 1B and 1C
library(cogena)
library(GSVA)
library(reshape2)

GOimmune <- gmt2list("~/Downloads/immune_10.18.gmt")
ssGSEA_BLBC_GOimmune <- gsva(expr = BLBCmat, gset.idx.list = GOimmune,
                             method = "ssgsea")
ssGSEA_BLBC_GOimmune <- t(ssGSEA_BLBC_GOimmune)
BLBC_immunesubsets <- read.csv("~/BLBC_immunesubsets.csv")
ssGSEA_BLBC_GOimmune <- merge(ssGSEA_BLBC_GOimmune,BLBC_immunesubsets, by.x =0, by.y = "BLBC.sample")

ssGSEA_BLBC_GOimmune <- subset(ssGSEA_BLBC_GOimmune, select = -c(GOBP_IMMUNE_RESPONSE_TO_TUMOR_CELL,                                         
                                                                 GOBP_LYMPHOCYTE_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE,                     
                                                                 GOBP_REGULATION_OF_T_CELL_MEDIATED_IMMUNE_RESPONSE_TO_TUMOR_CELL,           
                                                                 GOBP_T_CELL_MEDIATED_IMMUNE_RESPONSE_TO_TUMOR_CELL  ))

ssGSEA_BLBC_GOimmune <-  melt(ssGSEA_BLBC_GOimmune, measure.vars = colnames(ssGSEA_BLBC_GOimmune)[2:7],
                              variable.name = "GOterm", value.name = "ssGSEAenrichment")
ssGSEA_BLBC_GOimmune <- subset(ssGSEA_BLBC_GOimmune, category != "outlier")
View(compare_means(ssGSEAenrichment ~ category, data = ssGSEA_BLBC_GOimmune, group.by = "GOterm"))
BLBCImmstats <- compare_means(ssGSEAenrichment ~ category, data = ssGSEA_BLBC_GOimmune, group.by = "GOterm")
ssGSEA_BLBC_GOimmune$GOterm <- factor(ssGSEA_BLBC_GOimmune$GOterm,levels = rev(unique(ssGSEA_BLBC_GOimmune$GOterm)))



b <- ggplot(ssGSEA_BLBC_GOimmune, aes(x=GOterm, y=ssGSEAenrichment,fill = category, color = category)) + 
  geom_boxplot( outlier.size = 0.1,alpha = 0,width=0.75, lwd=0.25)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
              size=.1)+
  theme_classic()+
  theme(
    plot.margin = margin(0.1, 1, 0, 1.65, "cm"),
    axis.text.y = element_text(colour=c("black"),size=11),
    axis.text.x = element_text(colour=c("black"), size=11),
    axis.title = element_text(colour=c("black"), size=11),
    axis.title.x = element_text(colour=c("black"), size=11),
    legend.text = element_text(colour=c("black"), size=11),
    legend.title = element_text(colour=c("black"), size=11)
  )+coord_flip()+ylim(c(-0.36,0.73))+
  labs(title = "TCGA BLBC", y = "ssGSEA")+
  guides(fill = guide_legend(reverse = TRUE),col = guide_legend(reverse = TRUE))+
  scale_fill_manual(values = c("blue","red"))+
  scale_color_manual(values = c("blue","red"))


#10/20
C1020 <- gmt2list("~/Downloads/10.20.2023_1cGeneset.gmt")
ssGSEA_BLBC_C1020 <- gsva(expr = BLBCmat, gset.idx.list = C1020,
                          method = "ssgsea")
ssGSEA_BLBC_C1020 <- t(ssGSEA_BLBC_C1020)
BLBC_immunesubsets <- read.csv("~/BLBC_immunesubsets.csv")
ssGSEA_BLBC_C1020 <- merge(ssGSEA_BLBC_C1020,BLBC_immunesubsets, by.x =0, by.y = "BLBC.sample")

ssGSEA_BLBC_C1020 <-  melt(ssGSEA_BLBC_C1020, measure.vars = colnames(ssGSEA_BLBC_C1020)[2:9],
                           variable.name = "GOterm", value.name = "ssGSEAenrichment")
ssGSEA_BLBC_C1020 <- subset(ssGSEA_BLBC_C1020, category != "outlier")
View(compare_means(ssGSEAenrichment ~ category, data = ssGSEA_BLBC_C1020, group.by = "GOterm"))
BLBCGostats <- compare_means(ssGSEAenrichment ~ category, data = ssGSEA_BLBC_C1020, group.by = "GOterm")

BLBCstat <- rbind(BLBCImmstats, BLBCGostats)
write.csv(BLBCstat, file = "BLBCwilcoxstat.csv", quote = F)

ssGSEA_BLBC_C1020$GOterm <- factor(ssGSEA_BLBC_C1020$GOterm,levels = rev(unique(ssGSEA_BLBC_C1020$GOterm)))

c <- ggplot(ssGSEA_BLBC_C1020, aes(x=GOterm, y=ssGSEAenrichment,fill = category, color = category)) + 
  geom_boxplot( outlier.size = 0.1,alpha = 0,width=0.75, lwd=0.25)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
              size=.1)+
  theme_classic()+
  theme(
    plot.margin = margin(0.1, 1, 0, 1.65, "cm"),
    axis.text.y = element_text(colour=c("black"),size=11),
    axis.text.x = element_text(colour=c("black"), size=11),
    axis.title = element_text(colour=c("black"), size=11),
    axis.title.x = element_text(colour=c("black"), size=11),
    legend.text = element_text(colour=c("black"), size=11),
    legend.title = element_text(colour=c("black"), size=11)
  )+coord_flip()+ylim(c(-0.36,0.73))+
  labs(title = "TCGA BLBC", y = "ssGSEA")+
  guides(fill = guide_legend(reverse = TRUE),col = guide_legend(reverse = TRUE))+
  scale_fill_manual(values = c("blue","red"))+
  scale_color_manual(values = c("blue","red"))
#10/20
pdf("1BCraw_4.5.pdf",width = 12,  height = 6)
ggarrange(b,c, nrow = 2, common.legend = T, align = "v", legend = "right", heights = c(0.8,1))

dev.off()


#MGAT1 Figure 1D
#correlation plot in BLBC- immune genes and n-glycosylation genes
#First Batch with TCGA

immune <- read.table("Downloads/Immunegenes_Em.txt",header = T)
meta <- read_excel(path = "~/Downloads/11.8.2023 function category MGAT1 gene list.xlsx", sheet = "Sheet1")

immune_BLBC <- subset(BLBC_TCGA, TCGA %in% immune$ImmuneGenes)
immune_BLBC <- data.frame(immune_BLBC, row.names = 1)
library(ggcorrplot)

View(cor(t(immuneBLBC2), t(NglycoBLBC2), method = "pearson"))
F1D_TCGApearson <- cor(t(immuneBLBC2), t(NglycoBLBC2), method = "pearson")
ppearson <- cor_pmat(cbind(t(immuneBLBC2), t(NglycoBLBC2)), method = "pearson")
ppearson <- ppearson[rownames(F1D_TCGApearson),colnames(F1D_TCGApearson)]
write.csv(F1D_TCGApearson, file = "F1D_TCGA_Correlations.csv", quote = F)
write.csv(ppearson, file = "F1D_TCGA_Pvals.csv", quote = F)


F1D_TCGA <- read.csv("~/F1D_TCGA_Correlations.csv", header = T, row.names = 1)
pc <- read.csv("~/F1D_TCGA_Pvals.csv", header =T, row.names = 1)

d <- ggcorrplot(F1D_TCGA,
                p.mat = pc,colors = c("navy","white","red"),
                sig.level = 0.05,
                legend.title = "Correlation",
                tl.cex = 8,
                tl.srt = 90,
                #title = "Spearman correlation",
                tl.col = "black",
                insig = "blank")
pdf("Downloads/1draw.pdf", width = 13, height = 15)
d
dev.off()

