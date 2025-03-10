#11.17.2024
##follow up on MGAT1
##suuplemental figure scipt

#Divvy up all other subtypes
library(readxl)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
#will need to acquire rawTPM data from TCGA GDC for BRCA
#https://portal.gdc.cancer.gov/projects/TCGA-BRCA
TCGA <- read_excel("~/Downloads/TCGA_BRCA_FPKMandTPM.xlsx", sheet = "rawTPM")

LumA_TCGA <- TCGA[,c(1,which(TCGA[1,] == "LumA"))]
LumB_TCGA <- TCGA[,c(1,which(TCGA[1,] == "LumB"))]
Her2_TCGA <- TCGA[,c(1,which(TCGA[1,] == "Her2"))]
Normal_TCGA <- TCGA[,c(1,which(TCGA[1,] == "Normal"))]



immune <- read.table("Downloads/Immunegenes_Em.txt",header = T)
#immune <- read.table("~/Downloads/Emory_70immunes.txt",header = T)
meta <- read_excel(path = "~/Downloads/11.8.2023 function category MGAT1 gene list.xlsx", sheet = "Sheet1")
#change LumA_TCGA and TCGA
immune_LumA <- subset(LumA_TCGA, TCGA %in% immune$ImmuneGenes)
immune_LumA <- data.frame(immune_LumA, row.names = 1)
matimmune  <- as.matrix(immune_LumA)
matimmune <- apply(matimmune, 2, as.numeric)
matimmune <- t(scale(t(matimmune)))
rownames(matimmune) <- rownames(immune_LumA)
#change LumA_TCGA and TCGA
immune_LumAmeta <- data.frame(t(data.frame(subset(LumA_TCGA, TCGA %in% c("PAM50", "MGAT1")),row.names = 1)))
immune_LumAmeta$MGAT1 <- as.numeric(immune_LumAmeta$MGAT1)

MGAT1_color_breaks <- c(16, 130, 242)
MGAT1_color_panel <- colorRamps::blue2yellow(length(MGAT1_color_breaks))
MGAT1_colors <- colorRamp2(rev(MGAT1_color_breaks), rev(MGAT1_color_panel));

PAM50 <- c("blue","skyblue","red","hotpink","green")
names(PAM50) <- unique(immune_LumAmeta$PAM50)
colors <- list(MGAT1 = MGAT1_colors,
               PAM50 = PAM50)
immune_LumAmeta$PAM50 <- factor(immune_LumAmeta$PAM50, levels = c("Her2","Normal","LumB","Basal","LumA"
))
a0 <- HeatmapAnnotation(#MGAT1= immune_LumAmeta$MGAT1,
  PAM50 = immune_LumAmeta$PAM50,
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
a2 <- rowAnnotation(#MGAT1= immune_LumAmeta$MGAT1,
  Function = meta$`function`,
  col =list(
    Function = colors2
  ),
  show_annotation_name = TRUE)
matimmune <- matimmune[match(meta$Gene, rownames(matimmune)),]
matimmune[which(matimmune > 3)] <-3

col_fun = colorRamp2(c(-2, 0, 3), c("blue", "white", "red"))
col_fun(seq(-2, 3))


library(dendextend)
hclust_LumA_all129 <- hclust(dist(t(matimmune)))
treeeeLumA <- cutree(hclust_LumA_all129, k = 2)
table(treeeeLumA)
dendro <- as.dendrogram(hclust_LumA_all129) 
dendrocol <- color_branches(dendro, k=2)

a <- Heatmap(matrix = matimmune, 
             #top_annotation = a0,
             left_annotation = a2,
             cluster_rows = FALSE,
             #row_order = match(meta$Gene, rownames(matimmune)),
             show_row_dend = FALSE,
             col = col_fun,
             cluster_columns = dendrocol,
             #column_split = immune_LumAmeta$PAM50,
             #column_order = order(immune_LumAmeta$MGAT1),
             show_row_names = TRUE,show_column_names = FALSE,
             show_heatmap_legend = TRUE,
             row_names_gp = gpar(fontsize = 9),
             heatmap_legend_param = list(title = "Z-score")
) 

#LumB
immune_LumB <- subset(LumB_TCGA, TCGA %in% immune$ImmuneGenes)
immune_LumB <- data.frame(immune_LumB, row.names = 1)
matimmune  <- as.matrix(immune_LumB)
matimmune <- apply(matimmune, 2, as.numeric)
matimmune <- t(scale(t(matimmune)))
rownames(matimmune) <- rownames(immune_LumB)
#change LumB_TCGA and TCGA
immune_LumBmeta <- data.frame(t(data.frame(subset(LumB_TCGA, TCGA %in% c("PAM50", "MGAT1")),row.names = 1)))
immune_LumBmeta$MGAT1 <- as.numeric(immune_LumBmeta$MGAT1)

MGAT1_color_breaks <- c(16, 130, 242)
MGAT1_color_panel <- colorRamps::blue2yellow(length(MGAT1_color_breaks))
MGAT1_colors <- colorRamp2(rev(MGAT1_color_breaks), rev(MGAT1_color_panel));

PAM50 <- c("blue","skyblue","red","hotpink","green")
names(PAM50) <- unique(immune_LumBmeta$PAM50)
colors <- list(MGAT1 = MGAT1_colors,
               PAM50 = PAM50)
immune_LumBmeta$PAM50 <- factor(immune_LumBmeta$PAM50, levels = c("Her2","Normal","LumA","Basal","LumB"
))
a0 <- HeatmapAnnotation(#MGAT1= immune_LumBmeta$MGAT1,
  PAM50 = immune_LumBmeta$PAM50,
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
a2 <- rowAnnotation(#MGAT1= immune_LumBmeta$MGAT1,
  Function = meta$`function`,
  col =list(
    Function = colors2
  ),
  show_annotation_name = TRUE)
matimmune <- matimmune[match(meta$Gene, rownames(matimmune)),]
matimmune[which(matimmune > 3)] <-3

col_fun = colorRamp2(c(-2, 0, 3), c("blue", "white", "red"))
col_fun(seq(-2, 3))


library(dendextend)
hclust_LumB_all129 <- hclust(dist(t(matimmune)))
treeeeLumB <- cutree(hclust_LumB_all129, k = 2)
table(treeeeLumB)
dendro <- as.dendrogram(hclust_LumB_all129) 
dendrocol <- color_branches(dendro, k=2)

a <- Heatmap(matrix = matimmune, 
             #top_annotation = a0,
             left_annotation = a2,
             cluster_rows = FALSE,
             #row_order = match(meta$Gene, rownames(matimmune)),
             show_row_dend = FALSE,
             col = col_fun,
             cluster_columns = dendrocol,
             #column_split = immune_LumBmeta$PAM50,
             #column_order = order(immune_LumBmeta$MGAT1),
             show_row_names = TRUE,show_column_names = FALSE,
             show_heatmap_legend = TRUE,
             row_names_gp = gpar(fontsize = 9),
             heatmap_legend_param = list(title = "Z-score")
) 
#Her2
immune_Her2 <- subset(Her2_TCGA, TCGA %in% immune$ImmuneGenes)
immune_Her2 <- data.frame(immune_Her2, row.names = 1)
matimmune  <- as.matrix(immune_Her2)
matimmune <- apply(matimmune, 2, as.numeric)
matimmune <- t(scale(t(matimmune)))
rownames(matimmune) <- rownames(immune_Her2)
#change Her2_TCGA and TCGA
immune_Her2meta <- data.frame(t(data.frame(subset(Her2_TCGA, TCGA %in% c("PAM50", "MGAT1")),row.names = 1)))
immune_Her2meta$MGAT1 <- as.numeric(immune_Her2meta$MGAT1)

MGAT1_color_breaks <- c(16, 130, 242)
MGAT1_color_panel <- colorRamps::blue2yellow(length(MGAT1_color_breaks))
MGAT1_colors <- colorRamp2(rev(MGAT1_color_breaks), rev(MGAT1_color_panel));

PAM50 <- c("blue","skyblue","red","hotpink","green")
names(PAM50) <- unique(immune_Her2meta$PAM50)
colors <- list(MGAT1 = MGAT1_colors,
               PAM50 = PAM50)
immune_Her2meta$PAM50 <- factor(immune_Her2meta$PAM50, levels = c("Her2","Normal","LumA","Basal","LumB"
))
a0 <- HeatmapAnnotation(#MGAT1= immune_Her2meta$MGAT1,
  PAM50 = immune_Her2meta$PAM50,
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
a2 <- rowAnnotation(#MGAT1= immune_Her2meta$MGAT1,
  Function = meta$`function`,
  col =list(
    Function = colors2
  ),
  show_annotation_name = TRUE)
matimmune <- matimmune[match(meta$Gene, rownames(matimmune)),]
matimmune[which(matimmune > 3)] <-3

col_fun = colorRamp2(c(-2, 0, 3), c("blue", "white", "red"))
col_fun(seq(-2, 3))


library(dendextend)
hclust_Her2_all129 <- hclust(dist(t(matimmune)))
treeeeHer2 <- cutree(hclust_Her2_all129, k = 2)
table(treeeeHer2)
dendro <- as.dendrogram(hclust_Her2_all129) 
dendrocol <- color_branches(dendro, k=2)

a <- Heatmap(matrix = matimmune, 
             #top_annotation = a0,
             left_annotation = a2,
             cluster_rows = FALSE,
             #row_order = match(meta$Gene, rownames(matimmune)),
             show_row_dend = FALSE,
             col = col_fun,
             cluster_columns = dendrocol,
             #column_split = immune_Her2meta$PAM50,
             #column_order = order(immune_Her2meta$MGAT1),
             show_row_names = TRUE,show_column_names = FALSE,
             show_heatmap_legend = TRUE,
             row_names_gp = gpar(fontsize = 9),
             heatmap_legend_param = list(title = "Z-score")
)
#Normal
immune_Normal <- subset(Normal_TCGA, TCGA %in% immune$ImmuneGenes)
immune_Normal <- data.frame(immune_Normal, row.names = 1)
matimmune  <- as.matrix(immune_Normal)
matimmune <- apply(matimmune, 2, as.numeric)
matimmune <- t(scale(t(matimmune)))
rownames(matimmune) <- rownames(immune_Normal)
#change Normal_TCGA and TCGA
immune_Normalmeta <- data.frame(t(data.frame(subset(Normal_TCGA, TCGA %in% c("PAM50", "MGAT1")),row.names = 1)))
immune_Normalmeta$MGAT1 <- as.numeric(immune_Normalmeta$MGAT1)

MGAT1_color_breaks <- c(16, 130, 242)
MGAT1_color_panel <- colorRamps::blue2yellow(length(MGAT1_color_breaks))
MGAT1_colors <- colorRamp2(rev(MGAT1_color_breaks), rev(MGAT1_color_panel));

PAM50 <- c("blue","skyblue","red","hotpink","green")
names(PAM50) <- unique(immune_Normalmeta$PAM50)
colors <- list(MGAT1 = MGAT1_colors,
               PAM50 = PAM50)
immune_Normalmeta$PAM50 <- factor(immune_Normalmeta$PAM50, levels = c("Her2","Normal","LumA","Basal","LumB"
))
a0 <- HeatmapAnnotation(#MGAT1= immune_Normalmeta$MGAT1,
  PAM50 = immune_Normalmeta$PAM50,
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
a2 <- rowAnnotation(#MGAT1= immune_Normalmeta$MGAT1,
  Function = meta$`function`,
  col =list(
    Function = colors2
  ),
  show_annotation_name = TRUE)
matimmune <- matimmune[match(meta$Gene, rownames(matimmune)),]
matimmune[which(matimmune > 3)] <-3

col_fun = colorRamp2(c(-2, 0, 3), c("blue", "white", "red"))
col_fun(seq(-2, 3))


library(dendextend)
hclust_Normal_all129 <- hclust(dist(t(matimmune)))
treeeeNormal <- cutree(hclust_Normal_all129, k = 2)
table(treeeeNormal)
dendro <- as.dendrogram(hclust_Normal_all129) 
dendrocol <- color_branches(dendro, k=2)

 Heatmap(matrix = matimmune, 
             #top_annotation = a0,
             left_annotation = a2,
             cluster_rows = FALSE,
             #row_order = match(meta$Gene, rownames(matimmune)),
             show_row_dend = FALSE,
             col = col_fun,
             cluster_columns = dendrocol,
             #column_split = immune_Normalmeta$PAM50,
             #column_order = order(immune_Normalmeta$MGAT1),
             show_row_names = TRUE,show_column_names = FALSE,
             show_heatmap_legend = TRUE,
             row_names_gp = gpar(fontsize = 9),
             heatmap_legend_param = list(title = "Z-score")
)
 
 
 
 
 
 
 
 
 library(cogena)
 library(GSVA)
 
 LumA_TCGA <- LumA_TCGA[!duplicated(LumA_TCGA$TCGA),]
 LumA <- data.frame(LumA_TCGA, row.names = 1)
 LumA <- LumA[-1,]
 LumAmat <- apply(LumA, 2, as.numeric)
 rownames(LumAmat) <- rownames(LumA)
 GOimmune <- gmt2list("~/Downloads/immune_10.18.gmt")
 ssGSEA_LumA_GOimmune <- gsva(expr = LumAmat, gset.idx.list = GOimmune,
                              method = "ssgsea")
 ssGSEA_LumA_GOimmune <- t(ssGSEA_LumA_GOimmune)
 LumA_immunesubsets <- read.csv("~/LumAcat.csv")
 ssGSEA_LumA_GOimmune <- merge(ssGSEA_LumA_GOimmune,LumA_immunesubsets, by.x =0, by.y = "X")
 library(reshape2)
 ssGSEA_LumA_GOimmune <- subset(ssGSEA_LumA_GOimmune, select = -c(GOBP_IMMUNE_RESPONSE_TO_TUMOR_CELL,                                         
                                                                  GOBP_LYMPHOCYTE_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE,                     
                                                                  GOBP_REGULATION_OF_T_CELL_MEDIATED_IMMUNE_RESPONSE_TO_TUMOR_CELL,           
                                                                  GOBP_T_CELL_MEDIATED_IMMUNE_RESPONSE_TO_TUMOR_CELL  ))
 
 ssGSEA_LumA_GOimmune <-  melt(ssGSEA_LumA_GOimmune, measure.vars = colnames(ssGSEA_LumA_GOimmune)[2:7],
                               variable.name = "GOterm", value.name = "ssGSEAenrichment")
 #ssGSEA_LumA_GOimmune <- subset(ssGSEA_LumA_GOimmune, category != "outlier")
 View(compare_means(ssGSEAenrichment ~ x, data = ssGSEA_LumA_GOimmune, group.by = "GOterm"))
 LumAimmstat <- compare_means(ssGSEAenrichment ~ x, data = ssGSEA_LumA_GOimmune, group.by = "GOterm")
 ssGSEA_LumA_GOimmune$GOterm <- factor(ssGSEA_LumA_GOimmune$GOterm,levels = rev(unique(ssGSEA_LumA_GOimmune$GOterm)))
 ssGSEA_LumA_GOimmune$x <- factor(ssGSEA_LumA_GOimmune$x, levels = unique(ssGSEA_LumA_GOimmune$x))
 
 
 b <- ggplot(ssGSEA_LumA_GOimmune, aes(x=GOterm, y=ssGSEAenrichment,fill = x, color = x)) + 
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
   labs(title = "TCGA LumA", y = "ssGSEA")+
   guides(fill = guide_legend(reverse = TRUE),col = guide_legend(reverse = TRUE))+
   scale_fill_manual(values = c("blue","red"))+
   scale_color_manual(values = c("blue","red"))
 
 
 #10/20
 C1020 <- gmt2list("~/Downloads/10.20.2023_1cGeneset.gmt")
 ssGSEA_LumA_C1020 <- gsva(expr = LumAmat, gset.idx.list = C1020,
                           method = "ssgsea")
 ssGSEA_LumA_C1020 <- t(ssGSEA_LumA_C1020)
 #LumA_immunesubsets <- read.csv("~/BLBC_immunesubsets.csv")
 ssGSEA_LumA_C1020 <- merge(ssGSEA_LumA_C1020,LumA_immunesubsets, by.x =0, by.y = "X")
 
 ssGSEA_LumA_C1020 <-  melt(ssGSEA_LumA_C1020, measure.vars = colnames(ssGSEA_LumA_C1020)[2:9],
                            variable.name = "GOterm", value.name = "ssGSEAenrichment")
 #ssGSEA_LumA_C1020 <- subset(ssGSEA_LumA_C1020, category != "outlier")
 View(compare_means(ssGSEAenrichment ~ x, data = ssGSEA_LumA_C1020, group.by = "GOterm"))
 
 ssGSEA_LumA_C1020$GOterm <- factor(ssGSEA_LumA_C1020$GOterm,levels = rev(unique(ssGSEA_LumA_C1020$GOterm)))
 ssGSEA_LumA_C1020$x <- factor(ssGSEA_LumA_C1020$x, levels = unique(ssGSEA_LumA_C1020$x))
 
 LumAimmstat <- compare_means(ssGSEAenrichment ~ x, data = ssGSEA_LumA_GOimmune, group.by = "GOterm")
 LumAGostat <- compare_means(ssGSEAenrichment ~ x, data = ssGSEA_LumA_C1020, group.by = "GOterm")
 LumAstats <- rbind(LumAimmstat, LumAGostat)
 write.csv(LumAstats, file = "LumAwilcoxstat.csv", quote = F)
 
 
 c <- ggplot(ssGSEA_LumA_C1020, aes(x=GOterm, y=ssGSEAenrichment,fill = x, color = x)) + 
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
   labs(title = "TCGA LumA", y = "ssGSEA")+
   guides(fill = guide_legend(reverse = TRUE),col = guide_legend(reverse = TRUE))+
   scale_fill_manual(values = c("blue","red"))+
   scale_color_manual(values = c("blue","red"))
 #11/25
 pdf("1BCraw_LumA.pdf",width = 12,  height = 6)
 ggarrange(b,c, nrow = 2, common.legend = T, align = "v", legend = "right", heights = c(0.8,1))
 
 dev.off()
 
 
 
 
 
 LumB_TCGA <- LumB_TCGA[!duplicated(LumB_TCGA$TCGA),]
 LumB <- data.frame(LumB_TCGA, row.names = 1)
 LumB <- LumB[-1,]
 LumBmat <- apply(LumB, 2, as.numeric)
 rownames(LumBmat) <- rownames(LumB)
 GOimmune <- gmt2list("~/Downloads/immune_10.18.gmt")
 ssGSEA_LumB_GOimmune <- gsva(expr = LumBmat, gset.idx.list = GOimmune,
                              method = "ssgsea")
 ssGSEA_LumB_GOimmune <- t(ssGSEA_LumB_GOimmune)
 LumB_immunesubsets <- read.csv("~/LumBcat.csv")
 ssGSEA_LumB_GOimmune <- merge(ssGSEA_LumB_GOimmune,LumB_immunesubsets, by.x =0, by.y = "X")
 library(reshape2)
 ssGSEA_LumB_GOimmune <- subset(ssGSEA_LumB_GOimmune, select = -c(GOBP_IMMUNE_RESPONSE_TO_TUMOR_CELL,                                         
                                                                  GOBP_LYMPHOCYTE_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE,                     
                                                                  GOBP_REGULATION_OF_T_CELL_MEDIATED_IMMUNE_RESPONSE_TO_TUMOR_CELL,           
                                                                  GOBP_T_CELL_MEDIATED_IMMUNE_RESPONSE_TO_TUMOR_CELL  ))
 
 ssGSEA_LumB_GOimmune <-  melt(ssGSEA_LumB_GOimmune, measure.vars = colnames(ssGSEA_LumB_GOimmune)[2:7],
                               variable.name = "GOterm", value.name = "ssGSEAenrichment")
 #ssGSEA_LumB_GOimmune <- subset(ssGSEA_LumB_GOimmune, category != "outlier")
 View(compare_means(ssGSEAenrichment ~ x, data = ssGSEA_LumB_GOimmune, group.by = "GOterm"))
 
 ssGSEA_LumB_GOimmune$GOterm <- factor(ssGSEA_LumB_GOimmune$GOterm,levels = rev(unique(ssGSEA_LumB_GOimmune$GOterm)))
 ssGSEA_LumB_GOimmune$x <- factor(ssGSEA_LumB_GOimmune$x, levels = unique(ssGSEA_LumB_GOimmune$x))
 
 
 b <- ggplot(ssGSEA_LumB_GOimmune, aes(x=GOterm, y=ssGSEAenrichment,fill = x, color = x)) + 
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
   )+coord_flip()+ylim(c(-0.38,0.85))+
   labs(title = "TCGA LumB", y = "ssGSEA")+
   guides(fill = guide_legend(reverse = TRUE),col = guide_legend(reverse = TRUE))+
   scale_fill_manual(values = c("blue","red"))+
   scale_color_manual(values = c("blue","red"))
 
 
 #10/20
 C1020 <- gmt2list("~/Downloads/10.20.2023_1cGeneset.gmt")
 ssGSEA_LumB_C1020 <- gsva(expr = LumBmat, gset.idx.list = C1020,
                           method = "ssgsea")
 ssGSEA_LumB_C1020 <- t(ssGSEA_LumB_C1020)
 #LumB_immunesubsets <- read.csv("~/BLBC_immunesubsets.csv")
 ssGSEA_LumB_C1020 <- merge(ssGSEA_LumB_C1020,LumB_immunesubsets, by.x =0, by.y = "X")
 
 ssGSEA_LumB_C1020 <-  melt(ssGSEA_LumB_C1020, measure.vars = colnames(ssGSEA_LumB_C1020)[2:9],
                            variable.name = "GOterm", value.name = "ssGSEAenrichment")
 #ssGSEA_LumB_C1020 <- subset(ssGSEA_LumB_C1020, category != "outlier")
 View(compare_means(ssGSEAenrichment ~ x, data = ssGSEA_LumB_C1020, group.by = "GOterm"))
 
 ssGSEA_LumB_C1020$GOterm <- factor(ssGSEA_LumB_C1020$GOterm,levels = rev(unique(ssGSEA_LumB_C1020$GOterm)))
 ssGSEA_LumB_C1020$x <- factor(ssGSEA_LumB_C1020$x, levels = unique(ssGSEA_LumB_C1020$x))
 
 c <- ggplot(ssGSEA_LumB_C1020, aes(x=GOterm, y=ssGSEAenrichment,fill = x, color = x)) + 
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
   )+coord_flip()+ylim(c(-0.38,0.85))+
   labs(title = "TCGA LumB", y = "ssGSEA")+
   guides(fill = guide_legend(reverse = TRUE),col = guide_legend(reverse = TRUE))+
   scale_fill_manual(values = c("blue","red"))+
   scale_color_manual(values = c("blue","red"))
 
 LumBimmstat <- compare_means(ssGSEAenrichment ~ x, data = ssGSEA_LumB_GOimmune, group.by = "GOterm")
 LumBGostat <- compare_means(ssGSEAenrichment ~ x, data = ssGSEA_LumB_C1020, group.by = "GOterm")
 LumBstats <- rbind(LumBimmstat, LumBGostat)
 write.csv(LumBstats, file = "LumBwilcoxstat.csv", quote = F)
 
 #11/25
 pdf("1BCraw_LumB.pdf",width = 12,  height = 6)
 ggarrange(b,c, nrow = 2, common.legend = T, align = "v", legend = "right", heights = c(0.8,1))
 
 dev.off()
 
 
 
 #HER2 
 Her2_TCGA <- Her2_TCGA[!duplicated(Her2_TCGA$TCGA),]
 Her2 <- data.frame(Her2_TCGA, row.names = 1)
 Her2 <- Her2[-1,]
 Her2mat <- apply(Her2, 2, as.numeric)
 rownames(Her2mat) <- rownames(Her2)
 GOimmune <- gmt2list("~/Downloads/immune_10.18.gmt")
 ssGSEA_Her2_GOimmune <- gsva(expr = Her2mat, gset.idx.list = GOimmune,
                              method = "ssgsea")
 ssGSEA_Her2_GOimmune <- t(ssGSEA_Her2_GOimmune)
 Her2_immunesubsets <- read.csv("~/Her2cat.csv")
 ssGSEA_Her2_GOimmune <- merge(ssGSEA_Her2_GOimmune,Her2_immunesubsets, by.x =0, by.y = "X")
 library(reshape2)
 ssGSEA_Her2_GOimmune <- subset(ssGSEA_Her2_GOimmune, select = -c(GOBP_IMMUNE_RESPONSE_TO_TUMOR_CELL,                                         
                                                                  GOBP_LYMPHOCYTE_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE,                     
                                                                  GOBP_REGULATION_OF_T_CELL_MEDIATED_IMMUNE_RESPONSE_TO_TUMOR_CELL,           
                                                                  GOBP_T_CELL_MEDIATED_IMMUNE_RESPONSE_TO_TUMOR_CELL  ))
 
 ssGSEA_Her2_GOimmune <-  melt(ssGSEA_Her2_GOimmune, measure.vars = colnames(ssGSEA_Her2_GOimmune)[2:7],
                               variable.name = "GOterm", value.name = "ssGSEAenrichment")
 #ssGSEA_Her2_GOimmune <- subset(ssGSEA_Her2_GOimmune, category != "outlier")
 View(compare_means(ssGSEAenrichment ~ x, data = ssGSEA_Her2_GOimmune, group.by = "GOterm"))
 
 ssGSEA_Her2_GOimmune$GOterm <- factor(ssGSEA_Her2_GOimmune$GOterm,levels = rev(unique(ssGSEA_Her2_GOimmune$GOterm)))
 ssGSEA_Her2_GOimmune$x <- factor(ssGSEA_Her2_GOimmune$x, levels = unique(ssGSEA_Her2_GOimmune$x))
 
 
 b <- ggplot(ssGSEA_Her2_GOimmune, aes(x=GOterm, y=ssGSEAenrichment,fill = x, color = x)) + 
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
   )+coord_flip()+ylim(c(-0.36,0.75))+
   labs(title = "TCGA Her2", y = "ssGSEA")+
   guides(fill = guide_legend(reverse = TRUE),col = guide_legend(reverse = TRUE))+
   scale_fill_manual(values = c("blue","red"))+
   scale_color_manual(values = c("blue","red"))
 
 
 #10/20
 C1020 <- gmt2list("~/Downloads/10.20.2023_1cGeneset.gmt")
 ssGSEA_Her2_C1020 <- gsva(expr = Her2mat, gset.idx.list = C1020,
                           method = "ssgsea")
 ssGSEA_Her2_C1020 <- t(ssGSEA_Her2_C1020)
 #Her2_immunesubsets <- read.csv("~/BLBC_immunesubsets.csv")
 ssGSEA_Her2_C1020 <- merge(ssGSEA_Her2_C1020,Her2_immunesubsets, by.x =0, by.y = "X")
 
 ssGSEA_Her2_C1020 <-  melt(ssGSEA_Her2_C1020, measure.vars = colnames(ssGSEA_Her2_C1020)[2:9],
                            variable.name = "GOterm", value.name = "ssGSEAenrichment")
 #ssGSEA_Her2_C1020 <- subset(ssGSEA_Her2_C1020, category != "outlier")
 View(compare_means(ssGSEAenrichment ~ x, data = ssGSEA_Her2_C1020, group.by = "GOterm"))
 
 ssGSEA_Her2_C1020$GOterm <- factor(ssGSEA_Her2_C1020$GOterm,levels = rev(unique(ssGSEA_Her2_C1020$GOterm)))
 ssGSEA_Her2_C1020$x <- factor(ssGSEA_Her2_C1020$x, levels = unique(ssGSEA_Her2_C1020$x))
 
 c <- ggplot(ssGSEA_Her2_C1020, aes(x=GOterm, y=ssGSEAenrichment,fill = x, color = x)) + 
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
   )+coord_flip()+ylim(c(-0.36,0.75))+
   labs(title = "TCGA Her2", y = "ssGSEA")+
   guides(fill = guide_legend(reverse = TRUE),col = guide_legend(reverse = TRUE))+
   scale_fill_manual(values = c("blue","red"))+
   scale_color_manual(values = c("blue","red"))
 
 Her2immstat <- compare_means(ssGSEAenrichment ~ x, data = ssGSEA_Her2_GOimmune, group.by = "GOterm")
 Her2Gostat <- compare_means(ssGSEAenrichment ~ x, data = ssGSEA_Her2_C1020, group.by = "GOterm")
 Her2stats <- rbind(Her2immstat, Her2Gostat)
 write.csv(Her2stats, file = "Her2wilcoxstat.csv", quote = F)
 
 #11/25
 pdf("1BCraw_Her2.pdf",width = 12,  height = 6)
 ggarrange(b,c, nrow = 2, common.legend = T, align = "v", legend = "right", heights = c(0.8,1))
 
 dev.off()
 
 #normies, GET OUT!
 Normal_TCGA <- Normal_TCGA[!duplicated(Normal_TCGA$TCGA),]
 Normal <- data.frame(Normal_TCGA, row.names = 1)
 Normal <- Normal[-1,]
 Normalmat <- apply(Normal, 2, as.numeric)
 rownames(Normalmat) <- rownames(Normal)
 GOimmune <- gmt2list("~/Downloads/immune_10.18.gmt")
 ssGSEA_Normal_GOimmune <- gsva(expr = Normalmat, gset.idx.list = GOimmune,
                              method = "ssgsea")
 ssGSEA_Normal_GOimmune <- t(ssGSEA_Normal_GOimmune)
 Normal_immunesubsets <- read.csv("~/Normalcat.csv")
 ssGSEA_Normal_GOimmune <- merge(ssGSEA_Normal_GOimmune,Normal_immunesubsets, by.x =0, by.y = "X")
 library(reshape2)
 ssGSEA_Normal_GOimmune <- subset(ssGSEA_Normal_GOimmune, select = -c(GOBP_IMMUNE_RESPONSE_TO_TUMOR_CELL,                                         
                                                                  GOBP_LYMPHOCYTE_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE,                     
                                                                  GOBP_REGULATION_OF_T_CELL_MEDIATED_IMMUNE_RESPONSE_TO_TUMOR_CELL,           
                                                                  GOBP_T_CELL_MEDIATED_IMMUNE_RESPONSE_TO_TUMOR_CELL  ))
 
 ssGSEA_Normal_GOimmune <-  melt(ssGSEA_Normal_GOimmune, measure.vars = colnames(ssGSEA_Normal_GOimmune)[2:7],
                               variable.name = "GOterm", value.name = "ssGSEAenrichment")
 #ssGSEA_Normal_GOimmune <- subset(ssGSEA_Normal_GOimmune, category != "outlier")
 View(compare_means(ssGSEAenrichment ~ x, data = ssGSEA_Normal_GOimmune, group.by = "GOterm"))
 
 ssGSEA_Normal_GOimmune$GOterm <- factor(ssGSEA_Normal_GOimmune$GOterm,levels = rev(unique(ssGSEA_Normal_GOimmune$GOterm)))
 ssGSEA_Normal_GOimmune$x <- factor(ssGSEA_Normal_GOimmune$x, levels = unique(ssGSEA_Normal_GOimmune$x))
 
 
 b <- ggplot(ssGSEA_Normal_GOimmune, aes(x=GOterm, y=ssGSEAenrichment,fill = x, color = x)) + 
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
   )+coord_flip()+ylim(c(-0.36,0.9))+
   labs(title = "TCGA Normal", y = "ssGSEA")+
   guides(fill = guide_legend(reverse = TRUE),col = guide_legend(reverse = TRUE))+
   scale_fill_manual(values = c("blue","red"))+
   scale_color_manual(values = c("blue","red"))
 
 
 #10/20
 C1020 <- gmt2list("~/Downloads/10.20.2023_1cGeneset.gmt")
 ssGSEA_Normal_C1020 <- gsva(expr = Normalmat, gset.idx.list = C1020,
                           method = "ssgsea")
 ssGSEA_Normal_C1020 <- t(ssGSEA_Normal_C1020)
 #Normal_immunesubsets <- read.csv("~/BLBC_immunesubsets.csv")
 ssGSEA_Normal_C1020 <- merge(ssGSEA_Normal_C1020,Normal_immunesubsets, by.x =0, by.y = "X")
 
 ssGSEA_Normal_C1020 <-  melt(ssGSEA_Normal_C1020, measure.vars = colnames(ssGSEA_Normal_C1020)[2:9],
                            variable.name = "GOterm", value.name = "ssGSEAenrichment")
 #ssGSEA_Normal_C1020 <- subset(ssGSEA_Normal_C1020, category != "outlier")
 View(compare_means(ssGSEAenrichment ~ x, data = ssGSEA_Normal_C1020, group.by = "GOterm"))
 
 ssGSEA_Normal_C1020$GOterm <- factor(ssGSEA_Normal_C1020$GOterm,levels = rev(unique(ssGSEA_Normal_C1020$GOterm)))
 ssGSEA_Normal_C1020$x <- factor(ssGSEA_Normal_C1020$x, levels = unique(ssGSEA_Normal_C1020$x))
 
 c <- ggplot(ssGSEA_Normal_C1020, aes(x=GOterm, y=ssGSEAenrichment,fill = x, color = x)) + 
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
   )+coord_flip()+ylim(c(-0.36,0.9))+
   labs(title = "TCGA Normal", y = "ssGSEA")+
   guides(fill = guide_legend(reverse = TRUE),col = guide_legend(reverse = TRUE))+
   scale_fill_manual(values = c("blue","red"))+
   scale_color_manual(values = c("blue","red"))
 
 
 Normalimmstat <- compare_means(ssGSEAenrichment ~ x, data = ssGSEA_Normal_GOimmune, group.by = "GOterm")
 NormalGostat <- compare_means(ssGSEAenrichment ~ x, data = ssGSEA_Normal_C1020, group.by = "GOterm")
 Normalstats <- rbind(Normalimmstat, NormalGostat)
 write.csv(Normalstats, file = "Normalwilcoxstat.csv", quote = F)
 
 
 #11/25
 pdf("1BCraw_Normal.pdf",width = 12,  height = 6)
 ggarrange(b,c, nrow = 2, common.legend = T, align = "v", legend = "right", heights = c(0.8,1))
 
 dev.off()
 
