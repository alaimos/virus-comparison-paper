library(readxl)
ISG15 <- read_excel("../excels/ISG15.xlsx")
classes <- ISG15$`Node Class`
class(ISG15) <- "data.frame"
rownames(ISG15) <- make.names(ISG15$`Node Name`, unique = TRUE)
orig.names <- ISG15$`Node Name`
ISG15 <- ISG15[,-(1:2)]
ISG15 <- data.matrix(ISG15)

library(ComplexHeatmap)
library(circlize)
s <- classes != "JNK-AP1"
col_fun = colorRamp2(c(-6, 0, 6), c("blue", "white", "red"))
pdf("../figures/ISG15.pdf", height = 20)
ht <- Heatmap(ISG15[s,c(1,3,2,4,5)], row_names_side = "left", row_dend_side = "right", 
              column_names_side = "top", column_dend_side = "bottom", 
              row_split = classes[s], row_gap = unit(0, "mm"), border = TRUE, 
              row_title_rot = 0,
              row_labels = orig.names[s], cluster_columns = TRUE, name="Activity", 
              width = unit(3, "cm"), height = unit(20, "cm"), col=col_fun, 
              use_raster = TRUE, raster_device = "png")
draw(ht)
dev.off()

IFNg <- read_excel("../excels/IFNgamma.xlsx")
class(IFNg) <- "data.frame"
rownames(IFNg) <- make.names(IFNg$`Node Name`, unique = TRUE)
orig.names.ifn <- IFNg$`Node Name`
IFNg <- IFNg[,-1]
IFNg <- data.matrix(IFNg)

col_fun = colorRamp2(c(-6, 0, 6), c("blue", "white", "red"))
pdf("../figures/IFNg.pdf", height = 10)
ht <-Heatmap(IFNg[,c(1,3,2,4,5)], row_names_side = "left", row_dend_side = "right", 
             column_names_side = "top", column_dend_side = "bottom", 
             border = TRUE, row_labels = orig.names.ifn, cluster_columns = TRUE, 
             name="Activity", 
             width = unit(3, "cm"), height = unit(7, "cm"), col=col_fun, 
             use_raster = TRUE, raster_device = "png")
draw(ht)
dev.off()


library(readxl)
innate_ifn <- read_excel("../excels/virus_activity.xlsx", sheet = "Innate_IFN")
classes <- innate_ifn$pathway
class(innate_ifn) <- "data.frame"
rownames(innate_ifn) <- make.names(innate_ifn$gene, unique = TRUE)
orig.names <- sapply(strsplit(innate_ifn$gene, " "), function(x)(x[1]))
innate_ifn <- innate_ifn[,-(1:2)]
innate_ifn <- data.matrix(innate_ifn)

library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-6, 0, 6), c("blue", "white", "red"))
pdf("../figures/innate_inf.pdf", height = 10)
ht <- Heatmap(innate_ifn[,c(1,3,2,4,5)], row_names_side = "left", row_dend_side = "right", 
              column_names_side = "top", column_dend_side = "bottom", 
              row_split = classes, row_gap = unit(0, "mm"), border = TRUE, 
              row_labels = orig.names, cluster_columns = TRUE, name="Activity", 
              width = unit(3, "cm"), height = unit(20, "cm"), col=col_fun, 
              use_raster = TRUE, raster_device = "png")
draw(ht)
dev.off()


library(readxl)
nucleo <- read_excel("../excels/virus_activity.xlsx", sheet = "Nucleoporins")
class(nucleo) <- "data.frame"
position <- nucleo$position
rownames(nucleo) <- make.names(nucleo$gene, unique = TRUE)
orig.names <- sapply(strsplit(nucleo$gene, " "), function(x)(x[1]))
nucleo <- nucleo[,-(1:2)]
nucleo <- data.matrix(nucleo)

library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-6, 0, 6), c("blue", "white", "red"))
pdf("../figures/nucleoporins.pdf", height = 10)
ht <- Heatmap(nucleo[,c(1,3,2,4,5)], row_names_side = "left", row_dend_side = "right", 
              column_names_side = "top", column_dend_side = "bottom", 
              row_split = position, row_gap = unit(0, "mm"), row_title_rot = 0, 
              border = TRUE, row_labels = orig.names, cluster_columns = TRUE, 
              name="Activity", width = unit(3, "cm"), height = unit(20, "cm"), 
              col=col_fun, use_raster = TRUE, raster_device = "png")
draw(ht)
dev.off()

library(readxl)
muci <- read_excel("../excels/virus_activity.xlsx", sheet = "Mucins")
class(muci) <- "data.frame"
rownames(muci) <- make.names(muci$gene, unique = TRUE)
orig.names <- sapply(strsplit(muci$gene, " "), function(x)(x[1]))
muci <- muci[,-1]
muci <- data.matrix(muci)

library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-6, 0, 6), c("blue", "white", "red"))
pdf("../figures/mucins.pdf", height = 10)
ht <- Heatmap(muci[,c(1,3,2,4,5)], row_names_side = "left", row_dend_side = "right", 
              column_names_side = "top", column_dend_side = "bottom", 
              border = TRUE, row_labels = orig.names, cluster_columns = TRUE, 
              name="Activity", width = unit(3, "cm"), height = unit(20, "cm"), 
              col=col_fun, use_raster = TRUE, raster_device = "png")
draw(ht)
dev.off()
