library(VennDiagram)
cov2  <- readr::read_tsv("../simulations/A549_COV2_MOI_2.tsv")
hpiv3 <- readr::read_tsv("../simulations/A549_HPIV3.tsv")
hrv   <- readr::read_tsv("../simulations/A549_HRV_24hpi.tsv")
iav   <- readr::read_tsv("../simulations/A549_IAV.tsv")
rsv   <- readr::read_tsv("../simulations/A549_RSV.tsv")

get.pathways <- function (d, up) {
  dd <- unique(d[,c(1, 11)])
  print(nrow(dd))
  if (up) {
    dd <- dd[[1]][dd[[2]] > 0]
  } else {
    dd <- dd[[1]][dd[[2]] < 0]
  }
  return (dd)
}

get.genes <- function (d, up) {
  dd <- unique(d[,c(3, 7)])
  print(nrow(dd))
  if (up) {
    dd <- dd[[1]][dd[[2]] > 0]
  } else {
    dd <- dd[[1]][dd[[2]] < 0]
  }
  return (dd)
}

cov2_p_up <- get.pathways(cov2, TRUE)
hpiv3_p_up <- get.pathways(hpiv3, TRUE)
hrv_p_up <- get.pathways(hrv, TRUE)
iav_p_up <- get.pathways(iav, TRUE)
rsv_p_up <- get.pathways(rsv, TRUE)

cov2_p_down <- get.pathways(cov2, FALSE)
hpiv3_p_down <- get.pathways(hpiv3, FALSE)
hrv_p_down <- get.pathways(hrv, FALSE)
iav_p_down <- get.pathways(iav, FALSE)
rsv_p_down <- get.pathways(rsv, FALSE)

cov2_g_up <- get.genes(cov2, TRUE)
hpiv3_g_up <- get.genes(hpiv3, TRUE)
hrv_g_up <- get.genes(hrv, TRUE)
iav_g_up <- get.genes(iav, TRUE)
rsv_g_up <- get.genes(rsv, TRUE)

cov2_g_down <- get.genes(cov2, FALSE)
hpiv3_g_down <- get.genes(hpiv3, FALSE)
hrv_g_down <- get.genes(hrv, FALSE)
iav_g_down <- get.genes(iav, FALSE)
rsv_g_down <- get.genes(rsv, FALSE)

p_names <- cov2[,c(1:2)] %>% unique()
colnames(p_names) <- c("id", "name")

g_names <- cov2[,c(3,4)] %>% unique()
colnames(g_names) <- c("id", "name")

write_overlap <- function (p, p_names, fname) {
  nms <- p[[1]]
  for (x in p) {
    nms <- union(nms, x)
  }
  df <- data.frame(id=nms, matrix("No", nrow=length(nms), ncol=length(p)))
  colnames(df)[2:(length(p)+1)] <- names(p)
  for (c in seq_along(p)) {
    df[df$id %in% p[[c]],c+1] <- "Yes"
  }
  df <- df %>% inner_join(p_names)
  df <- df[,c(1, ncol(df), 2:(ncol(df)-1))]
  openxlsx::write.xlsx(df, fname, asTable=TRUE, rowNames=FALSE)
}

p <- list(
  "SARS-CoV2"=cov2_p_up,
  "HPIV3"=hpiv3_p_up,
  "HRV"=hrv_p_up,
  "IAV"=iav_p_up,
  "RSV"=rsv_p_up
)
venn.diagram(p, filename = "../figures/pathway_comparison_upregulated.png", na = "remove", resolution = 300, width = 10, height = 10, units = "in", disable.logging = TRUE, cat.pos=c(0, 0, 140, 100, 0), fill = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#9E2A2B"))
write_overlap(p, p_names, "../tables/pathway_comparison_upregulated.xlsx")

p <- list(
  "SARS-CoV2"=cov2_p_down,
  "HPIV3"=hpiv3_p_down,
  "HRV"=hrv_p_down,
  "IAV"=iav_p_down,
  "RSV"=rsv_p_down
)
venn.diagram(p, filename = "../figures/pathway_comparison_downregulated.png", na = "remove", resolution = 300, width = 10, height = 10, units = "in", disable.logging = TRUE, cat.pos=c(0, 0, 140, 100, 0), fill = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#9E2A2B"))
write_overlap(p, p_names, "../tables/pathway_comparison_downregulated.xlsx")

p <- list(
  "SARS-CoV2"=cov2_g_up,
  "HPIV3"=hpiv3_g_up,
  "HRV"=hrv_g_up,
  "IAV"=iav_g_up,
  "RSV"=rsv_g_up
)
venn.diagram(p, filename = "../figures/genes_comparison_upregulated.png", na = "remove", resolution = 300, width = 10, height = 10, units = "in", disable.logging = TRUE, cat.pos=c(0, 0, 140, 100, 0), fill = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#9E2A2B"))
write_overlap(p, g_names, "../tables/genes_comparison_upregulated.xlsx")

p <- list(
  "SARS-CoV2"=cov2_g_down,
  "HPIV3"=hpiv3_g_down,
  "HRV"=hrv_g_down,
  "IAV"=iav_g_down,
  "RSV"=rsv_g_down
)
venn.diagram(p, filename = "../figures/genes_comparison_downregulated.png", na = "remove", resolution = 300, width = 10, height = 10, units = "in", disable.logging = TRUE, cat.pos=c(0, 0, 140, 100, 0), fill = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#9E2A2B"))
write_overlap(p, g_names, "../tables/genes_comparison_downregulated.xlsx")


get.activities <- function (d, sel=NULL) {
  dd <- unique(d[,c(4, 7)])
  if (!is.null(sel)) {
    dd <- dd[dd[[1]] %in% sel,]
  }
  print(nrow(dd))
  return (dd)
}

pol3_genes_o <- c("RPAC1","RPAC2","RPC1","RPC2","RPC3","RPC4","RPC5","RPC6",
                "RPC7","RPC8","RPC9","RPC10","RPABC1","RPABC2","RPABC3",
                "RPABC4","RPABC5")
pol3_genes <- c("POLR1C","POLR1D","POLR3A","POLR3B","POLR3C","POLR3D","POLR3E",
                "POLR3F","POLR3G","POLR3H","CRCP","POLR3K","POLR2E","POLR2F",
                "POLR2H","POLR2K","POLR2L")
cov2_a  <- get.activities(cov2, pol3_genes)
hpiv3_a <- get.activities(hpiv3, pol3_genes)
hrv_a   <- get.activities(hrv, pol3_genes)
iav_a   <- get.activities(iav, pol3_genes)
rsv_a   <- get.activities(rsv, pol3_genes)

tmp <- cbind(cov2_a, hpiv3_a, hrv_a, iav_a, rsv_a)
rownames(tmp) <- tmp[,1]
tmp <- tmp[,c(2,4,6,8,10)]
colnames(tmp) <- c("SARS-CoV2 MOI 2.0", "HPIV3", "Rhinovirus", "IAV", "RSV")
tmp <- data.matrix(tmp)
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-4.83, 0, 4.83), c("blue", "white", "red"))
pdf("../figures/pol3.pdf", height = 10)
ht <- Heatmap(tmp, row_names_side = "left", row_dend_side = "right", 
              column_names_side = "top", column_dend_side = "bottom", 
              border = TRUE, row_labels = rownames(tmp), cluster_columns = TRUE, 
              name="Activity", width = unit(3, "cm"), height = unit(20, "cm"), 
              col=col_fun, use_raster = TRUE, raster_device = "png")
draw(ht)
dev.off()
