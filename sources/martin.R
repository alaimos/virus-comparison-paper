library(readr)
library(readxl)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
genes <- read_excel("../external/1-s2.0-S1097276521003130-mmc3.xlsx", 
                    sheet = "Overexpression screen", skip = 4)

colnames(genes)[1:2] <- c("gene_name", "gene_id")
genes$gene_id        <- as.character(genes$gene_id)
genes                <- genes[,1:2]

cov2  <- read_tsv("../simulations/A549_COV2_MOI_2.tsv")[,c(3,7)]
hpiv3 <- read_tsv("../simulations/A549_HPIV3.tsv")[,c(3,7)]
hrv   <- read_tsv("../simulations/A549_HRV_24hpi.tsv")[,c(3,7)]
iav   <- read_tsv("../simulations/A549_IAV.tsv")[,c(3,7)]
rsv   <- read_tsv("../simulations/A549_RSV.tsv")[,c(3,7)]
colnames(cov2)  <- c("gene_id", "activity")
colnames(hpiv3) <- c("gene_id", "activity")
colnames(hrv)   <- c("gene_id", "activity")
colnames(iav)   <- c("gene_id", "activity")
colnames(rsv)   <- c("gene_id", "activity")

cov2_ne  <- read_tsv("../simulations/A549_COV2_MOI_2_non_exp.txt", col_names = FALSE)
hpiv3_ne <- read_tsv("../simulations/A549_HPIV3_non_exp.txt", col_names = FALSE)
hrv_ne   <- read_tsv("../simulations/A549_HRV_24hpi_non_exp.txt", col_names = FALSE)
iav_ne   <- read_tsv("../simulations/A549_IAV_non_exp.txt", col_names = FALSE)
rsv_ne   <- read_tsv("../simulations/A549_RSV_non_exp.txt", col_names = FALSE)
colnames(cov2_ne)  <- c("gene_id")
colnames(hpiv3_ne) <- c("gene_id")
colnames(hrv_ne)   <- c("gene_id")
colnames(iav_ne)   <- c("gene_id")
colnames(rsv_ne)   <- c("gene_id")
cov2_ne$gene_id  <- as.character(cov2_ne$gene_id)
hpiv3_ne$gene_id <- as.character(hpiv3_ne$gene_id)
hrv_ne$gene_id   <- as.character(hrv_ne$gene_id)
iav_ne$gene_id   <- as.character(iav_ne$gene_id)
rsv_ne$gene_id   <- as.character(rsv_ne$gene_id)
cov2_ne$NE  <- TRUE
hpiv3_ne$NE <- TRUE
hrv_ne$NE   <- TRUE
iav_ne$NE   <- TRUE
rsv_ne$NE   <- TRUE

cov2  <- cov2 %>% inner_join(genes) %>% unique() %>% left_join(cov2_ne) %>% mutate(NE = ifelse(is.na(NE), FALSE, NE))
hpiv3 <- hpiv3 %>% inner_join(genes) %>% unique() %>% left_join(hpiv3_ne) %>% mutate(NE = ifelse(is.na(NE), FALSE, NE))
hrv   <- hrv %>% inner_join(genes) %>% unique() %>% left_join(hrv_ne) %>% mutate(NE = ifelse(is.na(NE), FALSE, NE)) 
iav   <- iav %>% inner_join(genes) %>% unique() %>% left_join(iav_ne) %>% mutate(NE = ifelse(is.na(NE), FALSE, NE))
rsv   <- rsv %>% inner_join(genes) %>% unique() %>% left_join(rsv_ne) %>% mutate(NE = ifelse(is.na(NE), FALSE, NE))
colnames(cov2)[2:4]  <- c("SARS-CoV2 MOI 2.0", "gene_name", "NE_cov2")
colnames(hpiv3)[2:4] <- c("HPIV3", "gene_name", "NE_hpiv3")
colnames(hrv)[2:4]   <- c("Rhinovirus", "gene_name", "NE_hrv")
colnames(iav)[2:4]   <- c("IAV", "gene_name", "NE_iav")
colnames(rsv)[2:4]   <- c("RSV", "gene_name", "NE_rsv")

data <- cov2 %>% 
  inner_join(hpiv3) %>% 
  inner_join(hrv) %>% 
  inner_join(iav) %>% 
  inner_join(rsv) %>%
  select(gene_name, `SARS-CoV2 MOI 2.0`, HPIV3, Rhinovirus, IAV, RSV)
data_m <- data.matrix(data[,-1])
rownames(data_m) <- data$gene_name
missing <- genes$gene_name[!(genes$gene_name %in% rownames(data_m))]
missing_m <- matrix(NA, nrow = length(missing), ncol = ncol(data_m))
colnames(missing_m) <- colnames(data_m)
rownames(missing_m) <- missing
data_m <- rbind(data_m, missing_m)
pdf("../figures/martin_sancho.pdf", height = 15)
col_fun = colorRamp2(c(-6, 0, 6), c("blue", "white", "red"))
ht <- Heatmap(data_m, row_names_side = "left", row_dend_side = "right", 
              column_names_side = "top", column_dend_side = "bottom", 
              border = TRUE, row_labels = rownames(data_m), cluster_columns = TRUE, 
              name="Activity", width = unit(3, "cm"), height = unit(25, "cm"), 
              col=col_fun, use_raster = TRUE, raster_device = "png",
              clustering_distance_rows = function (m) {
                tmp <- dist(m, method = "euclidean")
                tmp[is.na(tmp)] <- max(tmp, na.rm = TRUE) + 1
                return(tmp)
              })
draw(ht)
dev.off()




