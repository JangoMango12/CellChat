library(arrow)
pos <- read_parquet("tissue_positions.parquet")
names(pos)
sapply(pos[sapply(pos, is.numeric)], range, na.rm = TRUE)

install.packages("arrow")
getwd()
list.files(pattern = "parquet")
setwd("C:/Users/<you>/MAPK_spatial")
list.files()
setwd("C:/Users/write/OneDrive/Documents/MapKSeries_Segmented_MouseDuodenum")
list.files()
library(arrow)
pos <- read_parquet("tissue_positions.parquet")
names(pos)
sapply(pos[sapply(pos, is.numeric)], range, na.rm = TRUE)
library(rhdf5)
h5ls("filtered_feature_bc_matrix.h5")
library(rhdf5)
h5ls("filtered_feature_bc_matrix.h5")
library(rhdf5)
h5ls("filtered_feature_bc_matrix.h5")
library(Seurat)
mat <- Read10X_h5("filtered_feature_bc_matrix.h5")
head(colnames(mat))
install.packages("Seurat")
library(Seurat)
mat <- Read10X_h5("filtered_feature_bc_matrix.h5")
lab <- read.csv("Graph-Based.csv")

mean(lab$Barcode %in% colnames(mat))

mat <- Read10X_h5("filtered_feature_bc_matrix.h5")
lab <- read.csv("Graph-Based.csv")

mean(lab$Barcode %in% colnames(mat))
library(jsonlite)
fromJSON("scalefactors_json.json")library(arrow)
pos <- read_parquet("tissue_positions.parquet")
names(pos)
head(pos$barcode, 3)
library(arrow)
pos <- read_parquet("tissue_positions.parquet")
names(pos)
head(pos$barcode, 3)
library(Seurat)
mat2 <- Read10X_h5("filtered_feature_cell_matrix.h5")
dim(mat2)
head(colnames(mat2), 3)

lab <- read.csv("Graph-Based.csv")
mean(lab$Barcode %in% colnames(mat2))
proj <- read.csv("Default Pipeline Analysis-Spatial-Projection.csv")
gen  <- read.csv("Genotypes.csv")

names(proj)
head(proj, 3)
sapply(proj[sapply(proj, is.numeric)], range)
mean(proj$Barcode %in% colnames(mat2))

names(gen)
head(gen, 3)
mean(gen$Barcode %in% colnames(mat2))

library(Seurat); library(dplyr)

names(lab)
obj$celltype <- lab$Graph.based[match(colnames(obj), lab$Barcode)]

table(obj$celltype, useNA = "ifany")
table(obj$genotype, useNA = "ifany")
summary(obj$nCount_Spatial)

epi <- read.csv("ControlEpithelialCells.csv")
names(epi)
head(epi, 3)
nrow(epi)

library(Seurat)
set.seed(42)

# load full matrix, subset to your control epithelial cells
mat <- Read10X_h5("filtered_feature_cell_matrix.h5")
epi_bc <- epi$Barcode
epi_bc <- epi_bc[epi_bc %in% colnames(mat)]
cat("matched", length(epi_bc), "of", nrow(epi), "barcodes\n")   # want ~17777

# build object on just those cells
e <- CreateSeuratObject(mat[, epi_bc], assay = "Spatial")
e <- NormalizeData(e)
e <- FindVariableFeatures(e, nfeatures = 2000)
e <- ScaleData(e)
e <- RunPCA(e, npcs = 20)

# k-means = 5 on the PCA embedding (the thing Loupe wouldn't do)
km <- kmeans(Embeddings(e, "pca")[, 1:15], centers = 5, nstart = 25, iter.max = 100)
e$zone_kmeans <- factor(km$cluster)

table(e$zone_kmeans)

# write barcode -> cluster for Loupe import
out <- data.frame(Barcode = colnames(e), Zone = paste0("kZone_", km$cluster))
write.csv(out, "epi_kmeans5.csv", row.names = FALSE)
cat("wrote epi_kmeans5.csv\n")

proj <- read.csv("Default Pipeline Analysis-Spatial-Projection.csv")
km_out <- read.csv("epi_kmeans5.csv")

proj$Zone <- km_out$Zone[match(proj$Barcode, km_out$Barcode)]
proj$Zone[is.na(proj$Zone)] <- "non_epithelial"   # cells not in the epi subset

write.csv(proj, "projection_with_zones.csv", row.names = FALSE)
cat("wrote projection_with_zones.csv\n")

set.seed(42)
km6 <- kmeans(Embeddings(e, "pca")[, 1:15], centers = 6, nstart = 25, iter.max = 100)
e$zone_k6 <- factor(km6$cluster)
table(e$zone_k6)

out6 <- data.frame(Barcode = colnames(e), Zone = paste0("k6_", km6$cluster))
write.csv(out6, "epi_kmeans6.csv", row.names = FALSE)

library(dplyr)
q@meta.data %>%
  filter(genotype == "Villin CreERT2; BrafV600E/V600E") %>%
  group_by(zone_named) %>%
  summarise(median_score = median(score), n = n())

# control cells are in object e, with their scores? no - e was the reference, no scores.
# instead compare mutant villus to a genotype that should be normal, e.g. Erk DKO only:
library(dplyr)
q@meta.data %>%
  filter(genotype %in% c("Villin CreERT2; Erk1/Erk2 DKO",
                         "Villin CreERT2; BrafV600E/V600E")) %>%
  group_by(genotype, zone_named) %>%
  summarise(median_score = median(score), n = n(), .groups="drop")

sc <- read.csv("Stromal Clusters Features.csv")
rownames(sc) <- make.unique(sc$FeatureName)

# marker panel (cited sources + canonical immune)
panel <- list(
  Macrophage      = c("C1qa","C1qb","C1qc","Csf1r","Adgre1","Lyz2","Fcgr1"),
  Mac_SPP1        = c("Spp1"),
  Mac_SELENOP     = c("Selenop"),
  DC              = c("Itgax","Flt3","Xcr1","Batf3"),
  T_cell          = c("Cd3e","Cd3d","Cd3g"),
  CD8_T           = c("Cd8a","Cd8b1"),
  CD4_T           = c("Cd4"),
  Treg            = c("Foxp3","Il2ra"),
  B_cell          = c("Cd19","Cd79a","Ms4a1"),
  Plasma          = c("Jchain","Mzb1","Xbp1","Igha"),
  NK_ILC          = c("Ncr1","Klrb1c","Gzmb"),
  Mast            = c("Cpa3","Kit","Mcpt4","Cma1"),
  Fibro_pan       = c("Pdgfra","Col1a1","Col3a1","Dcn"),
  Telocyte        = c("Foxl1","Gli1"),
  Trophocyte      = c("Rspo3","Grem1"),
  Myofibroblast   = c("Acta2","Tagln","Myh11"),
  Pericyte        = c("Pdgfrb","Notch3","Rgs5"),
  Endothelial     = c("Pecam1","Cdh5","Vwf"),
  Lymphatic_endo  = c("Prox1","Lyve1")
)

avg_cols <- grep("Average", names(sc), value = TRUE)
clusters <- sub("\\.Average","", avg_cols)

# z-score each gene across clusters so highly-expressed genes don't dominate
M <- as.matrix(sc[, avg_cols]); colnames(M) <- clusters
M <- M[rowSums(M) > 0, ]
Z <- t(scale(t(M)))   # per-gene z across clusters
Z[is.na(Z)] <- 0

score <- sapply(panel, function(genes){
  g <- intersect(genes, rownames(Z))
  if(length(g)==0) return(rep(NA, ncol(Z)))
  colMeans(Z[g,,drop=FALSE])
})
rownames(score) <- clusters

# best call per cluster
for(cl in clusters){
  s <- sort(score[cl,], decreasing=TRUE)
  cat(sprintf("%-12s -> %-14s (%.2f) | next: %s (%.2f)\n",
              cl, names(s)[1], s[1], names(s)[2], s[2]))
}


# ============ CITED STROMAL/IMMUNE ANNOTATION ============
sc <- read.csv("Stromal Clusters Features.csv")
rownames(sc) <- make.unique(sc$FeatureName)

# panel — each entry cited in comments
panel <- list(
  Macrophage      = c("C1qa","C1qb","C1qc","Csf1r","Adgre1","Apoe","Mpeg1","Aif1"), # Oliveira25, Sathe23
  Mac_SPP1        = c("Spp1"),                          # Sathe23, Oliveira25
  Mac_SELENOP     = c("Selenop"),                       # Oliveira25
  DC              = c("Itgax","Flt3","Xcr1","Batf3"),   # canonical
  Tcell           = c("Cd3e","Cd3d","Cd3g"),            # canonical
  CD8_T           = c("Cd8a","Cd8b1"),
  CD4_T           = c("Cd4"),
  Treg            = c("Foxp3","Il2ra"),
  Cytotoxic_NK    = c("Gzma","Gzmb","Nkg7","Ncr1","Klrb1c","Klrd1"), # canonical
  Bcell           = c("Cd79a","Cd79b","Ms4a1"),
  Plasma          = c("Jchain","Mzb1","Ighm","Igha"),
  Mast            = c("Cpa3","Kit","Mcpt4","Cma1"),
  Fibro_pan       = c("Pdgfra","Col1a1","Col3a1","Dcn"),   # McCarthy20
  Telocyte        = c("Foxl1","Gli1","Bmp4","Bmp7"),       # McCarthy20 (Pdgfra-hi, BMP source)
  Trophocyte      = c("Cd81","Grem1","Rspo3"),             # McCarthy20 (Pdgfra-lo, BMP-antagonist)
  Myofibroblast   = c("Acta2","Tagln","Myh11"),
  Pericyte        = c("Pdgfrb","Notch3","Rgs5"),
  Endothelial     = c("Pecam1","Cdh5","Cldn5","Vwf"),
  Lymphatic_endo  = c("Prox1","Lyve1","Mmrn1")
)

avg_cols <- grep("Average", names(sc), value = TRUE)
clusters <- sub("\\.Average","", avg_cols)

# z-score per gene across clusters (relative specificity, not raw abundance)
M <- as.matrix(sc[, avg_cols]); colnames(M) <- clusters
M <- M[rowSums(M) > 0, ]
Z <- t(scale(t(M))); Z[is.na(Z)] <- 0

score <- sapply(panel, function(g){
  g <- intersect(g, rownames(Z))
  if(length(g)==0) return(rep(NA, ncol(Z)))
  colMeans(Z[g,,drop=FALSE])
})
rownames(score) <- clusters

# ---- OUTPUT 1: best call + margin (confidence) per cluster ----
cat("\n===== CALLS (score = relative specificity; margin = confidence) =====\n")
for(cl in clusters){
  s <- sort(score[cl,], decreasing=TRUE)
  margin <- s[1]-s[2]
  flag <- if(s[1]<0.5) "  <<LOW-no signature>>" else if(margin<0.3) "  <<AMBIGUOUS-tied>>" else ""
  cat(sprintf("%-11s -> %-14s (%.2f)  2nd:%-12s(%.2f)  margin=%.2f%s\n",
              cl, names(s)[1], s[1], names(s)[2], s[2], margin, flag))
}

# ---- OUTPUT 2: raw evidence — top expressed panel genes per cluster ----
all_panel <- unique(unlist(panel))
present <- all_panel[all_panel %in% sc$FeatureName]
tab <- sapply(avg_cols, function(a) sc[[a]][match(present, sc$FeatureName)])
rownames(tab) <- present; colnames(tab) <- clusters

cat("\n===== RAW MARKER EVIDENCE (top expressed lineage genes per cluster) =====\n")
for(cl in clusters){
  v <- sort(tab[,cl], decreasing=TRUE); v <- v[v>0.05][1:6]
  v <- v[!is.na(v)]
  cat(sprintf("%-11s: %s\n", cl, paste(sprintf("%s=%.2f",names(v),v), collapse="  ")))
}

library(Matrix)

strom <- strom[strom$Barcode %in% colnames(mat2), ]
cl_vec <- strom[[2]]                    # "Cluster 1", "Cluster 2", ...
names(cl_vec) <- strom$Barcode

sub <- mat2[, strom$Barcode]
sub <- NormalizeData(CreateSeuratObject(sub, assay="Spatial"))
expr <- GetAssayData(sub, layer="data")

clusters <- unique(cl_vec)
cluster_mean <- sapply(clusters, function(k){
  cells <- names(cl_vec)[cl_vec == k]
  Matrix::rowMeans(expr[, cells, drop=FALSE])
})
colnames(cluster_mean) <- clusters

Z <- t(scale(t(cluster_mean))); Z[is.na(Z)] <- 0

panel <- list(
  Macrophage=c("C1qa","C1qb","C1qc","Csf1r","Adgre1","Mpeg1","Aif1"),
  Mac_SPP1=c("Spp1"), Mac_SELENOP=c("Selenop"),
  DC=c("Itgax","Flt3","Xcr1","Batf3"),
  Tcell=c("Cd3e","Cd3d","Cd3g"), CD8_T=c("Cd8a","Cd8b1"), CD4_T=c("Cd4"),
  Treg=c("Foxp3","Il2ra"), Cytotoxic_NK=c("Gzma","Gzmb","Nkg7","Ncr1","Klrb1c","Klrd1"),
  Bcell=c("Cd79a","Cd79b","Ms4a1"), Plasma=c("Jchain","Mzb1","Ighm","Igha"),
  Mast=c("Cpa3","Kit","Mcpt4","Cma1"),
  Fibro_pan=c("Pdgfra","Col1a1","Col3a1","Dcn"),
  Telocyte=c("Foxl1","Gli1","Bmp4","Bmp7"),
  Trophocyte=c("Cd81","Grem1","Rspo3"),
  Myofibroblast=c("Acta2","Tagln","Myh11"),
  Pericyte=c("Pdgfrb","Notch3","Rgs5"),
  Endothelial=c("Pecam1","Cdh5","Cldn5","Vwf"),
  Lymphatic_endo=c("Prox1","Lyve1","Mmrn1")
)

score <- sapply(panel, function(g){
  g <- intersect(g, rownames(Z)); if(!length(g)) return(rep(NA,ncol(Z)))
  colMeans(Z[g,,drop=FALSE])
})
rownames(score) <- colnames(Z)

# order output Cluster 1..19
ord <- paste("Cluster", 1:19)
cat("\n===== CALLS (full matrix) =====\n")
for(cl in ord){
  s <- sort(score[cl,], decreasing=TRUE); m <- s[1]-s[2]
  flag <- if(s[1]<0.5)"  <<LOW>>" else if(m<0.3)"  <<TIED>>" else ""
  cat(sprintf("%-11s -> %-14s(%.2f) 2nd:%-12s(%.2f) margin=%.2f%s\n",
              cl,names(s)[1],s[1],names(s)[2],s[2],m,flag))
}

cat("\n===== KEY FIBROBLAST MARKERS (raw mean expr) =====\n")
fib <- c("Pdgfra","Col1a1","Col3a1","Cd81","Grem1","Foxl1","Gli1","Rspo3","Bmp4","Acta2","Tagln","Pdgfrb")
print(round(cluster_mean[intersect(fib,rownames(cluster_mean)), ord], 2))



library(Seurat)

# k-means = 10 on stromal cells
ss <- CreateSeuratObject(mat2[, strom$Barcode], assay="Spatial")
ss <- NormalizeData(ss)
ss <- FindVariableFeatures(ss, nfeatures=2000)
ss <- ScaleData(ss)
ss <- RunPCA(ss, npcs=20)

set.seed(42)
km10 <- kmeans(Embeddings(ss,"pca")[,1:15], centers=10, nstart=25, iter.max=100)
strom$kmeans <- paste0("kS_", km10$cluster)
table(strom$kmeans)

library(Seurat)

# make sure the k-means labels are on the object
ss$kmeans <- strom$kmeans[match(colnames(ss), strom$Barcode)]
Idents(ss) <- "kmeans"

# DE: each cluster vs all others (Wilcoxon, Seurat default)
markers <- FindAllMarkers(ss,
                          only.pos = TRUE,        # only genes UP in the cluster
                          min.pct = 0.10,         # expressed in >=10% of cluster cells
                          logfc.threshold = 0.25) # meaningful fold change

# top 15 markers per cluster by fold change
library(dplyr)
top <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 15) %>%
  ungroup()

# print cleanly, per cluster
for(k in sort(unique(top$cluster))){
  cat("\n=== ", k, " ===\n")
  d <- top[top$cluster==k, c("gene","avg_log2FC","pct.1","pct.2","p_val_adj")]
  print(as.data.frame(d), row.names=FALSE)
}

# also save to CSV so you can scan it / put in slides
write.csv(markers, "stromal_kmeans_markers.csv", row.names=FALSE)













library(Seurat)
ss <- CreateSeuratObject(mat2[, strom$Barcode], assay="Spatial")
ss <- NormalizeData(ss)
ss$cl <- strom$Stromal.Clusters[match(colnames(ss), strom$Barcode)]
Idents(ss) <- "cl"

genes <- c("Ptprc",              # immune (CD45) - splits immune vs stroma
           "C1qa","Csf1r","Lyz2", # MACROPHAGE
           "Cd3e","Cd3d",         # T CELL
           "Cd79a","Ms4a1",       # B CELL
           "Jchain","Mzb1",       # PLASMA
           "Cpa3","Kit",          # MAST
           "Itgax","Flt3",        # DC
           "Col1a1","Col3a1","Dcn", # FIBROBLAST
           "Acta2","Myh11",       # SMOOTH MUSCLE/MYOFIBRO
           "Pecam1","Cldn5",      # ENDOTHELIAL
           "Lyve1")               # LYMPHATIC
avg <- AverageExpression(ss, features=genes, assays="Spatial")$Spatial
round







































library(Seurat)
set.seed(42)

# k-means = 5 on stromal cells (same as epithelial)
ss <- CreateSeuratObject(mat2[, strom$Barcode], assay="Spatial")
ss <- NormalizeData(ss)
ss <- FindVariableFeatures(ss, nfeatures=2000)
ss <- ScaleData(ss)
ss <- RunPCA(ss, npcs=20)

km5 <- kmeans(Embeddings(ss,"pca")[,1:15], centers=5, nstart=25, iter.max=100)
ss$kmeans5 <- factor(km5$cluster)
table(ss$kmeans5)

# write barcode -> cluster
out <- data.frame(Barcode = colnames(ss), Zone = paste0("kStroma_", km5$cluster))
write.csv(out, "stroma_kmeans5.csv", row.names=FALSE)

# merge into the projection file for Loupe (same as epithelial)
proj <- read.csv("Default Pipeline Analysis-Spatial-Projection.csv")
km_out <- read.csv("stroma_kmeans5.csv")
proj$StromaZone <- km_out$Zone[match(proj$Barcode, km_out$Barcode)]
proj$StromaZone[is.na(proj$StromaZone)] <- "non_stromal"
write.csv(proj, "projection_with_stroma.csv", row.names=FALSE)
cat("wrote projection_with_stroma.csv\n")
















set.seed(42)
km8 <- kmeans(Embeddings(ss,"pca")[,1:15], centers=8, nstart=25, iter.max=100)
ss$k8 <- factor(km8$cluster)
table(ss$k8)

# export for loupe
out <- data.frame(Barcode = colnames(ss), Zone = paste0("kStroma8_", km8$cluster))
write.csv(out, "stroma_kmeans8.csv", row.names=FALSE)
proj <- read.csv("Default Pipeline Analysis-Spatial-Projection.csv")
proj$Stroma8 <- out$Zone[match(proj$Barcode, out$Barcode)]
proj$Stroma8[is.na(proj$Stroma8)] <- "non_stromal"
write.csv(proj, "projection_with_stroma8.csv", row.names=FALSE)

# markers per cluster
Idents(ss) <- "k8"
genes <- c("C1qa","Csf1r","Lyz2","Cd3e","Cd3d","Cd8a",
           "Jchain","Mzb1","Cd79a","Ms4a1","Cpa3","Kit",
           "Itgax","Flt3","Col1a1","Col3a1","Dcn","Acta2","Myh11",
           "Pecam1","Cldn5","Lyve1")
avg <- as.matrix(AverageExpression(ss, features=genes, assays="Spatial")$Spatial)
cat("\n=== 8 stromal clusters: top markers ===\n")
for(cl in colnames(avg)){
  v <- sort(avg[,cl], decreasing=TRUE)[1:6]
  cat(sprintf("%-6s: %s\n", cl, paste(sprintf("%s=%.1f",names(v),v),collapse="  ")))
}





ss$cl <- strom$Stromal.Clusters[match(colnames(ss), strom$Barcode)]
Idents(ss) <- "cl"
genes <- c("Jchain","Mzb1","Col1a1","Col3a1","Acta2","Dcn","C1qa","Csf1r","Lyz2",
           "Cd3e","Cd8a","Ccl5","Pecam1","Cldn5","Cpa3","Kit")
m <- round(as.matrix(AverageExpression(ss, features=genes, assays="Spatial")$Spatial), 2)
print(m)

type_map <- c(
  "Cluster 3"="Fibroblast","Cluster 19"="Fibroblast","Cluster 14"="Fibroblast",
  "Cluster 4"="Fibroblast","Cluster 6"="Fibroblast",
  "Cluster 17"="Macrophage","Cluster 18"="Macrophage","Cluster 7"="Macrophage",
  "Cluster 13"="Tcell","Cluster 5"="Tcell","Cluster 11"="Tcell",
  "Cluster 1"="Plasma","Cluster 2"="Plasma","Cluster 10"="Plasma",
  "Cluster 15"="Plasma","Cluster 12"="Plasma","Cluster 16"="Plasma",
  "Cluster 8"="Plasma","Cluster 9"="Plasma"
)
strom$celltype <- type_map[strom$Stromal.Clusters]
tab <- table(strom$celltype)
tab
round(100*tab/sum(tab), 1)

strom$celltype <- type_map[strom$Stromal.Clusters]
tab <- table(strom$celltype)
round(100*tab/sum(tab), 1)   # percentages
tab                          # raw counts



type_map <- c(
  # CONFIDENT
  "Cluster 3"="Fibroblast","Cluster 19"="Fibroblast",
  "Cluster 17"="Macrophage","Cluster 7"="Macrophage","Cluster 18"="Macrophage",
  "Cluster 13"="Tcell","Cluster 5"="Tcell","Cluster 11"="Tcell",
  "Cluster 1"="Plasma","Cluster 2"="Plasma","Cluster 10"="Plasma",
  "Cluster 12"="Plasma","Cluster 15"="Plasma",
  # UNCERTAIN -> excluded
  "Cluster 14"="Mixed_unsure","Cluster 4"="Mixed_unsure","Cluster 6"="Mixed_unsure",
  "Cluster 8"="Mixed_unsure","Cluster 9"="Mixed_unsure","Cluster 16"="Mixed_unsure"
)
strom$celltype <- type_map[strom$Stromal.Clusters]








library(Seurat); library(dplyr)

# ---- load everything ----
strom <- read.csv("Stromal Clusters.csv")            # Barcode, Stromal.Clusters
plas  <- read.csv("Plasma.csv")                       # Barcode, Plasma
macs  <- read.csv("Macrophages.csv")                  # Barcode, Macrophages
mixed <- read.csv("MixedReanalysis.csv")              # Barcode, MixedReanalysis
epiz  <- read.csv("epi_zones_all_genotypes.csv")      # Barcode, zone, genotype
gen   <- read.csv("Genotypes.csv")

# ---- 1. base broad-type map from original 19 stromal clusters ----
base_map <- c(
  "Cluster 1"="Plasma","Cluster 2"="Plasma","Cluster 10"="Plasma","Cluster 12"="Plasma",
  "Cluster 3"="Fibroblast","Cluster 4"="Fibroblast","Cluster 6"="Fibroblast",
  "Cluster 9"="Fibroblast","Cluster 19"="Fibroblast",
  "Cluster 7"="Macrophage","Cluster 15"="Macrophage",
  "Cluster 5"="Tcell","Cluster 13"="Tcell",
  "Cluster 16"="Endothelial","Cluster 11"="Endothelial",
  "Cluster 17"="DC",
  "Cluster 14"="Myeloid",
  "Cluster 18"="Mixed_unsure","Cluster 8"="Mixed_unsure"
)
strom$celltype <- base_map[strom$Stromal.Clusters]

# ---- 2. apply reanalysis overrides (these take precedence) ----
# helper: normalize + map, returns named vector barcode->final label ("DROP" or "EPI" special)
lab <- setNames(strom$celltype, strom$Barcode)

# Macrophages reanalysis
mac_map <- c("Macrophages"="Macrophage","Plasma Cells"="Plasma","Fibroblast"="Fibroblast",
             "Mixed"="DROP","Discard"="DROP","Discard_"="DROP")
lab[macs$Barcode] <- mac_map[macs$Macrophages]

# Mixed reanalysis
mix_map <- c("Macrophage"="Macrophage","Plasma"="Plasma","Cluster 3"="Mixed_c3",
             "Epithelium"="EPI","Epithelium 1"="EPI","Lymphoid Follicle(discard)"="DROP")
lab[mixed$Barcode] <- mix_map[mixed$MixedReanalysis]

# Plasma reanalysis (normalize epithelium capitalization)
plas$lab2 <- ifelse(grepl("^epithelium|^Epithelium", plas$Plasma, ignore.case=TRUE),
                    "EPI", "Plasma")
lab[plas$Barcode] <- plas$lab2

# ---- 3. project EPI-flagged cells onto the 4 epithelial zones ----
epi_bc <- names(lab)[lab == "EPI" & !is.na(lab)]
epi_bc <- epi_bc[epi_bc %in% colnames(mat2)]
cat("projecting", length(epi_bc), "epithelial-contaminant cells onto zones\n")

# reference = your named control epithelial object 'e' (must be in session)
q2 <- CreateSeuratObject(mat2[, epi_bc], assay="Spatial")
q2 <- NormalizeData(q2)
anch <- FindTransferAnchors(reference=e, query=q2, dims=1:15, reference.reduction="pca")
pr   <- TransferData(anchorset=anch, refdata=e$zone_named, dims=1:15)
proj_zone <- setNames(as.character(pr$predicted.id), colnames(q2))
lab[names(proj_zone)] <- proj_zone   # now Crypt/Junction/Villus/Villus_tip

# ---- 4. label Mixed_c3 (compute markers so we can name it) ----
c3_bc <- names(lab)[lab == "Mixed_c3" & !is.na(lab)]
c3_bc <- c3_bc[c3_bc %in% colnames(mat2)]
cat("\nMixed cluster 3 =", length(c3_bc), "cells. Top markers:\n")
if(length(c3_bc) > 20){
  tmp <- CreateSeuratObject(mat2[, c3_bc], assay="Spatial"); tmp <- NormalizeData(tmp)
  rest <- sample(setdiff(colnames(mat2), c3_bc), min(5000, ncol(mat2)-length(c3_bc)))
  comb <- CreateSeuratObject(mat2[, c(c3_bc, rest)], assay="Spatial"); comb <- NormalizeData(comb)
  comb$grp <- ifelse(colnames(comb) %in% c3_bc, "c3", "rest"); Idents(comb) <- "grp"
  m3 <- FindMarkers(comb, ident.1="c3", only.pos=TRUE, min.pct=0.1, logfc.threshold=0.25)
  print(head(m3[order(-m3$avg_log2FC),], 20))
}

# ---- 5. assemble the FINAL table: epithelial zones + stromal types ----
# start with all epithelial cells from the zonation
final <- data.frame(Barcode = epiz$Barcode, celltype = epiz$zone, stringsAsFactors=FALSE)

# add stromal cells (drop DROP and NA)
stro_final <- data.frame(Barcode = names(lab), celltype = unname(lab), stringsAsFactors=FALSE)
stro_final <- stro_final[!is.na(stro_final$celltype) & !(stro_final$celltype %in% c("DROP","EPI")), ]

# the projected epithelial-contaminant cells are already Crypt/Junction/etc in lab,
# so they're included in stro_final with zone labels — good.
# combine, dropping any barcode already in the epithelial zonation to avoid dupes
stro_final <- stro_final[!(stro_final$Barcode %in% final$Barcode), ]
final <- rbind(final, stro_final)

# attach genotype
final$genotype <- gen$Genotypes[match(final$Barcode, gen$Barcode)]
final <- final[!is.na(final$genotype), ]

cat("\n=== FINAL cell type counts ===\n")
print(table(final$celltype))
cat("\ntotal cells:", nrow(final), "\n")

write.csv(final, "FINAL_all_celltypes.csv", row.names=FALSE)
cat("wrote FINAL_all_celltypes.csv\n")




# ---- install CellChat if needed (one time) ----
if (!requireNamespace("CellChat", quietly=TRUE)) {
  if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
  BiocManager::install(c("ComplexHeatmap","BiocNeighbors"))
  install.packages("devtools")
  devtools::install_github("jinworks/CellChat")
}

library(CellChat); library(Seurat); library(dplyr)

# ================= SETUP =================
final <- read.csv("FINAL_all_celltypes.csv")   # Barcode, celltype, genotype
proj  <- read.csv("Default Pipeline Analysis-Spatial-Projection.csv")
MPP   <- 0.3443589

# cell types we want: 4 epithelial zones + 4 stromal
keep_types <- c("Crypt","Junction","Villus","Villus_tip",
                "Macrophage","Tcell","Fibroblast","DC")
final <- final[final$celltype %in% keep_types, ]

# coordinates in MICRONS, matched to cells
final$x <- proj$X.Coordinate[match(final$Barcode, proj$Barcode)] * MPP
final$y <- proj$Y.Coordinate[match(final$Barcode, proj$Barcode)] * MPP
final <- final[!is.na(final$x) & final$Barcode %in% colnames(mat2), ]

# ================= FUNCTION: run CellChat for ONE genotype =================
run_cellchat_geno <- function(geno) {
  cat("\n==== CellChat:", geno, "====\n")
  sub <- final[final$genotype == geno, ]
  cat("cells:", nrow(sub), " types:", paste(table(sub$celltype), collapse="/"), "\n")
  
  # expression (normalized) for these cells
  bc  <- sub$Barcode
  obj <- CreateSeuratObject(mat2[, bc], assay="Spatial")
  obj <- NormalizeData(obj)
  data.input <- GetAssayData(obj, layer="data")
  
  meta <- data.frame(labels = sub$celltype, row.names = bc)
  coords <- as.matrix(sub[, c("x","y")]); rownames(coords) <- bc
  
  # spatial factors (Visium HD: ratio and tol; conversion already in microns)
  spatial.factors <- data.frame(ratio = 1, tol = 5)
  
  cc <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                       datatype = "spatial", coordinates = coords,
                       spatial.factors = spatial.factors)
  cc@DB <- CellChatDB.mouse
  
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  
  # Rachel's parameters: truncatedMean trim 0.1%, spatial ranges
  cc <- computeCommunProb(cc, type = "truncatedMean", trim = 0.001,
                          distance.use = TRUE,
                          interaction.range = 200, scale.distance = 1,
                          contact.range = 20)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc
}

# ---- run CONTROL first (fastest to sanity-check) ----
cc_ctrl <- run_cellchat_geno("Villin CreERT2(control)")

# save it
saveRDS(cc_ctrl, "cellchat_control.rds")

# ---- quick look at results ----
cat("\n=== Pathways detected (control) ===\n")
print(cc_ctrl@netP$pathways)

# overall interaction strength between cell types
cat("\n=== Number of interactions between types ===\n")
print(cc_ctrl@net$count)

BiocManager::install(c("BiocNeighbors","ComplexHeatmap"))
remotes::install_github("jinworks/CellChat")


pkgbuild::has_build_tools(debug=TRUE)

remotes::install_github("jinworks/CellChat")

pkgbuild::has_build_tools(debug=TRUE)
library(CellChat)