library(Seurat)

epi_bc_all <- epi_all$Barcode
epi_bc_all <- epi_bc_all[epi_bc_all %in% colnames(mat)]

control_bc <- colnames(e)                      # 17777 control reference cells
query_bc   <- setdiff(epi_bc_all, control_bc)  # mutant epithelium to classify
cat("reference:", length(control_bc), " query:", length(query_bc), "\n")

q <- CreateSeuratObject(mat[, query_bc], assay = "Spatial")
q <- NormalizeData(q)

anchors <- FindTransferAnchors(reference = e, query = q,
                               dims = 1:15, reference.reduction = "pca")
pred <- TransferData(anchorset = anchors, refdata = e$zone_named, dims = 1:15)
q$zone_named <- pred$predicted.id
table(q$zone_named)

all_zones <- rbind(
  data.frame(Barcode = colnames(e), zone = as.character(e$zone_named)),
  data.frame(Barcode = colnames(q), zone = as.character(q$zone_named))
)
gen <- read.csv("Genotypes.csv")
all_zones$genotype <- gen$Genotypes[match(all_zones$Barcode, gen$Barcode)]

write.csv(all_zones, "epi_zones_all_genotypes.csv", row.names = FALSE)
table(all_zones$zone, all_zones$genotype)

# overall distribution of confidence
summary(pred$prediction.score.max)

# attach scores to the query cells with their genotype, then compare by genotype
q$score <- pred$prediction.score.max
q$genotype <- gen$Genotypes[match(colnames(q), gen$Barcode)]

# median confidence per genotype
tapply(q$score, q$genotype, median)

# and per genotype x predicted zone
aggregate(score ~ genotype + zone_named, data = q@meta.data, FUN = median)


remotes::install_github("jinworks/CellChat")

library(CellChat)
packageVersion("CellChat")

library(dplyr); library(ggplot2); library(tidyr)

# ---- pull zone proportions from the live object ----
az <- read.csv("epi_zones_all_genotypes.csv")   # Barcode, zone, genotype

# clean genotype names + set Braf-dose order (Erk arms after)
geno_map <- c(
  "Villin CreERT2(control)"                        = "Control",
  "VillinCreERT2; BrafV600E/+"                      = "Braf V600E/+",
  "Villin CreERT2; BrafV600E/V600E"                 = "Braf V600E/V600E",
  "Villin CreERT2; Erk1/Erk2 DKO"                   = "Erk DKO",
  "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"      = "Braf V600E/+;\nErk DKO",
  "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"  = "Braf V600E/V600E;\nErk DKO"
)
az$g <- geno_map[az$genotype]
geno_order <- c("Control","Braf V600E/+","Braf V600E/V600E",
                "Erk DKO","Braf V600E/+;\nErk DKO","Braf V600E/V600E;\nErk DKO")
zone_order <- c("Crypt","Junction","Villus","Villus_tip")

prop <- az %>%
  filter(!is.na(g)) %>%
  count(g, zone) %>%
  group_by(g) %>% mutate(pct = 100*n/sum(n)) %>% ungroup() %>%
  mutate(g = factor(g, levels=geno_order),
         zone = factor(zone, levels=zone_order))

# ---- Nature-style palette (crypt=deep blue base -> tip=warm red apex) ----
pal <- c("Crypt"="#2C5F8A", "Junction"="#7BA7C4",
         "Villus"="#E8C46B", "Villus_tip"="#C0392B")
zone_labs <- c("Crypt","Junction","Villus","Villus tip")

# ---- stacked bar, minimal aesthetic ----
p <- ggplot(prop, aes(g, pct, fill=zone)) +
  geom_col(width=0.72, color="white", linewidth=0.35) +
  scale_fill_manual(values=pal, labels=zone_labs, name="Epithelial zone") +
  scale_y_continuous(expand=expansion(mult=c(0,0.02)),
                     breaks=seq(0,100,25), labels=function(x) paste0(x,"%")) +
  labs(x=NULL, y="Proportion of epithelium") +
  # dose bracket annotation
  annotate("segment", x=0.6, xend=3.4, y=-8, yend=-8, linewidth=0.5) +
  annotate("text", x=2, y=-12, label="Increasing Braf dose", size=3.4, fontface="italic") +
  annotate("segment", x=3.6, xend=6.4, y=-8, yend=-8, linewidth=0.5, color="grey45") +
  annotate("text", x=5, y=-12, label="+ Erk1/2 deletion", size=3.4, fontface="italic", color="grey35") +
  coord_cartesian(ylim=c(0,100), clip="off") +
  theme_classic(base_size=12) +
  theme(
    axis.text.x = element_text(size=9.5, color="grey15", lineheight=0.9),
    axis.text.y = element_text(size=10, color="grey15"),
    axis.title.y= element_text(size=11, margin=margin(r=8)),
    axis.line   = element_line(linewidth=0.4, color="grey30"),
    axis.ticks  = element_line(linewidth=0.4, color="grey30"),
    legend.position="right",
    legend.key.size=unit(0.5,"cm"),
    legend.title=element_text(size=10, face="bold"),
    legend.text=element_text(size=9.5),
    plot.margin=margin(t=10,r=10,b=34,l=12)
  )

ggsave("zone_proportions_figure.pdf", p, width=8.2, height=4.6, device=cairo_pdf)
ggsave("zone_proportions_figure.png", p, width=8.2, height=4.6, dpi=400)
cat("saved zone_proportions_figure.pdf/.png\n")
print(prop %>% select(g,zone,pct) %>% arrange(g,zone), n=24)




library(dplyr); library(ggplot2); library(tidyr)

az <- read.csv("epi_zones_all_genotypes.csv")

geno_map <- c(
  "Villin CreERT2(control)"                        = "Control",
  "VillinCreERT2; BrafV600E/+"                      = "Braf V600E/+",
  "Villin CreERT2; BrafV600E/V600E"                 = "Braf V600E/V600E",
  "Villin CreERT2; Erk1/Erk2 DKO"                   = "Erk DKO",
  "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"      = "Braf V600E/+;\nErk DKO",
  "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"  = "Braf V600E/V600E;\nErk DKO"
)
az$g <- geno_map[az$genotype]
geno_order <- c("Control","Braf V600E/+","Braf V600E/V600E",
                "Erk DKO","Braf V600E/+;\nErk DKO","Braf V600E/V600E;\nErk DKO")
zone_order <- c("Crypt","Junction","Villus","Villus_tip")

prop <- az %>%
  filter(!is.na(g)) %>%
  count(g, zone) %>%
  group_by(g) %>% mutate(pct = 100*n/sum(n)) %>% ungroup() %>%
  mutate(g = factor(g, levels=geno_order),
         zone = factor(zone, levels=zone_order))

# label positions: midpoint of each stacked segment; hide labels < 4%
prop <- prop %>%
  group_by(g) %>%
  arrange(zone, .by_group=TRUE) %>%
  mutate(ypos = cumsum(pct) - pct/2,
         lab = ifelse(pct >= 4, paste0(round(pct), "%"), "")) %>%
  ungroup()

pal <- c("Crypt"="#2C5F8A", "Junction"="#7BA7C4",
         "Villus"="#E8C46B", "Villus_tip"="#C0392B")
# white text on dark segments, dark text on light — pick per zone for contrast
txt_col <- c("Crypt"="white","Junction"="grey15","Villus"="grey15","Villus_tip"="white")

p <- ggplot(prop, aes(g, pct, fill=zone)) +
  geom_col(width=0.72, color="white", linewidth=0.35) +
  geom_text(aes(y=ypos, label=lab, color=zone), size=3.1, fontface="bold", show.legend=FALSE) +
  scale_fill_manual(values=pal, labels=c("Crypt","Junction","Villus","Villus tip"), name="Epithelial zone") +
  scale_color_manual(values=txt_col, guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0,0.02)),
                     breaks=seq(0,100,25), labels=function(x) paste0(x,"%")) +
  labs(x=NULL, y="Proportion of epithelium") +
  annotate("segment", x=0.6, xend=3.4, y=-8, yend=-8, linewidth=0.5) +
  annotate("text", x=2, y=-12, label="Increasing Braf dose", size=3.4, fontface="italic") +
  annotate("segment", x=3.6, xend=6.4, y=-8, yend=-8, linewidth=0.5, color="grey45") +
  annotate("text", x=5, y=-12, label="+ Erk1/2 deletion", size=3.4, fontface="italic", color="grey35") +
  coord_cartesian(ylim=c(0,100), clip="off") +
  theme_classic(base_size=12) +
  theme(
    axis.text.x = element_text(size=9.5, color="grey15", lineheight=0.9),
    axis.text.y = element_text(size=10, color="grey15"),
    axis.title.y= element_text(size=11, margin=margin(r=8)),
    axis.line   = element_line(linewidth=0.4, color="grey30"),
    axis.ticks  = element_line(linewidth=0.4, color="grey30"),
    legend.position="right",
    legend.key.size=unit(0.5,"cm"),
    legend.title=element_text(size=10, face="bold"),
    legend.text=element_text(size=9.5),
    plot.margin=margin(t=10,r=10,b=34,l=12)
  )

ggsave("zone_proportions_figure.pdf", p, width=8.2, height=4.6, device=cairo_pdf)
ggsave("zone_proportions_figure.png", p, width=8.2, height=4.6, dpi=400)
cat("saved with percentage labels\n")






library(dplyr); library(ggplot2)

az <- read.csv("epi_zones_all_genotypes.csv")

geno_map <- c(
  "Villin CreERT2(control)"                        = "Control",
  "VillinCreERT2; BrafV600E/+"                      = "Braf V600E/+",
  "Villin CreERT2; BrafV600E/V600E"                 = "Braf V600E/V600E",
  "Villin CreERT2; Erk1/Erk2 DKO"                   = "Erk DKO",
  "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"      = "Braf V600E/+;\nErk DKO",
  "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"  = "Braf V600E/V600E;\nErk DKO"
)
az$g <- geno_map[az$genotype]
geno_order <- c("Control","Braf V600E/+","Braf V600E/V600E",
                "Erk DKO","Braf V600E/+;\nErk DKO","Braf V600E/V600E;\nErk DKO")
zone_order <- c("Crypt","Junction","Villus","Villus_tip")

prop <- az %>%
  filter(!is.na(g)) %>%
  count(g, zone) %>%
  group_by(g) %>% mutate(pct = 100*n/sum(n)) %>% ungroup() %>%
  mutate(g = factor(g, levels=geno_order),
         zone = factor(zone, levels=zone_order),
         lab = ifelse(pct >= 4, paste0(round(pct), "%"), ""))

pal <- c("Crypt"="#2C5F8A","Junction"="#7BA7C4","Villus"="#E8C46B","Villus_tip"="#C0392B")

p <- ggplot(prop, aes(g, pct, fill=zone)) +
  geom_col(width=0.72, color="white", linewidth=0.35) +
  geom_text(aes(label=lab), position=position_stack(vjust=0.5),
            size=3.1, fontface="bold", color="grey10") +
  scale_fill_manual(values=pal, labels=c("Crypt","Junction","Villus","Villus tip"),
                    name="Epithelial zone") +
  scale_y_continuous(expand=expansion(mult=c(0,0.02)),
                     breaks=seq(0,100,25), labels=function(x) paste0(x,"%")) +
  labs(x=NULL, y="Proportion of epithelium") +
  annotate("segment", x=0.6, xend=3.4, y=-8, yend=-8, linewidth=0.5) +
  annotate("text", x=2, y=-12, label="Increasing Braf dose", size=3.4, fontface="italic") +
  annotate("segment", x=3.6, xend=6.4, y=-8, yend=-8, linewidth=0.5, color="grey45") +
  annotate("text", x=5, y=-12, label="+ Erk1/2 deletion", size=3.4, fontface="italic", color="grey35") +
  coord_cartesian(ylim=c(0,100), clip="off") +
  theme_classic(base_size=12) +
  theme(
    axis.text.x = element_text(size=9.5, color="grey15", lineheight=0.9),
    axis.text.y = element_text(size=10, color="grey15"),
    axis.title.y= element_text(size=11, margin=margin(r=8)),
    axis.line   = element_line(linewidth=0.4, color="grey30"),
    axis.ticks  = element_line(linewidth=0.4, color="grey30"),
    legend.position="right", legend.key.size=unit(0.5,"cm"),
    legend.title=element_text(size=10, face="bold"), legend.text=element_text(size=9.5),
    plot.margin=margin(t=10,r=10,b=34,l=12)
  )

ggsave("zone_proportions_figure.pdf", p, width=8.2, height=4.6, device=cairo_pdf)
ggsave("zone_proportions_figure.png", p, width=8.2, height=4.6, dpi=400)
cat("done\n")

BiocManager::install("Biobase")

remotes::install_github("jinworks/CellChat")

library(Biobase)
library(CellChat
        )

packageVersion("CellChat")


library(Seurat); library(dplyr); library(jsonlite)

final <- read.csv("FINAL_all_celltypes.csv")
proj  <- read.csv("Default Pipeline Analysis-Spatial-Projection.csv")
scalef <- fromJSON("scalefactors_json.json")

keep_types <- c("Crypt","Junction","Villus","Villus_tip","Macrophage","Tcell","Fibroblast","DC")
final <- final[final$celltype %in% keep_types, ]
final$x <- proj$X.Coordinate[match(final$Barcode, proj$Barcode)]
final$y <- proj$Y.Coordinate[match(final$Barcode, proj$Barcode)]
final <- final[!is.na(final$x) & final$Barcode %in% colnames(mat2), ]

run_geno <- function(geno) {
  sub <- final[final$genotype == geno, ]
  cat("\n====", geno, "| cells:", nrow(sub), "====\n"); print(table(sub$celltype))
  bc <- sub$Barcode
  obj <- NormalizeData(CreateSeuratObject(mat2[, bc], assay="Spatial"))
  data.input <- GetAssayData(obj, layer="data")
  meta <- data.frame(labels=sub$celltype, samples=factor("s1"), row.names=bc)
  locs <- as.matrix(sub[, c("x","y")]); rownames(locs) <- bc
  
  spot.size <- 8
  conv <- spot.size / scalef$spot_diameter_fullres
  spatial.factors <- data.frame(ratio=conv, tol=spot.size/2)
  
  cc <- createCellChat(object=data.input, meta=meta, group.by="labels",
                       datatype="spatial", coordinates=locs, spatial.factors=spatial.factors)
  cc@DB <- subsetDB(CellChatDB.mouse, search=c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"), key="annotation")
  cc <- subsetData(cc)
  options(future.globals.maxSize=13000*1024^2)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc, variable.both=FALSE)
  cc <- computeCommunProb(cc, type="truncatedMean", trim=0.001,
                          distance.use=TRUE, interaction.range=200,
                          scale.distance=conv, contact.dependent=TRUE, contact.range=20)
  cc <- filterCommunication(cc, min.cells=10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc
}

cc_ctrl <- run_geno("Villin CreERT2(control)")
saveRDS(cc_ctrl, "cellchat_control.rds")
cat("\n=== pathways ===\n"); print(cc_ctrl@netP$pathways)
cat("\n=== interaction counts (rows=sender, cols=receiver) ===\n"); print(cc_ctrl@net$count)



run_geno <- function(geno) {
  sub <- final[final$genotype == geno, ]
  cat("\n====", geno, "| cells:", nrow(sub), "====\n"); print(table(sub$celltype))
  bc <- sub$Barcode
  obj <- NormalizeData(CreateSeuratObject(mat2[, bc], assay="Spatial"))
  data.input <- GetAssayData(obj, layer="data")
  meta <- data.frame(labels=sub$celltype, samples=factor("s1"), row.names=bc)
  locs <- as.matrix(sub[, c("x","y")]); rownames(locs) <- bc
  
  MPP <- 0.3443589                       # microns per pixel (from your scalefactors)
  spatial.factors <- data.frame(ratio = MPP, tol = 8/2)   # ratio=MPP, tol=half spot
  
  cc <- createCellChat(object=data.input, meta=meta, group.by="labels",
                       datatype="spatial", coordinates=locs, spatial.factors=spatial.factors)
  cc@DB <- subsetDB(CellChatDB.mouse, search=c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"), key="annotation")
  cc <- subsetData(cc)
  options(future.globals.maxSize=13000*1024^2)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc, variable.both=FALSE)
  cc <- computeCommunProb(cc, type="truncatedMean", trim=0.001,
                          distance.use=TRUE, interaction.range=200,
                          scale.distance=MPP, contact.dependent=TRUE, contact.range=20)
  cc <- filterCommunication(cc, min.cells=10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc
}

cc_ctrl <- run_geno("Villin CreERT2(control)")
saveRDS(cc_ctrl, "cellchat_control.rds")
cat("\n=== pathways ===\n"); print(cc_ctrl@netP$pathways)
cat("\n=== interaction counts ===\n"); print(cc_ctrl@net$count)








run_geno <- function(geno) {
  sub <- final[final$genotype == geno, ]
  cat("\n====", geno, "| cells:", nrow(sub), "====\n"); print(table(sub$celltype))
  bc <- sub$Barcode
  obj <- NormalizeData(CreateSeuratObject(mat2[, bc], assay="Spatial"))
  data.input <- GetAssayData(obj, layer="data")
  meta <- data.frame(labels=sub$celltype, samples=factor("s1"), row.names=bc)
  locs <- as.matrix(sub[, c("x","y")]); rownames(locs) <- bc
  
  MPP <- 0.3443589
  spatial.factors <- data.frame(ratio = MPP, tol = 4)
  cc <- computeCommunProb(cc, type="truncatedMean", trim=0.001,
                        distance.use=TRUE, interaction.range=100,   # was 200
                        scale.distance=0.3443589, contact.dependent=TRUE,
                        contact.range=20, nboot=20)
  
  
  cc <- createCellChat(object=data.input, meta=meta, group.by="labels",
                       datatype="spatial", coordinates=locs, spatial.factors=spatial.factors)
  cc@DB <- subsetDB(CellChatDB.mouse, search=c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"), key="annotation")
  cc <- subsetData(cc)
  options(future.globals.maxSize=13000*1024^2)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc, variable.both=FALSE)
  cc <- computeCommunProb(cc, type="truncatedMean", trim=0.001,
                          distance.use=TRUE, interaction.range=200,
                          scale.distance=MPP, contact.dependent=TRUE, contact.range=20)
  cc <- filterCommunication(cc, min.cells=10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc
}




run_geno <- function(geno) {
  sub <- final[final$genotype == geno, ]
  cat("\n====", geno, "| cells:", nrow(sub), "====\n"); print(table(sub$celltype))
  bc <- sub$Barcode
  obj <- NormalizeData(CreateSeuratObject(mat2[, bc], assay="Spatial"))
  data.input <- GetAssayData(obj, layer="data")
  meta <- data.frame(labels=sub$celltype, samples=factor("s1"), row.names=bc)
  locs <- as.matrix(sub[, c("x","y")]); rownames(locs) <- bc

  MPP <- 0.3443589
  spatial.factors <- data.frame(ratio = MPP, tol = 4)

  cc <- createCellChat(object=data.input, meta=meta, group.by="labels",
                       datatype="spatial", coordinates=locs, spatial.factors=spatial.factors)
  cc@DB <- subsetDB(CellChatDB.mouse, search=c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"), key="annotation")
  cc <- subsetData(cc)
  options(future.globals.maxSize=13000*1024^2)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc, variable.both=FALSE)
  cc <- computeCommunProb(cc, type="truncatedMean", trim=0.001,
                          distance.use=TRUE, interaction.range=100,
                          scale.distance=MPP, contact.dependent=TRUE,
                          contact.range=20, nboot=20)
  cc <- filterCommunication(cc, min.cells=10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc
}

cc_ctrl <- run_geno("Villin CreERT2(control)")
saveRDS(cc_ctrl, "cellchat_control.rds")
cat("\n=== pathways ===\n"); print(cc_ctrl@netP$pathways)
cat("\n=== interaction counts ===\n"); print(cc_ctrl@net$count)






run_geno <- function(geno) {
  sub <- final[final$genotype == geno, ]
  cat("\n====", geno, "| cells:", nrow(sub), "====\n"); print(table(sub$celltype))
  bc <- sub$Barcode
  obj <- NormalizeData(CreateSeuratObject(mat2[, bc], assay="Spatial"))
  data.input <- GetAssayData(obj, layer="data")
  meta <- data.frame(labels=sub$celltype, samples=factor("s1"), row.names=bc)
  locs <- as.matrix(sub[, c("x","y")]); rownames(locs) <- bc
  
  MPP <- 0.3443589
  spatial.factors <- data.frame(ratio = MPP, tol = 4)
  
  cc <- createCellChat(object=data.input, meta=meta, group.by="labels",
                       datatype="spatial", coordinates=locs, spatial.factors=spatial.factors)
  cc@DB <- subsetDB(CellChatDB.mouse, search=c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"), key="annotation")
  cc <- subsetData(cc)
  options(future.globals.maxSize=13000*1024^2)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc, variable.both=FALSE)
  cc <- computeCommunProb(cc, type="truncatedMean", trim=0.001,
                          distance.use=TRUE, interaction.range=100,
                          scale.distance=MPP, contact.dependent=TRUE,
                          contact.range=20, nboot=20)
  cc <- filterCommunication(cc, min.cells=10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc
}

cc_ctrl <- run_geno("Villin CreERT2(control)")
saveRDS(cc_ctrl, "cellchat_control.rds")
cat("\n=== pathways ===\n"); print(cc_ctrl@netP$pathways)
cat("\n=== interaction counts ===\n"); print(cc_ctrl@net$count)




library(CellChat)
cc_ctrl <- readRDS("cellchat_control.rds")
cc_ctrl <- netAnalysis_computeCentrality(cc_ctrl, slot.name="netP")

grp_cols <- c("Crypt"="#6A3D9A","DC"="#E7298A","Fibroblast"="#1B9E77",
              "Junction"="#F4C430","Macrophage"="#D95F02","Tcell"="#7570B3",
              "Villus"="#2C6FBB","Villus_tip"="#C0392B")
grp_cols <- grp_cols[levels(cc_ctrl@idents)]
grp_sizes <- as.numeric(table(cc_ctrl@idents))

# rank pathways by strength, keep top 25 for readable heatmaps
pw_strength <- sort(sapply(seq_along(cc_ctrl@netP$pathways), function(i)
  sum(cc_ctrl@netP$prob[,,i])), decreasing=TRUE)
names(pw_strength) <- cc_ctrl@netP$pathways[order(-sapply(seq_along(cc_ctrl@netP$pathways),
                                                          function(i) sum(cc_ctrl@netP$prob[,,i])))]
top_pw <- names(pw_strength)[1:25]

# ============ FIG 1: circle networks (count + strength) ============
pdf("fig1_networks.pdf", width=13, height=6.5)
par(mfrow=c(1,2), xpd=TRUE, mar=c(1,1,3,1))
netVisual_circle(cc_ctrl@net$count, vertex.weight=grp_sizes, weight.scale=TRUE,
                 label.edge=FALSE, color.use=grp_cols, title.name="Number of interactions")
netVisual_circle(cc_ctrl@net$weight, vertex.weight=grp_sizes, weight.scale=TRUE,
                 label.edge=FALSE, color.use=grp_cols, title.name="Interaction strength")
dev.off()

# ============ FIG 2: incoming/outgoing role heatmaps (top 25, readable) ============
pdf("fig2_signaling_roles.pdf", width=9, height=10)
ph_out <- netAnalysis_signalingRole_heatmap(cc_ctrl, pattern="outgoing",
                                            signaling=top_pw, color.use=grp_cols, width=9, height=16, font.size=7)
ph_in  <- netAnalysis_signalingRole_heatmap(cc_ctrl, pattern="incoming",
                                            signaling=top_pw, color.use=grp_cols, width=9, height=16, font.size=7)
print(ph_out); print(ph_in)
dev.off()

# ============ FIG 3: epithelium -> fibroblast/macrophage (bubble, L-R pairs) ============
pdf("fig3_epi_to_stroma.pdf", width=7, height=9)
netVisual_bubble(cc_ctrl,
                 sources.use=c("Crypt","Junction","Villus","Villus_tip"),
                 targets.use=c("Fibroblast","Macrophage"),
                 remove.isolate=TRUE, font.size=8)
dev.off()

# ============ FIG 4: fibroblast/macrophage -> epithelium (reverse) ============
pdf("fig4_stroma_to_epi.pdf", width=7, height=9)
netVisual_bubble(cc_ctrl,
                 sources.use=c("Fibroblast","Macrophage"),
                 targets.use=c("Crypt","Junction","Villus","Villus_tip"),
                 remove.isolate=TRUE, font.size=8)
dev.off()

# ============ FIG 5: key pathway heatmaps (one per page, tied to your story) ============
key_paths <- intersect(c("BMP","TGFb","WNT","SPP1","PDGF","FGF","GRN","MK"),
                       cc_ctrl@netP$pathways)
pdf("fig5_key_pathways.pdf", width=6.5, height=5.5)
for(pw in key_paths){
  print(netVisual_heatmap(cc_ctrl, signaling=pw, color.heatmap="Reds",
                          color.use=grp_cols, title.name=paste(pw, "signaling")))
}
dev.off()

cat("saved fig1-fig5\n")
cat("key pathways found:", paste(key_paths, collapse=", "), "\n")





# reconstruct: original base map (19 clusters -> broad type) vs final assignment
strom <- read.csv("Stromal Clusters.csv")

base_map <- c(
  "Cluster 1"="Plasma","Cluster 2"="Plasma","Cluster 10"="Plasma","Cluster 12"="Plasma",
  "Cluster 3"="Fibroblast","Cluster 4"="Fibroblast","Cluster 6"="Fibroblast",
  "Cluster 9"="Fibroblast","Cluster 19"="Fibroblast",
  "Cluster 7"="Macrophage","Cluster 15"="Macrophage",
  "Cluster 5"="Tcell","Cluster 13"="Tcell",
  "Cluster 16"="Endothelial","Cluster 11"="Endothelial",
  "Cluster 17"="DC","Cluster 14"="Myeloid",
  "Cluster 18"="Mixed_unsure","Cluster 8"="Mixed_unsure"
)
strom$original <- base_map[strom$Stromal.Clusters]

# final labels
final <- read.csv("FINAL_all_celltypes.csv")
strom$final <- final$celltype[match(strom$Barcode, final$Barcode)]
strom$final[is.na(strom$final)] <- "removed/epithelial"

# how many cells changed label?
cat("total stromal cells:", nrow(strom), "\n")
cat("unchanged:", sum(strom$original == strom$final, na.rm=TRUE), "\n")
cat("reassigned:", sum(strom$original != strom$final, na.rm=TRUE), "\n\n")

# the full transition table (original -> final)
print(table(strom$original, strom$final))










final <- read.csv("FINAL_all_celltypes.csv")

# what's in the final annotation
cat("=== FINAL annotation (all labeled cells) ===\n")
print(table(final$celltype))

# what CellChat actually used
kept <- c("Crypt","Junction","Villus","Villus_tip","Macrophage","Tcell","Fibroblast","DC")
final$in_cellchat <- ifelse(final$celltype %in% kept, "USED in CellChat", "excluded from CellChat")

cat("\n=== used vs excluded in CellChat ===\n")
print(table(final$in_cellchat))

cat("\n=== which cell types were excluded ===\n")
print(table(final$celltype[final$in_cellchat=="excluded from CellChat"]))

# and the original stromal cells that never made it to final at all (dropped entirely)
strom <- read.csv("Stromal Clusters.csv")
dropped <- sum(!(strom$Barcode %in% final$Barcode))
cat("\noriginal stromal cells dropped entirely (not in final):", dropped,
    "of", nrow(strom), "\n")





# ---- CellChat parameters (applied identically to every genotype) ----
MPP <- 0.3443589   # microns per pixel

cc@DB <- subsetDB(CellChatDB.mouse,
                  search = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"),
                  key = "annotation")

cc <- computeCommunProb(cc,
                        type            = "truncatedMean",
                        trim            = 0.001,
                        distance.use    = TRUE,
                        interaction.range = 100,     # secreted-signal diffusion range (µm)
                        contact.range   = 20,        # contact-dependent range (µm)
                        contact.dependent = TRUE,
                        scale.distance  = MPP,
                        nboot           = 20)

cc <- filterCommunication(cc, min.cells = 10)


# ---- Run every genotype through the identical pipeline ----
genotypes <- c(
  "Villin CreERT2(control)",
  "VillinCreERT2; BrafV600E/+",
  "Villin CreERT2; BrafV600E/V600E",
  "Villin CreERT2; Erk1/Erk2 DKO",
  "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO",
  "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"
)

cc_list <- list()
for (g in genotypes) {
  cc_list[[g]] <- run_geno(g)   # run_geno holds the parameter block above
  saveRDS(cc_list[[g]], paste0("cellchat_", gsub("[^A-Za-z0-9]","_",g), ".rds"))
  cat("done:", g, "\n")
}
saveRDS(cc_list, "cellchat_all_genotypes.rds")











library(CellChat); library(Seurat); library(dplyr)

# ---- reload inputs (in case session was restarted) ----
final <- read.csv("FINAL_all_celltypes.csv")
proj  <- read.csv("Default Pipeline Analysis-Spatial-Projection.csv")
keep_types <- c("Crypt","Junction","Villus","Villus_tip","Macrophage","Tcell","Fibroblast","DC")
final <- final[final$celltype %in% keep_types, ]
final$x <- proj$X.Coordinate[match(final$Barcode, proj$Barcode)]
final$y <- proj$Y.Coordinate[match(final$Barcode, proj$Barcode)]
final <- final[!is.na(final$x) & final$Barcode %in% colnames(mat2), ]

# ---- identical pipeline used for control ----
run_geno <- function(geno) {
  sub <- final[final$genotype == geno, ]
  cat("\n====", geno, "| cells:", nrow(sub), "====\n"); print(table(sub$celltype))
  bc <- sub$Barcode
  obj <- NormalizeData(CreateSeuratObject(mat2[, bc], assay="Spatial"))
  data.input <- GetAssayData(obj, layer="data")
  meta <- data.frame(labels=sub$celltype, samples=factor("s1"), row.names=bc)
  locs <- as.matrix(sub[, c("x","y")]); rownames(locs) <- bc
  
  MPP <- 0.3443589
  spatial.factors <- data.frame(ratio = MPP, tol = 4)
  
  cc <- createCellChat(object=data.input, meta=meta, group.by="labels",
                       datatype="spatial", coordinates=locs, spatial.factors=spatial.factors)
  cc@DB <- subsetDB(CellChatDB.mouse,
                    search=c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"), key="annotation")
  cc <- subsetData(cc)
  options(future.globals.maxSize=13000*1024^2)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc, variable.both=FALSE)
  cc <- computeCommunProb(cc, type="truncatedMean", trim=0.001,
                          distance.use=TRUE, interaction.range=100,
                          scale.distance=MPP, contact.dependent=TRUE,
                          contact.range=20, nboot=20)
  cc <- filterCommunication(cc, min.cells=10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc
}

# ---- run the 5 remaining genotypes ----
remaining <- c(
  "VillinCreERT2; BrafV600E/+",
  "Villin CreERT2; BrafV600E/V600E",
  "Villin CreERT2; Erk1/Erk2 DKO",
  "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO",
  "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"
)

for (g in remaining) {
  cc <- run_geno(g)
  fn <- paste0("cellchat_", gsub("[^A-Za-z0-9]","_",g), ".rds")
  saveRDS(cc, fn)
  cat("SAVED:", fn, "\n")
}
cat("\nAll remaining genotypes done.\n")

library(CellChat); library(dplyr)
cc_list <- readRDS("cellchat_all_genotypes.rds")

# ---- extract significant L-R interactions per genotype ----
sig_all <- lapply(names(cc_list), function(g){
  df <- subsetCommunication(cc_list[[g]])   # all inferred, already p<0.05 filtered
  df$genotype <- g
  df
})
sig_all <- do.call(rbind, sig_all)
write.csv(sig_all, "significant_interactions_all_genotypes.csv", row.names=FALSE)

# ---- focus: epithelium -> fibroblast/macrophage, key pathways ----
epi <- c("Crypt","Junction","Villus","Villus_tip")
key <- sig_all %>%
  filter(source %in% epi, target %in% c("Fibroblast","Macrophage")) %>%
  filter(pathway_name %in% c("BMP","TGFb","WNT","ncWNT","SPP1","PDGF","FGF"))

# ---- pathway strength per genotype (for the dose-response) ----
library(tidyr)
pw_by_geno <- key %>%
  group_by(genotype, pathway_name) %>%
  summarise(total_prob = sum(prob), n_pairs = n(), .groups="drop")
print(pw_by_geno)

# merge objects for CellChat's built-in comparison
cc_merged <- mergeCellChat(cc_list, add.names = names(cc_list))




library(CellChat); library(Seurat); library(dplyr)

final <- read.csv("FINAL_all_celltypes.csv")
proj  <- read.csv("Default Pipeline Analysis-Spatial-Projection.csv")

# NOW INCLUDE PLASMA
keep_types <- c("Crypt","Junction","Villus","Villus_tip",
                "Macrophage","Tcell","Fibroblast","DC","Plasma")
final <- final[final$celltype %in% keep_types, ]
final$x <- proj$X.Coordinate[match(final$Barcode, proj$Barcode)]
final$y <- proj$Y.Coordinate[match(final$Barcode, proj$Barcode)]
final <- final[!is.na(final$x) & final$Barcode %in% colnames(mat2), ]

run_geno <- function(geno) {
  sub <- final[final$genotype == geno, ]
  cat("\n====", geno, "| cells:", nrow(sub), "====\n"); print(table(sub$celltype))
  bc <- sub$Barcode
  obj <- NormalizeData(CreateSeuratObject(mat2[, bc], assay="Spatial"))
  data.input <- GetAssayData(obj, layer="data")
  meta <- data.frame(labels=sub$celltype, samples=factor("s1"), row.names=bc)
  locs <- as.matrix(sub[, c("x","y")]); rownames(locs) <- bc
  MPP <- 0.3443589
  spatial.factors <- data.frame(ratio=MPP, tol=4)
  cc <- createCellChat(object=data.input, meta=meta, group.by="labels",
                       datatype="spatial", coordinates=locs, spatial.factors=spatial.factors)
  cc@DB <- subsetDB(CellChatDB.mouse,
                    search=c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"), key="annotation")
  cc <- subsetData(cc)
  options(future.globals.maxSize=13000*1024^2)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc, variable.both=FALSE)
  cc <- computeCommunProb(cc, type="truncatedMean", trim=0.001,
                          distance.use=TRUE, interaction.range=100,
                          scale.distance=MPP, contact.dependent=TRUE,
                          contact.range=20, nboot=20)
  cc <- filterCommunication(cc, min.cells=10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc, slot.name="netP")   # for role analysis
  cc
}

genotypes <- c(
  "Villin CreERT2(control)","VillinCreERT2; BrafV600E/+",
  "Villin CreERT2; BrafV600E/V600E","Villin CreERT2; Erk1/Erk2 DKO",
  "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO",
  "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO")

cc_list <- list()
for (g in genotypes) {
  cc_list[[g]] <- run_geno(g)
  saveRDS(cc_list[[g]], paste0("cc_plasma_", gsub("[^A-Za-z0-9]","_",g), ".rds"))
  cat("SAVED:", g, "\n")
}
saveRDS(cc_list, "cc_plasma_all.rds")
cat("\nSTAGE 1 DONE\n")


library(CellChat); library(dplyr); library(tidyr)
cc_list <- readRDS("cc_plasma_all.rds")

# clean short names + dose order
short <- c("Villin CreERT2(control)"="Control",
           "VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo",
           "Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
names(cc_list) <- short[names(cc_list)]
dose_order <- c("Control","Braf_het","Braf_homo","ErkDKO","Braf_het_ErkDKO","Braf_homo_ErkDKO")

# ---- TABLE 1: total interactions & strength per genotype ----
t1 <- data.frame(
  genotype = names(cc_list),
  n_interactions = sapply(cc_list, function(c) sum(c@net$count)),
  total_strength = sapply(cc_list, function(c) sum(c@net$weight)),
  n_pathways = sapply(cc_list, function(c) length(c@netP$pathways)))
t1 <- t1[match(dose_order, t1$genotype),]
write.csv(t1, "T1_interactions_per_genotype.csv", row.names=FALSE)
cat("\n=== TABLE 1: interactions per genotype ===\n"); print(t1, row.names=FALSE)

# ---- TABLE 2: all significant L-R pairs, all genotypes ----
lr_all <- bind_rows(Map(function(cc,g){
  df <- subsetCommunication(cc); df$genotype <- g; df
}, cc_list, names(cc_list)))
write.csv(lr_all, "T2_all_LR_pairs.csv", row.names=FALSE)
cat("\n=== TABLE 2: significant L-R pairs (n rows) ===\n"); print(nrow(lr_all))

# ---- TABLE 3: pathway ranking by total strength ----
t3 <- lr_all %>% group_by(pathway_name) %>%
  summarise(total_prob=sum(prob), n_pairs=n(), n_genotypes=n_distinct(genotype),
            .groups="drop") %>% arrange(desc(total_prob))
write.csv(t3, "T3_pathway_ranking.csv", row.names=FALSE)
cat("\n=== TABLE 3: top pathways ===\n"); print(as.data.frame(head(t3,25)), row.names=FALSE)

# ---- TABLE 4: pathway strength x genotype (dose-response) ----
t4 <- lr_all %>% group_by(pathway_name, genotype) %>%
  summarise(prob=sum(prob), .groups="drop") %>%
  pivot_wider(names_from=genotype, values_from=prob, values_fill=0)
t4 <- t4[, c("pathway_name", intersect(dose_order, names(t4)))]
t4 <- t4[order(-rowSums(t4[,-1])),]
write.csv(t4, "T4_pathway_by_genotype.csv", row.names=FALSE)
cat("\n=== TABLE 4: pathway x genotype ===\n"); print(as.data.frame(head(t4,25)), row.names=FALSE)

# ---- TABLE 5: epithelium->stroma focus, top L-R ----
epi <- c("Crypt","Junction","Villus","Villus_tip")
stro <- c("Fibroblast","Macrophage","Plasma","Tcell","DC")
t5 <- lr_all %>%
  filter(source %in% epi, target %in% stro) %>%
  group_by(genotype, pathway_name, ligand, receptor) %>%
  summarise(prob=sum(prob), .groups="drop") %>% arrange(desc(prob))
write.csv(t5, "T5_epi_to_stroma_LR.csv", row.names=FALSE)
cat("\n=== TABLE 5: top epi->stroma L-R (head) ===\n"); print(as.data.frame(head(t5,25)), row.names=FALSE)

# ---- TABLE 6: top ligands & receptors overall ----
t6l <- lr_all %>% group_by(ligand) %>% summarise(prob=sum(prob), n=n(), .groups="drop") %>%
  arrange(desc(prob)) %>% head(20)
t6r <- lr_all %>% group_by(receptor) %>% summarise(prob=sum(prob), n=n(), .groups="drop") %>%
  arrange(desc(prob)) %>% head(20)
write.csv(t6l,"T6_top_ligands.csv",row.names=FALSE); write.csv(t6r,"T6_top_receptors.csv",row.names=FALSE)
cat("\n=== TOP LIGANDS ===\n"); print(as.data.frame(t6l),row.names=FALSE)
cat("\n=== TOP RECEPTORS ===\n"); print(as.data.frame(t6r),row.names=FALSE)
cat("\nSTAGE 2 DONE\n")







library(CellChat); library(ggplot2); library(patchwork)
cc_list <- readRDS("cc_plasma_all.rds")
names(cc_list) <- short[names(cc_list)]
cc_list <- cc_list[dose_order]   # enforce dose order

grp_cols <- c("Crypt"="#6A3D9A","Junction"="#F4C430","Villus"="#2C6FBB","Villus_tip"="#C0392B",
              "Fibroblast"="#1B9E77","Macrophage"="#D95F02","Tcell"="#7570B3",
              "DC"="#E7298A","Plasma"="#66A61E")

# ---- FIG 1: circle network per genotype (interaction strength), one page ----
pdf("F1_networks_per_genotype.pdf", width=15, height=10)
par(mfrow=c(2,3), xpd=TRUE, mar=c(1,1,3,1))
for(g in names(cc_list)){
  cc <- cc_list[[g]]
  gc <- grp_cols[levels(cc@idents)]
  netVisual_circle(cc@net$weight, vertex.weight=as.numeric(table(cc@idents)),
                   weight.scale=TRUE, label.edge=FALSE, color.use=gc, title.name=g)
}
dev.off()

# ---- FIG 2: merged comparison — interactions & strength per genotype ----
cc_merged <- mergeCellChat(cc_list, add.names=names(cc_list))
pdf("F2_compare_interactions.pdf", width=10, height=5)
p1 <- compareInteractions(cc_merged, show.legend=FALSE, group=1:length(cc_list), measure="count")
p2 <- compareInteractions(cc_merged, show.legend=FALSE, group=1:length(cc_list), measure="weight")
print(p1 + p2)
dev.off()

# ---- FIG 3: pathway dose-response (from Table 4) ----
t4 <- read.csv("T4_pathway_by_genotype.csv", check.names=FALSE)
key <- c("BMP","TGFb","WNT","ncWNT","SPP1","PDGF","FGF","GRN","COLLAGEN","LAMININ")
t4k <- t4[t4$pathway_name %in% key,]
long <- tidyr::pivot_longer(t4k, -pathway_name, names_to="genotype", values_to="prob")
long$genotype <- factor(long$genotype, levels=dose_order)
pdf("F3_pathway_dose_response.pdf", width=10, height=6)
ggplot(long, aes(genotype, prob, group=pathway_name, color=pathway_name)) +
  geom_line(linewidth=1) + geom_point(size=2) +
  labs(x=NULL, y="Total communication probability", color="Pathway",
       title="Signaling pathway strength across the allelic series") +
  theme_classic() + theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off()

# ---- FIG 4: differential interaction heatmap (control vs Braf_homo) ----
cc_pair <- mergeCellChat(list(Control=cc_list[["Control"]], Braf_homo=cc_list[["Braf_homo"]]),
                         add.names=c("Control","Braf_homo"))
pdf("F4_differential_heatmap.pdf", width=7, height=6)
netVisual_heatmap(cc_pair, measure="weight",
                  title.name="Differential interaction strength (Braf_homo vs Control)")
dev.off()

# ---- FIG 5: signaling role (sender/receiver) per genotype, key genotypes ----
pdf("F5_signaling_roles.pdf", width=12, height=5)
par(mfrow=c(1,3))
for(g in c("Control","Braf_homo","Braf_homo_ErkDKO")){
  netAnalysis_signalingRole_scatter(cc_list[[g]], title=g)
}
dev.off()

# ---- FIG 6: key pathway bubble across genotypes (BMP) ----
pdf("F6_BMP_across_genotypes.pdf", width=13, height=6)
netVisual_bubble(cc_merged,
                 sources.use=c("Crypt","Junction","Villus","Villus_tip"),
                 targets.use=c("Fibroblast","Macrophage","Plasma"),
                 signaling=intersect("BMP", unlist(lapply(cc_list,function(c)c@netP$pathways))),
                 comparison=1:length(cc_list), angle.x=45)
dev.off()

cat("\nSTAGE 3 DONE — figures F1-F6 saved\n")






library(CellChat); library(dplyr); library(tidyr)
cc_list <- readRDS("cc_plasma_all.rds")
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
names(cc_list) <- short[names(cc_list)]
dose_order <- c("Control","Braf_het","Braf_homo","ErkDKO","Braf_het_ErkDKO","Braf_homo_ErkDKO")

# all significant interactions, all genotypes
lr_all <- bind_rows(Map(function(cc,g){d<-subsetCommunication(cc); d$genotype<-g; d}, cc_list, names(cc_list)))

# ===== TABLE A: DATA-DRIVEN pathway ranking (how "top" is defined) =====
# rank = total communication probability across all cell-pairs & genotypes
tabA <- lr_all %>% group_by(pathway_name) %>%
  summarise(total_prob = round(sum(prob),4),
            n_LR_pairs = n(),
            mean_pval  = round(mean(pval),4),
            n_genotypes_present = n_distinct(genotype),
            .groups="drop") %>%
  arrange(desc(total_prob)) %>%
  mutate(rank = row_number())
write.csv(tabA, "TABLE_A_pathway_ranking_method.csv", row.names=FALSE)
cat("\n=== TABLE A: DATA-DRIVEN PATHWAY RANKING (top 20) ===\n")
print(as.data.frame(head(tabA,20)), row.names=FALSE)

# ===== TABLE B: TOP 10 pathways x genotype (dose-response) =====
top10_pw <- head(tabA$pathway_name, 10)
tabB <- lr_all %>% filter(pathway_name %in% top10_pw) %>%
  group_by(pathway_name, genotype) %>% summarise(prob=round(sum(prob),4),.groups="drop") %>%
  pivot_wider(names_from=genotype, values_from=prob, values_fill=0)
tabB <- tabB[, c("pathway_name", intersect(dose_order,names(tabB)))]
tabB <- tabB[match(top10_pw, tabB$pathway_name),]
write.csv(tabB,"TABLE_B_top10_by_genotype.csv",row.names=FALSE)
cat("\n=== TABLE B: TOP 10 PATHWAYS x GENOTYPE ===\n"); print(as.data.frame(tabB),row.names=FALSE)

# ===== TABLE C: HYPOTHESIZED pathways x genotype (BMP/WNT/FGF/etc) =====
hyp <- c("BMP","WNT","ncWNT","FGF","TGFb","SPP1","EGF","PDGF","NOTCH","IGF")
tabC <- lr_all %>% filter(pathway_name %in% hyp) %>%
  group_by(pathway_name, genotype) %>% summarise(prob=round(sum(prob),4),.groups="drop") %>%
  pivot_wider(names_from=genotype, values_from=prob, values_fill=0)
tabC <- tabC[, c("pathway_name", intersect(dose_order,names(tabC)))]
write.csv(tabC,"TABLE_C_hypothesized_by_genotype.csv",row.names=FALSE)
cat("\n=== TABLE C: HYPOTHESIZED PATHWAYS x GENOTYPE ===\n"); print(as.data.frame(tabC),row.names=FALSE)

# ===== FIGURE: heatmap of pathway x genotype (top + hypothesized) =====
library(ggplot2)
plot_pw <- unique(c(top10_pw, intersect(hyp, tabA$pathway_name)))
hm <- lr_all %>% filter(pathway_name %in% plot_pw) %>%
  group_by(pathway_name, genotype) %>% summarise(prob=sum(prob),.groups="drop")
# normalize each pathway to its own max (shows dose-trend, not absolute magnitude)
hm <- hm %>% group_by(pathway_name) %>% mutate(rel = prob/max(prob)) %>% ungroup()
hm$genotype <- factor(hm$genotype, levels=dose_order)
hm$pathway_name <- factor(hm$pathway_name, levels=rev(plot_pw))
pdf("FIG_pathway_heatmap.pdf", width=8, height=8)
ggplot(hm, aes(genotype, pathway_name, fill=rel)) +
  geom_tile(color="white") +
  scale_fill_gradient(low="#F7FBFF", high="#08306B", name="Relative\nstrength") +
  labs(x=NULL,y=NULL,title="Pathway signaling across the allelic series",
       subtitle="Each pathway normalized to its own maximum") +
  theme_minimal() + theme(axis.text.x=element_text(angle=45,hjust=1))
dev.off()
cat("\nDONE — tables A/B/C + heatmap\n")







library(CellChat); library(dplyr); library(tidyr); library(ggplot2)
cc_list <- readRDS("cc_plasma_all.rds")
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
names(cc_list) <- short[names(cc_list)]
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")
cc_list <- cc_list[dose_order]

lr_all <- bind_rows(Map(function(cc,g){d<-subsetCommunication(cc);d$genotype<-g;d}, cc_list, names(cc_list)))
lr_all$genotype <- factor(lr_all$genotype, levels=dose_order)

interesting <- intersect(c("BMP","WNT","ncWNT","FGF","TGFb","SPP1","EGF","PDGF","IGF","NOTCH"),
                         unique(lr_all$pathway_name))

# ---- TABLE: interesting pathways x genotype ----
tab <- lr_all %>% filter(pathway_name %in% interesting) %>%
  group_by(pathway_name, genotype) %>% summarise(prob=round(sum(prob),3),.groups="drop") %>%
  pivot_wider(names_from=genotype, values_from=prob, values_fill=0)
tab <- tab[, c("pathway_name", dose_order)]
write.csv(tab,"TABLE_interesting_ordered.csv",row.names=FALSE)
print(as.data.frame(tab), row.names=FALSE)

# ---- BAR GRAPHS per pathway (ordered) ----
bar_df <- lr_all %>% filter(pathway_name %in% interesting) %>%
  group_by(pathway_name, genotype) %>% summarise(prob=sum(prob),.groups="drop")
bar_df$genotype <- factor(bar_df$genotype, levels=dose_order)
geno_cols <- c("Braf_homo"="#B2182B","Braf_het"="#F4A582","Control"="grey60",
               "ErkDKO"="#92C5DE","Braf_homo_ErkDKO"="#2166AC","Braf_het_ErkDKO"="#4393C3")
pdf("FIG_bars_ordered.pdf", width=12, height=8)
ggplot(bar_df, aes(genotype, prob, fill=genotype)) +
  geom_col() + facet_wrap(~pathway_name, scales="free_y", ncol=3) +
  scale_fill_manual(values=geno_cols) +
  labs(x=NULL, y="Communication probability", title="Pathway strength across genotypes (MAPK-dose order)") +
  theme_bw() + theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="none",
                     strip.text=element_text(face="bold"))
dev.off()

# ---- HEATMAP (ordered, normalized per pathway) ----
hm <- bar_df %>% group_by(pathway_name) %>% mutate(rel=prob/max(prob)) %>% ungroup()
hm$genotype <- factor(hm$genotype, levels=dose_order)
hm$pathway_name <- factor(hm$pathway_name, levels=rev(interesting))
pdf("FIG_heatmap_ordered.pdf", width=7, height=6)
ggplot(hm, aes(genotype, pathway_name, fill=rel)) +
  geom_tile(color="white") + geom_text(aes(label=round(prob,2)), size=2.6) +
  scale_fill_gradient(low="#FFF7EC", high="#7F0000", name="Relative") +
  labs(x=NULL,y=NULL,title="Pathway signaling (MAPK-dose order)") +
  theme_minimal() + theme(axis.text.x=element_text(angle=45,hjust=1))
dev.off()

# ---- circle networks ordered ----
grp_cols <- c("Crypt"="#6A3D9A","Junction"="#F4C430","Villus"="#2C6FBB","Villus_tip"="#C0392B",
              "Fibroblast"="#1B9E77","Macrophage"="#D95F02","Tcell"="#7570B3","DC"="#E7298A","Plasma"="#66A61E")
pdf("FIG_networks_ordered.pdf", width=15, height=10)
par(mfrow=c(2,3), xpd=TRUE, mar=c(1,1,3,1))
for(g in dose_order){
  cc <- cc_list[[g]]; gc <- grp_cols[levels(cc@idents)]
  netVisual_circle(cc@net$weight, vertex.weight=as.numeric(table(cc@idents)),
                   weight.scale=TRUE, color.use=gc, title.name=g)
}
dev.off()
cat("DONE — ordered figures\n")


cc_merged <- mergeCellChat(cc_list, add.names=dose_order)
pdf("FIG_bubble_nature_style.pdf", width=14, height=7)
netVisual_bubble(cc_merged,
                 sources.use=c("Macrophage"),
                 targets.use=c("Tcell","Crypt","Villus","Villus_tip"),
                 comparison=1:length(cc_list),
                 angle.x=45, remove.isolate=TRUE)
dev.off()















library(CellChat); library(dplyr); library(tidyr); library(ggplot2)
cc_list <- readRDS("cc_plasma_all.rds")
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
names(cc_list) <- short[names(cc_list)]
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")
cc_list <- cc_list[dose_order]

lr_all <- bind_rows(Map(function(cc,g){d<-subsetCommunication(cc);d$genotype<-g;d}, cc_list, names(cc_list)))

# map your pathway names to CellChat DB names
# (CellChat labels: BMP, WNT/ncWNT, FGF, no direct Hippo, EGF, IL22? check, TGFb, NOTCH, NRG, no LXR)
wanted <- c(BMP="BMP", WNT="WNT", ncWNT="ncWNT", FGF="FGF", EGF="EGF",
            IL22="IL22", TGFb="TGFb", NOTCH="NOTCH", NRG="NRG")
# check what's actually present
present <- unique(lr_all$pathway_name)
cat("=== pathways present in YOUR data ===\n"); print(sort(present))
cat("\n=== of your wanted list, found: ===\n")
print(intersect(unname(wanted), present))
cat("\n=== NOT found (no signal / different name / absent): ===\n")
print(setdiff(unname(wanted), present))

# build table + bars for the ones present
myp <- intersect(unname(wanted), present)
tab <- lr_all %>% filter(pathway_name %in% myp) %>%
  group_by(pathway_name, genotype) %>% summarise(prob=round(sum(prob),3),.groups="drop") %>%
  pivot_wider(names_from=genotype, values_from=prob, values_fill=0)
tab <- tab[, c("pathway_name", dose_order)]
write.csv(tab,"TABLE_10pathways_ordered.csv",row.names=FALSE)
cat("\n=== YOUR 10-PATHWAY TABLE (MAPK-dose order) ===\n"); print(as.data.frame(tab),row.names=FALSE)

bar_df <- lr_all %>% filter(pathway_name %in% myp) %>%
  group_by(pathway_name, genotype) %>% summarise(prob=sum(prob),.groups="drop")
bar_df$genotype <- factor(bar_df$genotype, levels=dose_order)
geno_cols <- c("Braf_homo"="#B2182B","Braf_het"="#F4A582","Control"="grey60",
               "ErkDKO"="#92C5DE","Braf_homo_ErkDKO"="#2166AC","Braf_het_ErkDKO"="#4393C3")
pdf("FIG_10pathways_bars.pdf", width=13, height=8)
ggplot(bar_df, aes(genotype, prob, fill=genotype)) +
  geom_col() + facet_wrap(~pathway_name, scales="free_y", ncol=3) +
  scale_fill_manual(values=geno_cols) +
  labs(x=NULL,y="Communication probability",
       title="Signaling pathways across the allelic series (MAPK-dose order)") +
  theme_bw() + theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="none",
                     strip.text=element_text(face="bold", size=11))
dev.off()
cat("\nDONE\n")



lig <- c("Tgfb1","Wnt2b","Wnt5a","Nrg1","Dll1","Jag1","Areg","Ereg","Hbegf")  # EGF ligands too
# mean expr in Fibroblast + Macrophage per genotype












library(CellChat); library(dplyr); library(tidyr); library(ggplot2)
cc_list <- readRDS("cc_plasma_all.rds")
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
names(cc_list) <- short[names(cc_list)]
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")

lr_all <- bind_rows(Map(function(cc,g){d<-subsetCommunication(cc);d$genotype<-g;d}, cc_list, names(cc_list)))

# ===== restrict to FIBROBLAST or MACROPHAGE as SENDER, epithelium as receiver =====
epi <- c("Crypt","Junction","Villus","Villus_tip")
drivers <- c("Fibroblast","Macrophage")
fm <- lr_all %>% filter(source %in% drivers, target %in% epi)

# pathway strength per genotype (fibroblast/mac -> epithelium only)
pw <- fm %>% group_by(pathway_name, genotype) %>%
  summarise(prob=sum(prob), .groups="drop") %>%
  pivot_wider(names_from=genotype, values_from=prob, values_fill=0)

# ===== define Erk-dependence score =====
# Erk-DKO arms = ErkDKO, Braf_homo_ErkDKO, Braf_het_ErkDKO
# Erk-intact   = Braf_homo, Braf_het, Control
erk_arms <- c("ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")
int_arms <- c("Braf_homo","Braf_het","Control")
pw$mean_ErkDKO   <- rowMeans(pw[,erk_arms])
pw$mean_Erkintact<- rowMeans(pw[,int_arms])
pw$erk_effect    <- pw$mean_ErkDKO - pw$mean_Erkintact   # + = UP when Erk deleted
pw$fold          <- (pw$mean_ErkDKO+1e-6)/(pw$mean_Erkintact+1e-6)

# pathways most ELEVATED by Erk deletion (candidate ERK-suppressed niche signals)
erk_up <- pw %>% filter(mean_ErkDKO+mean_Erkintact > 0) %>%
  arrange(desc(erk_effect)) %>%
  select(pathway_name, all_of(dose_order), mean_Erkintact, mean_ErkDKO, erk_effect, fold)
write.csv(erk_up, "TABLE_erk_dependent_fibro_mac.csv", row.names=FALSE)
cat("\n=== FIBRO/MAC->EPI pathways MOST ELEVATED by Erk deletion ===\n")
print(as.data.frame(head(erk_up, 15)), row.names=FALSE)

# ===== TREND TEST across the ordered series (descriptive) =====
# Spearman correlation of pathway strength vs an MAPK-activity rank
# rank: Braf_homo=3(highest MAPK), Braf_het=2, Control=1, Erk arms=0 (ERK removed)
mapk_rank <- c(Braf_homo=3, Braf_het=2, Control=1, ErkDKO=0,
               Braf_homo_ErkDKO=0, Braf_het_ErkDKO=0)
trend <- fm %>% group_by(pathway_name, genotype) %>% summarise(prob=sum(prob),.groups="drop")
trend_stat <- trend %>% group_by(pathway_name) %>%
  filter(n()>=4) %>%
  summarise(rho = suppressWarnings(cor(prob, mapk_rank[genotype], method="spearman")),
            .groups="drop") %>%
  arrange(rho)   # most NEGATIVE rho = suppressed by MAPK / released by Erk-DKO
write.csv(trend_stat, "TABLE_mapk_trend_correlation.csv", row.names=FALSE)
cat("\n=== pathways NEGATIVELY correlated with MAPK activity (released by Erk loss) ===\n")
print(as.data.frame(head(trend_stat, 15)), row.names=FALSE)







# ===== FIG: top Erk-released niche pathways, bars across genotypes =====
top_erk <- head(erk_up$pathway_name, 9)
fig_df <- fm %>% filter(pathway_name %in% top_erk) %>%
  group_by(pathway_name, genotype) %>% summarise(prob=sum(prob),.groups="drop")
fig_df$genotype <- factor(fig_df$genotype, levels=dose_order)
geno_cols <- c("Braf_homo"="#B2182B","Braf_het"="#F4A582","Control"="grey60",
               "ErkDKO"="#92C5DE","Braf_homo_ErkDKO"="#2166AC","Braf_het_ErkDKO"="#4393C3")
pdf("FIG_erk_released_pathways.pdf", width=13, height=8)
ggplot(fig_df, aes(genotype, prob, fill=genotype)) +
  geom_col() + facet_wrap(~pathway_name, scales="free_y", ncol=3) +
  scale_fill_manual(values=geno_cols) +
  labs(x=NULL, y="Fibroblast/Macrophage -> Epithelium signaling",
       title="Niche pathways released by Erk1/2 deletion") +
  theme_bw() + theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="none",
                     strip.text=element_text(face="bold"))
dev.off()

# ===== FIG: Erk-intact vs Erk-DKO paired comparison (the mechanism in one plot) =====
comp <- erk_up %>% filter(mean_Erkintact+mean_ErkDKO > 0.002) %>%   # drop near-zero
  select(pathway_name, mean_Erkintact, mean_ErkDKO) %>%
  pivot_longer(-pathway_name, names_to="Erk_status", values_to="strength")
pdf("FIG_erk_mechanism_paired.pdf", width=7, height=6)
ggplot(comp, aes(Erk_status, strength, group=pathway_name, color=pathway_name)) +
  geom_line(linewidth=1) + geom_point(size=2.5) +
  scale_x_discrete(labels=c("mean_Erkintact"="Erk intact","mean_ErkDKO"="Erk DKO")) +
  labs(x=NULL, y="Fibro/Mac -> Epi signaling strength",
       title="Niche signaling is released by Erk1/2 deletion", color="Pathway") +
  theme_classic()
dev.off()
cat("\nDONE — mechanism analysis\n")

























library(CellChat); library(dplyr); library(tidyr); library(ggplot2)
cc_list <- readRDS("cc_plasma_all.rds")
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
names(cc_list) <- short[names(cc_list)]
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")
geno_cols <- c("Braf_homo"="#B2182B","Braf_het"="#F4A582","Control"="grey60",
               "ErkDKO"="#92C5DE","Braf_homo_ErkDKO"="#2166AC","Braf_het_ErkDKO"="#4393C3")

lr_all <- bind_rows(Map(function(cc,g){d<-subsetCommunication(cc);d$genotype<-g;d}, cc_list, names(cc_list)))

drivers <- c("Fibroblast","Macrophage","Tcell","Plasma")
zones   <- c("Villus_tip","Villus","Crypt","Junction")

# ---- rescue analysis: driver compartment -> single zone ----
analyze_pair <- function(driver, zone){
  fm <- lr_all %>% filter(source == driver, target == zone)
  if(nrow(fm)==0) return(NULL)
  pw <- fm %>% group_by(pathway_name, genotype) %>%
    summarise(prob=sum(prob), .groups="drop") %>%
    pivot_wider(names_from=genotype, values_from=prob, values_fill=0)
  for(g in dose_order) if(!g %in% names(pw)) pw[[g]] <- 0
  pw <- pw %>% mutate(
    braf_effect_homo = Braf_homo - Control,
    rescue_homo      = Braf_homo - Braf_homo_ErkDKO,
    braf_effect_het  = Braf_het  - Control,
    rescue_het       = Braf_het  - Braf_het_ErkDKO,
    driver = driver, zone = zone)
  pw %>% filter(braf_effect_homo>0, rescue_homo>0) %>% arrange(desc(rescue_homo))
}

# ---- run all driver x zone combos, collect ----
all_res <- list()
for(d in drivers) for(z in zones){
  r <- analyze_pair(d, z)
  if(!is.null(r) && nrow(r)>0) all_res[[paste(d,z,sep="_")]] <- r
}
combined <- bind_rows(all_res)
write.csv(combined, "FIGGIE_rescue_driver_x_zone.csv", row.names=FALSE)

cat("\n=== Braf-up & Erk-rescued signals, by DRIVER -> ZONE ===\n")
cat("(each row: a pathway where Braf raises signaling & Erk-DKO reverses it)\n\n")
summary_tbl <- combined %>%
  group_by(driver, zone) %>%
  summarise(n_rescued_pathways = n(),
            top_pathways = paste(head(pathway_name,4), collapse=", "),
            .groups="drop") %>%
  arrange(desc(n_rescued_pathways))
print(as.data.frame(summary_tbl), row.names=FALSE)

# ---- focus figure: signaling INTO Villus_tip (the expanded zone) ----
tip <- combined %>% filter(zone=="Villus_tip") %>% arrange(desc(rescue_homo))
cat("\n=== signals INTO VILLUS_TIP that are Braf-up / Erk-rescued ===\n")
print(as.data.frame(head(tip,15)), row.names=FALSE)
write.csv(tip, "FIGGIE_rescue_into_villustip.csv", row.names=FALSE)

# bar figure: top tip-targeting rescued pathways, per driver
if(nrow(tip)>0){
  top_tip <- head(unique(tip$pathway_name), 6)
  figdf <- lr_all %>%
    filter(target=="Villus_tip", source %in% drivers, pathway_name %in% top_tip) %>%
    group_by(pathway_name, source, genotype) %>% summarise(prob=sum(prob),.groups="drop")
  figdf$genotype <- factor(figdf$genotype, levels=dose_order)
  pdf("FIGGIE_villustip_rescue.pdf", width=14, height=9)
  print(ggplot(figdf, aes(genotype, prob, fill=genotype)) +
          geom_col() + facet_grid(source ~ pathway_name, scales="free_y") +
          scale_fill_manual(values=geno_cols) +
          labs(x=NULL, y="Signaling -> Villus_tip",
               title="Braf-induced, Erk-reversed signaling into Villus_tip, by source") +
          theme_bw() + theme(axis.text.x=element_text(angle=45,hjust=1,size=6),
                             legend.position="none", strip.text=element_text(face="bold",size=9)))
  dev.off()
}
cat("\nDONE — driver x zone rescue analysis\n")







library(CellChat); library(dplyr); library(tidyr); library(ggplot2)
cc_list <- readRDS("cc_plasma_all.rds")
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
names(cc_list) <- short[names(cc_list)]
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")
geno_cols <- c("Braf_homo"="#B2182B","Braf_het"="#F4A582","Control"="grey60",
               "ErkDKO"="#92C5DE","Braf_homo_ErkDKO"="#2166AC","Braf_het_ErkDKO"="#4393C3")

lr_all <- bind_rows(Map(function(cc,g){d<-subsetCommunication(cc);d$genotype<-g;d}, cc_list, names(cc_list)))

# ---- INSPECT: what's actually in the EGF pathway? ----
egf <- lr_all %>% filter(pathway_name=="EGF")
cat("=== EGF ligand-receptor pairs present ===\n")
print(egf %>% distinct(interaction_name, ligand, receptor))
cat("\n=== EGF sources (senders) ===\n"); print(sort(unique(egf$source)))
cat("\n=== EGF targets (receivers) ===\n"); print(sort(unique(egf$target)))
cat("\n=== EGF total by genotype ===\n")
print(egf %>% group_by(genotype) %>% summarise(prob=sum(prob), n=n()) %>%
        arrange(match(genotype,dose_order)))












library(CellChat); library(dplyr); library(tidyr); library(ggplot2)
egf$genotype <- factor(egf$genotype, levels=dose_order)

# ---- FIGGIE 1: EGF dose-response (headline) ----
f1 <- egf %>% group_by(genotype) %>% summarise(prob=sum(prob),.groups="drop")
f1$genotype <- factor(f1$genotype, levels=dose_order)
pdf("FIGGIE_EGF_1_doseresponse.pdf", width=7, height=5)
print(ggplot(f1, aes(genotype, prob, fill=genotype)) + geom_col(width=0.7) +
        scale_fill_manual(values=geno_cols) +
        labs(x=NULL,y="Total EGF communication probability",
             title="EGF/EGFR signaling scales with Braf dose, reversed by Erk1/2 deletion") +
        theme_classic() + theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="none"))
dev.off()

# ---- FIGGIE 2: WHO SENDS EGF (the mechanism-defining figure) ----
f2 <- egf %>% group_by(source, genotype) %>% summarise(prob=sum(prob),.groups="drop")
f2$genotype <- factor(f2$genotype, levels=dose_order)
src_cols <- c("Crypt"="#6A3D9A","Junction"="#F4C430","Villus"="#2C6FBB","Villus_tip"="#C0392B",
              "Fibroblast"="#1B9E77","Macrophage"="#D95F02","Tcell"="#7570B3","DC"="#E7298A","Plasma"="#66A61E")
pdf("FIGGIE_EGF_2_who_sends.pdf", width=9, height=5)
print(ggplot(f2, aes(genotype, prob, fill=source)) + geom_col(position="stack") +
        scale_fill_manual(values=src_cols, name="EGF sender") +
        labs(x=NULL,y="EGF communication probability",
             title="Source of EGF signaling across genotypes") +
        theme_classic() + theme(axis.text.x=element_text(angle=45,hjust=1)))
dev.off()

# ---- FIGGIE 3: ligand breakdown (which EGF ligand drives it) ----
f3 <- egf %>% group_by(ligand, genotype) %>% summarise(prob=sum(prob),.groups="drop")
f3$genotype <- factor(f3$genotype, levels=dose_order)
pdf("FIGGIE_EGF_3_ligands.pdf", width=11, height=6)
print(ggplot(f3, aes(genotype, prob, fill=genotype)) + geom_col() +
        facet_wrap(~ligand, scales="free_y") + scale_fill_manual(values=geno_cols) +
        labs(x=NULL,y="Communication probability", title="EGF-family ligands across genotypes") +
        theme_bw() + theme(axis.text.x=element_text(angle=45,hjust=1,size=6), legend.position="none",
                           strip.text=element_text(face="bold")))
dev.off()

# ---- FIGGIE 4: source->target heatmap, Control vs Braf_homo ----
f4 <- egf %>% filter(genotype %in% c("Control","Braf_homo")) %>%
  group_by(genotype, source, target) %>% summarise(prob=sum(prob),.groups="drop")
pdf("FIGGIE_EGF_4_source_target.pdf", width=11, height=5)
print(ggplot(f4, aes(target, source, fill=prob)) + geom_tile(color="white") +
        facet_wrap(~genotype) +
        scale_fill_gradient(low="#FFF5F0", high="#B2182B", name="EGF prob") +
        labs(x="Receiver", y="Sender", title="EGF source->target: Control vs Braf homozygous") +
        theme_minimal() + theme(axis.text.x=element_text(angle=45,hjust=1)))
dev.off()
cat("\nDONE — 4 EGF figures\n")

# print the sender breakdown so we can read the mechanism now
cat("\n=== EGF sent BY each cell type, Control vs Braf_homo ===\n")
print(egf %>% filter(genotype %in% c("Control","Braf_homo")) %>%
        group_by(genotype, source) %>% summarise(prob=round(sum(prob),4),.groups="drop") %>%
        pivot_wider(names_from=genotype, values_from=prob, values_fill=0) %>%
        arrange(desc(Braf_homo)))






egf_lig <- egf %>% group_by(ligand, genotype) %>% summarise(prob=round(sum(prob),4),.groups="drop") %>%
  pivot_wider(names_from=genotype, values_from=prob, values_fill=0)
egf_lig <- egf_lig[, c("ligand", dose_order)]
print(as.data.frame(egf_lig), row.names=FALSE)



egf %>% filter(ligand=="Egf") %>% group_by(source, genotype) %>%
  summarise(prob=round(sum(prob),4),.groups="drop") %>%
  pivot_wider(names_from=genotype, values_from=prob, values_fill=0)










library(CellChat); library(dplyr); library(tidyr); library(ggplot2)
cc_list <- readRDS("cc_plasma_all.rds")
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
names(cc_list) <- short[names(cc_list)]
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")
geno_cols <- c("Braf_homo"="#B2182B","Braf_het"="#F4A582","Control"="grey60",
               "ErkDKO"="#92C5DE","Braf_homo_ErkDKO"="#2166AC","Braf_het_ErkDKO"="#4393C3")

lr_all <- bind_rows(Map(function(cc,g){d<-subsetCommunication(cc);d$genotype<-g;d}, cc_list, names(cc_list)))

stroma <- c("Fibroblast","Macrophage","Plasma","Tcell","DC")
tip <- "Villus_tip"

# both directions: stroma <-> villus tip
vt <- lr_all %>% filter((source %in% stroma & target==tip) | (source==tip & target %in% stroma))
vt$direction <- ifelse(vt$target==tip, "to_tip", "from_tip")
vt$partner   <- ifelse(vt$target==tip, vt$source, vt$target)

# ---- TABLE 1: pathway x genotype, STROMA -> villus tip ----
to_tip <- vt %>% filter(direction=="to_tip") %>%
  group_by(pathway_name, genotype) %>% summarise(prob=sum(prob),.groups="drop") %>%
  pivot_wider(names_from=genotype, values_from=prob, values_fill=0)
to_tip <- to_tip[, c("pathway_name", dose_order)]
to_tip <- to_tip[order(-rowSums(to_tip[,-1])),]
write.csv(to_tip, "FIGGIE_stroma_to_tip.csv", row.names=FALSE)
cat("\n=== STROMA -> VILLUS_TIP pathways (top 20) ===\n")
print(as.data.frame(head(to_tip,20)), row.names=FALSE)

# ---- TABLE 2: which stromal cell sends most to tip, per genotype ----
by_sender <- vt %>% filter(direction=="to_tip") %>%
  group_by(partner, genotype) %>% summarise(prob=sum(prob),.groups="drop") %>%
  pivot_wider(names_from=genotype, values_from=prob, values_fill=0)
by_sender <- by_sender[, c("partner", dose_order)]
cat("\n=== signaling TO tip, by stromal sender x genotype ===\n")
print(as.data.frame(by_sender), row.names=FALSE)

# ---- rescue filter: Braf-up AND Erk-reversed, stroma->tip ----
resc <- to_tip %>%
  mutate(braf_eff = Braf_homo - Control, rescue = Braf_homo - Braf_homo_ErkDKO) %>%
  filter(braf_eff > 0, rescue > 0) %>% arrange(desc(rescue))
cat("\n=== STROMA->TIP pathways: Braf-up & Erk-rescued ===\n")
print(as.data.frame(head(resc,15)), row.names=FALSE)
write.csv(resc, "FIGGIE_stroma_to_tip_rescued.csv", row.names=FALSE)

# ---- FIGURE: top stroma->tip pathways across genotypes ----
top_pw <- head(to_tip$pathway_name, 9)
fig <- vt %>% filter(direction=="to_tip", pathway_name %in% top_pw) %>%
  group_by(pathway_name, genotype) %>% summarise(prob=sum(prob),.groups="drop")
fig$genotype <- factor(fig$genotype, levels=dose_order)
pdf("FIGGIE_stroma_to_tip_bars.pdf", width=13, height=8)
print(ggplot(fig, aes(genotype, prob, fill=genotype)) +
        geom_col() + facet_wrap(~pathway_name, scales="free_y", ncol=3) +
        scale_fill_manual(values=geno_cols) +
        labs(x=NULL, y="Signaling -> Villus_tip",
             title="Top stromal->villus tip signaling pathways across genotypes") +
        theme_bw() + theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="none",
                           strip.text=element_text(face="bold")))
dev.off()

# ---- FIGURE: heatmap of stromal sender x pathway INTO tip (Braf_homo) ----
hm <- vt %>% filter(direction=="to_tip", genotype=="Braf_homo", pathway_name %in% top_pw) %>%
  group_by(partner, pathway_name) %>% summarise(prob=sum(prob),.groups="drop")
pdf("FIGGIE_tip_sender_pathway_heatmap.pdf", width=8, height=5)
print(ggplot(hm, aes(pathway_name, partner, fill=prob)) +
        geom_tile(color="white") +
        scale_fill_gradient(low="#FFF5F0", high="#B2182B", name="prob") +
        labs(x="Pathway", y="Stromal sender", title="Stromal->villus tip signaling in Braf homozygous") +
        theme_minimal() + theme(axis.text.x=element_text(angle=45,hjust=1)))
dev.off()

cat("\nDONE\n")















library(CellChat); library(dplyr); library(tidyr); library(ggplot2)
cc_list <- readRDS("cc_plasma_all.rds")
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
names(cc_list) <- short[names(cc_list)]
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")
geno_cols <- c("Braf_homo"="#B2182B","Braf_het"="#F4A582","Control"="grey60",
               "ErkDKO"="#92C5DE","Braf_homo_ErkDKO"="#2166AC","Braf_het_ErkDKO"="#4393C3")
lr_all <- bind_rows(Map(function(cc,g){d<-subsetCommunication(cc);d$genotype<-g;d}, cc_list, names(cc_list)))

# ---- figure for ONE directed pair: sender -> receiver ----
pair_fig <- function(sender, receiver){
  d <- lr_all %>% filter(source==sender, target==receiver)
  if(nrow(d)==0){ cat("\n",sender,"->",receiver,": no interactions\n"); return(invisible()) }
  pw <- d %>% group_by(pathway_name, genotype) %>% summarise(prob=sum(prob),.groups="drop") %>%
    pivot_wider(names_from=genotype, values_from=prob, values_fill=0)
  for(g in dose_order) if(!g %in% names(pw)) pw[[g]] <- 0
  pw <- pw[, c("pathway_name", dose_order)]; pw <- pw[order(-rowSums(pw[,-1])),]
  cat("\n===============", sender, "->", receiver, "===============\n")
  print(as.data.frame(head(pw,12)), row.names=FALSE)
  write.csv(pw, paste0("FIGGIE_", sender, "_to_", receiver, ".csv"), row.names=FALSE)
  
  top <- head(pw$pathway_name, 9)
  fig <- d %>% filter(pathway_name %in% top) %>%
    group_by(pathway_name, genotype) %>% summarise(prob=sum(prob),.groups="drop")
  fig$genotype <- factor(fig$genotype, levels=dose_order)
  pdf(paste0("FIGGIE_", sender, "_to_", receiver, ".pdf"), width=13, height=8)
  print(ggplot(fig, aes(genotype, prob, fill=genotype)) +
          geom_col() + facet_wrap(~pathway_name, scales="free_y", ncol=3) +
          scale_fill_manual(values=geno_cols) +
          labs(x=NULL, y="Communication probability",
               title=paste0(sender, " \u2192 ", receiver, " signaling across genotypes")) +
          theme_bw() + theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="none",
                             strip.text=element_text(face="bold")))
  dev.off()
  invisible(pw)
}

# ---- run the stroma -> villus tip pairs, each separately ----
pair_fig("Fibroblast", "Villus_tip")
pair_fig("Macrophage", "Villus_tip")
pair_fig("Plasma",     "Villus_tip")
pair_fig("Tcell",      "Villus_tip")
pair_fig("DC",         "Villus_tip")

cat("\nDONE — one figure per directed cell-type pair\n")







library(Seurat); library(dplyr); library(tidyr); library(ggplot2)
final <- read.csv("FINAL_all_celltypes.csv")
cells <- final$Barcode[final$Barcode %in% colnames(mat2)]
obj <- CreateSeuratObject(mat2[, cells], assay="Spatial"); obj <- NormalizeData(obj)
obj$genotype <- final$genotype[match(colnames(obj), final$Barcode)]
obj$celltype <- final$celltype[match(colnames(obj), final$Barcode)]
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
obj$geno <- short[obj$genotype]
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")

# modules for each hypothesis
modules <- list(
  Purinergic = c("Nt5e","Ada","Adora1","Adora2a","Adora2b","Entpd1","Entpd2","Panx1"),
  AKT_survival = c("Akt3","Akt1","Akt2","Pik3ca","Pdk1","Mtor","Bcl2","Bcl2l1","Xiap"),
  TightJunction = c("Cldn4","Cldn3","Cldn7","Cldn2","Tjp1","Ocln","F11r","Cgn"),
  BMP_target = c("Id1","Id2","Id3","Smad6","Smad7","Bambi","Msx1","Msx2")  # BMP pathway ACTIVITY readout
)
for(m in names(modules)) modules[[m]] <- modules[[m]][modules[[m]] %in% rownames(obj)]
obj <- AddModuleScore(obj, features=modules, name=names(modules))
# rename score columns
sc <- grep("^(Purinergic|AKT_survival|TightJunction|BMP_target)[0-9]", colnames(obj@meta.data), value=TRUE)

# tip only, module scores per genotype
tip <- subset(obj, celltype=="Villus_tip")
res <- lapply(seq_along(modules), function(i){
  col <- paste0(names(modules)[i], i)
  data.frame(module=names(modules)[i],
             tip@meta.data %>% mutate(geno=factor(geno,levels=dose_order)) %>%
               group_by(geno) %>% summarise(mean=round(mean(.data[[col]]),3),.groups="drop"))
})
cat("=== MODULE SCORES in VILLUS_TIP across genotypes ===\n")
for(r in res){ cat("\n--", r$module[1], "--\n"); print(r %>% select(geno,mean), row.names=FALSE) }

# individual key genes in tip
key <- c("Nt5e","Ada","Akt3","Cldn4")
key <- key[key %in% rownames(obj)]
Idents(tip) <- "geno"
cat("\n=== KEY TIP GENES (raw mean expr) x genotype ===\n")
print(round(as.data.frame(as.matrix(AverageExpression(tip, features=key, assays="Spatial")$Spatial))[,dose_order],3))












library(Seurat); library(dplyr)
# (obj already built: normalized, celltype + geno assigned)
epi_machinery <- list(
  Polycomb    = c("Ezh2","Eed","Suz12","Rnf2","Cbx7"),
  DNAmeth     = c("Dnmt1","Dnmt3a","Dnmt3b","Uhrf1"),
  DNAdemeth   = c("Tet1","Tet2","Tet3"),
  H3K27demeth = c("Kdm6a","Kdm6b"),
  HDAC        = c("Hdac1","Hdac2","Hdac3"),
  Acetyl      = c("Ep300","Crebbp","Kat2b")
)
for(m in names(epi_machinery)) epi_machinery[[m]] <- epi_machinery[[m]][epi_machinery[[m]] %in% rownames(obj)]
obj <- AddModuleScore(obj, features=epi_machinery, name=names(epi_machinery))

tip <- subset(obj, celltype=="Villus_tip")
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")
sc_cols <- grep(paste0("^(",paste(names(epi_machinery),collapse="|"),")[0-9]"), colnames(tip@meta.data), value=TRUE)
cat("=== epigenetic-machinery module scores in VILLUS_TIP by genotype ===\n")
for(i in seq_along(epi_machinery)){
  col <- paste0(names(epi_machinery)[i], i)
  cat("\n--", names(epi_machinery)[i], "--\n")
  print(tip@meta.data %>% mutate(geno=factor(geno,levels=dose_order)) %>%
          group_by(geno) %>% summarise(mean=round(mean(.data[[col]]),3),.groups="drop"), row.names=FALSE)
}
# key individual genes
key <- intersect(c("Ezh2","Dnmt3b","Dnmt1","Tet1","Kdm6b"), rownames(obj))
Idents(tip) <- "geno"
cat("\n=== key epigenetic genes (raw) x genotype ===\n")
print(round(as.data.frame(as.matrix(AverageExpression(tip, features=key, assays="Spatial")$Spatial))[,dose_order],3))




library(Seurat); library(dplyr)

# ---- rebuild obj cleanly from scratch ----
final <- read.csv("FINAL_all_celltypes.csv")
cells <- final$Barcode[final$Barcode %in% colnames(mat2)]
obj <- CreateSeuratObject(mat2[, cells], assay="Spatial")
obj <- NormalizeData(obj)

# attach metadata by matching barcodes to colnames(obj) explicitly
m <- final[match(colnames(obj), final$Barcode), ]
obj$celltype <- m$celltype
obj$genotype <- m$genotype
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
obj$geno <- short[obj$genotype]
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")

cat("cells:", ncol(obj), " | tip cells:", sum(obj$celltype=="Villus_tip"), "\n")

# ---- epigenetic machinery modules ----
epi_machinery <- list(
  Polycomb    = c("Ezh2","Eed","Suz12","Rnf2","Cbx7"),
  DNAmeth     = c("Dnmt1","Dnmt3a","Dnmt3b","Uhrf1"),
  DNAdemeth   = c("Tet1","Tet2","Tet3"),
  H3K27demeth = c("Kdm6a","Kdm6b"),
  HDAC        = c("Hdac1","Hdac2","Hdac3"),
  Acetyl      = c("Ep300","Crebbp","Kat2b")
)
epi_machinery <- lapply(epi_machinery, function(g) g[g %in% rownames(obj)])
epi_machinery <- epi_machinery[sapply(epi_machinery, length) > 0]  # drop empty
obj <- AddModuleScore(obj, features=epi_machinery, name=names(epi_machinery))

# ---- report scores in villus tip ----
tip <- subset(obj, celltype=="Villus_tip")
cat("\n=== epigenetic-machinery module scores in VILLUS_TIP by genotype ===\n")
for(i in seq_along(epi_machinery)){
  col <- paste0(names(epi_machinery)[i], i)
  if(col %in% colnames(tip@meta.data)){
    cat("\n--", names(epi_machinery)[i], "(",paste(epi_machinery[[i]],collapse=","),") --\n")
    print(tip@meta.data %>% mutate(geno=factor(geno,levels=dose_order)) %>%
            group_by(geno) %>% summarise(mean=round(mean(.data[[col]]),3),.groups="drop") %>%
            as.data.frame(), row.names=FALSE)
  }
}

# ---- key individual genes, raw ----
key <- intersect(c("Ezh2","Dnmt3b","Dnmt1","Tet1","Tet2","Kdm6b","Uhrf1"), rownames(obj))
if(length(key)>0){
  Idents(tip) <- "geno"
  cat("\n=== key epigenetic genes (raw mean expr) x genotype ===\n")
  av <- AverageExpression(tip, features=key, assays="Spatial")$Spatial
  print(round(as.data.frame(as.matrix(av))[, dose_order], 3))
}




library(Seurat); library(dplyr)

final <- read.csv("FINAL_all_celltypes.csv")
cells <- final$Barcode[final$Barcode %in% colnames(mat2)]
obj <- CreateSeuratObject(mat2[, cells], assay="Spatial")
obj <- NormalizeData(obj)

m <- final[match(colnames(obj), final$Barcode), ]
obj$celltype <- m$celltype
obj$genotype <- m$genotype

short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")

# unname() strips the vector names so Seurat assigns by position, not by name
obj$geno <- unname(short[obj$genotype])
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")

cat("cells:", ncol(obj), " | tip cells:", sum(obj$celltype=="Villus_tip", na.rm=TRUE), "\n")
cat("geno values:", paste(names(table(obj$geno)), collapse=", "), "\n")

# ---- epigenetic machinery modules ----
epi_machinery <- list(
  Polycomb    = c("Ezh2","Eed","Suz12","Rnf2","Cbx7"),
  DNAmeth     = c("Dnmt1","Dnmt3a","Dnmt3b","Uhrf1"),
  DNAdemeth   = c("Tet1","Tet2","Tet3"),
  H3K27demeth = c("Kdm6a","Kdm6b"),
  HDAC        = c("Hdac1","Hdac2","Hdac3"),
  Acetyl      = c("Ep300","Crebbp","Kat2b")
)
epi_machinery <- lapply(epi_machinery, function(g) g[g %in% rownames(obj)])
epi_machinery <- epi_machinery[sapply(epi_machinery, length) > 0]
# report which genes were found
cat("\ngenes found per module:\n")
for(nm in names(epi_machinery)) cat(" ", nm, ":", paste(epi_machinery[[nm]], collapse=", "), "\n")

obj <- AddModuleScore(obj, features=epi_machinery, name=names(epi_machinery))

tip <- subset(obj, celltype=="Villus_tip")
cat("\n=== epigenetic-machinery module scores in VILLUS_TIP by genotype ===\n")
for(i in seq_along(epi_machinery)){
  col <- paste0(names(epi_machinery)[i], i)
  if(col %in% colnames(tip@meta.data)){
    cat("\n--", names(epi_machinery)[i], "--\n")
    d <- tip@meta.data
    d$geno <- factor(d$geno, levels=dose_order)
    print(d %>% group_by(geno) %>% summarise(mean=round(mean(.data[[col]]),3),.groups="drop") %>%
            as.data.frame(), row.names=FALSE)
  }
}

key <- intersect(c("Ezh2","Dnmt3b","Dnmt1","Tet1","Tet2","Uhrf1"), rownames(obj))
if(length(key)>0){
  Idents(tip) <- "geno"
  cat("\n=== key epigenetic genes (raw) x genotype ===\n")
  av <- AverageExpression(tip, features=key, assays="Spatial")$Spatial
  print(round(as.data.frame(as.matrix(av))[, dose_order], 3))
  
  
  
  
  
  
  
  
  
  
}





library(Seurat); library(dplyr); library(tidyr)
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")

# ---- hypothesis-specific gene modules ----
mods <- list(
  # A: methyl-donor / SAM-SAH machinery
  MethylCycle = c("Ahcy","Mat2a","Mat2b","Mtr","Mthfr","Ada","Nt5e","Adk"),
  # B: ERK->EZH2 / Polycomb
  Polycomb    = c("Ezh2","Eed","Suz12"),
  # C: YAP targets vs BMP targets (score separately)
  YAP_targets = c("Ctgf","Cyr61","Amotl2","Ankrd1","Ccn1","Ccn2"),
  BMP_targets = c("Id1","Id2","Id3","Smad6","Smad7","Bambi"),
  # D: anoikis / survival
  Anoikis_resist = c("Akt3","Bcl2l1","Xiap","Birc5"),
  Proapoptotic   = c("Bcl2l11","Bmf","Bax","Casp3")
)
mods <- lapply(mods, function(g) g[g %in% rownames(obj)])
mods <- mods[sapply(mods,length)>0]
cat("genes found per module:\n"); for(n in names(mods)) cat(" ",n,":",paste(mods[[n]],collapse=", "),"\n")

obj <- AddModuleScore(obj, features=mods, name=names(mods))

# ---- report each module in VILLUS_TIP across genotypes ----
tip <- subset(obj, celltype=="Villus_tip")
d <- tip@meta.data; d$geno <- factor(d$geno, levels=dose_order)
cat("\n===== MODULE SCORES IN VILLUS_TIP (by genotype, dose order) =====\n")
for(i in seq_along(mods)){
  col <- paste0(names(mods)[i], i)
  if(col %in% colnames(d)){
    cat("\n--", names(mods)[i], "--\n")
    print(d %>% group_by(geno) %>% summarise(mean=round(mean(.data[[col]]),3),.groups="drop") %>%
            as.data.frame(), row.names=FALSE)
  }
}

# ---- key individual genes raw, in tip ----
key <- intersect(c("Akt3","Ezh2","Ahcy","Ada","Nt5e","Ctgf","Cyr61","Id1","Bcl2l11","Cldn4"), rownames(obj))
Idents(tip) <- "geno"
cat("\n===== KEY GENES (raw mean expr) IN TIP x genotype =====\n")
av <- AverageExpression(tip, features=key, assays="Spatial")$Spatial
print(round(as.data.frame(as.matrix(av))[, dose_order], 3))

# ---- Hypothesis C test: YAP/BMP ratio (divergence) ----
yc <- paste0("YAP_targets", which(names(mods)=="YAP_targets"))
bc <- paste0("BMP_targets", which(names(mods)=="BMP_targets"))
if(yc %in% colnames(d) & bc %in% colnames(d)){
  cat("\n===== YAP vs BMP target activity in tip (Hypothesis C) =====\n")
  print(d %>% group_by(geno) %>%
          summarise(YAP=round(mean(.data[[yc]]),3), BMP=round(mean(.data[[bc]]),3),
                    YAP_minus_BMP=round(mean(.data[[yc]])-mean(.data[[bc]]),3), .groups="drop") %>%
          as.data.frame(), row.names=FALSE)
}





library(Seurat); library(dplyr)
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")

key <- intersect(c("Akt3","Ezh2","Ahcy","Ada","Nt5e","Ctgf","Cyr61","Id1","Bcl2l11","Cldn4"), rownames(obj))
Idents(tip) <- "geno"
av <- AverageExpression(tip, features=key, assays="Spatial")$Spatial
av <- as.data.frame(as.matrix(av))

cat("columns actually returned:\n"); print(colnames(av))

# reorder to dose_order but only keep columns that exist
cols <- intersect(dose_order, colnames(av))
cat("\n===== KEY GENES (raw mean expr) IN TIP =====\n")
print(round(av[, cols, drop=FALSE], 3))

library(Seurat); library(dplyr); library(tidyr)
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")
key <- intersect(c("Akt3","Ezh2","Ahcy","Ada","Nt5e","Ctgf","Cyr61","Id1","Bcl2l11","Cldn4","Bmp8a"), rownames(obj))

# pull expression + geno directly, compute means by hand
expr <- FetchData(tip, vars = c(key, "geno"))
cat("geno levels in fetched data:\n"); print(table(expr$geno))

tab <- expr %>%
  group_by(geno) %>%
  summarise(across(all_of(key), ~round(mean(.x),3)), .groups="drop")

# transpose to genes x genotypes, ordered
tab_t <- as.data.frame(t(tab[,-1]))
colnames(tab_t) <- tab$geno
cols <- intersect(dose_order, colnames(tab_t))
cat("\n===== KEY GENES (mean expr) IN TIP, all genotypes =====\n")
print(tab_t[, cols, drop=FALSE])



key2 <- intersect(c("Atf3","Il22ra1","Stat3","Socs3","Cflar","Reg3b","Reg3g"), rownames(obj))
# Il22ra1 = IL-22 receptor, Socs3/Reg3 = STAT3 targets, Cflar = the ATF3 death-suppression target
expr2 <- FetchData(tip, vars=c(key2,"geno"))
expr2 %>% group_by(geno) %>% summarise(across(all_of(key2), ~round(mean(.x),3))) %>%
  as.data.frame() %>% {rownames(.)<-.$geno; t(.[,-1])} %>% .[, dose_order]






library(Seurat); library(dplyr); library(tidyr)

# ---- reload the count matrix ----
mat2 <- Read10X_h5("filtered_feature_cell_matrix.h5")   # the segmented matrix
# if that errors, use whatever load line you originally used for mat2

# ---- rebuild the annotated Seurat object ----
final <- read.csv("FINAL_all_celltypes.csv")
cells <- final$Barcode[final$Barcode %in% colnames(mat2)]
obj <- CreateSeuratObject(mat2[, cells], assay="Spatial")
obj <- NormalizeData(obj)

# attach metadata by barcode
m <- final[match(colnames(obj), final$Barcode), ]
obj$celltype <- m$celltype
obj$genotype <- m$genotype
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
obj$geno <- unname(short[obj$genotype])
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")

# tip subset
tip <- subset(obj, celltype=="Villus_tip")
cat("tip cells per genotype:\n"); print(table(tip$geno))




setwd("C:/Users/write/OneDrive/Documents/MapKSeries_Segmented_MouseDuodenum")
list.files(pattern="\\.h5$")   # confirm the h5 file is here


library(Seurat); library(dplyr); library(tidyr)

mat2 <- Read10X_h5("filtered_feature_cell_matrix.h5")

final <- read.csv("FINAL_all_celltypes.csv")
cells <- final$Barcode[final$Barcode %in% colnames(mat2)]
obj <- CreateSeuratObject(mat2[, cells], assay="Spatial")
obj <- NormalizeData(obj)

m <- final[match(colnames(obj), final$Barcode), ]
obj$celltype <- m$celltype
obj$genotype <- m$genotype
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
obj$geno <- unname(short[obj$genotype])
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")

tip <- subset(obj, celltype=="Villus_tip")
cat("tip cells per genotype:\n"); print(table(tip$geno))

# save so a restart never costs you this again
saveRDS(obj, "obj_annotated.rds")



library(Seurat); library(dplyr); library(tidyr)
# obj + tip already loaded (obj_annotated.rds); if not: obj <- readRDS("obj_annotated.rds"); tip <- subset(obj, celltype=="Villus_tip")
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")

# the mechanism's genes:
#  ATF3 = pro-apoptotic brake (predict DOWN with Braf, restored by ErkDKO)
#  pro-death effectors ATF3 drives (predict DOWN with Braf)
#  anti-death (cFLIP/survival) (predict UP with Braf) -- the reciprocal arm
mech <- intersect(c(
  "Atf3",                          # the brake
  "Ddit3","Gdf15","Bbc3","Bax","Bcl2l11",  # pro-apoptotic (CHOP=Ddit3, NAG-1=Gdf15, PUMA=Bbc3, Bim=Bcl2l11)
  "Cflar","Bcl2l1","Birc5"         # anti-apoptotic / survival (cFLIP, Bcl-xL, survivin)
), rownames(obj))

expr <- FetchData(tip, vars=c(mech,"geno"))
tab <- expr %>% group_by(geno) %>% summarise(across(all_of(mech), ~round(mean(.x),3)),.groups="drop")
tab_t <- as.data.frame(t(tab[,-1])); colnames(tab_t) <- tab$geno
tab_t <- tab_t[, intersect(dose_order, colnames(tab_t)), drop=FALSE]

# annotate each gene's role + whether the direction fits the mechanism
role <- c(Atf3="brake(down?)", Ddit3="pro-death", Gdf15="pro-death", Bbc3="pro-death",
          Bax="pro-death", Bcl2l11="pro-death", Cflar="survival(up?)",
          Bcl2l1="survival", Birc5="survival")
tab_t$role <- role[rownames(tab_t)]

cat("===== ATF3 MECHANISM — villus tip, MAPK-dose order =====\n")
cat("Prediction: Atf3 & pro-death genes DOWN in Braf_homo; survival genes UP; all reverse in ErkDKO arms\n\n")
print(tab_t)










library(CellChat); library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
cc_list <- readRDS("cc_plasma_all.rds")
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
names(cc_list) <- short[names(cc_list)]
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")
cc_list <- cc_list[dose_order]
# compute centrality for all (needed for role scatter)
cc_list <- lapply(cc_list, function(cc) netAnalysis_computeCentrality(cc, slot.name="netP"))

# ============ PANEL G: signaling role scatter, per genotype ============
pdf("PANEL_G_role_scatter.pdf", width=15, height=10)
library(patchwork)
plots <- lapply(names(cc_list), function(g)
  netAnalysis_signalingRole_scatter(cc_list[[g]], title=g) )
print(wrap_plots(plots, ncol=3))
dev.off()

# ============ PANEL H: differential interaction heatmap (pairwise) ============
# pick the key contrast: Braf_homo vs Control
cc_pair <- mergeCellChat(list(Control=cc_list[["Control"]], Braf_homo=cc_list[["Braf_homo"]]),
                         add.names=c("Control","Braf_homo"))
pdf("PANEL_H_diff_interaction.pdf", width=7, height=6)
netVisual_heatmap(cc_pair, measure="weight",
                  title.name="Differential interaction strength (Braf_homo vs Control)")
dev.off()

# ============ PANEL I: differential pathway scatter (Braf_homo vs Control) ============
# outgoing vs incoming signaling strength per pathway, two conditions
get_pw_strength <- function(cc){
  cc <- netAnalysis_computeCentrality(cc, slot.name="netP")
  out <- rowSums(sapply(cc@netP$pathways, function(p){
    i <- which(cc@netP$pathways==p); c(out=sum(cc@netP$prob[,,i]))
  }))
  data.frame(pathway=cc@netP$pathways,
             outgoing=sapply(seq_along(cc@netP$pathways), function(i) sum(cc@netP$prob[,,i][,])),
             incoming=sapply(seq_along(cc@netP$pathways), function(i) sum(cc@netP$prob[,,i])))
}
# simpler: total strength per pathway per condition, compare
pw_ctrl <- sapply(seq_along(cc_list[["Control"]]@netP$pathways), function(i) sum(cc_list[["Control"]]@netP$prob[,,i]))
names(pw_ctrl) <- cc_list[["Control"]]@netP$pathways
pw_homo <- sapply(seq_along(cc_list[["Braf_homo"]]@netP$pathways), function(i) sum(cc_list[["Braf_homo"]]@netP$prob[,,i]))
names(pw_homo) <- cc_list[["Braf_homo"]]@netP$pathways
allpw <- union(names(pw_ctrl), names(pw_homo))
dfI <- data.frame(pathway=allpw,
                  Control=pw_ctrl[allpw], Braf_homo=pw_homo[allpw])
dfI[is.na(dfI)] <- 0
dfI$diff <- dfI$Braf_homo - dfI$Control
pdf("PANEL_I_pathway_scatter.pdf", width=7, height=7)
print(ggplot(dfI, aes(Control, Braf_homo, label=pathway)) +
        geom_abline(slope=1, linetype="dashed", color="grey") +
        geom_point(aes(color=diff), size=2) +
        scale_color_gradient2(low="#2166AC", mid="grey80", high="#B2182B", name="Δ") +
        ggrepel::geom_text_repel(data=subset(dfI, abs(diff)>quantile(abs(dfI$diff),0.8)), size=3) +
        labs(title="Pathway signaling: Braf_homo vs Control",
             x="Control signaling strength", y="Braf_homo signaling strength") +
        theme_classic())
dev.off()

# ============ PANEL J: relative signaling strength heatmap (pathway x genotype) ============
lr_all <- bind_rows(Map(function(cc,g){d<-subsetCommunication(cc);d$genotype<-g;d}, cc_list, names(cc_list)))
top_pw <- lr_all %>% group_by(pathway_name) %>% summarise(s=sum(prob),.groups="drop") %>%
  arrange(desc(s)) %>% head(18) %>% pull(pathway_name)
hmJ <- lr_all %>% filter(pathway_name %in% top_pw) %>%
  group_by(pathway_name, genotype) %>% summarise(prob=sum(prob),.groups="drop") %>%
  group_by(pathway_name) %>% mutate(rel=prob/max(prob)) %>% ungroup()
hmJ$genotype <- factor(hmJ$genotype, levels=dose_order)
hmJ$pathway_name <- factor(hmJ$pathway_name, levels=rev(top_pw))
pdf("PANEL_J_relative_strength.pdf", width=7, height=7)
print(ggplot(hmJ, aes(genotype, pathway_name, fill=rel)) +
        geom_tile(color="white") +
        scale_fill_gradient(low="#FFF7EC", high="#7F0000", name="Relative\nstrength") +
        labs(x=NULL, y=NULL, title="Relative signaling strength across genotypes") +
        theme_minimal() + theme(axis.text.x=element_text(angle=45,hjust=1)))
dev.off()

cat("panels G, H, I, J saved\n")









library(CellChat); library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
cc_list <- readRDS("cc_plasma_all.rds")
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
names(cc_list) <- short[names(cc_list)]
cc_list <- lapply(cc_list, function(cc) netAnalysis_computeCentrality(cc, slot.name="netP"))

# ---- function to make the full panel set for ONE pairwise contrast ----
make_contrast <- function(condA, condB, tag){
  # condA = Braf (mutant), condB = Braf+ErkDKO (rescue); positive = up in Braf
  pair <- mergeCellChat(list(cc_list[[condA]], cc_list[[condB]]), add.names=c(condA, condB))
  
  # PANEL 1: differential interaction heatmap (sender x receiver)
  pdf(paste0("CONTRAST_",tag,"_1_diff_heatmap.pdf"), width=7, height=6)
  print(netVisual_heatmap(pair, measure="weight",
                          title.name=paste0("Diff interaction strength: ",condA," vs ",condB)))
  dev.off()
  
  # PANEL 2: differential interaction COUNT circle
  pdf(paste0("CONTRAST_",tag,"_2_diff_circle.pdf"), width=8, height=8)
  netVisual_diffInteraction(pair, weight.scale=TRUE, measure="weight",
                            title.name=paste0(condA," vs ",condB))
  dev.off()
  
  # PANEL 3: role scatter side by side
  pdf(paste0("CONTRAST_",tag,"_3_role_scatter.pdf"), width=11, height=5)
  p1 <- netAnalysis_signalingRole_scatter(cc_list[[condA]], title=condA)
  p2 <- netAnalysis_signalingRole_scatter(cc_list[[condB]], title=condB)
  print(p1 + p2)
  dev.off()
  
  # PANEL 4: pathway strength scatter (A vs B), labeled
  pwA <- setNames(sapply(seq_along(cc_list[[condA]]@netP$pathways),
                         function(i) sum(cc_list[[condA]]@netP$prob[,,i])), cc_list[[condA]]@netP$pathways)
  pwB <- setNames(sapply(seq_along(cc_list[[condB]]@netP$pathways),
                         function(i) sum(cc_list[[condB]]@netP$prob[,,i])), cc_list[[condB]]@netP$pathways)
  allpw <- union(names(pwA), names(pwB))
  df <- data.frame(pathway=allpw, A=pwA[allpw], B=pwB[allpw]); df[is.na(df)] <- 0
  df$diff <- df$A - df$B
  pdf(paste0("CONTRAST_",tag,"_4_pathway_scatter.pdf"), width=7, height=7)
  print(ggplot(df, aes(B, A, label=pathway)) +
          geom_abline(slope=1, linetype="dashed", color="grey60") +
          geom_point(aes(color=diff), size=2.5) +
          scale_color_gradient2(low="#2166AC", mid="grey85", high="#B2182B", name="Δ (Braf-rescue)") +
          ggrepel::geom_text_repel(data=subset(df, abs(diff) > quantile(abs(df$diff),0.75)), size=3, max.overlaps=20) +
          labs(title=paste0("Pathway signaling: ",condA," (y) vs ",condB," (x)"),
               subtitle="above diagonal = higher in Braf mutant (rescued by ErkDKO)",
               x=paste0(condB," strength"), y=paste0(condA," strength")) +
          theme_classic())
  dev.off()
  
  # PANEL 5: which pathways change most (bar of the difference)
  df2 <- df[order(-df$diff),]
  df2 <- rbind(head(df2,12), tail(df2,8))   # top up + top down
  df2$pathway <- factor(df2$pathway, levels=df2$pathway)
  pdf(paste0("CONTRAST_",tag,"_5_pathway_diff_bar.pdf"), width=7, height=7)
  print(ggplot(df2, aes(diff, pathway, fill=diff>0)) + geom_col() +
          scale_fill_manual(values=c("TRUE"="#B2182B","FALSE"="#2166AC"),
                            labels=c("TRUE"="up in Braf (rescued)","FALSE"="up in ErkDKO"), name=NULL) +
          labs(x=paste0("Δ signaling (",condA," − ",condB,")"), y=NULL,
               title=paste0("Pathways changed by Erk rescue: ",tag)) +
          theme_classic())
  dev.off()
  
  cat("contrast", tag, "done: 5 panels\n")
  df  # return the pathway diff table
}

# ---- run BOTH contrasts ----
homo_diff <- make_contrast("Braf_homo", "Braf_homo_ErkDKO", "HOMO")
het_diff  <- make_contrast("Braf_het",  "Braf_het_ErkDKO",  "HET")

# ---- compare the two rescues: is the same pathway rescued at both doses? ----
merged <- merge(homo_diff[,c("pathway","diff")], het_diff[,c("pathway","diff")],
                by="pathway", suffixes=c("_homo","_het"), all=TRUE)
merged[is.na(merged)] <- 0
merged <- merged[order(-(merged$diff_homo + merged$diff_het)),]
write.csv(merged, "CONTRAST_rescue_both_doses.csv", row.names=FALSE)
cat("\n=== pathways rescued by ErkDKO at BOTH Braf doses (top 15) ===\n")
print(head(merged, 15), row.names=FALSE)

# figure: rescue at homo vs rescue at het (concordance)
pdf("CONTRAST_rescue_concordance.pdf", width=7, height=7)
print(ggplot(merged, aes(diff_het, diff_homo, label=pathway)) +
        geom_hline(yintercept=0, color="grey70") + geom_vline(xintercept=0, color="grey70") +
        geom_point(color="#B2182B", size=2) +
        ggrepel::geom_text_repel(data=subset(merged, abs(diff_homo)+abs(diff_het) >
                                               quantile(abs(merged$diff_homo)+abs(merged$diff_het),0.8)), size=3, max.overlaps=20) +
        labs(title="Erk rescue concordance across Braf doses",
             subtitle="upper-right = rescued (Braf-up) at both doses = robust mechanism",
             x="Δ at Braf_het (het − het+ErkDKO)", y="Δ at Braf_homo (homo − homo+ErkDKO)") +
        theme_classic())
dev.off()
cat("\nALL DONE\n")








library(CellChat); library(dplyr); library(tidyr); library(ggplot2)
cc_list <- readRDS("cc_plasma_all.rds")
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
names(cc_list) <- short[names(cc_list)]
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")
cc_list <- cc_list[dose_order]

# ---- build pathway -> category map from CellChatDB ----
db <- CellChatDB.mouse$interaction
# each pathway's annotation category: "Secreted Signaling","ECM-Receptor","Cell-Cell Contact"
pw_cat <- db %>% distinct(pathway_name, annotation) %>%
  mutate(category = case_when(
    annotation=="Secreted Signaling" ~ "Secreted",
    annotation=="ECM-Receptor"       ~ "ECM",
    annotation=="Cell-Cell Contact"  ~ "Non-protein/Contact",
    TRUE ~ annotation))
# some pathways appear in >1; take first
pw_cat <- pw_cat %>% group_by(pathway_name) %>% slice(1) %>% ungroup()

lr_all <- bind_rows(Map(function(cc,g){d<-subsetCommunication(cc);d$genotype<-g;d}, cc_list, names(cc_list)))
lr_all$category <- pw_cat$category[match(lr_all$pathway_name, pw_cat$pathway_name)]
lr_all$category[is.na(lr_all$category)] <- "Other"

cat_order <- c("Secreted","ECM","Non-protein/Contact")

# ============ PANEL I: pathway scatter, FACETED by category ============
# contrast Braf_homo vs Braf_homo_ErkDKO (rescue)
condA <- "Braf_homo"; condB <- "Braf_homo_ErkDKO"
pwstr <- lr_all %>% filter(genotype %in% c(condA,condB)) %>%
  group_by(pathway_name, category, genotype) %>% summarise(prob=sum(prob),.groups="drop") %>%
  pivot_wider(names_from=genotype, values_from=prob, values_fill=0)
names(pwstr)[names(pwstr)==condA] <- "A"; names(pwstr)[names(pwstr)==condB] <- "B"
pwstr$diff <- pwstr$A - pwstr$B
pwstr$category <- factor(pwstr$category, levels=cat_order)
pwstr <- pwstr[!is.na(pwstr$category),]

pdf("PANEL_I_scatter_by_category.pdf", width=6, height=11)
print(ggplot(pwstr, aes(B, A, label=pathway_name)) +
        geom_abline(slope=1, linetype="dashed", color="grey60") +
        geom_point(aes(color=diff), size=2.5) +
        scale_color_gradient2(low="#2166AC", mid="grey85", high="#B2182B", name="Δ") +
        ggrepel::geom_text_repel(size=2.8, max.overlaps=15) +
        facet_wrap(~category, ncol=1, scales="free") +
        labs(title=paste0("Pathway signaling by category: ",condA," vs ",condB),
             subtitle="above diagonal = higher in Braf (rescued by ErkDKO)",
             x=paste0(condB," strength"), y=paste0(condA," strength")) +
        theme_bw() + theme(strip.text=element_text(face="bold", size=11)))
dev.off()

# ============ PANEL J: relative strength heatmap, Y-AXIS GROUPED by category ============
hmJ <- lr_all %>%
  group_by(pathway_name, category, genotype) %>% summarise(prob=sum(prob),.groups="drop") %>%
  group_by(pathway_name) %>% mutate(rel=prob/max(prob)) %>% ungroup()
# keep top pathways per category so it's readable
top_per_cat <- hmJ %>% group_by(category, pathway_name) %>% summarise(s=sum(prob),.groups="drop") %>%
  group_by(category) %>% slice_max(s, n=8) %>% ungroup()
hmJ <- hmJ %>% filter(pathway_name %in% top_per_cat$pathway_name)
hmJ$category <- factor(hmJ$category, levels=cat_order)
hmJ$genotype <- factor(hmJ$genotype, levels=dose_order)
# order pathways within category by total strength
pw_ord <- hmJ %>% group_by(category, pathway_name) %>% summarise(s=sum(prob),.groups="drop") %>%
  arrange(category, s) %>% pull(pathway_name)
hmJ$pathway_name <- factor(hmJ$pathway_name, levels=unique(pw_ord))

pdf("PANEL_J_heatmap_by_category.pdf", width=7, height=9)
print(ggplot(hmJ, aes(genotype, pathway_name, fill=rel)) +
        geom_tile(color="white") +
        facet_grid(category ~ ., scales="free_y", space="free_y") +
        scale_fill_gradient(low="#FFF7EC", high="#7F0000", name="Relative\nstrength") +
        labs(x=NULL, y=NULL, title="Relative signaling strength by category") +
        theme_minimal() +
        theme(axis.text.x=element_text(angle=45,hjust=1),
              strip.text.y=element_text(face="bold", angle=0),
              panel.spacing=unit(0.4,"lines")))
dev.off()

cat("done — panels I and J split by Secreted / ECM / Non-protein\n")
# print the category assignments so you can verify
cat("\npathway categories:\n")
print(pw_cat %>% arrange(category, pathway_name) %>% as.data.frame())



library(Seurat)
obj <- readRDS("obj_annotated.rds")
proj <- read.csv("Default Pipeline Analysis-Spatial-Projection.csv")
md <- data.frame(Barcode=colnames(obj), celltype=obj$celltype, geno=obj$geno)
md$x <- proj$X.Coordinate[match(md$Barcode, proj$Barcode)]
md$y <- proj$Y.Coordinate[match(md$Barcode, proj$Barcode)]
md2 <- md %>% filter(celltype %in% c("Tcell","Villus_tip"), !is.na(x))
md2$geno <- factor(md2$geno, levels=dose_order)
pdf("Tcell_tip_spatial.pdf", width=14, height=9)
print(ggplot(md2, aes(x, y, color=celltype)) +
        geom_point(size=0.3, alpha=0.6) +
        scale_color_manual(values=c("Tcell"="#7570B3","Villus_tip"="#C0392B")) +
        facet_wrap(~geno, ncol=3) + coord_fixed() +
        labs(title="Spatial distribution: T cells and villus tip across genotypes") +
        theme_void() + theme(legend.position="bottom"))
dev.off()





library(CellChat); library(dplyr); library(tidyr)
cc_list <- readRDS("cc_plasma_all.rds")
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
names(cc_list) <- short[names(cc_list)]
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")
cc_list <- cc_list[dose_order]

lr_all <- bind_rows(Map(function(cc,g){d<-subsetCommunication(cc);d$genotype<-g;d}, cc_list, names(cc_list)))

# the regeneration panel (map to CellChat DB names)
regen <- c("TGFb","MIF","NRG","COMPLEMENT","GUCA","EGF","VEGF","ANGPTL",
           "IGF","TNF","GALECTIN","COLLAGEN","TENASCIN","PROSTAGLANDIN","RA")
# check which are present, and print CellChat's actual naming
present <- sort(unique(lr_all$pathway_name))
cat("=== regeneration pathways FOUND in your data ===\n")
print(intersect(regen, present))
cat("\n=== requested but NOT found (may have different name) ===\n")
print(setdiff(regen, present))
cat("\n=== all pathways containing related terms (to catch renamed ones) ===\n")
print(grep("RA|RETINOIC|PROSTAG|PTGS|GUCA|GUCY", present, value=TRUE, ignore.case=TRUE))




regen_present <- c("TGFb","MIF","NRG","COMPLEMENT","GUCA","EGF","VEGF","ANGPTL",
                   "IGF","TNF","GALECTIN","COLLAGEN","TENASCIN")



# ===== ANGLE 1: total strength per pathway per genotype (dose-response) =====
A1 <- lr_all %>% filter(pathway_name %in% regen_present) %>%
  group_by(pathway_name, genotype) %>% summarise(prob=round(sum(prob),4),.groups="drop") %>%
  pivot_wider(names_from=genotype, values_from=prob, values_fill=0)
A1 <- A1[, c("pathway_name", dose_order)]
# add trend metrics
A1$braf_effect <- A1$Braf_homo - A1$Control                    # + = Braf raises
A1$erk_rescue  <- A1$Braf_homo - A1$Braf_homo_ErkDKO           # + = ErkDKO reverses
A1$erkdko_alone<- A1$ErkDKO - A1$Control                        # ERK deletion effect w/o Braf
cat("\n===== ANGLE 1: pathway strength x genotype + trend metrics =====\n")
print(as.data.frame(A1), row.names=FALSE)

# ===== ANGLE 2: classify each pathway's behavior =====
A2 <- A1 %>% mutate(
  pattern = case_when(
    braf_effect > 0.001 & erk_rescue > 0.001 ~ "Braf-UP / Erk-rescued",
    braf_effect < -0.001 & erk_rescue < -0.001 ~ "Braf-DOWN / Erk-restored",
    erkdko_alone > 0.001 & braf_effect < 0.001 ~ "ERK-suppressed (up in DKO)",
    TRUE ~ "flat/other")
) %>% select(pathway_name, braf_effect, erk_rescue, erkdko_alone, pattern)
cat("\n===== ANGLE 2: behavior classification =====\n")
print(as.data.frame(A2[order(A2$pattern, -A2$braf_effect),]), row.names=FALSE)

# ===== ANGLE 3: WHO sends these pathways (dominant source per pathway) =====
A3 <- lr_all %>% filter(pathway_name %in% regen_present) %>%
  group_by(pathway_name, source) %>% summarise(prob=sum(prob),.groups="drop") %>%
  group_by(pathway_name) %>% slice_max(prob, n=2) %>%
  summarise(top_sources=paste(source, collapse=", "),.groups="drop")
cat("\n===== ANGLE 3: dominant SENDER cell types per pathway =====\n")
print(as.data.frame(A3), row.names=FALSE)

# ===== ANGLE 4: WHO receives (dominant target per pathway) =====
A4 <- lr_all %>% filter(pathway_name %in% regen_present) %>%
  group_by(pathway_name, target) %>% summarise(prob=sum(prob),.groups="drop") %>%
  group_by(pathway_name) %>% slice_max(prob, n=2) %>%
  summarise(top_targets=paste(target, collapse=", "),.groups="drop")
cat("\n===== ANGLE 4: dominant RECEIVER cell types per pathway =====\n")
print(as.data.frame(A4), row.names=FALSE)

# ===== ANGLE 5: signaling specifically INTO villus tip =====
A5 <- lr_all %>% filter(pathway_name %in% regen_present, target=="Villus_tip") %>%
  group_by(pathway_name, genotype) %>% summarise(prob=round(sum(prob),4),.groups="drop") %>%
  pivot_wider(names_from=genotype, values_from=prob, values_fill=0)
if(nrow(A5)>0){ A5 <- A5[, c("pathway_name", intersect(dose_order,names(A5)))]
cat("\n===== ANGLE 5: regeneration pathways INTO villus tip =====\n")
print(as.data.frame(A5), row.names=FALSE) }

# ===== ANGLE 6: a "regeneration score" — sum of all regen pathways per genotype =====
A6 <- lr_all %>% filter(pathway_name %in% regen_present) %>%
  group_by(genotype) %>% summarise(total_regen=round(sum(prob),3), n_interactions=n(),.groups="drop")
A6 <- A6[match(dose_order, A6$genotype),]
cat("\n===== ANGLE 6: TOTAL regeneration signaling per genotype =====\n")
print(as.data.frame(A6), row.names=FALSE)




library(CellChat); library(Seurat); library(dplyr); library(tidyr)
cc_list <- readRDS("cc_plasma_all.rds")
obj <- readRDS("obj_annotated.rds")
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
names(cc_list) <- short[names(cc_list)]
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")
cc_list <- cc_list[dose_order]
lr_all <- bind_rows(Map(function(cc,g){d<-subsetCommunication(cc);d$genotype<-g;d}, cc_list, names(cc_list)))

# ===== 1. GUCA ligand-receptor pairs =====
guca <- lr_all %>% filter(pathway_name=="GUCA")
cat("=== GUCA L-R pairs ===\n")
print(guca %>% distinct(interaction_name, ligand, receptor))

# ===== 2. RAW EXPRESSION of GUCA genes across genotypes (capture check) =====
guca_genes <- intersect(c("Guca2a","Guca2b","Gucy2c"), rownames(obj))
cat("\nGUCA genes present:", paste(guca_genes, collapse=", "), "\n")
Idents(obj) <- "geno"
cat("\n=== raw mean expression, WHOLE tissue x genotype ===\n")
av <- AverageExpression(obj, features=guca_genes, assays="Spatial")$Spatial
avdf <- as.data.frame(as.matrix(av)); cols <- intersect(dose_order, colnames(avdf))
print(round(avdf[, cols, drop=FALSE], 3))

# % of cells expressing (capture reliability)
cat("\n=== % of ALL cells expressing each GUCA gene ===\n")
ex <- FetchData(obj, vars=guca_genes)
print(sapply(ex, function(x) round(100*mean(x>0),1)))

# ===== 3. GUCA genes in VILLUS TIP specifically, across genotypes =====
tip <- subset(obj, celltype=="Villus_tip")
Idents(tip) <- "geno"
cat("\n=== GUCA genes in VILLUS_TIP x genotype (raw) ===\n")
avt <- as.data.frame(as.matrix(AverageExpression(tip, features=guca_genes, assays="Spatial")$Spatial))
colst <- intersect(dose_order, colnames(avt))
print(round(avt[, colst, drop=FALSE], 3))

# ===== 4. GUCA by cell type (who expresses the ligands vs receptor) =====
Idents(obj) <- "celltype"
cat("\n=== GUCA genes by cell type (which cells make ligand vs receptor) ===\n")
avc <- as.data.frame(as.matrix(AverageExpression(obj, features=guca_genes, assays="Spatial")$Spatial))
print(round(avc, 3))

# ===== 5. GUCA signaling into tip, per genotype (from CellChat) =====
cat("\n=== GUCA CellChat signaling INTO villus tip x genotype ===\n")
print(guca %>% filter(target=="Villus_tip") %>%
        group_by(genotype, source) %>% summarise(prob=round(sum(prob),4),.groups="drop") %>%
        pivot_wider(names_from=genotype, values_from=prob, values_fill=0) %>% as.data.frame())




library(dplyr); library(tidyr)
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")
guca_genes <- c("Guca2a","Guca2b","Gucy2c")

# whole tissue, all six genotypes
e_all <- FetchData(obj, vars=c(guca_genes,"geno"))
cat("=== GUCA genes, WHOLE tissue, all genotypes ===\n")
t_all <- e_all %>% group_by(geno) %>% summarise(across(all_of(guca_genes), ~round(mean(.x),3)),.groups="drop")
tt <- as.data.frame(t(t_all[,-1])); colnames(tt) <- t_all$geno
print(tt[, intersect(dose_order, colnames(tt))])

# villus tip only, all six genotypes
e_tip <- FetchData(tip, vars=c(guca_genes,"geno"))
cat("\n=== GUCA genes, VILLUS TIP, all genotypes ===\n")
t_tip <- e_tip %>% group_by(geno) %>% summarise(across(all_of(guca_genes), ~round(mean(.x),3)),.groups="drop")
tt2 <- as.data.frame(t(t_tip[,-1])); colnames(tt2) <- t_tip$geno
print(tt2[, intersect(dose_order, colnames(tt2))])


for(ct in c("Crypt","Villus","Junction","Villus_tip")){
  sub <- subset(obj, celltype==ct)
  e <- FetchData(sub, vars=c("Guca2a","Guca2b","Gucy2c","geno"))
  cat("\n---", ct, "---\n")
  t <- e %>% group_by(geno) %>% summarise(across(c(Guca2a,Guca2b,Gucy2c), ~round(mean(.x),3)),.groups="drop")
  tt <- as.data.frame(t(t[,-1])); colnames(tt) <- t$geno
  print(tt[, intersect(dose_order, colnames(tt))])
}







library(Seurat); library(dplyr); library(tidyr); library(ggplot2)
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")
geno_cols <- c("Braf_homo"="#B2182B","Braf_het"="#F4A582","Control"="grey60",
               "ErkDKO"="#92C5DE","Braf_homo_ErkDKO"="#2166AC","Braf_het_ErkDKO"="#4393C3")
guca_genes <- c("Guca2a","Guca2b","Gucy2c")

# gather expression by zone x genotype
zones <- c("Crypt","Junction","Villus","Villus_tip")
dat <- lapply(zones, function(z){
  sub <- subset(obj, celltype==z)
  e <- FetchData(sub, vars=c(guca_genes,"geno"))
  e %>% group_by(geno) %>% summarise(across(all_of(guca_genes), ~mean(.x)),.groups="drop") %>%
    mutate(zone=z)
}) %>% bind_rows() %>%
  pivot_longer(all_of(guca_genes), names_to="gene", values_to="expr")
dat$geno <- factor(dat$geno, levels=dose_order)
dat$zone <- factor(dat$zone, levels=zones)

# ---- FIGURE 1: Guca2a across zones and genotypes (the clean dose-response) ----
pdf("GUCA_fig1_Guca2a_zones.pdf", width=11, height=4)
print(ggplot(subset(dat, gene=="Guca2a"), aes(geno, expr, fill=geno)) +
        geom_col() + facet_wrap(~zone, nrow=1) +
        scale_fill_manual(values=geno_cols) +
        labs(x=NULL, y="Guca2a (guanylin) expression",
             title="Guanylin (Guca2a) induced by Braf in differentiated zones, reversed by Erk deletion") +
        theme_bw() + theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="none",
                           strip.text=element_text(face="bold")))
dev.off()

# ---- FIGURE 2: all three GUCA genes in the villus tip ----
pdf("GUCA_fig2_tip_allgenes.pdf", width=10, height=4)
print(ggplot(subset(dat, zone=="Villus_tip"), aes(geno, expr, fill=geno)) +
        geom_col() + facet_wrap(~gene, scales="free_y", nrow=1) +
        scale_fill_manual(values=geno_cols) +
        labs(x=NULL, y="Expression",
             title="Guanylin-GUCY2C axis in the villus tip across the allelic series") +
        theme_bw() + theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="none",
                           strip.text=element_text(face="bold")))
dev.off()

# ---- FIGURE 3: heatmap, gene x zone in Braf_homo vs Control (spatial induction) ----
hm <- dat %>% filter(geno %in% c("Control","Braf_homo")) %>%
  select(gene, zone, geno, expr) %>%
  pivot_wider(names_from=geno, values_from=expr) %>%
  mutate(log2FC = log2((Braf_homo+0.01)/(Control+0.01)))
pdf("GUCA_fig3_induction_heatmap.pdf", width=6, height=4)
print(ggplot(hm, aes(zone, gene, fill=log2FC)) +
        geom_tile(color="white") +
        geom_text(aes(label=round(log2FC,1)), size=3) +
        scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B", name="log2FC\nBraf/Ctrl") +
        labs(x=NULL, y=NULL, title="GUCA induction (Braf_homo vs Control) by zone") +
        theme_minimal() + theme(axis.text.x=element_text(angle=45,hjust=1)))
dev.off()
cat("GUCA figures saved\n")





library(Seurat); library(dplyr); library(tidyr)
obj <- readRDS("obj_annotated.rds")   # if not loaded
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")

# SREBP-2 / mevalonate-cholesterol biosynthesis program
srebp <- list(mevalonate = c(
  "Srebf2",           # SREBP-2 itself (master TF)
  "Hmgcr","Hmgcs1",   # rate-limiting mevalonate enzymes (statin target = Hmgcr)
  "Mvk","Mvd","Pmvk", # mevalonate kinase pathway
  "Fdps","Fdft1",     # farnesyl/squalene
  "Sqle",             # squalene epoxidase
  "Lss","Cyp51",      # lanosterol
  "Dhcr7","Dhcr24",   # distal cholesterol synthesis
  "Idi1","Insig1",    # isopentenyl / SREBP regulation
  "Ldlr"))            # LDL receptor (SREBP-2 target)
srebp$mevalonate <- srebp$mevalonate[srebp$mevalonate %in% rownames(obj)]
cat("mevalonate genes found:", paste(srebp$mevalonate, collapse=", "), "\n")

# crypt cells only
crypt <- subset(obj, celltype=="Crypt")
cat("crypt cells per genotype:\n"); print(table(crypt$geno))

# ---- module score in crypt across genotypes ----
obj <- AddModuleScore(obj, features=srebp, name="Mevalonate")
crypt <- subset(obj, celltype=="Crypt")
e <- FetchData(crypt, vars=c("Mevalonate1","geno"))
cat("\n=== SREBP-2/mevalonate MODULE SCORE in crypt x genotype ===\n")
print(e %>% group_by(geno) %>% summarise(mean=round(mean(Mevalonate1),3),.groups="drop") %>%
        arrange(match(geno,dose_order)) %>% as.data.frame(), row.names=FALSE)

# ---- key individual genes (raw) in crypt, all genotypes ----
key <- srebp$mevalonate
ec <- FetchData(crypt, vars=c(key,"geno"))
tab <- ec %>% group_by(geno) %>% summarise(across(all_of(key), ~round(mean(.x),3)),.groups="drop")
tt <- as.data.frame(t(tab[,-1])); colnames(tt) <- tab$geno
cat("\n=== key mevalonate genes (raw) in crypt x genotype ===\n")
print(tt[, intersect(dose_order, colnames(tt)), drop=FALSE])



library(ggplot2)
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")
geno_cols <- c("Braf_homo"="#B2182B","Braf_het"="#F4A582","Control"="grey60",
               "ErkDKO"="#92C5DE","Braf_homo_ErkDKO"="#2166AC","Braf_het_ErkDKO"="#4393C3")
e$geno <- factor(e$geno, levels=dose_order)
pdf("SREBP2_crypt_module.pdf", width=6, height=5)
print(ggplot(e, aes(geno, Mevalonate1, fill=geno)) +
        geom_violin(scale="width", alpha=0.5) +
        stat_summary(fun=mean, geom="crossbar", width=0.5) +
        scale_fill_manual(values=geno_cols) +
        labs(x=NULL, y="SREBP-2/mevalonate module score",
             title="Cholesterol-biosynthesis program in Braf crypts (dose-dependent, Erk-reversed)") +
        theme_classic() + theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="none"))
dev.off()





library(CellChat); library(dplyr); library(tidyr); library(ggplot2)

# Ensure cc_list and dose_order are loaded in your environment
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")
geno_cols <- c("Braf_homo"="#B2182B","Braf_het"="#F4A582","Control"="grey60",
               "ErkDKO"="#92C5DE","Braf_homo_ErkDKO"="#2166AC","Braf_het_ErkDKO"="#4393C3")

# 1. Bind all interactions into one master dataframe
lr_all <- bind_rows(Map(function(cc, g) {
  d <- subsetCommunication(cc)
  d$genotype <- g
  d
}, cc_list, names(cc_list)))
lr_all$genotype <- factor(lr_all$genotype, levels = dose_order)

# 2. Extract Data for Each Hypothesis

# Hypothesis A: Raft Sensitization (INCOMING to Crypt)
hypA_paths <- c("WNT", "ncWNT", "EGF", "NOTCH")
df_A <- lr_all %>% 
  filter(target == "Crypt", pathway_name %in% hypA_paths) %>%
  mutate(hypothesis = "A: Raft Sensitization (Crypt Incoming)")

# Hypothesis B: Lipid-Modified Secretion (OUTGOING from Crypt)
hypB_paths <- c("HH", "WNT")
df_B <- lr_all %>% 
  filter(source == "Crypt", pathway_name %in% hypB_paths) %>%
  mutate(hypothesis = "B: Lipid Secretion (Crypt Outgoing)")

# Hypothesis C: Positional Boundary Breakdown (Ephrin signaling at the boundary)
hypC_paths <- c("EPHA", "EPHB")
df_C <- lr_all %>% 
  filter((source %in% c("Crypt", "Junction") | target %in% c("Crypt", "Junction")), 
         pathway_name %in% hypC_paths) %>%
  mutate(hypothesis = "C: Positional Slippage (Ephrin Breakdown)")

# Combine into one analysis frame
hyp_all <- bind_rows(df_A, df_B, df_C)

# 3. Build the Summary Table
tab <- hyp_all %>% 
  group_by(hypothesis, pathway_name, genotype) %>% 
  summarise(prob = round(sum(prob), 4), .groups = "drop") %>%
  pivot_wider(names_from = genotype, values_from = prob, values_fill = 0)

# Reorder columns to match your MAPK dose order
tab <- tab[, c("hypothesis", "pathway_name", dose_order)]

cat("\n=== CHOLESTEROL MEMBRANE HYPOTHESES: CELLCHAT RESULTS ===\n")
print(as.data.frame(tab), row.names = FALSE)
write.csv(tab, "TABLE_cholesterol_membrane_hypotheses.csv", row.names = FALSE)

# 4. Generate the Visualization
plot_df <- hyp_all %>% 
  group_by(hypothesis, pathway_name, genotype) %>% 
  summarise(prob = sum(prob), .groups = "drop")

plot_df$genotype <- factor(plot_df$genotype, levels = dose_order)

pdf("FIGGIE_cholesterol_membrane_hypotheses.pdf", width = 14, height = 8)
p <- ggplot(plot_df, aes(x = genotype, y = prob, fill = genotype)) +
  geom_col() +
  facet_wrap(~ hypothesis + pathway_name, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = geno_cols) +
  labs(x = NULL, y = "Communication Probability",
       title = "Membrane Mechanics: Raft Sensitization, Secretion, and Positional Slippage",
       subtitle = "Hypothesis A & B should rise with Braf dose. Hypothesis C (Ephrins) should collapse with Braf dose.") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        strip.text = element_text(face = "bold", size = 9),
        legend.position = "none")
print(p)
dev.off()

cat("\nSaved FIGGIE_cholesterol_membrane_hypotheses.pdf\n")






library(dplyr)
library(tidyr)

# Isolate Ephrin signaling strictly involving the Crypt (removing Junction)
hypC_crypt_paths <- c("EPHA", "EPHB")

df_C_crypt <- lr_all %>% 
  filter((source == "Crypt" | target == "Crypt") & 
           !(source == "Junction" | target == "Junction"), 
         pathway_name %in% hypC_crypt_paths)

# Summarize the Crypt-isolated Ephrin data
tab_crypt_ephrin <- df_C_crypt %>% 
  group_by(pathway_name, genotype) %>% 
  summarise(prob = round(sum(prob), 4), .groups = "drop") %>%
  pivot_wider(names_from = genotype, values_from = prob, values_fill = 0)

# Order columns by the MAPK dose
tab_crypt_ephrin <- tab_crypt_ephrin[, c("pathway_name", dose_order)]

cat("\n=== ISOLATED CRYPT EPHRIN SIGNALING ===\n")
print(as.data.frame(tab_crypt_ephrin), row.names = FALSE)




library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Fibroblast -> Crypt (The Stromal Response)
# Looking specifically at EGF (which we know went up) and WNT (which went down), plus BMP
fibro_to_crypt_paths <- c("EGF", "WNT", "ncWNT", "BMP", "FGF")

df_fibro_crypt <- lr_all %>% 
  filter(source == "Fibroblast", target == "Crypt", pathway_name %in% fibro_to_crypt_paths) %>%
  group_by(pathway_name, genotype) %>% 
  summarise(prob = round(sum(prob), 4), .groups = "drop") %>%
  pivot_wider(names_from = genotype, values_from = prob, values_fill = 0)

# Reorder columns
df_fibro_crypt <- df_fibro_crypt[, c("pathway_name", intersect(dose_order, names(df_fibro_crypt)))]

cat("\n=== FIBROBLAST -> CRYPT SIGNALING ===\n")
print(as.data.frame(df_fibro_crypt), row.names = FALSE)


# 2. Crypt -> Fibroblast (The Metabolic Exhaust / Epithelial Distress)
# Looking at remodeling and inflammatory signals
crypt_to_fibro_paths <- c("TGFb", "SPP1", "MIF", "COLLAGEN", "PDGF")

df_crypt_fibro <- lr_all %>% 
  filter(source == "Crypt", target == "Fibroblast", pathway_name %in% crypt_to_fibro_paths) %>%
  group_by(pathway_name, genotype) %>% 
  summarise(prob = round(sum(prob), 4), .groups = "drop") %>%
  pivot_wider(names_from = genotype, values_from = prob, values_fill = 0)

# Reorder columns
df_crypt_fibro <- df_crypt_fibro[, c("pathway_name", intersect(dose_order, names(df_crypt_fibro)))]

cat("\n=== CRYPT -> FIBROBLAST SIGNALING ===\n")
print(as.data.frame(df_crypt_fibro), row.names = FALSE)





library(Seurat)
library(dplyr)
library(tidyr)

# 1. The Autocrine LXR-AREG Loop (The Breakthrough)
# Nr1h3 (LXRalpha), Nr1h2 (LXRbeta), Abca1/Abcg1 (LXR targets/efflux), Areg (EGF ligand)
lxr_areg_genes <- c("Nr1h3", "Nr1h2", "Abca1", "Abcg1", "Areg")

# 2. Ephrin Boundary Tightening (Positional Slippage check)
ephrin_genes <- c("Ephb2", "Ephb3", "Efnb1", "Efnb2")

# 3. Inflammatory Distress 
distress_genes <- c("Mif")

# 4. Stromal Starvation (Fibroblasts turning off life support)
stromal_starvation_genes <- c("Wnt2b", "Bmp4", "Bmp5")

# Assuming obj and dose_order are already in your environment
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")

# ==========================================
# PART A: CRYPT ANALYSIS 
# Testing the LXR-AREG loop, Ephrin tightening, and Distress
# ==========================================
crypt <- subset(obj, celltype == "Crypt")
Idents(crypt) <- "geno"
crypt_genes <- c(lxr_areg_genes, ephrin_genes, distress_genes)
crypt_genes_present <- intersect(crypt_genes, rownames(crypt))

cat("\n=== RAW EXPRESSION IN CRYPT (LXR-AREG Loop, Ephrins, Distress) ===\n")
crypt_expr <- AverageExpression(crypt, features = crypt_genes_present, assays = "Spatial")$Spatial
crypt_df <- round(as.data.frame(as.matrix(crypt_expr))[, intersect(dose_order, colnames(crypt_expr))], 3)
print(crypt_df)

# ==========================================
# PART B: FIBROBLAST ANALYSIS 
# Testing the Stromal Starvation response
# ==========================================
fibro <- subset(obj, celltype == "Fibroblast")
Idents(fibro) <- "geno"
fibro_genes_present <- intersect(stromal_starvation_genes, rownames(fibro))

cat("\n=== RAW EXPRESSION IN FIBROBLAST (Stromal Starvation) ===\n")
fibro_expr <- AverageExpression(fibro, features = fibro_genes_present, assays = "Spatial")$Spatial
fibro_df <- round(as.data.frame(as.matrix(fibro_expr))[, intersect(dose_order, colnames(fibro_expr))], 3)
print(fibro_df)


# Co-expression mapping for Autocrine Loop
areg_egfr <- FetchData(crypt, vars = c("Areg", "Egfr"))
cor_val <- cor(areg_egfr$Areg, areg_egfr$Egfr, method = "spearman")
cat("\n=== Spearman correlation of Areg and Egfr in Crypts: ", round(cor_val, 3), " ===\n")

# Spatial plot of the autocrine potential
obj$Areg_Egfr_Score <- (obj@assays$Spatial@data["Areg",] > 0 & obj@assays$Spatial@data["Egfr",] > 0)
pdf("SPATIAL_Autocrine_Loop.pdf", width=8, height=6)
SpatialPlot(obj, features = "Areg_Egfr_Score", pt.size.factor = 1.5) + ggtitle("Autocrine Loop Potential")
dev.off()



# Correct way to pull expression for Assay5 objects
areg_vals <- GetAssayData(obj, assay = "Spatial", layer = "data")["Areg", ]
egfr_vals <- GetAssayData(obj, assay = "Spatial", layer = "data")["Egfr", ]

# Create the score
obj$Areg_Egfr_Score <- (areg_vals > 0 & egfr_vals > 0)

# Spatial Plot
library(Seurat)
SpatialPlot(obj, features = "Areg_Egfr_Score", pt.size.factor = 1.5) + 
  ggtitle("Spatial Co-localization of Areg and Egfr")

# 1. Define sets of barcodes where genes are detected
areg_pos <- colnames(obj)[GetAssayData(obj, assay="Spatial", layer="counts")["Areg",] > 0]
egfr_pos <- colnames(obj)[GetAssayData(obj, assay="Spatial", layer="counts")["Egfr",] > 0]

# 2. Calculate distance matrix from your projection coordinates
coords <- obj@meta.data[, c("x", "y")] # Ensure these are the projection coordinates
dist_mat <- as.matrix(dist(coords))

# 3. Test if Areg spots are closer to Egfr spots than they are to random spots
# Extract distances between all Areg+ and Egfr+ pairs
relevant_dists <- dist_mat[areg_pos, egfr_pos]

# Compare to the distance of all spots in the tissue
null_dists <- dist_mat[sample(colnames(obj), 500), sample(colnames(obj), 500)]

cat("\n=== Spatial Proximity Analysis (Areg vs Egfr) ===\n")
cat("Mean distance between Areg and Egfr spots:", mean(relevant_dists), "\n")
cat("Mean distance between random spots:", mean(null_dists), "\n")
library(ggplot2)

plot_df <- obj@meta.data
plot_df$Areg <- GetAssayData(obj, assay="Spatial", layer="counts")["Areg",] > 0
plot_df$Egfr <- GetAssayData(obj, assay="Spatial", layer="counts")["Egfr",] > 0

# Create a combined status column
plot_df$Status <- "Neither"
plot_df$Status[plot_df$Areg] <- "Areg+"
plot_df$Status[plot_df$Egfr] <- "Egfr+"
plot_df$Status[plot_df$Areg & plot_df$Egfr] <- "Both"

pdf("SPATIAL_Areg_Egfr_Map.pdf", width=8, height=7)
ggplot(plot_df, aes(x=x, y=y, color=Status)) +
  geom_point(size=0.8, alpha=0.7) +
  scale_color_manual(values=c("Neither"="grey90", "Areg+"="blue", "Egfr+"="orange", "Both"="purple")) +
  theme_void() + 
  ggtitle("Spatial Distribution: Areg and Egfr")
dev.off()



# 1. Add coordinates to meta.data if they aren't there
proj <- read.csv("Default Pipeline Analysis-Spatial-Projection.csv")
# Ensure the barcodes match
rownames(proj) <- proj$Barcode
obj@meta.data$x <- proj[colnames(obj), "X.Coordinate"]
obj@meta.data$y <- proj[colnames(obj), "Y.Coordinate"]

# 2. Re-run Proximity Analysis
areg_pos <- colnames(obj)[GetAssayData(obj, assay="Spatial", layer="counts")["Areg",] > 0]
egfr_pos <- colnames(obj)[GetAssayData(obj, assay="Spatial", layer="counts")["Egfr",] > 0]

# Subset to valid spots only
coords <- obj@meta.data[, c("x", "y")]
dist_mat <- as.matrix(dist(coords))

relevant_dists <- dist_mat[areg_pos, egfr_pos]
mean_dist <- mean(relevant_dists)


library(Seurat); library(dplyr); library(tidyr)
obj <- readRDS("obj_annotated.rds")
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")

# ============ Q1: is it LXR-driven? SREBP-2 vs LXR are DIFFERENT programs ============
# SREBP-2 = cholesterol SYNTHESIS (what we found). LXR = cholesterol EFFLUX/oxysterol sensing.
# These are opposing arms — important to distinguish!
programs <- list(
  SREBP2_synthesis = c("Srebf2","Hmgcr","Hmgcs1","Sqle","Fdft1","Cyp51","Dhcr24","Ldlr","Insig1"),
  LXR_program      = c("Nr1h3","Nr1h2","Abca1","Abcg1","Abcg5","Abcg8","Srebf1","Fasn","Scd1","Apoe"),
  Oxysterol_enzymes= c("Cyp27a1","Cyp7a1","Cyp7b1","Ch25h","Cyp46a1"),   # make LXR ligands
  AREG_EGF         = c("Areg","Ereg","Hbegf","Egfr"),
  Ephrin           = c("Efna1","Efnb1","Efnb2","Epha2","Ephb2","Ephb3","Ephb4")
)
programs <- lapply(programs, function(g) g[g %in% rownames(obj)])
for(n in names(programs)) cat(n, "found:", paste(programs[[n]], collapse=", "), "\n")

obj <- AddModuleScore(obj, features=programs, name=names(programs))

# ============ Q2 + Q3: WHERE (which zone) and WHICH CELLS drive each program ============
score_cols <- paste0(names(programs), seq_along(programs))
# by cell type
cat("\n===== module scores BY CELL TYPE =====\n")
md <- FetchData(obj, vars=c(score_cols, "celltype"))
by_ct <- md %>% group_by(celltype) %>% summarise(across(all_of(score_cols), ~round(mean(.x),3)),.groups="drop")
colnames(by_ct) <- c("celltype", names(programs))
print(as.data.frame(by_ct), row.names=FALSE)
    











library(Seurat); library(dplyr)
# proliferation module
prolif <- list(prolif=intersect(c("Mki67","Pcna","Top2a","Ccnb1","Cdk1","Birc5","Mcm2","Mcm6"), rownames(obj)))
obj <- AddModuleScore(obj, features=prolif, name="Prolif")

# in crypt: does cholesterol synthesis track proliferation? per genotype
crypt <- subset(obj, celltype=="Crypt")
cd <- FetchData(crypt, vars=c("SREBP2_synthesis1","Prolif1","geno"))
cat("=== crypt: SREBP2-cholesterol vs proliferation correlation, per genotype ===\n")
for(g in dose_order){
  sub <- cd[cd$geno==g,]
  if(nrow(sub)>30) cat(sprintf("%-18s r = %+.3f (n=%d)\n", g,
                               cor(sub$SREBP2_synthesis1, sub$Prolif1, method="spearman"), nrow(sub)))
}

# and does proliferation itself go up with Braf dose in crypt?
cat("\n=== crypt proliferation score by genotype ===\n")
print(cd %>% group_by(geno) %>% summarise(prolif=round(mean(Prolif1),3), chol=round(mean(SREBP2_synthesis1),3),.groups="drop") %>%
        arrange(match(geno,dose_order)) %>% as.data.frame(), row.names=FALSE)







library(Seurat); library(dplyr); library(tidyr)
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")

# 1. Do oxysterol-PRODUCING enzymes rise with Braf in epithelium? (the ligand source)
oxy <- intersect(c("Cyp27a1","Ch25h","Cyp7b1","Cyp46a1","Hsd11b1"), rownames(obj))
epi <- subset(obj, celltype %in% c("Crypt","Junction","Villus","Villus_tip"))
eo <- FetchData(epi, vars=c(oxy,"geno"))
cat("=== oxysterol-producing enzymes in EPITHELIUM x genotype ===\n")
t <- eo %>% group_by(geno) %>% summarise(across(all_of(oxy), ~round(mean(.x),3)),.groups="drop")
print(as.data.frame(t(t[,-1]) %>% `colnames<-`(t$geno))[,intersect(dose_order,t$geno)])

# 2. Immune cell COMPOSITION per genotype (is there T-cell depletion with Braf?)
cat("\n=== immune cell counts + fractions per genotype ===\n")
comp <- FetchData(obj, vars=c("celltype","geno")) %>%
  filter(celltype %in% c("Tcell","Macrophage","DC","Plasma","Myeloid")) %>%
  group_by(geno, celltype) %>% summarise(n=n(),.groups="drop") %>%
  group_by(geno) %>% mutate(frac=round(n/sum(n),3)) %>% ungroup()
print(comp %>% select(-n) %>% pivot_wider(names_from=celltype, values_from=frac) %>%
        arrange(match(geno,dose_order)) %>% as.data.frame(), row.names=FALSE)

# 3. LXR TARGET activity in macrophages (the abundant immune cell) across genotypes
mac <- subset(obj, celltype=="Macrophage")
lxr_t <- intersect(c("Abca1","Abcg1","Apoe","Srebf1"), rownames(obj))
em <- FetchData(mac, vars=c(lxr_t,"geno"))
cat("\n=== LXR target genes in MACROPHAGES x genotype ===\n")
tm <- em %>% group_by(geno) %>% summarise(across(all_of(lxr_t), ~round(mean(.x),3)),.groups="drop")
print(as.data.frame(t(tm[,-1]) %>% `colnames<-`(tm$geno))[,intersect(dose_order,tm$geno)])




library(Seurat); library(dplyr); library(tidyr)
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")

# apoptosis / cell-death modules (both pro-death and survival)
death_mods <- list(
  Apoptosis   = c("Bax","Bak1","Bcl2l11","Bbc3","Pmaip1","Casp3","Casp8","Casp9","Cycs","Apaf1","Bid"),
  Survival    = c("Bcl2","Bcl2l1","Mcl1","Bcl2a1a","Birc5","Xiap","Cflar","Bcl2l2"),
  Anoikis     = c("Bcl2l11","Bmf","Bak1","Bax"),   # detachment death (crypt cells at base)
  ISR_stress  = c("Ddit3","Atf4","Atf3","Trib3","Chac1")  # integrated stress response
)
death_mods <- lapply(death_mods, function(g) g[g %in% rownames(obj)])
for(n in names(death_mods)) cat(n, ":", paste(death_mods[[n]], collapse=", "), "\n")
obj <- AddModuleScore(obj, features=death_mods, name=names(death_mods))
dcols <- paste0(names(death_mods), seq_along(death_mods))

# ============ TEST 1: apoptosis/survival in CRYPT across genotypes ============
crypt <- subset(obj, celltype=="Crypt")
cd <- FetchData(crypt, vars=c("SREBP2_synthesis1", dcols, "geno"))
colnames(cd) <- c("Cholesterol", names(death_mods), "geno")
cat("\n=== CRYPT: cholesterol + death/survival modules by genotype ===\n")
print(cd %>% group_by(geno) %>%
        summarise(Cholesterol=round(mean(Cholesterol),3),
                  Apoptosis=round(mean(Apoptosis),3),
                  Survival=round(mean(Survival),3),
                  Anoikis=round(mean(Anoikis),3),
                  ISR_stress=round(mean(ISR_stress),3),.groups="drop") %>%
        arrange(match(geno,dose_order)) %>% as.data.frame(), row.names=FALSE)

# ============ TEST 2: single-cell correlation cholesterol vs apoptosis (the shield test) ============
cat("\n=== CRYPT: correlation cholesterol vs apoptosis, per genotype ===\n")
cat("(shield hypothesis predicts NEGATIVE: high cholesterol = low apoptosis)\n")
for(g in dose_order){
  s <- cd[cd$geno==g,]
  if(nrow(s)>30) cat(sprintf("%-18s  chol-vs-Apoptosis r = %+.3f | chol-vs-Survival r = %+.3f (n=%d)\n",
                             g, cor(s$Cholesterol, s$Apoptosis, method="spearman"),
                             cor(s$Cholesterol, s$Survival, method="spearman"), nrow(s)))
}

# ============ TEST 3: split crypt cells into cholesterol-high vs -low, compare apoptosis ============
cat("\n=== apoptosis in cholesterol-HIGH vs -LOW crypt cells (within Braf_homo) ===\n")
bh <- cd[cd$geno=="Braf_homo",]
bh$chol_group <- ifelse(bh$Cholesterol > median(bh$Cholesterol), "chol_HIGH", "chol_LOW")
print(bh %>% group_by(chol_group) %>%
        summarise(Apoptosis=round(mean(Apoptosis),3), Survival=round(mean(Survival),3),
                  Anoikis=round(mean(Anoikis),3), n=n(),.groups="drop") %>% as.data.frame(), row.names=FALSE)






library(Seurat); library(dplyr); library(tidyr)
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")

# ============ modules: BMP, Wnt, Myc, E2F, Fetal ============
mods <- list(
  BMP_activity = c("Id1","Id2","Id3","Smad6","Smad7","Bambi","Smad4"),
  Wnt_targets  = c("Axin2","Lgr5","Ascl2","Ccnd1","Sox9","Ephb2","Rnf43","Tcf7"),
  Myc_targets  = c("Myc","Npm1","Ncl","Nop56","Ppan","Srm","Odc1"),
  E2F_targets  = c("E2f1","Mki67","Pcna","Mcm2","Mcm6","Ccne1","Cdc6","Rrm2"),
  Fetal_revert = c("Ly6a","Anxa1","Anxa3","Clu","Tacstd2","Areg","Il33","Cldn4","Spp1")  # fetal/regen markers
)
mods <- lapply(mods, function(g) g[g %in% rownames(obj)])
for(n in names(mods)) cat(n,":",paste(mods[[n]],collapse=", "),"\n")
obj <- AddModuleScore(obj, features=mods, name=names(mods))
mcols <- paste0(names(mods), seq_along(mods))

# ============ TEST A: all modules in CRYPT across genotypes (co-enrichment with cholesterol) ============
crypt <- subset(obj, celltype=="Crypt")
cd <- FetchData(crypt, vars=c("SREBP2_synthesis1", mcols, "geno"))
colnames(cd) <- c("Cholesterol", names(mods), "geno")
cat("\n=== CRYPT: cholesterol + parallel pathways by genotype ===\n")
print(cd %>% group_by(geno) %>% summarise(across(c("Cholesterol",names(mods)), ~round(mean(.x),3)),.groups="drop") %>%
        arrange(match(geno,dose_order)) %>% as.data.frame(), row.names=FALSE)

# ============ TEST B: single-cell correlation of cholesterol with each pathway (crypt, Braf_homo) ============
cat("\n=== correlation of cholesterol with each pathway (Braf_homo crypt cells) ===\n")
bh <- cd[cd$geno=="Braf_homo",]
for(m in names(mods)){
  r <- cor(bh$Cholesterol, bh[[m]], method="spearman")
  cat(sprintf("  Cholesterol vs %-14s r = %+.3f\n", m, r))
}

# ============ TEST C: the 6 master TFs from the paper — are they Braf-up in crypt? ============
tfs <- intersect(c("Cebpb","Creb3","Nr2f1","Klf16","Sp6","Fosl1"), rownames(obj))
cat("\n=== the 6 cholesterol master-TFs (Rzasa) in CRYPT x genotype ===\n")
et <- FetchData(crypt, vars=c(tfs,"geno"))
tt <- et %>% group_by(geno) %>% summarise(across(all_of(tfs), ~round(mean(.x),3)),.groups="drop")
print(as.data.frame(t(tt[,-1]) %>% `colnames<-`(tt$geno))[,intersect(dose_order,tt$geno)])

# ============ TEST D: SMAD4/BMP-cholesterol vulnerability link ============
# do cells with HIGH BMP activity also depend on (co-express) cholesterol?
cat("\n=== cholesterol in BMP-high vs BMP-low crypt cells (Braf_homo) ===\n")
bh$bmp_group <- ifelse(bh$BMP_activity > median(bh$BMP_activity), "BMP_high","BMP_low")
print(bh %>% group_by(bmp_group) %>% summarise(Cholesterol=round(mean(Cholesterol),3), n=n(),.groups="drop") %>%
        as.data.frame(), row.names=FALSE)






library(Seurat); library(CellChat); library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
obj <- readRDS("obj_annotated.rds")
cc_list <- readRDS("cc_plasma_all.rds")
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
names(cc_list) <- short[names(cc_list)]
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")
geno_cols <- c("Braf_homo"="#B2182B","Braf_het"="#F4A582","Control"="grey60",
               "ErkDKO"="#92C5DE","Braf_homo_ErkDKO"="#2166AC","Braf_het_ErkDKO"="#4393C3")
obj$geno <- factor(obj$geno, levels=dose_order)

# ===== GRAPH 1: master summary — all key modules x genotype, one heatmap =====
# (need the modules computed; recompute the key ones)
mods <- list(
  Cholesterol = c("Srebf2","Hmgcr","Hmgcs1","Sqle","Fdft1","Cyp51","Dhcr24","Ldlr","Insig1"),
  Fetal       = c("Ly6a","Anxa1","Anxa3","Clu","Tacstd2","Areg","Il33","Cldn4","Spp1"),
  EGF_ligands = c("Areg","Hbegf","Tgfa","Ereg"),
  STAT3_axis  = c("Stat3","Socs3","Cflar","Il22ra1"),
  Guanylin    = c("Guca2a","Guca2b","Gucy2c")
)
mods <- lapply(mods, function(g) g[g %in% rownames(obj)])
obj <- AddModuleScore(obj, features=mods, name=names(mods))
mcols <- paste0(names(mods), seq_along(mods))

# scale each module 0-1 across genotypes for the heatmap
summ <- FetchData(obj, vars=c(mcols,"geno")) %>%
  group_by(geno) %>% summarise(across(all_of(mcols), mean),.groups="drop")
colnames(summ) <- c("geno", names(mods))
sl <- summ %>% pivot_longer(-geno, names_to="module", values_to="score") %>%
  group_by(module) %>% mutate(scaled=(score-min(score))/(max(score)-min(score))) %>% ungroup()
sl$geno <- factor(sl$geno, levels=dose_order)

pdf("SUMMARY_1_module_heatmap.pdf", width=7, height=4)
print(ggplot(sl, aes(geno, module, fill=scaled)) +
        geom_tile(color="white") +
        geom_text(aes(label=round(score,2)), size=3) +
        scale_fill_gradient(low="#F7FBFF", high="#B2182B", name="scaled") +
        labs(x=NULL, y=NULL, title="Key programs across the allelic series (MAPK-dose order)") +
        theme_minimal() + theme(axis.text.x=element_text(angle=45,hjust=1)))
dev.off()

# ===== GRAPH 2: FOSL1 — the driver — clean dose-response =====
crypt <- subset(obj, celltype=="Crypt")
fd <- FetchData(crypt, vars=c("Fosl1","geno"))
pdf("SUMMARY_2_FOSL1.pdf", width=6, height=4)
print(ggplot(fd %>% group_by(geno) %>% summarise(Fosl1=mean(Fosl1),.groups="drop"),
             aes(geno, Fosl1, fill=geno)) + geom_col() +
        scale_fill_manual(values=geno_cols) +
        labs(x=NULL, y="Fosl1 expression (crypt)",
             title="FOSL1 (AP-1/ERK cholesterol driver): strictly Braf-dependent, ERK-reversed") +
        theme_classic() + theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="none"))
dev.off()

# ===== GRAPH 3: which CELL TYPES are the signaling hubs? outgoing strength per cell type x genotype =====
# aggregate CellChat outgoing strength by cell type
hub <- bind_rows(lapply(dose_order, function(g){
  cc <- cc_list[[g]]
  s <- rowSums(cc@net$weight)   # total outgoing per cell type
  data.frame(celltype=names(s), outgoing=s, genotype=g)
}))
hub$genotype <- factor(hub$genotype, levels=dose_order)
pdf("SUMMARY_3_signaling_hubs.pdf", width=10, height=5)
print(ggplot(hub, aes(genotype, outgoing, fill=celltype)) +
        geom_col(position="stack") +
        labs(x=NULL, y="Total outgoing signaling", fill="Cell type",
             title="Which cell populations drive signaling across genotypes") +
        theme_classic() + theme(axis.text.x=element_text(angle=45,hjust=1)))
dev.off()

# ===== GRAPH 4: T-cell depletion (absolute density) =====
comp <- FetchData(obj, vars=c("celltype","geno"))
dens <- comp %>% group_by(geno) %>%
  summarise(Tcell_per_1k = round(1000*sum(celltype=="Tcell")/n(),1),
            Macro_per_1k = round(1000*sum(celltype=="Macrophage")/n(),1),.groups="drop")
dens$geno <- factor(dens$geno, levels=dose_order)
pdf("SUMMARY_4_Tcell_depletion.pdf", width=6, height=4)
print(ggplot(dens, aes(geno, Tcell_per_1k, fill=geno)) + geom_col() +
        scale_fill_manual(values=geno_cols) +
        labs(x=NULL, y="T cells per 1000 cells",
             title="T-cell depletion in Braf_homo, restored by Erk deletion") +
        theme_classic() + theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="none"))
dev.off()

cat("4 summary graphs saved\n")
print(as.data.frame(dens), row.names=FALSE)  # show T-cell density numbers






library(Seurat); library(dplyr); library(tidyr); library(ggplot2)
obj <- readRDS("obj_annotated.rds")
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")
geno_cols <- c("Braf_homo"="#B2182B","Braf_het"="#F4A582","Control"="grey60",
               "ErkDKO"="#92C5DE","Braf_homo_ErkDKO"="#2166AC","Braf_het_ErkDKO"="#4393C3")
obj$geno <- factor(obj$geno, levels=dose_order)

# cholesterol module
chol <- list(Cholesterol=intersect(c("Srebf2","Hmgcr","Hmgcs1","Mvk","Mvd","Pmvk","Fdps","Fdft1",
                                     "Sqle","Lss","Cyp51","Dhcr7","Dhcr24","Idi1","Insig1","Ldlr"), rownames(obj)))
obj <- AddModuleScore(obj, features=chol, name="Chol")

# ===== GRAPH 1: cholesterol by CELL TYPE x genotype (WHICH cells carry it) =====
d1 <- FetchData(obj, vars=c("Chol1","celltype","geno"))
by_ct <- d1 %>% group_by(celltype, geno) %>% summarise(m=mean(Chol1),.groups="drop")
by_ct$geno <- factor(by_ct$geno, levels=dose_order)
# order cell types by overall cholesterol
ct_ord <- by_ct %>% group_by(celltype) %>% summarise(s=mean(m),.groups="drop") %>% arrange(desc(s)) %>% pull(celltype)
by_ct$celltype <- factor(by_ct$celltype, levels=ct_ord)
pdf("CHOL_1_by_celltype.pdf", width=9, height=5)
print(ggplot(by_ct, aes(celltype, m, fill=geno)) +
        geom_col(position="dodge") +
        scale_fill_manual(values=geno_cols) +
        labs(x=NULL, y="Cholesterol module score", fill="Genotype",
             title="Cholesterol biosynthesis by cell type across genotypes") +
        theme_classic() + theme(axis.text.x=element_text(angle=45,hjust=1)))
dev.off()

# ===== GRAPH 2: cholesterol by ZONE (epithelial) x genotype — the spatial gradient =====
zones <- c("Crypt","Junction","Villus","Villus_tip")
d2 <- d1 %>% filter(celltype %in% zones)
d2$celltype <- factor(d2$celltype, levels=zones)
d2$geno <- factor(d2$geno, levels=dose_order)
z_summ <- d2 %>% group_by(celltype, geno) %>% summarise(m=mean(Chol1),.groups="drop")
pdf("CHOL_2_zone_gradient.pdf", width=9, height=5)
print(ggplot(z_summ, aes(celltype, m, group=geno, color=geno)) +
        geom_line(linewidth=1) + geom_point(size=2.5) +
        scale_color_manual(values=geno_cols) +
        labs(x="Epithelial zone (crypt→tip)", y="Cholesterol module score", color="Genotype",
             title="Cholesterol biosynthesis along the crypt-villus axis") +
        theme_classic())
dev.off()

# ===== GRAPH 3: the pathway genes as a heatmap, crypt only, x genotype =====
key <- chol$Cholesterol
crypt <- subset(obj, celltype=="Crypt")
ec <- FetchData(crypt, vars=c(key,"geno"))
gm <- ec %>% group_by(geno) %>% summarise(across(all_of(key), mean),.groups="drop") %>%
  pivot_longer(-geno, names_to="gene", values_to="expr") %>%
  group_by(gene) %>% mutate(scaled=(expr-min(expr))/(max(expr)-min(expr)+1e-9)) %>% ungroup()
gm$geno <- factor(gm$geno, levels=dose_order)
# order genes by pathway position (roughly)
gene_ord <- c("Srebf2","Insig1","Hmgcs1","Hmgcr","Mvk","Pmvk","Mvd","Idi1","Fdps","Fdft1","Sqle","Lss","Cyp51","Dhcr24","Dhcr7","Ldlr")
gm$gene <- factor(gm$gene, levels=rev(intersect(gene_ord, unique(gm$gene))))
pdf("CHOL_3_gene_heatmap.pdf", width=7, height=6)
print(ggplot(gm, aes(geno, gene, fill=scaled)) +
        geom_tile(color="white") +
        scale_fill_gradient(low="#FFF7EC", high="#7F0000", name="scaled\nexpr") +
        labs(x=NULL, y=NULL, title="Mevalonate/cholesterol pathway genes in crypt (dose order)") +
        theme_minimal() + theme(axis.text.x=element_text(angle=45,hjust=1)))
dev.off()

cat("=== cholesterol by cell type (mean across genotypes) ===\n")
print(by_ct %>% group_by(celltype) %>% summarise(mean=round(mean(m),3),.groups="drop") %>%
        arrange(desc(mean)) %>% as.data.frame(), row.names=FALSE)







library(CellChat)
cc_list <- readRDS("cc_plasma_all.rds")
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
names(cc_list) <- short[names(cc_list)]
cc_list <- cc_list[dose_order]
lr_all <- bind_rows(Map(function(cc,g){d<-subsetCommunication(cc);d$genotype<-g;d}, cc_list, names(cc_list)))

# what signals are RECEIVED by the Crypt (where cholesterol is high)?
# does any incoming pathway to crypt track the cholesterol dose-response?
crypt_in <- lr_all %>% filter(target=="Crypt") %>%
  group_by(pathway_name, genotype) %>% summarise(prob=sum(prob),.groups="drop") %>%
  pivot_wider(names_from=genotype, values_from=prob, values_fill=0)
crypt_in <- crypt_in[, c("pathway_name", dose_order)]
crypt_in$braf_effect <- crypt_in$Braf_homo - crypt_in$Control
crypt_in <- crypt_in[order(-crypt_in$braf_effect),]
cat("\n=== pathways INTO crypt, ranked by Braf-induction (potential cholesterol drivers) ===\n")
print(as.data.frame(head(crypt_in, 12)), row.names=FALSE)




library(Seurat); library(CellChat); library(dplyr); library(tidyr)
# obj, cc_list, lr_all, dose_order loaded
obj <- readRDS("obj_annotated.rds"); dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")
chol <- list(Chol=intersect(c("Srebf2","Hmgcr","Hmgcs1","Sqle","Fdft1","Cyp51","Dhcr24","Ldlr","Insig1"),rownames(obj)))
obj <- AddModuleScore(obj, features=chol, name="Chol")

# ===== TEST 1: which stromal/immune signals go INTO crypt from FIBRO/MAC specifically, tracking cholesterol =====
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
cc_list <- readRDS("cc_plasma_all.rds"); names(cc_list) <- short[names(cc_list)]; cc_list <- cc_list[dose_order]
lr_all <- bind_rows(Map(function(cc,g){d<-subsetCommunication(cc);d$genotype<-g;d}, cc_list, names(cc_list)))

# ONLY fibroblast/macrophage -> crypt, ranked by tracking cholesterol dose-response
fm_crypt <- lr_all %>% filter(source %in% c("Fibroblast","Macrophage"), target=="Crypt") %>%
  group_by(pathway_name, source, genotype) %>% summarise(prob=sum(prob),.groups="drop") %>%
  pivot_wider(names_from=genotype, values_from=prob, values_fill=0)
for(g in dose_order) if(!g %in% names(fm_crypt)) fm_crypt[[g]] <- 0
fm_crypt$braf_up   <- fm_crypt$Braf_homo - fm_crypt$Control
fm_crypt$erk_rev   <- fm_crypt$Braf_homo - fm_crypt$Braf_homo_ErkDKO
# a real catalyst candidate: braf_up>0 AND erk_rev>0 (matches cholesterol pattern)
cand <- fm_crypt %>% filter(braf_up>0, erk_rev>0) %>% arrange(desc(braf_up)) %>%
  select(pathway_name, source, Control, Braf_homo, Braf_homo_ErkDKO, braf_up, erk_rev)
cat("=== FIBRO/MAC -> crypt signals that MATCH cholesterol pattern (Braf-up, Erk-rev) ===\n")
print(as.data.frame(head(cand,15)), row.names=FALSE)

# ===== TEST 2: inflammatory signals specifically (the "inflammatory catalyst" idea) =====
inflam <- c("TNF","IL1","IL6","IL17","IFN-II","CXCL","CCL","COMPLEMENT","OSM","TGFb")
inf_crypt <- lr_all %>% filter(pathway_name %in% inflam, target=="Crypt") %>%
  group_by(pathway_name, genotype) %>% summarise(prob=round(sum(prob),4),.groups="drop") %>%
  pivot_wider(names_from=genotype, values_from=prob, values_fill=0)
inf_crypt <- inf_crypt[, c("pathway_name", intersect(dose_order,names(inf_crypt)))]
cat("\n=== inflammatory signals INTO crypt x genotype ===\n")
print(as.data.frame(inf_crypt), row.names=FALSE)

# ===== TEST 3: SPATIAL — are FOSL1-high/cholesterol-high crypt cells physically near fibroblasts/macrophages? =====
# proximity test: for each crypt cell, is cholesterol higher when near a fibroblast/macrophage?
proj <- read.csv("Default Pipeline Analysis-Spatial-Projection.csv")
md <- FetchData(obj, vars=c("Chol1","Fosl1","celltype","geno"))
md$Barcode <- rownames(md)
md$x <- proj$X.Coordinate[match(md$Barcode, proj$Barcode)]
md$y <- proj$Y.Coordinate[match(md$Barcode, proj$Barcode)]
md <- md[!is.na(md$x),]

# for Braf_homo: compute each crypt cell's distance to nearest fibroblast/macrophage
library(FNN)
bh <- md[md$geno=="Braf_homo",]
crypt_bh <- bh[bh$celltype=="Crypt",]
stroma_bh <- bh[bh$celltype %in% c("Fibroblast","Macrophage"),]
if(nrow(stroma_bh)>5 & nrow(crypt_bh)>5){
  nn <- get.knnx(stroma_bh[,c("x","y")], crypt_bh[,c("x","y")], k=1)
  crypt_bh$dist_to_stroma <- nn$nn.dist[,1]
  # split crypt cells: near vs far from stroma, compare cholesterol
  crypt_bh$prox <- ifelse(crypt_bh$dist_to_stroma < median(crypt_bh$dist_to_stroma), "NEAR_stroma","FAR_stroma")
  cat("\n=== Braf_homo crypt: cholesterol NEAR vs FAR from fibroblast/macrophage ===\n")
  print(crypt_bh %>% group_by(prox) %>%
          summarise(Cholesterol=round(mean(Chol1),3), Fosl1=round(mean(Fosl1),3), n=n(),.groups="drop") %>%
          as.data.frame(), row.names=FALSE)
  cat("\ncorrelation: distance-to-stroma vs cholesterol r =",
      round(cor(crypt_bh$dist_to_stroma, crypt_bh$Chol1, method="spearman"),3), "\n")
  cat("(NEGATIVE r = closer to stroma → higher cholesterol = supports stromal catalyst)\n")
}











library(Seurat); library(dplyr); library(tidyr)
obj <- readRDS("obj_annotated.rds")
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")

# the two opposing arms
mods <- list(
  SREBP2_synth = c("Srebf2","Hmgcr","Hmgcs1","Mvk","Fdps","Fdft1","Sqle","Cyp51","Dhcr24","Ldlr","Insig1"),
  LXR_efflux   = c("Nr1h3","Nr1h2","Abca1","Abcg1","Abcg5","Abcg8","Apoe","Srebf1")
)
mods <- lapply(mods, function(g) g[g %in% rownames(obj)])
obj <- AddModuleScore(obj, features=mods, name=names(mods))

# ===== TEST 1: correlation between the two arms, per genotype, in CRYPT =====
crypt <- subset(obj, celltype=="Crypt")
cd <- FetchData(crypt, vars=c("SREBP2_synth1","LXR_efflux2","geno"))
cat("=== SREBP2 vs LXR correlation in CRYPT, per genotype ===\n")
cat("(inverse seesaw predicts NEGATIVE r)\n")
for(g in dose_order){
  s <- cd[cd$geno==g,]
  if(nrow(s)>30) cat(sprintf("%-18s r = %+.3f (n=%d)\n", g,
                             cor(s$SREBP2_synth1, s$LXR_efflux2, method="spearman"), nrow(s)))
}

# ===== TEST 2: do the two arms move oppositely ACROSS genotypes? (crypt means) =====
cat("\n=== both arms across genotypes (crypt) — do they diverge? ===\n")
print(cd %>% group_by(geno) %>%
        summarise(SREBP2_synth=round(mean(SREBP2_synth1),3),
                  LXR_efflux=round(mean(LXR_efflux2),3),.groups="drop") %>%
        arrange(match(geno,dose_order)) %>% as.data.frame(), row.names=FALSE)

# ===== TEST 3: same correlation but across ALL epithelial cells (more power) =====
epi <- subset(obj, celltype %in% c("Crypt","Junction","Villus","Villus_tip"))
ed <- FetchData(epi, vars=c("SREBP2_synth1","LXR_efflux2","geno"))
cat("\n=== SREBP2 vs LXR correlation, ALL epithelium, per genotype ===\n")
for(g in dose_order){
  s <- ed[ed$geno==g,]
  if(nrow(s)>30) cat(sprintf("%-18s r = %+.3f (n=%d)\n", g,
                             cor(s$SREBP2_synth1, s$LXR_efflux2, method="spearman"), nrow(s)))
}

# ===== TEST 4: split cells by SREBP2 level, look at LXR (the clean seesaw test) =====
cat("\n=== LXR in SREBP2-HIGH vs SREBP2-LOW epithelial cells (Braf_homo) ===\n")
bh <- ed[ed$geno=="Braf_homo",]
bh$grp <- ifelse(bh$SREBP2_synth1 > median(bh$SREBP2_synth1), "SREBP2_high","SREBP2_low")
print(bh %>% group_by(grp) %>%
        summarise(LXR_efflux=round(mean(LXR_efflux2),3), n=n(),.groups="drop") %>%
        as.data.frame(), row.names=FALSE)






library(dplyr)
library(tidyr)

# The extracellular pathways controlled by LXR
lxr_echo_paths <- c("EGF", "APOE", "TNF", "CXCL")

df_lxr_echo <- lr_all %>% 
  filter(source == "Crypt", pathway_name %in% lxr_echo_paths) %>%
  group_by(pathway_name, ligand, genotype) %>% 
  summarise(prob = round(sum(prob), 4), .groups = "drop") %>%
  pivot_wider(names_from = genotype, values_from = prob, values_fill = 0)

# Reorder to match MAPK dose
df_lxr_echo <- df_lxr_echo[, c("pathway_name", "ligand", dose_order)]

cat("\n=== THE LXR EXTRACELLULAR ECHO (Crypt Outgoing) ===\n")
print(as.data.frame(df_lxr_echo), row.names = FALSE)



library(dplyr)
library(tidyr)

macrophages <- subset(obj, celltype == "Macrophage")
mac_ligands <- intersect(c("Nrg1", "Ereg", "Hbegf", "Areg", "Egf", "Tgfa", "Btc"), rownames(macrophages))

# Fetch raw data and manually calculate the means to bypass the Seurat bug
mac_expr_raw <- FetchData(macrophages, vars = c(mac_ligands, "geno"))

mac_tab <- mac_expr_raw %>% 
  group_by(geno) %>% 
  summarise(across(all_of(mac_ligands), ~round(mean(.x), 3)), .groups = "drop")

mac_tab_t <- as.data.frame(t(mac_tab[,-1]))
colnames(mac_tab_t) <- mac_tab$geno

# Reorder to match your standard dose
mac_tab_t <- mac_tab_t[, intersect(dose_order, colnames(mac_tab_t)), drop = FALSE]

cat("\n=== TRUE MACROPHAGE LIGAND EXPRESSION (All Genotypes) ===\n")
print(mac_tab_t)










# =================================================================
# SPATIAL COUPLING (KNN Approach: Scale-Invariant)
# =================================================================

# We will look at the 5 closest macrophages to each crypt spot
k_macs <- 5
crypt_df$Local_Mac_EGF_Load <- 0

cat("Calculating KNN spatial proximity...\n")

for(i in 1:nrow(crypt_df)) {
  # Calculate distances to all macrophages
  dists <- sqrt((mac_df$x - crypt_df$x[i])^2 + (mac_df$y - crypt_df$y[i])^2)
  
  # Find the indices of the 'k' closest macrophages
  closest_idx <- order(dists)[1:k_macs]
  
  # Sum the EGF output of those specific closest macrophages
  crypt_df$Local_Mac_EGF_Load[i] <- sum(mac_df$Mac_EGF_Output[closest_idx])
}

# The Moment of Truth: Correlation
cor_result <- cor(crypt_df$Local_Mac_EGF_Load, crypt_df$SREBP2_synth1, method = "spearman")

cat("\n=== SPATIAL COUPLING RESULTS (Braf_homo) ===\n")
cat("Spearman Correlation (Local Macrophage EGF vs Crypt SREBP-2):", round(cor_result, 3), "\n")

# Split crypts by Median Macrophage EGF Load
med_load <- median(crypt_df$Local_Mac_EGF_Load)
crypt_df$Exposure <- ifelse(crypt_df$Local_Mac_EGF_Load > med_load, "Near_HIGH_EGF_Macs", "Near_LOW_EGF_Macs")

exposure_summary <- aggregate(SREBP2_synth1 ~ Exposure, data = crypt_df, FUN = mean)
print(exposure_summary)









library(CellChat); library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
cc_list <- readRDS("cc_plasma_all.rds")
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
names(cc_list) <- short[names(cc_list)]

# ---- extract just the two we want, compute centrality ----
ctrl <- netAnalysis_computeCentrality(cc_list[["Control"]], slot.name="netP")
homo <- netAnalysis_computeCentrality(cc_list[["Braf_homo"]], slot.name="netP")
pair <- mergeCellChat(list(Control=ctrl, Braf_homo=homo), add.names=c("Control","Braf_homo"))

# consistent cell-type colors
grp_cols <- c("Crypt"="#6A3D9A","Junction"="#F4C430","Villus"="#2C6FBB","Villus_tip"="#C0392B",
              "Fibroblast"="#1B9E77","Macrophage"="#D95F02","Tcell"="#7570B3","DC"="#E7298A","Plasma"="#66A61E")


top <- pw %>% arrange(desc(diff))
top <- rbind(head(top,10), tail(top,8))
top$pathway_name <- factor(top$pathway_name, levels=top$pathway_name)
pdf("BvC_6_pathway_bar.pdf", width=7, height=7)
print(ggplot(top, aes(diff, pathway_name, fill=diff>0)) +
        geom_col() +
        scale_fill_manual(values=c("TRUE"="#B2182B","FALSE"="#2166AC"),
                          labels=c("TRUE"="up in Braf_homo","FALSE"="up in Control"), name=NULL) +
        labs(x="Δ signaling strength (Braf_homo − Control)", y=NULL,
             title="Pathways most changed in Braf_homo vs Control") +
        theme_classic(base_size=12) + theme(plot.title=element_text(face="bold")))
dev.off()

cat("saved 6 Braf_homo-vs-Control figures\n")








library(CellChat); library(Seurat); library(dplyr); library(tidyr); library(ggplot2)
cc_list <- readRDS("cc_plasma_all.rds")
obj <- readRDS("obj_annotated.rds")
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
names(cc_list) <- short[names(cc_list)]
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")
cc_list <- cc_list[dose_order]
lr_all <- bind_rows(Map(function(cc,g){d<-subsetCommunication(cc);d$genotype<-g;d}, cc_list, names(cc_list)))

# ===== 1. what L-R pairs make up the JAM pathway? =====
jam <- lr_all %>% filter(pathway_name=="JAM")
cat("=== JAM ligand-receptor pairs ===\n")
print(jam %>% distinct(interaction_name, ligand, receptor))
cat("\n=== the genes involved ===\n")
jam_genes <- unique(c(jam$ligand, jam$receptor))
# CellChat stores complexes; split any underscores
jam_genes <- unique(unlist(strsplit(jam_genes, "_")))
print(jam_genes)

# ===== 2. JAM L-R pairs into villus tip, across genotypes =====
cat("\n=== JAM interactions INTO villus tip x genotype ===\n")
print(jam %>% filter(target=="Villus_tip") %>%
        group_by(interaction_name, genotype) %>% summarise(prob=round(sum(prob),4),.groups="drop") %>%
        pivot_wider(names_from=genotype, values_from=prob, values_fill=0) %>%
        as.data.frame(), row.names=FALSE)

# ===== 3. which cell pairs drive JAM at the tip (Braf_homo) =====
cat("\n=== JAM source->target at tip (Braf_homo) ===\n")
print(jam %>% filter(genotype=="Braf_homo", target=="Villus_tip") %>%
        group_by(source, interaction_name) %>% summarise(prob=round(sum(prob),4),.groups="drop") %>%
        arrange(desc(prob)) %>% head(15) %>% as.data.frame(), row.names=FALSE)

# ===== 4. RAW EXPRESSION of the JAM genes in villus tip across genotypes =====
jam_genes <- intersect(jam_genes, rownames(obj))
cat("\nJAM genes in data:", paste(jam_genes, collapse=", "), "\n")
tip <- subset(obj, celltype=="Villus_tip")
e <- FetchData(tip, vars=c(jam_genes,"geno"))
t <- e %>% group_by(geno) %>% summarise(across(all_of(jam_genes), ~round(mean(.x),3)),.groups="drop")
tt <- as.data.frame(t(t[,-1])); colnames(tt) <- t$geno
cat("\n=== JAM gene raw expression in VILLUS TIP x genotype ===\n")
print(tt[, intersect(dose_order, colnames(tt)), drop=FALSE])

# ===== 5. % of tip cells expressing each JAM gene (capture check) =====
cat("\n=== % of tip cells expressing each JAM gene ===\n")
print(sapply(FetchData(tip, vars=jam_genes), function(x) round(100*mean(x>0),1)))








library(CellChat); library(Seurat); library(dplyr); library(tidyr)
cc_list <- readRDS("cc_plasma_all.rds")
obj <- readRDS("obj_annotated.rds")
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
names(cc_list) <- short[names(cc_list)]
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")
cc_list <- cc_list[dose_order]
obj$geno <- factor(obj$geno, levels=dose_order)
lr_all <- bind_rows(Map(function(cc,g){d<-subsetCommunication(cc);d$genotype<-g;d}, cc_list, names(cc_list)))
db <- CellChatDB.mouse$interaction
comp <- CellChatDB.mouse$complex

pathways <- c("JAM","APP","MHC-II","LAMININ","GALECTIN","WNT","VEGF","HSPG")
tip <- subset(obj, celltype=="Villus_tip")

for(pw in pathways){
  cat("\n\n########################  ", pw, "  ########################\n")
  
  # ---- PART A: L-R pairs ranked by total communication probability ----
  d <- lr_all %>% filter(pathway_name==pw)
  if(nrow(d)>0){
    lr <- d %>% group_by(interaction_name, ligand, receptor, genotype) %>%
      summarise(prob=sum(prob),.groups="drop") %>%
      pivot_wider(names_from=genotype, values_from=prob, values_fill=0)
    for(g in dose_order) if(!g %in% names(lr)) lr[[g]] <- 0
    lr$total <- rowSums(lr[,dose_order])
    lr <- lr %>% arrange(desc(total))
    lr <- lr[, c("interaction_name","ligand","receptor", dose_order)]
    lr[dose_order] <- round(lr[dose_order],4)
    cat("\n--- L-R PAIRS (ranked by total signaling) ---\n")
    print(as.data.frame(head(lr,10)), row.names=FALSE)
    write.csv(lr, paste0("LR_", pw, ".csv"), row.names=FALSE)
  } else cat("\n[no L-R pairs found]\n")
  
  # ---- PART B: individual genes (from DB) ranked by tip expression ----
  rows <- db[db$pathway_name==pw,]
  genes <- unique(c(rows$ligand, rows$receptor))
  genes <- unlist(lapply(genes, function(g) if(g %in% rownames(comp)) unlist(comp[g,],use.names=FALSE) else g))
  genes <- unique(genes[!is.na(genes) & genes!=""])
  genes <- intersect(genes, rownames(obj))
  if(length(genes)>0){
    e <- FetchData(tip, vars=c(genes,"geno"))
    gt <- e %>% group_by(geno) %>% summarise(across(all_of(genes), ~mean(.x)),.groups="drop")
    gtt <- as.data.frame(t(gt[,-1])); colnames(gtt) <- gt$geno
    gtt <- gtt[, intersect(dose_order, colnames(gtt)), drop=FALSE]
    gtt$overall <- round(rowMeans(gtt),3)
    gtt <- gtt[order(-gtt$overall),]
    gtt <- round(gtt,3)
    gtt$pct_cells <- sapply(FetchData(tip, vars=rownames(gtt)), function(x) round(100*mean(x>0),1))[rownames(gtt)]
    cat("\n--- GENES (ranked by tip expression, with % capture) ---\n")
    print(head(gtt,10))
    write.csv(head(gtt,10), paste0("GENES_", pw, ".csv"))
  } else cat("\n[no genes found in data]\n")
}
cat("\n\nDONE — LR_<pathway>.csv and GENES_<pathway>.csv saved for each\n")






library(CellChat); library(Seurat); library(dplyr); library(tidyr)
cc_list <- readRDS("cc_plasma_all.rds")
obj <- readRDS("obj_annotated.rds")
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
obj$geno <- unname(short[obj$genotype])   # unname() is the fix
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")
db <- CellChatDB.mouse$interaction; comp <- CellChatDB.mouse$complex
CAPTURE_MIN <- 10
tip <- subset(obj, celltype=="Villus_tip")

rank_by_fc <- function(pw){
  rows <- db[db$pathway_name==pw,]
  genes <- unique(c(rows$ligand, rows$receptor))
  genes <- unlist(lapply(genes, function(g) if(g %in% rownames(comp)) unlist(comp[g,],use.names=FALSE) else g))
  genes <- unique(genes[!is.na(genes) & genes!=""]); genes <- intersect(genes, rownames(obj))
  if(length(genes)==0){ cat("\n[",pw,"] no genes\n"); return(invisible()) }
  e <- FetchData(tip, vars=c(genes,"geno"))
  m <- e %>% group_by(geno) %>% summarise(across(all_of(genes), ~mean(.x)),.groups="drop")
  mm <- as.data.frame(t(m[,-1])); colnames(mm) <- m$geno
  pct <- sapply(FetchData(tip, vars=genes), function(x) 100*mean(x>0))
  mm$pct_cells <- round(pct[rownames(mm)],1)
  mm$Braf_homo <- round(mm$Braf_homo,3); mm$Control <- round(mm$Control,3)
  mm$log2FC <- round(log2((mm$Braf_homo+0.01)/(mm$Control+0.01)),2)
  keep <- mm[mm$pct_cells >= CAPTURE_MIN, ]
  keep <- keep[order(-abs(keep$log2FC)), ]
  out <- keep[, c("pct_cells","Control","Braf_homo","log2FC")]
  cat("\n==================", pw, "— well-captured genes (>=",CAPTURE_MIN,"% cells), ranked by |log2FC| ==================\n")
  if(nrow(out)==0) cat("  (no genes pass capture threshold)\n") else print(out)
  write.csv(out, paste0("FC_", pw, ".csv"))
  invisible(out)
}

for(pw in c("JAM","APP","MHC-II","LAMININ","GALECTIN","WNT","VEGF","HSPG")){
  rank_by_fc(pw)
}
cat("\nDONE\n")











library(CellChat); library(Seurat); library(dplyr); library(tidyr)
obj <- readRDS("obj_annotated.rds")
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
obj$geno <- unname(short[obj$genotype])
db <- CellChatDB.mouse$interaction; comp <- CellChatDB.mouse$complex
tip <- subset(obj, celltype=="Villus_tip")

all_genes_fc <- function(pw){
  rows <- db[db$pathway_name==pw,]
  genes <- unique(c(rows$ligand, rows$receptor))
  genes <- unlist(lapply(genes, function(g) if(g %in% rownames(comp)) unlist(comp[g,],use.names=FALSE) else g))
  genes <- unique(genes[!is.na(genes) & genes!=""]); genes <- intersect(genes, rownames(obj))
  if(length(genes)==0){ cat("\n[",pw,"] no genes\n"); return(invisible()) }
  e <- FetchData(tip, vars=c(genes,"geno"))
  m <- e %>% group_by(geno) %>% summarise(across(all_of(genes), ~mean(.x)),.groups="drop")
  mm <- as.data.frame(t(m[,-1])); colnames(mm) <- m$geno
  mm$pct_cells <- round(sapply(FetchData(tip, vars=genes), function(x) 100*mean(x>0))[rownames(mm)],1)
  mm$role <- ifelse(rownames(mm) %in% rows$ligand, "ligand",
                    ifelse(rownames(mm) %in% rows$receptor, "receptor", "both/complex"))
  mm$Control <- round(mm$Control,3); mm$Braf_homo <- round(mm$Braf_homo,3)
  mm$log2FC <- round(log2((mm$Braf_homo+0.01)/(mm$Control+0.01)),2)
  out <- mm[order(-mm$log2FC), c("role","pct_cells","Control","Braf_homo","log2FC")]
  cat("\n==================", pw, "— ALL genes, sorted by log2FC (Braf_homo vs Control) ==================\n")
  print(out)
  write.csv(out, paste0("ALLGENES_", pw, ".csv"))
  invisible(out)
}

for(pw in c("JAM","APP","MHC-II","LAMININ","GALECTIN","WNT","VEGF","HSPG")) all_genes_fc(pw)
cat("\nDONE — ALLGENES_<pathway>.csv saved\n")







library(CellChat)
library(Seurat)
library(dplyr)
library(tidyr)

# Function to extract only the genes DRIVING the pathway differences
get_pathway_drivers <- function(obj, pathways, logfc_thresh = 0.5, min_pct = 5) {
  db <- CellChatDB.mouse$interaction
  comp <- CellChatDB.mouse$complex
  
  # Ensure we are looking at the Villus tip as in your previous code
  target_cells <- subset(obj, celltype == "Villus_tip")
  
  results_list <- list()
  
  for (pw in pathways) {
    rows <- db[db$pathway_name == pw, ]
    if(nrow(rows) == 0) next
    
    # Extract all unique genes for this pathway from CellChat DB
    genes <- unique(c(rows$ligand, rows$receptor))
    genes <- unlist(lapply(genes, function(g) if(g %in% rownames(comp)) unlist(comp[g,], use.names=FALSE) else g))
    genes <- unique(genes[!is.na(genes) & genes != ""])
    
    # Reality check: only keep genes actually present in your Seurat object
    genes <- intersect(genes, rownames(target_cells))
    if(length(genes) == 0) next
    
    # Fetch expression data
    e <- FetchData(target_cells, vars = c(genes, "geno"))
    
    # Calculate Mean Expression
    m <- e %>% 
      filter(geno %in% c("Control", "Braf_homo")) %>%
      group_by(geno) %>% 
      summarise(across(all_of(genes), ~mean(.x)), .groups="drop")
    
    if(nrow(m) < 2) next # Skip if a genotype is missing
    
    mm <- as.data.frame(t(m[,-1]))
    colnames(mm) <- m$geno
    
    # Calculate percent capture (Max between Control and Braf_homo)
    pct_b_homo <- sapply(FetchData(subset(target_cells, geno == "Braf_homo"), vars = genes), function(x) 100 * mean(x > 0))
    pct_ctrl <- sapply(FetchData(subset(target_cells, geno == "Control"), vars = genes), function(x) 100 * mean(x > 0))
    mm$pct_cells_max <- round(pmax(pct_b_homo[rownames(mm)], pct_ctrl[rownames(mm)]), 1)
    
    # Calculate Log2FC (using a small pseudocount to prevent infinity)
    mm$Control <- round(mm$Control, 3)
    mm$Braf_homo <- round(mm$Braf_homo, 3)
    mm$log2FC <- round(log2((mm$Braf_homo + 0.01) / (mm$Control + 0.01)), 2)
    
    # Annotate metadata
    mm$Pathway <- pw
    mm$Gene <- rownames(mm)
    mm$Role <- ifelse(mm$Gene %in% rows$ligand, "ligand", 
                      ifelse(mm$Gene %in% rows$receptor, "receptor", "both/complex"))
    
    # THE FILTER: Keep only meaningful drivers
    drivers <- mm %>%
      filter(abs(log2FC) >= logfc_thresh & pct_cells_max >= min_pct) %>%
      arrange(desc(log2FC)) %>%
      select(Pathway, Gene, Role, pct_cells_max, Control, Braf_homo, log2FC)
    
    if(nrow(drivers) > 0) {
      results_list[[pw]] <- drivers
    }
  }
  
  # Combine all pathways into one clean dataframe
  final_df <- bind_rows(results_list)
  return(final_df)
}

# Run the function on your pathways of interest
pathways_to_test <- c("JAM", "APP", "MHC-II", "LAMININ", "GALECTIN", "WNT", "VEGF", "HSPG")

# You can adjust logfc_thresh or min_pct here if you want more or fewer genes
driver_genes_df <- get_pathway_drivers(obj, pathways_to_test, logfc_thresh = 0.5, min_pct = 5.0)

# Display in console
cat("\n=== PATHWAY DRIVERS (Braf_homo vs Control) ===\n")
print(as.data.frame(driver_genes_df), row.names = FALSE)

# Save to a single clean CSV for the PI
write.csv(driver_genes_df, "Pathway_Drivers_Braf_vs_Control.csv", row.names = FALSE)
cat("\nSaved to 'Pathway_Drivers_Braf_vs_Control.csv'\n")




library(CellChat)
library(Seurat)
library(dplyr)
library(tidyr)

# Prepare the data
db <- CellChatDB.mouse$interaction
comp <- CellChatDB.mouse$complex
tip <- subset(obj, celltype == "Villus_tip")

pathways_to_test <- c("JAM", "APP", "MHC-II", "LAMININ", "GALECTIN", "WNT", "VEGF", "HSPG")
results_list <- list()

for (pw in pathways_to_test) {
  rows <- db[db$pathway_name == pw, ]
  if(nrow(rows) == 0) next
  
  # Extract genes and unpack complexes
  genes <- unique(c(rows$ligand, rows$receptor))
  genes <- unlist(lapply(genes, function(g) if(g %in% rownames(comp)) unlist(comp[g,],use.names=FALSE) else g))
  genes <- unique(genes[!is.na(genes) & genes != ""])
  genes <- intersect(genes, rownames(obj))
  
  if(length(genes) == 0) next
  
  # Fetch data and calculate means
  e <- FetchData(tip, vars = c(genes, "geno"))
  m <- e %>% group_by(geno) %>% summarise(across(all_of(genes), ~mean(.x)), .groups="drop")
  
  mm <- as.data.frame(t(m[,-1]))
  colnames(mm) <- m$geno
  
  # Calculate % detection (pct_cells)
  mm$pct_cells <- round(sapply(FetchData(tip, vars=genes), function(x) 100*mean(x>0))[rownames(mm)], 1)
  
  # Annotate roles
  mm$role <- ifelse(rownames(mm) %in% rows$ligand, "ligand", 
                    ifelse(rownames(mm) %in% rows$receptor, "receptor", "both/complex"))
  
  # Format calculations and Log2FC
  mm$Control <- round(mm$Control, 3)
  mm$Braf_homo <- round(mm$Braf_homo, 3)
  mm$log2FC <- round(log2((mm$Braf_homo + 0.01) / (mm$Control + 0.01)), 2)
  
  # Add Pathway and Gene columns, sort by fold-change
  out <- mm %>%
    mutate(Pathway = pw, Gene = rownames(mm)) %>%
    arrange(desc(log2FC)) %>%
    select(Pathway, Gene, role, pct_cells, Control, Braf_homo, log2FC)
  
  results_list[[pw]] <- out
}

# Merge everything into ONE master dataframe
master_unfiltered_df <- bind_rows(results_list)

# Print a preview and save the single file
cat("\n=== PREVIEW OF MASTER FILE ===\n")
print(head(master_unfiltered_df, 10))

write.csv(master_unfiltered_df, "MASTER_ALL_GENES_Unfiltered.csv", row.names = FALSE)
cat("\n=======================================================\n")
cat("SUCCESS: All pathways saved to 'MASTER_ALL_GENES_Unfiltered.csv'\n")
cat("=======================================================\n")





library(Seurat)
library(dplyr)
library(CellChat)

# 1. Define the 8 pathways and isolate the exact target genes
pathways_to_test <- c("JAM", "APP", "MHC-II", "LAMININ", "GALECTIN", "WNT", "VEGF", "HSPG")
db <- CellChatDB.mouse$interaction
comp <- CellChatDB.mouse$complex

pw_genes <- unique(c(db$ligand[db$pathway_name %in% pathways_to_test],
                     db$receptor[db$pathway_name %in% pathways_to_test]))
pw_genes <- unlist(lapply(pw_genes, function(g) if(g %in% rownames(comp)) unlist(comp[g,],use.names=FALSE) else g))
pw_genes <- unique(pw_genes[!is.na(pw_genes) & pw_genes != ""])
pw_genes <- intersect(pw_genes, rownames(obj))

# 2. Subset the Seurat object strictly to the two groups we are comparing
obj_subset <- subset(obj, geno %in% c("Control", "Braf_homo"))
Idents(obj_subset) <- "geno"

celltypes <- unique(obj_subset$celltype)
results_list <- list()

cat("Calculating True Statistical DEGs (Braf_homo vs Control) across ALL cell types...\n")

for(ct in celltypes) {
  cat("Testing:", ct, "...\n")
  
  # Isolate this specific cell type
  ct_obj <- subset(obj_subset, celltype == ct)
  
  # Safety check: Ensure there are enough cells to run statistics
  if(sum(ct_obj$geno == "Braf_homo") < 3 | sum(ct_obj$geno == "Control") < 3) {
    cat("  -> Skipping (Not enough cells)\n")
    next
  }
  
  # 3. EXACT CONSTRAINT ENFORCED HERE:
  # The Wilcoxon Rank Sum test is strictly run on Braf_homo vs Control
  degs <- FindMarkers(ct_obj, 
                      ident.1 = "Braf_homo", 
                      ident.2 = "Control", 
                      features = pw_genes,
                      logfc.threshold = 0, # Keep at 0 to capture everything
                      min.pct = 0,         # Keep at 0 to capture everything
                      verbose = FALSE)
  
  if(nrow(degs) > 0) {
    degs$Gene <- rownames(degs)
    degs$CellType <- ct
    results_list[[ct]] <- degs
  }
}

# 4. Combine into one master presentation table
master_stats <- bind_rows(results_list) %>%
  select(CellType, Gene, p_val, p_val_adj, avg_log2FC, pct.1, pct.2) %>%
  arrange(CellType, p_val_adj, desc(abs(avg_log2FC)))

cat("\n=== TRUE STATISTICAL DEGs (Example Preview: Tnfrsf21 across tissue) ===\n")
print(master_stats %>% filter(Gene == "Tnfrsf21") %>% arrange(p_val_adj))

write.csv(master_stats, "ALL_CellTypes_Statistical_DEGs_Braf_vs_Control.csv", row.names = FALSE)
cat("\nSUCCESS: Master table saved to 'ALL_CellTypes_Statistical_DEGs_Braf_vs_Control.csv'\n")





library(Seurat)
library(dplyr)
library(CellChat)

pathways_to_test <- c("JAM", "APP", "MHC-II", "LAMININ", "GALECTIN", "WNT", "VEGF", "HSPG")
db <- CellChatDB.mouse$interaction
comp <- CellChatDB.mouse$complex

# 1. Build a Master Map of Genes to Pathways
pathway_map_list <- list()

for(pw in pathways_to_test) {
  rows <- db[db$pathway_name == pw, ]
  if(nrow(rows) == 0) next
  
  genes <- unique(c(rows$ligand, rows$receptor))
  genes <- unlist(lapply(genes, function(g) if(g %in% rownames(comp)) unlist(comp[g,],use.names=FALSE) else g))
  genes <- unique(genes[!is.na(genes) & genes != ""])
  genes <- intersect(genes, rownames(obj))
  
  if(length(genes) > 0) {
    pathway_map_list[[pw]] <- data.frame(Pathway = pw, Gene = genes)
  }
}

# Combine and handle genes that exist in multiple pathways (e.g., Itga3)
pathway_map <- bind_rows(pathway_map_list) %>%
  group_by(Gene) %>%
  summarise(Pathway = paste(Pathway, collapse = " / "), .groups = "drop")

pw_genes <- unique(pathway_map$Gene)

# 2. Subset object
obj_subset <- subset(obj, geno %in% c("Control", "Braf_homo"))
Idents(obj_subset) <- "geno"

celltypes <- unique(obj_subset$celltype)
results_list <- list()

cat("Calculating Statistical DEGs (Braf_homo vs Control) WITH Pathway Mapping...\n")

for(ct in celltypes) {
  cat("Testing:", ct, "...\n")
  
  ct_obj <- subset(obj_subset, celltype == ct)
  
  if(sum(ct_obj$geno == "Braf_homo") < 3 | sum(ct_obj$geno == "Control") < 3) {
    cat("  -> Skipping (Not enough cells)\n")
    next
  }
  
  degs <- FindMarkers(ct_obj, 
                      ident.1 = "Braf_homo", 
                      ident.2 = "Control", 
                      features = pw_genes,
                      logfc.threshold = 0,
                      min.pct = 0,
                      verbose = FALSE)
  
  if(nrow(degs) > 0) {
    degs$Gene <- rownames(degs)
    degs$CellType <- ct
    
    # 3. Glue the Pathway name back onto the statistical results
    degs <- left_join(degs, pathway_map, by = "Gene")
    results_list[[ct]] <- degs
  }
}

# 4. Master presentation table: Filtered and grouped by Pathway
master_stats <- bind_rows(results_list) %>%
  select(Pathway, CellType, Gene, p_val, p_val_adj, avg_log2FC, pct.1, pct.2) %>%
  # THE CUTOFF: Keep only statistically significant results (FDR < 0.05)
  filter(p_val_adj < 0.05) %>%
  # THE SORTING: Group by Pathway first, then rank by significance
  arrange(Pathway, p_val_adj, desc(abs(avg_log2FC)))

cat("\n=== SIGNIFICANT DEGs (Ranked by p-value within each Pathway) ===\n")
print(head(master_stats, 15))

write.csv(master_stats, "Significant_Pathway_DEGs_Ranked.csv", row.names = FALSE)
cat("\nSUCCESS: Filtered and ranked table saved to 'Significant_Pathway_DEGs_Ranked.csv'\n")





library(CellChat); library(dplyr); library(tidyr)
cc_list <- readRDS("cc_plasma_all.rds")
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
names(cc_list) <- short[names(cc_list)]
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")
cc_list <- cc_list[dose_order]
lr_all <- bind_rows(Map(function(cc,g){d<-subsetCommunication(cc);d$genotype<-g;d}, cc_list, names(cc_list)))

# for each pathway: top sender -> receiver pairs, in Braf_homo and Control (to see the change)
who_talks <- function(pw){
  d <- lr_all %>% filter(pathway_name==pw)
  if(nrow(d)==0){ cat("\n[",pw,"] none\n"); return(invisible()) }
  # aggregate by sender->receiver per genotype
  st <- d %>% group_by(source, target, genotype) %>% summarise(prob=sum(prob),.groups="drop") %>%
    pivot_wider(names_from=genotype, values_from=prob, values_fill=0)
  for(g in dose_order) if(!g %in% names(st)) st[[g]] <- 0
  st$total <- rowSums(st[,dose_order])
  st <- st %>% arrange(desc(total)) %>%
    mutate(pair = paste(source, "->", target))
  st <- st[, c("pair", dose_order)]
  st[dose_order] <- round(st[dose_order],4)
  cat("\n================== ", pw, " — top sender -> receiver pairs (dose order) ==================\n")
  print(as.data.frame(head(st, 10)), row.names=FALSE)
  write.csv(st, paste0("WHO_", pw, ".csv"), row.names=FALSE)
  invisible(st)
}

for(pw in c("JAM","APP","MHC-II","LAMININ","GALECTIN","WNT","VEGF")){
  who_talks(pw)
}
cat("\nDONE — WHO_<pathway>.csv saved\n")






library(CellChat); library(Seurat); library(dplyr); library(tidyr)
cc_list <- readRDS("cc_plasma_all.rds")
obj <- readRDS("obj_annotated.rds")
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
names(cc_list) <- short[names(cc_list)]
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")
cc_list <- cc_list[dose_order]
lr_all <- bind_rows(Map(function(cc,g){d<-subsetCommunication(cc);d$genotype<-g;d}, cc_list, names(cc_list)))

# ===== 1. ALL pathways Fibroblast -> Crypt, ranked, across genotypes =====
fc <- lr_all %>% filter(source=="Fibroblast", target=="Crypt") %>%
  group_by(pathway_name, genotype) %>% summarise(prob=sum(prob),.groups="drop") %>%
  pivot_wider(names_from=genotype, values_from=prob, values_fill=0)
for(g in dose_order) if(!g %in% names(fc)) fc[[g]] <- 0
fc$total <- rowSums(fc[,dose_order])
fc <- fc %>% arrange(desc(total))
fc[dose_order] <- round(fc[dose_order],4)
cat("=== ALL Fibroblast -> Crypt pathways (ranked, dose order) ===\n")
print(as.data.frame(fc[, c("pathway_name", dose_order)]), row.names=FALSE)

# ===== 2. the specific L-R pairs for the top fibroblast->crypt pathways =====
cat("\n=== top Fibroblast -> Crypt L-R pairs (Braf_homo) ===\n")
print(lr_all %>% filter(source=="Fibroblast", target=="Crypt", genotype=="Braf_homo") %>%
        group_by(interaction_name, ligand, receptor) %>% summarise(prob=round(sum(prob),4),.groups="drop") %>%
        arrange(desc(prob)) %>% head(12) %>% as.data.frame(), row.names=FALSE)

# ===== 3. is Fibroblast->Crypt specifically Braf-driven or just tissue-wide? =====
# compare fibroblast->crypt total vs fibroblast->ALL, to see if crypt is special
cat("\n=== Fibroblast total OUTGOING signaling by target, per genotype ===\n")
print(lr_all %>% filter(source=="Fibroblast") %>%
        group_by(target, genotype) %>% summarise(prob=sum(prob),.groups="drop") %>%
        pivot_wider(names_from=genotype, values_from=prob, values_fill=0) %>%
        mutate(total=rowSums(across(all_of(dose_order)))) %>% arrange(desc(total)) %>%
        select(target, all_of(dose_order)) %>% as.data.frame() %>% head(10), row.names=FALSE)





library(Seurat); library(CellChat); library(dplyr); library(tidyr)
obj <- readRDS("obj_annotated.rds")
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
obj$geno <- unname(short[obj$genotype])
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")

# cholesterol module
chol <- intersect(c("Srebf2","Hmgcr","Hmgcs1","Mvk","Mvd","Pmvk","Fdps","Fdft1","Sqle","Lss",
                    "Cyp51","Dhcr7","Dhcr24","Idi1","Insig1","Ldlr"), rownames(obj))
obj <- AddModuleScore(obj, features=list(chol), name="Chol")

# ===== APPROACH A: per-cell correlation of cholesterol vs each pathway's RECEPTOR module, in crypt =====
# build receptor-gene modules per pathway from CellChatDB (receptors = what the crypt cell "receives")
db <- CellChatDB.mouse$interaction; comp <- CellChatDB.mouse$complex
pathways <- unique(db$pathway_name)
get_receptors <- function(pw){
  r <- db$receptor[db$pathway_name==pw]
  g <- unlist(lapply(r, function(x) if(x %in% rownames(comp)) unlist(comp[x,],use.names=FALSE) else x))
  unique(intersect(g[!is.na(g)&g!=""], rownames(obj)))
}

crypt <- subset(obj, celltype=="Crypt")
crypt_bh <- subset(crypt, geno=="Braf_homo")   # test within Braf_homo (where cholesterol is high)

results <- data.frame()
for(pw in pathways){
  recs <- get_receptors(pw)
  if(length(recs)<1) next
  # only use receptors captured in >5% of crypt cells (avoid noise)
  pct <- sapply(FetchData(crypt_bh, vars=recs), function(x) mean(x>0))
  recs <- recs[pct>0.05]
  if(length(recs)<1) next
  crypt_bh <- AddModuleScore(crypt_bh, features=list(recs), name="tmpRec")
  d <- FetchData(crypt_bh, vars=c("Chol1","tmpRec1"))
  r <- cor(d$Chol1, d$tmpRec1, method="spearman")
  results <- rbind(results, data.frame(pathway=pw, n_receptors=length(recs), r=round(r,3)))
  crypt_bh$tmpRec1 <- NULL
}
results <- results %>% arrange(desc(r))
cat("=== pathways whose RECEPTOR expression correlates with cholesterol in Braf_homo crypt cells ===\n")
cat("(top = positively correlated with cholesterol program)\n\n")
print(head(results, 15), row.names=FALSE)
cat("\n--- most NEGATIVELY correlated ---\n")
print(tail(results, 8), row.names=FALSE)
write.csv(results, "chol_pathway_correlations.csv", row.names=FALSE)



library(Seurat); library(CellChat); library(dplyr); library(tidyr)
obj <- readRDS("obj_annotated.rds")
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
obj$geno <- unname(short[obj$genotype])
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")

chol <- intersect(c("Srebf2","Hmgcr","Hmgcs1","Mvk","Mvd","Pmvk","Fdps","Fdft1","Sqle","Lss",
                    "Cyp51","Dhcr7","Dhcr24","Idi1","Insig1","Ldlr"), rownames(obj))
obj <- AddModuleScore(obj, features=list(chol), name="Chol")

db <- CellChatDB.mouse$interaction; comp <- CellChatDB.mouse$complex
pathways <- unique(db$pathway_name)
get_receptors <- function(pw){
  r <- db$receptor[db$pathway_name==pw]
  g <- unlist(lapply(r, function(x) if(x %in% rownames(comp)) unlist(comp[x,],use.names=FALSE) else x))
  unique(intersect(g[!is.na(g)&g!=""], rownames(obj)))
}

# screen within a given cell subset, correlate cholesterol vs each pathway's receptor module
screen_corr <- function(cells_obj, label){
  res <- data.frame()
  for(pw in pathways){
    recs <- get_receptors(pw)
    if(length(recs)<1) next
    pct <- sapply(FetchData(cells_obj, vars=recs), function(x) mean(x>0))
    recs <- recs[pct>0.05]
    if(length(recs)<1) next
    tmp <- AddModuleScore(cells_obj, features=list(recs), name="tmpRec")
    d <- FetchData(tmp, vars=c("Chol1","tmpRec1"))
    res <- rbind(res, data.frame(pathway=pw, n_rec=length(recs), r=round(cor(d$Chol1,d$tmpRec1,method="spearman"),3)))
  }
  res <- res %>% arrange(desc(r))
  cat("\n===== ", label, " (top 12 positive) =====\n")
  print(head(res,12), row.names=FALSE)
  cat("--- (top 5 negative) ---\n")
  print(tail(res,5), row.names=FALSE)
  res$subset <- label
  res
}

# run across: all epithelium, and each zone, all within Braf_homo
all_res <- list()
epi_bh <- subset(obj, celltype %in% c("Crypt","Junction","Villus","Villus_tip") & geno=="Braf_homo")
all_res[["All_epithelium"]] <- screen_corr(epi_bh, "ALL EPITHELIUM (Braf_homo)")
for(z in c("Crypt","Junction","Villus","Villus_tip")){
  sub <- subset(obj, celltype==z & geno=="Braf_homo")
  if(ncol(sub)>50) all_res[[z]] <- screen_corr(sub, paste0(z, " (Braf_homo)"))
}

# combine: which pathways are consistently correlated across zones?
combined <- bind_rows(all_res) %>%
  select(pathway, subset, r) %>%
  pivot_wider(names_from=subset, values_from=r)
write.csv(combined, "chol_corr_by_zone.csv", row.names=FALSE)
cat("\n===== pathways x zone correlation table (saved) =====\n")
# show pathways with high mean correlation across zones
combined$mean_r <- round(rowMeans(combined[,-1], na.rm=TRUE),3)
print(as.data.frame(combined %>% arrange(desc(mean_r)) %>% head(15)), row.names=FALSE)






library(Seurat); library(CellChat); library(dplyr); library(tidyr)
obj <- readRDS("obj_annotated.rds")
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
obj$geno <- unname(short[obj$genotype])
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")

# YAP/TAZ target-gene signature (canonical YAP/TEAD targets) + the components
yap_targets <- intersect(c("Ctgf","Ccn2","Cyr61","Ccn1","Ankrd1","Amotl2","Cav1","Tead1",
                           "Thbs1","Axl","Igfbp3","Serpine1","Cdc6"), rownames(obj))
yap_core <- intersect(c("Yap1","Wwtr1","Tead1","Tead2","Tead4"), rownames(obj))  # Wwtr1=TAZ
chol <- intersect(c("Srebf2","Hmgcr","Hmgcs1","Mvk","Mvd","Pmvk","Fdps","Fdft1","Sqle","Lss",
                    "Cyp51","Dhcr7","Dhcr24","Idi1","Insig1","Ldlr"), rownames(obj))

obj <- AddModuleScore(obj, features=list(YAP=yap_targets, Chol=chol), name=c("YAP","Chol"))

# ===== 1. YAP target score by genotype, per zone (where is YAP up/down?) =====
cat("=== YAP/TAZ target signature by zone x genotype ===\n")
for(z in c("Crypt","Junction","Villus","Villus_tip")){
  sub <- subset(obj, celltype==z)
  d <- FetchData(sub, vars=c("YAP1","geno"))
  r <- d %>% group_by(geno) %>% summarise(YAP=round(mean(YAP1),3),.groups="drop")
  cat("\n--", z, "--\n"); print(as.data.frame(r[match(dose_order,r$geno),]), row.names=FALSE)
}

# ===== 2. YAP core + target genes raw, in crypt, across genotypes =====
crypt <- subset(obj, celltype=="Crypt"); crypt$geno <- factor(crypt$geno, levels=dose_order)
g <- c(yap_core, yap_targets)
e <- FetchData(crypt, vars=c(g,"geno"))
tab <- e %>% group_by(geno) %>% summarise(across(all_of(g), ~round(mean(.x),3)),.groups="drop")
tt <- as.data.frame(t(tab[,-1])); colnames(tt) <- tab$geno
tt$pct <- round(sapply(FetchData(crypt, vars=g), function(x) 100*mean(x>0)),1)[rownames(tt)]
cat("\n=== YAP core + target genes in CRYPT x genotype (+capture) ===\n")
print(tt[, c(intersect(dose_order,colnames(tt)),"pct")])

# 1. Recalculate the module scores separately so they are named correctly
obj <- AddModuleScore(obj, features=list(yap_targets), name="YAP_Score_")
obj <- AddModuleScore(obj, features=list(chol), name="Chol_Score_")

# 2. Re-run the correlation loop
cat("\n=== YAP vs cholesterol correlation, per zone (Braf_homo) ===\n")
for(z in c("Crypt","Junction","Villus","Villus_tip")){
  
  # Subset to the specific zone and genotype
  sub <- subset(obj, celltype==z & geno=="Braf_homo")
  
  # Fetch the newly generated score columns (Seurat automatically appends '1')
  d <- FetchData(sub, vars=c("YAP_Score_1", "Chol_Score_1"))
  
  # Run the spearman correlation
  r_val <- cor(d$YAP_Score_1, d$Chol_Score_1, method="spearman")
  
  # Print the results cleanly
  cat(sprintf("%-12s r = %+.3f (n=%d)\n", z, r_val, nrow(d)))
}



# Load the pre-processed Seurat object
seurat_obj <- readRDS("path/to/your/file.rds")


# 1. Define the downstream targets of the YAP/TAZ pathway
# (Replace with highly specific tissue-context targets if necessary)
yap_taz_targets <- list(c("Ctgf", "Cyr61", "Birc5", "Axl"))

# 2. Score for YAP/TAZ pathway stabilization
seurat_obj <- AddModuleScore(
  object = seurat_obj,
  features = yap_taz_targets,
  name = "YAP_TAZ_Score"
)

# 3. Evaluate the correlation between the Cholesterol score and TAZ upregulation
FeatureScatter(
  seurat_obj,
  feature1 = "Cholesterol_Biosynthesis_Score1",
  feature2 = "YAP_TAZ_Score1",
  group.by = "seurat_clusters"
)

# 4. Assess standard stem cell homeostasis subversion alongside the pathway
FeaturePlot(
  seurat_obj, 
  features = c("Lgr5", "Ascl2", "YAP_TAZ_Score1"), 
  reduction = "umap",
  keep.scale = "all"
)


# 1. Check which folder R is currently looking inside
current_folder <- getwd()
message("R is currently looking in: ", current_folder)

# 2. Get a list of all files in this folder
all_my_files <- list.files()

# 3. Print the list to the screen
message("\nHere are the files in this folder:")
print(all_my_files)


# 1. Check the distribution of your cell types
message("--- CELL TYPE COUNTS ---")
print(table(seurat_obj$celltype))

# 2. Check if the target genes are actually in the dataset
tfmc_apkc_markers <- c("Prkci", "Prkcz", "Anxa10")
present_genes <- intersect(tfmc_apkc_markers, rownames(seurat_obj))
message("\n--- GENES FOUND ---")
print(present_genes)

present_genes <- c("Prkci", "Prkcz", "Anxa10")

table(seurat_obj$genotype)




library(Seurat)
library(DESeq2)

Idents(seurat_obj) <- "celltype"

# --- JUNCTION ZONE ---
junction_de <- FindMarkers(seurat_obj, ident.1 = "VillinCreERT2; BrafV600E/+", ident.2 = "Villin CreERT2(control)", group.by = "genotype", subset.ident = "Junction", test.use = "DESeq2", slot = "counts", logfc.threshold = 0, min.pct = 0)

junction_de$gene <- rownames(junction_de)
junction_de$wald_stat <- sign(junction_de$avg_log2FC) * -log10(junction_de$p_val + 1e-300)
junction_de <- junction_de[order(-junction_de$wald_stat), ]

gene_list_junction <- junction_de$wald_stat
names(gene_list_junction) <- junction_de$gene

# --- CRYPT ZONE ---
crypt_de <- FindMarkers(seurat_obj, ident.1 = "VillinCreERT2; BrafV600E/+", ident.2 = "Villin CreERT2(control)", group.by = "genotype", subset.ident = "Crypt", test.use = "DESeq2", slot = "counts", logfc.threshold = 0, min.pct = 0)

crypt_de$gene <- rownames(crypt_de)
crypt_de$wald_stat <- sign(crypt_de$avg_log2FC) * -log10(crypt_de$p_val + 1e-300)
crypt_de <- crypt_de[order(-crypt_de$wald_stat), ]

gene_list_crypt <- crypt_de$wald_stat
names(gene_list_crypt) <- crypt_de$gene







library(CellChat); library(dplyr); library(tidyr)
cc_list <- readRDS("cc_plasma_all.rds")
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
names(cc_list) <- short[names(cc_list)]
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")
cc_list <- cc_list[dose_order]
lr_all <- bind_rows(Map(function(cc,g){d<-subsetCommunication(cc);d$genotype<-g;d}, cc_list, names(cc_list)))

pathways <- c("JAM","APP","MHC-II","LAMININ","GALECTIN","WNT","VEGF","HSPG")

layer1 <- function(pw){
  d <- lr_all %>% filter(pathway_name==pw)
  if(nrow(d)==0){ cat("\n########", pw, ": not found ########\n"); return(invisible()) }
  cat("\n\n############################", pw, "############################\n")
  
  # sender->receiver strength per genotype
  sr <- d %>% group_by(source, target, genotype) %>% summarise(prob=sum(prob),.groups="drop") %>%
    pivot_wider(names_from=genotype, values_from=prob, values_fill=0)
  for(g in dose_order) if(!g %in% names(sr)) sr[[g]] <- 0
  sr <- sr %>% mutate(pair=paste0(source,"->",target))
  
  # ---- A: top pairs PER GENOTYPE ----
  cat("\n--- TOP 5 sender->receiver pairs, per genotype ---\n")
  for(g in dose_order){
    top <- sr %>% arrange(desc(.data[[g]])) %>% filter(.data[[g]]>0) %>% head(5)
    cat("\n[", g, "]\n")
    if(nrow(top)>0) print(as.data.frame(top %>% transmute(pair, prob=round(.data[[g]],4))), row.names=FALSE)
    else cat("  (no signal)\n")
  }
  
  # ---- B: WHOLE (summed across genotypes), top pairs overall ----
  sr$total <- rowSums(sr[,dose_order])
  cat("\n--- TOP 8 pairs OVERALL (summed all genotypes) ---\n")
  print(as.data.frame(sr %>% arrange(desc(total)) %>% head(8) %>%
                        transmute(pair, !!!setNames(lapply(dose_order,function(g) round(sr %>% arrange(desc(total)) %>% head(8) %>% pull(g),4)), dose_order))), row.names=FALSE)
  
  # ---- C: BRAF_HOMO vs CONTROL contrast ----
  sr$diff <- sr$Braf_homo - sr$Control
  cat("\n--- Braf_homo vs Control: pairs GAINED in Braf_homo (top 6) ---\n")
  print(as.data.frame(sr %>% arrange(desc(diff)) %>% head(6) %>%
                        transmute(pair, Control=round(Control,4), Braf_homo=round(Braf_homo,4), diff=round(diff,4))), row.names=FALSE)
  cat("\n--- pairs LOST in Braf_homo (top 6) ---\n")
  print(as.data.frame(sr %>% arrange(diff) %>% head(6) %>%
                        transmute(pair, Control=round(Control,4), Braf_homo=round(Braf_homo,4), diff=round(diff,4))), row.names=FALSE)
  
  write.csv(sr[,c("pair",dose_order,"total","diff")], paste0("L1_", gsub("-","_",pw), ".csv"), row.names=FALSE)
}

for(pw in pathways) layer1(pw)
cat("\n\nLAYER 1 DONE — L1_<pathway>.csv saved for each\n")








library(CellChat); library(Seurat); library(dplyr); library(tidyr)
cc_list <- readRDS("cc_plasma_all.rds"); obj <- readRDS("obj_annotated.rds")
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
names(cc_list) <- short[names(cc_list)]; obj$geno <- unname(short[obj$genotype])
lr_all <- bind_rows(Map(function(cc,g){d<-subsetCommunication(cc);d$genotype<-g;d}, cc_list, names(cc_list)))

# for a pathway: for its top pairs, show ligand expr in SENDER + receptor expr in RECEIVER,
# in Braf_homo vs Control — so we see if a change is ligand-driven, receptor-driven, or new pair
db <- CellChatDB.mouse$interaction; comp <- CellChatDB.mouse$complex
expand <- function(x) unique(unlist(lapply(x, function(g) if(g %in% rownames(comp)) unlist(comp[g,],use.names=FALSE) else g)))

# mean expr of a gene set in a cell type, in a genotype
gexpr <- function(ct, gt, genes){
  genes <- intersect(genes, rownames(obj))
  if(length(genes)==0) return(NA)
  cells <- colnames(obj)[obj$celltype==ct & obj$geno==gt]
  if(length(cells)<10) return(NA)
  mean(as.matrix(GetAssayData(obj,layer="data")[genes, cells, drop=FALSE]))
}

layer2 <- function(pw, pairs){
  rows <- db[db$pathway_name==pw,]
  ligs <- expand(rows$ligand); recs <- expand(rows$receptor)
  cat("\n\n############", pw, "— ligand vs receptor decomposition ############\n")
  cat("ligand genes:", paste(intersect(ligs,rownames(obj)),collapse=","), "\n")
  cat("receptor genes:", paste(intersect(recs,rownames(obj)),collapse=","), "\n\n")
  out <- data.frame()
  for(p in pairs){
    s <- strsplit(p,"->")[[1]][1]; t <- strsplit(p,"->")[[1]][2]
    out <- rbind(out, data.frame(
      pair=p,
      ligand_Control=round(gexpr(s,"Control",ligs),3),
      ligand_Brafhomo=round(gexpr(s,"Braf_homo",ligs),3),
      receptor_Control=round(gexpr(t,"Control",recs),3),
      receptor_Brafhomo=round(gexpr(t,"Braf_homo",recs),3)))
  }
  out$lig_FC <- round(log2((out$ligand_Brafhomo+.01)/(out$ligand_Control+.01)),2)
  out$rec_FC <- round(log2((out$receptor_Brafhomo+.01)/(out$receptor_Control+.01)),2)
  print(out, row.names=FALSE)
  invisible(out)
}

# feed the key pairs from Layer 1 for each pathway
layer2("JAM", c("Villus_tip->Villus_tip","Junction->Junction","Villus_tip->Villus","Tcell->Tcell"))
layer2("LAMININ", c("Junction->Junction","Crypt->Crypt","Fibroblast->Crypt","Macrophage->Tcell"))
layer2("MHC-II", c("Macrophage->Macrophage","Tcell->Macrophage","DC->DC","Macrophage->DC"))
layer2("APP", c("Macrophage->Macrophage","Tcell->Macrophage","Villus_tip->Villus_tip","DC->DC"))
layer2("GALECTIN", c("Villus_tip->Villus_tip","Junction->Junction","Macrophage->Villus_tip","Crypt->Crypt"))
layer2("VEGF", c("Macrophage->Macrophage","Villus_tip->Macrophage","Tcell->Tcell","Tcell->Macrophage"))
cat("\nLAYER 2 DONE\n")








###############################################################################
##  REPRODUCIBLE CellChat DIFFERENTIAL PIPELINE  —  Braf_homo vs Control
##  Every threshold is a named constant at the top. No "top N" selections.
##  Layers run sequentially and write outputs; nothing needs pasting back.
##
##  METHOD BASIS: CellChat permutation test (Jin et al. 2021), benchmarked as a
##  top CCC method by ESICCC (Luo et al. 2023, Genome Res, doi:10.1101/gr.278001.123).
###############################################################################

suppressPackageStartupMessages({
  library(CellChat); library(Seurat); library(dplyr); library(tidyr); library(ggplot2)
})

## ======================= PARAMETERS (all cutoffs, justified) =================
P_SIG        <- 0.05    # CellChat permutation p-value cutoff for a "significant" interaction
MIN_CELLS    <- 10      # min cells per group for a valid interaction (CellChat default)
FC_PAIR      <- 2       # fold-change to call a sender->receiver pair "gained"/"lost"
PSEUDO       <- 1e-4    # pseudocount for log2FC on communication probabilities
CONTRIB_MIN  <- 0.01    # a pair must contribute >=1% of the pathway total to be reported
WILCOX_P     <- 0.05    # Wilcoxon p cutoff for ligand/receptor expression change
LOG2FC_EXPR  <- 0.58    # |log2FC| >= 0.58 (=1.5x) to call an expression change meaningful
G1 <- "Braf_homo"; G2 <- "Control"   # the two anchor genotypes for the contrast
OUTDIR <- "cellchat_pipeline_out"; dir.create(OUTDIR, showWarnings=FALSE)

## ======================= LOAD =================================================
cc_list <- readRDS("cc_plasma_all.rds")
obj     <- readRDS("obj_annotated.rds")
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
names(cc_list) <- short[names(cc_list)]
obj$geno <- unname(short[obj$genotype])
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")
cc_list <- cc_list[dose_order]

## subsetCommunication returns per-interaction: source,target,ligand,receptor,
## interaction_name, pathway_name, prob, pval  (pval = CellChat permutation p-value)
lr_all <- bind_rows(Map(function(cc,g){
  d <- subsetCommunication(cc); d$genotype <- g; d
}, cc_list, names(cc_list)))

## keep ONLY significant interactions (Step 1: p<P_SIG). min.cells already applied
## upstream in filterCommunication during object construction.
lr_sig <- lr_all %>% filter(pval < P_SIG)
write.csv(lr_sig, file.path(OUTDIR,"00_significant_interactions_all_genotypes.csv"), row.names=FALSE)

###############################################################################
## LAYER 0 — DIFFERENTIAL PATHWAYS (the bar graph), fully reproducible
##  - pathway total prob = sum of SIGNIFICANT L-R probs in that genotype
##  - include pathway if significant in >=1 of the two anchor genotypes
##  - rank by log2FC(Braf_homo / Control)
###############################################################################
pw_geno <- lr_sig %>%
  filter(genotype %in% c(G1,G2)) %>%
  group_by(pathway_name, genotype) %>%
  summarise(prob = sum(prob), n_interactions = n(), .groups="drop") %>%
  pivot_wider(names_from=genotype, values_from=c(prob,n_interactions), values_fill=0)

# ensure both columns exist
for(cn in c(paste0("prob_",G1), paste0("prob_",G2))) if(!cn %in% names(pw_geno)) pw_geno[[cn]] <- 0
pw_geno <- pw_geno %>%
  mutate(prob_G1 = .data[[paste0("prob_",G1)]],
         prob_G2 = .data[[paste0("prob_",G2)]],
         log2FC  = log2((prob_G1+PSEUDO)/(prob_G2+PSEUDO)),
         status  = case_when(prob_G1>0 & prob_G2==0 ~ paste0(G1,"-only"),
                             prob_G2>0 & prob_G1==0 ~ paste0(G2,"-only"),
                             TRUE ~ "both")) %>%
  filter(prob_G1>0 | prob_G2>0) %>%          # EXCLUSION: significant in neither -> dropped
  arrange(desc(log2FC))
write.csv(pw_geno, file.path(OUTDIR,"01_differential_pathways.csv"), row.names=FALSE)

# the bar graph (this IS the reproducible version of the heatmap/bar)
pw_geno$pathway_name <- factor(pw_geno$pathway_name, levels=pw_geno$pathway_name)
p <- ggplot(pw_geno, aes(log2FC, pathway_name, fill=log2FC>0)) +
  geom_col() +
  scale_fill_manual(values=c("TRUE"="#B2182B","FALSE"="#2166AC"),
                    labels=c("TRUE"=paste0("up in ",G1),"FALSE"=paste0("up in ",G2)), name=NULL) +
  labs(x=paste0("log2 FC signaling probability (",G1," / ",G2,")"), y=NULL,
       title=paste0("Differential signaling pathways: ",G1," vs ",G2),
       subtitle=paste0("significant interactions only (perm. p<",P_SIG,"); pathways significant in >=1 genotype")) +
  theme_classic(base_size=11)
ggsave(file.path(OUTDIR,"01_differential_pathways_bar.pdf"), p, width=7, height=9)

cat("LAYER 0 done:", nrow(pw_geno), "pathways pass inclusion (significant in >=1 of",G1,"/",G2,")\n")
cat("  up in",G1,":", sum(pw_geno$log2FC>0), " | up in",G2,":", sum(pw_geno$log2FC<0), "\n")

## PATHWAYS CARRIED FORWARD: those that pass inclusion (no manual picking)
pathways_keep <- as.character(pw_geno$pathway_name)

###############################################################################
## LAYER 1 — sender->receiver structure per pathway (significance-based, not top-N)
##  For each kept pathway, report every sender->receiver pair that is
##  (a) significant (already filtered) AND (b) contributes >= CONTRIB_MIN of the
##  pathway's total prob in that genotype. Classify gained/lost by FC_PAIR.
###############################################################################
layer1_all <- list()
for(pw in pathways_keep){
  d <- lr_sig %>% filter(pathway_name==pw)
  sr <- d %>% group_by(source,target,genotype) %>% summarise(prob=sum(prob),.groups="drop") %>%
    pivot_wider(names_from=genotype, values_from=prob, values_fill=0)
  for(g in dose_order) if(!g %in% names(sr)) sr[[g]] <- 0
  # contribution filter within the two anchor genotypes
  totG1 <- sum(sr[[G1]]); totG2 <- sum(sr[[G2]])
  sr$contrib_G1 <- if(totG1>0) sr[[G1]]/totG1 else 0
  sr$contrib_G2 <- if(totG2>0) sr[[G2]]/totG2 else 0
  sr$pair <- paste0(sr$source,"->",sr$target)
  sr$FC   <- log2((sr[[G1]]+PSEUDO)/(sr[[G2]]+PSEUDO))
  sr$call <- case_when(
    (sr$contrib_G1>=CONTRIB_MIN | sr$contrib_G2>=CONTRIB_MIN) & sr$FC >=  log2(FC_PAIR) ~ "gained_in_Braf_homo",
    (sr$contrib_G1>=CONTRIB_MIN | sr$contrib_G2>=CONTRIB_MIN) & sr$FC <= -log2(FC_PAIR) ~ "lost_in_Braf_homo",
    (sr$contrib_G1>=CONTRIB_MIN | sr$contrib_G2>=CONTRIB_MIN) ~ "stable",
    TRUE ~ "below_contribution_floor")
  sr$pathway <- pw
  layer1_all[[pw]] <- sr %>% select(pathway,pair,source,target,all_of(dose_order),
                                    contrib_G1,contrib_G2,FC,call)
}
layer1_df <- bind_rows(layer1_all)
write.csv(layer1_df, file.path(OUTDIR,"02_sender_receiver_pairs.csv"), row.names=FALSE)
cat("LAYER 1 done: sender->receiver pairs classified for", length(pathways_keep), "pathways\n")

###############################################################################
## LAYER 2 — rewiring vs volume: decompose each GAINED/LOST pair into
##  ligand (in sender) and receptor (in receiver) expression change, with a
##  Wilcoxon test on the actual cells. Classification by WILCOX_P + LOG2FC_EXPR.
###############################################################################
db <- CellChatDB.mouse$interaction; comp <- CellChatDB.mouse$complex
expand_genes <- function(x){
  unique(unlist(lapply(x, function(g) if(g %in% rownames(comp)) unlist(comp[g,],use.names=FALSE) else g)))
}
DefaultAssay(obj) <- "Spatial"
expr_mat <- GetAssayData(obj, layer="data")

# module-mean expression of a gene set in given cells
mean_expr <- function(genes, cells){
  genes <- intersect(genes, rownames(expr_mat))
  if(length(genes)==0 || length(cells)==0) return(NA_real_)
  mean(expr_mat[genes, cells, drop=FALSE])
}
# Wilcoxon on per-cell mean of a gene set, G1 vs G2 cells
wilcox_set <- function(genes, cellsG1, cellsG2){
  genes <- intersect(genes, rownames(expr_mat))
  if(length(genes)<1 || length(cellsG1)<10 || length(cellsG2)<10) return(c(p=NA, l2fc=NA))
  v1 <- colMeans(expr_mat[genes, cellsG1, drop=FALSE])
  v2 <- colMeans(expr_mat[genes, cellsG2, drop=FALSE])
  p  <- tryCatch(wilcox.test(v1, v2)$p.value, error=function(e) NA)
  l2 <- log2((mean(v1)+PSEUDO)/(mean(v2)+PSEUDO))
  c(p=p, l2fc=l2)
}

changed <- layer1_df %>% filter(call %in% c("gained_in_Braf_homo","lost_in_Braf_homo"))
layer2_rows <- list()
for(i in seq_len(nrow(changed))){
  pw <- changed$pathway[i]; s <- changed$source[i]; t <- changed$target[i]
  rows <- db[db$pathway_name==pw,]
  ligs <- expand_genes(rows$ligand); recs <- expand_genes(rows$receptor)
  cs_G1 <- colnames(obj)[obj$celltype==s & obj$geno==G1]  # sender cells, Braf_homo
  cs_G2 <- colnames(obj)[obj$celltype==s & obj$geno==G2]
  ct_G1 <- colnames(obj)[obj$celltype==t & obj$geno==G1]  # receiver cells, Braf_homo
  ct_G2 <- colnames(obj)[obj$celltype==t & obj$geno==G2]
  lig <- wilcox_set(ligs, cs_G1, cs_G2)   # ligand change in sender
  rec <- wilcox_set(recs, ct_G1, ct_G2)   # receptor change in receiver
  lig_sig <- !is.na(lig["p"]) && lig["p"]<WILCOX_P && abs(lig["l2fc"])>=LOG2FC_EXPR
  rec_sig <- !is.na(rec["p"]) && rec["p"]<WILCOX_P && abs(rec["l2fc"])>=LOG2FC_EXPR
  mech <- case_when(
    lig_sig &  rec_sig ~ "ligand+receptor (both change)",
    lig_sig & !rec_sig ~ "ligand-driven (sender makes more)",
    !lig_sig &  rec_sig ~ "receptor-driven (receiver expresses more)",
    TRUE               ~ "neither -> spatial/proximity (test in Layer 3)")
  layer2_rows[[i]] <- data.frame(pathway=pw, pair=changed$pair[i], call=changed$call[i],
                                 ligand_log2FC=round(unname(lig["l2fc"]),2), ligand_p=signif(unname(lig["p"]),3),
                                 receptor_log2FC=round(unname(rec["l2fc"]),2), receptor_p=signif(unname(rec["p"]),3),
                                 mechanism=mech)
}
layer2_df <- bind_rows(layer2_rows)
write.csv(layer2_df, file.path(OUTDIR,"03_rewiring_vs_volume.csv"), row.names=FALSE)
cat("LAYER 2 done: decomposed", nrow(layer2_df), "changed pairs into ligand/receptor/spatial\n")
print(table(layer2_df$mechanism))

###############################################################################
## LAYER 3 — spatial/proximity test for pairs classified "neither -> spatial"
##  Question: are sender & receiver cells physically closer in Braf_homo than
##  Control? Uses spatial coordinates; compares nearest-neighbour distance.
##  THRESHOLD: Wilcoxon on per-receiver-cell nearest-sender distance, p<WILCOX_P.
###############################################################################
proj <- tryCatch(read.csv("Default Pipeline Analysis-Spatial-Projection.csv"), error=function(e) NULL)
if(!is.null(proj) && requireNamespace("FNN", quietly=TRUE)){
  library(FNN)
  md <- data.frame(Barcode=colnames(obj), celltype=obj$celltype, geno=obj$geno)
  md$x <- proj$X.Coordinate[match(md$Barcode, proj$Barcode)]
  md$y <- proj$Y.Coordinate[match(md$Barcode, proj$Barcode)]
  md <- md[!is.na(md$x),]
  spatial_pairs <- layer2_df %>% filter(grepl("spatial", mechanism)) %>%
    tidyr::separate(pair, into=c("s","t"), sep="->", remove=FALSE)
  l3 <- list()
  for(i in seq_len(nrow(spatial_pairs))){
    s <- spatial_pairs$s[i]; t <- spatial_pairs$t[i]
    dist_geno <- function(gt){
      sc <- md[md$geno==gt & md$celltype==s, c("x","y")]
      tc <- md[md$geno==gt & md$celltype==t, c("x","y")]
      if(nrow(sc)<5 || nrow(tc)<5) return(NULL)
      get.knnx(sc, tc, k=1)$nn.dist[,1]   # each receiver's distance to nearest sender
    }
    dG1 <- dist_geno(G1); dG2 <- dist_geno(G2)
    if(is.null(dG1) || is.null(dG2)) next
    p <- tryCatch(wilcox.test(dG1,dG2)$p.value, error=function(e) NA)
    l3[[i]] <- data.frame(pathway=spatial_pairs$pathway[i], pair=spatial_pairs$pair[i],
                          median_dist_G1=round(median(dG1),1), median_dist_G2=round(median(dG2),1),
                          closer_in=ifelse(median(dG1)<median(dG2), G1, G2), wilcox_p=signif(p,3))
  }
  if(length(l3)>0){
    layer3_df <- bind_rows(l3)
    write.csv(layer3_df, file.path(OUTDIR,"04_spatial_proximity.csv"), row.names=FALSE)
    cat("LAYER 3 done: spatial proximity tested for", nrow(layer3_df), "pairs\n")
  } else cat("LAYER 3: no pairs required spatial testing (or too few cells)\n")
} else {
  cat("LAYER 3 skipped: spatial projection or FNN package unavailable\n")
}

cat("\n=========== PIPELINE COMPLETE — outputs in", OUTDIR, "===========\n")
cat("  00_significant_interactions_all_genotypes.csv\n")
cat("  01_differential_pathways.csv  + _bar.pdf   (the reproducible bar graph)\n")
cat("  02_sender_receiver_pairs.csv               (who talks to whom, classified)\n")
cat("  03_rewiring_vs_volume.csv                  (ligand vs receptor vs spatial)\n")
cat("  04_spatial_proximity.csv                   (proximity test for spatial-class pairs)\n")









library(CellChat); library(dplyr); library(ggplot2)
setwd("C:/Users/write/OneDrive/Documents/MapKSeries_Segmented_MouseDuodenum")

# ---- load the two anchor objects (your existing per-genotype CellChat objects) ----
cc_list <- readRDS("cc_plasma_all.rds")
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
names(cc_list) <- short[names(cc_list)]

# make sure centrality is computed (needed for rankNet), then merge
cc_control <- netAnalysis_computeCentrality(cc_list[["Control"]],   slot.name="netP")
cc_homo    <- netAnalysis_computeCentrality(cc_list[["Braf_homo"]], slot.name="netP")
# object.list ORDER matters: rankNet compares (2nd / 1st) i.e. Braf_homo relative to Control
object.list <- list(Control = cc_control, Braf_homo = cc_homo)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# ---- rankNet: CellChat-native differential pathway ranking with significance ----
# mode="comparison" does the two-condition test; returns per-pathway flow + p-values.
gg <- rankNet(cellchat, mode = "comparison", measure = "weight",
              comparison = c(1, 2),          # Control vs Braf_homo
              stacked = FALSE, do.stat = TRUE, return.data = TRUE)
rn <- gg$signaling.contribution              # data frame: name, group, contribution, pvalues
write.csv(rn, "cellchat_pipeline_out/rankNet_braf_vs_control.csv", row.names = FALSE)

# ---- reshape to one row per pathway: Control vs Braf_homo strength + significance ----
library(tidyr)
wide <- rn %>%
  select(name, group, contribution) %>%
  pivot_wider(names_from = group, values_from = contribution, values_fill = 0)
# p-value is per pathway (same across the two groups in rankNet output)
pvals <- rn %>% group_by(name) %>% summarise(pval = min(pvalues, na.rm = TRUE), .groups="drop")
wide <- left_join(wide, pvals, by = "name")

# log2 fold-change Braf_homo vs Control (pseudocount for zeros)
PSEUDO <- 1e-4
wide <- wide %>%
  mutate(log2FC = log2((Braf_homo + PSEUDO) / (Control + PSEUDO))) %>%
  filter(Braf_homo > 0 | Control > 0)        # INCLUSION: present in >=1 condition

# ---- SIGNIFICANCE FILTER + rank by log2FC, take top 25 ----
P_SIG  <- 0.05
N_SHOW <- 25
sig <- wide %>% filter(pval < P_SIG)          # keep only significantly-different pathways
cat("pathways significant (rankNet p<", P_SIG, "):", nrow(sig), "of", nrow(wide), "\n")

plotdf <- sig %>% arrange(desc(abs(log2FC))) %>% head(N_SHOW) %>% arrange(desc(log2FC))
plotdf$name <- factor(plotdf$name, levels = plotdf$name)

# ---- the bar graph ----
p <- ggplot(plotdf, aes(log2FC, name, fill = log2FC > 0)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE"="#B2182B","FALSE"="#2166AC"),
                    labels = c("TRUE"="up in Braf_homo","FALSE"="up in Control"), name = NULL) +
  labs(x = "log2 FC signaling strength (Braf_homo / Control)", y = NULL,
       title = "Differential signaling pathways: Braf_homo vs Control",
       subtitle = paste0("rankNet (weight), paired Wilcoxon p<", P_SIG,
                         "; top ", N_SHOW, " by |log2FC|; nboot=20")) +
  theme_classic(base_size = 12) +
  theme(axis.text.y = element_text(size = 9))

ggsave("cellchat_pipeline_out/differential_pathways_rankNet_top25.pdf", p, width = 7, height = 7)
cat("saved differential_pathways_rankNet_top25.pdf\n")
print(plotdf[, c("name","Control","Braf_homo","log2FC","pval")], row.names = FALSE)




library(CellChat); library(dplyr)
setwd("C:/Users/write/OneDrive/Documents/MapKSeries_Segmented_MouseDuodenum")

future::plan("sequential")   # single-threaded: avoids the memory-driven worker crashes

cc_list <- readRDS("cc_plasma_all.rds")
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
names(cc_list) <- short[names(cc_list)]
dose_order <- c("Braf_homo","Braf_het","Control","ErkDKO","Braf_homo_ErkDKO","Braf_het_ErkDKO")
cc_list <- cc_list[dose_order]

options(future.globals.maxSize = 30000 * 1024^2)

cc_list_100 <- list()
for(g in dose_order){
  cat("\n==== recomputing", g, "with nboot=100 ====", format(Sys.time()), "\n")
  cc <- cc_list[[g]]
  # Tyler's exact params — ONLY nboot changes (20 -> 100)
  cc <- computeCommunProb(cc, type = "truncatedMean", trim = 0.001,
                          distance.use = TRUE, interaction.range = 100, scale.distance = 0.344,
                          contact.dependent = TRUE, contact.range = 20,
                          nboot = 100)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")
  cc_list_100[[g]] <- cc
  saveRDS(cc_list_100, "cc_plasma_all_nboot100.rds")   # save after EACH genotype
  cat("   done:", length(cc@netP$pathways), "pathways —", g, "saved\n")
}
cat("\n\nALL 6 DONE — cc_plasma_all_nboot100.rds\n")

library(CellChat); library(Seurat); library(patchwork); library(NMF); library(ggalluvial); library(dplyr)
setwd("C:/Users/write/OneDrive/Documents/MapKSeries_Segmented_MouseDuodenum")
dir.create("plots/control", recursive=TRUE, showWarnings=FALSE)
dir.create("plots/comparison", recursive=TRUE, showWarnings=FALSE)

cc_list <- readRDS("cc_plasma_all_nboot100.rds")   # your rebuilt nboot=100 objects
# names should already be the short labels; if not, they are in dose_order used at build
cellchat_control <- cc_list[["Control"]]
cellchat_braf    <- cc_list[["Braf_homo"]]

# centrality is already computed, but re-run to be safe (cheap)
cellchat_control <- netAnalysis_computeCentrality(cellchat_control, slot.name="netP")
cellchat_braf    <- netAnalysis_computeCentrality(cellchat_braf,    slot.name="netP")

groupSize <- as.numeric(table(cellchat_control@idents))
png("plots/control/01a_interaction_count.png", width=700, height=700, res=150)
netVisual_circle(cellchat_control@net$count, vertex.weight=groupSize, weight.scale=TRUE,
                 label.edge=FALSE, title.name="Number of Interactions — Control")
dev.off()
png("plots/control/01b_interaction_strength.png", width=700, height=700, res=150)
netVisual_circle(cellchat_control@net$weight, vertex.weight=groupSize, weight.scale=TRUE,
                 label.edge=FALSE, title.name="Interaction Strength — Control")
dev.off()


set.seed(42)
sig_pathways <- cellchat_control@netP$pathways          # significant pathways in Control
two_random <- sample(sig_pathways, 2)
cat("randomly chosen pathways:", paste(two_random, collapse=", "), "\n")

for(pw in two_random){
  gg <- netAnalysis_contribution(cellchat_control, signaling = pw)
  ggsave(paste0("plots/control/02_LR_contribution_", pw, ".png"), gg, width=6, height=4, dpi=150)
  # visualize the top-contributing L-R pair for that pathway
  pairLR <- extractEnrichedLR(cellchat_control, signaling=pw, geneLR.return=FALSE)
  png(paste0("plots/control/03_top_LR_", pw, ".png"), width=800, height=800, res=150)
  netVisual_individual(cellchat_control, signaling=pw, pairLR.use=pairLR[1,], layout="circle")
  dev.off()
}

png("plots/control/04_bubble_all.png", width=1650, height=1350, res=150)
netVisual_bubble(cellchat_control, remove.isolate=FALSE, angle.x=45)
dev.off()


# use your Seurat object; subset to Control cells to match
obj <- readRDS("obj_annotated.rds")
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
obj$geno <- unname(short[obj$genotype])
seurat_control <- subset(obj, geno=="Control")
Idents(seurat_control) <- "celltype"

for(pw in two_random){
  lr <- subsetCommunication(cellchat_control, signaling=pw)
  genes <- unique(c(lr$ligand, lr$receptor)); genes <- intersect(genes, rownames(seurat_control))
  if(length(genes)>0){
    h <- 750*ceiling(length(genes)/3)
    png(paste0("plots/control/05_genes_", pw, ".png"), width=1500, height=h, res=150)
    print(VlnPlot(seurat_control, features=genes, group.by="celltype", pt.size=0, ncol=3) &
            theme(axis.text.x=element_text(angle=45, hjust=1)))
    dev.off()
  }
}





png("plots/control/06_role_heatmap_outgoing.png", width=1200, height=1000, res=150)
netAnalysis_signalingRole_heatmap(cellchat_control, pattern="outgoing", height=12)
dev.off()
png("plots/control/06_role_heatmap_incoming.png", width=1200, height=1000, res=150)
netAnalysis_signalingRole_heatmap(cellchat_control, pattern="incoming", height=12)
dev.off()
png("plots/control/07_role_scatter.png", width=900, height=750, res=150)
print(netAnalysis_signalingRole_scatter(cellchat_control))
dev.off()


saveRDS(cellchat_control, "cellchat_Control.rds")
saveRDS(cellchat_braf,    "cellchat_Braf_homo.rds")




object.list <- list(Control = cellchat_control, Braf_homo = cellchat_braf)
# sanity check labels match (tutorial's warning)
stopifnot(identical(levels(cellchat_control@idents), levels(cellchat_braf@idents)))
cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list))



gg1 <- compareInteractions(cellchat_merged, show.legend=FALSE, group=c(1,2))
gg2 <- compareInteractions(cellchat_merged, show.legend=FALSE, group=c(1,2), measure="weight")
png("plots/comparison/01_totals.png", width=1200, height=750, res=150)
print(gg1 + gg2)
dev.off()


gg1 <- netVisual_heatmap(cellchat_merged)                    # count; red=up in Braf_homo
gg2 <- netVisual_heatmap(cellchat_merged, measure="weight")  # strength
png("plots/comparison/04_diff_heatmaps.png", width=1800, height=750, res=150)
print(gg1 + gg2)
dev.off()


png("plots/comparison/05_rankNet_stacked.png", width=1050, height=1200, res=150)
print(rankNet(cellchat_merged, mode="comparison", stacked=TRUE,  do.stat=TRUE))
dev.off()
png("plots/comparison/06_rankNet_sidebyside.png", width=1050, height=1200, res=150)
print(rankNet(cellchat_merged, mode="comparison", stacked=FALSE, do.stat=TRUE))
dev.off()


png("plots/comparison/07_outgoing_Control.png", width=1200, height=1000, res=150)
netAnalysis_signalingRole_heatmap(object.list[[1]], pattern="outgoing",
                                  title="Control — Outgoing", height=12, color.heatmap="GnBu")
dev.off()
png("plots/comparison/08_outgoing_Braf_homo.png", width=1200, height=1000, res=150)
netAnalysis_signalingRole_heatmap(object.list[[2]], pattern="outgoing",
                                  title="Braf_homo — Outgoing", height=12, color.heatmap="GnBu")
dev.off()
png("plots/comparison/09_incoming_Control.png", width=1200, height=1000, res=150)
netAnalysis_signalingRole_heatmap(object.list[[1]], pattern="incoming",
                                  title="Control — Incoming", height=12, color.heatmap="OrRd")
dev.off()
png("plots/comparison/10_incoming_Braf_homo.png", width=1200, height=1000, res=150)
netAnalysis_signalingRole_heatmap(object.list[[2]], pattern="incoming",
                                  title="Braf_homo — Incoming", height=12, color.heatmap="OrRd")
dev.off()


library(CellChat); library(Seurat); library(patchwork); library(dplyr); library(ggplot2)
setwd("C:/Users/write/OneDrive/Documents/MapKSeries_Segmented_MouseDuodenum")

# ===================== SETUP: folders + objects =====================
unlink("plots/control",    recursive = TRUE)
unlink("plots/comparison", recursive = TRUE)
dir.create("plots/control",    recursive = TRUE, showWarnings = FALSE)
dir.create("plots/comparison", recursive = TRUE, showWarnings = FALSE)

cc_list <- readRDS("cc_plasma_all_nboot100.rds")
cellchat_control <- netAnalysis_computeCentrality(cc_list[["Control"]],   slot.name="netP")
cellchat_braf    <- netAnalysis_computeCentrality(cc_list[["Braf_homo"]], slot.name="netP")

obj <- readRDS("obj_annotated.rds")
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
obj$geno <- unname(short[obj$genotype])

# ############################################################################
# ##########################   PART 1 — CONTROL   ############################
# ############################################################################

# ---- 1.8.1 overview circle plots (count + strength) ----
groupSize <- as.numeric(table(cellchat_control@idents))
png("plots/control/01a_interaction_count.png", width=900, height=900, res=150)
netVisual_circle(cellchat_control@net$count, vertex.weight=groupSize, weight.scale=TRUE,
                 label.edge=FALSE, title.name="Number of Interactions — Control")
dev.off()
png("plots/control/01b_interaction_strength.png", width=900, height=900, res=150)
netVisual_circle(cellchat_control@net$weight, vertex.weight=groupSize, weight.scale=TRUE,
                 label.edge=FALSE, title.name="Interaction Strength — Control")
dev.off()

# ---- 1.9.3 L-R contribution for 2 random significant pathways + top L-R pair ----
set.seed(42)
two_random <- sample(cellchat_control@netP$pathways, 2)
cat("randomly chosen pathways:", paste(two_random, collapse=", "), "\n")
for(pw in two_random){
  gg <- netAnalysis_contribution(cellchat_control, signaling=pw)
  ggsave(paste0("plots/control/02_LR_contribution_", pw, ".png"), gg, width=7, height=5, dpi=150)
  pairLR <- extractEnrichedLR(cellchat_control, signaling=pw, geneLR.return=FALSE)
  png(paste0("plots/control/03_top_LR_", pw, ".png"), width=900, height=900, res=150)
  netVisual_individual(cellchat_control, signaling=pw, pairLR.use=pairLR[1,], layout="circle")
  dev.off()
}

# ---- 1.9.4 bubble plot: all significant interactions (BIG canvas, small font) ----
png("plots/control/04_bubble_all.png", width=2200, height=4000, res=150)
print(netVisual_bubble(cellchat_control, remove.isolate=FALSE, angle.x=45) +
        theme(axis.text=element_text(size=6)))
dev.off()

# ---- 1.9.5 gene expression violins for the 2 random pathways ----
seurat_control <- subset(obj, geno=="Control"); Idents(seurat_control) <- "celltype"
for(pw in two_random){
  lr <- subsetCommunication(cellchat_control, signaling=pw)
  genes <- intersect(unique(c(lr$ligand, lr$receptor)), rownames(seurat_control))
  if(length(genes)>0){
    h <- 700*ceiling(length(genes)/3)
    png(paste0("plots/control/05_genes_", pw, ".png"), width=1600, height=h, res=150)
    print(VlnPlot(seurat_control, features=genes, group.by="celltype", pt.size=0, ncol=3) &
            theme(axis.text.x=element_text(angle=45, hjust=1, size=8)))
    dev.off()
  }
}

# ---- 1.9 systems-level: role heatmaps (TALL + small font) + scatter ----
png("plots/control/06_role_heatmap_outgoing.png", width=1400, height=3500, res=150)
print(netAnalysis_signalingRole_heatmap(cellchat_control, pattern="outgoing", height=28, font.size=6))
dev.off()
png("plots/control/06_role_heatmap_incoming.png", width=1400, height=3500, res=150)
print(netAnalysis_signalingRole_heatmap(cellchat_control, pattern="incoming", height=28, font.size=6))
dev.off()
png("plots/control/07_role_scatter.png", width=1000, height=850, res=150)
print(netAnalysis_signalingRole_scatter(cellchat_control))
dev.off()

# ---- 1.12 save ----
saveRDS(cellchat_control, "cellchat_Control.rds")
saveRDS(cellchat_braf,    "cellchat_Braf_homo.rds")

# ############################################################################
# ####################   PART 3 — CONTROL vs BRAF_HOMO   #####################
# ####  direction: Control = index 1 (ref), Braf_homo = index 2          #####
# ####  red/positive = up in Braf_homo ; blue/negative = up in Control   #####
# ############################################################################
object.list <- list(Control = cellchat_control, Braf_homo = cellchat_braf)
stopifnot(identical(levels(cellchat_control@idents), levels(cellchat_braf@idents)))
cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list))

# ---- 3.2 compareInteractions: total count + strength ----
gg1 <- compareInteractions(cellchat_merged, show.legend=FALSE, group=c(1,2))
gg2 <- compareInteractions(cellchat_merged, show.legend=FALSE, group=c(1,2), measure="weight")
png("plots/comparison/01_totals.png", width=1200, height=750, res=150)
print(gg1 + gg2)
dev.off()

# ---- 3.4 differential interaction heatmaps: count + strength ----
gg1 <- netVisual_heatmap(cellchat_merged)                    # count; red=up in Braf_homo
gg2 <- netVisual_heatmap(cellchat_merged, measure="weight")  # strength
png("plots/comparison/04_diff_heatmaps.png", width=2000, height=900, res=150)
print(gg1 + gg2)
dev.off()

# ---- 3.5 rankNet: differential pathways (TALL + small font, all pathways) ----
png("plots/comparison/05_rankNet_stacked.png", width=1400, height=4000, res=150)
print(rankNet(cellchat_merged, mode="comparison", stacked=TRUE,  do.stat=TRUE, font.size=6))
dev.off()
png("plots/comparison/06_rankNet_sidebyside.png", width=1400, height=4000, res=150)
print(rankNet(cellchat_merged, mode="comparison", stacked=FALSE, do.stat=TRUE, font.size=6))
dev.off()

# ---- 3.6 outgoing/incoming role heatmaps per condition (TALL + small font) ----
png("plots/comparison/07_outgoing_Control.png", width=1400, height=3500, res=150)
print(netAnalysis_signalingRole_heatmap(object.list[[1]], pattern="outgoing",
                                        title="Control — Outgoing", height=28, font.size=6, color.heatmap="GnBu"))
dev.off()
png("plots/comparison/08_outgoing_Braf_homo.png", width=1400, height=3500, res=150)
print(netAnalysis_signalingRole_heatmap(object.list[[2]], pattern="outgoing",
                                        title="Braf_homo — Outgoing", height=28, font.size=6, color.heatmap="GnBu"))
dev.off()
png("plots/comparison/09_incoming_Control.png", width=1400, height=3500, res=150)
print(netAnalysis_signalingRole_heatmap(object.list[[1]], pattern="incoming",
                                        title="Control — Incoming", height=28, font.size=6, color.heatmap="OrRd"))
dev.off()
png("plots/comparison/10_incoming_Braf_homo.png", width=1400, height=3500, res=150)
print(netAnalysis_signalingRole_heatmap(object.list[[2]], pattern="incoming",
                                        title="Braf_homo — Incoming", height=28, font.size=6, color.heatmap="OrRd"))
dev.off()

# ---- 3.9 rankSimilarity: condition-specific pathways ----
cellchat_merged <- computeNetSimilarityPairwise(cellchat_merged, type="functional")
cellchat_merged <- netEmbedding(cellchat_merged, type="functional", umap.method="uwot")
options(future.rng.onMisuse="ignore")
cellchat_merged <- netClustering(cellchat_merged, type="functional")
png("plots/comparison/16_rankSimilarity.png", width=1200, height=3500, res=150)
print(rankSimilarity(cellchat_merged, type="functional", font.size=6))
dev.off()
saveRDS(cellchat_merged, "cellchat_merged_Control_vs_Braf_homo.rds")

cat("\nDONE — all figures in plots/control and plots/comparison\n")





library(CellChat); library(ComplexHeatmap); library(circlize); library(png); library(grid)
setwd("C:/Users/write/OneDrive/Documents/MapKSeries_Segmented_MouseDuodenum")

unlink(list.files("plots/comparison", pattern="outgoing|incoming", full.names=TRUE))

cc_list <- readRDS("cc_plasma_all_nboot100.rds")
cc_c <- netAnalysis_computeCentrality(cc_list[["Control"]],   slot.name="netP")
cc_b <- netAnalysis_computeCentrality(cc_list[["Braf_homo"]], slot.name="netP")
pw_common <- sort(union(cc_c@netP$pathways, cc_b@netP$pathways))
celltypes <- levels(cc_c@idents)

get_mat <- function(cc, pattern, pws){
  slot <- if(pattern=="outgoing") "outdeg" else "indeg"
  m <- matrix(0, length(pws), length(celltypes), dimnames=list(pws, celltypes))
  for(pw in pws){ ce <- cc@netP$centr[[pw]]
  if(!is.null(ce) && !is.null(ce[[slot]])) m[pw, names(ce[[slot]])] <- ce[[slot]] }
  m
}

aligned_chunked <- function(pattern, cmap, outfile, chunk=10){
  mc <- get_mat(cc_c, pattern, pw_common); mb <- get_mat(cc_b, pattern, pw_common)
  rmax <- pmax(apply(mc,1,max), apply(mb,1,max)); rmax[rmax==0] <- 1
  mc <- mc/rmax; mb <- mb/rmax
  cols <- if(cmap=="GnBu") c("#F7FCF0","#7BCCC4","#0868AC") else c("#FFF7EC","#FC8D59","#B30000")
  cf <- colorRamp2(c(0,0.5,1), cols)
  
  # right-bar = per-pathway total; find the REAL max across both conditions -> set axis to it
  tot_c <- rowSums(mc); tot_b <- rowSums(mb)
  BARMAX <- ceiling(max(c(tot_c,tot_b))*10)/10      # true ceiling, no clipping
  cat(pattern, "- max pathway total strength =", round(max(c(tot_c,tot_b)),2),
      "-> bar axis 0 to", BARMAX, "\n")
  
  chunks <- split(pw_common, ceiling(seq_along(pw_common)/chunk))
  tmp <- character(length(chunks))
  for(i in seq_along(chunks)){
    p <- chunks[[i]]; tmp[i] <- tempfile(fileext=".png")
    # Control strength bar on its LEFT (far-left, visible); Braf on its RIGHT (far-right, visible)
    la_c <- rowAnnotation(strength=anno_barplot(tot_c[p], ylim=c(0,BARMAX),
                                                gp=gpar(fill="grey40"), width=unit(2.2,"cm"),
                                                axis_param=list(at=round(c(0,BARMAX/2,BARMAX),1))))
    ra_b <- rowAnnotation(strength=anno_barplot(tot_b[p], ylim=c(0,BARMAX),
                                                gp=gpar(fill="grey40"), width=unit(2.2,"cm"),
                                                axis_param=list(at=round(c(0,BARMAX/2,BARMAX),1))))
    
    hc <- Heatmap(mc[p,,drop=FALSE], name="importance", col=cf, cluster_rows=FALSE,
                  cluster_columns=FALSE, show_column_names=TRUE, column_title="Control",
                  row_names_side="left", left_annotation=la_c,   # bar on FAR LEFT = visible
                  column_names_rot=45, row_names_gp=gpar(fontsize=10), column_names_gp=gpar(fontsize=10))
    hb <- Heatmap(mb[p,,drop=FALSE], name="importance", col=cf, cluster_rows=FALSE,
                  cluster_columns=FALSE, show_column_names=TRUE, column_title="Braf_homo",
                  show_row_names=FALSE, right_annotation=ra_b,    # bar on FAR RIGHT = visible
                  column_names_rot=45, column_names_gp=gpar(fontsize=10))
    
    png(tmp[i], width=2400, height=length(p)*70+320, res=150)
    draw(hc + hb, ht_gap=unit(6,"mm"),
         column_title=if(i==1) paste0(tools::toTitleCase(pattern),
                                      " — Control vs Braf_homo (outer bars = pathway strength, fixed 0–", BARMAX, ")") else NULL,
         column_title_gp=gpar(fontsize=13, fontface="bold"),
         padding=unit(c(18,4,4,4),"mm"))
    dev.off()
  }
  imgs <- lapply(tmp, readPNG); W <- max(sapply(imgs,function(x)dim(x)[2]))
  gap <- 30; totalH <- sum(sapply(imgs,function(x)dim(x)[1])) + gap*(length(imgs)-1)
  png(outfile, width=W, height=totalH)
  par(mar=c(0,0,0,0)); plot.new(); plot.window(c(0,W),c(0,totalH), xaxs="i", yaxs="i")
  yt <- totalH; for(im in imgs){ih<-dim(im)[1];iw<-dim(im)[2];rasterImage(im,0,yt-ih,iw,yt);yt<-yt-ih-gap}
  dev.off(); unlink(tmp); cat("wrote", outfile, "\n")
}

aligned_chunked("outgoing", "GnBu", "plots/comparison/07_08_outgoing_Control_vs_Braf.png")
aligned_chunked("incoming", "OrRd", "plots/comparison/09_10_incoming_Control_vs_Braf.png")
cat("\nDONE\n")




library(CellChat); library(ggplot2)
setwd("C:/Users/write/OneDrive/Documents/MapKSeries_Segmented_MouseDuodenum")
dir.create("plots/comparison/bubbles", recursive=TRUE, showWarnings=FALSE)

cc_list <- readRDS("cc_plasma_all_nboot100.rds")
cc_c <- netAnalysis_computeCentrality(cc_list[["Control"]],   slot.name="netP")
cc_b <- netAnalysis_computeCentrality(cc_list[["Braf_homo"]], slot.name="netP")

# direction: Control = index 1 (ref), Braf_homo = index 2
# red/positive = up in Braf_homo ; blue = up in Control
object.list <- list(Control = cc_c, Braf_homo = cc_b)
cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list))

celltypes <- levels(cc_c@idents)
cat("cell types:", paste(celltypes, collapse=", "), "\n")


library(CellChat); library(ggplot2)
setwd("C:/Users/write/OneDrive/Documents/MapKSeries_Segmented_MouseDuodenum")
dir.create("plots/comparison/bubbles", recursive=TRUE, showWarnings=FALSE)

cc_list <- readRDS("cc_plasma_all_nboot100.rds")
cc_c <- netAnalysis_computeCentrality(cc_list[["Control"]],   slot.name="netP")
cc_b <- netAnalysis_computeCentrality(cc_list[["Braf_homo"]], slot.name="netP")

# Control = index 1 (ref), Braf_homo = index 2
# red/positive = up in Braf_homo ; blue = up in Control
object.list <- list(Control = cc_c, Braf_homo = cc_b)
cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list))

# the four specific sender -> receiver pairs you want
pairs <- list(
  "Crypt_to_Fibroblast"  = list(src="Crypt",  tgt="Fibroblast"),
  "Crypt_to_Macrophage"  = list(src="Crypt",  tgt="Macrophage"),
  "Villus_to_Plasma"     = list(src="Villus", tgt="Plasma"),
  "Villus_to_Tcell"      = list(src="Villus", tgt="Tcell")
)

for(nm in names(pairs)){
  p <- pairs[[nm]]
  gg <- tryCatch(
    netVisual_bubble(cellchat_merged, sources.use=p$src, targets.use=p$tgt,
                     comparison=c(1,2), angle.x=45, remove.isolate=TRUE,
                     title.name=paste0(gsub("_"," ",nm), " (Control vs Braf_homo)")),
    error=function(e){ cat("  no signal for", nm, ":", conditionMessage(e), "\n"); NULL })
  if(!is.null(gg)){
    n_rows <- length(unique(gg$data$interaction_name))
    ggsave(paste0("plots/comparison/bubbles/", nm, ".png"),
           gg + theme(axis.text=element_text(size=10),
                      axis.text.x=element_text(angle=45, hjust=1)),
           width = 7, height = max(4, n_rows*0.32), dpi=150, limitsize=FALSE)
    cat("wrote", nm, "(", n_rows, "L-R pairs )\n")
  }
}
cat("\nDONE — 4 bubble plots in plots/comparison/bubbles/\n")