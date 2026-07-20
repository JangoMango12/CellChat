###############################################################################
##  VILLUS-TIP BRAF SIGNATURE — FULL PIPELINE (consolidated, clean)
##  Steps: 1 DE  ->  2 GSEA  ->  3 ORA  ->  4 allelic-series modules
##         5 broad signature  ->  6 human serrated-lesion validation
##  Run top to bottom. Files on disk (obj_annotated.rds, GSE76987 xlsx) survive
##  an R session reset; only the in-memory objects were lost.
###############################################################################

setwd("C:/Users/write/OneDrive/Documents/MapKSeries_Segmented_MouseDuodenum")

## ---- one-time package install (safe to re-run; skips if present) ----------
if(!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
bioc <- c("fgsea","clusterProfiler","org.Mm.eg.db","org.Hs.eg.db","enrichplot",
          "UCell","GSVA","edgeR","GEOquery")
for(p in bioc) if(!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, update=FALSE, ask=FALSE)
cran <- c("dplyr","ggplot2","msigdbr","clinfun","babelgene","readxl","pROC","tidyr")
for(p in cran) if(!requireNamespace(p, quietly=TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(tidyr)
})

## genotype label map + colors (used throughout) ----------------------------
short <- c("Villin CreERT2(control)"="Control","VillinCreERT2; BrafV600E/+"="Braf_het",
           "Villin CreERT2; BrafV600E/V600E"="Braf_homo","Villin CreERT2; Erk1/Erk2 DKO"="ErkDKO",
           "Villin CreERT2; BrafV600E/+; Erk1/Erk2 DKO"="Braf_het_ErkDKO",
           "Villin CreERT2; BrafV600E/V600E; Erk1/Erk2 DKO"="Braf_homo_ErkDKO")
geno_cols <- c(Braf_homo="#B2182B", Braf_het="#F4A582", Control="grey60",
               ErkDKO="#92C5DE", Braf_homo_ErkDKO="#2166AC", Braf_het_ErkDKO="#4393C3")

###############################################################################
## STEP 1 — DE: Braf_homo tip vs Control tip (Wilcoxon, Control = reference)
###############################################################################
dir.create("tip_DE", showWarnings=FALSE)
obj <- readRDS("obj_annotated.rds")
obj$geno <- unname(short[obj$genotype])

tip <- subset(obj, celltype=="Villus_tip")
tip$geno <- factor(tip$geno); Idents(tip) <- "geno"
cat("tip cells  Braf_homo:", sum(tip$geno=="Braf_homo"),
    " Control:", sum(tip$geno=="Control"), "\n")

de <- FindMarkers(tip, ident.1="Braf_homo", ident.2="Control",
                  test.use="wilcox", logfc.threshold=0, min.pct=0)  # test ALL genes
de$gene <- rownames(de)
write.csv(de, "tip_DE/tip_DE_full.csv", row.names=FALSE)

# filtered lists (log2FC>=0.58, pct>=0.3, padj<0.05) — the chosen thresholds
UP <- de %>% filter(p_val_adj<0.05, avg_log2FC>= 0.58, pct.1>=0.3) %>% arrange(desc(avg_log2FC))
DN <- de %>% filter(p_val_adj<0.05, avg_log2FC<=-0.58, pct.2>=0.3) %>% arrange(avg_log2FC)
write.csv(UP, "tip_DE/tip_UP_loose.csv", row.names=FALSE)
write.csv(DN, "tip_DE/tip_DOWN_loose.csv", row.names=FALSE)
cat("DE done — UP:", nrow(UP), " DOWN:", nrow(DN), "\n")

###############################################################################
## STEP 2 — GSEA (pre-ranked, signed significance; mouse MSigDB)
###############################################################################
suppressPackageStartupMessages({ library(fgsea); library(msigdbr) })
dir.create("tip_GSEA", showWarnings=FALSE)

de <- read.csv("tip_DE/tip_DE_full.csv")
pv <- de$p_val; pv[pv==0] <- .Machine$double.xmin
de$rm <- sign(de$avg_log2FC) * -log10(pv)
fm <- max(abs(de$rm[is.finite(de$rm)])); de$rm[de$rm==Inf]<-fm; de$rm[de$rm==-Inf]<--fm
ranks <- sort(setNames(de$rm, de$gene), decreasing=TRUE)

get_sets <- function(coll, sub=NULL){
  m <- msigdbr(db_species="MM", species="mouse", collection=coll, subcollection=sub)
  split(m$gene_symbol, m$gs_name)
}
collections <- list(
  Hallmark = get_sets("MH"),
  GO_BP    = get_sets("M5","GO:BP"),
  GO_MF    = get_sets("M5","GO:MF"),
  Reactome = get_sets("M2","CP:REACTOME"))

run_gsea <- function(sets, label){
  set.seed(42)
  res <- fgsea(pathways=sets, stats=ranks, minSize=10, maxSize=500, eps=0)
  res <- res[order(res$NES, decreasing=TRUE), ]
  out <- as.data.frame(res); out$leadingEdge <- sapply(out$leadingEdge, paste, collapse=";")
  write.csv(out, file.path("tip_GSEA", paste0("GSEA_", label, ".csv")), row.names=FALSE)
  s <- res[res$padj<0.25, ]
  cat("\n[",label,"] FDR<0.25:", nrow(s), " (FDR<0.05:", sum(res$padj<0.05),
      ") UP:", sum(s$NES>0), " DOWN:", sum(s$NES<0), "\n")
}
for(nm in names(collections)) run_gsea(collections[[nm]], nm)

## GSEA figures: NES barplot + dotplot (FDR<0.05, slide-clean) ---------------
dir.create("tip_GSEA/figures", showWarnings=FALSE)
gf <- list(Hallmark="tip_GSEA/GSEA_Hallmark.csv",GO_BP="tip_GSEA/GSEA_GO_BP.csv",
           GO_MF="tip_GSEA/GSEA_GO_MF.csv",Reactome="tip_GSEA/GSEA_Reactome.csv")
g <- bind_rows(lapply(names(gf), function(n){d<-read.csv(gf[[n]]);d$collection<-n;d}))
sig <- g %>% filter(padj<0.05) %>% arrange(desc(NES))
if(nrow(sig)>0){
  bardf <- sig %>% mutate(dir=ifelse(NES>0,"UP in Braf tip","DOWN in Braf tip"))
  p1 <- ggplot(bardf, aes(reorder(pathway,NES), NES, fill=dir)) + geom_col(width=0.7) + coord_flip() +
    scale_fill_manual(values=c("UP in Braf tip"="#B2182B","DOWN in Braf tip"="#2166AC")) +
    theme_classic(base_size=11) + theme(axis.text.y=element_text(size=8), legend.position="top") +
    labs(title="GSEA significant pathways (FDR<0.05)", x=NULL, y="NES", fill=NULL)
  ggsave("tip_GSEA/figures/NES_barplot_FDR05.png", p1, width=9, height=max(4,nrow(sig)*0.4), dpi=150)
}

###############################################################################
## STEP 3 — ORA (GO over-representation; clusterProfiler; measured-gene bkg)
###############################################################################
suppressPackageStartupMessages({ library(clusterProfiler); library(org.Mm.eg.db) })
dir.create("tip_ORA", showWarnings=FALSE)

de  <- read.csv("tip_DE/tip_DE_full.csv")
UP  <- read.csv("tip_DE/tip_UP_loose.csv")$gene
DN  <- read.csv("tip_DE/tip_DOWN_loose.csv")$gene
sym2ent <- function(s){ m<-mapIds(org.Mm.eg.db, keys=s, column="ENTREZID", keytype="SYMBOL", multiVals="first"); m[!is.na(m)] }
up_e<-sym2ent(UP); dn_e<-sym2ent(DN); univ_e<-sym2ent(de$gene)

run_go <- function(genes, ont, tag){
  ego <- enrichGO(gene=genes, universe=univ_e, OrgDb=org.Mm.eg.db, ont=ont,
                  keyType="ENTREZID", pvalueCutoff=0.05, qvalueCutoff=0.05,
                  minGSSize=10, maxGSSize=500, readable=TRUE)
  if(is.null(ego)||nrow(as.data.frame(ego))==0){ cat("  ",tag,ont,": none\n"); return() }
  d<-as.data.frame(ego); write.csv(d, paste0("tip_ORA/ORA_",tag,"_",ont,".csv"), row.names=FALSE)
  cat("  ",tag,ont,":",nrow(d),"terms\n")
  n<-min(15,nrow(d))
  ggsave(paste0("tip_ORA/dot_",tag,"_",ont,".png"),
         dotplot(ego, showCategory=n)+ggtitle(paste0(tag," GO:",ont)),
         width=9, height=max(4,n*0.4), dpi=150)
}
cat("\nORA UP:\n");   for(o in c("BP","MF","CC")) run_go(up_e,o,"UP")
cat("ORA DOWN:\n"); for(o in c("BP","MF","CC")) run_go(dn_e,o,"DOWN")

###############################################################################
## STEP 4 — allelic-series module scores (UCell): dose-response + ERK reversal
###############################################################################
suppressPackageStartupMessages({ library(UCell); library(AnnotationDbi); library(clinfun) })
dir.create("tip_modules", showWarnings=FALSE)

obj <- readRDS("obj_annotated.rds"); obj$geno <- unname(short[obj$genotype])
tipA <- subset(obj, celltype=="Villus_tip")

go_terms <- c(Energy_respiration="GO:0045333", BrushBorder="GO:0005903", MHCII_antigen="GO:0019886")
detected <- rownames(tipA)
modules <- lapply(go_terms, function(id){
  g <- unique(AnnotationDbi::select(org.Mm.eg.db, keys=id, keytype="GOALL", columns="SYMBOL")$SYMBOL)
  intersect(g[!is.na(g)], detected)
}); names(modules) <- names(go_terms)

tipA <- AddModuleScore_UCell(tipA, features=modules)
dose_order <- c("Braf_homo_ErkDKO","Braf_het_ErkDKO","ErkDKO","Control","Braf_het","Braf_homo")
tipA$geno <- factor(tipA$geno, levels=dose_order)
md <- tipA@meta.data

for(m in names(modules)){
  col <- paste0(m,"_UCell")
  p <- ggplot(md, aes(geno, .data[[col]], fill=geno)) +
    geom_violin(scale="width", trim=TRUE) + geom_boxplot(width=0.12, fill="white", outlier.size=0.2) +
    scale_fill_manual(values=geno_cols) + theme_classic(base_size=11) +
    theme(axis.text.x=element_text(angle=35,hjust=1), legend.position="none") +
    labs(title=paste0(m," module across allelic series"),
         subtitle="MAPK dose low->high; ErkDKO arms = reversal test", x=NULL, y="UCell score")
  ggsave(paste0("tip_modules/module_",m,".png"), p, width=8, height=5, dpi=150)
  
  ba <- md %>% filter(geno %in% c("Control","Braf_het","Braf_homo"))
  ba$dose <- as.integer(factor(ba$geno, levels=c("Control","Braf_het","Braf_homo")))
  jt <- jonckheere.test(ba[[col]], ba$dose, alternative="two.sided", nperm=1000)
  rs <- wilcox.test(md[[col]][md$geno=="Braf_homo"], md[[col]][md$geno=="Braf_homo_ErkDKO"])
  cat("\n[",m,"] dose-response p=",signif(jt$p.value,3),
      " | ERK-rescue p=",signif(rs$p.value,3),"\n")
  print(round(tapply(md[[col]], md$geno, median),4))
}

###############################################################################
## STEP 5 — broad UP signature (union of all sig UP GO terms) + mouse check
###############################################################################
up_files <- list.files("tip_ORA", pattern="ORA_UP_(BP|MF|CC).csv", full.names=TRUE)
up_sig <- unique(unlist(lapply(up_files, function(f) unlist(strsplit(read.csv(f)$geneID,"/")))))
up_sig <- up_sig[!is.na(up_sig) & up_sig!=""]
writeLines(up_sig, "tip_modules/broad_UP_signature_mouse.txt")
cat("\nbroad UP signature:", length(up_sig), "genes\n")

tipA <- AddModuleScore_UCell(tipA, features=list(BroadUP=intersect(up_sig, rownames(tipA))))
cat("mouse tip BroadUP median by genotype:\n")
print(round(tapply(tipA$BroadUP_UCell, tipA$geno, median),4))

###############################################################################
## STEP 6 — human serrated-lesion validation (GSE76987 right colon)
##   ASSOCIATION test only (serrated vs control, serrated vs conventional).
##   NOT a diagnostic/ML biomarker — n too small, no validation cohort.
###############################################################################
suppressPackageStartupMessages({
  library(babelgene); library(readxl); library(org.Hs.eg.db)
  library(GSVA); library(edgeR); library(pROC)
})
dir.create("tip_modules/serrated_validation", showWarnings=FALSE)

## mouse signature -> human orthologs
up_sig <- readLines("tip_modules/broad_UP_signature_mouse.txt")
orth <- orthologs(genes=up_sig, species="mouse", human=FALSE)   # human=FALSE: input is mouse
hcol <- intersect(c("human_symbol","symbol"), colnames(orth))[1]
human_sig <- unique(orth[[hcol]]); human_sig <- human_sig[!is.na(human_sig) & human_sig!=""]
writeLines(human_sig, "tip_modules/broad_UP_signature_human.txt")
cat("\nmouse->human:", length(up_sig), "->", length(human_sig), "genes\n")

## load right-colon counts (serrated SSA/P + control CR + conventional AP live here)
r <- as.data.frame(read_excel("tip_modules/GSE76987/GSE76987_RightColonProcessed.xlsx", sheet=1))
cc <- grep("^Counts ", colnames(r), value=TRUE)
mat <- as.matrix(r[,cc]); rownames(mat) <- r$Ensembl_ID; colnames(mat) <- gsub("^Counts ","",cc)

## normalize (log2-CPM, low-count filter)
mat_f <- mat[rowSums(mat>=10)>=5, ]
cpm <- edgeR::cpm(mat_f, log=TRUE, prior.count=1)

## signature symbols -> Ensembl to match matrix
sig_ens <- mapIds(org.Hs.eg.db, keys=human_sig, column="ENSEMBL", keytype="SYMBOL", multiVals="first")
sig_in <- intersect(unique(sig_ens[!is.na(sig_ens)]), rownames(cpm))
cat("signature genes in dataset:", length(sig_in), "of", length(human_sig), "\n")

## ssGSEA score per sample (new GSVA API)
scores <- gsva(ssgseaParam(cpm, list(BRAF_tip_sig=sig_in)))
sdf <- data.frame(sample=colnames(cpm), score=as.numeric(scores[1,]),
                  group=gsub("-[0-9]+$","",colnames(cpm))) %>%
  filter(group %in% c("SSA/P","CR","AP"))
sdf$group <- factor(sdf$group, levels=c("CR","AP","SSA/P"),
                    labels=c("Control","Conventional adenoma","Serrated (SSA/P)"))
write.csv(sdf, "tip_modules/serrated_validation/scores.csv", row.names=FALSE)

t1 <- wilcox.test(score~group, data=sdf %>% filter(group %in% c("Control","Serrated (SSA/P)")))
t2 <- wilcox.test(score~group, data=sdf %>% filter(group %in% c("Conventional adenoma","Serrated (SSA/P)")))
cat("\nTEST1 Serrated vs Control     p=",signif(t1$p.value,3),"\n")
cat(  "TEST2 Serrated vs Conventional p=",signif(t2$p.value,3),"\n")
print(round(tapply(sdf$score, sdf$group, median),4))

roc_df <- sdf %>% filter(group %in% c("Control","Serrated (SSA/P)"))
ro <- roc(response=roc_df$group, predictor=roc_df$score, levels=c("Control","Serrated (SSA/P)"))
cat("AUROC serrated vs control (this cohort):", round(as.numeric(auc(ro)),3), "\n")

p <- ggplot(sdf, aes(group, score, fill=group)) +
  geom_boxplot(outlier.shape=NA, width=0.6) + geom_jitter(width=0.15, size=1.5, alpha=0.7) +
  scale_fill_manual(values=c("Control"="grey70","Conventional adenoma"="#F4A582","Serrated (SSA/P)"="#B2182B")) +
  theme_classic(base_size=12) + theme(legend.position="none", axis.text.x=element_text(angle=20,hjust=1)) +
  labs(title="Mouse BRAF-tip signature in human colon lesions (GSE76987)",
       subtitle=paste0("Serrated vs Control p=",signif(t1$p.value,2),
                       " | vs Conventional p=",signif(t2$p.value,2)),
       x=NULL, y="ssGSEA signature score")
ggsave("tip_modules/serrated_validation/signature_boxplot.png", p, width=7, height=5, dpi=150)

cat("\n=== PIPELINE COMPLETE ===\n")



library(dplyr)
setwd("C:/Users/write/OneDrive/Documents/MapKSeries_Segmented_MouseDuodenum")

sdf <- read.csv("tip_modules/serrated_validation/scores.csv")
sdf$serrated <- factor(ifelse(sdf$group=="Serrated (SSA/P)","serrated","other"),
                       levels=c("other","serrated"))
cat("class balance:\n"); print(table(sdf$serrated))

# leave-one-out CV predictions
set.seed(1)
loo_pred <- numeric(nrow(sdf))
for(i in 1:nrow(sdf)){
  fit <- glm(serrated ~ score, data=sdf[-i,], family=binomial)
  loo_pred[i] <- predict(fit, newdata=sdf[i,], type="response")
}

# explicit pROC:: to dodge any masking
roc_cv <- pROC::roc(sdf$serrated, loo_pred, levels=c("other","serrated"), quiet=TRUE)
auc_cv <- as.numeric(pROC::auc(roc_cv))
ci_cv  <- pROC::ci.auc(roc_cv)
cat("\nCross-validated AUROC (serrated vs everything else):", round(auc_cv,3), "\n")
cat("95% CI:", round(ci_cv[1],3), "-", round(ci_cv[3],3), "\n")
cat("serrated prevalence (baseline):", round(mean(sdf$serrated=="serrated"),3), "\n")

# ROC curve via ggplot (avoids the masked base plot method entirely)
roc_pts <- data.frame(fpr = 1 - roc_cv$specificities, tpr = roc_cv$sensitivities)
roc_pts <- roc_pts[order(roc_pts$fpr), ]
library(ggplot2)
g <- ggplot(roc_pts, aes(fpr, tpr)) +
  geom_line(color="#B2182B", linewidth=1) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey60") +
  coord_equal() + theme_classic(base_size=12) +
  labs(title=paste0("Serrated vs everything else (LOO-CV)  AUROC=", round(auc_cv,3)),
       x="False positive rate (1 - specificity)", y="True positive rate (sensitivity)")
ggsave("tip_modules/serrated_validation/ROC_serrated_vs_other.png", g, width=6, height=6, dpi=150)
cat("wrote ROC figure\n")