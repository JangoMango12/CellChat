# CellChat
CellChat Analysis Code for MapK Allelic Series Mouse Duodenum Day 2 from Spatial Transcriptomic VisiumHD.
# MAPK Allelic Series — Spatial Transcriptomics of Mouse Duodenum

Visium HD spatial transcriptomics analysis of a 6-genotype BRAF/ERK (MAPK)
allelic series in segmented mouse duodenum (day 2 post-tamoxifen, pre-neoplastic).

## Genotypes
Control, Braf V600E/+, Braf V600E/V600E, Erk1/2 DKO,
Braf V600E/+ ; Erk DKO, Braf V600E/V600E ; Erk DKO

## Analysis overview

**1. Cell-type annotation & epithelial zonation**
Epithelial cells assigned to four crypt–villus zones (Crypt, Junction, Villus,
Villus_tip) via k-means on PCA + label transfer; stromal/immune populations
annotated from marker panels.

**2. Villus-tip differential expression & pathway analysis**
- DE: Braf_homo tip vs Control tip (Wilcoxon).
- GO over-representation (ORA, clusterProfiler) on the DE gene list, measured-gene
  background. Primary method given the single-section design.
- GSEA (fgsea, mouse MSigDB) as a complementary pre-ranked analysis.
- Allelic-series module scores (UCell) testing dose-dependence and ERK-reversal.

**3. Human validation**
Mouse tip signature mapped to human orthologs and scored (ssGSEA) in a public
serrated-lesion cohort (GSE76987). Association test only, not a diagnostic biomarker.

**4. Cell–cell communication (CellChat, spatial, nboot=100)**
Per-genotype spatial CellChat objects. Differential pathway analysis
(Braf_homo vs Control and vs Braf_homo_ErkDKO), per-pathway composites, and a
tip→immune antigen-presentation analysis (per-cell normalized) plus a
permutation-based spatial neighborhood-enrichment test.

## Repository structure
- `R/` — analysis scripts (see order below)
- `README.md` — this file
- Data files are not tracked (see `.gitignore`); available on request.

## Script order
1. Cell-type annotation & zonation
2. Villus-tip DE → ORA/GSEA → allelic-series modules
3. Human serrated-lesion validation
4. CellChat build (nboot=100) → comparison figures → tip-immune analysis

## Notes / limitations
- One tissue section per genotype (n=1); the allelic-series dose-response is the
  primary rigor given this design.
- CellChat infers communication from co-expression and spatial proximity
  (hypothesis-generating, not direct measurement).

_Rutgers Cancer Institute of New Jersey._
