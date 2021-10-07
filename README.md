
-   [****1**
    Blastoid](https://data.bioinfo.vbc.ac.at/rivron.grp/blastoid/blastoid.html)
    -   [**Read
        preprocessing](https://data.bioinfo.vbc.ac.at/rivron.grp/blastoid/blastoid.html#read-preprocessing)
    -   [**Loading packages and preprocessed
        data](https://data.bioinfo.vbc.ac.at/rivron.grp/blastoid/blastoid.html#loading-packages-and-preprocessed-data "Loading packages and preprocessed data")
    -   [**Filtering](https://data.bioinfo.vbc.ac.at/rivron.grp/blastoid/blastoid.html#filtering)
    -   [**Data
        Analysis](https://data.bioinfo.vbc.ac.at/rivron.grp/blastoid/blastoid.html#data-analysis)
    -   [**Data
        Visualization](https://data.bioinfo.vbc.ac.at/rivron.grp/blastoid/blastoid.html#data-visualization)
        -   [**UMAP/Sampletime](https://data.bioinfo.vbc.ac.at/rivron.grp/blastoid/blastoid.html#umapsampletime)
        -   [**UMAP/Clusters](https://data.bioinfo.vbc.ac.at/rivron.grp/blastoid/blastoid.html#umapclusters)
        -   [**UMAP/Marker](https://data.bioinfo.vbc.ac.at/rivron.grp/blastoid/blastoid.html#umapmarker)
        -   [**Marker
            Heatmap](https://data.bioinfo.vbc.ac.at/rivron.grp/blastoid/blastoid.html#marker-heatmap)
-   [****2** Blastoid Data
    Integration](https://data.bioinfo.vbc.ac.at/rivron.grp/blastoid/blastoid-data-integration.html)
    -   [**Read
        preprocessing](https://data.bioinfo.vbc.ac.at/rivron.grp/blastoid/blastoid-data-integration.html#read-preprocessing-1)
    -   [**Loading packages and preprocessed
        data](https://data.bioinfo.vbc.ac.at/rivron.grp/blastoid/blastoid-data-integration.html#loading-packages-and-preprocessed-data-1 "Loading packages and preprocessed data")
    -   [**Filtering](https://data.bioinfo.vbc.ac.at/rivron.grp/blastoid/blastoid-data-integration.html#filtering-1)
    -   [**Library
        normalization](https://data.bioinfo.vbc.ac.at/rivron.grp/blastoid/blastoid-data-integration.html#library-normalization)
    -   [**FastMNN](https://data.bioinfo.vbc.ac.at/rivron.grp/blastoid/blastoid-data-integration.html#fastmnn)

[**](https://data.bioinfo.vbc.ac.at/rivron.grp/blastoid/# "Toggle Sidebar")[**](https://data.bioinfo.vbc.ac.at/rivron.grp/blastoid/# "Search")

[**](https://data.bioinfo.vbc.ac.at/rivron.grp/blastoid/# "Share")

Facebook

Twitter

LinkedIn

Weibo

Instapaper

[**](https://data.bioinfo.vbc.ac.at/rivron.grp/blastoid/# "Font Settings")

A

A

Serif

Sans

White

Sepia

Night

[**](https://data.bioinfo.vbc.ac.at/rivron.grp/blastoid/# "Information about the toolbar")

**[Blastoid scRNAseq processing](https://data.bioinfo.vbc.ac.at/rivron.grp/blastoid/)
=====================================================================================

1 Blastoid[](https://data.bioinfo.vbc.ac.at/rivron.grp/blastoid/#blastoid)
==========================================================================

Read preprocessing[](https://data.bioinfo.vbc.ac.at/rivron.grp/blastoid/#read-preprocessing) {.hasAnchor}
--------------------------------------------------------------------------------------------

Raw reads can be obtained from
[GSE177689](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE177689).

Smart-Seq transcriptome sequencing experiments were analyzed using
genome sequence and gene annotation from Ensembl GRCh38 release 103 as
reference. For gene expression quantification RNA-seq reads were first
trimmed using trim-galore v0.6.6 and thereafter aligned to the genome
using hisat2 v2.2.1. Uniquely mapping reads in genes were quantified
using htseq-count v0.13.5 with parameter -s no. TPM estimates were
obtained using RSEM v1.3.3 with parameter –single-cell-prior

Loading packages and preprocessed data[](https://data.bioinfo.vbc.ac.at/rivron.grp/blastoid/#loading-packages-and-preprocessed-data) {.hasAnchor}
------------------------------------------------------------------------------------------------------------------------------------

-   Raw count data and metadata files can be obtained from NCBI GEO
    [GSE177616](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE177616)
    or from here:
    -   [blastoid.H9.okae\_bts5\_\_scRNAseq.unfiltered\_\_counts.RDS](https://data.bioinfo.vbc.ac.at/rivron.grp/blastoid/blastoid.H9.okae_bts5__scRNAseq.unfiltered__counts.RDS)
    -   [blastoid.H9.okae\_bts5\_\_scRNAseq.unfiltered\_\_metadata.tsv](https://data.bioinfo.vbc.ac.at/rivron.grp/blastoid/blastoid.H9.okae_bts5__scRNAseq.unfiltered__metadata.tsv)

**

``` {.sourceCode .r}
suppressPackageStartupMessages({
    library(RColorBrewer)
    library(ggplot2)
    library(dplyr)
    library(scran)
    library(Seurat)
    library(kableExtra)
})

counts.all <- readRDS("blastoid.H9.okae_bts5__scRNAseq.unfiltered__counts.RDS") %>% tibble::column_to_rownames("gene_id")
meta.all <- readr::read_tsv("blastoid.H9.okae_bts5__scRNAseq.unfiltered__metadata.tsv") %>%
        dplyr::mutate(time = case_when(
                      grepl("24h",samplename) ~ "24h",
                      grepl("60h",samplename) ~ "60h",
                      grepl("96h",samplename) ~ "96h",
                      TRUE ~ gsub("-.*","",samplename)
                  )) %>% data.frame()
```

**

``` {.sourceCode .r}
rownames(meta.all) <- meta.all$sampleid

mtname="^ENSG[[:digit:]]+-MT-"
```

Filtering[](https://data.bioinfo.vbc.ac.at/rivron.grp/blastoid/#filtering) {.hasAnchor}
--------------------------------------------------------------------------

Based on initial evaluation of per-cell quality control metrices and
outlier identification using the median absolute deviation algorithm,
cells with \<= 2000 detected genes or \>= 12.5% mitochondrial gene
percentage were filtered out. Only genes detected in at least 5 cells
were retained.

**

``` {.sourceCode .r}
meta.filter <- subset(meta.all,
    nFeature_RNA > 2000 & percent.mt < 12.5)
counts.filter <- counts.all[,rownames(meta.filter)]

sobj <- CreateSeuratObject(counts.filter, min.cells = 5, meta.data = meta.filter)
```

Data Analysis[](https://data.bioinfo.vbc.ac.at/rivron.grp/blastoid/#data-analysis) {.hasAnchor}
----------------------------------------------------------------------------------

Count-data were log-normalized, top 3000 highly variable were selected,
and standardization of per gene expression values across cells was
performed using NormalizeData, FindVariableFeatures and ScaleData data
functions in Seurat. Principal component analysis (PCA) based on the
standardized highly variable features was used for linear dimension
reduction, a shared nearest neighbor (SNN) graph was constructed on the
dimensionally reduced data, and the graph was partitioned using a SNN
modularity optimization based clustering algorithm at a range of
resolutions using RunPCA, FindNeighbors and FindClusters from Seurat
with default settings. Cluster marker genes were identified with the
Wilcox likelihood-ratio test using the FindAllMarkers function. Uniform
Manifold Approximation and Projection (UMAP) was used for visualization.

**

``` {.sourceCode .r}
set.seed(1234)
varGene = 3000
pc = 20
resolutionwanted=c(0.02,1)
sobj <- sobj %>% Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(nfeatures = varGene) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    FindNeighbors(verbose = FALSE) %>%
    FindClusters(resolution = resolutionwanted, verbose=FALSE) %>%
    RunUMAP(dims = 1:pc, n.components = 3L, min.dist = 0.5, verbose=FALSE)
```

    ## Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    ## To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    ## This message will be shown once per session

Data Visualization[](https://data.bioinfo.vbc.ac.at/rivron.grp/blastoid/#data-visualization) {.hasAnchor}
--------------------------------------------------------------------------------------------

### UMAP/Sampletime[](https://data.bioinfo.vbc.ac.at/rivron.grp/blastoid/#umapsampletime) {.hasAnchor}

UMAP projection colored according to sample time

**

``` {.sourceCode .r}
idents <- c("naive_H9","24h","60h","96h","primed_H9","okae_bts5")
colors <- rev(c("#ff7f00","#4d076a","#ea7bc0","#619CFF","#fdbf6f","#33a02c"))
Idents(sobj) <- "time"
p <- DimPlot(sobj, order = idents, cols = colors, pt.size=0.4)
p[[1]]$layers[[1]]$aes_params$alpha = .7
p <- p + theme(aspect.ratio = 1, legend.position = "bottom", legend.title = element_blank()) +
    guides(color = guide_legend(override.aes = list(size=3)))
print(p)
```

![](./Blastoid%20scRNAseq%20processing_files/blastoid-data-analysis-umap-sampletime-1.png)

### UMAP/Clusters[](https://data.bioinfo.vbc.ac.at/rivron.grp/blastoid/#umapclusters) {.hasAnchor}

UMAP projection colored according to cluster identity at various
resolutions

**

``` {.sourceCode .r}
ResolutionList <- grep("_snn_res", colnames(sobj@meta.data), value = TRUE)
for(Resolution in ResolutionList){
    print(DimPlot(sobj, label = TRUE, group.by = Resolution, order = T, pt.size = 0.4) +
          theme(aspect.ratio=1) +
          theme(legend.title=element_blank())
          )
          
}
```

![](./Blastoid%20scRNAseq%20processing_files/blastoid-data-analysis-umap-cluster-1.png)![](./Blastoid%20scRNAseq%20processing_files/blastoid-data-analysis-umap-cluster-2.png)

### UMAP/Marker[](https://data.bioinfo.vbc.ac.at/rivron.grp/blastoid/#umapmarker) {.hasAnchor}

UMAP projection colored according to expression of selected EPI, TE, and
PrE markers

**

``` {.sourceCode .r}
markers.list<-list(
    "TE" = c("GATA2","GATA3"),
    "EPI" = c("POU5F1","KLF17"),
    "PRE" = c("GATA4","SOX17")
)
for(set in names(markers.list)){
    ens.gn <- gs2ensid(markers.list[[set]])
    print(FeaturePlot(object = sobj, features = ens.gn, order = T, pt.size = 0.4, max.cutoff = "q99",cols = c("#FEE0D2","#67000D")) + theme(aspect.ratio=1))
}
```

![](./Blastoid%20scRNAseq%20processing_files/blastoid-data-analysis-umap-markers-1.png)![](./Blastoid%20scRNAseq%20processing_files/blastoid-data-analysis-umap-markers-2.png)![](./Blastoid%20scRNAseq%20processing_files/blastoid-data-analysis-umap-markers-3.png)

### Marker Heatmap[](https://data.bioinfo.vbc.ac.at/rivron.grp/blastoid/#marker-heatmap) {.hasAnchor}

Identify top markers based on avg\_log2FC at resolution 0.02 and display
their expression in a row-scaled heatmap sorted by cluster and sample
time.

**

``` {.sourceCode .r}
#wanted.resolution <- "RNA_snn_res.0.02" 
wanted.resolution <- "blastoid.Fig2b.lowres"
## Identify top markers at wanted resolution

Idents(sobj) <- wanted.resolution
top.markers <- FindAllMarkers(sobj, only.pos = TRUE, verbose = FALSE) %>%
    dplyr::filter(p_val_adj < 0.01) %>%
    group_by(cluster) %>%
    top_n(n = 30, wt = avg_log2FC) %>%
    dplyr::filter(cluster %in% c(0,1,3))


## Order time points of interest

times.order <- c("naive_H9","24h","60h","96h")
sobjT <- subset(sobj, subset = time %in% times.order)
sobjT[["time"]] <- factor(sobjT@meta.data$time, levels = times.order)


## Get CB ordered by cluster id and time

sorted.cells <- subset(sobjT, subset = time %in% times.order)@meta.data %>%
                                                      arrange_at(c(wanted.resolution,"time")) %>%
                                                      pull(sampleid) %>% as.character()

# Get their normalized expr values for the identified top markers

lognormExp <- as.matrix(GetAssayData(sobjT, slot = "data"))[ top.markers$gene, sorted.cells]
lognormExp.scaled <- t(apply(lognormExp,1,scale)) %>% data.frame() %>%
    mutate(across(everything(), ~ ifelse(. > 2.5, 2.5, .))) %>%
    mutate(across(everything(), ~ ifelse(. < -2.5, -2.5, .))) %>%
    setNames(colnames(lognormExp))

colanno <- sobjT@meta.data[match(colnames(lognormExp.scaled),rownames(sobjT@meta.data)),c(wanted.resolution,"time")]

heat.col <- colorRampPalette(c("#0D0887FF","#0D0887FF","#0D0887FF","#0D0887FF","#7E03A8FF","#CC4678FF","#F89441FF","#F0F921FF","#F0F921FF"))(100)

pheatmap::pheatmap(lognormExp.scaled,
                   width = 500,
                   height = 300,
                   annotation = colanno,
                   cluster_rows = F,
                   cluster_cols = F,
                   scale = "none",
                   show_colnames = F,
                   show_rownames = T,
                   color = heat.col,
                   border_color = "NA",
                   fontsize = 4,
                   fontsize_row = 2)
```

![](./Blastoid%20scRNAseq%20processing_files/blastoid-data-analysis-heatmap-markers-1.png)

[**](https://data.bioinfo.vbc.ac.at/rivron.grp/blastoid/blastoid-data-integration.html)
