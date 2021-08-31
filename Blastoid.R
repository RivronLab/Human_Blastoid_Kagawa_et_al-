library(tidyr)



mapT<-readr::read_tsv("mapTable.tsv")

gs2ensid<-function(fgenes){

   gn <- data.frame(Name=fgenes) %>%
        left_join(mapT, by=c("Name"="gene_name")) %>%
        dplyr::filter(!is.na(gene_id)) %>%
        unite(Name, c("gene_id","Name"), sep="-") %>%
        pull(Name)
    return(gn)
}

## Loading packages and preprocessed data
suppressPackageStartupMessages({
    library(RColorBrewer)
    library(ggplot2)
    library(dplyr)
    library(scran)
    library(Seurat)
    library(kableExtra)
})

options(digits = 4)
options(future.globals.maxSize= 3001289600)
numCores <- 10
library(doParallel)
registerDoParallel(numCores)

rename <- dplyr::rename
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate

## import counts.all,meta.all,mt.gene,ribo.gene
load("preprocessed.GSE177689.rda")


## Filtering 

meta.filter <- meta.all %>%
    filter(nFeature_RNA > 2000 & percent.mt < 12.5)
counts.filter <- counts.all[,rownames(meta.filter)]

sobj <- CreateSeuratObject(counts.filter, min.cells = 5, meta.data = meta.filter)

## Data Analysis

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

## UMAP projection colored according to sample time

idents <- c("Naive_H9","48h","84h","120h","Primed_H9","OKAE_bTS5")
colors <- rev(c("#ff7f00","#4d076a","#ea7bc0","#619CFF","#fdbf6f","#33a02c"))
Idents(sobj) <- "time"
p <- DimPlot(sobj, order = idents, cols = colors, pt.size=0.4)
p[[1]]$layers[[1]]$aes_params$alpha = .7
p <- p + theme(aspect.ratio = 1, legend.position = "bottom", legend.title = element_blank()) +
    guides(color = guide_legend(override.aes = list(size=3)))


## UMAP projection colored according to cluster identity at various resolutions

ResolutionList <- grep("_snn_res", colnames(sobj@meta.data), value = TRUE)
for(Resolution in ResolutionList){
    print(DimPlot(sobj, label = TRUE, group.by = Resolution, order = T, pt.size = 0.4) +
          theme(aspect.ratio=1) +
          theme(legend.title=element_blank())
          )
          
}

## UMAP projection colored according to expression of selected EPI, TE, and PrE markers

markers.list<-list(
    "TE" = c("GATA2","GATA3"),
    "EPI" = c("POU5F1","KLF17"),
    "PRE" = c("GATA4","SOX17")
)
for(set in names(markers.list)){
    ens.gn <- gs2ensid(markers.list[[set]])
    print(FeaturePlot(object = sobj, features = ens.gn, order = T, pt.size = 0.4, max.cutoff = "q99", cols = c("#FEE0D2","#67000D")) + theme(aspect.ratio=1))
}

## Identify top markers based on avg_log2FC at resolution 0.02 and display their expression in a row-scaled heatmap sorted by cluster and sample time.

wanted.resolution <- "RNA_snn_res.0.02" 

## Identify top markers at wanted resolution

Idents(sobj) <- wanted.resolution
top.markers <- FindAllMarkers(sobj, only.pos = TRUE, verbose = FALSE) %>%
    dplyr::filter(p_val_adj < 0.01) %>%
    group_by(cluster) %>%
    top_n(n = 30, wt = avg_log2FC) %>%
    dplyr::filter(cluster %in% c(0,1,3))


## Order time points of interest

times.order <- c("Naive_H9","48h","84h","120h")
sobjT <- subset(sobj, subset = time %in% times.order)
sobjT[["time"]] <- factor(sobjT@meta.data$time, levels = times.order)


## Get CB ordered by cluster id and time

sorted.cells <- subset(sobjT, subset = time %in% times.order)@meta.data %>%
                                                      arrange_at(c(wanted.resolution,"time")) %>%
                                                      pull(cell)

## Get their normalized expr values for the identified top markers

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


# Blastoid Data Integration 


## Loading packages and preprocessed data 

suppressPackageStartupMessages({
    library(dplyr)
    library(scran)
    library(batchelor)
    library(Seurat)
    library(scales)
    library(kableExtra)
})

## import counts.all,meta.all,mt.gene,ribo.gene
load("preprocessed.E-MTAB-3929.GSE109555.Tyser2020.GSE177689.rda")

## Filtering

meta.filter <- meta.all %>%
    filter(! EML %in% c("Hemogenic_Endothelial_Progenitors","Erythroblasts")) %>%
    filter(
        !(pj == "GSE109555") | 
        (pj == "GSE109555" & zhou.sample.1000 == 1)
    ) %>%
    filter(nGene > 2000 & mt.perc < 0.125)


counts.filter <- counts.all[setdiff(rownames(counts.all), mt.gene), meta.filter$cell]

meta.filter <- meta.filter %>%
    mutate(EML=gsub("Primitive_Streak","PriS",EML)) %>%
    mutate(EML=gsub("Emergent_Mesoderm","EmMes",EML)) %>%
    mutate(EML=gsub("YS_Mesoderm","YsMes",EML)) %>%
    mutate(EML=gsub("Nascent_Mesoderm","NasMes",EML)) %>%
    mutate(EML=gsub("Axial_Mesoderm","AxMes",EML)) %>%
    mutate(EML=gsub("Advanced_Mesoderm","AdvMes",EML)) %>%
    mutate(EML=gsub("Non-Neural_Ectoderm","Ectoderm",EML)) %>%
    mutate(EML=gsub("ExE_Mesoderm","ExE_Mes",EML)) %>%
    dplyr::mutate(blastoid_cells = case_when(
                      pj == "GSE177689" & EML == "primed_H9" ~ EML,
                      pj == "GSE177689" & EML == "naive_H9" ~ EML,
                      pj == "GSE177689" ~ "blastoid",
                      TRUE ~ "other"
                  )) %>%
    dplyr::mutate(embryo_cells = case_when(
                      time == "E3" | time == "E4" | time == "E5" ~ "preBl",
                      time == "E6" | time == "E7" | time == "D6" ~ "Bl",
                      pj == "Tyser2020"  ~ "Tyser2020",
                      ! pj == "GSE177689"  ~ "postImpl",
                      TRUE ~ "other"
                  )) %>%
    dplyr::mutate(zhou_time = case_when(
                      pj == "GSE109555" ~ time,
                      TRUE ~ "other"
                  )) %>%
    dplyr::mutate(petro_time = case_when(
                      pj == "E-MTAB-3929" ~ time,
                      TRUE ~ "other"
                  )) %>%
    dplyr::mutate(tyser_tissue = case_when(
                      pj == "Tyser2020" ~ EML,
                      TRUE ~ "other"
                  ))
rownames(counts.filter) <- gsub("_","-",rownames(counts.filter))

## Library normalization 

experiments <- c("E-MTAB-3929","GSE177689","Tyser2020","GSE109555")

expG.set <- list()
for (b in experiments){ 
    temp.cell <- meta.filter %>% filter(pj==b) %>% pull(cell)
    expG.set[[b]] <- rownames(counts.filter )[rowSums(counts.filter[,temp.cell] >=1) >=5]
}
sel.expG <- unlist(expG.set) %>% unique() %>% as.vector()

sce.ob <- list()
for (b in experiments){ 
    temp.M <- meta.filter %>% filter(pj == b) 
    temp.sce <-  SingleCellExperiment(list(counts = as.matrix(counts.filter[sel.expG,temp.M$cell])), colData = temp.M) %>%
        computeSumFactors()
    sce.ob[[b]] <- temp.sce
}

mBN.sce.ob <- multiBatchNorm(sce.ob$`E-MTAB-3929`, sce.ob$GSE177689, sce.ob$Tyser2020, sce.ob$GSE109555)
names(mBN.sce.ob) <- experiments

lognormExp.mBN <- mBN.sce.ob %>%
    lapply(function(x){logcounts(x) %>% as.data.frame()  %>% return()}) %>%
    do.call("bind_cols",.)


## FastMNN 

varGene = 2000
pc = 20

sobj <-
    CreateSeuratObject(
        counts.filter[rownames(lognormExp.mBN),meta.filter$cell], meta.data = meta.filter
    ) %>%
    NormalizeData(verbose = FALSE)

sobj@assays$RNA@data <- as.matrix(lognormExp.mBN[rownames(lognormExp.mBN),colnames(sobj)])

sobj.list <- SplitObject(sobj, split.by = "pj") %>%
    lapply(function(x){ x = FindVariableFeatures(x, verbose=F, nfeatures = varGene)})
sobj.list <- sobj.list[experiments]

set.seed(123)
sobj <- SeuratWrappers::RunFastMNN(sobj.list, verbose = F, features = varGene) %>%
    RunUMAP(reduction = "mnn", dims = 1:pc, verbose = F, n.components = 3L) %>%
    FindNeighbors(reduction = "mnn", dims = 1:pc, verbose = F)

for(plotby in c("blastoid_cells","embryo_cells","zhou_time","petro_time","tyser_tissue")){
    Idents(sobj) <- plotby

    if(plotby == "blastoid_cells"){
        orderidents <- c("naive_H9","primed_H9","blastoid","other")
        cols<-rev(c("#f0590e","#cca95e","#960000","gray90"))
    } else if(plotby == "embryo_cells"){
        orderidents <- c("Bl","postImpl","preBl","CS7","other")
        cols<-rev(c("#5ecc74","#cca95e","#5e69cc","#b81212","gray90"))
    } else if(plotby == "zhou_time"){
        orderidents <- c("D6","D12","D10","D8","other")
        cols<-rev(c("#d73027","#4d076a","#ea7bc0","#619CFF","gray90"))
    } else if(plotby == "petro_time"){ 
        orderidents <- c("E7","E6","E5","E4","E3","other")
        cols<-rev(c("#a50026","#f46d43","#c6e868","#41ab5d","#005824","gray90"))
    }else{
        orderidents <- c(setdiff(levels(Idents(sobj)),"other"),"other")
        cols <- c("gray90", hue_pal()(length(orderidents)-1))
    }
    
    sobj[[plotby]] <- factor(Idents(sobj), levels = orderidents)
    ncol <- max(3, floor(length(orderidents)/6))

    print(
        DimPlot(sobj, order=orderidents, cols=cols, pt.size=0.3) +
        theme(aspect.ratio=1) +
        theme(legend.position="bottom", legend.title=element_blank()) +
        guides(color = guide_legend(override.aes=list(size=3), ncol=ncol))
    )

}
