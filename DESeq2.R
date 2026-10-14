#----------------------------------------------------------------------------
#Statistical analysis
#----------------------------------------------------------------------------
library(DESeq2)
#----------------------------------------------------------------------

FDR <- 0.1 #setear

#We make a *DESeqDataSet* (dds) from a count matrix and column data
dds <- DESeqDataSetFromMatrix( countData = count.matrix,
                               colData = sample.table,
                               design = ~ group )

#ANALISIS ESTADISTICO CON DESeq2

# add analysis to dds object
dds <- DESeq(dds)

#results
res.DESeq2 <- results(dds)

#vemos las primeras 5 filas de results
head(res.DESeq2)

#resumen estandar de los resultados que arroja DESeq2
summary(res.DESeq2)

#conteos de resultados
table(res.DESeq2$padj < FDR)

table(res.DESeq2$padj < FDR & res.DESeq2$log2FoldChange< -1) # downregulated
table(res.DESeq2$padj < FDR & res.DESeq2$log2FoldChange> 1) # upregulated

#convertir ensembl a symbol y agregar a results
library(org.Hs.eg.db) #si trabajo en raton usar org.Mm.eg.db

#agrego a results la comlumna symbol
res.DESeq2$symbol <- mapIds(org.Hs.eg.db, #Hs homo sapiens
                     keys = rownames(res.DESeq2), 
                     keytype = "ENSEMBL", 
                     column="SYMBOL")

head(res.DESeq2)

#seleccionamos genes UP y DOWN
UP = res.DESeq2[which(res.DESeq2$padj < 0.1 & res.DESeq2$log2FoldChange > 1),] %>% na.omit()
DOWN = res.DESeq2[which(res.DESeq2$padj < 0.1 & res.DESeq2$log2FoldChange< -1),] %>% na.omit()

#Unimos UP y DOWN en una sola tabla
DEGs.DESeq2 = rbind(UP,DOWN)
DEGs.DESeq2 <- as.data.frame(DEGs.DESeq2)

head(DEGs.DESeq2)

#search genes of interest
#res.DESeq2[which(res.DESeq2$symbol=="ACTA2"),]

#export results
write.csv2(DEGs.DESeq2, "results/DEGs PD vs Control DESeq2.csv")

#-------------------------------------------------------------------------------
#visualizaciones
#-------------------------------------------------------------------------------
#aplicamos transformaciones para mejorar la visualizacion
#covert to regularized logarithm
rld <- rlog(dds)

#convert to variance stabilizing transformation
vst <- vst(dds)

#A PCA plot and a heatmap of the top genes: we need to use the rld object

#Transformación logarítmica regularizada
plotPCA(rld, intgroup="group")

#Transformación estabilizadora de varianza
plotPCA(vst, intgroup="group")

library(EnhancedVolcano)

#hacemos un volcano plot de los resultados
EnhancedVolcano(
  res.DESeq2,
  lab = res.DESeq2$symbol,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = FDR,
  FCcutoff = 1,
  title = "PD vs Control"
)

#heatmap de DEGs
library(pheatmap)

indexDEGs <- which(res.DESeq2$symbol %in% DEGs.DESeq2$symbol)

#heatmap for DEGs
mat.DEseq2 <- assay(rld)[indexDEGs,] #opcion 1

#mat.DEseq2 <- assay(vst)[indexDEGs,] #opcion 2

mat.DEseq2 <- t(scale(t(mat.DEseq2)))  # z score

rownames(mat.DEseq2) <- mapIds(org.Hs.eg.db, 
                               keys = rownames(mat.DEseq2), 
                               keytype = "ENSEMBL", column="SYMBOL")

head(mat.DEseq2)

# Renombrar las muestras para el heatmap
colnames(mat.DEseq2) <- paste(  colnames(mat.DEseq2), "_",
                                sample.table$group)

colnames(mat.DEseq2)

pheatmap(mat.DEseq2, 
         fontsize = 7.5, #tamanio de la fuente 
         cutree_rows = 2) #cortes en las filas

#guardamos el grafico en pdf
pheatmap(mat.DEseq2, 
         fontsize = 7.5, #tamanio de la fuente 
         cutree_rows = 2,
         filename = "plots/heatmap_DESeq2.pdf")

graphics.off()
