#----------------------------------------------------------------------------
# Statistical analysis with limma-voom
#----------------------------------------------------------------------------
library(edgeR)
library(limma)
library(dplyr)
library(org.Hs.eg.db)

FDR <- 0.1 #setear

#------------------------------------------------------------
# Design matrix
#------------------------------------------------------------

# Simple model
design <- model.matrix(~ labels)

design

#------------------------------------------------------------
# Create DGEList
#------------------------------------------------------------

y <- DGEList(counts = count.matrix)

# TMM normalization
y <- calcNormFactors(y)

#------------------------------------------------------------
# voom transformation
#------------------------------------------------------------

#v <- voom(y, design,  plot = TRUE)

v <- voomWithQualityWeights(y, design,  plot = TRUE)  #otra opcion

#------------------------------------------------------------
# Linear model
#------------------------------------------------------------

fit <- lmFit( v,  design )

fit <- eBayes(fit)

# PD vs Control
res.voom <- topTable(
  fit,
  coef = paste0("labels", levels(labels)[2]),
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

head(res.voom)

#------------------------------------------------------------
# Significant genes
#------------------------------------------------------------

table(res.voom$adj.P.Val < FDR)

table(res.voom$adj.P.Val < FDR & res.voom$logFC < -1 )

table( res.voom$adj.P.Val < FDR & res.voom$logFC > 1 )

#------------------------------------------------------------
# Add symbols
#------------------------------------------------------------

res.voom$symbol <- mapIds(
  org.Hs.eg.db,
  keys = rownames(res.voom),
  keytype = "ENSEMBL",
  column = "SYMBOL"
)

UP.voom <- res.voom[which(res.voom$adj.P.Val < FDR & res.voom$logFC > 1),] %>% na.omit()

DOWN.voom <- res.voom[which(res.voom$adj.P.Val < FDR & res.voom$logFC < -1),] %>% na.omit()

DEGs.voom <- rbind( UP.voom, DOWN.voom )

# Export
write.csv2( DEGs.voom, "results/DEGs_limma_voom_PD_vs_Control.csv")

# Gene of interest
#res.voom[  which(res.voom$symbol == "ACTA2"),]

#----------------------------------------------------------------------------
# Visualizations limma-voom
#----------------------------------------------------------------------------

library(EnhancedVolcano)
library(pheatmap)

# logCPM normalized expression
logCPM <- cpm( y,  log = TRUE,  prior.count = 1)

# PCA using voom expression matrix
pca <- prcomp( t(v$E) )

plot(
  pca$x[,1],
  pca$x[,2],
  col = as.numeric(labels),
  pch = 19,
  xlab = "PC1",
  ylab = "PC2"
)

legend(
  "topright",
  legend = levels(labels),
  col = 1:length(levels(labels)),
  pch = 19
)

# Volcano
EnhancedVolcano(
  res.voom,
  lab = res.voom$symbol,
  x = "logFC",
  y = "adj.P.Val",
  pCutoff = FDR,
  FCcutoff = 1,
  title = "PD vs Control - limma-voom"
)

# Heatmap of DEGs
indexDEGs <- which( rownames(logCPM) %in% rownames(DEGs.voom) )

mat.voom <- logCPM[indexDEGs, ]

mat.voom <- t(scale(t(mat.voom))) # z score

rownames(mat.voom) <- mat.symbol <- mapIds( org.Hs.eg.db,
                                            keys = rownames(mat.voom),
                                            keytype = "ENSEMBL",
                                            column = "SYMBOL")

colnames(mat.voom) <- paste0( colnames(mat.voom), "_", labels )

pheatmap(  mat.voom,
           fontsize = 7.5,
           cutree_rows = 2)

#guardamos el grafico en pdf
pheatmap(mat.voom, 
         fontsize = 7.5, #tamanio de la fuente 
         cutree_rows = 2,
         filename = "plots/heatmap_edgeR.pdf")

graphics.off()
