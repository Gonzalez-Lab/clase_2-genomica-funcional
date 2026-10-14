#----------------------------------------------------------------------------
# Statistical analysis with edgeR
#----------------------------------------------------------------------------

library(edgeR)
library(dplyr)
library(org.Hs.eg.db)

FDR <- 0.1 #setear

#------------------------------------------------------------
# Design matrix
#------------------------------------------------------------

# Simple model: PD vs Control
design <- model.matrix(~ labels)

design

# Create DGEList object
y <- DGEList(counts = count.matrix)

# TMM normalization
y <- calcNormFactors(y)

#------------------------------------------------------------
# Estimate dispersion and fit model
#------------------------------------------------------------
y <- estimateDisp(y, design)

fit <- glmQLFit(y, design)

# groupPD coefficient = PD vs Control
qlf <- glmQLFTest(fit, 
                  coef = paste0("labels", levels(labels)[2]) )

# Results
res.edgeR <- topTags(qlf, n = Inf, sort.by = "PValue")$table

head(res.edgeR)

# Add FDR significance summaries
table(res.edgeR$FDR < FDR)

table(res.edgeR$FDR < FDR & res.edgeR$logFC < -1)

table(res.edgeR$FDR < FDR &  res.edgeR$logFC > 1)

#------------------------------------------------------------
# Add gene symbols
#------------------------------------------------------------

res.edgeR$symbol <- mapIds(
  org.Hs.eg.db,
  keys = rownames(res.edgeR),
  keytype = "ENSEMBL",
  column = "SYMBOL"
)

# DEGs
UP.edgeR <- res.edgeR[which(res.edgeR$FDR < FDR & res.edgeR$logFC > 1),] %>% na.omit()

DOWN.edgeR <- res.edgeR[which(res.edgeR$FDR < FDR & res.edgeR$logFC < -1),] %>% na.omit()

DEGs.edgeR <- rbind(UP.edgeR, DOWN.edgeR)

# Search gene of interest
#res.edgeR[which(res.edgeR$symbol == "ACTA2"), ]

# Export
write.csv2( DEGs.edgeR, "results/DEGs_edgeR_PD_vs_Control.csv")

head(DEGs.edgeR)

#----------------------------------------------------------------------------
# Visualizations edgeR
#----------------------------------------------------------------------------

library(EnhancedVolcano)
library(pheatmap)

# logCPM normalized expression
logCPM <- cpm( y , log = TRUE , prior.count = 1)

# PCA
pca <- prcomp(t(logCPM))

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
  res.edgeR,
  lab = res.edgeR$symbol,
  x = "logFC",
  y = "FDR",
  pCutoff = FDR,
  FCcutoff = 1,
  title = "PD vs Control - edgeR"
)

# Heatmap of DEGs
indexDEGs <- which( rownames(logCPM) %in% rownames(DEGs.edgeR) )

mat.edgeR <- logCPM[indexDEGs, ]

mat.edgeR <- t(scale(t(mat.edgeR))) # z score

rownames(mat.edgeR) <- mapIds(org.Hs.eg.db,
                       keys = rownames(mat.edgeR),
                       keytype = "ENSEMBL",
                       column = "SYMBOL")

colnames(mat.edgeR) <- paste0( colnames(mat.edgeR), "_",
                               sample.table$group )

pheatmap(
  mat.edgeR,
  fontsize = 7.5,
  cutree_rows = 2
)

#guardamos el grafico en pdf
pheatmap(mat.edgeR, 
         fontsize = 7.5, #tamanio de la fuente 
         cutree_rows = 2,
         filename = "plots/heatmap_edgeR.pdf")

graphics.off()
