###############################################################################
# -----------------------------------------------------------------------------
# Enriquecimiento funcional
# -----------------------------------------------------------------------------

# Evaluamos si los DEGs se encuentran enriquecidos en procesos biologicos especificos.
#
# Se utiliza Gene Ontology (GO) - Biological Process (BP)
#                                 Cellular component (CC)
#                                 Molecular function (MF)

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ggplot2)
library(dplyr)

# -----------------------------------------------------------------------------
# 1 Definir universo de genes y DEGs
# -----------------------------------------------------------------------------

# Como universo utilizamos todos los genes evaluados (estan en ensembl)
all_genes <- rownames(count.matrix)

#DEGs a evaluar
DEGs <- rownames(DEGs.DESeq2) #cargar aca los DEGs
#DEGs <- rownames(DEGs.edgeR)
#DEGs <- rownames(DEGs.voom)

# -----------------------------------------------------------------------------
# 2 Gene Ontology con clusterProfiler
# -----------------------------------------------------------------------------

go_bp <- clusterProfiler::enrichGO(
  gene          = DEGs,
  universe      = all_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENSEMBL",         #ENTREZID
  ont           = "BP",             # BP, CC, MF
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 1,
  readable      = TRUE
)

# Convertir resultados a tabla
go_bp_df <- as.data.frame(go_bp)

# Guardar resultados
write.csv(go_bp_df, "results/ORA_DEGs_GO_BP.voom.csv")

# -----------------------------------------------------------------------------
# 3 Visualizacion del enriquecimiento funcional
# -----------------------------------------------------------------------------
# El dotplot muestra los principales procesos biologicos enriquecidos.
# Eje X:GeneRatio = proporcion de genes de nuestra lista asociados al termino.
# Tamaño del punto: numero de genes asociados al termino.
# Color: significancia estadistica ajustada (p.adjust).

top_terms <- 10 #numero maximo de terminos que quiero ver en el grafico

p_go <- enrichplot::dotplot(
    go_bp,
    showCategory = top_terms, #cuantos terminos quiero mostrar
    x = "GeneRatio") +
    labs(
      title = "GO enrichment of DEGs",
      subtitle = "RNA-seq subtantia nigra control vs PD",
      x = "Gene ratio",
      y = NULL
    ) +
    theme_classic(base_size = 13) +
    theme( plot.title = element_text(face = "bold"),
           axis.text.y = element_text(size = 10)  )
  
  
# Mostrar grafico
p_go

ggsave(
  "plots/ORA_GO_BP_dotplot.png",
  plot = p_go,
  width = 10,
  height = 7,
  dpi = 300
)

#################################################################################
# Gene set enrichment analysis GSEA 
################################################################################

# Resultados COMPLETOS, no solo los DEGs
# para DESeq2
geneList <- res.DESeq2$stat #generamos un vector con los valores del estadistico
names(geneList) <- rownames(res.DESeq2) #nombramos el vector con los geneid

#para edgeR
#geneList <- res.edgeR$F
#names(geneList) <- rownames(res.edgeR)

#para limma voom
geneList <- res.voom$t 
names(geneList) <- rownames(res.voom)

#limpiamos los valores
geneList <- geneList[!is.na(geneList) & is.finite(geneList)]

# Ordenar de mayor a menor
geneList <- sort(geneList, decreasing = TRUE)

# Verificar
head(geneList)
length(geneList)

gsea_GO <- gseGO(
  geneList = geneList,
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  seed = TRUE,
  verbose = FALSE
)

# Ver resultados
gsea_results <- as.data.frame(gsea_GO)

head(gsea_results[,1:7])

# Número de procesos enriquecidos
nrow(gsea_results)

write.csv(gsea_results,"results/gsea bp functional analysis.csv")

#VISUALIZACION
#defino numero maximo de terminos que quiero graficar
top_terms <- 10

#dot plot
gsea_plot <- dotplot( gsea_GO,
                      showCategory = top_terms) +
                      ggtitle("GSEA - Parkinson vs controles") +
                      theme_bw()

gsea_plot

ggsave( "plots/GSEA_GO_BP_dotplot.png",
          plot = gsea_plot,
          width = 10,
          height = 7,
          dpi = 300)

#visualizacion UP vs DOWN
gsea_plot_up.down <- dotplot(gsea_GO,
                             showCategory = top_terms,
                             split = ".sign") +
                             facet_grid(. ~ .sign) +
                             ggtitle("GSEA - Parkinson vs controles") +
                             theme_bw()

gsea_plot_up.down

ggsave("plots/GSEA.up.down_GO_BP_dotplot.png",
        plot = gsea_plot_up.down,
        width = 10,
        height = 7,
        dpi = 300)

# Graficar el proceso más significativo
gseaplot2(
    gsea_GO,
    geneSetID = gsea_results$ID[1],
    title = gsea_results$Description[1])
