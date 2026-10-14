# Instalar si hace falta
#install.packages("ggVennDiagram")

library(ggVennDiagram)
library(ggplot2)

# Extraer identificadores de genes
# Usamos los rownames si contienen los IDs de genes
genes_deseq2 <- rownames(DEGs.DESeq2)
genes_edger  <- rownames(DEGs.edgeR)
genes_voom   <- rownames(DEGs.voom)

# Crear lista de genes por método
gene_lists <- list(
  DESeq2 = unique(genes_deseq2),
  edgeR = unique(genes_edger),
  `limma-voom` = unique(genes_voom)
)

# Diagrama de Venn
ggVennDiagram(gene_lists, label = "count") +
  scale_fill_gradient(low = "white", high = "#4a63e2") +
  labs(
    title = "Genes diferencialmente expresados",
    subtitle = "Comparación DESeq2, edgeR y limma-voom"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 11)
  )

# Crear tabla de correspondencia entre Ensembl ID y symbol
anotacion <- rbind(
  data.frame(gene_id = rownames(DEGs.DESeq2),
             symbol = DEGs.DESeq2$symbol),
  data.frame(gene_id = rownames(DEGs.edgeR),
             symbol = DEGs.edgeR$symbol),
  data.frame(gene_id = rownames(DEGs.voom),
             symbol = DEGs.voom$symbol)
)

anotacion <- anotacion[!duplicated(anotacion$gene_id), ]

# Función para recuperar los símbolos
get_symbols <- function(ids) {
  anotacion$symbol[match(ids, anotacion$gene_id)]
}

# Genes compartidos por los tres métodos
genes_comunes <- Reduce(intersect, gene_lists)
get_symbols(genes_comunes)

# Genes compartidos por pares
get_symbols(intersect(genes_deseq2, genes_edger))

get_symbols(intersect(genes_deseq2, genes_voom))

get_symbols(intersect(genes_edger, genes_voom))
