# Clase 2 | RNA-seq: de los counts a la interpretación biológica

**Curso:** Descifrando el código epigenético: de la genómica funcional al análisis bioinformático de datos ómicos  
**Modalidad:** práctico en R / RStudio

## Objetivo

Partiendo de una matriz de *counts* obtenida tras el procesamiento de RNA-seq para identificar genes diferencialmente expresados (DEGs), comparar tres métodos estadísticos y explorar su significado biológico mediante enriquecimiento funcional.

**Caso de estudio:** RNA-seq de sustancia negra humana, con muestras de enfermedad de Parkinson (PD) y controles (dataset **GSE168496**, trabajado en la clase anterior). El objetivo es comparar **PD vs. Control**. Los resultados del ejercicio no deben interpretarse como una validación clínica.

## Material del repositorio

| Script | Contenido |
| --- | --- |
| `RNAseq data preparation.R` | Importación de counts y metadata, correspondencia entre muestras, definición de grupos y filtrado de genes de baja expresión con `filterByExpr`. |
| `DESeq2.R` | Análisis diferencial con DESeq2; anotación de genes, selección de DEGs, PCA, volcano plot y heatmap. |
| `edgeR.R` | Normalización TMM, modelo binomial negativo con prueba quasi-likelihood, DEGs y visualizaciones. |
| `limma voom.R` | Normalización TMM, transformación voom con pesos de calidad, modelos lineales moderados, DEGs y visualizaciones. |
| `comparacion metodos.R` | Diagrama de Venn, genes compartidos entre métodos y recuperación de símbolos génicos. |
| `enriquecimiento funcional.R` | Enriquecimiento GO Biological Process mediante ORA y GSEA, tablas y gráficos. |

## Archivos de entrada

descomprimir data.rar dentro del proyecto:

- `counts_matrix.tabular`: matriz de counts; genes en filas y muestras en columnas. La primera columna contiene los identificadores de genes.
- `SraRunTable.csv`: metadata de las muestras, incluida la columna `Run` y la variable `disease_state`.

La matriz de counts procede del procesamiento de la **Clase 1** (alineamiento y cuantificación).

Estructura:

```text
proyecto/
├── data/
│   ├── counts_matrix.tabular
│   └── SraRunTable.csv
├── results/
├── plots/
├── RNAseq data preparation.R
├── DESeq2.R
├── edgeR.R
├── limma voom.R
├── comparacion metodos.R
└── enriquecimiento funcional.R
```

Crear las carpetas de salida antes de ejecutar los análisis manualmente, o con este codigo:

```r
dir.create("results", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)
```

## Orden de ejecución

### 1. Preparación de datos

Ejecutar **`RNAseq data preparation.R`**. El script:

1. Importa la matriz de counts y la metadata.
2. Limpia los identificadores Ensembl y ajusta, cuando hace falta, los nombres de las muestras.
3. Reordena la metadata para que coincida con las columnas de la matriz.
4. Define los grupos `Control` y `PD`.
5. Filtra genes de baja expresión mediante `edgeR::filterByExpr()`.

**Objetos principales:** `count.matrix`, `labels`, `sample.table`.


### 2. Expresión diferencial: tres métodos

Ejecutar, a partir de los mismos counts filtrados y la misma metadata:

- **`DESeq2.R`**: modelo binomial negativo, normalización por *median of ratios* y prueba de Wald; utiliza `rlog` y VST para PCA y heatmap.
- **`edgeR.R`**: normalización TMM, estimación de dispersión y prueba quasi-likelihood.
- **`limma voom.R`**: normalización TMM, `voomWithQualityWeights`, modelo lineal y moderación *empirical Bayes*.

**Objetos de resultados:** `res.DESeq2`, `res.edgeR`, `res.voom`; listas DEGs seleccionados: `DEGs.DESeq2`, `DEGs.edgeR`, `DEGs.voom`.

**Visualizaciones:** PCA, volcano plots y heatmaps. Los valores de un heatmap escalado por fila (Z-score) representan expresión relativa dentro de cada gen.

### 3. Comparación entre métodos

Ejecutar **`comparacion metodos.R`** después de los tres análisis. El script construye un diagrama de Venn y permite recuperar símbolos de los genes compartidos. Compara identificadores Ensembl presentes en las listas de DEGs.

### 4. Enriquecimiento funcional

Ejecutar **`enriquecimiento funcional.R`** una vez generados los objetos de expresión diferencial.

- **ORA (`enrichGO`)**: analiza una lista de DEGs frente a un universo de genes evaluados, utilizando términos GO Biological Process.
- **GSEA (`gseGO`)**: analiza el ranking completo de genes, sin limitarse a los DEGs significativos. Para DESeq2 se usa el estadístico Wald (`stat`); para limma-voom, el estadístico `t`.

El script genera tablas de resultados, dot plots y gráficos separados según la dirección del enriquecimiento. Un enriquecimiento positivo o negativo se interpreta según el contraste y el estadístico usados para ordenar los genes; **no equivale automáticamente a activación o inhibición funcional**.

## Adaptación a otros experimentos

Para reutilizar el práctico con otro dataset, revisar:

1. Rutas y formatos de la matriz de counts y la metadata.
2. Nombre de la columna que identifica cada muestra (`Run` en este ejemplo) y de la columna con los grupos (`disease_state`).
3. Etiquetas y niveles del factor experimental: el orden define la referencia y el signo del log2FC.
4. Tipo de identificador génico y organismo: `org.Hs.eg.db` corresponde a humano; para ratón hay que utilizar `org.Mm.eg.db`.
5. Diseño estadístico: los scripts actuales suponen una comparación sencilla entre **dos grupos**; experimentos con covariables, lotes o diseños pareados requieren modificar el diseño.

## Preguntas para discutir

- ¿Cuántos DEGs detecta cada método y cuántos comparten?
- ¿Las muestras se separan por condición en el PCA? ¿Qué otros factores podrían explicar la variación?
- ¿Qué patrones se observan en los heatmaps de DEGs?
- ¿En qué se diferencian ORA y GSEA? ¿Qué aporta cada enfoque?
- ¿Los cambios de expresión en tejido completo reflejan regulación dentro de las células, cambios en composición celular o ambos?

---
