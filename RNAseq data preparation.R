#RNAseq analysis 
library(dplyr)

#cargamos las counts
data <- read.delim("data/counts_matrix.tabular")

#cargamos la informacion de las muestras
metadata <- read.csv("data/SraRunTable.csv")

#chequear nombres de columnas de data
colnames(data) #vemos que esta duplicado el run
#ejecutar esto SI ES NECESARIO
colnames(data)[-1] <- substr(colnames(data)[-1], 1, 11)

#------------------------------------------------------
#para adaptar a otro experimento
#data <- read.csv2("data/counts_matrix_be.csv")
#metadata <- read.csv2("data/metadata.csv")
#------------------------------------------------------

#------------------------------------------------------------
# Preparamos la matriz de counts
#------------------------------------------------------------

# en las columnas solo tienen que quedar las counts
count.matrix <- data[, -1]

ensembl.id <- substr(data[,1],1,15) #limpiamos el codigo ensembl

rownames(count.matrix) <- ensembl.id #reasignamos el nombre a las filas

#limpiamos los datos y convertimos explicitamente a matriz
count.matrix = na.omit(count.matrix) # sacamos NAs

#convertimos a matriz
count.matrix <- as.matrix(count.matrix)

#chequeamos la matriz final
dim(count.matrix)
head(count.matrix)

#------------------------------------------------------------
# Ordenamos metadata según las columnas de la matriz
#------------------------------------------------------------
idx <- match(colnames(count.matrix), metadata$Run)

metadata.ordered <- metadata[idx, ]

# Chequeos
all(!is.na(idx))
identical(colnames(count.matrix), metadata.ordered$Run)

#------------------------------------------------------------
# creamos etiquetas de grupos y tabla de metadata de muestras
#------------------------------------------------------------

#etiquetas de los grupos
labels <- metadata.ordered$disease_state
#labels <- metadata.ordered$labels

# Simplificamos nombres si es necesario
labels <- ifelse(labels == "Non-demented control", "Control", "PD")

#convertimos a factor
labels <- factor(labels, 
                 levels = c("Control", "PD")) #ACA HAY QUE CAMBIAR SEGUN EL DATASET

#labels <- factor(labels, levels = c("control", "mrp4"))

# Chequeo
table(labels)

# creamos sample.table
sample.table <- data.frame(row.names = metadata.ordered$Run,
                           group = labels)

# Verificamos correspondencia
stopifnot(identical(rownames(sample.table), colnames(count.matrix)))

sample.table

#---------------------------------------------------------------------------
# Filtramos genes de baja expresion
#----------------------------------------------------------------------------

#keep <- rowSums(count.matrix) >= 30 #aplicamos un filtro basico

## Filter low-expression genes with edgeR
design <- model.matrix(~ labels)

keep <- edgeR::filterByExpr(edgeR::DGEList(counts = count.matrix), 
                                           design = design)

count.matrix <- count.matrix[keep, ]

dim(count.matrix)