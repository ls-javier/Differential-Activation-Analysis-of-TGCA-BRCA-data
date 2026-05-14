
#############-------------------################

## Analisis de Activacion Diferencial - TCGA-BRCA 
## Javier Luque Serrano 
## Master en Bioinformatica - Universidad Europea de Madrid

#############-------------------################

### Limpieza de los datos

#setwd("D:/Académico/Universitario/Master en Bioinformática/Cursos/
#Módulo 4. Metodología de la Investigación I/Actividades/
#Actividad 4. Análisis de activación diferencial")

DATA_DIR <- '.' #ruta relativa ('.' quiere decir "de aquí en adelante)
OUT_DIR  <- file.path(DATA_DIR, 'outputs') # cambiar path
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE) # crear directorio si no existe

# Cargamos el archivo
exp_mat <- read.delim("./BRCA_exp_matrix.tsv", header=T, sep="\t", fill=T,
                      stringsAsFactors = F)
dim(exp_mat)

# Hay valores NA en los nombres de las muestras
clinical <- read.delim("./clinical_info_TCGA-BRCA.tsv", header = TRUE, 
                       sep = "\t", check.names = T, fill = TRUE,
                       stringsAsFactors = F, row.names = NULL)

# Quitar valores NA - columna row.names tiene errores
clinical <- clinical[!is.na(clinical$bcr_patient_uuid), ]
# Eliminar duplicados mediante una fila que se repite entre muestras
clinical <- clinical[!duplicated(clinical$snames), ]
# Establecer nombre de las filas por orden
rownames(clinical) <- clinical[, 1]
clinical <- clinical[, -1]
rownames(clinical) <- clinical$sample

cat("Dimensiones (muestras x variables)", dim(clinical))

# Buscamos sólo muestras comunes
sum(colnames(exp_mat)%in%clinical$sample)
snames <- intersect(colnames(exp_mat), clinical$sample)

# Escogemos sólo las muestras comunes
exp_mat_clean <- exp_mat[ ,snames]
clinical_clean <- clinical[snames, ]

# Comprobamos la columna er_status_by_ihc 
cat("Grupos de ER:", table(clinical_clean$er_status_by_ihc))

# Eliminamos muestras sin valores en er
clinical_cleaner <- clinical_clean[
  !clinical_clean$er_status_by_ihc=="[Not Evaluated]" 
  &!clinical_clean$er_status_by_ihc=="Indeterminate", ]

table(clinical_cleaner$er_status_by_ihc)

# Nos quedamos con las mismas muestras en la matriz de expresión
snames_clean <- intersect(colnames(exp_mat_clean), clinical_cleaner$sample)
exp_mat_cleaner <- exp_mat_clean[, snames_clean]

# Comprobamos que todas las muestras coinciden en orden
all(rownames(clinical_cleaner)==colnames(exp_mat_cleaner))

# Convertimos a factor la columna de er_status
ER <- factor(clinical_cleaner$er_status_by_ihc, levels = c("Positive", "Negative"))

# Eliminar matrices crudas y renombrar las nuevas
clinical <- clinical_cleaner
exp_mat <- exp_mat_cleaner
rm(clinical_clean,clinical_cleaner, exp_mat_clean, exp_mat_cleaner)

### Normalizacion

# limma es requerido por edgeR
library(limma)
# Normalización y regresión
library(edgeR)

# Crear un objeto DGEList que edgeR entienda
y <- DGEList(exp_mat)
class(y)
str(y)

# Normalizamos en función de la cantidad de lecturas que tenemos
y <- calcNormFactors(y, method = "TMM")
str(y)

# Matriz de diseño para el ajuste de regresión lineal
design <- model.matrix(~ER)

# Quitamos los genes que son estables en la distribución de ambos grupos (caso y control), 
# ya que no aportan informacion
keep <- filterByExpr(y, design)
cat("Filtramos", sum(!keep), "genes") # Nos quitamos 3008 genes

# Nos quedamos con los genes de keep y todas las muestras, 
#manteniendo el tamaño de las librerias
y <- y[keep,,keep.lib.sizes=TRUE]
str(y)

# Hacemos una estimacion 
y <- estimateCommonDisp(y)

# Ajustamos la distribución al modelo
fit <- glmQLFit(y, design)

# Reproducimos las características que mejor se adaptan al modelo
qlf <- glmQLFTest(fit, coef = "ERNegative")

# Obtenemos las variables que mejor predicen el modelo
pvalues <- topTags(qlf, n=Inf, adjust.method ="BH", sort.by = "PValue", p.value=1)
# Convertir en dataframe
p_values <- as.data.frame(pvalues)

# Obtenemos los genes significativos
alpha <- 0.05
sig <- p_values[p_values$FDR<alpha, ]

# Sacamos una lista de los genes significativos
listGenes <- rownames(sig)
# Normalizamos los datos para tenerlos en una distribucion continua
data.norm <- cpm(y, log=TRUE, prior.count = 3)

# Observar dataset
cat("Dimensiones del dataset normalizado (genes x muestras): ", dim(data.norm))

# Hacemos un histograma con la expresion de los genes normalizados
hist(data.norm, 1000, main = paste("Histograma de la matriz de conteos normalizada "))

# Guardamos las tablas
write.table(data.norm, "./outputs/01_data_norm.tsv", col.names = T, 
            row.names = T, sep = "\t", quote = F)
write.table(listGenes, "./outputs/01_listGenes.tsv", col.names = F, 
            row.names = F, sep="\t", quote=F)

### Exploracion de los datos obtenidos

# Obtenemos la matriz de conteos sin normalizar
data <- cpm(exp_mat, log = TRUE, prior.count = 3)

# Tomar los indices de las 50 muestras elegidas al azar
set.seed(246)
samples_indices = sample(1:ncol(y), 50)
samples_indices

subset_norm = data.norm[, samples_indices]
subset = data[, samples_indices]

# Boxplot de datos no normalizados vs. normalizados
#png(file.path(OUT_DIR, '01_boxplots.png'), width=1400, height=700)
par(mfrow=c(1,2), oma = c(0, 0, 2, 0)) # oma añade margen externo 
boxplot(subset, las=2,  main="")
title(main="A. Datos crudos",ylab="Log-cpm")
boxplot(subset_norm, las=2,  main="")
title(main="B. Datos normalizados",ylab="Log-cpm")
mtext("Efecto de la normalizacion TMM sobre 50 muestras aleatorias", 
      side = 3, line = 0, outer = TRUE, cex = 1.3, font = 2) # Añade titulo principal
par(mfrow=c(1,1))

######################################
### Implementacion del metodo Hipathia
######################################

### Exploracion y control de calidad

# Previamente, instalamos hipathia desde el repo de GitHub para evitar problemas con load_pathways()
suppressPackageStartupMessages({
  library(hipathia)  # carga dependencias (SummarizedExperiment, etc.)
})
# Para convertir los IDs
library(AnnotationHub)

set.seed(246) # Semilla de reproducibilidad, para evitar que el 
# componente de aleatoriedad cambie los resultados
options(stringsAsFactors = FALSE) # Cambiar los parámetros con los que vamos a trabajar

# Convertir la matriz de expresion en un objeto matriz
exp_matrix <- as.matrix(data.norm)
mode(exp_matrix) <- "numeric"

# Construir una matriz de diseño con las muestras y su grupo
design2 <- data.frame(clinical$sample, ER)
colnames(design2) <- gsub("clinical.sample", "sample",
                          colnames(design2), fixed = T)
colnames(design2)
dim(design2)

# Crear subset de genes para graficar el boxplot
set.seed(246)
sel_genes <- sample(seq_len(nrow(exp_matrix)), 50)
boxplot(t(exp_matrix[sel_genes, ]), las=2, main='Distribucion sin normalizar - 50 genes',
        ylab='Expresion', outline=FALSE)

### Escalado y normalizado

# Construir un objeto SummarizedExperiment
se_raw <- SummarizedExperiment(
  assays = S4Vectors::SimpleList(raw = exp_matrix),
  colData = design2
)
str(se_raw)
save(se_raw, file = file.path(OUT_DIR, "02_SumExp_ER"))

# Traducir IDs
se_entrez <- translate_data(se_raw, species = "hsa")

# Normalizar a escala [0,1]
# Escalado global (by_gene=F) / by_quantiles=T puede estabilizar la distribucion
se_norm <- normalize_data(se_entrez, by_quantiles = TRUE, by_gene = FALSE)
str(se_norm)

## Boxplot tras normalización ##
# Matriz de expresion normalizada con los entrezIDs para grafico
norm_mat <- assay(se_norm)
# Definir el subset con los mismos genes que anteriormente
set.seed(246)
sel_genes2 <- sample(seq_len(nrow(norm_mat)), 50)
boxplot(t(norm_mat[sel_genes2, ]), las=2, main='Distribucion normalizada [0-1] - 50 genes',
        ylab='Expresion normalizada', outline=FALSE)

### Esta celda esta siendo ejecutada desde Docker -> hipathia 3.5 -> Bioconductor 3.18 -> R 4.3
# bioconductor/bioconductor_docker:RELEASE_3_18
library(hipathia)

# Cargar rutas de KEGG
#pathways <- load_pathways(species = "hsa")

# Ejecutar Hipathia: devuelve un objeto con nodos efectores y circuitos
#res <- hipathia(se_norm, pathways)

# Extraer la matriz de activacion de los circuitos
#path_mat <- get_paths_data(res, matrix = TRUE)
# Guardar objetos para que sean reutilizables
#saveRDS(pathways,      file = file.path(OUT_DIR, '03_pathways_object.rds'))
#saveRDS(res,      file = file.path(OUT_DIR, '03_results_object.rds'))
#saveRDS(path_mat, file = file.path(OUT_DIR, '03_path_activity_matrix.rds'))

pathways <- readRDS("D:/Académico/Universitario/Master en Bioinformática/Cursos/Módulo 4. Metodología de la Investigación I/Actividades/Actividad 4. Análisis de activación diferencial/outputs/03_pathways_object.rds")
path_mat <- readRDS("D:/Académico/Universitario/Master en Bioinformática/Cursos/Módulo 4. Metodología de la Investigación I/Actividades/Actividad 4. Análisis de activación diferencial/outputs/03_path_activity_matrix.rds")
res <- readRDS("D:/Académico/Universitario/Master en Bioinformática/Cursos/Módulo 4. Metodología de la Investigación I/Actividades/Actividad 4. Análisis de activación diferencial/outputs/03_results_object.rds")

dim(path_mat)

### Análisis de activación diferencial: Positive vs Negative

# Convertir la matriz de diseño en un vector
ERvec <- as.character(design2$ER)
names(ERvec) <- design2$sample
head(ERvec, 5)
table(ERvec)

# Wilcoxon por circuito (no parametrico)
# Ajustar los datos por multiples comparaciones
comp_paths <- do_wilcoxon(path_mat, ERvec[colnames(path_mat)], g1='Positive', g2='Negative',
                          paired=FALSE, adjust=TRUE, order=TRUE)  

write.csv(comp_paths, file = file.path(OUT_DIR, '04_wilcoxon_circuitos.csv'))

# Seleccionar circuitos significativos (FDR < 0.05 si existe esa columna)
fdr_col <- if ('FDRp.value' %in% colnames(comp_paths)) 'FDRp.value' else 'p.value'
sig_comp <- comp_paths[comp_paths[[fdr_col]] < 0.05, , drop=FALSE]
# Traducir a nombres legibles "nombre de ruta hsa: nombre gen efector"
sig_comp$path_name <- get_path_names(pathways, rownames(sig_comp))
write.csv(sig_comp, file = file.path(OUT_DIR, '04_wilcoxon_circuitos_significativos.csv'))

cat("Hay", nrow(sig_comp), "circuitos significativos (", fdr_col, " < 0.05). ")
cat("De los cuales,", table(sig_comp$`UP/DOWN`)[1] ,"son DOWNREGULADOS y",
    table(sig_comp$`UP/DOWN`)[2], "son UPREGULADOS")

# Rutas de señalizacion de los 10 circuitos mas significativos
sig_comp[1:10, ]

# Calculamos la varianza de cada circuito (filas en tu path_mat original)
varianzas <- apply(path_mat, 1, var)

# Filtramos la matriz para conservar SOLO los circuitos que varían (varianza > 0)
path_mat_limpia <- path_mat[varianzas > 0, ]
dim(path_mat_limpia)

# Calculamos el PCA creando un objeto princomp
pca_model <- do_pca(path_mat_limpia[seq_len(ncol(path_mat_limpia)),])

# Hacemos la grafica usando la columna "scores" que se creo en pca_model
pca_plot(pca_model, 
         group = ERvec[rownames(pca_model$scores)],
         main='PCA (actividad de circuitos) Positive vs Negative')

# Heatmap completo de los top circuitos por FDR/pvalue
library(pheatmap)

# Ordenar tabla de resultados por la columna del FDR (de menor a mayor)
comp_paths_sorted <- comp_paths[order(comp_paths$FDRp.value), ]
# Extraer el top 20
n_top <- 20 
sel_fdr <- head(rownames(comp_paths_sorted), n_top)
# Nos aseguramos de cruzar solo los que existen en la matriz
sel_fdr_validos <- sel_fdr[sel_fdr %in% rownames(path_mat)]

# Preparar la anotacion clinica 
# pheatmap necesita que la informacion de los pacientes sea un data.frame
# y que los nombres de las filas coincidan EXACTAMENTE con las columnas de la matriz
annotation_df <- data.frame(Fenotipo = ERvec)
rownames(annotation_df) <- colnames(path_mat)

png(file.path(OUT_DIR, '05_heatmap_paths.png'), width=700, height=500) # tambien se puede definir la resolucion
par(mar=c(7,7,7,7))

# Personalizar los colores de la leyenda clinica
ann_colors <- list(
  Fenotipo = c("Positive" = "#008080", "Negative" = "#87CEEB") 
)

# Renderizar el heatmap definitivo
pheatmap(
  mat = path_mat[sel_fdr_validos, ],  # Matriz filtrada con los top 20 circuitos
  annotation_col = annotation_df,  # Barra de pacientes (ER+ vs ER-)
  annotation_colors = ann_colors,         
  scale = "row",   # Convierte los valores a Z-score para máximo contraste visual
  color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(50), 
  cluster_rows = TRUE,  # Agrupa los circuitos biológicamente similares
  cluster_cols = TRUE,   # Agrupa a los pacientes "a ciegas" (tu validación)
  show_colnames = FALSE, # Ocultamos los IDs de los pacientes para que quede limpio
  fontsize_row = 10,
  main = "Firma Mecanística: Top 20 Circuitos Activados (ER- vs ER+)"
)
dev.off()

# Diferencia nodos
colors.de <- node_color_per_de(res, pathways, ER, "Positive","Negative") 

# Creamos y visualizamos el report, debemos abrir un navegador y pegar la 
# url que nos indique la función visualize_report()
url<-create_report(comp_paths, pathways, "NEG_vs_POS", node_colors = colors.de)
visualize_report(url, port=4054) #El número de puerto podemos ir cambiándolo
servr::daemon_stop(1) #para cerrar el puerto y poderlo usar de nuevo, comentado para knitr

message('FIN. Copia y pega la url proporcionada en un navegador para visualizar las rutas \n
        y los circuitos diferencialmente activados de manera interactiva.\n
        Importante cerrar la conexión al servidor antes de abrir otro puerto.')

sessionInfo()