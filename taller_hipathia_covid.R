# ===============================================================
# Clase síncrona 7 - 23/01/2026 - 17:30-21:00h
# Taller  de análisis mecanístico con Hipathia
#
# María Peña    Fri Jan 23 03:51:48 2026
#
# ===============================================================

# ------------------------------
# 0) Setup (paquetes y carpetas)
# ------------------------------

# Estructura del directorio:
#   ./matiz4hipathia.txt
#   ./condition4hipathia.txt
#   ./phisiological_paths.tsv
#   ./taller_hipathia_covid.R   (este script)
#   ./outputs/                 (se crea automáticamente)


# Nota: Hipathia se instala desde Bioconductor.
# Si NO está instalado, se recomienda usar la versión más reciente:

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.21", force = T) # por defecto, R no va a downgradear paquetes ya instalados

BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))
BiocManager::install(c("BiocFileCache", "AnnotationHub"), dependencies = TRUE)

if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install('hipathia', force = T)

suppressPackageStartupMessages({
  library(hipathia)  # carga dependencias (SummarizedExperiment, etc.)
})

set.seed(246) # Semilla de reproducibilidad, para evitar que el componente de aleatoriedad cambie los resultados
options(stringsAsFactors = FALSE) # Cambiar los parámetros con los que vamos a trabajar

setwd('D:\\Académico\\Universitario\\Master en Bioinformática\\Cursos\\Módulo 4. Metodología de la Investigación I\\Unidad 6. Biología de sistemas y análisis mecanístico de datos ómicos\\Taller MPA')

DATA_DIR <- '.' #ruta relativa ('.' quiere decir "de aquí en adelante)
OUT_DIR  <- file.path(DATA_DIR, 'outputs') # cambiar path
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE) # crear directorio si no existe

expr_file <- file.path(DATA_DIR, 'matiz4hipathia.txt')
cond_file <- file.path(DATA_DIR, 'condition4hipathia.txt')

# ------------------------------
# 1) Cargar datos
# ------------------------------

# Matriz de expresión: filas = genes, columnas = muestras.
expr_df <- read.delim(expr_file, sep = '\t', header = TRUE, row.names = 1) # a veces funciona mejor que read_table. Cuantos mas argumentos, mejor leera la tabla
# Si tengo una RAM pequeña, es útil elegir un nrows bajito
# Si no pongo header = TRUE, R va a elegir nombres para las cabeceras

expr_mat <- as.matrix(expr_df) # posteriores funciones necesitan convertirlo en matriz (todos los elementos del mismo tipo)
mode(expr_mat) <- 'numeric'

# Diseño experimental: POS vs NEG
design <- read.delim(cond_file, sep = '\t', header = FALSE)

colnames(design) <- c('sample', 'group')
design$sample <- trimws(design$sample) #elimina posibles espacios (causan errores)
design$group  <- trimws(design$group)

# Comprobación de que los nombres de muestra están en la matriz de diseño
all(colnames(expr_mat) %in% design$sample) # Si es FALSE no podemos continuar con el analisis
all(colnames(expr_mat) == design$sample) # Esto ademas comprueba que estan en el mismo orden

# Reordenar el diseño para que coincida con el orden de columnas de la matriz
ord <- match(colnames(expr_mat), design$sample)
design <- design[ord, ]
stopifnot(all(design$sample == colnames(expr_mat)))

condition<-as.factor(design$group)
design$group <- factor(design$group, levels = c('NEG','POS')) # importante asignar los niveles, la comparacion puede no ser como queremos

# Estas 3 funciones muestran informacion en la consola. Interesante para un informe en RMarkDown
cat('Dimensiones (genes x muestras): ', nrow(expr_mat), ' x ', ncol(expr_mat))
message('Grupos:')
print(table(design$group))

# Guardar un log
write.csv(data.frame(n_genes=nrow(expr_mat), n_samples=ncol(expr_mat),
                     POS=sum(design$group=='POS'), NEG=sum(design$group=='NEG')),
          file = file.path(OUT_DIR, '00_resumen_dataset.csv'), row.names = FALSE)

# ----------------------------------
# 2) Exploración  / control QC
# ----------------------------------

# Boxplot de una submuestra de genes (para no saturar el gráfico)
sel_genes <- sample(seq_len(nrow(expr_mat)), 100)
png(file.path(OUT_DIR, '01_boxplot_raw.png'), width=1400, height=700) # tambien se puede definir la resolucion
par(mar=c(7,4,2,1))
boxplot(t(expr_mat[sel_genes, ]), las=2, main='Distribución (raw) - 100 genes aleatorios',
        ylab='Expresión (escala original)', outline=FALSE)
dev.off() # cerrar el archivo

# -------------------------------------------------------------
# 3) Preprocesado para Hipathia: escalado y normalización
# -------------------------------------------------------------

# Hipathia requiere entrez IDs en rownames y valores escalados [0,1]
# si no están en entrez ids, se usa translate_data 

# (A) Construir SummarizedExperiment. 
# Contiene muestras en bioconductor. Se puede hacer directamente sobre la matriz, pero este paso es interesante
se_raw <- SummarizedExperiment(
  assays  = S4Vectors::SimpleList(raw = expr_mat),
  colData = S4Vectors::DataFrame(group = design$group, row.names = design$sample)
)

str(se_raw) # ver estructura del objeto
# Si quiero acceder a cierto tipo de datos dentro del objeto
se_raw@colData@nrows

# (B) Traducir IDs (HGNC -> Entrez)
# IMPORTANTE: requiere acceso a recursos de anotación (AnnotationHub) la primera vez.
se_entrez <- translate_data(se_raw, species = 'hsa')

# (C) Normalizar/escala a [0,1]
# Por defecto: scaling global (by_gene = FALSE) recomendado por el vignette.
# by_quantiles = TRUE puede estabilizar distribuciones (opcional).
se_norm <- normalize_data(se_entrez, by_quantiles = TRUE, by_gene = FALSE)

# Boxplot tras normalización. Importante definir la escala para que R no la elija a su modo
norm_mat <- assay(se_norm)
sel_genes2 <- sample(seq_len(nrow(norm_mat)), 100)
png(file.path(OUT_DIR, '02_boxplot_norm01.png'), width=1400, height=700)
par(mar=c(7,4,2,1))
boxplot(t(norm_mat[sel_genes2, ]), las=2, main='Distribución normalizada [0,1] - 100 genes',
        ylab='Expresión normalizada', outline=FALSE)
dev.off()

# -------------------------------------------
# 4) Cálculo de activación de circuitos 
# -------------------------------------------

# Cargar vías KEGG fisiológicas disponibles (hsa)

path_list<-read.table("./physiological_paths.tsv",sep = "\t",stringsAsFactors = F,header = F,quote = "") # read.delim funciona igual
pathways<-load_pathways("hsa",pathways_list = path_list[,2]) #carga las rutas fisiológicas solo
#pathways <- load_pathways(species = 'hsa') #carga todas las rutas

#Si load_pathways falla por resource AH69111
data("pathways", package="hipathia") #carga un subset muy pequeño de circuitos integrados en el paquete sin llamar a igraph
   


# Ejecutar Hipathia: devuelve un objeto "results" con nodos efectores y circuitos
res <- hipathia(se_norm, pathways)

# Extraer la matriz de activación de los circuitos
path_mat <- get_paths_data(res, matrix = TRUE)

# Guardar por si se quiere reutilizar sin recalcular o continuar con el análisis más adelante
#OJO! si se ha usado data("pathways") la carga y cálculo será rápido, PERO si se usa KEGG entero puede tardar, por lo que es recomendable guardar
saveRDS(res,      file = file.path(OUT_DIR, '03_results_object.rds'))
saveRDS(path_mat, file = file.path(OUT_DIR, '03_path_activity_matrix.rds'))

# -------------------------------------------
# 5) Análisis mecanístico: POS vs NEG
#   (análisis de activación diferencial)
# -------------------------------------------

group_vec <- as.character(design$group)
names(group_vec) <- design$sample

# Wilcoxon por circuito (no paramétrico; robusto). No todos los conjuntos de datos son normales, hay casos que son binomiales
comp_paths <- do_wilcoxon(path_mat, group_vec[colnames(path_mat)], g1='POS', g2='NEG',
                          paired=FALSE, adjust=TRUE, order=TRUE)  # Ajustar los datos por multiples comparaciones

write.csv(comp_paths, file = file.path(OUT_DIR, '04_wilcoxon_circuitos.csv'))

# Seleccionar circuitos significativos (FDR < 0.05 si existe esa columna)
fdr_col <- if ('FDRp.value' %in% colnames(comp_paths)) 'FDRp.value' else 'p.value'
sig_comp <- comp_paths[comp_paths[[fdr_col]] < 0.05, , drop=FALSE]
write.csv(sig_comp, file = file.path(OUT_DIR, '04_wilcoxon_circuitos_significativos.csv'))

message('Circuitos significativos (', fdr_col, ' < 0.05): ', nrow(sig_comp))
# dim(sig_comp)
# nrow(sig_comp)

# -------------------------------------------
# 6) Visualizaciones: PCA + heatmap (top circuitos)
# -------------------------------------------

# PCA de muestras usando actividad de circuitos
pca_model <- do_pca(path_mat)

png(file.path(OUT_DIR, '05_pca_circuitos.png'), width=900, height=700)
pca_plot(pca_model, group = group_vec[rownames(pca_model$scores)],
         main='PCA (actividad de circuitos) POS vs NEG')
dev.off()

# Heatmap de los top circuitos por FDR/pvalue
n_top <- 18 #si hemos usado el subset de rutas, el número de significativos es menor, por lo que se puede reducir para que el heatmap sea visualmente más interpretable
sel <- head(rownames(comp_paths), n_top)

heatmap_plot(path_mat[sel, , drop=FALSE], group = group_vec[colnames(path_mat)],
             colors='hipathia', sample_clust=TRUE, variable_clust=TRUE,
             save_png = file.path(OUT_DIR, '06_heatmap_top50_circuitos.png'),
             main = paste0('Top ', n_top, ' circuitos (ordenados por ', fdr_col, ')'))


# Los IDs de circuitos codifican pares ruta/efector
# get_path_names traduce a nombres legibles "nombre de ruta hsa: nombre gen efector"
# Hacer antes del heatmap
comp_paths$path_name <- get_path_names(pathways, rownames(comp_paths))

write.csv(comp_paths, file = file.path(OUT_DIR, '07_wilcoxon_circuitos_con_nombres.csv'))


# -------------------------------------------
# 7) (Opcional) Activación de funciones (GO / Uniprot)
# -------------------------------------------

# Este paso nos permite agrupar la señal de los circuitos por funciones anotadas,
# calculando la activación o el nivel de actividad de funciones concretas, 
# puede ser más fácil de interpretar, pero proporciona información complementaria a los circuitos.

# Equivalente a un analisis de enriquecimiento mecanistico

# Podemos anotar vías a términos GO

# GO
fun_mat_go <- quantify_terms(res, pathways, "GO", out_matrix = TRUE)

# o Uniprot
fun_mat_uniprot <- quantify_terms(res, pathways, "uniprot", out_matrix = TRUE)


comp_fun <- do_wilcoxon(fun_mat_go, group_vec[colnames(fun_mat)], g1='POS', g2='NEG',
                        paired=FALSE, adjust=TRUE, order=TRUE)

# Traducir GO IDs a nombres (si aplica)
if ('GO' %in% colnames(comp_fun)) {
  comp_fun$GO_name <- get_go_names(comp_fun$GO)
}

write.csv(comp_fun, file = file.path(OUT_DIR, '08_wilcoxon_funciones_GO.csv'))

message('FIN. Revisa la carpeta outputs/ para resultados y figuras.')


# ------------------------------------------------
# 8) (Opcional) Cálculo de la activación con limma. Así ajustamos los datos a una distribucion lineal
# --------------------------------------------------

#Si no está instalada instalamos limma
#BiocManager::install("limma")

#Cargamos limma
library(limma)

#Volvemos a generar el diseño para modelo lineal (debe seguir la fórmula ~ y tener Intercept)
design=model.matrix(~condition) # si tuvieramos otra variable a tener en cuenta, habria que añadirla en la formula. Clave para que el analisis este bien hecho

#Ajustamos a un modelo lineal
fit <- lmFit (path_mat, design=design) #Matriz de disenyo
res.limma <- eBayes(fit) # transofrmacion bayesiana

#Obtenemos las tags significativas
comp <- topTable(res.limma, number = Inf,  coef="conditionPOS", adjust="fdr",p.value = 1)

# -------------------------------------------------------
# 9) (Opcional) Creación y visualización de report interactivo
# --------------------------------------------------------

#Para que la visualización sea correcta debemos ajustar aquellos circuitos sin variabilidad
equals <- apply(path_mat,1,function(x) all(x==x[1]))
names(equals) <- rownames(path_mat)

#Ajustamos el objeto de resultados para visualizar, esto debe hacerse porque hemos usado una función custom (limma) y no la incluida en el paquete hipathia
comp$P.Value[equals[comp$ID]] <- 1
comp$adj.P.Val[equals[comp$ID]] <- 1
comp$t[equals[comp$ID]] <- 0
comp$logFC[equals[comp$ID]] <- 0
comp <- comp[order(comp$P.Value,decreasing=F),]
path.names <- get_path_names(pathways, rownames(comp)) #Obtener los nombres de path
colnames(comp) <- gsub("adj.P.Val","FDRp.value",colnames(comp))
colnames(comp) <- gsub("P.Value","p.value",colnames(comp))
comp <- cbind(path.names, comp)##change path id to path name
comp$FDRp.value<-comp$p.value
comp <- comp[,c("path.names","t","logFC","p.value","FDRp.value")]
colnames(comp)[1]<-"path"


comp$"UP/DOWN" <- "UP"
comp$"UP/DOWN"[comp$logFC<0] <- "DOWN"
colors.de <- node_color_per_de(res, pathways, condition, "POS","NEG") #Diferencia nodos

comp_sig<-comp[comp$p.value<0.05,]
path_ef<-strsplit(":",comp_sig$path) #divide cadenas de caracteres en un vector
#Guardamos la tabla de resultados significativos
write.table(comp_sig,"09.pathways_sig_05_GEIS32_ALL_PFS.txt",col.names = T,row.names = T)

#Creamos y visualizamos el report, debemos abrir un navegador y pegar la url que nos indique la función visualize_report()
url<-create_report(comp, pathways, "NEG_vs_POS", node_colors = colors.de)
visualize_report(url, port=4000) #El número de puerto podemos ir cambiándolo
servr::daemon_stop(1) #para cerrar el puerto y poderlo usar de nuevo, el número depende de cuantos puertos hayamos usado

message('FIN. Copia y pega la url proporcionada en un navegador para visualizar las rutas y los \ncircuitos diferencialmente activados de manera interactiva, acuérdate de cerrar la conexión \nal servidor antes de abrir otro puerto.')
