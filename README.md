# Differential Activation Analysis of TCGA-BRCA Data

---

## English

### Overview

This repository contains the code, data processing scripts, and results for a **Differential Activation Analysis** of breast cancer samples from the **TCGA-BRCA** dataset. The project applies a mechanistic approach using the **hiPathia** (*High Throughput Pathway Interpretation and Analysis*) method to infer functional cell context by integrating gene expression data with signaling pathway topology from **KEGG**.

The analysis compares **ER-positive (ER+)** versus **ER-negative (ER-)** luminal breast cancer samples to identify differentially activated signaling circuits, providing insights into hormonal, metabolic, and stress-related cross-talk mechanisms.

### Author & Affiliation

- **Author:** Javier Luque Serrano
- **Institution:** European University of Madrid
- **Degree:** Master's in Bioinformatics
- **Module:** Research Methodology I (Module 4)

### Objectives

- Normalize raw gene expression data using the **TMM** method via `edgeR`.
- Perform differential expression analysis to select significant genes.
- Translate gene identifiers to **Entrez IDs** and scale expression values to [0, 1].
- Compute pathway activation scores for **146 KEGG signaling pathways** using hiPathia.
- Identify significantly differentially activated circuits between ER+ and ER- phenotypes using the **Wilcoxon rank-sum test**.
- Visualize results through **PCA**, **heatmap clustering**, and interactive pathway maps.

### Repository Structure

```
.
├── app/                                         # Interactive Shiny web application
│   ├── app.R                                    # Main Shiny app (UI + server + data loading)
│   ├── www/                                     # Static assets
│   │   ├── UEM.jfif                             # University logo
│   │   ├── estrogenos.png                       # Estrogen pathway map
│   │   ├── tiroidea.png                         # Thyroid pathway map
│   │   ├── diabetes.png                         # Diabetes pathway map
│   │   └── endocrinos.png                       # Endocrine pathway map
│   └── README.md                                # App documentation
├── bin/                                         # Scripts and intermediate/output files
│   ├── Analisis_activacion_diferencial.R         # Main analysis script (R)
│   ├── analisis_activacion_diferencial_JLS.Rmd   # R Markdown report (Spanish)
│   ├── taller_hipathia_covid.R                   # Supplementary Hipathia workshop script
│   ├── image.RData                               # Large R workspace image (>100 MB, ignored)
│   └── outputs/                                  # Results and intermediate files
│       ├── 01_data_norm.tsv                      # Normalized expression matrix
│       ├── 01_listGenes.tsv                      # Significant gene list
│       ├── 02_SumExp_ER                          # SummarizedExperiment object
│       ├── 03_pathways_object.rds                # KEGG pathways (hiPathia)
│       ├── 03_path_activity_matrix.rds           # Circuit activation matrix
│       ├── 03_results_object.rds               # hiPathia results object
│       ├── 04_wilcoxon_circuitos.csv             # Wilcoxon test results (all circuits)
│       ├── 04_wilcoxon_circuitos_significativos.csv # Significant circuits only
│       ├── 04_heatmap_paths.png                  # Heatmap of top 20 circuits
│       ├── 01_boxplots.png                       # QC: raw vs normalized boxplots
│       ├── 02_boxplot_norm01.png                 # QC: normalized distribution
│       ├── data_norm.tsv                         # Alternative normalized matrix
│       └── listGenes.tsv                         # Alternative gene list
├── data/                                         # Raw input datasets
│   ├── BRCA_exp_matrix.tsv                       # Gene expression matrix (TCGA-BRCA)
│   └── clinical_info_TCGA-BRCA.tsv               # Clinical metadata
├── pdf/                                          # Reports and pathway figures
│   ├── analisis_activacion_diferencial_JLS.pdf   # Compiled PDF report
│   ├── U6-actividad-R-HIPATHIA.pdf               # Activity statement (ignored)
│   └── rutas_señalizacion/                       # Signaling pathway map screenshots
│       ├── diabetes.png
│       ├── endocrinos.png
│       ├── estrogenos.png
│       └── tiroidea.png
├── LICENSE                                       # MIT License
└── README.md                                     # This file
```

### Interactive Shiny App

An interactive web application is available in the `app/` directory for exploring the analysis results:

**Features:**
- **Data Quality Tab:** Boxplots with adjustable gene count (slider 10-200) and toggle between raw/normalized views
- **Differential Expression Tab:** Top 10 up/downregulated circuits with sortable, color-coded DT tables
- **PCA Tab:** Interactive PCA plot with zoom/pan, colored by ER status (teal=ER+, sky blue=ER-)
- **Heatmap Tab:** Interactive heatmap of top 20 circuits by FDR with zoom/pan, sample labels, and Z-score scaling
- **Pathway Maps Tab:** Browse all 1,876 circuits from 146 KEGG pathways in a searchable table, plus featured pathway maps with direct KEGG links

**Launch the app:**
```bash
# Windows
launch_app.bat

# Or from R
shiny::runApp("app", port = 7890, launch.browser = TRUE)
```

See `app/README.md` for detailed documentation.

### Methods Summary

1. **Data Cleaning:** Intersection of expression and clinical matrices; removal of unannotated/duplicated ER-status samples.
2. **Normalization:** TMM normalization (`edgeR::calcNormFactors`) and conversion to log-CPM.
3. **Differential Expression:** `glmQLFTest` with BH correction (FDR < 0.05).
4. **Hipathia Preprocessing:** ID translation to Entrez, quantile scaling [0,1], and `SummarizedExperiment` creation.
5. **Pathway Activation:** Signal propagation through 146 KEGG pathways to compute effector-circuit activation values.
6. **Statistical Testing:** Non-parametric Wilcoxon test with FDR correction comparing ER+ vs ER-.
7. **Visualization:** Unsupervised clustering heatmap and PCA of circuit activation profiles.

### Requirements

- R >= 4.0
- Bioconductor packages: `edgeR`, `limma`, `hipathia`, `SummarizedExperiment`, `AnnotationHub`
- CRAN packages: `shiny`, `ggplot2`, `plotly`, `DT`, `pheatmap`, `data.table`, `bslib`

> **Note:** The original hiPathia KEGG pathway download was performed inside a **Docker** container (`bioconductor/bioconductor_docker:RELEASE_3_18`) due to compatibility issues between modern R/Bioconductor versions and the required AnnotationHub resource.

---

## Español

### Descripción General

Este repositorio contiene el código, los scripts de procesamiento de datos y los resultados de un **Análisis de Activación Diferencial** de muestras de cáncer de mama del conjunto de datos **TCGA-BRCA**. El proyecto aplica un enfoque mecanicista mediante el método **hiPathia** (*High Throughput Pathway Interpretation and Analysis*) para inferir el contexto celular funcional integrando datos de expresión génica con la topología de las vías de señalización de **KEGG**.

El análisis compara muestras **luminales receptor de estrógenos positivo (ER+)** versus **receptor de estrógenos negativo (ER-)** para identificar circuitos de señalización diferencialmente activados, aportando información sobre los mecanismos de comunicación cruzada hormonal, metabólica y de estrés.

### Autor y Afiliación

- **Autor:** Javier Luque Serrano
- **Institución:** European University of Madrid (Universidad Europea de Madrid)
- **Título:** Máster en Bioinformática
- **Módulo:** Metodología de la Investigación I (Módulo 4)

### Objetivos

- Normalizar los datos de expresión génica mediante el método **TMM** con `edgeR`.
- Realizar un análisis de expresión diferencial para seleccionar genes significativos.
- Traducir los identificadores de genes a **Entrez IDs** y escalar los valores de expresión a [0, 1].
- Calcular los valores de activación de las vías para **146 vías de señalización de KEGG** utilizando hiPathia.
- Identificar circuitos significativamente diferencialmente activados entre los fenotipos ER+ y ER- mediante el **test de Wilcoxon**.
- Visualizar los resultados mediante **PCA**, **clustering heatmap** e informes interactivos de mapas de vías.

### Estructura del Repositorio

```
.
├── app/                                         # Aplicación web interactiva Shiny
│   ├── app.R                                    # Aplicación Shiny principal (UI + servidor + carga de datos)
│   ├── www/                                     # Recursos estáticos
│   │   ├── UEM.jfif                             # Logo de la universidad
│   │   ├── estrogenos.png                       # Mapa de ruta de estrógenos
│   │   ├── tiroidea.png                         # Mapa de ruta tiroidea
│   │   ├── diabetes.png                         # Mapa de ruta de diabetes
│   │   └── endocrinos.png                       # Mapa de ruta endocrina
│   └── README.md                                # Documentación de la app
├── bin/                                         # Scripts y archivos intermedios/finales
│   ├── Analisis_activacion_diferencial.R         # Script principal de análisis (R)
│   ├── analisis_activacion_diferencial_JLS.Rmd   # Informe en R Markdown (español)
│   ├── taller_hipathia_covid.R                   # Script complementario taller Hipathia
│   ├── image.RData                               # Imagen grande del espacio de trabajo R (>100 MB, ignorado)
│   └── outputs/                                  # Resultados y archivos intermedios
│       ├── 01_data_norm.tsv                      # Matriz de expresión normalizada
│       ├── 01_listGenes.tsv                      # Lista de genes significativos
│       ├── 02_SumExp_ER                          # Objeto SummarizedExperiment
│       ├── 03_pathways_object.rds                # Vías KEGG (hiPathia)
│       ├── 03_path_activity_matrix.rds           # Matriz de activación de circuitos
│       ├── 03_results_object.rds               # Objeto de resultados hiPathia
│       ├── 04_wilcoxon_circuitos.csv             # Resultados Wilcoxon (todos los circuitos)
│       ├── 04_wilcoxon_circuitos_significativos.csv # Circuitos significativos únicamente
│       ├── 04_heatmap_paths.png                  # Heatmap de los 20 principales circuitos
│       ├── 01_boxplots.png                       # QC: boxplots datos crudos vs normalizados
│       ├── 02_boxplot_norm01.png                 # QC: distribución normalizada
│       ├── data_norm.tsv                         # Matriz normalizada alternativa
│       └── listGenes.tsv                         # Lista de genes alternativa
├── data/                                         # Conjuntos de datos de entrada
│   ├── BRCA_exp_matrix.tsv                       # Matriz de expresión génica (TCGA-BRCA)
│   └── clinical_info_TCGA-BRCA.tsv               # Metadatos clínicos
├── pdf/                                          # Informes y figuras de vías
│   ├── analisis_activacion_diferencial_JLS.pdf   # Informe PDF compilado
│   ├── U6-actividad-R-HIPATHIA.pdf               # Enunciado de la actividad (ignorado)
│   └── rutas_señalizacion/                       # Capturas de pantalla de mapas de señalización
│       ├── diabetes.png
│       ├── endocrinos.png
│       ├── estrogenos.png
│       └── tiroidea.png
├── LICENSE                                       # Licencia MIT
└── README.md                                     # Este archivo
```

### Resumen de Métodos

1. **Limpieza de datos:** Intersección de matrices de expresión y clínicas; eliminación de muestras sin anotación o duplicadas según estado de ER.
2. **Normalización:** Normalización TMM (`edgeR::calcNormFactors`) y conversión a log-CPM.
3. **Expresión Diferencial:** `glmQLFTest` con corrección BH (FDR < 0,05).
4. **Preprocesado Hipathia:** Traducción de IDs a Entrez, escalado cuantílico [0,1] y creación de `SummarizedExperiment`.
5. **Activación de Vías:** Propagación de la señal a través de 146 vías de KEGG para calcular valores de activación de circuitos efectores.
6. **Test Estadístico:** Test no paramétrico de Wilcoxon con corrección FDR comparando ER+ vs ER-.
7. **Visualización:** Heatmap de clustering no supervisado y PCA de perfiles de activación de circuitos.

### Requisitos

- R >= 4.0
- Paquetes de Bioconductor: `edgeR`, `limma`, `hipathia`, `SummarizedExperiment`, `AnnotationHub`
- Paquetes de CRAN: `shiny`, `ggplot2`, `plotly`, `DT`, `pheatmap`, `data.table`, `bslib`

> **Nota:** La descarga original de las vías KEGG para hiPathia se realizó dentro de un contenedor **Docker** (`bioconductor/bioconductor_docker:RELEASE_3_18`) debido a problemas de compatibilidad entre versiones modernas de R/Bioconductor y el recurso de AnnotationHub requerido.
