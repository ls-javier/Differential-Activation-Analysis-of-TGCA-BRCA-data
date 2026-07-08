# TCGA-BRCA Differential Activation Analysis - Interactive Shiny App

## Overview

Interactive web application for exploring the results of the differential activation analysis pipeline for TCGA-BRCA breast cancer samples (ER+ vs ER- phenotypes) using the hiPathia mechanistic approach.

## File Structure

```
app/
├── app.R          # Self-contained Shiny app (data loading + UI + server)
├── www/           # Static assets
│   ├── UEM.jfif   # University logo
│   ├── estrogenos.png   # Estrogen signaling pathway map
│   ├── tiroidea.png     # Thyroid hormone signaling pathway map
│   ├── diabetes.png     # Type II diabetes mellitus pathway map
│   └── endocrinos.png   # Endocrine system - calcium reabsorption map
└── README.md      # This file
```

## Launch the App

### Windows
Double-click `launch_app.bat` in the project root directory.

### From R Console
```r
setwd("path/to/Actividad 4. Análisis de activación diferencial")
shiny::runApp("app", port = 7890, launch.browser = TRUE)
```

### Command Line
```bash
"C:\Program Files\R\R-4.5.2\bin\Rscript.exe" -e "shiny::runApp('app', port = 7890, launch.browser = TRUE)"
```

The app will open automatically in your default browser at `http://localhost:7890`.

## App Features

### Tab 1: Data Quality
- **Boxplots** showing gene expression distributions
- **Slider** to adjust number of genes displayed (10-200, default: 50)
- **Toggle button** to switch between raw (log-CPM) and normalized [0-1] views
- Reproducible random sampling (seed = 246)

### Tab 2: Differential Expression
- **Top 10 upregulated** circuits (ER- vs ER+)
- **Top 10 downregulated** circuits (ER- vs ER+)
- Sortable DT tables with pathway name, effector gene, statistic, and FDR p-value
- Color-coded direction indicators (red=UP, blue=DOWN)

### Tab 3: PCA Analysis
- **Interactive PCA plot** of circuit activation profiles
- Points colored by ER status (teal = ER+, sky blue = ER-)
- **Zoom in/out** with mouse wheel or toolbar buttons
- **Pan/drag** to explore clusters
- Hover tooltips showing sample ID and ER status
- Axes show variance explained (PC1: 12.6%, PC2: 10.8%)

### Tab 4: Heatmap
- **Interactive heatmap** of top 20 circuits by FDR
- Row-scaled Z-scores for maximum contrast
- Color palette: teal (inhibited) → white → red (activated)
- **Zoom in/out** and **pan/drag** capabilities
- Sample names displayed with angled labels
- Hierarchical clustering of rows and columns

### Tab 5: Pathway Maps
- **Pathway Browser Table:** All 1,876 circuits from 146 KEGG pathways with search, sort, and color-coding
- **Featured Pathway Maps:**
  - Estrogen Signaling Pathway (hsa04915)
  - Thyroid Hormone Signaling Pathway (hsa04919)
  - Type II Diabetes Mellitus Pathway (hsa04930)
  - Endocrine System - Calcium Reabsorption (hsa04912)
- Clickable images for full-size view
- Direct links to KEGG pathway pages (`https://www.genome.jp/kegg-bin/show_pathway?map=XXXXX`)
- Biological significance descriptions

## Required R Packages

- shiny
- ggplot2
- plotly
- DT
- pheatmap
- data.table
- edgeR
- bslib

All packages must be installed before running the app.

## Data Sources

The app uses output files from the analysis pipeline:
- `bin/outputs/01_data_norm.tsv` - Normalized expression matrix
- `bin/outputs/03_path_activity_matrix.rds` - Circuit activation matrix
- `bin/outputs/04_wilcoxon_circuitos.csv` - Wilcoxon test results
- `bin/outputs/04_wilcoxon_circuitos_significativos.csv` - Significant circuits
- `app/www/*.png` - Pathway map images

## Troubleshooting

### App won't start
- Ensure all required R packages are installed
- Check that all output files exist in the correct locations
- Verify R version >= 4.0

### Plots not rendering
- Check browser console for JavaScript errors
- Try refreshing the page
- Ensure plotly package is installed

### Pathway images not loading
- Verify images exist in `app/www/` directory
- KEGG links will still work even if local images are missing

### KEGG links return 404
- The correct URL format is: `https://www.genome.jp/kegg-bin/show_pathway?map=XXXXX`
- KEGG changed their URL structure; old `genome.jp/kegg/pathway/` links no longer work

## Author

Javier Luque Serrano  
Master's in Bioinformatics  
European University of Madrid  
Module 4: Research Methodology I
