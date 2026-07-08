# app.R - Main Shiny application (self-contained)

# Load required packages first (before UI definition)
suppressPackageStartupMessages({
  library(shiny)
  library(ggplot2)
  library(plotly)
  library(DT)
  library(pheatmap)
  library(data.table)
  library(edgeR)
  library(bslib)
})

# ============================================================
# DATA LOADING SECTION
# ============================================================

# Set working directory to project root
proj_root <- normalizePath(file.path(getwd(), ".."), mustWork = FALSE)
if (!dir.exists(file.path(proj_root, "bin"))) {
  proj_root <- getwd()
}

OUT_DIR <- file.path(proj_root, "bin", "outputs")
DATA_DIR <- file.path(proj_root, "data")
PATHWAY_DIR <- file.path(proj_root, "pdf", "rutas_señalizacion")

# Check required files exist
required_files <- c(
  file.path(OUT_DIR, "01_data_norm.tsv"),
  file.path(OUT_DIR, "01_listGenes.tsv"),
  file.path(OUT_DIR, "03_path_activity_matrix.rds"),
  file.path(OUT_DIR, "04_wilcoxon_circuitos.csv"),
  file.path(OUT_DIR, "04_wilcoxon_circuitos_significativos.csv"),
  file.path(DATA_DIR, "clinical_info_TCGA-BRCA.tsv")
)

for (f in required_files) {
  if (!file.exists(f)) {
    stop(paste("Required file not found:", f))
  }
}

# Load clinical data
clinical <- read.delim(file.path(DATA_DIR, "clinical_info_TCGA-BRCA.tsv"),
                       header = TRUE, sep = "\t", check.names = TRUE,
                       fill = TRUE, stringsAsFactors = FALSE, row.names = NULL)

# Clean clinical data
clinical <- clinical[!is.na(clinical$bcr_patient_uuid), ]
clinical <- clinical[!duplicated(clinical$snames), ]
rownames(clinical) <- clinical[, 1]
clinical <- clinical[, -1]
rownames(clinical) <- clinical$sample

# Filter by ER status
clinical <- clinical[
  !clinical$er_status_by_ihc == "[Not Evaluated]" &
  !clinical$er_status_by_ihc == "Indeterminate", ]

ER <- factor(clinical$er_status_by_ihc, levels = c("Positive", "Negative"))
ERvec <- as.character(ER)
names(ERvec) <- clinical$sample

# Load normalized expression data
message("Loading normalized expression matrix...")
data_norm <- read.delim(file.path(OUT_DIR, "01_data_norm.tsv"), header = TRUE, sep = "\t", row.names = 1)
data_norm <- as.matrix(data_norm)
mode(data_norm) <- "numeric"

# Load raw expression data for boxplot comparison
message("Loading raw expression matrix...")
exp_mat <- read.delim(file.path(DATA_DIR, "BRCA_exp_matrix.tsv"),
                      header = TRUE, sep = "\t", fill = TRUE,
                      stringsAsFactors = FALSE)

# Match samples
snames <- intersect(colnames(exp_mat), clinical$sample)
exp_mat <- exp_mat[, snames]
clinical_matched <- clinical[snames, ]

# Compute log-CPM for raw data
y_raw <- DGEList(exp_mat)
data_raw <- cpm(y_raw, log = TRUE, prior.count = 3)

# Load pathway activation matrix
message("Loading pathway activation matrix...")
path_mat <- readRDS(file.path(OUT_DIR, "03_path_activity_matrix.rds"))

# Subset ERvec to match path_mat columns
common_samples <- intersect(names(ERvec), colnames(path_mat))
ERvec <- ERvec[common_samples]

# Load Wilcoxon results
comp_paths <- read.csv(file.path(OUT_DIR, "04_wilcoxon_circuitos.csv"), row.names = 1)
sig_comp <- read.csv(file.path(OUT_DIR, "04_wilcoxon_circuitos_significativos.csv"), row.names = 1)

# Extract effector gene from path_name (after the colon)
sig_comp$effector <- sapply(strsplit(sig_comp$path_name, ": "), function(x) {
  if (length(x) > 1) gsub("\\*", "", trimws(x[2])) else "Unknown"
})

# Create DE table with proper column names
de_table <- sig_comp[order(sig_comp$FDRp.value), ]
de_table$Direction <- ifelse(de_table$UP.DOWN == "UP", "Upregulated", "Downregulated")

# Prepare PCA data
message("Computing PCA...")
varianzas <- apply(path_mat, 1, var)
path_mat_limpia <- path_mat[varianzas > 0, ]

# Transpose for PCA (samples as rows)
pca_data <- t(path_mat_limpia)
pca_result <- prcomp(pca_data, scale. = TRUE)

# Create PCA dataframe for plotting
pca_df <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Sample = rownames(pca_result$x),
  ER_Status = ERvec[rownames(pca_result$x)]
)

var_explained <- summary(pca_result)$importance[2, ] * 100

# Prepare heatmap data
comp_paths_sorted <- comp_paths[order(comp_paths$FDRp.value), ]
n_top <- 20
sel_fdr <- head(rownames(comp_paths_sorted), n_top)
sel_fdr_validos <- sel_fdr[sel_fdr %in% rownames(path_mat)]

# Annotation for heatmap
annotation_df <- data.frame(Fenotipo = ERvec)
rownames(annotation_df) <- colnames(path_mat)

ann_colors <- list(
  Fenotipo = c("Positive" = "#008080", "Negative" = "#87CEEB")
)

# Pathway images and KEGG links
pathway_info <- list(
  estrogenos = list(
    file = "estrogenos.png",
    name = "Estrogen Signaling Pathway",
    kegg_id = "hsa04915",
    KEGG_URL = "https://www.genome.jp/kegg-bin/show_pathway?map=04915",
    description = "Estrogen receptor signaling cascade. Hyperactivation converges on ESR1 and CREB3, indicating coordinated hormonal response and endoplasmic reticulum stress in proliferating cells."
  ),
  tiroidea = list(
    file = "tiroidea.png",
    name = "Thyroid Hormone Signaling Pathway",
    kegg_id = "hsa04919",
    KEGG_URL = "https://www.genome.jp/kegg-bin/show_pathway?map=04919",
    description = "Thyroid hormone signaling pathway. Cascades converge on ESR1 and CREB3 nodes, showing cross-talk with estrogen signaling."
  ),
  diabetes = list(
    file = "diabetes.png",
    name = "Type II Diabetes Mellitus Pathway",
    kegg_id = "hsa04930",
    KEGG_URL = "https://www.genome.jp/kegg-bin/show_pathway?map=04930",
    description = "Metabolic reprogramming pathway. Overexcitation from INSR through IRS1 activating MAPK and PI3K kinases, with altered insulin sensitivity and bioenergetic machinery."
  ),
  endocrinos = list(
    file = "endocrinos.png",
    name = "Endocrine System - Calcium Reabsorption",
    kegg_id = "hsa04912",
    KEGG_URL = "https://www.genome.jp/kegg-bin/show_pathway?map=04912",
    description = "Calcium homeostasis pathway. Vitamin D and estradiol signals drive TRPV5 and CALB1, altering ionic homeostasis for cellular motility and secretion."
  )
)

message("Data loading complete!")

# Prepare pathway browser data
pathway_browser <- data.frame(
  Circuit = rownames(comp_paths),
  Direction = comp_paths$UP.DOWN,
  Statistic = round(comp_paths$statistic, 2),
  `P-value` = formatC(comp_paths$p.value, format = "e", digits = 2),
  `FDR P-value` = formatC(comp_paths$FDRp.value, format = "e", digits = 2),
  stringsAsFactors = FALSE
)

# Extract pathway ID and name from circuit ID (P-hsaXXXXX-NN)
pathway_browser$Pathway_ID <- sapply(strsplit(pathway_browser$Circuit, "-"), function(x) x[2])
pathway_browser$KEGG_URL <- paste0("https://www.genome.jp/kegg/pathway/", pathway_browser$Pathway_ID)

# Get unique pathways
unique_pathways <- unique(pathway_browser[, c("Pathway_ID", "KEGG_URL")])
unique_pathways <- unique_pathways[order(unique_pathways$Pathway_ID), ]

# ============================================================
# UI DEFINITION
# ============================================================

ui <- fluidPage(
  theme = bslib::bs_theme(version = 5, bootswatch = "cosmo"),
  
  titlePanel(
    title = div(
      style = "display: flex; align-items: center; gap: 15px;",
      img(src = "UEM.jfif", height = 60, style = "border-radius: 8px;"),
      div(
        h2("TCGA-BRCA Differential Activation Analysis", style = "margin: 0;"),
        p("HiPathia Mechanistic Analysis: ER+ vs ER- Phenotypes", 
          style = "margin: 5px 0 0 0; color: #6c757d;")
      )
    )
  ),
  
  hr(),
  
  tabsetPanel(
    id = "main_tabs",
    type = "tabs",
    
    # Tab 1: Boxplots
    tabPanel(
      title = div(icon("chart-bar"), " Data Quality"),
      fluidRow(
        column(3,
          wellPanel(
            h4("Controls"),
            sliderInput(
              "n_genes",
              "Number of Genes:",
              min = 10, max = 200, value = 50, step = 10
            ),
            actionButton(
              "toggle_data",
              label = "Switch to Normalized",
              icon = icon("exchange-alt"),
              class = "btn-primary",
              style = "width: 100%;"
            ),
            hr(),
            p(strong("Note:"), "Random sampling uses seed 246 for reproducibility.")
          )
        ),
        column(9,
          plotOutput("boxplot", height = "600px")
        )
      )
    ),
    
    # Tab 2: Differential Expression
    tabPanel(
      title = div(icon("table"), " Differential Expression"),
      fluidRow(
        column(12,
          h3("Top 10 Upregulated Circuits (ER- vs ER+)"),
          DTOutput("table_up"),
          hr(),
          h3("Top 10 Downregulated Circuits (ER- vs ER+)"),
          DTOutput("table_down")
        )
      )
    ),
    
    # Tab 3: PCA
    tabPanel(
      title = div(icon("project-diagram"), " PCA Analysis"),
      fluidRow(
        column(12,
          h3("PCA of Circuit Activation Profiles"),
          p("Principal Component Analysis of pathway circuit activity, colored by ER status."),
          plotlyOutput("pca_plot", height = "550px")
        )
      )
    ),
    
    # Tab 4: Heatmap
    tabPanel(
      title = div(icon("th"), " Heatmap"),
      fluidRow(
        column(12,
          h3("Top 20 Activated Circuits (ER- vs ER+)"),
          p("Row-scaled Z-scores with hierarchical clustering. Teal = inhibited, Red = activated."),
          plotlyOutput("heatmap_plot", height = "650px")
        )
      )
    ),
    
    # Tab 5: Pathway Maps
    tabPanel(
      title = div(icon("sitemap"), " Pathway Maps"),
      fluidRow(
        column(12,
          h3("KEGG Signaling Pathway Browser"),
          p("Browse all 146 KEGG pathways analyzed in this study. Click pathway IDs to view detailed maps on KEGG."),
          hr()
        )
      ),
      fluidRow(
        column(12,
          h4("All Pathways (146 total)"),
          DTOutput("pathway_table")
        )
      ),
      hr(),
      fluidRow(
        column(12,
          h3("Featured Pathway Maps"),
          p("Key pathways identified as differentially activated between ER+ and ER- phenotypes."),
          hr()
        )
      ),
      uiOutput("pathway_maps")
    )
  ),
  
  hr(),
  tags$footer(
    div(
      style = "text-align: center; padding: 20px; color: #6c757d;",
      p("Javier Luque Serrano | Master's in Bioinformatics | European University of Madrid"),
      p("Module 4: Research Methodology I - Differential Activation Analysis")
    )
  )
)

# ============================================================
# SERVER LOGIC
# ============================================================

server <- function(input, output, session) {
  
  # Reactive value for data toggle
  data_type <- reactiveVal("raw")
  
  observeEvent(input$toggle_data, {
    if (data_type() == "raw") {
      data_type("normalized")
      updateActionButton(session, "toggle_data", label = "Switch to Raw")
    } else {
      data_type("raw")
      updateActionButton(session, "toggle_data", label = "Switch to Normalized")
    }
  })
  
  # Tab 1: Boxplot
  output$boxplot <- renderPlot({
    set.seed(246)
    
    if (data_type() == "raw") {
      sel_genes <- sample(seq_len(nrow(data_raw)), input$n_genes)
      plot_data <- t(data_raw[sel_genes, ])
      ylab_text <- "Log-CPM"
      main_text <- "Raw Gene Expression Distribution"
    } else {
      sel_genes <- sample(seq_len(nrow(data_norm)), input$n_genes)
      plot_data <- t(data_norm[sel_genes, ])
      ylab_text <- "Normalized Expression [0-1]"
      main_text <- "Normalized Gene Expression Distribution"
    }
    
    boxplot(
      plot_data,
      las = 2,
      main = main_text,
      ylab = ylab_text,
      outline = FALSE,
      col = ifelse(data_type() == "raw", "#FF6B6B", "#4ECDC4")
    )
  }, res = 150)
  
  # Tab 2: Differential Expression Tables
    output$table_up <- renderDT({
    up_genes <- de_table[de_table$UP.DOWN == "UP", ][1:10, ]
    
    datatable(
      up_genes[, c("path_name", "effector", "statistic", "FDRp.value", "UP.DOWN")],
      options = list(
        pageLength = 10,
        ordering = TRUE,
        columnDefs = list(
          list(className = 'dt-center', targets = '_all')
        )
      ),
      rownames = FALSE,
      colnames = c("Pathway", "Effector Gene", "Statistic", "FDR p-value", "Direction")
    ) %>%
      formatStyle(
        'UP.DOWN',
        backgroundColor = styleEqual('UP', '#FF6B6B'),
        color = 'white',
        fontWeight = 'bold'
      )
  })
  
  output$table_down <- renderDT({
    down_genes <- de_table[de_table$UP.DOWN == "DOWN", ][1:10, ]
    
    datatable(
      down_genes[, c("path_name", "effector", "statistic", "FDRp.value", "UP.DOWN")],
      options = list(
        pageLength = 10,
        ordering = TRUE,
        columnDefs = list(
          list(className = 'dt-center', targets = '_all')
        )
      ),
      rownames = FALSE,
      colnames = c("Pathway", "Effector Gene", "Statistic", "FDR p-value", "Direction")
    ) %>%
      formatStyle(
        'UP.DOWN',
        backgroundColor = styleEqual('DOWN', '#4DBBD5'),
        color = 'white',
        fontWeight = 'bold'
      )
  })
  
  # Tab 3: PCA Plot
  output$pca_plot <- renderPlotly({
    p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = ER_Status)) +
      geom_point(size = 4, alpha = 0.8) +
      scale_color_manual(values = c("Positive" = "#008080", "Negative" = "#87CEEB")) +
      labs(
        title = "PCA (Circuit Activity) ER+ vs ER-",
        x = paste0("PC1 (", round(var_explained[1], 1), "% variance)"),
        y = paste0("PC2 (", round(var_explained[2], 1), "% variance)"),
        color = "ER Status"
      ) +
      theme_minimal(base_size = 16) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right",
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_line(color = "grey95")
      )
    
    ggplotly(p, tooltip = c("Sample", "ER_Status")) %>%
      layout(
        dragmode = "zoom",
        xaxis = list(title = paste0("PC1 (", round(var_explained[1], 1), "%)")),
        yaxis = list(title = paste0("PC2 (", round(var_explained[2], 1), "%)"))
      )
  })
  
  # Tab 4: Heatmap
  output$heatmap_plot <- renderPlotly({
    # Get heatmap matrix
    hm_mat <- path_mat[sel_fdr_validos, ]
    
    # Row scale (Z-score)
    hm_mat_scaled <- t(scale(t(hm_mat)))
    
    # Create plotly heatmap
    p <- plot_ly(
      z = hm_mat_scaled,
      x = colnames(hm_mat),
      y = rownames(hm_mat),
      type = "heatmap",
      colorscale = list(
        c(0, "#4DBBD5"),
        c(0.5, "white"),
        c(1, "#E64B35")
      ),
      colorbar = list(
        title = "Z-score",
        tickvals = c(-2, -1, 0, 1, 2),
        ticktext = c("-2", "-1", "0", "1", "2")
      ),
      hovertemplate = "Circuit: %{y}<br>Sample: %{x}<br>Z-score: %{z:.2f}<extra></extra>"
    ) %>%
      layout(
        xaxis = list(
          title = "Samples",
          tickfont = list(size = 6),
          tickangle = -45
        ),
        yaxis = list(
          title = "Circuits",
          autorange = "reversed",
          tickfont = list(size = 9)
        ),
        margin = list(l = 150, r = 20, b = 120, t = 30)
      )
    
    # Add annotation bar for ER status
    er_colors <- ifelse(ERvec[colnames(hm_mat)] == "Positive", "#008080", "#87CEEB")
    
    p <- p %>%
      add_trace(
        type = "bar",
        x = seq_along(er_colors),
        y = rep(1, length(er_colors)),
        marker = list(color = er_colors),
        yaxis = "y2",
        showlegend = FALSE,
        hovertemplate = "Sample: %{customdata}<br>ER Status: %{marker.color}<extra></extra>",
        customdata = colnames(hm_mat)
      ) %>%
      layout(
        yaxis2 = list(
          overlaying = "y",
          side = "right",
          showticklabels = FALSE,
          showgrid = FALSE,
          zeroline = FALSE
        )
      )
    
    p
  })
  
  # Tab 5: Pathway Browser Table
  output$pathway_table <- renderDT({
    datatable(
      pathway_browser[, c("Circuit", "Pathway_ID", "Direction", "Statistic", "FDR P-value")],
      options = list(
        pageLength = 25,
        ordering = TRUE,
        search = list(regex = TRUE),
        columnDefs = list(
          list(className = 'dt-center', targets = '_all')
        )
      ),
      rownames = FALSE,
      colnames = c("Circuit", "Pathway ID", "Direction", "Statistic", "FDR P-value"),
      escape = FALSE
    ) %>%
      formatStyle(
        'Direction',
        backgroundColor = styleEqual(c('UP', 'DOWN'), c('#FF6B6B', '#4DBBD5')),
        color = 'white',
        fontWeight = 'bold'
      ) %>%
      formatStyle(
        'Pathway_ID',
        target = 'row',
        cursor = 'pointer'
      )
  })
  
  # Tab 5: Pathway Maps
  output$pathway_maps <- renderUI({
    lapply(names(pathway_info), function(name) {
      info <- pathway_info[[name]]
      kegg_url <- info$KEGG_URL
      
      fluidRow(
        column(12,
          wellPanel(
            h3(info$name),
            fluidRow(
              column(6,
                tags$a(
                  href = info$file,
                  target = "_blank",
                  tags$img(
                    src = info$file,
                    style = "width: 100%; border-radius: 8px; cursor: pointer;",
                    alt = info$name
                  )
                )
              ),
              column(6,
                h4("Description"),
                p(info$description),
                hr(),
                h4("KEGG Pathway"),
                tags$a(
                  href = kegg_url,
                  target = "_blank",
                  class = "btn btn-primary",
                  icon("external-link-alt"),
                  paste("View", info$name, "on KEGG")
                ),
                br(),
                br(),
                p(strong("KEGG ID:"), code(info$kegg_id))
              )
            )
          )
        )
      )
    })
  })
}

# Run the app
shinyApp(ui = ui, server = server)
