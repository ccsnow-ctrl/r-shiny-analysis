library(shiny)
library(bslib)
library(tidyverse)
library(DT)
library(pheatmap)
library(ggrepel)
library(igraph)
library(showtext)

# load Press Start 2P for plots
font_add_google("Press Start 2P", "pressstart")
showtext_auto()

#  pre functions
make_summary <- function(df) {
  data.frame(
    Column = colnames(df),
    Type = sapply(df, class),
    Value = sapply(df, function(col) {
      if (is.numeric(col)) {
        paste0(round(mean(col), 2), " (+/- ", round(sd(col), 2), ")")
      } else {
        paste(unique(col), collapse = ", ")
      }
    })
  )
}

options(shiny.maxRequestSize = 32 * 1024^2) #counts was too big

# ui
ui <- page_fluid(
  theme = bs_theme(
    bg = "#cccccc", fg = "#0d0c0c", primary = "#dd2020",
    base_font = font_google("Press Start 2P"),
    code_font = font_google("Press Start 2P"),
    "font-size-base" = "0.7rem", "enable-rounded" = FALSE  # smaller font
  ) %>%
    bs_add_rules(
      list(
        sass::sass_file("nes.min.css"),
        "
        .nav-link { font-size: 0.6rem !important; }
        .btn { font-size: 0.6rem !important; }
        .btn { font-size: 0.65rem !important; }
        input[type='file']::file-selector-button { font-size: 0.65rem !important; }
        input[type='file'] { font-size: 0.65rem !important; }
        .form-control { font-size: 0.65rem !important; }
        .app-header { padding: 1rem 0; border-bottom: 3px solid #dd2020; margin-bottom: 1rem; }
        .tab-description { color: #555; font-size: 0.45rem; margin-bottom: 0.5rem; padding: 0.4rem; border-left: 3px solid #dd2020; }
    "
      )
    ),
  
  # header
  div(class = "app-header",
      p("Christine Snow | Final Project | BF530", style = "margin: 0; color: #dd2020;"),
      h2("RNA-seq Expression Analysis: Huntington's Disease vs. Neurologically Normal Controls (GSE64810)",
         style = "margin: 0.25rem 0 0 0; font-size: 0.55rem;")
  ),
  
  navset_tab(
    
    # tab 1, sample info
    nav_panel("Sample Info",
              sidebarLayout(
                sidebarPanel(
                  fileInput("file1", "Choose Sample Info CSV File", accept = ".csv"), 
                  selectInput("plot_col", "Select column to plot", choices = NULL)
                ),
                mainPanel(
                  # tab description
                  div(class = "tab-description",
                      p("Upload the sample information CSV to explore metadata across samples.
                      The Summary tab shows column types and values. The Data tab shows a sortable
                      table. The Plot tab shows distributions of continuous variables.")
                  ),
                  navset_card_underline(
                    nav_panel("Plot", plotOutput("plot")), #tab1a
                    nav_panel("Summary", tableOutput("summary")), #tab1b
                    nav_panel("Data", DT::dataTableOutput("data")) #tab1c
                  )
                )
              )
    ),
    
    # tab 2 counts
    nav_panel("Counts Exploration",
              sidebarLayout(
                sidebarPanel(
                  fileInput("file2", "Choose Counts CSV File", accept = ".csv"), #side bar
                  sliderInput("variance_slider", "Minimum variance percentile",
                              min = 0, max = 100, value = 50),
                  sliderInput("min_nonzero", "Minimum non-zero samples",
                              min = 0, max = 100, value = 10)
                ),
                mainPanel(
                  # tab description
                  div(class = "tab-description",
                      p("Upload the normalized counts CSV. Use the sliders to filter genes by
                      variance percentile and minimum non-zero samples. Plots update reactively
                      as filters change.")
                  ),
                  navset_card_underline(
                    nav_panel("Summary", tableOutput("summary2")), #tab 2a
                    nav_panel("Scatter Plot", plotOutput("scatterplot")),  # tab2b; median vs variance
                    nav_panel("Heatmap", plotOutput("heatmap")),       #tab 2c
                    nav_panel("PCA", plotOutput("pca")),          #tab2d
                    nav_panel("Data", DT::dataTableOutput("data2")) #tab2e
                  )
                )
              )
    ),
    
    # tab 3 diff exp
    nav_panel("Differential Expression",
              sidebarLayout(
                sidebarPanel(
                  fileInput("file3", "Choose Differential Expression CSV File", accept = ".csv"), #side bar
                ),
                mainPanel(
                  # tab description
                  div(class = "tab-description",
                      p("Upload the DESeq2 results CSV. The Summary tab shows a sortable table of
                      all DE results. The Volcano Plot shows log2 fold change vs -log10 adjusted
                      p-value, colored by significance (padj < 0.05).")
                  ),
                  navset_card_underline(
                    nav_panel("Summary", DT::dataTableOutput("summary3")), #tab 2a
                    nav_panel("Volcano Plot", plotOutput("volcanoplot")),  # tab2b; median vs variance #tab2e
                  )
                )
              )
    ),
    
    # tab 4 network analysis
    nav_panel("Network Analysis",       
              sidebarLayout(
                sidebarPanel(
                  fileInput("file4", "Choose Normalized Counts CSV File", accept = ".csv"),
                  textAreaInput("text1", "Enter Gene Name(s)", width = "100%", height = "150px", 
                                placeholder = "pls dont leave me empty :("),
                  sliderInput("corrslider", "correlation", min = 0, max = 1, value = .5),
                  actionButton("button1", "send it", icon = shiny::icon("bomb"))
                ),
                mainPanel(
                  # tab description
                  div(class = "tab-description",
                      p("Upload normalized counts, enter gene names one per line, set a minimum
                      correlation threshold, and hit send it. The network shows pairwise gene
                      correlations as edges. The stats table shows degree, closeness, and
                      betweenness centrality per gene.")
                  ),
                  navset_card_underline(
                    nav_panel("Stats", DT::dataTableOutput("data4")), 
                    nav_panel("Heatmap", plotOutput("heatmap4")), 
                    nav_panel("Network Plot", plotOutput("networkplot"))
                  )
                )
              )
    )                                 
  )                                    
)                                      


# server
server <- function(input, output, session) {
  
  #data
  
  # Tab 1: load sample info CSV
  data <- reactive({
    file <- input$file1
    req(file)
    ext <- tools::file_ext(file$datapath)
    validate(need(ext == "csv", "Please upload a csv file"))
    read.csv(file$datapath) %>% rename(accession = X)
  })
  
  # tab 2 reactives
  data2 <- reactive({
    file <- input$file2
    req(file)
    ext <- tools::file_ext(file$datapath)
    validate(need(ext == "csv", "Please upload a csv file"))
    read.csv(file$datapath)
  })
  
  counts_filtered2 <- reactive({
    req(data2())
    counts_only <- data2() %>% select(where(is.numeric)) %>% select(-GeneID)
    gene_var <- apply(counts_only, 1, var)
    nonzero <- rowSums(counts_only > 0)
    var_threshold <- quantile(gene_var, input$variance_slider / 100)
    keep <- gene_var >= var_threshold & nonzero >= input$min_nonzero
    counts_only[keep, ]
  })
  
  output$scatterplot <- renderPlot({
    req(counts_filtered2())
    gene_stats <- data.frame(
      median = apply(counts_filtered2(), 1, median),
      variance = apply(counts_filtered2(), 1, var)
    )
    ggplot(gene_stats, aes(x = median, y = variance)) +
      geom_point(alpha = 0.3, color = "#dd2020") +
      scale_x_log10() +
      scale_y_log10() +
      theme_bw() +
      theme(text = element_text(family = "pressstart")) + # match font
      labs(title = "Median Count vs Variance", x = "Median Count", y = "Variance")
  })
  
  #tab 3 reactives
  
  data3 <- reactive({
    file <- input$file3
    req(file)
    ext <- tools::file_ext(file$datapath)
    validate(need(ext == "csv", "Please upload a csv file"))
    read.csv(file$datapath)
  })
  
  output$volcanoplot <- renderPlot({
    req(data3())
    volcplot <- function(data, padj_threshold = 0.05, fc = 1, plot_title = 'Volcano Plot', plot_subtitle = NULL, genelist_vector = NULL, genelist_filter = FALSE) {
      
      # Set the fold-change thresholds
      neg_log2fc <- -log2(fc)
      pos_log2fc <- log2(fc)
      
      # Make a dataset for plotting, add the status as a new column
      plot_ready_data <- data %>%
        mutate_at('padj', ~replace(.x, is.na(.x), 1)) %>%
        mutate_at('log2FoldChange', ~replace(.x, is.na(.x), 0)) %>%
        mutate(
          log2fc_threshold = ifelse(log2FoldChange >= pos_log2fc & padj <= padj_threshold, 'up',
                                    ifelse(log2FoldChange <= neg_log2fc & padj <= padj_threshold, 'down', 'ns')
          )
        ) %>%
        mutate(Symbol = replace_na(Symbol, 'none'))
      
      if (genelist_filter) {
        plot_ready_data <- plot_ready_data %>% filter(Symbol %in% genelist_vector)
      }
      
      if(!is.null(genelist_vector)) {
        plot_ready_data <- plot_ready_data %>% mutate(Symbol = ifelse(Symbol %in% genelist_vector & padj < padj_threshold & log2fc_threshold != 'ns',Symbol, ''))
      }
      
      # Get the number of up, down, and unchanged genes
      up_genes <- plot_ready_data %>% filter(log2fc_threshold == 'up') %>% nrow()
      down_genes <- plot_ready_data %>% filter(log2fc_threshold == 'down') %>% nrow()
      unchanged_genes <- plot_ready_data %>% filter(log2fc_threshold == 'ns') %>% nrow()
      
      # Make the labels for the legend
      legend_labels <- c(
        str_c('Up: ', up_genes),
        str_c('NS: ', unchanged_genes),
        str_c('Down: ', down_genes)
      )
      
      # Set the x axis limits, rounded to the next even number
      x_axis_limits <- DescTools::RoundTo(
        max(abs(plot_ready_data$log2FoldChange)),
        2,
        ceiling
      )
      
      # updated plot colors to match theme
      plot_colors <- c(
        'up' = '#dd2020',      # matches primary red
        'ns' = '#c5c5c8',      # matches bg grey
        'down' = '#4a90d9'     # softer blue
      )
      
      
      # Make the plot, these options are a reasonable strting point
      plot <- ggplot(plot_ready_data) +
        geom_point(
          alpha = 0.25,
          size = 1.5
        ) +
        aes(
          x = log2FoldChange,
          y = -log10(padj),
          color = log2fc_threshold,
          label = Symbol
        ) +
        geom_vline(
          xintercept = c(neg_log2fc, pos_log2fc),
          linetype = 'dashed'
        ) +
        geom_hline(
          yintercept = -log10(padj_threshold),
          linetype = 'dashed'
        ) +
        scale_x_continuous(
          'log2(FC)',
          limits = c(-x_axis_limits, x_axis_limits)
        ) +
        scale_color_manual(
          values = plot_colors,
          labels = legend_labels
        ) +
        labs(
          color = str_c(fc, '-fold, padj ≤', padj_threshold),
          title = plot_title,
          subtitle = plot_subtitle
        ) +
        theme_bw(base_size = 24) +
        theme(
          text = element_text(family = "pressstart"), # match font
          axis.text = element_text(color = 'black'),
          legend.margin = margin(0, 0, 0, 0),
          legend.spacing.x = unit(0.2, 'cm')
        )
      
      # Add gene labels if needed
      if (!is.null(genelist_vector)) {
        plot <- plot +
          geom_label_repel(
            size = 6,
            force = 0.1,
            max.overlaps = 100000,
            nudge_x = 1,
            segment.color = 'black',
            min.segment.length = 0,
            show.legend = FALSE
          )
      }
      plot
    }
    
    volcplot(data3())
  })
  
  
  # tab 4 reactives
  
  data4 <- reactive({
    file <- input$file4
    req(file)
    ext <- tools::file_ext(file$datapath)
    validate(need(ext == "csv", "Please upload a csv file"))
    read.csv(file$datapath)
  })                                        
  
  network_data <- eventReactive(input$button1, {
    req(data4(), input$text1)
    genes <- unlist(strsplit(input$text1, "\n")) %>% trimws()
    counts_mat <- data4() %>% 
      column_to_rownames("Symbol") %>%
      select(where(is.numeric)) %>%
      select(-GeneID) 
    found <- genes[genes %in% rownames(counts_mat)]
    not_found <- genes[!genes %in% rownames(counts_mat)]
    subset_mat <- counts_mat[found, ]       
    cor_mat <- cor(t(subset_mat))            
    list(subset = subset_mat, cor_mat = cor_mat, not_found = not_found)
  })                                       
  
  # tab1 output
  
  # update plot column dropdown when file loads
  observe({
    updateSelectInput(session, "plot_col",
                      choices = colnames(data())[sapply(data(), is.numeric)])
  })
  
  # histogram of selected column
  output$plot <- renderPlot({
    req(input$plot_col)
    ggplot(data(), aes(x = .data[[input$plot_col]])) +
      geom_histogram(fill = "#dd2020", color = "#c5c5c8", bins = 30) +
      labs(title = paste("Distribution of", input$plot_col),
           x = input$plot_col, y = "Count") +
      theme_bw() +
      theme(text = element_text(family = "pressstart")) # match font
  })
  
  output$data <- DT::renderDataTable({ data() })
  output$summary <- renderTable({ make_summary(data()) })
  
  # tab 2 outputs
  
  # heatmap output 
  output$heatmap <- renderPlot({
    req(counts_filtered2())
    pheatmap(counts_filtered2(),
             show_rownames = FALSE,
             show_colnames = TRUE,
             scale = "row",
             color = colorRampPalette(c("#4a90d9", "#cccccc", "#dd2020"))(100))
  })
  
  # t2 pca output
  
  output$pca <- renderPlot({
    req(counts_filtered2())
    # run PCA - need to transpose so samples are rows
    pca_result <- prcomp(t(counts_filtered2()), scale. = TRUE)
    # get variance explained
    var_explained <- round(summary(pca_result)$importance[2,] * 100, 1)
    # build dataframe for plotting
    pca_df <- data.frame(
      PC1 = pca_result$x[,1],
      PC2 = pca_result$x[,2],
      sample = rownames(pca_result$x)
    )
    ggplot(pca_df, aes(x = PC1, y = PC2, label = sample)) +
      geom_point(size = 3, color = "#dd2020") +
      labs(
        title = "PCA",
        x = paste0("PC1 (", var_explained[1], "%)"),
        y = paste0("PC2 (", var_explained[2], "%)")
      ) +
      theme_bw() +
      theme(text = element_text(family = "pressstart")) # match font
  })
  
  # t2 data table output
  
  output$data2 <- DT::renderDataTable({ data2() })
  
  # t2 summary table output
  
  output$summary2 <- renderTable({
    req(data2(), counts_filtered2())
    counts_only <- data2() %>% select(where(is.numeric)) %>% select(-GeneID)
    data.frame(
      Metric = c("Number of samples", "Total genes", "Genes passing filter", "Genes not passing filter"),
      Value = c(
        ncol(counts_only),
        nrow(counts_only),
        nrow(counts_filtered2()),
        nrow(counts_only) - nrow(counts_filtered2())
      )
    )
  })
  
  #tab 3 output
  output$summary3 <- DT::renderDataTable({ data3() })
  
  #tab 4 outputs
  output$heatmap4 <- renderPlot({
    pheatmap(network_data()$subset,
             show_rownames = TRUE,
             show_colnames = TRUE,
             scale = "row",
             color = colorRampPalette(c("#4a90d9", "#cccccc", "#dd2020"))(100))
  })
  # network plot
  output$networkplot <- renderPlot({
    req(network_data())
    cor_mat <- network_data()$cor_mat
    
    # apply correlation threshold from slider
    adj_mat <- cor_mat
    adj_mat[abs(adj_mat) < input$corrslider] <- 0
    diag(adj_mat) <- 0
    
    # build igraph object
    g <- graph_from_adjacency_matrix(adj_mat, 
                                     mode = "undirected", 
                                     weighted = TRUE)
    par(bg = "#c0c0c0")
    
    # plot network
    plot(g, 
         vertex.label = V(g)$name,
         vertex.color = "#dd2020",
         vertex.label.color = "white",
         vertex.label.cex = 0.8,
         vertex.size = 35,
         edge.width = E(g)$weight * 2,
         edge.color = "#0d0c0c",        # add this - dark edge color
         main = "Gene Correlation Network",
         bg = "#c0c0c0")
  })
  
  # stats table
  output$data4 <- DT::renderDataTable({
    req(network_data())
    cor_mat <- network_data()$cor_mat
    adj_mat <- cor_mat
    adj_mat[abs(adj_mat) < input$corrslider] <- 0
    diag(adj_mat) <- 0
    g <- graph_from_adjacency_matrix(adj_mat, mode = "undirected", weighted = TRUE)
    
    data.frame(
      Gene = V(g)$name,
      Degree = degree(g),
      Closeness = round(closeness(g), 4),
      Betweenness = round(betweenness(g), 4)
    )
  })
}

# ---- RUN APP ----
shinyApp(ui = ui, server = server)