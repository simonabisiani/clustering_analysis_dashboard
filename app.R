library(shiny)
library(bslib)
library(DT)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(googlesheets4)
library(colorspace)
library(RColorBrewer)
library(ggbeeswarm)
library(clValid)
library(cluster)
library(tibble)
library(ggradar)
library(scales)
library(purrr)

format_algorithm_name <- function(algorithm) {
  case_when(
    algorithm == "hc_ward" ~ "hierarchical with Ward linkage",
    algorithm == "hc_complete" ~ "hierarchical with complete linkage",
    algorithm == "kmeans" ~ "k-means",
    algorithm == "diana" ~ "DIANA",
    algorithm == "hdbscan" ~ "HDBSCAN",
    TRUE ~ algorithm
  )
}

# Calculate additional clustering metrics
calculate_metrics <- function(data, clusters) {
  # Remove outliers for calculations
  valid_idx <- clusters != 0
  data_clean <- data[valid_idx, ]
  clusters_clean <- clusters[valid_idx]

  # Within-cluster sum of squares
  wcss <- sum(sapply(unique(clusters_clean), function(k) {
    cluster_data <- data_clean[clusters_clean == k, ]
    if (nrow(cluster_data) > 1) {
      sum(dist(cluster_data)^2) / (2 * nrow(cluster_data))
    } else {
      0
    }
  }))

  list(wcss = wcss)
}

# Load all data at startup
load_all_data <- function() {
  cat("Loading data...\n")
  gs4_deauth()

  df <- read.csv("outlet_cluster_directory.csv", stringsAsFactors = FALSE)

  cluster_stats <- read_sheet(
    "https://docs.google.com/spreadsheets/d/1Y1xGVuFqMzaKnbQAQFgaIdgVJ1Ciik3UV5ougUzxXXE/edit?gid=1256911426#gid=1256911426",
    sheet = "clusters_statistics_transposed"
  )

  feature_cols <- c(
    "Area",
    "Radius",
    "Districts",
    "Entropy",
    "Gini",
    "MoranI",
    "DistCV",
    "Pct10km"
  )

  experiments <- list()
  unique_experiments <- unique(cluster_stats$Experiment)
  cluster_cols <- grep("^(X)?\\d{2}_", colnames(df), value = TRUE)

  exp_to_col <- list()
  for (col in cluster_cols) {
    col_clean <- gsub("^X", "", col)
    if (grepl("^\\d+_PCA", col_clean)) {
      # For PCA columns like "01_PCA_80_kmeans", we want "PCA_80_kmeans"
      exp_name <- sub("^\\d+_", "", col_clean)
    } else {
      # For regular columns like "01_Minimal_kmeans", we want "Minimal_kmeans"
      parts <- str_split(col_clean, "_", n = 4)[[1]] # Increased from 3 to 4
      if (length(parts) >= 3) {
        exp_name <- paste(parts[2:length(parts)], collapse = "_")
      } else {
        exp_name <- paste(parts[2], parts[3], sep = "_")
      }
    }
    exp_to_col[[exp_name]] <- col
  }

  for (exp_name in unique_experiments) {
    col_name <- exp_to_col[[exp_name]]
    if (is.null(col_name)) {
      next
    }

    exp_stats <- cluster_stats %>% filter(Experiment == exp_name)
    col_clean <- gsub("^X", "", col_name)
    parts <- str_split(col_clean, "_")[[1]]

    if (grepl("^\\d+_PCA", col_clean)) {
      dataset <- paste(parts[2], parts[3], sep = "_") # This captures "PCA_80"
      method <- parts[4]
    } else {
      dataset <- parts[2]
      method <- parts[3]
    }

    # Handle hierarchical clustering variants
    if (
      method == "hc" &&
        length(parts) > (if (grepl("^\\d+_PCA", col_clean)) 4 else 3)
    ) {
      linkage_method <- parts[if (grepl("^\\d+_PCA", col_clean)) 5 else 4]
      method <- paste0("hc_", linkage_method)
    }

    # Calculate additional metrics
    metrics <- calculate_metrics(df[, feature_cols], df[[col_name]])

    experiments[[col_name]] <- list(
      data = df[, feature_cols],
      clusters = df[[col_name]],
      algorithm = method,
      dataset = dataset,
      silhouette = exp_stats$Silhouette[1],
      n_clusters = exp_stats$N_Clusters[1],
      total_clustered = exp_stats$Total_Clustered[1],
      outliers = exp_stats$Outliers[1],
      cluster_sizes = exp_stats$Size,
      cluster_stats = exp_stats,
      full_name = exp_name,
      col_name = col_name,
      domain = df$Domain,
      wcss = metrics$wcss
    )
  }

  experiments <- experiments[order(sapply(experiments, function(x) x$dataset))]
  cat(paste("Loaded", length(experiments), "experiments\n"))

  return(list(
    experiments = experiments,
    df = df,
    cluster_stats = cluster_stats
  ))
}

perform_pca <- function(data) {
  scaled_data <- scale(data)
  pca <- prcomp(scaled_data, scale. = FALSE)
  list(pca = pca, data = as.data.frame(pca$x[, 1:2]))
}

CLUSTERING_DATA <- load_all_data()

dark_theme <- function() {
  theme_minimal() +
    theme(
      text = element_text(color = "#E43D12"),
      axis.text = element_text(color = "#E43D12", size = rel(1.1)),
      axis.title = element_text(color = "#E43D12", size = rel(1.1)),
      legend.text = element_text(color = "#E43D12"),
      legend.title = element_text(color = "#E43D12", face = "bold"),
      panel.grid.major = element_line(color = "#e1ded2", size = 0.5),
      strip.text = element_text(
        color = "#E43D12",
        face = "bold",
        size = rel(1.1),
        family = "Manrope"
      ),
      plot.title = element_text(
        color = "#E43D12",
        face = "bold",
        hjust = 0.5,
        size = rel(1.5)
      )
    )
}

thematic::thematic_shiny(
  font = "auto",
  bg = "auto",
  fg = "auto",
  qualitative = RColorBrewer::brewer.pal(8, "Set2")
)

ui <- page_navbar(
  title = "Clustering Analysis Dashboard",
  theme = bslib::bs_theme(
    bg = "#EBE9E1",
    fg = "#E43D12",
    primary = "#E33529",
    base_font = bslib::font_google("Manrope"),
    heading_font = bslib::font_google("Manrope"),
  ),
  bg = "#E33529",
  inverse = TRUE,

  # Add custom CSS for hover effects
  tags$head(
    tags$style(HTML("
      .experiment-button {
        transition: background-color 0.2s ease;
      }
      .experiment-button:hover {
        background-color: #f4a5a5 !important;
      }
      .experiment-button.selected:hover {
        background-color: #cc2a1e !important;
      }
    "))
  ),

  nav_panel(
    title = "Home",
    div(
      class = "container-fluid mt-4",
      div(
        h4("About this Analysis"),
        p(
          "This interactive dashboard accompanies a study exploring the spatial dimensions of local news coverage. We present the results of 25 clustering experiments applied to a set of spatial features derived from the locations mentioned in three continuous years (2020-2022) of news from 358 subnational media outlets from the UK. Each experiment represents a different configuration of input features and clustering algorithms (e.g., k-means, HDBSCAN), and aims to identify distinct types of local news outlets based on the spatial patterns of their reported locations."
        ),

        h5("Selected Experiment: Minimal Features + K-means"),
        p(
          "For the main study findings, we selected the",
          strong("'Minimal - kmeans'"),
          "experiment, which applies k-means clustering to a reduced set of minimally correlated spatial features. This configuration was chosen for its optimal balance between statistical quality and analytical interpretability. While other experiments achieved higher silhouette scores (particularly DIANA clustering with PCA-reduced features), they often produced severely imbalanced cluster solutions with most outlets concentrated in single dominant groups, limiting their analytical utility."
        ),

        p(
          "The minimal feature approach addresses multicollinearity issues inherent in spatial data while preserving the substantive meaning of geographic coverage patterns. This method successfully identified",
          strong("six distinct outlet typologies"),
          "spanning five orders of magnitude in coverage area: Hyperlocal, Local-Regional, Rural, Metropolitan, Sub-Regional, and National outlets. The k-means algorithm's centroid-based approach proved particularly suitable for capturing the hierarchical nature of spatial news organization, from neighborhood-focused hyperlocals (∼33 km²) to national broadcasters (∼79,000 km²)."
        ),

        h5("Interactive Features"),
        p("Users can:"),
        tags$ul(
          tags$li(
            "Select and explore different experiments, comparing PCA-based and feature-based visualisations of the cluster structure."
          ),
          tags$li(
            "Inspect cluster profiles, showing aggregate feature statistics for each group."
          ),
          tags$li(
            "Browse cluster members, filtering by individual outlet names."
          ),
          tags$li(
            "View feature distributions, to assess the variability of spatial indicators within and across clusters."
          )
        ),

        p(
          "This tool is designed to support transparency and interpretation in computational analyses of local media, offering researchers a way to investigate emergent patterns without enforcing a predefined typology. By comparing across experiments, users can assess the robustness of clustering solutions and understand how methodological choices influence the identification of media outlet spatial typologies."
        ),
      )
    )
  ),

  # Compare Experiments Tab
  nav_panel(
    title = "Compare Experiments",
      h2("Experiment Comparison"),
      p(
        "This table shows the performance metrics for all 25 clustering experiments conducted in the study."
      ),
      DTOutput("experiments_table")
    ),

# Explore One Experiment Tab
nav_panel(
  title = "Explore One Experiment",
  div(
    class = "container-fluid mt-4",
    
    # Experiment selection and summary side by side
    fluidRow(
      # Left column - Experiment selector grid
      column(
        7,
        card(
          card_header(
            h4("Select an Experiment", style = "margin: 0;")
          ),
          card_body(
            p("Click on a cell in the table below to select an experiment:", 
              style = "margin-bottom: 20px; color: #666;"),
            uiOutput("experiment_selector_table")
          )
        )
      ),
      
      # Right column - Experiment summary
      column(
        3,
        div(style = "margin-left: 20px;",
          uiOutput("exp_metrics") 
        )
      )
    ),

    # Add spacing before the analysis tabs
    div(style = "margin-top: 20px;",
      navset_card_underline(
        nav_panel(
          "Visualization",
          card_body(
            fluidRow(
              column(
                3,
                radioButtons(
                  "viz_type",
                  "Plot Type:",
                  choices = list(
                    "PCA Space" = "pca",
                    "Original Features" = "original"
                  ),
                  selected = "pca",
                  inline = TRUE
                )
              )
            ),
            conditionalPanel(
              condition = "input.viz_type == 'original'",
              fluidRow(
                column(
                  2,
                  div(
                    style = "text-align: left;",
                    selectInput(
                      "x_var",
                      "X Variable:",
                      choices = c(
                        "Area",
                        "Radius",
                        "Districts",
                        "Entropy",
                        "Gini",
                        "MoranI",
                        "DistCV",
                        "Pct10km"
                      ),
                      selected = "Area"
                    )
                  )
                ),
                column(
                  2,
                  div(
                    style = "text-align: left;",
                    selectInput(
                      "y_var",
                      "Y Variable:",
                      choices = c(
                        "Area",
                        "Radius",
                        "Districts",
                        "Entropy",
                        "Gini",
                        "MoranI",
                        "DistCV",
                        "Pct10km"
                      ),
                      selected = "Districts"
                    )
                  )
                ),
                column(
                  2,
                  div(
                    style = "text-align: left; padding-top: 25px;",
                    checkboxInput("log_scale", "Use log scale", FALSE)
                  )
                )
              )
            ),
            conditionalPanel(
              condition = "input.viz_type != 'dendrogram'",
              plotOutput("cluster_plot", height = "600px")
            ),
            conditionalPanel(
              condition = "input.viz_type == 'pca'",
              hr(),
              h5("PCA Feature Contributions"),
              plotOutput("pca_loadings", height = "400px")
            )
          )
        ),

        nav_panel(
          "Cluster Members",
          card_body(
            fluidRow(
              column(
                3,
                selectInput(
                  "selected_cluster",
                  "Select Cluster:",
                  choices = NULL,
                  width = "100%"
                )
              ),
              column(
                3,
                textInput(
                  "outlet_search",
                  "Search Outlets:",
                  placeholder = "Type outlet name...",
                  width = "100%"
                )
              )
            ),
            hr(),
            DTOutput("cluster_members")
          )
        ),

        nav_panel(
          "Cluster Profiles",
          card_body(
            h5("Cluster Characteristics"),
            DTOutput("cluster_profiles")
          )
        ),

        nav_panel(
          "Feature Distributions",
          card_body(
            h5("Distribution of Features by Cluster"),
            plotOutput("boxplots_plots", height = "375px")
          )
        )
      ) # Close navset_card_underline
    ) # Close div for spacing
  ) # Close main div
), # Close nav_panel


  # Methodology Tab
  nav_panel(
    title = "Methodology",
    div(
      class = "container-fluid mt-4",
      h4("Clustering Approach"),
      p(
        "This app presents the results of clustering experiments applied to UK local news outlets. Each experiment groups outlets based on spatial characteristics of their coverage, derived from place mentions in news content."
      ),

      h4("Input Variables"),
      tags$ul(
        tags$li(
          strong("Area:"),
          " Total area covered by outlet's reported locations (sq km)."
        ),
        tags$li(strong("Radius:"), " Average radius of coverage."),
        tags$li(
          strong("Districts:"),
          " Number of unique local authority districts mentioned."
        ),
        tags$li(
          strong("Entropy:"),
          " Distributional uniformity of location mentions."
        ),
        tags$li(strong("Gini:"), " Spatial inequality index."),
        tags$li(strong("Moran's I:"), " Spatial autocorrelation of coverage."),
        tags$li(
          strong("DistCV:"),
          " Coefficient of variation in distances between locations."
        ),
        tags$li(
          strong("Pct10km:"),
          " Percentage of mentions within 10 km of outlet centroid."
        )
      ),

      h4("Clustering Algorithms"),
      tags$ul(
        tags$li(
          strong("K-means:"),
          " Partitions outlets into clusters by minimizing within-cluster variance."
        ),
        tags$li(
          strong("Agglomerative Hierarchical:"),
          " Builds a tree of clusters using bottom-up merging."
        ),
        tags$li(
          strong("Diana:"),
          " A divisive hierarchical method starting with one large cluster split recursively."
        ),
        tags$li(
          strong("HDBSCAN:"),
          " A density-based clustering method that detects clusters of varying shapes and densities, and can identify outliers."
        )
      ),

      h4("Experimental Design"),
      p(
        "Each experiment represents a unique combination of input variables and clustering algorithms, optionally preceded by PCA for dimensionality reduction. The goal is to explore how varying these conditions affects the structure of resulting clusters."
      )
    )
  ),

  # Navigation spacer and links
  nav_spacer(),
  nav_menu(
    title = "About Me",
    align = "right",
    nav_item(tags$a(
      "Website",
      href = "https://simonabisiani.github.io",
      target = "_blank"
    )),
    nav_item(tags$a(
      "Bluesky",
      href = "https://bsky.app/profile/simonabisiani.bsky.social",
      target = "_blank"
    ))
  )
)

server <- function(input, output, session) {
  output$experiments_table <- renderDT({
    experiments <- CLUSTERING_DATA$experiments

    # Create summary table
    summary_data <- map_dfr(experiments, function(exp) {
      data.frame(
        Dataset = exp$dataset,
        Algorithm = format_algorithm_name(exp$algorithm),
        Silhouette = round(exp$silhouette, 3),
        N_Clusters = exp$n_clusters,
        Total_Clustered = exp$total_clustered,
        Outliers = exp$outliers,
        stringsAsFactors = FALSE
      )
    })

    # Sort by silhouette score descending
    summary_data <- summary_data[order(-summary_data$Silhouette), ]

    datatable(
      summary_data,
      options = list(pageLength = 15, scrollX = TRUE),
      class = 'cell-border stripe',
      rownames = FALSE
    ) %>%
      formatStyle(
        "Silhouette",
        background = styleColorBar(summary_data$Silhouette, '#E3B11D'),
        backgroundSize = '100% 90%',
        backgroundRepeat = 'no-repeat',
        backgroundPosition = 'center'
      )
  })

  # Add this reactive value to track the selected experiment
  values <- reactiveValues(selected_experiment = NULL)

output$experiment_selector_table <- renderUI({
  # Define cell dimensions once
  cell_width <- "120px"   # Make wider to fit text
  cell_height <- "40px"   # Make shorter/more rectangular
  
  experiments <- CLUSTERING_DATA$experiments
  
  # Get unique datasets and algorithms
  all_datasets <- unique(sapply(experiments, function(x) x$dataset))
  all_algorithms <- unique(sapply(experiments, function(x) x$algorithm))
  
  # Sort them for consistent ordering
  all_datasets <- sort(all_datasets)
  all_algorithms <- sort(all_algorithms)
  
  # Create a matrix to track which combinations exist
  experiment_matrix <- matrix(FALSE, nrow = length(all_datasets), ncol = length(all_algorithms))
  rownames(experiment_matrix) <- all_datasets
  colnames(experiment_matrix) <- all_algorithms
  
  # Fill in available experiments and store their keys
  experiment_keys <- matrix("", nrow = length(all_datasets), ncol = length(all_algorithms))
  silhouette_scores <- matrix(NA, nrow = length(all_datasets), ncol = length(all_algorithms))
  
  for(exp_key in names(experiments)) {
    exp <- experiments[[exp_key]]
    dataset_idx <- which(all_datasets == exp$dataset)
    algorithm_idx <- which(all_algorithms == exp$algorithm)
    
    if(length(dataset_idx) > 0 && length(algorithm_idx) > 0) {
      experiment_matrix[dataset_idx, algorithm_idx] <- TRUE
      experiment_keys[dataset_idx, algorithm_idx] <- exp_key
      silhouette_scores[dataset_idx, algorithm_idx] <- exp$silhouette
    }
  }
  
  # Create the grid HTML - using separate width/height
  grid_style <- paste0("display: grid; gap: 2px; grid-template-columns: 200px repeat(%d, ", cell_width, "); grid-template-rows: ", cell_height, " repeat(%d, ", cell_height, "); margin: 20px 0;")
  
  # Create header row
  header_cells <- list(
    # Corner cell with NO background color (transparent)
    div(style = "padding: 5px; font-weight: bold; text-align: center; font-size: 12px;", "Dataset / Algorithm")
  )
  
  for(algo in all_algorithms) {
    header_cells <- append(header_cells, list(
      div(style = "padding: 5px; font-weight: bold; text-align: center; font-size: 10px; word-wrap: break-word; overflow-wrap: break-word;", 
          format_algorithm_name(algo))
    ))
  }
  
  # Create data rows
  data_cells <- list()
  for(i in 1:length(all_datasets)) {
    dataset <- all_datasets[i]
    
    # Dataset name cell with NO background color (transparent)
    data_cells <- append(data_cells, list(
      div(style = "padding: 5px; font-weight: bold; text-align: right; font-size: 12px;", 
          dataset)
    ))

    # Algorithm cells
    for(j in 1:length(all_algorithms)) {
      if(experiment_matrix[i, j]) {
        # Available experiment - create clickable button
        button_id <- paste0("exp_", i, "_", j)
        exp_key <- experiment_keys[i, j]
        
        # Check if this is the currently selected experiment
        is_selected <- !is.null(values$selected_experiment) && values$selected_experiment == exp_key
        bg_color <- if(is_selected) "#E33529" else "#e1ded2"  # Red if selected, grey if not
        css_class <- if(is_selected) "experiment-button selected" else "experiment-button"
        
        data_cells <- append(data_cells, list(
          actionButton(
            button_id,
            label = "",
            class = css_class,
            style = paste0("width: ", cell_width, "; height: ", cell_height, "; border: 1px solid #ccc; background: ", bg_color, "; margin: 0; padding: 0;"),
            onclick = paste0("Shiny.setInputValue('selected_experiment_key', '", exp_key, "', {priority: 'event'});")
          )
        ))
      } else {
        # Unavailable combination - keep light grey
        data_cells <- append(data_cells, list(
          div(style = paste0("width: ", cell_width, "; height: ", cell_height, "; background: #f5f5f5; border: 1px solid #ddd; margin: 0;"))
        ))
      }
    }
  } # Close the i loop
  
  # Combine all cells
  all_cells <- append(header_cells, data_cells)
  
  div(
    style = sprintf(grid_style, length(all_algorithms), length(all_datasets)),
    all_cells
  )
})
# Handle button clicks
observeEvent(input$selected_experiment_key, {
  req(input$selected_experiment_key)
  values$selected_experiment <- input$selected_experiment_key
}, priority = 100)  # Higher priority to ensure it fires
  
  # Update grid colors when selection changes
observe({
  req(values$selected_experiment)
  # This will trigger a re-render of the grid with updated colors
  # The grid will automatically update because it's reactive to values$selected_experiment
})
  
# Initialize with the first available experiment
observe({
  if (is.null(values$selected_experiment)) {
    experiments <- CLUSTERING_DATA$experiments
    if (length(experiments) > 0) {
      values$selected_experiment <- names(experiments)[1]
    }
  }
})

  # Update the current_experiment reactive to use the new selection method
  current_experiment <- reactive({
    req(values$selected_experiment)
    CLUSTERING_DATA$experiments[[values$selected_experiment]]
  })

  # Current prototypes reactive
  current_prototypes <- reactive({
    req(current_experiment())
    exp <- current_experiment()

    exp$data[exp$clusters != 0, ] %>%
      mutate(cluster = exp$clusters[exp$clusters != 0]) %>%
      group_by(cluster) %>%
      summarise(across(everything(), median)) %>%
      mutate(cluster = paste("Cluster", cluster))
  })

  output$exp_metrics <- renderUI({
    exp <- current_experiment()
    req(exp)

    # Use the formatted algorithm name
    algorithm_name <- format_algorithm_name(exp$algorithm)

    dataset_footnote <- if (exp$dataset == "Minimal") {
      HTML(
        "<br><span style='font-size: 0.9em; color: #E43D12;'>Note: The 'Minimal' dataset includes the following variables: <b>Area</b>, <b>Districts</b>, <b>DistCV</b>, and <b>MoranI</b>.</span>"
      )
    } else {
      NULL
    }

    div(
      class = "p-3 rounded mt-3",
      style = "background-color: #e1ded2; width: 100%;",
      div(
        style = "display: flex; flex-direction: column;",
        h4(
          "Experiment Summary",
          class = "text-left",
          style = "color: #e33529;"
        ),
        div(
          class = "text-left",
          HTML(sprintf(
            "<p style='color: #e33529;'>Here below are the clustering results for an experiment which used the <b>%s</b> dataset and the <b>%s</b> algorithm. The model produced <b>%d</b> clusters, with a silhouette score of <b>%.3f</b>. A total of <b>%d</b> outlets were clustered, while <b>%d</b> were identified as outliers.</p>",
            exp$dataset,
            algorithm_name,
            exp$n_clusters,
            exp$silhouette,
            exp$total_clustered,
            exp$outliers
          ))
        ),
        dataset_footnote
      )
    )
  })

  observe({
    exp <- current_experiment()
    req(exp)

    clusters <- unique(exp$clusters[exp$clusters != 0])
    cluster_choices <- setNames(
      clusters,
      paste(
        "Cluster",
        clusters,
        "(n=",
        table(exp$clusters)[as.character(clusters)],
        ")"
      )
    )

    updateSelectInput(
      session,
      "selected_cluster",
      choices = cluster_choices,
      selected = clusters[1]
    )
  })

  output$cluster_plot <- renderPlot({
    exp <- current_experiment()
    req(exp)

    if (input$viz_type == "pca") {
      pca_result <- perform_pca(exp$data)
      plot_data <- pca_result$data %>%
        mutate(
          cluster = as.factor(exp$clusters),
          domain = exp$domain
        ) |>
        group_by(cluster) %>%
        mutate(
          cluster_label = paste0("Cluster ", cluster, " (n=", n(), ")")
        ) %>%
        ungroup()

      p <- ggplot(plot_data, aes(x = PC1, y = PC2, color = cluster)) +
        geom_point(size = 3, alpha = 0.7) +
        labs(
          x = paste0(
            "PC1 (",
            round(summary(pca_result$pca)$importance[2, 1] * 100, 1),
            "% var)"
          ),
          y = paste0(
            "PC2 (",
            round(summary(pca_result$pca)$importance[2, 2] * 100, 1),
            "% var)"
          ),
          color = "Cluster"
        ) +
        dark_theme()
    } else {
      plot_data <- exp$data %>%
        mutate(
          cluster = as.factor(exp$clusters),
          domain = exp$domain,
          x_val = .data[[input$x_var]],
          y_val = .data[[input$y_var]]
        )

      if (input$log_scale) {
        plot_data <- plot_data %>%
          mutate(
            x_val = log10(x_val + 1),
            y_val = log10(y_val + 1)
          )
      }

      x_label <- ifelse(
        input$log_scale,
        paste0("log10(", input$x_var, " + 1)"),
        input$x_var
      )
      y_label <- ifelse(
        input$log_scale,
        paste0("log10(", input$y_var, " + 1)"),
        input$y_var
      )

      p <- ggplot(plot_data, aes(x = x_val, y = y_val, color = cluster)) +
        geom_point(size = 3, alpha = 0.7) +
        labs(
          x = x_label,
          y = y_label,
          color = "Cluster"
        ) +
        dark_theme()
    }

    p
  })

  output$pca_loadings <- renderPlot({
    exp <- current_experiment()
    req(exp)

    pca_result <- perform_pca(exp$data)
    loadings <- as.data.frame(pca_result$pca$rotation[, 1:2])
    loadings$feature <- rownames(loadings)

    ggplot(loadings, aes(x = PC1, y = PC2, label = feature)) +
      geom_segment(
        aes(x = 0, y = 0, xend = PC1, yend = PC2),
        arrow = arrow(length = unit(0.3, "cm")),
        color = "#E3B11D",
        size = 1
      ) +
      geom_text(color = "#e33529", nudge_x = 0.02, nudge_y = 0.02, size = 5) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
      geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
      labs(x = "PC1 Loading", y = "PC2 Loading") +
      dark_theme()
  })

  output$cluster_members <- renderDT({
    exp <- current_experiment()
    req(exp, input$selected_cluster)

    cluster_num <- as.numeric(input$selected_cluster)

    members_data <- data.frame(
      Domain = exp$domain[exp$clusters == cluster_num],
      stringsAsFactors = FALSE
    )

    feature_data <- exp$data[exp$clusters == cluster_num, ]
    members_data <- cbind(members_data, round(feature_data, 2))

    if (!is.null(input$outlet_search) && input$outlet_search != "") {
      members_data <- members_data %>%
        filter(grepl(input$outlet_search, Domain, ignore.case = TRUE))
    }

    datatable(
      members_data,
      options = list(pageLength = 10, scrollX = TRUE),
      class = 'cell-border stripe',
      rownames = FALSE
    )
  })

  output$cluster_profiles <- renderDT({
    exp <- current_experiment()
    req(exp)

    profile_data <- exp$cluster_stats %>%
      select(
        Cluster,
        Size,
        Area,
        Radius,
        Districts,
        Entropy,
        Gini,
        MoranI,
        DistCV,
        Pct10km
      )

    if (
      "Label" %in%
        names(exp$cluster_stats) &&
        !all(is.na(exp$cluster_stats$Label))
    ) {
      profile_data <- profile_data %>%
        mutate(Label = exp$cluster_stats$Label) %>%
        select(Cluster, Label, everything())
    }

    numeric_cols <- c(
      "Area",
      "Radius",
      "Districts",
      "Entropy",
      "Gini",
      "MoranI",
      "DistCV",
      "Pct10km"
    )
    profile_data[numeric_cols] <- round(profile_data[numeric_cols], 3)

    datatable(
      profile_data,
      options = list(pageLength = 10, scrollX = TRUE),
      class = 'cell-border stripe',
      rownames = FALSE
    )
  })

  output$boxplots_plots <- renderPlot({
    exp <- current_experiment()
    req(exp)

    plot_data <- exp$data %>%
      mutate(Cluster = as.factor(exp$clusters)) %>%
      filter(Cluster != "0")

    # Standardize the features (z-score normalization)
    plot_data_scaled <- plot_data %>%
      mutate(across(-Cluster, ~ scale(.x)[, 1])) %>%
      pivot_longer(cols = -Cluster, names_to = "Feature", values_to = "Value")

    ggplot(plot_data_scaled, aes(x = Feature, y = Value)) +
      geom_boxplot(alpha = 0.4, outlier.size = 0.8, outlier.alpha = 0.3) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) + # Exactly 5 breaks
      coord_flip() +
      facet_wrap(~Cluster, scales = "free_x", ncol = 6) +
      labs(x = "", y = "Standardised Value") +
      dark_theme() +
      theme(
        panel.spacing = unit(2, "lines")
      )
  })
}

shinyApp(ui = ui, server = server)

