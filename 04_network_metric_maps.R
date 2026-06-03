################################################################################
# 04_network_metric_maps.R
#
# Purpose:
#   Convert species ensemble maps to binary presence maps, then calculate spatial
#   maps of network metrics for each climate scenario.
################################################################################

project_root <- normalizePath(".", mustWork = FALSE)
setwd(project_root)

paths <- list(
  sdm_database_dir = file.path("outputs", "04_sdm_database"),
  ensemble_map_dir = file.path("outputs", "06_sdm_maps"),
  network_map_dir = file.path("outputs", "07_network_maps"),
  interaction_csv = file.path("data", "species_interaction2.csv"),
  threshold_csv = file.path("outputs", "07_network_maps", "species_thresholds_current.csv")
)

csv_encoding <- "UTF-8"
recalculate_thresholds <- FALSE
current_scenario_id <- "current_1991_2020"
network_metrics_to_map <- c("LD", "Generality", "Modularity", "pBasal", "Nest")

required_packages <- c("terra", "dismo", "igraph", "NetIndices", "bipartite")

load_required_packages <- function(package_names) {
  missing_packages <- package_names[!vapply(package_names, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop(
      "Install required packages before running this script: ",
      paste(missing_packages, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(lapply(package_names, library, character.only = TRUE))
}

load_required_packages(required_packages)
dir.create(paths$network_map_dir, recursive = TRUE, showWarnings = FALSE)

safe_mean <- function(values) {
  if (all(is.na(values))) {
    return(NA_real_)
  }
  mean(values, na.rm = TRUE)
}

safe_max <- function(values) {
  if (all(is.na(values))) {
    return(NA_real_)
  }
  max(values, na.rm = TRUE)
}

make_species_node_name <- function(species_name) {
  name_parts <- strsplit(species_name, " ")[[1]]
  paste0(substr(name_parts[1], 1, 1), ". ", name_parts[2])
}

calculate_species_thresholds <- function(paths, current_scenario_id) {
  threshold_table <- NULL
  database_files <- list.files(
    paths$sdm_database_dir,
    pattern = "^db_sdm_k[0-9]+_sp[0-9]{2}\\.RData$",
    full.names = TRUE
  )
  database_files <- database_files[order(database_files)]

  for (database_file in database_files) {
    load(database_file)

    current_map_file <- file.path(paths$ensemble_map_dir, paste0("ensemble_", species_id, "_", current_scenario_id, ".tif"))
    if (!file.exists(current_map_file)) {
      stop("Missing current ensemble map: ", current_map_file, call. = FALSE)
    }

    current_map <- terra::rast(current_map_file)
    extracted_values <- terra::extract(current_map, as.matrix(pa_data_original[, c("longitude", "latitude")]))
    prediction_column <- names(extracted_values)[2]

    evaluation_object <- dismo::evaluate(
      p = extracted_values[[prediction_column]][pa_data_original$occurrence == 1],
      a = extracted_values[[prediction_column]][pa_data_original$occurrence == 0]
    )

    species_threshold <- dismo::threshold(evaluation_object, "spec_sens")

    threshold_table <- rbind(
      threshold_table,
      data.frame(
        species_id = species_id,
        species = unique(pa_data_original$species),
        threshold = species_threshold,
        stringsAsFactors = FALSE
      )
    )
  }

  threshold_table$species_node <- vapply(threshold_table$species, make_species_node_name, character(1))
  threshold_table
}

calculate_network_metrics <- function(graph) {
  node_count <- igraph::gorder(graph)
  edge_count <- igraph::gsize(graph)

  if (node_count == 0) {
    return(data.frame(
      Size = 0,
      Omnivory = NA_real_,
      Links = 0,
      LD = NA_real_,
      Connectance = NA_real_,
      PathLength = NA_real_,
      Clustering = NA_real_,
      TLmean = NA_real_,
      TLmax = NA_real_,
      Vulnerability = NA_real_,
      VulSD = NA_real_,
      Generality = NA_real_,
      GenSD = NA_real_,
      pBasal = NA_real_,
      pTop = NA_real_,
      pInter = NA_real_,
      pCanni = NA_real_,
      Modularity = NA_real_,
      Nest = NA_real_
    ))
  }

  adjacency_matrix <- as.matrix(igraph::as_adjacency_matrix(graph, sparse = FALSE))
  trophic_levels <- tryCatch(NetIndices::TrophInd(adjacency_matrix)$TL, error = function(error) rep(NA_real_, node_count))

  vulnerability_values <- colSums(adjacency_matrix)
  generality_values <- rowSums(adjacency_matrix)
  diagonal_values <- diag(adjacency_matrix)

  data.frame(
    Size = node_count,
    Omnivory = safe_mean(trophic_levels - floor(trophic_levels)),
    Links = edge_count,
    LD = ifelse(node_count > 0, edge_count / node_count, NA_real_),
    Connectance = tryCatch(igraph::edge_density(graph), error = function(error) NA_real_),
    PathLength = tryCatch(igraph::mean_distance(graph, directed = FALSE, unconnected = TRUE), error = function(error) NA_real_),
    Clustering = tryCatch(igraph::transitivity(graph, type = "global"), error = function(error) NA_real_),
    TLmean = safe_mean(trophic_levels),
    TLmax = safe_max(trophic_levels),
    Vulnerability = safe_mean(vulnerability_values),
    VulSD = stats::sd(vulnerability_values),
    Generality = safe_mean(generality_values),
    GenSD = stats::sd(generality_values),
    pBasal = sum(rowSums(adjacency_matrix) == 0) / node_count,
    pTop = sum(colSums(adjacency_matrix) == 0) / node_count,
    pInter = 1 - (sum(rowSums(adjacency_matrix) == 0) / node_count) - (sum(colSums(adjacency_matrix) == 0) / node_count),
    pCanni = sum(diagonal_values > 0) / node_count,
    Modularity = tryCatch(igraph::modularity(igraph::cluster_louvain(graph)), error = function(error) NA_real_),
    Nest = tryCatch(bipartite::nested(adjacency_matrix, "NODF"), error = function(error) NA_real_)
  )
}

apply_species_thresholds <- function(species_rasters, thresholds) {
  binary_rasters <- species_rasters

  for (layer_index in seq_len(terra::nlyr(species_rasters))) {
    binary_rasters[[layer_index]] <- terra::ifel(species_rasters[[layer_index]] >= thresholds[layer_index], 1, 0)
  }

  binary_rasters
}

make_network_metric_map <- function(binary_maps, metric_name, reference_network, output_file) {
  calculate_cell_metric <- function(cell_values) {
    if (all(is.na(cell_values))) {
      return(NA_real_)
    }

    if (sum(cell_values, na.rm = TRUE) == 0) {
      return(0)
    }

    absent_species <- names(binary_maps)[!is.na(cell_values) & cell_values == 0]
    absent_species <- intersect(absent_species, igraph::V(reference_network)$name)
    local_network <- igraph::delete_vertices(reference_network, absent_species)
    as.numeric(calculate_network_metrics(local_network)[[metric_name]])
  }

  terra::app(binary_maps, fun = calculate_cell_metric, filename = output_file, overwrite = TRUE)
}

interaction_table <- utils::read.csv(paths$interaction_csv, fileEncoding = csv_encoding)
if (!all(c("Predator", "Prey") %in% names(interaction_table))) {
  stop("The interaction CSV must contain 'Predator' and 'Prey' columns.", call. = FALSE)
}

reference_network <- igraph::graph_from_data_frame(interaction_table, directed = FALSE)

if (recalculate_thresholds || !file.exists(paths$threshold_csv)) {
  threshold_table <- calculate_species_thresholds(paths, current_scenario_id)
  utils::write.csv(threshold_table, paths$threshold_csv, row.names = FALSE, fileEncoding = csv_encoding)
} else {
  threshold_table <- utils::read.csv(paths$threshold_csv, fileEncoding = csv_encoding)
}

if (!"species_node" %in% names(threshold_table)) {
  threshold_table$species_node <- vapply(threshold_table$species, make_species_node_name, character(1))
}

network_species <- sort(unique(c(interaction_table$Predator, interaction_table$Prey)))
threshold_species <- sort(threshold_table$species_node)
if (!identical(network_species, threshold_species)) {
  warning("Species names in the threshold table and network table do not fully match.")
}

model_count_files <- list.files(
  paths$ensemble_map_dir,
  pattern = "^model_counts_by_species_.*\\.csv$",
  full.names = TRUE
)
model_count_files <- model_count_files[order(model_count_files)]
scenario_ids <- sub("^model_counts_by_species_(.*)\\.csv$", "\\1", basename(model_count_files))

for (scenario_index in seq_along(scenario_ids)) {
  scenario_id <- scenario_ids[scenario_index]
  message("Calculating network maps for scenario: ", scenario_id)

  ensemble_files <- list.files(
    paths$ensemble_map_dir,
    pattern = paste0("^ensemble_sp[0-9]{2}_", scenario_id, "\\.tif$"),
    full.names = TRUE
  )
  ensemble_files <- ensemble_files[order(ensemble_files)]

  if (length(ensemble_files) == 0) {
    warning("No ensemble maps found for scenario: ", scenario_id)
    next
  }

  map_species_ids <- sub("^ensemble_(sp[0-9]{2})_.*", "\\1", basename(ensemble_files))
  threshold_match <- match(map_species_ids, threshold_table$species_id)

  if (any(is.na(threshold_match))) {
    stop("Missing threshold values for one or more species in scenario: ", scenario_id, call. = FALSE)
  }

  species_maps <- terra::rast(ensemble_files)
  species_thresholds <- threshold_table$threshold[threshold_match]
  binary_species_maps <- apply_species_thresholds(species_maps, species_thresholds)
  names(binary_species_maps) <- threshold_table$species_node[threshold_match]

  terra::writeRaster(
    binary_species_maps,
    filename = file.path(paths$network_map_dir, paste0("binary_species_maps_", scenario_id, ".tif")),
    overwrite = TRUE
  )

  for (metric_name in network_metrics_to_map) {
    output_file <- file.path(paths$network_map_dir, paste0("network_", scenario_id, "_", metric_name, ".tif"))
    make_network_metric_map(
      binary_maps = binary_species_maps,
      metric_name = metric_name,
      reference_network = reference_network,
      output_file = output_file
    )
  }
}
