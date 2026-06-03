################################################################################
# 01_data_preprocessing.R
#
# Purpose:
#   1. Load environmental raster layers.
#   2. Build species-by-environment occurrence tables.
#   3. Create species buffers and background points.
#   4. Build presence-background databases.
#   5. Split each database into repeated train/test sets for SDM modeling.
################################################################################

project_root <- normalizePath(".", mustWork = FALSE)
setwd(project_root)

paths <- list(
  current_bioclim_dir = file.path("data", "worldclim_current"),
  current_bioclim_cropped = file.path("data", "worldclim_current", "cropped_bc_30ys_1991_2020.tif"),
  landcover_dir = file.path("data", "landcover"),
  gbif_rdata = file.path("data", "gbif", "gbif_records.RData"),
  buffer_dir = file.path("outputs", "01_buffers"),
  background_dir = file.path("outputs", "02_background_points"),
  presence_absence_dir = file.path("outputs", "03_presence_absence_db"),
  sdm_database_dir = file.path("outputs", "04_sdm_database")
)

csv_encoding <- "UTF-8"
run_landcover_percent_conversion <- FALSE
run_presence_extraction <- TRUE
run_buffer_generation <- TRUE
run_background_sampling <- TRUE
run_presence_absence_build <- TRUE
run_train_test_split <- TRUE

background_point_count <- 15000
selected_background_count <- 10000
spatial_weight <- 0.30
repeat_count <- 10
test_fraction <- 0.20
base_seed <- 910520

required_packages <- c(
  "terra",
  "rnaturalearth",
  "sf",
  "megaSDM",
  "dplyr",
  "usdm"
)

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
invisible(lapply(paths[c("buffer_dir", "background_dir", "presence_absence_dir", "sdm_database_dir")], dir.create, recursive = TRUE, showWarnings = FALSE))

target_species <- c(
  "Ciconia episcopus",
  "Leptoptilos dubius",
  "Duttaphrynus melanostictus",
  "Hoplobatrachus tigerinus",
  "Lethocerus indicus",
  "Channa punctata",
  "Fowlea piscator",
  "Gryllus bimaculatus",
  "Calotes versicolor",
  "Pomacea canaliculata"
)

valid_basis_of_record <- c(
  "HUMAN_OBSERVATION",
  "MACHINE_OBSERVATION",
  "OCCURRENCE",
  "OBSERVATION"
)

make_species_id <- function(species_index) {
  paste0("sp", sprintf("%02d", species_index))
}

make_species_slug <- function(species_name) {
  gsub("_+$", "", gsub("[^A-Za-z0-9]+", "_", species_name))
}

create_landcover_percent_rasters <- function(raw_landcover_dir, reference_raster, output_dir = raw_landcover_dir) {
  raw_landcover_files <- list.files(raw_landcover_dir, pattern = "\\.tif$", full.names = TRUE)
  raw_landcover_files <- raw_landcover_files[!grepl("Percent", basename(raw_landcover_files), ignore.case = TRUE)]

  if (length(raw_landcover_files) == 0) {
    stop("No raw land-cover .tif files were found.", call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  for (landcover_file in raw_landcover_files) {
    message("Converting land-cover raster to percentage layer: ", basename(landcover_file))
    landcover_raster <- terra::rast(landcover_file)
    aggregation_factor <- terra::res(reference_raster)[1] / terra::res(landcover_raster)[1]

    percent_raster <- terra::aggregate(
      landcover_raster,
      fact = aggregation_factor,
      fun = function(values, ...) mean(values == 1, na.rm = TRUE) * 100
    )

    names(percent_raster) <- paste0("Percent_", names(landcover_raster))
    terra::writeRaster(
      percent_raster,
      filename = file.path(output_dir, paste0("Percent_", basename(landcover_file))),
      overwrite = TRUE
    )
  }
}

load_current_predictors <- function(paths) {
  bioclim_rasters <- terra::rast(paths$current_bioclim_cropped)

  landcover_percent_files <- list.files(
    paths$landcover_dir,
    pattern = "Percent.*\\.tif$",
    full.names = TRUE
  )

  if (length(landcover_percent_files) == 0) {
    stop("No percentage land-cover .tif files were found.", call. = FALSE)
  }

  landcover_percent_rasters <- terra::rast(landcover_percent_files)
  c(bioclim_rasters, landcover_percent_rasters)
}

remove_african_records <- function(occurrence_records) {
  world_polygons <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  africa_polygons <- world_polygons[world_polygons$continent == "Africa", ]

  occurrence_points <- sf::st_as_sf(
    occurrence_records,
    coords = c("decimallongitude", "decimallatitude"),
    crs = 4326
  )

  africa_polygons <- sf::st_transform(africa_polygons, sf::st_crs(occurrence_points))
  africa_polygons <- sf::st_make_valid(africa_polygons)
  africa_buffer <- sf::st_buffer(africa_polygons, dist = 100000)

  overlap_list <- sf::st_within(occurrence_points, africa_buffer, sparse = TRUE)
  records_in_africa <- lengths(overlap_list) > 0

  occurrence_records[!records_in_africa, ]
}

build_species_environment_tables <- function(occurrence_raw, target_species, predictor_rasters) {
  occurrence_environment_tables <- vector("list", length(target_species))

  for (species_index in seq_along(target_species)) {
    species_name <- target_species[species_index]
    message("Extracting predictors for ", species_name)

    species_records <- occurrence_raw[occurrence_raw$species == species_name, ]
    species_records <- species_records[species_records$year >= 1990, ]
    species_records <- species_records[species_records$basisofrecord %in% valid_basis_of_record, ]

    if (species_name == "Ciconia episcopus") {
      species_records <- remove_african_records(species_records)
    }

    occurrence_table <- species_records[, c(
      "species",
      "year",
      "month",
      "decimallongitude",
      "decimallatitude"
    )]

    occurrence_points <- terra::vect(
      occurrence_table,
      geom = c("decimallongitude", "decimallatitude"),
      crs = "epsg:4326",
      keepgeom = FALSE
    )

    extracted_predictors <- terra::extract(predictor_rasters, occurrence_points)

    occurrence_environment_tables[[species_index]] <- cbind(
      extracted_predictors[1],
      occurrence_table,
      extracted_predictors[-1]
    )
  }

  names(occurrence_environment_tables) <- target_species
  occurrence_environment_tables
}

calculate_buffer_width <- function(occurrence_points) {
  if (nrow(occurrence_points) < 2) {
    stop("At least two occurrence points are required to calculate a buffer width.", call. = FALSE)
  }

  distance_matrix <- as.matrix(terra::distance(occurrence_points))
  minimum_distances <- numeric(0)

  for (point_index in seq_len(ncol(distance_matrix))) {
    positive_distances <- distance_matrix[distance_matrix[, point_index] > 0, point_index]
    minimum_distances <- c(minimum_distances, min(positive_distances))
  }

  2 * stats::quantile(minimum_distances, 0.95, na.rm = TRUE)
}

build_species_buffers <- function(occurrence_environment_tables, output_dir) {
  complete_occurrence_tables <- vector("list", length(occurrence_environment_tables))

  for (species_index in seq_along(occurrence_environment_tables)) {
    species_name <- names(occurrence_environment_tables)[species_index]
    species_id <- make_species_id(species_index)
    message("Building sampling buffer for ", species_name)

    species_table <- occurrence_environment_tables[[species_index]]
    complete_species_table <- species_table[stats::complete.cases(species_table), ]
    complete_occurrence_tables[[species_index]] <- complete_species_table

    occurrence_points <- terra::vect(
      complete_species_table,
      geom = c("decimallongitude", "decimallatitude"),
      crs = "epsg:4326",
      keepgeom = FALSE
    )

    buffer_width_m <- calculate_buffer_width(occurrence_points)
    sampling_buffer <- terra::buffer(occurrence_points, buffer_width_m)
    sampling_buffer <- terra::aggregate(sampling_buffer)
    sampling_buffer <- terra::project(sampling_buffer, terra::crs(occurrence_points))

    terra::writeVector(
      sampling_buffer,
      filename = file.path(output_dir, paste0("buffer_", species_id, ".shp")),
      filetype = "ESRI Shapefile",
      overwrite = TRUE
    )

    message("Buffer width for ", species_name, ": ", round(buffer_width_m, 2), " m")
  }

  names(complete_occurrence_tables) <- names(occurrence_environment_tables)
  complete_occurrence_tables
}

sample_background_points <- function(target_species, predictor_rasters, buffer_dir, output_dir) {
  buffer_files <- list.files(buffer_dir, pattern = "\\.shp$", full.names = TRUE)
  buffer_files <- buffer_files[order(buffer_files)]

  if (length(buffer_files) != length(target_species)) {
    stop("The number of buffer files does not match the number of target species.", call. = FALSE)
  }

  megaSDM::BackgroundPoints(
    spplist = target_species,
    envdata = predictor_rasters,
    output = output_dir,
    spatial_weights = spatial_weight,
    nbg = background_point_count,
    buffers = buffer_files,
    method = "Varela",
    ncores = 1
  )
}

varela_environmental_filter <- function(environment_occurrences, bin_count, use_pca = TRUE, pca_axes = "auto") {
  clean_occurrences <- environment_occurrences[stats::complete.cases(environment_occurrences), ]

  if (use_pca) {
    predictor_matrix <- clean_occurrences[, 3:ncol(clean_occurrences), drop = FALSE]
    pca_environment <- stats::prcomp(predictor_matrix, scale. = TRUE)
    pca_importance <- summary(pca_environment)$importance

    if (is.numeric(pca_axes)) {
      number_axes <- pca_axes
    } else {
      candidate_axes <- which(pca_importance[3, ] > 0.95)
      number_axes <- if (length(candidate_axes) == 0) 2 else max(2, min(candidate_axes))
    }

    environment_occurrences <- data.frame(
      clean_occurrences[, 1:2, drop = FALSE],
      pca_environment$x[, 1:number_axes, drop = FALSE]
    )
  } else {
    environment_occurrences <- clean_occurrences
  }

  sample_counts <- numeric(0)

  for (bin_index in seq_along(bin_count)) {
    bin_membership <- environment_occurrences[, 1:2, drop = FALSE]

    for (predictor_index in 3:ncol(environment_occurrences)) {
      predictor_values <- environment_occurrences[, predictor_index]
      value_range <- range(predictor_values, na.rm = TRUE)
      bin_resolution <- (value_range[2] - value_range[1]) / bin_count[bin_index]
      scaled_values <- (predictor_values - value_range[1]) / bin_resolution
      predictor_bins <- ceiling(scaled_values)
      predictor_bins[predictor_bins == 0] <- 1
      names(predictor_bins) <- names(environment_occurrences)[predictor_index]
      bin_membership <- cbind(bin_membership, predictor_bins)
      names(bin_membership)[ncol(bin_membership)] <- names(environment_occurrences)[predictor_index]
    }

    unique_bins <- dplyr::distinct(bin_membership[, -(1:2), drop = FALSE])
    unique_bins$group_id <- seq_len(nrow(unique_bins))

    join_columns <- setdiff(names(unique_bins), "group_id")
    bin_membership <- suppressMessages(dplyr::left_join(bin_membership, unique_bins, by = join_columns))

    filtered_points <- data.frame(x = numeric(), y = numeric())

    for (group_index in seq_len(nrow(unique_bins))) {
      group_members <- bin_membership[bin_membership$group_id == group_index, 1:2, drop = FALSE]
      selected_member <- group_members[sample(seq_len(nrow(group_members)), 1), ]
      filtered_points <- rbind(filtered_points, selected_member)
    }

    filtered_points <- data.frame(x = filtered_points[, 1], y = filtered_points[, 2])
    filtered_points <- merge(filtered_points, clean_occurrences, by = c("x", "y"), all.x = TRUE)
    sample_counts <- c(sample_counts, nrow(filtered_points))
  }

  if (length(bin_count) == 1) {
    filtered_points
  } else {
    data.frame(sample_count = sample_counts, climate_bin_count = bin_count)
  }
}

build_presence_absence_tables <- function(complete_occurrence_tables, target_species, background_dir, output_dir) {
  for (species_index in seq_along(target_species)) {
    species_name <- target_species[species_index]
    species_id <- make_species_id(species_index)
    species_slug <- make_species_slug(species_name)
    message("Building presence-background table for ", species_name)

    presence_table <- complete_occurrence_tables[[species_name]]
    background_file <- file.path(background_dir, paste0(gsub(" ", "_", species_name), "_background.csv"))
    background_table <- utils::read.csv(background_file, fileEncoding = csv_encoding)

    names(presence_table)[5:6] <- c("longitude", "latitude")
    names(background_table)[2:3] <- c("longitude", "latitude")

    presence_predictors <- presence_table[, -(1:6), drop = FALSE]
    predictor_variance <- vapply(presence_predictors, stats::var, numeric(1), na.rm = TRUE)
    usable_predictors <- names(predictor_variance)[predictor_variance != 0]

    set.seed(base_seed + species_index)
    filtered_presence_index <- varela_environmental_filter(
      environment_occurrences = data.frame(
        x = seq_len(nrow(presence_predictors)),
        y = seq_len(nrow(presence_predictors)),
        presence_predictors[, usable_predictors, drop = FALSE]
      ),
      bin_count = 25,
      use_pca = TRUE,
      pca_axes = "auto"
    )

    filtered_presence <- presence_table[filtered_presence_index$x, ]

    presence_data <- filtered_presence[, c("longitude", "latitude", usable_predictors), drop = FALSE]
    background_data <- background_table[, c("longitude", "latitude", usable_predictors), drop = FALSE]

    presence_data$occurrence <- 1
    background_data$occurrence <- 0

    set.seed(base_seed * species_index + 919)
    selected_background_rows <- sample(seq_len(nrow(background_data)), size = min(selected_background_count, nrow(background_data)))

    species_database <- rbind(presence_data, background_data[selected_background_rows, ])
    species_database$species <- unique(filtered_presence$species)

    utils::write.csv(
      species_database,
      file = file.path(output_dir, paste0("pa_predictor_", species_id, "_", species_slug, ".csv")),
      row.names = FALSE,
      fileEncoding = csv_encoding
    )

    message("Presence-background counts for ", species_name, ":")
    print(table(species_database$occurrence))
  }
}

split_sdm_databases <- function(presence_absence_dir, output_dir) {
  presence_absence_files <- list.files(
    presence_absence_dir,
    pattern = "^pa_predictor_sp[0-9]{2}_.*\\.csv$",
    full.names = TRUE
  )
  presence_absence_files <- presence_absence_files[order(presence_absence_files)]

  for (presence_absence_file in presence_absence_files) {
    message("Splitting SDM database: ", basename(presence_absence_file))

    species_id <- sub("^pa_predictor_(sp[0-9]{2})_.*", "\\1", basename(presence_absence_file))
    species_database <- utils::read.csv(presence_absence_file, fileEncoding = csv_encoding)
    pa_data_original <- species_database
    species_name <- unique(pa_data_original$species)

    predictor_columns <- setdiff(names(species_database), c("longitude", "latitude", "occurrence", "species"))
    predictor_data <- species_database[, predictor_columns, drop = FALSE]

    vif_correlation <- usdm::vifcor(predictor_data, th = 0.90)
    predictors_after_correlation <- usdm::exclude(predictor_data, vif_correlation)

    vif_stepwise <- usdm::vifstep(predictors_after_correlation, th = 5)
    predictors_after_vif <- usdm::exclude(predictors_after_correlation, vif_stepwise)

    selected_predictors <- names(predictors_after_vif)
    sdm_formula <- stats::as.formula(paste("occurrence ~", paste(selected_predictors, collapse = " + ")))

    modeling_database <- cbind(
      species_database["occurrence"],
      species_database[, selected_predictors, drop = FALSE]
    )

    presence_rows <- modeling_database[modeling_database$occurrence == 1, , drop = FALSE]
    background_rows <- modeling_database[modeling_database$occurrence == 0, , drop = FALSE]

    train_sets <- vector("list", repeat_count)
    test_sets <- vector("list", repeat_count)

    for (repeat_id in seq_len(repeat_count)) {
      set.seed(base_seed + repeat_id^2)
      presence_test_count <- max(1, round(test_fraction * nrow(presence_rows)))
      background_test_count <- max(1, round(test_fraction * nrow(background_rows)))

      presence_test_rows <- sample(seq_len(nrow(presence_rows)), size = presence_test_count)
      background_test_rows <- sample(seq_len(nrow(background_rows)), size = background_test_count)

      test_sets[[repeat_id]] <- rbind(
        presence_rows[presence_test_rows, , drop = FALSE],
        background_rows[background_test_rows, , drop = FALSE]
      )

      train_sets[[repeat_id]] <- rbind(
        presence_rows[-presence_test_rows, , drop = FALSE],
        background_rows[-background_test_rows, , drop = FALSE]
      )
    }

    save(
      species_name,
      species_id,
      selected_predictors,
      sdm_formula,
      pa_data_original,
      train_sets,
      test_sets,
      file = file.path(output_dir, paste0("db_sdm_k", repeat_count, "_", species_id, ".RData"))
    )
  }
}

if (run_landcover_percent_conversion) {
  reference_bioclim <- terra::rast(paths$current_bioclim_cropped)
  create_landcover_percent_rasters(paths$landcover_dir, reference_bioclim, paths$landcover_dir)
}

predictor_rasters <- load_current_predictors(paths)

if (run_presence_extraction) {
  loaded_objects <- load(paths$gbif_rdata)
  if (!"data" %in% loaded_objects) {
    stop("The GBIF RData file must contain an object named 'data'.", call. = FALSE)
  }

  occurrence_raw <- get("data")
  occurrence_environment_by_species <- build_species_environment_tables(
    occurrence_raw = occurrence_raw,
    target_species = target_species,
    predictor_rasters = predictor_rasters
  )

  save(
    target_species,
    occurrence_environment_by_species,
    file = file.path(paths$presence_absence_dir, "presence_environment_by_species.RData")
  )
}

if (run_buffer_generation) {
  load(file.path(paths$presence_absence_dir, "presence_environment_by_species.RData"))
  complete_occurrence_by_species <- build_species_buffers(
    occurrence_environment_tables = occurrence_environment_by_species,
    output_dir = paths$buffer_dir
  )

  save(
    target_species,
    complete_occurrence_by_species,
    file = file.path(paths$presence_absence_dir, "presence_environment_complete.RData")
  )
}

if (run_background_sampling) {
  sample_background_points(
    target_species = target_species,
    predictor_rasters = predictor_rasters,
    buffer_dir = paths$buffer_dir,
    output_dir = paths$background_dir
  )
}

if (run_presence_absence_build) {
  load(file.path(paths$presence_absence_dir, "presence_environment_complete.RData"))
  build_presence_absence_tables(
    complete_occurrence_tables = complete_occurrence_by_species,
    target_species = target_species,
    background_dir = paths$background_dir,
    output_dir = paths$presence_absence_dir
  )
}

if (run_train_test_split) {
  split_sdm_databases(paths$presence_absence_dir, paths$sdm_database_dir)
}
