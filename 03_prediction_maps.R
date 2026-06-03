################################################################################
# 03_prediction_maps.R
#
# Purpose:
#   Use selected SDM repeats to create current and future prediction maps, then
#   create species-level ensemble maps across available SDM methods.
################################################################################

project_root <- normalizePath(".", mustWork = FALSE)
setwd(project_root)

paths <- list(
  model_output_dir = file.path("outputs", "05_sdm_results"),
  map_output_dir = file.path("outputs", "06_sdm_maps"),
  landcover_dir = file.path("data", "landcover"),
  current_bioclim_cropped = file.path("data", "worldclim_current", "cropped_bc_30ys_1991_2020.tif"),
  future_bioclim_dir = file.path("data", "worldclim_future")
)

csv_encoding <- "UTF-8"
raster_core_count <- 5
minimum_test_tss <- 0.50
minimum_test_auc <- 0.70

required_packages <- c("ranger", "maxnet", "terra")

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
dir.create(paths$map_output_dir, recursive = TRUE, showWarnings = FALSE)

make_clean_id <- function(text_value) {
  gsub("_+$", "", gsub("[^A-Za-z0-9]+", "_", text_value))
}

load_landcover_rasters <- function(landcover_dir) {
  landcover_files <- list.files(landcover_dir, pattern = "Percent.*\\.tif$", full.names = TRUE)
  landcover_files <- landcover_files[order(landcover_files)]

  if (length(landcover_files) == 0) {
    stop("No percentage land-cover .tif files were found.", call. = FALSE)
  }

  terra::rast(landcover_files)
}

maxent_raster_predict <- function(model, newdata, ...) {
  library(maxnet)
  predict(model, newdata, ..., type = "cloglog")
}

rf_raster_predict <- function(model, newdata, ...) {
  library(ranger)
  predict(model, newdata, ...)$predictions
}

mean_raster_stack <- function(raster_stack) {
  terra::app(raster_stack, fun = mean, na.rm = TRUE)
}

model_files <- list.files(
  paths$model_output_dir,
  pattern = "^sdm_output_sp[0-9]{2}_(RF|MaxEnt)\\.RData$",
  full.names = TRUE
)
model_files <- model_files[order(model_files)]

future_files <- list.files(paths$future_bioclim_dir, pattern = "\\.tif$", full.names = TRUE)
future_files <- future_files[order(future_files)]

scenario_table <- data.frame(
  scenario_id = "current_1991_2020",
  bioclim_file = NA_character_,
  is_current = TRUE,
  stringsAsFactors = FALSE
)

if (length(future_files) > 0) {
  scenario_table <- rbind(
    scenario_table,
    data.frame(
      scenario_id = vapply(tools::file_path_sans_ext(basename(future_files)), make_clean_id, character(1)),
      bioclim_file = future_files,
      is_current = FALSE,
      stringsAsFactors = FALSE
    )
  )
}

landcover_rasters <- load_landcover_rasters(paths$landcover_dir)
current_bioclim_names <- NULL
pipeline_timer <- proc.time()

for (scenario_index in seq_len(nrow(scenario_table))) {
  scenario_id <- scenario_table$scenario_id[scenario_index]
  message("Creating prediction maps for scenario: ", scenario_id)

  if (isTRUE(scenario_table$is_current[scenario_index])) {
    bioclim_rasters <- terra::rast(paths$current_bioclim_cropped)
    current_bioclim_names <- names(bioclim_rasters)
  } else {
    raw_future_rasters <- terra::rast(scenario_table$bioclim_file[scenario_index])
    bioclim_rasters <- terra::crop(raw_future_rasters, landcover_rasters)

    if (!is.null(current_bioclim_names) && length(current_bioclim_names) == terra::nlyr(bioclim_rasters)) {
      names(bioclim_rasters) <- current_bioclim_names
    }
  }

  predictor_rasters <- c(bioclim_rasters, landcover_rasters)
  model_count_records <- list()
  record_index <- 0

  for (model_index in seq_along(model_files)) {
    load(model_files[model_index])

    model_method <- unique(evaluation_table$model_method)
    selected_repeat_ids <- evaluation_table$repeat_id[
      !is.na(evaluation_table$test_tss) &
        evaluation_table$test_tss >= minimum_test_tss &
        evaluation_table$test_auc >= minimum_test_auc
    ]

    record_index <- record_index + 1
    model_count_records[[record_index]] <- data.frame(
      species_id = species_id,
      species = species_name,
      model_method = model_method,
      selected_model_count = length(selected_repeat_ids),
      selected_repeat_ids = paste(selected_repeat_ids, collapse = ","),
      stringsAsFactors = FALSE
    )

    if (length(selected_repeat_ids) == 0) {
      message("No selected repeats for ", species_id, " ", model_method)
      next
    }

    predicted_maps <- vector("list", length(selected_repeat_ids))

    for (repeat_position in seq_along(selected_repeat_ids)) {
      repeat_id <- selected_repeat_ids[repeat_position]
      message(
        "Predicting ", species_id, " ", model_method,
        " repeat ", repeat_id, " for ", scenario_id
      )

      if (model_method == "MaxEnt") {
        predicted_maps[[repeat_position]] <- terra::predict(
          predictor_rasters,
          sdm_models[[repeat_id]],
          fun = maxent_raster_predict,
          na.rm = TRUE,
          cores = raster_core_count
        )
      } else if (model_method == "RF") {
        predicted_maps[[repeat_position]] <- terra::predict(
          predictor_rasters,
          sdm_models[[repeat_id]],
          fun = rf_raster_predict,
          na.rm = TRUE,
          cores = raster_core_count
        )
      } else {
        stop("Unsupported model method: ", model_method, call. = FALSE)
      }
    }

    selected_map_stack <- terra::rast(predicted_maps)
    mean_prediction_map <- mean_raster_stack(selected_map_stack)

    terra::writeRaster(
      mean_prediction_map,
      filename = file.path(paths$map_output_dir, paste0("map_", species_id, "_", model_method, "_", scenario_id, ".tif")),
      overwrite = TRUE
    )

    rm(predicted_maps, selected_map_stack, mean_prediction_map)
    gc()
  }

  model_count_table <- do.call(rbind, model_count_records)
  utils::write.csv(
    model_count_table,
    file = file.path(paths$map_output_dir, paste0("model_counts_by_species_", scenario_id, ".csv")),
    row.names = FALSE,
    fileEncoding = csv_encoding
  )

  selected_species_ids <- sort(unique(model_count_table$species_id[model_count_table$selected_model_count > 0]))

  for (species_id_for_ensemble in selected_species_ids) {
    species_map_files <- list.files(
      paths$map_output_dir,
      pattern = paste0("^map_", species_id_for_ensemble, "_.*_", scenario_id, "\\.tif$"),
      full.names = TRUE
    )
    species_map_files <- species_map_files[order(species_map_files)]

    if (length(species_map_files) == 0) {
      next
    }

    species_map_stack <- terra::rast(species_map_files)
    ensemble_map <- mean_raster_stack(species_map_stack)

    terra::writeRaster(
      ensemble_map,
      filename = file.path(paths$map_output_dir, paste0("ensemble_", species_id_for_ensemble, "_", scenario_id, ".tif")),
      overwrite = TRUE
    )

    rm(species_map_stack, ensemble_map)
    gc()
  }

  message("Scenario complete: ", scenario_id, "; elapsed time: ", round((proc.time() - pipeline_timer)[3], 1), " sec")
}
