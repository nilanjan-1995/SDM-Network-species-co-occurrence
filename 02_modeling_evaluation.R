################################################################################
# 02_modeling_evaluation.R
#
# Purpose:
#   Fit Random Forest and MaxEnt SDMs, tune model parameters with Bayesian
#   optimization, evaluate train/test performance, and save model outputs.
################################################################################

project_root <- normalizePath(".", mustWork = FALSE)
setwd(project_root)

paths <- list(
  sdm_database_dir = file.path("outputs", "04_sdm_database"),
  model_output_dir = file.path("outputs", "05_sdm_results")
)

csv_encoding <- "UTF-8"
model_methods <- c("RF", "MaxEnt")
base_seed <- 910520
cv_fold_count <- 5
bayes_init_points <- 5
bayes_iteration_count <- 5
bayes_exploitation_count <- 5
use_parallel_optimization <- FALSE
detected_core_count <- parallel::detectCores()
parallel_core_count <- if (is.na(detected_core_count)) 1 else max(1, detected_core_count - 3)

maxent_param_bounds <- list(
  regmult = c(0.1, 2.0)
)

rf_param_bounds <- list(
  num.trees = c(500L, 1500L),
  mtry = c(3L, 6L),
  min.node.size = c(3L, 7L)
)

required_packages <- c(
  "maxnet",
  "ranger",
  "ParBayesianOptimization",
  "dismo",
  "caret",
  "MLmetrics",
  "doParallel"
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
dir.create(paths$model_output_dir, recursive = TRUE, showWarnings = FALSE)

safe_metric <- function(metric_expression) {
  value <- tryCatch(metric_expression, error = function(error) NA_real_)
  ifelse(is.nan(value), NA_real_, value)
}

scale_score_vector <- function(score_values) {
  if (all(is.na(score_values)) || length(unique(score_values[!is.na(score_values)])) <= 1) {
    return(rep(0, length(score_values)))
  }
  (score_values - min(score_values, na.rm = TRUE)) /
    (max(score_values, na.rm = TRUE) - min(score_values, na.rm = TRUE))
}

calculate_optimization_scores <- function(actual, probability, predicted_class) {
  sensitivity <- safe_metric(MLmetrics::Sensitivity(y_true = actual, y_pred = predicted_class, positive = 1))
  specificity <- safe_metric(MLmetrics::Specificity(y_true = actual, y_pred = predicted_class, positive = 1))
  tss <- sensitivity + specificity - 1
  pr_auc <- safe_metric(MLmetrics::PRAUC(y_pred = probability, y_true = actual))

  data.frame(
    tss = ifelse(is.na(tss), 0, tss),
    pr_auc = ifelse(is.na(pr_auc), 0, pr_auc)
  )
}

calculate_mcc <- function(actual, predicted_class) {
  true_positive <- as.numeric(sum(predicted_class == 1 & actual == 1))
  true_negative <- as.numeric(sum(predicted_class == 0 & actual == 0))
  false_positive <- as.numeric(sum(predicted_class == 1 & actual == 0))
  false_negative <- as.numeric(sum(predicted_class == 0 & actual == 1))

  numerator <- (true_positive * true_negative) - (false_positive * false_negative)
  denominator <- sqrt(
    (true_positive + false_positive) *
      (true_positive + false_negative) *
      (true_negative + false_positive) *
      (true_negative + false_negative)
  )

  if (is.na(denominator) || denominator == 0) {
    return(NA_real_)
  }

  numerator / denominator
}

calculate_evaluation_metrics <- function(actual, probability, predicted_class) {
  sensitivity <- safe_metric(MLmetrics::Sensitivity(y_true = actual, y_pred = predicted_class, positive = 1))
  specificity <- safe_metric(MLmetrics::Specificity(y_true = actual, y_pred = predicted_class, positive = 1))

  data.frame(
    accuracy = safe_metric(MLmetrics::Accuracy(y_pred = predicted_class, y_true = actual)),
    auc = safe_metric(MLmetrics::AUC(y_pred = probability, y_true = actual)),
    pr_auc = safe_metric(MLmetrics::PRAUC(y_pred = probability, y_true = actual)),
    tss = sensitivity + specificity - 1,
    f1_presence = safe_metric(MLmetrics::F1_Score(y_true = actual, y_pred = predicted_class, positive = 1)),
    sensitivity_presence = sensitivity,
    specificity_presence = specificity,
    precision_presence = safe_metric(MLmetrics::Precision(y_true = actual, y_pred = predicted_class, positive = 1)),
    recall_presence = safe_metric(MLmetrics::Recall(y_true = actual, y_pred = predicted_class, positive = 1)),
    mcc = calculate_mcc(actual = actual, predicted_class = predicted_class)
  )
}

fit_maxent_model <- function(training_data, selected_predictors, regmult) {
  occurrence_vector <- training_data$occurrence
  predictor_data <- training_data[, selected_predictors, drop = FALSE]

  maxnet::maxnet(
    p = occurrence_vector,
    data = predictor_data,
    f = maxnet::maxnet.formula(occurrence_vector, predictor_data, classes = "default"),
    regmult = regmult,
    addsamplestobackground = TRUE
  )
}

predict_maxent_model <- function(model, new_data) {
  predict(model, new_data, type = "cloglog")
}

fit_rf_model <- function(training_data, sdm_formula, selected_predictors, num_trees, mtry, min_node_size) {
  mtry <- max(1, min(length(selected_predictors), as.integer(round(mtry))))

  ranger::ranger(
    formula = sdm_formula,
    data = training_data,
    num.trees = as.integer(round(num_trees)),
    mtry = mtry,
    min.node.size = as.integer(round(min_node_size))
  )
}

predict_rf_model <- function(model, new_data) {
  predict(model, new_data)$predictions
}

calculate_spec_sens_threshold <- function(probability, occurrence) {
  evaluation_object <- dismo::evaluate(
    p = probability[occurrence == 1],
    a = probability[occurrence == 0]
  )
  dismo::threshold(evaluation_object, "spec_sens")
}

run_bayes_optimization <- function(optimization_function, parameter_bounds, worker_packages = character(0)) {
  optimization_result <- NULL
  cluster_object <- NULL

  if (use_parallel_optimization) {
    cluster_object <- parallel::makeCluster(parallel_core_count)
    doParallel::registerDoParallel(cluster_object)
    parallel::clusterCall(
      cluster_object,
      function(package_names) invisible(lapply(package_names, library, character.only = TRUE)),
      c("ParBayesianOptimization", worker_packages)
    )
  }

  optimization_result <- tryCatch(
    ParBayesianOptimization::bayesOpt(
      FUN = optimization_function,
      bounds = parameter_bounds,
      initPoints = bayes_init_points,
      iters.n = bayes_iteration_count,
      iters.k = bayes_exploitation_count,
      parallel = use_parallel_optimization
    ),
    error = function(error) {
      message("Bayesian optimization failed. Default parameters will be used. Message: ", error$message)
      NULL
    }
  )

  if (!is.null(cluster_object)) {
    parallel::stopCluster(cluster_object)
    doParallel::registerDoSEQ()
  }

  optimization_result
}

extract_best_parameters <- function(optimization_result, model_method, selected_predictors) {
  if (!is.null(optimization_result)) {
    best_parameters <- as.data.frame(as.list(ParBayesianOptimization::getBestPars(optimization_result)))
    best_parameters$parameter_status <- NA_character_
    return(best_parameters)
  }

  if (model_method == "MaxEnt") {
    return(data.frame(regmult = 1.0, parameter_status = "default_after_error"))
  }

  data.frame(
    num.trees = 500L,
    mtry = max(1, round(sqrt(length(selected_predictors)))),
    min.node.size = 5L,
    parameter_status = "default_after_error"
  )
}

sdm_database_files <- list.files(
  paths$sdm_database_dir,
  pattern = "^db_sdm_k[0-9]+_sp[0-9]{2}\\.RData$",
  full.names = TRUE
)
sdm_database_files <- sdm_database_files[order(sdm_database_files)]

for (database_index in seq_along(sdm_database_files)) {
  species_timer <- proc.time()
  load(sdm_database_files[database_index])

  message("Processing ", species_id, " - ", species_name)

  evaluation_by_method <- list()
  prediction_by_method <- list()

  for (method_index in seq_along(model_methods)) {
    current_method <- model_methods[method_index]
    method_timer <- proc.time()
    message("Fitting method: ", current_method)

    sdm_models <- vector("list", length(test_sets))
    evaluation_table <- NULL
    prediction_table <- NULL
    parameter_table <- NULL

    for (repeat_id in seq_along(test_sets)) {
      training_data <- train_sets[[repeat_id]]
      testing_data <- test_sets[[repeat_id]]

      set.seed(base_seed * cv_fold_count + database_index)
      folds <- caret::createFolds(training_data$occurrence, k = cv_fold_count, returnTrain = FALSE)

      maxent_optimization_function <- function(regmult) {
        fold_tss <- numeric(0)
        fold_pr_auc <- numeric(0)

        for (fold_id in seq_len(cv_fold_count)) {
          training_fold <- training_data[-folds[[fold_id]], , drop = FALSE]
          validation_fold <- training_data[folds[[fold_id]], , drop = FALSE]

          set.seed(base_seed + repeat_id * method_index * fold_id)
          candidate_model <- fit_maxent_model(training_fold, selected_predictors, regmult)
          validation_probability <- predict_maxent_model(candidate_model, validation_fold)
          validation_threshold <- calculate_spec_sens_threshold(
            probability = validation_probability,
            occurrence = validation_fold$occurrence
          )
          validation_class <- ifelse(validation_probability >= validation_threshold, 1, 0)

          fold_scores <- calculate_optimization_scores(
            actual = validation_fold$occurrence,
            probability = validation_probability,
            predicted_class = validation_class
          )

          fold_tss <- c(fold_tss, fold_scores$tss)
          fold_pr_auc <- c(fold_pr_auc, fold_scores$pr_auc)
        }

        combined_score <- scale_score_vector(fold_tss) + scale_score_vector(fold_pr_auc)
        list(Score = mean(combined_score, na.rm = TRUE))
      }

      rf_optimization_function <- function(num.trees, mtry, min.node.size) {
        fold_tss <- numeric(0)
        fold_pr_auc <- numeric(0)

        for (fold_id in seq_len(cv_fold_count)) {
          training_fold <- training_data[-folds[[fold_id]], , drop = FALSE]
          validation_fold <- training_data[folds[[fold_id]], , drop = FALSE]

          set.seed(base_seed + repeat_id * method_index * fold_id)
          candidate_model <- fit_rf_model(
            training_data = training_fold,
            sdm_formula = sdm_formula,
            selected_predictors = selected_predictors,
            num_trees = num.trees,
            mtry = mtry,
            min_node_size = min.node.size
          )

          validation_probability <- predict_rf_model(candidate_model, validation_fold)
          validation_threshold <- calculate_spec_sens_threshold(
            probability = validation_probability,
            occurrence = validation_fold$occurrence
          )
          validation_class <- ifelse(validation_probability >= validation_threshold, 1, 0)

          fold_scores <- calculate_optimization_scores(
            actual = validation_fold$occurrence,
            probability = validation_probability,
            predicted_class = validation_class
          )

          fold_tss <- c(fold_tss, fold_scores$tss)
          fold_pr_auc <- c(fold_pr_auc, fold_scores$pr_auc)
        }

        combined_score <- scale_score_vector(fold_tss) + scale_score_vector(fold_pr_auc)
        list(Score = mean(combined_score, na.rm = TRUE))
      }

      if (current_method == "MaxEnt") {
        optimization_result <- run_bayes_optimization(
          optimization_function = maxent_optimization_function,
          parameter_bounds = maxent_param_bounds,
          worker_packages = c("maxnet", "dismo", "MLmetrics")
        )
      } else {
        optimization_result <- run_bayes_optimization(
          optimization_function = rf_optimization_function,
          parameter_bounds = rf_param_bounds,
          worker_packages = c("ranger", "dismo", "MLmetrics")
        )
      }

      best_parameters <- extract_best_parameters(
        optimization_result = optimization_result,
        model_method = current_method,
        selected_predictors = selected_predictors
      )

      if (current_method == "MaxEnt") {
        set.seed(base_seed + repeat_id * method_index)
        final_model <- NULL
        maxent_attempt <- 0

        while (is.null(final_model) && maxent_attempt < 20) {
          final_model <- tryCatch(
            fit_maxent_model(training_data, selected_predictors, best_parameters$regmult),
            error = function(error) NULL
          )

          if (is.null(final_model)) {
            best_parameters$regmult <- best_parameters$regmult + 0.1
            maxent_attempt <- maxent_attempt + 1
          }
        }

        if (is.null(final_model)) {
          stop("MaxEnt failed after repeated attempts for ", species_id, " repeat ", repeat_id, call. = FALSE)
        }

        train_probability <- predict_maxent_model(final_model, training_data)
        test_probability <- predict_maxent_model(final_model, testing_data)
      } else {
        set.seed(base_seed + repeat_id * method_index)
        final_model <- fit_rf_model(
          training_data = training_data,
          sdm_formula = sdm_formula,
          selected_predictors = selected_predictors,
          num_trees = best_parameters$num.trees,
          mtry = best_parameters$mtry,
          min_node_size = best_parameters$min.node.size
        )

        train_probability <- predict_rf_model(final_model, training_data)
        test_probability <- predict_rf_model(final_model, testing_data)
      }

      sdm_models[[repeat_id]] <- final_model

      train_threshold <- calculate_spec_sens_threshold(
        probability = train_probability,
        occurrence = training_data$occurrence
      )
      test_threshold <- calculate_spec_sens_threshold(
        probability = test_probability,
        occurrence = testing_data$occurrence
      )

      train_class <- ifelse(train_probability >= test_threshold, 1, 0)
      test_class <- ifelse(test_probability >= test_threshold, 1, 0)

      train_metrics <- calculate_evaluation_metrics(
        actual = training_data$occurrence,
        probability = train_probability,
        predicted_class = train_class
      )
      names(train_metrics) <- paste0("train_", names(train_metrics))

      test_metrics <- calculate_evaluation_metrics(
        actual = testing_data$occurrence,
        probability = test_probability,
        predicted_class = test_class
      )
      names(test_metrics) <- paste0("test_", names(test_metrics))

      repeat_evaluation <- cbind(
        data.frame(
          species = species_name,
          species_id = species_id,
          model_method = current_method,
          method_index = method_index,
          repeat_id = repeat_id,
          parameter_status = best_parameters$parameter_status,
          n_train = nrow(training_data),
          n_train_background = sum(training_data$occurrence == 0),
          n_train_presence = sum(training_data$occurrence == 1),
          n_test = nrow(testing_data),
          n_test_background = sum(testing_data$occurrence == 0),
          n_test_presence = sum(testing_data$occurrence == 1),
          threshold_train = train_threshold,
          threshold_test = test_threshold
        ),
        train_metrics,
        test_metrics
      )

      repeat_parameters <- cbind(
        data.frame(
          species = species_name,
          species_id = species_id,
          model_method = current_method,
          method_index = method_index,
          repeat_id = repeat_id
        ),
        best_parameters
      )

      prediction_inputs <- rbind(
        cbind(data_split = "train", training_data[, c("occurrence", selected_predictors), drop = FALSE]),
        cbind(data_split = "test", testing_data[, c("occurrence", selected_predictors), drop = FALSE])
      )

      repeat_predictions <- cbind(
        data.frame(
          species = species_name,
          species_id = species_id,
          model_method = current_method,
          method_index = method_index,
          repeat_id = repeat_id,
          threshold_train = train_threshold,
          threshold_test = test_threshold
        ),
        prediction_inputs
      )

      repeat_predictions$prediction_probability <- c(train_probability, test_probability)
      repeat_predictions$prediction_class <- c(train_class, test_class)

      evaluation_table <- rbind(evaluation_table, repeat_evaluation)
      prediction_table <- rbind(prediction_table, repeat_predictions)
      parameter_table <- rbind(parameter_table, repeat_parameters)

      message(
        database_index, "/", length(sdm_database_files), " files, ",
        method_index, "/", length(model_methods), " methods, repeat ",
        repeat_id, "/", length(test_sets), "; elapsed method time: ",
        round((proc.time() - method_timer)[3], 1), " sec"
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
      sdm_models,
      evaluation_table,
      prediction_table,
      parameter_table,
      file = file.path(paths$model_output_dir, paste0("sdm_output_", species_id, "_", current_method, ".RData"))
    )

    evaluation_by_method[[current_method]] <- evaluation_table
    prediction_by_method[[current_method]] <- prediction_table
  }

  evaluation_all_methods <- do.call(rbind, evaluation_by_method)
  prediction_all_methods <- do.call(rbind, prediction_by_method)

  utils::write.csv(
    evaluation_all_methods,
    file = file.path(paths$model_output_dir, paste0("evaluation_", species_id, ".csv")),
    row.names = FALSE,
    fileEncoding = csv_encoding
  )

  utils::write.csv(
    prediction_all_methods,
    file = file.path(paths$model_output_dir, paste0("predictions_", species_id, ".csv")),
    row.names = FALSE,
    fileEncoding = csv_encoding
  )

  message(
    database_index, "/", length(sdm_database_files), " complete; elapsed species time: ",
    round((proc.time() - species_timer)[3], 1), " sec"
  )
}
