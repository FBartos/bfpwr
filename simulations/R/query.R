bfpwr_sim_query_fixed_power <- function(analysis, design, trajectories, n) {
    if (!n %in% trajectories$n) {
        stop("requested n is not present in trajectory grid: ", n)
    }
    subset <- trajectories[trajectories$n == n, , drop = FALSE]
    log_bf01 <- bfpwr_sim_log_bf(analysis, subset)
    pH1 <- mean(log_bf01 <= log(analysis$k1))
    pH0 <- mean(log_bf01 >= log(analysis$k0))
    ref <- bfpwr_sim_reference_fixed_power(analysis, design, n)
    data.frame(
        analysis_case_id = analysis$analysis_case_id,
        design_case_id = design$design_case_id,
        n = n,
        nsim = length(log_bf01),
        pH1 = pH1,
        pH0 = pH0,
        mcse_pH1 = sqrt(pH1 * (1 - pH1) / length(log_bf01)),
        mcse_pH0 = sqrt(pH0 * (1 - pH0) / length(log_bf01)),
        reference_pH1 = ref$reference_pH1,
        reference_pH0 = ref$reference_pH0
    )
}

bfpwr_sim_query_fixed_power_from_chunks <- function(analysis, design,
                                                    corpus_root, n) {
    if (!n %in% design$look_grid) {
        stop("requested n is not present in design look_grid: ", n)
    }
    chunk_files <- bfpwr_sim_design_chunk_files(corpus_root, design)
    missing <- chunk_files[!file.exists(chunk_files)]
    if (length(missing) > 0) {
        bfpwr_sim_stop_missing_corpus_files(
            corpus_root,
            missing = missing,
            asset_set = "chunk-validation",
            label = "design chunk files")
    }

    counts <- c(H1 = 0L, H0 = 0L, total = 0L)
    for (file in chunk_files) {
        trajectories <- readRDS(file)
        subset <- trajectories[trajectories$n == n, , drop = FALSE]
        log_bf01 <- bfpwr_sim_log_bf(analysis, subset)
        counts[["H1"]] <- counts[["H1"]] + sum(log_bf01 <= log(analysis$k1))
        counts[["H0"]] <- counts[["H0"]] + sum(log_bf01 >= log(analysis$k0))
        counts[["total"]] <- counts[["total"]] + length(log_bf01)
    }

    pH1 <- counts[["H1"]] / counts[["total"]]
    pH0 <- counts[["H0"]] / counts[["total"]]
    ref <- bfpwr_sim_reference_fixed_power(analysis, design, n)
    data.frame(
        analysis_case_id = analysis$analysis_case_id,
        design_case_id = design$design_case_id,
        n = n,
        nsim = counts[["total"]],
        pH1 = pH1,
        pH0 = pH0,
        mcse_pH1 = sqrt(pH1 * (1 - pH1) / counts[["total"]]),
        mcse_pH0 = sqrt(pH0 * (1 - pH0) / counts[["total"]]),
        reference_pH1 = ref$reference_pH1,
        reference_pH0 = ref$reference_pH0
    )
}
