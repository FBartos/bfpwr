bfpwr_sim_package_function <- function(name) {
    fn <- get(name, mode = "function", inherits = TRUE)
    if (!is.function(fn)) {
        stop("required bfpwr function is not available: ", name)
    }
    fn
}

bfpwr_sim_validate_bf_evaluation_case <- function(case) {
    if (!is.null(case$bf_prior_id)) {
        return(bfpwr_sim_validate_bf_prior_case(case))
    }
    bfpwr_sim_validate_analysis_case(case)
}

bfpwr_sim_bf_evaluation_id <- function(case) {
    if (!is.null(case$bf_prior_id)) {
        return(case$bf_prior_id)
    }
    case$analysis_case_id
}

bfpwr_sim_log_bf <- function(analysis, trajectories) {
    bfpwr_sim_validate_bf_evaluation_case(analysis)
    prior <- analysis$analysis_prior

    if (analysis$test_family == "z") {
        if (analysis$bf_type == "normal") {
            return(bfpwr_sim_package_function("bf01")(
                estimate = trajectories$estimate,
                se = trajectories$se,
                null = prior$null,
                pm = prior$pm,
                psd = prior$psd,
                log = TRUE
            ))
        }
        if (analysis$bf_type == "directional") {
            return(bfpwr_sim_package_function("dirbf01")(
                estimate = trajectories$estimate,
                se = trajectories$se,
                null = prior$null,
                pm = prior$pm,
                psd = prior$psd,
                log = TRUE
            ))
        }
        if (analysis$bf_type == "moment") {
            return(bfpwr_sim_package_function("nmbf01")(
                estimate = trajectories$estimate,
                se = trajectories$se,
                null = prior$null,
                psd = prior$psd,
                log = TRUE
            ))
        }
    }

    if (analysis$test_family == "t") {
        neff <- if (prior$type == "two.sample") {
            1 / (1 / trajectories$n1 + 1 / trajectories$n2)
        } else {
            trajectories$n1
        }
        t_eval <- (trajectories$estimate - prior$null) * sqrt(neff)
        return(bfpwr_sim_package_function("tbf01")(
            t = t_eval,
            n = trajectories$n1,
            n1 = trajectories$n1,
            n2 = trajectories$n2,
            plocation = prior$plocation - prior$null,
            pscale = prior$pscale,
            pdf = prior$pdf,
            type = prior$type,
            alternative = prior$alternative,
            log = TRUE
        ))
    }

    if (analysis$test_family == "binomial") {
        return(bfpwr_sim_package_function("binbf01")(
            x = trajectories$x,
            n = trajectories$n,
            p0 = prior$p0,
            type = analysis$bf_type,
            a = prior$a,
            b = prior$b,
            log = TRUE
        ))
    }

    stop("unsupported analysis family/type")
}

bfpwr_sim_bf_rows <- function(analysis, trajectories) {
    log_bf01 <- bfpwr_sim_log_bf(analysis, trajectories)
    bf_rows <- trajectories[c("design_case_id", "replicate_id", "chunk_id",
                              "look", "n")]
    if ("n1" %in% names(trajectories)) bf_rows$n1 <- trajectories$n1
    if ("n2" %in% names(trajectories)) bf_rows$n2 <- trajectories$n2
    bf_rows$bf_prior_id <- bfpwr_sim_bf_evaluation_id(analysis)
    if (!is.null(analysis$analysis_case_id)) {
        bf_rows$analysis_case_id <- analysis$analysis_case_id
    }
    bf_rows$log_bf01 <- log_bf01
    bf_rows
}

bfpwr_sim_materialize_analysis <- function(analysis, design, trajectories,
                                           include_bf_trajectories = FALSE) {
    if (!identical(analysis$design_case_id, design$design_case_id)) {
        stop("analysis case does not point to supplied design case")
    }
    if (!identical(analysis$test_family, design$test_family)) {
        stop("analysis and design families differ")
    }

    bf_rows <- bfpwr_sim_bf_rows(analysis, trajectories)
    outcomes <- bfpwr_sim_stop_outcomes(analysis, bf_rows)
    summary <- bfpwr_sim_summarize_analysis(analysis, design, outcomes)

    out <- list(analysis = analysis, outcomes = outcomes, summary = summary)
    if (include_bf_trajectories) {
        out$bf_trajectories <- bf_rows
    }
    out
}

bfpwr_sim_materialize_analysis_from_chunks <- function(analysis, design,
                                                       corpus_root,
                                                       include_bf_trajectories = FALSE,
                                                       bf_chunk_dir = NULL) {
    if (!identical(analysis$design_case_id, design$design_case_id)) {
        stop("analysis case does not point to supplied design case")
    }
    if (!identical(analysis$test_family, design$test_family)) {
        stop("analysis and design families differ")
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
    if (!is.null(bf_chunk_dir)) {
        dir.create(bf_chunk_dir, recursive = TRUE, showWarnings = FALSE)
    }

    outcomes <- vector("list", length(chunk_files))
    bf_chunks <- if (include_bf_trajectories && is.null(bf_chunk_dir)) {
        vector("list", length(chunk_files))
    } else {
        NULL
    }
    bf_chunk_files <- character()

    for (i in seq_along(chunk_files)) {
        trajectories <- readRDS(chunk_files[[i]])
        bfpwr_sim_validate_design_trajectories(design, trajectories)
        bf_rows <- bfpwr_sim_bf_rows(analysis, trajectories)
        outcomes[[i]] <- bfpwr_sim_stop_outcomes(analysis, bf_rows)

        if (include_bf_trajectories) {
            if (is.null(bf_chunk_dir)) {
                bf_chunks[[i]] <- bf_rows
            } else {
                out_file <- file.path(bf_chunk_dir,
                                      sprintf("chunk-%04d.rds",
                                              design$chunks$chunk_id[[i]]))
                bfpwr_sim_write_rds(bf_rows, out_file)
                bf_chunk_files <- c(bf_chunk_files, out_file)
            }
        }
    }

    outcomes <- do.call(rbind, outcomes)
    summary <- bfpwr_sim_summarize_analysis(analysis, design, outcomes)
    out <- list(analysis = analysis, outcomes = outcomes, summary = summary)
    if (include_bf_trajectories) {
        if (is.null(bf_chunk_dir)) {
            out$bf_trajectories <- do.call(rbind, bf_chunks)
        } else {
            out$bf_chunk_files <- bf_chunk_files
        }
    }
    out
}

bfpwr_sim_stop_outcomes <- function(analysis, bf_rows) {
    log_k1 <- log(analysis$k1)
    log_k0 <- log(analysis$k0)
    bf_rows <- bf_rows[order(bf_rows$replicate_id, bf_rows$look), ]
    split_rows <- split(bf_rows, bf_rows$replicate_id)

    rows <- lapply(split_rows, function(x) {
        hit_h0 <- x$log_bf01 >= log_k0
        hit_h1 <- x$log_bf01 <= log_k1
        hit <- which(hit_h0 | hit_h1)
        if (length(hit) == 0) {
            j <- nrow(x)
            decision <- "inconclusive"
        } else {
            j <- hit[1]
            decision <- if (hit_h0[j]) "H0" else "H1"
        }
        data.frame(
            analysis_case_id = analysis$analysis_case_id,
            design_case_id = analysis$design_case_id,
            replicate_id = x$replicate_id[j],
            stop_look = x$look[j],
            stop_n = x$n[j],
            decision = decision,
            final_log_bf01 = x$log_bf01[j]
        )
    })
    do.call(rbind, rows)
}
