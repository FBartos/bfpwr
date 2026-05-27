bfpwr_sim_bf_prior_split_chunk_file <- function(corpus_root, bf_prior,
                                                chunk_id, subchunk_id) {
    file.path(
        bfpwr_sim_bf_prior_dir(corpus_root, bf_prior),
        "bf_split_chunks",
        sprintf("chunk-%04d", as.integer(chunk_id)),
        sprintf("part-%04d.rds", as.integer(subchunk_id))
    )
}

bfpwr_sim_expected_bf_chunk_file <- function(corpus_root, bf_prior,
                                             chunk_id) {
    file.path(
        bfpwr_sim_bf_prior_dir(corpus_root, bf_prior),
        "bf_chunks",
        sprintf("chunk-%04d.rds", as.integer(chunk_id))
    )
}

bfpwr_sim_validate_bf_file_range <- function(path, bf_prior, design,
                                             chunk_id, replicate_start,
                                             replicate_end) {
    bf_rows <- readRDS(path)
    bfpwr_sim_validate_bf_rows(bf_prior, design, bf_rows)
    expected <- seq.int(as.integer(replicate_start), as.integer(replicate_end))
    observed <- sort(unique(bf_rows$replicate_id))
    if (!identical(observed, expected)) {
        stop("BF rows in ", path, " do not contain expected replicate range")
    }
    if (!all(bf_rows$chunk_id == as.integer(chunk_id))) {
        stop("BF rows in ", path, " contain unexpected chunk_id values")
    }
    invisible(TRUE)
}

bfpwr_sim_bf_file_valid <- function(path, bf_prior, design, chunk_id,
                                    replicate_start, replicate_end) {
    if (!file.exists(path) || file.info(path)$size <= 0) return(FALSE)
    tryCatch({
        bfpwr_sim_validate_bf_file_range(
            path, bf_prior, design, chunk_id, replicate_start, replicate_end)
        TRUE
    }, error = function(e) FALSE)
}

bfpwr_sim_bf_file_present <- function(path) {
    file.exists(path) && file.info(path)$size > 0
}

bfpwr_sim_bf_prior_plan_rows <- function(case_set = "production",
                                         plan_index_file = NULL) {
    designs <- bfpwr_sim_design_case_set(case_set)
    bf_priors <- bfpwr_sim_bf_prior_case_set(case_set)
    plan <- bfpwr_sim_bf_prior_run_plan(bf_priors, designs)
    plan$plan_index <- seq_len(nrow(plan))

    if (is.null(plan_index_file)) {
        return(plan)
    }

    requested <- utils::read.csv(plan_index_file, stringsAsFactors = FALSE,
                                 check.names = FALSE)
    if ("plan_index" %in% names(requested)) {
        hit <- match(as.integer(requested$plan_index), plan$plan_index)
    } else {
        required <- c("bf_prior_id", "design_case_id", "chunk_id")
        missing <- setdiff(required, names(requested))
        if (length(missing) > 0) {
            stop("plan file must contain plan_index or columns: ",
                 paste(required, collapse = ", "))
        }
        key <- paste(plan$bf_prior_id, plan$design_case_id, plan$chunk_id,
                     sep = "\r")
        requested_key <- paste(requested$bf_prior_id,
                               requested$design_case_id,
                               requested$chunk_id, sep = "\r")
        hit <- match(requested_key, key)
    }
    if (any(is.na(hit))) {
        stop("could not map all rows in plan file to the production BF plan")
    }
    out <- plan[hit, , drop = FALSE]
    row.names(out) <- NULL
    out$source_plan_row <- seq_len(nrow(out))
    out
}

bfpwr_sim_bf_prior_slow_reason <- function(bf_prior, design,
                                           replicate_start = NULL,
                                           replicate_end = NULL,
                                           row_threshold = 200000L) {
    if (bf_prior$test_family != "t") return(NA_character_)

    prior <- bf_prior$analysis_prior
    reasons <- character()
    n_replicates <- if (!is.null(replicate_start) && !is.null(replicate_end)) {
        as.integer(replicate_end) - as.integer(replicate_start) + 1L
    } else {
        design$chunk_size
    }
    estimated_rows <- n_replicates * length(design$look_grid)
    one_sided <- prior$alternative %in% c("less", "greater")

    if (one_sided && estimated_rows >= row_threshold) {
        reasons <- c(reasons, "large one-sided t BF grid")
    }

    d <- bfpwr_sim_design_prior_mean_sd(design)
    if (one_sided) {
        opposite <- (prior$alternative == "less" && d$dpm > 0) ||
            (prior$alternative == "greater" && d$dpm < 0)
        if (opposite) {
            reasons <- c(reasons,
                         "one-sided t prior is opposite the design mean")
        }
    }

    if (one_sided && design$tier == "long" && !isTRUE(all.equal(prior$null, 0))) {
        reasons <- c(reasons, "long shifted-null one-sided t BF")
    }

    if (length(reasons) == 0) return(NA_character_)
    paste(unique(reasons), collapse = "; ")
}

bfpwr_sim_split_range <- function(replicate_start, replicate_end, split_size) {
    starts <- seq.int(as.integer(replicate_start), as.integer(replicate_end),
                      by = as.integer(split_size))
    ends <- pmin(starts + as.integer(split_size) - 1L,
                 as.integer(replicate_end))
    data.frame(
        subchunk_id = seq_along(starts),
        replicate_start = starts,
        replicate_end = ends
    )
}

bfpwr_sim_split_bf_prior_plan <- function(plan, split_size = 25L,
                                          case_set = "production",
                                          only_slow = FALSE,
                                          corpus_root = NULL,
                                          skip_present = FALSE,
                                          validate_present = FALSE) {
    designs <- bfpwr_sim_design_case_set(case_set)
    bf_priors <- bfpwr_sim_bf_prior_case_set(case_set)
    design_ids <- vapply(designs, function(x) x$design_case_id, character(1))
    bf_prior_ids <- vapply(bf_priors, function(x) x$bf_prior_id, character(1))

    rows <- vector("list", nrow(plan))
    used <- 0L
    for (i in seq_len(nrow(plan))) {
        design <- designs[[match(plan$design_case_id[[i]], design_ids)]]
        bf_prior <- bf_priors[[match(plan$bf_prior_id[[i]], bf_prior_ids)]]
        if (is.null(design) || is.null(bf_prior)) {
            stop("could not resolve design/BF prior for plan row ", i)
        }

        slow_reason <- bfpwr_sim_bf_prior_slow_reason(
            bf_prior, design, plan$replicate_start[[i]],
            plan$replicate_end[[i]])
        slow_task <- !is.na(slow_reason)
        if (only_slow && !slow_task) next

        if (skip_present) {
            if (is.null(corpus_root)) {
                stop("corpus_root is required when skip_present = TRUE")
            }
            chunk_file <- bfpwr_sim_expected_bf_chunk_file(
                corpus_root, bf_prior, plan$chunk_id[[i]])
            present <- if (validate_present) {
                bfpwr_sim_bf_file_valid(
                    chunk_file, bf_prior, design, plan$chunk_id[[i]],
                    plan$replicate_start[[i]], plan$replicate_end[[i]])
            } else {
                bfpwr_sim_bf_file_present(chunk_file)
            }
            if (present) {
                next
            }
        }

        split <- bfpwr_sim_split_range(plan$replicate_start[[i]],
                                       plan$replicate_end[[i]], split_size)
        used <- used + 1L
        rows[[used]] <- data.frame(
            plan_index = plan$plan_index[[i]],
            source_plan_row = if ("source_plan_row" %in% names(plan)) {
                plan$source_plan_row[[i]]
            } else {
                i
            },
            phase = plan$phase[[i]],
            tier = plan$tier[[i]],
            family = plan$family[[i]],
            design_case_id = plan$design_case_id[[i]],
            bf_prior_id = plan$bf_prior_id[[i]],
            chunk_id = plan$chunk_id[[i]],
            chunk_replicate_start = plan$replicate_start[[i]],
            chunk_replicate_end = plan$replicate_end[[i]],
            subchunk_id = split$subchunk_id,
            replicate_start = split$replicate_start,
            replicate_end = split$replicate_end,
            split_size = as.integer(split_size),
            slow_task = slow_task,
            slow_reason = if (slow_task) slow_reason else NA_character_,
            stringsAsFactors = FALSE
        )
    }

    out <- if (used == 0L) {
        data.frame(
            split_plan_index = integer(),
            plan_index = integer(),
            source_plan_row = integer(),
            phase = character(),
            tier = character(),
            family = character(),
            design_case_id = character(),
            bf_prior_id = character(),
            chunk_id = integer(),
            chunk_replicate_start = integer(),
            chunk_replicate_end = integer(),
            subchunk_id = integer(),
            replicate_start = integer(),
            replicate_end = integer(),
            split_size = integer(),
            slow_task = logical(),
            slow_reason = character()
        )
    } else {
        do.call(rbind, rows[seq_len(used)])
    }
    out$split_plan_index <- seq_len(nrow(out))
    out[, c("split_plan_index", setdiff(names(out), "split_plan_index"))]
}

bfpwr_sim_bf_prior_split_merge_plan <- function(split_plan) {
    if (nrow(split_plan) == 0L) {
        return(data.frame(
            merge_plan_index = integer(),
            plan_index = integer(),
            phase = character(),
            tier = character(),
            family = character(),
            design_case_id = character(),
            bf_prior_id = character(),
            chunk_id = integer(),
            replicate_start = integer(),
            replicate_end = integer(),
            n_subchunks = integer(),
            slow_task = logical(),
            slow_reason = character()
        ))
    }
    keys <- unique(split_plan$plan_index)
    rows <- lapply(seq_along(keys), function(i) {
        x <- split_plan[split_plan$plan_index == keys[[i]], , drop = FALSE]
        data.frame(
            merge_plan_index = i,
            plan_index = x$plan_index[[1]],
            phase = x$phase[[1]],
            tier = x$tier[[1]],
            family = x$family[[1]],
            design_case_id = x$design_case_id[[1]],
            bf_prior_id = x$bf_prior_id[[1]],
            chunk_id = x$chunk_id[[1]],
            replicate_start = x$chunk_replicate_start[[1]],
            replicate_end = x$chunk_replicate_end[[1]],
            n_subchunks = nrow(x),
            slow_task = any(x$slow_task),
            slow_reason = paste(unique(stats::na.omit(x$slow_reason)),
                                collapse = "; "),
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, rows)
}
