bfpwr_sim_set_chunk_seed <- function(seed) {
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
    invisible(seed)
}

bfpwr_sim_run_design_chunk <- function(design, chunk_id) {
    bfpwr_sim_validate_design_case(design)
    chunk <- design$chunks[design$chunks$chunk_id == chunk_id, , drop = FALSE]
    if (nrow(chunk) != 1) {
        stop("unknown chunk_id for design: ", chunk_id)
    }

    bfpwr_sim_set_chunk_seed(chunk$seed)
    replicate_ids <- seq(chunk$replicate_start, chunk$replicate_end)

    out <- switch(design$test_family,
                  z = bfpwr_sim_generate_z_chunk(design, replicate_ids, chunk_id),
                  t = bfpwr_sim_generate_t_chunk(design, replicate_ids, chunk_id),
                  binomial = bfpwr_sim_generate_binomial_chunk(design, replicate_ids, chunk_id))
    bfpwr_sim_validate_design_trajectories(design, out)
    out
}

bfpwr_sim_draw_continuous_effects <- function(design, nrep) {
    prior <- design$design_prior
    if (identical(prior$family, "point")) {
        rep(prior$mean, nrep)
    } else if (identical(prior$family, "normal")) {
        stats::rnorm(nrep, mean = prior$mean, sd = prior$sd)
    } else {
        stop("unsupported continuous design prior: ", prior$family)
    }
}

bfpwr_sim_draw_binomial_probs <- function(design, nrep) {
    prior <- design$design_prior
    if (identical(prior$family, "point")) {
        rep(prior$prob, nrep)
    } else if (identical(prior$family, "beta")) {
        draws <- numeric(nrep)
        filled <- 0L
        while (filled < nrep) {
            proposal <- stats::rbeta(max(1000L, 2L * (nrep - filled)),
                                     shape1 = prior$shape1,
                                     shape2 = prior$shape2)
            proposal <- proposal[proposal >= prior$lower & proposal <= prior$upper]
            take <- min(length(proposal), nrep - filled)
            if (take > 0) {
                draws[seq.int(filled + 1L, filled + take)] <- proposal[seq_len(take)]
                filled <- filled + take
            }
        }
        draws
    } else {
        stop("unsupported binomial design prior: ", prior$family)
    }
}

bfpwr_sim_generate_z_chunk <- function(design, replicate_ids, chunk_id) {
    look_grid <- design$look_grid
    nlook <- length(look_grid)
    usd <- design$generation$usd
    max_n <- max(look_grid)
    effects <- bfpwr_sim_draw_continuous_effects(design, length(replicate_ids))

    rows <- vector("list", length(replicate_ids))
    for (i in seq_along(replicate_ids)) {
        error <- stats::rnorm(max_n, mean = 0, sd = usd)
        estimate <- effects[i] + cumsum(error)[look_grid] / look_grid
        rows[[i]] <- data.frame(
            design_case_id = design$design_case_id,
            replicate_id = replicate_ids[i],
            chunk_id = chunk_id,
            look = seq_len(nlook),
            n = look_grid,
            se = usd / sqrt(look_grid),
            true_effect = effects[i],
            estimate = estimate
        )
    }
    do.call(rbind, rows)
}

bfpwr_sim_generate_t_chunk <- function(design, replicate_ids, chunk_id) {
    test_type <- match.arg(design$generation$type,
                           c("two.sample", "one.sample", "paired"))
    look_grid <- design$look_grid
    nlook <- length(look_grid)
    n1_grid <- look_grid
    n2_grid <- if (test_type == "two.sample") {
        pmax(2, round(look_grid * design$generation$n2_multiplier))
    } else {
        look_grid
    }
    effects <- bfpwr_sim_draw_continuous_effects(design, length(replicate_ids))

    rows <- vector("list", length(replicate_ids))
    for (i in seq_along(replicate_ids)) {
        if (test_type == "two.sample") {
            y1 <- stats::rnorm(max(n1_grid), mean = 0, sd = 1)
            y2 <- stats::rnorm(max(n2_grid), mean = effects[i], sd = 1)
            cs1 <- cumsum(y1)
            cs2 <- cumsum(y2)
            csq1 <- cumsum(y1^2)
            csq2 <- cumsum(y2^2)
            m1 <- cs1[n1_grid] / n1_grid
            m2 <- cs2[n2_grid] / n2_grid
            v1 <- (csq1[n1_grid] - cs1[n1_grid]^2 / n1_grid) / (n1_grid - 1)
            v2 <- (csq2[n2_grid] - cs2[n2_grid]^2 / n2_grid) / (n2_grid - 1)
            sp <- sqrt(((n1_grid - 1) * v1 + (n2_grid - 1) * v2) /
                           (n1_grid + n2_grid - 2))
            estimate <- (m2 - m1) / sp
            tstat <- (m2 - m1) / (sp * sqrt(1 / n1_grid + 1 / n2_grid))
        } else {
            y <- stats::rnorm(max(n1_grid), mean = effects[i], sd = 1)
            cs <- cumsum(y)
            csq <- cumsum(y^2)
            m <- cs[n1_grid] / n1_grid
            s <- sqrt((csq[n1_grid] - cs[n1_grid]^2 / n1_grid) / (n1_grid - 1))
            estimate <- m / s
            tstat <- sqrt(n1_grid) * estimate
        }
        rows[[i]] <- data.frame(
            design_case_id = design$design_case_id,
            replicate_id = replicate_ids[i],
            chunk_id = chunk_id,
            look = seq_len(nlook),
            n = n1_grid,
            n1 = n1_grid,
            n2 = n2_grid,
            true_effect = effects[i],
            estimate = estimate,
            t = tstat
        )
    }
    do.call(rbind, rows)
}

bfpwr_sim_generate_binomial_chunk <- function(design, replicate_ids, chunk_id) {
    look_grid <- design$look_grid
    nlook <- length(look_grid)
    max_n <- max(look_grid)
    probs <- bfpwr_sim_draw_binomial_probs(design, length(replicate_ids))

    rows <- vector("list", length(replicate_ids))
    for (i in seq_along(replicate_ids)) {
        y <- stats::rbinom(max_n, size = 1, prob = probs[i])
        rows[[i]] <- data.frame(
            design_case_id = design$design_case_id,
            replicate_id = replicate_ids[i],
            chunk_id = chunk_id,
            look = seq_len(nlook),
            n = look_grid,
            true_p = probs[i],
            x = cumsum(y)[look_grid]
        )
    }
    do.call(rbind, rows)
}
