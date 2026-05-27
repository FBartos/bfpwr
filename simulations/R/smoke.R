bfpwr_sim_smoke_design_cases <- function() {
    list(
        bfpwr_sim_design_case(
            design_case_id = "smoke-z-dnorm-m0p4-s0p1-gridshort",
            test_family = "z",
            look_grid_name = "short",
            nsim = 30,
            chunk_size = 10,
            master_seed = 20260524,
            design_prior = list(family = "normal", mean = 0.4, sd = 0.1),
            generation = list(usd = sqrt(2)),
            tier = "short",
            tags = c("developer-smoke"),
            rationale = "Developer-only architecture smoke case; not a stored fixture tier."
        ),
        bfpwr_sim_design_case(
            design_case_id = "smoke-t-dpoint-m0p4-gridshort",
            test_family = "t",
            look_grid_name = "short",
            nsim = 12,
            chunk_size = 6,
            master_seed = 20260525,
            design_prior = list(family = "point", mean = 0.4, sd = 0),
            generation = list(type = "two.sample", n2_multiplier = 1),
            tier = "short",
            tags = c("developer-smoke"),
            rationale = "Developer-only architecture smoke case; not a stored fixture tier."
        ),
        bfpwr_sim_design_case(
            design_case_id = "smoke-binomial-dpoint-p0p6-gridshort",
            test_family = "binomial",
            look_grid_name = "short",
            nsim = 30,
            chunk_size = 10,
            master_seed = 20260526,
            design_prior = list(family = "point", prob = 0.6),
            generation = list(),
            tier = "short",
            tags = c("developer-smoke"),
            rationale = "Developer-only architecture smoke case; not a stored fixture tier."
        )
    )
}

bfpwr_sim_smoke_analysis_cases <- function() {
    list(
        bfpwr_sim_analysis_case(
            analysis_case_id = "smoke-z-normal-pm0-psd1-k0p1-10",
            design_case_id = "smoke-z-dnorm-m0p4-s0p1-gridshort",
            test_family = "z",
            bf_type = "normal",
            k1 = 1/10,
            k0 = 10,
            analysis_prior = list(null = 0, pm = 0, psd = 1),
            strict = FALSE,
            tier = "short",
            tags = c("developer-smoke"),
            rationale = "Developer-only architecture smoke case; not a stored fixture tier."
        ),
        bfpwr_sim_analysis_case(
            analysis_case_id = "smoke-z-moment-psd0p5-k0p1-10",
            design_case_id = "smoke-z-dnorm-m0p4-s0p1-gridshort",
            test_family = "z",
            bf_type = "moment",
            k1 = 1/10,
            k0 = 10,
            analysis_prior = list(null = 0, psd = 0.5),
            strict = FALSE,
            tier = "short",
            tags = c("developer-smoke"),
            rationale = "Developer-only architecture smoke case; not a stored fixture tier."
        ),
        bfpwr_sim_analysis_case(
            analysis_case_id = "smoke-t-jzs-greater-k0p1-10",
            design_case_id = "smoke-t-dpoint-m0p4-gridshort",
            test_family = "t",
            bf_type = "t",
            k1 = 1/10,
            k0 = 10,
            analysis_prior = list(null = 0, plocation = 0, pscale = 1/sqrt(2),
                                  pdf = 1, type = "two.sample",
                                  alternative = "greater"),
            strict = FALSE,
            tier = "short",
            tags = c("developer-smoke"),
            rationale = "Developer-only architecture smoke case; not a stored fixture tier."
        ),
        bfpwr_sim_analysis_case(
            analysis_case_id = "smoke-binomial-point-beta1-1-k0p1-10",
            design_case_id = "smoke-binomial-dpoint-p0p6-gridshort",
            test_family = "binomial",
            bf_type = "point",
            k1 = 1/10,
            k0 = 10,
            analysis_prior = list(p0 = 0.5, a = 1, b = 1),
            strict = FALSE,
            tier = "short",
            tags = c("developer-smoke"),
            rationale = "Developer-only architecture smoke case; not a stored fixture tier."
        )
    )
}

bfpwr_sim_smoke_bf_prior_cases <- function() {
    designs <- bfpwr_sim_smoke_design_cases()
    ids <- vapply(designs, function(x) x$design_case_id, character(1))
    by_id <- function(id) designs[[match(id, ids)]]
    list(
        bfpwr_sim_bf_prior_z_bf01(
            by_id("smoke-z-dnorm-m0p4-s0p1-gridshort"),
            null = 0, pm = 0.2, psd = 0,
            tags = "developer-smoke",
            rationale = "Developer-only BF prior smoke case."),
        bfpwr_sim_bf_prior_z_moment(
            by_id("smoke-z-dnorm-m0p4-s0p1-gridshort"),
            null = 0, psd = 0.5,
            tags = "developer-smoke",
            rationale = "Developer-only BF prior smoke case."),
        bfpwr_sim_bf_prior_t(
            by_id("smoke-t-dpoint-m0p4-gridshort"),
            null = 0, plocation = 0, pscale = 1 / sqrt(2), pdf = 1,
            alternative = "greater",
            tags = "developer-smoke",
            rationale = "Developer-only BF prior smoke case."),
        bfpwr_sim_bf_prior_binomial(
            by_id("smoke-binomial-dpoint-p0p6-gridshort"),
            bf_type = "point", p0 = 0.5, a = 1, b = 1,
            tags = "developer-smoke",
            rationale = "Developer-only BF prior smoke case.")
    )
}

bfpwr_sim_smoke_run <- function(corpus_root = tempfile("bfpwr-sim-corpus-")) {
    designs <- bfpwr_sim_smoke_design_cases()
    analyses <- bfpwr_sim_smoke_analysis_cases()

    dir.create(corpus_root, recursive = TRUE, showWarnings = FALSE)
    design_trajectories <- list()

    for (design in designs) {
        design_dir <- bfpwr_sim_design_dir(corpus_root, design)
        bfpwr_sim_write_rds(design, file.path(design_dir, "design.rds"))
        chunks <- vector("list", nrow(design$chunks))
        for (chunk_id in design$chunks$chunk_id) {
            chunk <- bfpwr_sim_run_design_chunk(design, chunk_id)
            chunks[[chunk_id]] <- chunk
            bfpwr_sim_write_rds(chunk, file.path(design_dir, "chunks",
                                                 sprintf("chunk-%04d.rds", chunk_id)))
        }
        trajectories <- do.call(rbind, chunks)
        bfpwr_sim_validate_design_trajectories(design, trajectories)
        design_trajectories[[design$design_case_id]] <- trajectories
        bfpwr_sim_write_rds(trajectories, file.path(design_dir, "trajectories.rds"))
    }

    summaries <- list()
    fixed_queries <- list()
    for (analysis in analyses) {
        design <- designs[[which(vapply(designs, function(x) x$design_case_id,
                                        character(1)) == analysis$design_case_id)]]
        trajectories <- design_trajectories[[design$design_case_id]]
        materialized <- bfpwr_sim_materialize_analysis(
            analysis, design, trajectories, include_bf_trajectories = TRUE)
        analysis_dir <- bfpwr_sim_analysis_dir(corpus_root, analysis)
        bfpwr_sim_write_rds(analysis, file.path(analysis_dir, "analysis.rds"))
        bfpwr_sim_write_rds(materialized$outcomes, file.path(analysis_dir, "outcomes.rds"))
        bfpwr_sim_write_rds(materialized$summary, file.path(analysis_dir, "summary.rds"))
        utils::write.csv(materialized$summary, file.path(analysis_dir, "summary.csv"),
                         row.names = FALSE)
        summaries[[analysis$analysis_case_id]] <- materialized$summary
        fixed_queries[[analysis$analysis_case_id]] <- bfpwr_sim_query_fixed_power(
            analysis, design, trajectories, n = 100)
    }

    summary <- do.call(rbind, summaries)
    fixed <- do.call(rbind, fixed_queries)
    list(corpus_root = corpus_root, summary = summary, fixed = fixed,
         validation = bfpwr_sim_validate_summary(summary))
}
