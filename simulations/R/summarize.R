bfpwr_sim_summarize_analysis <- function(analysis, design, outcomes) {
    nsim <- nrow(outcomes)
    pH1 <- mean(outcomes$decision == "H1")
    pH0 <- mean(outcomes$decision == "H0")
    pInc <- mean(outcomes$decision == "inconclusive")
    EN <- mean(outcomes$stop_n)
    VarN <- mean(outcomes$stop_n^2) - EN^2
    ref <- bfpwr_sim_reference_summary(analysis, design)

    data.frame(
        analysis_case_id = analysis$analysis_case_id,
        design_case_id = design$design_case_id,
        test_family = analysis$test_family,
        bf_type = analysis$bf_type,
        nsim = nsim,
        look_grid_name = design$look_grid_name,
        k1 = analysis$k1,
        k0 = analysis$k0,
        pH1 = pH1,
        pH0 = pH0,
        pInc = pInc,
        EN = EN,
        VarN = VarN,
        mcse_pH1 = sqrt(pH1 * (1 - pH1) / nsim),
        mcse_pH0 = sqrt(pH0 * (1 - pH0) / nsim),
        mcse_pInc = sqrt(pInc * (1 - pInc) / nsim),
        reference_pH1 = ref$reference_pH1,
        reference_pH0 = ref$reference_pH0,
        reference_pInc = ref$reference_pInc,
        reference_EN = ref$reference_EN,
        reference_VarN = ref$reference_VarN,
        status = ref$status,
        notes = ref$notes,
        stringsAsFactors = FALSE
    )
}

bfpwr_sim_reference_summary <- function(analysis, design) {
    empty <- list(reference_pH1 = NA_real_, reference_pH0 = NA_real_,
                  reference_pInc = NA_real_, reference_EN = NA_real_,
                  reference_VarN = NA_real_, status = "no_reference",
                  notes = "no sequential deterministic reference implemented")

    if (analysis$test_family == "z") {
        dp <- bfpwr_sim_design_prior_mean_sd(design)
        prior <- analysis$analysis_prior
        null <- prior$null
        se <- design$generation$usd / sqrt(design$look_grid)
        args <- list(k1 = analysis$k1, k0 = analysis$k0, se = se,
                     n = design$look_grid, psd = prior$psd,
                     dpm = dp$dpm - null, dpsd = dp$dpsd,
                     type = analysis$bf_type, strict = analysis$strict)
        if (analysis$bf_type != "moment") {
            args$pm <- prior$pm - null
        }
        res <- try(do.call(bfpwr_sim_package_function("pbf01seq"), args),
                   silent = TRUE)
        if (inherits(res, "try-error")) {
            return(utils::modifyList(empty, list(status = "reference_error",
                                                 notes = as.character(res))))
        }
        return(list(reference_pH1 = tail(res$cumpH1, 1),
                    reference_pH0 = tail(res$cumpH0, 1),
                    reference_pInc = tail(res$cumpInc, 1),
                    reference_EN = res$EN,
                    reference_VarN = res$VarN,
                    status = "ok",
                    notes = "pbf01seq"))
    }

    if (analysis$test_family == "t") {
        dp <- bfpwr_sim_design_prior_mean_sd(design)
        prior <- analysis$analysis_prior
        n1 <- design$look_grid
        n2 <- if (prior$type == "two.sample") {
            pmax(2, round(n1 * design$generation$n2_multiplier))
        } else {
            n1
        }
        args <- list(k1 = analysis$k1, k0 = analysis$k0, n = n1,
                     n1 = n1, n2 = n2,
                     plocation = prior$plocation - prior$null,
                     pscale = prior$pscale, pdf = prior$pdf,
                     dpm = dp$dpm - prior$null, dpsd = dp$dpsd,
                     type = prior$type,
                     alternative = prior$alternative,
                     strict = analysis$strict, drange = analysis$drange)
        res <- try(do.call(bfpwr_sim_package_function("ptbf01seq"), args),
                   silent = TRUE)
        if (inherits(res, "try-error")) {
            return(utils::modifyList(empty, list(status = "reference_error",
                                                 notes = as.character(res))))
        }
        return(list(reference_pH1 = tail(res$cumpH1, 1),
                    reference_pH0 = tail(res$cumpH0, 1),
                    reference_pInc = tail(res$cumpInc, 1),
                    reference_EN = res$EN1,
                    reference_VarN = res$VarN1,
                    status = "ok",
                    notes = "ptbf01seq"))
    }

    if (analysis$test_family == "binomial") {
        return(list(reference_pH1 = NA_real_,
                    reference_pH0 = NA_real_,
                    reference_pInc = NA_real_,
                    reference_EN = NA_real_,
                    reference_VarN = NA_real_,
                    status = "no_sequential_reference",
                    notes = "use fixed-look queries for pbinbf01 references"))
    }

    empty
}

bfpwr_sim_reference_fixed_power <- function(analysis, design, n) {
    prior <- analysis$analysis_prior
    if (analysis$test_family == "z") {
        dp <- bfpwr_sim_design_prior_mean_sd(design)
        if (analysis$bf_type == "normal") {
            pH1 <- bfpwr_sim_package_function("pbf01")(
                k = analysis$k1, n = n, usd = design$generation$usd,
                null = prior$null, pm = prior$pm, psd = prior$psd,
                dpm = dp$dpm, dpsd = dp$dpsd, lower.tail = TRUE)
            pH0 <- bfpwr_sim_package_function("pbf01")(
                k = analysis$k0, n = n, usd = design$generation$usd,
                null = prior$null, pm = prior$pm, psd = prior$psd,
                dpm = dp$dpm, dpsd = dp$dpsd, lower.tail = FALSE)
        } else if (analysis$bf_type == "moment") {
            pH1 <- bfpwr_sim_package_function("pnmbf01")(
                k = analysis$k1, n = n, usd = design$generation$usd,
                null = prior$null, psd = prior$psd,
                dpm = dp$dpm, dpsd = dp$dpsd, lower.tail = TRUE)
            pH0 <- bfpwr_sim_package_function("pnmbf01")(
                k = analysis$k0, n = n, usd = design$generation$usd,
                null = prior$null, psd = prior$psd,
                dpm = dp$dpm, dpsd = dp$dpsd, lower.tail = FALSE)
        } else if (analysis$bf_type == "directional") {
            res <- bfpwr_sim_package_function("pbf01seq")(
                k1 = analysis$k1,
                k0 = analysis$k0,
                se = design$generation$usd / sqrt(n),
                n = n,
                pm = prior$pm - prior$null,
                psd = prior$psd,
                dpm = dp$dpm - prior$null,
                dpsd = dp$dpsd,
                type = "directional",
                strict = analysis$strict)
            pH1 <- tail(res$cumpH1, 1)
            pH0 <- tail(res$cumpH0, 1)
        } else {
            pH1 <- pH0 <- NA_real_
        }
    } else if (analysis$test_family == "t") {
        dp <- bfpwr_sim_design_prior_mean_sd(design)
        pH1 <- bfpwr_sim_package_function("ptbf01")(
            k = analysis$k1, n = n,
            n1 = n,
            n2 = if (prior$type == "two.sample") {
                pmax(2, round(n * design$generation$n2_multiplier))
            } else {
                n
            },
            null = 0,
            plocation = prior$plocation - prior$null,
            pscale = prior$pscale,
            pdf = prior$pdf,
            dpm = dp$dpm - prior$null,
            dpsd = dp$dpsd,
            type = prior$type, alternative = prior$alternative,
            lower.tail = TRUE, drange = analysis$drange)
        pH0 <- bfpwr_sim_package_function("ptbf01")(
            k = analysis$k0, n = n,
            n1 = n,
            n2 = if (prior$type == "two.sample") {
                pmax(2, round(n * design$generation$n2_multiplier))
            } else {
                n
            },
            null = 0,
            plocation = prior$plocation - prior$null,
            pscale = prior$pscale,
            pdf = prior$pdf,
            dpm = dp$dpm - prior$null,
            dpsd = dp$dpsd,
            type = prior$type, alternative = prior$alternative,
            lower.tail = FALSE, drange = analysis$drange)
    } else if (analysis$test_family == "binomial") {
        dp <- bfpwr_sim_binomial_design_args(design)
        pH1 <- bfpwr_sim_package_function("pbinbf01")(
            k = analysis$k1, n = n, p0 = prior$p0, type = analysis$bf_type,
            a = prior$a, b = prior$b, dp = dp$dp, da = dp$da,
            db = dp$db, dl = dp$dl, du = dp$du, lower.tail = TRUE)
        pH0 <- bfpwr_sim_package_function("pbinbf01")(
            k = analysis$k0, n = n, p0 = prior$p0, type = analysis$bf_type,
            a = prior$a, b = prior$b, dp = dp$dp, da = dp$da,
            db = dp$db, dl = dp$dl, du = dp$du, lower.tail = FALSE)
    } else {
        pH1 <- pH0 <- NA_real_
    }
    data.frame(n = n, reference_pH1 = pH1, reference_pH0 = pH0)
}
