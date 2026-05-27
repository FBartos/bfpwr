source(file.path("simulations", "scripts", "common.R"))

args <- bfpwr_sim_parse_args()
tiers <- if (is.null(args[["tiers"]])) {
    NULL
} else {
    strsplit(args[["tiers"]], ",", fixed = TRUE)[[1]]
}
families <- if (is.null(args[["families"]])) {
    NULL
} else {
    strsplit(args[["families"]], ",", fixed = TRUE)[[1]]
}

designs <- bfpwr_sim_design_cases(tiers = tiers, families = families)
analyses <- bfpwr_sim_analysis_cases(tiers = tiers, families = families,
                                     designs = bfpwr_sim_design_cases())
bf_priors <- bfpwr_sim_bf_prior_cases(tiers = tiers, families = families,
                                      designs = bfpwr_sim_design_cases())
bfpwr_sim_validate_registry(bfpwr_sim_design_cases(), bfpwr_sim_analysis_cases(),
                            bfpwr_sim_bf_prior_cases())

cat("registry ok\n")
cat("selected designs:", length(designs), "\n")
cat("selected analyses:", length(analyses), "\n")
cat("selected BF priors:", length(bf_priors), "\n")
if (length(designs) > 0) {
    print(table(vapply(designs, function(x) x$tier, character(1)),
                vapply(designs, function(x) x$test_family, character(1))))
}
if (length(bf_priors) > 0) {
    cat("\nBF priors by tier/family:\n")
    print(table(vapply(bf_priors, function(x) x$tier, character(1)),
                vapply(bf_priors, function(x) x$test_family, character(1))))
}
