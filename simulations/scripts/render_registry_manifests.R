source(file.path("simulations", "scripts", "common.R"))

args <- bfpwr_sim_parse_args()
repo_root <- bfpwr_sim_repo_root()
output_dir <- if (is.null(args[["output-dir"]])) {
    file.path(repo_root, "simulations", "registry", "manifests")
} else {
    args[["output-dir"]]
}

designs <- bfpwr_sim_design_cases()
analyses <- bfpwr_sim_analysis_cases(designs = designs)
bf_priors <- bfpwr_sim_bf_prior_cases(designs = designs)

bfpwr_sim_write_registry_manifests(output_dir, designs, analyses, bf_priors)
cat("wrote registry manifests to", output_dir, "\n")
