source(file.path("simulations", "scripts", "common.R"))

args <- bfpwr_sim_parse_args()
case_set <- if (is.null(args[["case-set"]])) "production" else args[["case-set"]]
plan_index_file <- args[["plan-index-file"]]
output_file <- args[["output-file"]]

plan <- bfpwr_sim_bf_prior_plan_rows(case_set = case_set,
                                     plan_index_file = plan_index_file)
designs <- bfpwr_sim_design_case_set(case_set)
bf_priors <- bfpwr_sim_bf_prior_case_set(case_set)
design_ids <- vapply(designs, function(x) x$design_case_id, character(1))
bf_prior_ids <- vapply(bf_priors, function(x) x$bf_prior_id, character(1))

reason <- character(nrow(plan))
for (i in seq_len(nrow(plan))) {
    design <- designs[[match(plan$design_case_id[[i]], design_ids)]]
    bf_prior <- bf_priors[[match(plan$bf_prior_id[[i]], bf_prior_ids)]]
    reason[[i]] <- bfpwr_sim_bf_prior_slow_reason(
        bf_prior, design, plan$replicate_start[[i]], plan$replicate_end[[i]])
}
plan$slow_task <- !is.na(reason)
plan$slow_reason <- reason

if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(plan, output_file, row.names = FALSE)
}

cat("plan rows:", nrow(plan), "\n")
cat("slow rows:", sum(plan$slow_task), "\n")
if (any(plan$slow_task)) {
    cat("\nslow-task reasons:\n")
    print(sort(table(plan$slow_reason[plan$slow_task]), decreasing = TRUE))
}
if (!is.null(output_file)) cat("wrote:", output_file, "\n")
