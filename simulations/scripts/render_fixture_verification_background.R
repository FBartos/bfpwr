source(file.path("simulations", "scripts", "common.R"))

corpus_root <- bfpwr_sim_require_corpus_root(
    "simulations/corpus/v1",
    asset_set = "fixture-tests",
    purpose = "fixture verification vignette rendering")

Sys.setenv(
    BFPWR_SIM_CORPUS = corpus_root,
    BFPWR_SIM_REPO = normalizePath(".", winslash = "/", mustWork = TRUE)
)

cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Repo:", getwd(), "\n")
cat("Corpus:", Sys.getenv("BFPWR_SIM_CORPUS"), "\n")

t0 <- proc.time()[["elapsed"]]
old <- setwd("simulations")
on.exit(setwd(old), add = TRUE)
knitr::knit2pdf("fixture-verification.Rnw", clean = TRUE)

cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Elapsed seconds:", proc.time()[["elapsed"]] - t0, "\n")
