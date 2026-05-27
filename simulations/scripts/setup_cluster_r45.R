parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
    out <- list()
    i <- 1L
    while (i <= length(args)) {
        key <- args[[i]]
        if (!startsWith(key, "--")) {
            stop("unexpected positional argument: ", key)
        }
        key <- sub("^--", "", key)
        if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
            out[[key]] <- TRUE
            i <- i + 1L
        } else {
            out[[key]] <- args[[i + 1L]]
            i <- i + 2L
        }
    }
    out
}

package_version_or_na <- function(package) {
    if (requireNamespace(package, quietly = TRUE)) {
        as.character(utils::packageVersion(package))
    } else {
        NA_character_
    }
}

args <- parse_args()
lib <- args[["lib"]]
if (is.null(lib)) {
    lib <- Sys.getenv("R_LIBS_USER", unset = "/auto/brno2/home/fbartos/Rpackages-45")
}
repo_root <- args[["repo-root"]]
if (is.null(repo_root)) {
    repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
package_dir <- args[["package-dir"]]
if (is.null(package_dir)) {
    package_dir <- file.path(repo_root, "package")
}
repos <- args[["repos"]]
if (is.null(repos)) {
    repos <- "https://cloud.r-project.org"
}
force <- isTRUE(args[["force"]])
ncpus <- as.integer(Sys.getenv("BFPWR_INSTALL_NCPUS", unset = "1"))
if (is.na(ncpus) || ncpus < 1) ncpus <- 1L

cran_packages <- c("lamW", "mvtnorm", "withr", "qrng", "tinytest")

dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))

cat("R:", R.version.string, "\n")
cat("library:", normalizePath(lib, winslash = "/", mustWork = FALSE), "\n")
cat("package dir:", normalizePath(package_dir, winslash = "/", mustWork = FALSE), "\n")
cat("repos:", repos, "\n")
print(.libPaths())

installed_ok <- vapply(cran_packages, requireNamespace, logical(1),
                       quietly = TRUE)
to_install <- if (force) cran_packages else cran_packages[!installed_ok]
if (length(to_install) > 0) {
    cat("installing CRAN packages:", paste(to_install, collapse = ", "), "\n")
    utils::install.packages(to_install, lib = lib, repos = repos,
                            dependencies = c("Depends", "Imports", "LinkingTo"),
                            Ncpus = ncpus)
} else {
    cat("all CRAN packages already available\n")
}

missing <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1),
                                 quietly = TRUE)]
if (length(missing) > 0) {
    stop("missing packages after install: ", paste(missing, collapse = ", "))
}

if (!file.exists(file.path(package_dir, "DESCRIPTION"))) {
    stop("package DESCRIPTION not found under ", package_dir)
}

install_args <- c("CMD", "INSTALL", "--no-multiarch", "--with-keep.source",
                  paste0("--library=", lib), package_dir)
cat("installing local package with: R", paste(install_args, collapse = " "), "\n")
status <- system2(file.path(R.home("bin"), "R"), install_args)
if (!identical(status, 0L)) {
    stop("local package install failed with status ", status)
}

if (!requireNamespace("bfpwr", quietly = TRUE)) {
    stop("installed bfpwr package cannot be loaded")
}

cat("installed package versions:\n")
print(setNames(vapply(c(cran_packages, "bfpwr"), package_version_or_na,
                      character(1)), c(cran_packages, "bfpwr")))

cat("setup complete\n")
