args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
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

safe_system <- function(command, args = character()) {
    out <- tryCatch(
        system2(command, args = args, stdout = TRUE, stderr = TRUE),
        error = function(e) conditionMessage(e)
    )
    status <- attr(out, "status")
    if (is.null(status)) status <- 0L
    list(status = status, output = out)
}

safe_git <- function(args, project_dir = NA_character_) {
    if (!is.na(project_dir) && dir.exists(project_dir)) {
        return(safe_system("git", c("-C", project_dir, args)))
    }
    safe_system("git", args)
}

source_snapshot <- function(root) {
    if (is.na(root) || !dir.exists(root)) {
        root <- getwd()
    }
    root <- normalizePath(root, winslash = "/", mustWork = FALSE)
    include_dirs <- c(
        "package",
        file.path("simulations", "R"),
        file.path("simulations", "scripts"),
        file.path("simulations", "registry"),
        file.path("simulations", "schema"),
        file.path("simulations", "jobs")
    )
    files <- unlist(lapply(file.path(root, include_dirs), function(dir) {
        if (!dir.exists(dir)) {
            return(character())
        }
        list.files(dir, recursive = TRUE, full.names = TRUE, all.files = TRUE,
                   no.. = TRUE)
    }), use.names = FALSE)
    if (!length(files)) {
        return(data.frame())
    }
    info <- file.info(files)
    files <- files[!info$isdir]
    info <- info[!info$isdir, , drop = FALSE]
    normalized <- normalizePath(files, winslash = "/", mustWork = FALSE)
    rel <- sub(paste0("^", root, "/?"), "", normalized)
    keep <- !grepl("(^|/)[.]Rproj[.]user(/|$)|(^|/)[.]git(/|$)", rel)
    files <- files[keep]
    info <- info[keep, , drop = FALSE]
    rel <- rel[keep]
    data.frame(
        path = rel,
        size = as.numeric(info$size),
        mtime = format(info$mtime, tz = "UTC", usetz = TRUE),
        md5 = unname(tools::md5sum(files)),
        stringsAsFactors = FALSE
    )
}

package_version_or_na <- function(package) {
    if (requireNamespace(package, quietly = TRUE)) {
        as.character(utils::packageVersion(package))
    } else {
        NA_character_
    }
}

args <- parse_args(args)
output_file <- args[["output-file"]]
if (is.null(output_file)) {
    output_file <- file.path("simulations", "corpus", "v1", "environment",
                             "environment.rds")
}

project_dir <- Sys.getenv("BFPWR_PROJECT_DIR", unset = NA_character_)

description <- if (file.exists(file.path("package", "DESCRIPTION"))) {
    as.list(read.dcf(file.path("package", "DESCRIPTION"))[1, ])
} else {
    list()
}

registry_counts <- NULL
if (file.exists(file.path("simulations", "scripts", "common.R"))) {
    registry_counts <- tryCatch({
        source(file.path("simulations", "scripts", "common.R"))
        designs <- bfpwr_sim_design_cases()
        data.frame(
            family = vapply(designs, function(x) x$test_family, character(1)),
            tier = vapply(designs, function(x) x$tier, character(1))
        )
    }, error = function(e) conditionMessage(e))
}

packages <- c("lamW", "mvtnorm", "withr", "qrng", "tinytest")
record <- list(
    recorded_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    working_directory = normalizePath(getwd(), winslash = "/", mustWork = FALSE),
    project_directory = normalizePath(project_dir, winslash = "/",
                                      mustWork = FALSE),
    command_args = commandArgs(FALSE),
    sys_info = as.list(Sys.info()),
    r = list(
        version_string = R.version.string,
        version = R.version,
        home = R.home(),
        lib_paths = .libPaths()
    ),
    environment = as.list(Sys.getenv(c(
        "R_LIBS_USER", "LOADEDMODULES", "MODULEPATH", "PBS_JOBID",
        "PBS_ARRAY_INDEX", "HOSTNAME", "SCRATCHDIR", "BFPWR_PROJECT_DIR",
        "BFPWR_SIM_CORPUS", "BFPWR_DESIGN_CASE", "BFPWR_R_MODULE",
        "BFPWR_R_LIBS_USER"
    ), unset = NA_character_)),
    git = list(
        rev_parse_head = safe_git(c("rev-parse", "HEAD"), project_dir),
        status_short = safe_git(c("status", "--short"), project_dir)
    ),
    source_snapshot = source_snapshot(project_dir),
    package_description = description,
    package_versions = setNames(vapply(packages, package_version_or_na, character(1)),
                                packages),
    registry_counts = registry_counts,
    session_info = utils::sessionInfo()
)

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(record, output_file, compress = "xz")

txt <- sub("[.]rds$", ".txt", output_file)
capture.output(str(record, max.level = 3), file = txt)
cat("wrote environment record to", normalizePath(output_file, winslash = "/"),
    "\n")
