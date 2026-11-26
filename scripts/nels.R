# Read files from NeLS (https://nels.elixir.no) directly in R, over SSH.
# Adapted from R/nels.R in the 3D-wgd repo.
#
# For interactive exploration and drafting. Pipelines (run_all.sh) should still
# fetch their inputs with analysis/<name>/download_data.sh, so that they keep
# working when the data moves to Zenodo (see CLAUDE.md).
#
# Setup: the `nels.storage` SSH host from NELS_SETUP.md. Run `ssh nels.storage ls`
# once in a terminal to accept the host keys.
#
# Usage:
#   source("scripts/nels.R")                         # from the repo root
#   nels_ls("data")                                  # list a folder
#   sp  <- nels_read("data/species.tsv")             # stream a table into R
#   f   <- nels_file("data/phylogeny/sco_cds_b500.date.nwk")  # download once, local path
#
# Relative paths are relative to the project folder on NeLS
# (Projects/NMBU_threeD_Sandve_2024/Mucor), which mirrors the repo root. So
# "data/species.tsv" means the same file on NeLS as in the repo, and nels_file()
# saves it to the same place as scripts/download_files.sh. Absolute paths are
# used as-is on NeLS and cached in the user cache folder.
#
# Settings (environment variables or options()):
#   NELS_HOST       / nels.host   SSH host alias       (default "nels.storage")
#   MUCOR_NELS_ROOT / nels.root   project folder on NeLS
#                                 (default "Projects/NMBU_threeD_Sandve_2024/Mucor")
#   NELS_CACHE      / nels.cache  cache for absolute paths
#                                 (default: user cache dir, outside Dropbox)
# MUCOR_NELS_ROOT is separate from the NELS_ROOT used by 3D-wgd, so both
# projects can be configured in the same ~/.Renviron.

.nels_opt <- function(opt, env, default) {
  val <- getOption(opt, Sys.getenv(env, ""))
  if (is.null(val) || identical(val, "")) default else val
}
.nels_host  <- function() .nels_opt("nels.host",  "NELS_HOST",       "nels.storage")
.nels_root  <- function() .nels_opt("nels.root",  "MUCOR_NELS_ROOT",
                                    "Projects/NMBU_threeD_Sandve_2024/Mucor")
.nels_cache <- function() .nels_opt("nels.cache", "NELS_CACHE",
                                    tools::R_user_dir("nels", which = "cache"))
.nels_ssh   <- function() getOption("nels.ssh", "ssh")  # overridable for testing

# Repo root: the folder above scripts/, found when this file is sourced;
# otherwise the git top level of the working directory.
.nels_repo <- local({
  ofile <- NULL
  for (f in rev(sys.frames())) if (!is.null(f$ofile)) { ofile <- f$ofile; break }
  if (!is.null(ofile)) {
    dirname(dirname(normalizePath(ofile)))
  } else {
    top <- suppressWarnings(tryCatch(
      system("git rev-parse --show-toplevel", intern = TRUE, ignore.stderr = TRUE),
      error = function(e) character()))
    if (length(top) == 1) top else getwd()
  }
})

# Full remote path: absolute paths are used as-is, others are relative to the root.
nels_path <- function(path = "") {
  root <- .nels_root()
  if (startsWith(path, "/") || root == "") return(path)
  if (path == "") root else paste(sub("/+$", "", root), path, sep = "/")
}

# Shell command (local) that runs `remote_cmd` on NeLS.
.nels_cmd <- function(remote_cmd) {
  paste(.nels_ssh(), "-o BatchMode=yes", shQuote(.nels_host()), shQuote(remote_cmd))
}

.nels_run <- function(remote_cmd) {
  out <- suppressWarnings(system(.nels_cmd(remote_cmd), intern = TRUE))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0)
    stop("NeLS command failed (exit ", status, "): ", remote_cmd, call. = FALSE)
  out
}

# Why reading `p` failed: a clear message for folders and missing paths
# (checked only after a failure, so successful reads cost no extra SSH call).
.nels_fail <- function(p, what) {
  kind <- tryCatch(.nels_run(paste(
    "if [ -d", shQuote(p), "]; then echo dir;",
    "elif [ -e", shQuote(p), "]; then echo file; else echo missing; fi")),
    error = function(e) "unreachable")
  msg <- switch(kind[1],
    dir = paste0("'", p, "' is a folder, not a file. ",
                 "Use nels_ls() to see its files and give the path of one of them."),
    missing = paste0("'", p, "' does not exist on NeLS. ",
                     "Use nels_ls() to check the folder."),
    unreachable = paste0("Could not connect to NeLS host '", .nels_host(),
                         "'. Check `ssh ", .nels_host(), " ls` in a terminal."),
    paste0(what, " failed for '", p, "' (see the SSH output above)."))
  stop(msg, call. = FALSE)
}

#' List files in a NeLS folder.
#' @param path folder, relative to the project root (default: the root itself)
#' @param long if TRUE, return `ls -l` style lines with sizes and dates
nels_ls <- function(path = "", long = FALSE) {
  p <- nels_path(path)
  .nels_run(paste(if (long) "ls -lh" else "ls -1", if (p == "") "" else shQuote(p)))
}

#' Read a text table (tsv/csv/bed/...; plain or .gz/.bz2) straight from NeLS.
#' Streams over SSH, so nothing is written to disk.
#' @param ... passed on to data.table::fread (e.g. header, select, col.names)
nels_read <- function(path, ...) {
  p <- nels_path(path)
  reader <- if (grepl("\\.gz$", p)) "gzip -dc" else if (grepl("\\.bz2$", p)) "bzip2 -dc" else "cat"
  tryCatch(data.table::fread(cmd = .nels_cmd(paste(reader, shQuote(p))), ...),
           error = function(e) {
             # fread reports any failed command as a possibly full temp disk
             if (grepl("External command failed", conditionMessage(e))) .nels_fail(p, "Reading")
             stop(e)
           })
}

# Local copy of a NeLS file: relative paths mirror into the repo (as
# scripts/download_files.sh does), absolute paths go to the cache.
.nels_local <- function(path) {
  if (startsWith(path, "/")) {
    file.path(.nels_cache(), .nels_host(), sub("^/+", "", path))
  } else {
    file.path(.nels_repo, path)
  }
}

#' Download a file from NeLS once and return its local path.
#' Use for files that need random access (.hic, .mcool, bigWig, BAM, RDS, ...).
#' Later calls reuse the local copy unless `refresh = TRUE`.
nels_file <- function(path, refresh = FALSE) {
  p <- nels_path(path)
  local <- .nels_local(path)
  if (refresh || !file.exists(local)) {
    dir.create(dirname(local), recursive = TRUE, showWarnings = FALSE)
    tmp <- paste0(local, ".part")
    status <- system(paste(.nels_cmd(paste("cat", shQuote(p))), ">", shQuote(tmp)))
    if (status != 0) {
      unlink(tmp)
      .nels_fail(p, "Download")
    }
    file.rename(tmp, local)
  }
  local
}

#' Read an .rds file from NeLS (downloaded once via nels_file()).
nels_readRDS <- function(path, refresh = FALSE) readRDS(nels_file(path, refresh))
