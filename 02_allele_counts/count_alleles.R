#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
parse_args <- function(x) {
  if (length(x) %% 2L != 0L) stop("Arguments must be --key value pairs")
  keys <- x[seq(1L, length(x), by = 2L)]
  if (any(!startsWith(keys, "--"))) stop("Expected --key value pairs")
  setNames(as.list(x[seq(2L, length(x), by = 2L)]), substring(keys, 3L))
}
`%||%` <- function(x, y) if (is.null(x)) y else x
opt <- parse_args(args)
required <- c("bam", "loci-dir", "output-dir", "sample-id")
missing <- required[!vapply(required, function(k) nzchar(opt[[k]] %||% ""), logical(1))]
if (length(missing)) stop("Missing required arguments: ", paste(missing, collapse = ", "))
bam <- normalizePath(opt[["bam"]], mustWork = TRUE)
loci_dir <- normalizePath(opt[["loci-dir"]], mustWork = TRUE)
out_root <- opt[["output-dir"]]
sample_id <- opt[["sample-id"]]
allele_counter <- opt[["allele-counter"]] %||% "alleleCounter"
min_baq <- as.integer(opt[["min-baq"]] %||% "20")
min_maq <- as.integer(opt[["min-maq"]] %||% "35")
dense <- tolower(opt[["dense-snps"]] %||% "true") == "true"
allelecount_only <- TRUE
chromosomes <- strsplit(opt[["chromosomes"]] %||% paste(paste0("chr", c(1:22, "X")), collapse = ","), ",", fixed = TRUE)[[1L]]
chromosomes <- trimws(chromosomes)
if (any(!grepl("^chr([1-9]|1[0-9]|2[0-2]|X)$", chromosomes))) stop("Invalid chromosome names")
if (anyNA(c(min_baq, min_maq))) stop("Invalid quality thresholds")
if (!nzchar(Sys.which(allele_counter)) && !file.exists(allele_counter)) stop("alleleCounter not found: ", allele_counter)
ac_dir <- file.path(out_root, "allele_count")
dir.create(ac_dir, recursive = TRUE, showWarnings = FALSE)

find_loci <- function(chr) {
  files <- list.files(loci_dir, pattern = "\\.txt$", full.names = TRUE)
  hits <- files[grepl(paste0("(^|[._-])", chr, "([._-]|$)"), basename(files))]
  if (length(hits) != 1L) stop("Expected one loci file for ", chr, "; found ", length(hits))
  hits
}

for (chr in chromosomes) {
  loci_file <- find_loci(chr)
  count_file <- file.path(ac_dir, paste0("allelecount_", sample_id, "_", chr, ".txt"))
  cmd <- c("-b", bam, "-l", loci_file, "-o", count_file,
           "-m", as.character(min_baq), "-q", as.character(min_maq), "-f", "0", "-F", "3844")
  if (dense) cmd <- c(cmd, "--dense-snps")
  status <- system2(allele_counter, args = shQuote(cmd))
  if (!identical(as.integer(status), 0L) || !file.exists(count_file)) stop("alleleCounter failed for ", chr)
  counts <- fread(count_file)
  if (!nrow(counts)) stop("Empty allele counts for ", chr)
  if ("Good_depth" %in% names(counts)) {
    counts <- counts[Good_depth != 0]
    fwrite(counts, count_file, sep = "\t")
  }
  message("Wrote ", count_file)
}
