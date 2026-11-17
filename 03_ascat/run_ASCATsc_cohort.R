suppressPackageStartupMessages({

  library(ASCAT.sc)

  library(Rsamtools)

})

task_id <- as.integer(Sys.getenv("SGE_TASK_ID", "1"))

tasks <- expand.grid(

  donor = trimws(strsplit(Sys.getenv("ASCAT_SAMPLE_IDS", "sample_01"), ",", fixed = TRUE)[[1]]),

  bin = c("30kb", "240kb", "500kb", "5mb"),

  stringsAsFactors = FALSE

)

if (is.na(task_id) || task_id < 1L || task_id > nrow(tasks)) {

  stop("Invalid SGE_TASK_ID: ", task_id)

}

donor <- tasks$donor[task_id]

bin_label <- tasks$bin[task_id]

bin_values <- c(

  "30kb" = 30000L,

  "240kb" = 240000L,

  "500kb" = 500000L,

  "5mb" = 5000000L

)

binsize <- unname(bin_values[[bin_label]])

root <- Sys.getenv("ASCAT_INPUT_ROOT", unset = "")
if (!nzchar(root)) stop("Set ASCAT_INPUT_ROOT")

output_root <- Sys.getenv("ASCAT_OUTPUT_ROOT", unset = "")
if (!nzchar(output_root)) stop("Set ASCAT_OUTPUT_ROOT")

sample_ids <- unique(tasks$donor)
resolve_sample_path <- function(sample_id, key, fallback) {
  specific <- Sys.getenv(paste0(key, "_", sample_id), unset = "")
  if (nzchar(specific)) return(specific)
  template <- Sys.getenv(key, unset = "")
  if (nzchar(template)) return(gsub("{sample}", sample_id, template, fixed = TRUE))
  fallback
}
settings <- setNames(lapply(sample_ids, function(sample_id) {
  list(
    bam_dir = resolve_sample_path(sample_id, "ASCAT_BAM_DIR", file.path(root, sample_id, "dup_marked")),
    ac_root = resolve_sample_path(sample_id, "ASCAT_CELL_AC_ROOT", file.path(root, sample_id, "bam_to_beagle_output")),
    sex = Sys.getenv(paste0("ASCAT_SEX_", sample_id), "female"),
    expected_cells = as.integer(Sys.getenv(paste0("ASCAT_NCELLS_", sample_id), "0"))
  )
}), sample_ids)

setting <- settings[[donor]]

new_phase_dir <- resolve_sample_path(
  donor,
  "ASCAT_PHASE_DIR",
  file.path(root, donor, "WGS/bam_to_beagle_output/beagle_confident_v1")
)

outdir <- file.path(

  output_root,

  paste0(

    donor, "_ascat_all_methods_", bin_label,

    "_ploidy1.5-2.5_confidentBeagle_chrXnonPAR_snpBB_NANB"

  )

)

if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)) > 0L) {

  stop("Output directory is not empty; refusing to overwrite: ", outdir)

}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

message("Donor: ", donor)

message("Bin size: ", binsize)

message("Ploidy range: 1.5-2.5")

message("Confident Beagle phasing (chr1-22, chrX): ", new_phase_dir)

message("Output: ", outdir)

allchr <- c(paste0("chr", 1:22), "chrX")

bam_files <- sort(list.files(

  setting$bam_dir,

  pattern = Sys.getenv("ASCAT_BAM_PATTERN", "\\.hg38\\.sorted\\.dedup\\.bam$"),

  full.names = TRUE

))

if (!length(bam_files)) stop("No BAM files found")

cell_ids <- sub(Sys.getenv("ASCAT_BAM_CELL_SUFFIX", "\\.dedup\\.bam$"), "", basename(bam_files))
if (anyDuplicated(cell_ids)) stop("Non-unique cell IDs derived from BAM filenames")

if (setting$expected_cells > 0L && length(bam_files) != setting$expected_cells) {

  stop(

    "Unexpected number of BAMs for ", donor,

    ": found ", length(bam_files),

    ", expected ", setting$expected_cells

  )

}

hdr <- scanBamHeader(bam_files[1])[[1]]$targets

bam_chrs <- names(hdr)

chrstring_bam <- if (

  all(sub("^chr", "", allchr) %in% bam_chrs)

) {

  ""

} else if (all(allchr %in% bam_chrs)) {

  "chr"

} else {

  stop("Chromosome naming mismatch in BAM")

}

vcf_vec <- setNames(
  vapply(allchr, function(chr) {
    chrom <- sub("^chr", "", chr)
    file.path(new_phase_dir, chr,
              paste0(donor, "_chr", chrom, ".beagleoutput.vcf.gz"))
  }, character(1)),
  allchr
)

if (!identical(names(vcf_vec), allchr)) {

  stop("Unexpected chromosome order in phasing VCFs")

}

if (any(!file.exists(vcf_vec))) {

  stop(

    "Missing Beagle VCFs: ",

    paste(vcf_vec[!file.exists(vcf_vec)], collapse=", ")

  )

}

if (any(!file.exists(paste0(vcf_vec, ".tbi")))) {

  stop(

    "Missing Beagle VCF indexes: ",

    paste(

      vcf_vec[!file.exists(paste0(vcf_vec, ".tbi"))],

      collapse=", "

    )

  )

}

message("Phasing VCFs verified: ", length(vcf_vec))

path_to_phases <- list(vcf_vec)

get_ac_vec <- function(cell_id) {

  ac_dir <- file.path(setting$ac_root, cell_id, "allele_count")

  if (!dir.exists(ac_dir)) {

    stop("Missing allele-count directory: ", ac_dir)

  }

  ac_files <- list.files(
    ac_dir,
    pattern = "_chr([0-9]+|X)\\.txt$",
    full.names = TRUE
  )
  get_chr <- function(x) {
    sub(".*_(chr([0-9]+|X))\\.txt$", "\\1", basename(x))
  }
  ac_names <- get_chr(ac_files)

  if (anyDuplicated(ac_names)) {

    stop("Duplicate allele-count chromosomes: ", cell_id)

  }

  ac_vec <- setNames(ac_files, ac_names)

  ac_vec <- ac_vec[allchr]

  if (length(ac_vec) != 23L || anyNA(ac_vec)) {

    stop("Incomplete allele counts: ", cell_id)

  }

  ac_vec

}

list_ac_counts_paths <- lapply(cell_ids, get_ac_vec)

message(

  "Input preflight passed: ",

  length(bam_files), " BAMs, ",

  length(list_ac_counts_paths), " allele-count sets, ",

  length(vcf_vec), " phase VCFs"

)

if (identical(Sys.getenv("ASCAT_PREFLIGHT_ONLY"), "1")) {

  message("PREFLIGHT PASSED: ", donor, " ", bin_label)

  quit(save="no", status=0)

}

purs <- replicate(

  length(bam_files),

  seq(0.999, 1.000, by = 0.001),

  simplify = FALSE

)

ploidies <- replicate(

  length(bam_files),

  seq(1.5, 2.5, by = 0.01),

  simplify = FALSE

)

cores <- 8L

message("Starting full ASCAT.sc analysis")

message("Number of cells: ", length(bam_files))

res <- run_sc_sequencing(

  tumour_bams = bam_files,

  sex = rep(setting$sex, length(bam_files)),

  allchr = allchr,

  chrstring_bam = chrstring_bam,

  build = "hg38",

  binsize = binsize,

  purs = purs,

  ploidies = ploidies,

  sc_filters = TRUE,

  multipcf = FALSE,

  predict_refit = FALSE,

  path_to_phases = path_to_phases,

  list_ac_counts_paths = list_ac_counts_paths,

  projectname = paste0(donor, "_", bin_label, "_noIMP"),

  MC.CORES = cores,

  outdir = outdir,

  betabinom = FALSE

)

if (is.null(res$allSolutions) || is.null(res$allProfiles)) {

  stop("Full ASCAT.sc analysis did not return the expected results")

}

message("Full ASCAT.sc analysis completed")

message("Running integrated ASCAT.sc SNP-level beta-binomial model")

# Use the getAS_CNA implementation provided by the installed ASCAT.sc package.
# Verify its beta-binomial API and output equivalence against the analysis version.

res <- getAS_CNA(

  res = res,

  path_to_phases = path_to_phases,

  list_ac_counts_paths = list_ac_counts_paths,

  purs = res$purs,

  ploidies = res$ploidies,

  chrstring = "chr",

  projectname = paste0(donor, "_", bin_label, "_noIMP_snpBB"),

  outdir = outdir,

  betabinom = TRUE,

  mc.cores = cores,

  steps = 5L,

  snp_bb_min_spacing = 1000L,

  snp_bb_include_chrX = TRUE

)

profiles <- res$allProfiles_AS

if (length(profiles) != length(bam_files)) {

  stop("Unexpected number of allele-specific profiles")

}

valid <- vapply(profiles, function(cell) {

  if (!is.list(cell)) return(FALSE)

  p <- cell$nprof.fixed

  fit <- cell$beta_bin_fit

  is.data.frame(p) &&

    all(c(

      "BAF", "BAF_beta_bin", "rho_beta_bin",

      "n_snps_beta_bin",

      "nA_best_overdispersion",

      "nB_best_overdispersion",

      "nA", "nB"

    ) %in% names(p)) &&

    is.list(fit) &&

    length(fit$rho) == 1L &&

    is.finite(fit$rho)

}, logical(1))

if (!all(valid)) {

  stop("Incomplete SNP beta-binomial results: ",

       sum(!valid), " cells")

}

output_rds <- file.path(

  outdir,

  paste0(donor, "_single_cells_ASCAT_sc_noIMP_snpBB.rds")

)

saveRDS(res, output_rds)

message("SUCCESS: ", donor, " ", bin_label)

message("Cells: ", length(profiles))

message("Valid beta-binomial profiles: ", sum(valid))

message("Saved: ", output_rds)

