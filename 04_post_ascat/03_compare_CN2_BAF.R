suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))

root <- Sys.getenv("ASCAT_INPUT_ROOT", unset = "")
if (!nzchar(root)) stop("Set ASCAT_INPUT_ROOT")

ascat_root <- Sys.getenv("ASCAT_OUTPUT_ROOT", unset = file.path(root, "ASCAT_Output"))

outdir <- file.path(
  ascat_root,
  "scWGS_CN2_BAF_comparison"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

donor <- Sys.getenv("ASCAT_SAMPLE_ID", unset = "")
if (!nzchar(donor)) stop("Set ASCAT_SAMPLE_ID")

target_cells <- trimws(strsplit(Sys.getenv("ASCAT_CELL_IDS", unset = ""), ",", fixed = TRUE)[[1]])
if (!length(target_cells) || any(!nzchar(target_cells))) stop("Set ASCAT_CELL_IDS")

bin_sizes <- c("240kb", "500kb", "5mb")

all_profiles <- list()
missing_cells <- list()

for (bin in bin_sizes) {
  
  result_dir <- file.path(
    ascat_root,
    paste0(
      donor,
      "_ascat_all_methods_",
      bin,
      "_ploidy1.5-2.5_confidentBeagle_chrXnonPAR_snpBB_NANB"
    )
  )
  
  rds_path <- file.path(
    result_dir,
    paste0(
      donor,
      "_single_cells_ASCAT_sc_noIMP_snpBB.rds"
    )
  )
  
  cat("\nChecking ", bin, "\n", sep = "")
  
  if (!file.exists(rds_path)) {
    warning("Missing ASCAT RDS: ", rds_path)
    
    missing_cells[[length(missing_cells) + 1L]] <- data.table(
      bin_size = bin,
      cell = target_cells,
      reason = "ASCAT RDS missing"
    )
    
    next
  }
  
  res <- readRDS(rds_path)
  
  track_names <- names(res$allTracks)
  
  if (is.null(track_names)) {
    stop("Unnamed allTracks in: ", rds_path)
  }
  
  for (cell_id in target_cells) {
    
    idx <- grep(
      paste0("^", cell_id, "([_-]|$)"),
      track_names,
      ignore.case = TRUE
    )
    
    if (length(idx) != 1L) {
      missing_cells[[length(missing_cells) + 1L]] <- data.table(
        bin_size = bin,
        cell = cell_id,
        reason = paste0("Expected one track, found ", length(idx))
      )
      next
    }
    
    if (
      is.null(res$allProfiles_AS) ||
      length(res$allProfiles_AS) < idx ||
      is.null(res$allProfiles_AS[[idx]]$nprof.fixed)
    ) {
      missing_cells[[length(missing_cells) + 1L]] <- data.table(
        bin_size = bin,
        cell = cell_id,
        reason = "Allele-specific profile missing"
      )
      next
    }
    
    d <- as.data.table(
      res$allProfiles_AS[[idx]]$nprof.fixed
    )
    
    required <- c(
      "chr", "startpos", "endpos",
      "total_copy_number", "BAF", "BAF_beta_bin"
    )
    
    absent <- setdiff(required, names(d))
    
    if (length(absent)) {
      missing_cells[[length(missing_cells) + 1L]] <- data.table(
        bin_size = bin,
        cell = cell_id,
        reason = paste0(
          "Missing columns: ",
          paste(absent, collapse = ", ")
        )
      )
      next
    }
    
    d <- d[
      total_copy_number == 2 &
        is.finite(BAF) &
        is.finite(BAF_beta_bin)
    ]
    
    if (!nrow(d)) {
      missing_cells[[length(missing_cells) + 1L]] <- data.table(
        bin_size = bin,
        cell = cell_id,
        reason = "No CN2 segments with both BAF estimates"
      )
      next
    }
    
    d[, cell := cell_id]
    d[, bin_size := bin]
    
    all_profiles[[length(all_profiles) + 1L]] <- copy(d)
    
    cat(
      cell_id, ": ",
      nrow(d),
      " matched CN2 segments\n",
      sep = ""
    )
  }
  
  rm(res)
  gc(verbose = FALSE)
}

if (!length(all_profiles)) {
  stop("No CN2 profiles found for the selected cells.")
}

cn2 <- rbindlist(all_profiles, fill = TRUE)

cn2[, cell := factor(cell, levels = target_cells)]
cn2[, bin_size := factor(bin_size, levels = bin_sizes)]

baf_long <- melt(
  cn2,
  id.vars = c(
    "cell", "bin_size",
    "chr", "startpos", "endpos"
  ),
  measure.vars = c("BAF", "BAF_beta_bin"),
  variable.name = "model",
  value.name = "BAF_value"
)

baf_long[
  ,
  model := factor(
    model,
    levels = c("BAF", "BAF_beta_bin"),
    labels = c("Binomial", "Beta-binomial")
  )
]

summary <- cn2[
  ,
  .(
    n_cn2_segments = .N,
    median_binomial = median(BAF),
    median_beta_binomial = median(BAF_beta_bin),
    median_absolute_deviation_binomial = median(abs(BAF - 0.5)),
    median_absolute_deviation_beta_binomial = median(abs(BAF_beta_bin - 0.5)),
    median_improvement = median(
      abs(BAF - 0.5) - abs(BAF_beta_bin - 0.5)
    ),
    proportion_segments_improved = mean(
      abs(BAF_beta_bin - 0.5) < abs(BAF - 0.5)
    )
  ),
  by = .(bin_size, cell)
][order(bin_size, cell)]

cat("\nCN2 BAF comparison summary:\n")
print(summary)

fwrite(
  cn2,
  file.path(outdir, "scWGS_matched_CN2_segments.tsv"),
  sep = "\t",
  na = "NA"
)

fwrite(
  summary,
  file.path(outdir, "scWGS_CN2_BAF_summary.tsv"),
  sep = "\t",
  na = "NA"
)

if (length(missing_cells)) {
  missing_dt <- rbindlist(missing_cells, fill = TRUE)
  
  cat("\nMissing cells / outputs:\n")
  print(missing_dt)
  
  fwrite(
    missing_dt,
    file.path(outdir, "missing_cells.tsv"),
    sep = "\t"
  )
}

model_colours <- c(
  "Binomial" = "#567D9A",
  "Beta-binomial" = "#C27A65"
)

for (bin in bin_sizes) {
  
  d <- baf_long[bin_size == bin]
  
  if (!nrow(d)) {
    cat("No data for ", bin, "\n", sep = "")
    next
  }
  
  counts <- summary[
    bin_size == bin,
    .(cell, n_cn2_segments)
  ]
  
  labels <- setNames(
    paste0(
      as.character(counts$cell),
      "\n(n=", counts$n_cn2_segments, ")"
    ),
    as.character(counts$cell)
  )
  
  p <- ggplot(
    d,
    aes(
      x = cell,
      y = BAF_value,
      fill = model
    )
  ) +
    geom_hline(
      yintercept = 0.5,
      linetype = "dashed",
      colour = "grey55",
      linewidth = 0.5
    ) +
    geom_boxplot(
      position = position_dodge(width = 0.75),
      width = 0.62,
      outlier.shape = NA,
      linewidth = 0.55
    ) +
    geom_point(
      aes(colour = model),
      position = position_jitterdodge(
        jitter.width = 0.09,
        dodge.width = 0.75,
        seed = 42
      ),
      size = 0.9,
      alpha = 0.3
    ) +
    scale_fill_manual(values = model_colours) +
    scale_colour_manual(values = model_colours) +
    scale_x_discrete(labels = labels, drop = FALSE) +
    scale_y_continuous(
      limits = c(0.5, NA),
      expand = expansion(mult = c(0.01, 0.05))
    ) +
    labs(
      title = paste0(
        "Single-cell — CN2 BAF comparison (",
        bin, ")"
      ),
      subtitle = "Each point represents one CN2 segment; matched segments for both models",
      x = NULL,
      y = "BAF",
      fill = NULL
    ) +
    guides(colour = "none") +
    theme_classic(base_size = 13) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold", size = 14),
      axis.text.x = element_text(size = 10),
      panel.grid = element_blank()
    )
  
  print(p)
  
  ggsave(
    filename = file.path(
      outdir,
      paste0("scWGS_CN2_BAF_", bin, ".pdf")
    ),
    plot = p,
    width = 11,
    height = 5.5
  )
}

cat(
  "\nFinished. Results saved to:\n",
  outdir,
  "\n",
  sep = ""
)