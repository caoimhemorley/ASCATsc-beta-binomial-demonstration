suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))

root <- Sys.getenv("ASCAT_INPUT_ROOT", unset = "")
if (!nzchar(root)) stop("Set ASCAT_INPUT_ROOT")

donor <- Sys.getenv("ASCAT_SAMPLE_ID", unset = "")
if (!nzchar(donor)) stop("Set ASCAT_SAMPLE_ID")
cell <- Sys.getenv("ASCAT_CELL_ID", unset = "")
if (!nzchar(cell)) stop("Set ASCAT_CELL_ID")
loss_selection <- "manual"
manual_chromosome <- Sys.getenv("ASCAT_REGION_CHROM", "chr10")
manual_start_mb <- as.numeric(Sys.getenv("ASCAT_REGION_START_MB", "6"))
manual_end_mb <- as.numeric(Sys.getenv("ASCAT_REGION_END_MB", "16"))
flank_mb <- 5
minimum_window_mb <- 25
bin_size <- "240kb"

point_colour <- "#666666"
line_colour <- "#222222"

loss_shade <- "#DCE8F0"
gain_shade <- "#EAD7CC"

loss_point_colour <- "#648BA4"
gain_point_colour <- "#A96855"

shade_alpha <- 0.50
axis_title_size <- 17       # Axis titles (pt)
axis_tick_size <- 14        # Axis tick labels (pt)
plot_title_size <- 20       # Main title (pt)
panel_border_width <- 0.9   # Border around each panel
dashed_line_width <- 0.8    # All dashed reference lines

logr_ylim <- c(-1.5, 1.1)
cn_ylim <- c(0, 3.5)
allele_ylim <- c(0, 1)
model_baf_ylim <- c(0.5, 1.0)

shade_colours <- c(
  Loss = loss_shade,
  Gain = gain_shade
)

outdir <- file.path(
  Sys.getenv("ASCAT_OUTPUT_ROOT", unset = file.path(root, "ASCAT_Output")),
  "regional_haplotype_BAF_plots"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

ascat_dir <- file.path(
  Sys.getenv("ASCAT_OUTPUT_ROOT", unset = file.path(root, "ASCAT_Output")),
  paste0(
    donor,
    "_ascat_all_methods_",
    bin_size,
    "_ploidy1.5-2.5_confidentBeagle_chrXnonPAR_snpBB_NANB"
  )
)

rds_path <- file.path(
  ascat_dir,
  paste0(donor, "_single_cells_ASCAT_sc_noIMP_snpBB.rds")
)

if (!file.exists(rds_path)) {
  stop("Missing ASCAT RDS: ", rds_path)
}

cat("Reading ASCAT results:\n", rds_path, "\n")

res <- readRDS(rds_path)

track_names <- names(res$allTracks)

idx <- grep(
  paste0("^", cell, "([_-]|$)"),
  track_names,
  ignore.case = TRUE
)

if (length(idx) != 1L) {
  stop(
    "Could not uniquely identify ASCAT track for ",
    cell,
    ". Matches: ",
    paste(track_names[idx], collapse = ", ")
  )
}

selected_track <- track_names[idx]

all_profile <- as.data.table(copy(res$allProfiles_AS[[idx]]$nprof.fixed))
required_loss_columns <- c("chr", "startpos", "endpos", "total_copy_number")
if (length(setdiff(required_loss_columns, names(all_profile)))) {
  stop("ASCAT profile is missing columns required for loss selection.")
}

loss_candidates <- all_profile[
  is.finite(total_copy_number) & total_copy_number < 2 &
    sub("^chr", "", as.character(chr), ignore.case = TRUE) %in% as.character(1:22) &
    is.finite(startpos) & is.finite(endpos) & endpos > startpos
]
loss_candidates[, width_bp := endpos - startpos]
setorder(loss_candidates, -width_bp)

if (loss_selection == "largest") {
  if (!nrow(loss_candidates)) {
    stop("No autosomal CN<2 segments in ", donor, " ", cell,
         ". Change cell or choose loss_selection = 'manual'.")
  }
  target <- loss_candidates[1L]
  chromosome <- paste0("chr", sub("^chr", "", as.character(target$chr), ignore.case = TRUE))
  midpoint_mb <- (target$startpos + target$endpos) / 2e6
  window_width_mb <- max(minimum_window_mb, target$width_bp / 1e6 + 2 * flank_mb)
  start_mb <- max(0, floor(midpoint_mb - window_width_mb / 2))
  end_mb <- ceiling(midpoint_mb + window_width_mb / 2)
  cat("Selected widest autosomal loss: ", chromosome, ":",
      target$startpos, "-", target$endpos, " CN=", target$total_copy_number, "\n", sep = "")
} else if (loss_selection == "manual") {
  chromosome <- manual_chromosome
  start_mb <- manual_start_mb
  end_mb <- manual_end_mb
} else {
  stop("loss_selection must be 'largest' or 'manual'.")
}

start_bp <- start_mb * 1e6
end_bp <- end_mb * 1e6
chr_number <- sub("^chr", "", chromosome, ignore.case = TRUE)
if (start_bp >= end_bp) stop("start_mb must be smaller than end_mb.")
cat("Plotting ", donor, " ", cell, " ", chromosome, ":", start_mb,
    "-", end_mb, " Mb\n", sep = "")

processed_idx <- match(
  selected_track,
  names(res$allTracks.processed)
)

if (is.na(processed_idx)) {
  stop("Could not match processed ASCAT track: ", selected_track)
}

processed <- res$allTracks.processed[[processed_idx]]

chr_idx <- match(chromosome, names(processed$lCTS))

if (is.na(chr_idx)) {
  stop("Chromosome not found in processed ASCAT tracks: ", chromosome)
}

if (is.null(processed$lCTS[[chromosome]])) {
  stop("No processed logR bins for ", chromosome)
}

logr <- as.data.table(
  copy(processed$lCTS[[chromosome]])
)

required_logr <- c("start", "end", "smoothed")

missing_logr <- setdiff(required_logr, names(logr))

if (length(missing_logr)) {
  stop(
    "Missing logR columns: ",
    paste(missing_logr, collapse = ", ")
  )
}

logr <- logr[
  end >= start_bp &
    start <= end_bp &
    is.finite(smoothed)
]

if ("space" %in% names(logr)) {
  logr <- logr[
    sub("^chr", "", as.character(space), ignore.case = TRUE) == chr_number
  ]
}

if (!nrow(logr)) {
  stop("No ASCAT logR bins in selected region.")
}

setorder(logr, start, end)

logr[, ascat_bin_id := .I]

logr[
  ,
  `:=`(
    position_mb = (start + end) / 2e6,
    logR = smoothed,
    plot_start = pmax(start, start_bp) / 1e6,
    plot_end = pmin(end, end_bp) / 1e6
  )
]

if (
  length(processed$lSegs) < chr_idx ||
  is.null(processed$lSegs[[chr_idx]]$output)
) {
  stop("No CBS logR segmentation for ", chromosome)
}

logr_segments <- as.data.table(
  copy(processed$lSegs[[chr_idx]]$output)
)

required_segments <- c(
  "chrom",
  "loc.start",
  "loc.end",
  "seg.mean"
)

missing_segments <- setdiff(
  required_segments,
  names(logr_segments)
)

if (length(missing_segments)) {
  stop(
    "Missing logR segmentation columns: ",
    paste(missing_segments, collapse = ", ")
  )
}

logr_segments <- logr_segments[
  sub("^chr", "", as.character(chrom), ignore.case = TRUE) == chr_number &
    loc.end >= start_bp &
    loc.start <= end_bp &
    is.finite(seg.mean)
]

logr_segments[
  ,
  `:=`(
    plot_start = pmax(loc.start, start_bp) / 1e6,
    plot_end = pmin(loc.end, end_bp) / 1e6
  )
]

logr_segments <- logr_segments[
  plot_end > plot_start
]

setorder(logr_segments, loc.start)

profile <- as.data.table(
  copy(res$allProfiles_AS[[idx]]$nprof.fixed)
)

required_profile <- c(
  "chr",
  "startpos",
  "endpos",
  "total_copy_number",
  "BAF",
  "BAF_beta_bin"
)

missing_profile <- setdiff(
  required_profile,
  names(profile)
)

if (length(missing_profile)) {
  stop(
    "Missing ASCAT profile columns: ",
    paste(missing_profile, collapse = ", ")
  )
}

profile <- profile[
  sub("^chr", "", as.character(chr), ignore.case = TRUE) == chr_number &
    endpos >= start_bp &
    startpos <= end_bp
]

profile[
  ,
  `:=`(
    plot_start = pmax(startpos, start_bp) / 1e6,
    plot_end = pmin(endpos, end_bp) / 1e6
  )
]

profile <- profile[
  plot_end > plot_start
]

if (!nrow(profile)) {
  stop("No ASCAT CN/BAF segments overlap region.")
}

setorder(profile, startpos, endpos)

profile[, segment_id := .I]

profile[
  ,
  cn_state := fifelse(
    is.na(total_copy_number),
    "Unknown",
    fifelse(
      total_copy_number < 2,
      "Loss",
      fifelse(
        total_copy_number > 2,
        "Gain",
        "CN2"
      )
    )
  )
]

profile_plot <- copy(profile)

profile_plot[
  ,
  `:=`(
    display_start = plot_start,
    display_end = plot_end
  )
]

if (nrow(profile_plot) > 1L) {
  for (i in seq_len(nrow(profile_plot) - 1L)) {
    midpoint <- (
      profile_plot$plot_end[i] +
        profile_plot$plot_start[i + 1L]
    ) / 2
    
    profile_plot$display_end[i] <- midpoint
    profile_plot$display_start[i + 1L] <- midpoint
  }
}

profile_plot[
  ,
  `:=`(
    display_start = pmax(display_start, start_mb),
    display_end = pmin(display_end, end_mb)
  )
]

non_cn2_plot <- profile_plot[
  cn_state == "Loss",
  .(
    segment_id,
    xmin = display_start,
    xmax = display_end,
    ymin = -Inf,
    ymax = Inf,
    cn_state
  )
]

logr_segments_plot <- copy(logr_segments)

logr_segments_plot[
  ,
  `:=`(
    display_start = plot_start,
    display_end = plot_end
  )
]

if (nrow(logr_segments_plot) > 1L) {
  for (i in seq_len(nrow(logr_segments_plot) - 1L)) {
    midpoint <- (
      logr_segments_plot$plot_end[i] +
        logr_segments_plot$plot_start[i + 1L]
    ) / 2
    
    logr_segments_plot$display_end[i] <- midpoint
    logr_segments_plot$display_start[i + 1L] <- midpoint
  }
}

logr_segments_plot[
  ,
  `:=`(
    display_start = pmax(display_start, start_mb),
    display_end = pmin(display_end, end_mb)
  )
]

ac_root <- Sys.getenv(
  "ASCAT_PLOT_AC_ROOT",
  unset = file.path(root, donor, "allele_counts_single_cells")
)

folders <- list.dirs(
  ac_root,
  recursive = FALSE,
  full.names = FALSE
)

matches <- folders[
  grepl(
    paste0("^", cell, "([^0-9]|$)"),
    folders,
    ignore.case = TRUE
  )
]

if (length(matches) != 1L) {
  stop(
    "Expected exactly one allele-count folder for ",
    cell,
    ". Found: ",
    paste(matches, collapse = ", ")
  )
}

ac_dir <- file.path(
  ac_root,
  matches,
  "allele_count"
)

ac_files <- list.files(
  ac_dir,
  pattern = paste0("_chr", chr_number, "\\.txt$"),
  full.names = TRUE
)

if (length(ac_files) != 1L) {
  stop("Expected exactly one allele-count file for ", chromosome)
}

cat("\nReading allele counts:\n", ac_files, "\n")

ac <- fread(ac_files)

if ("#CHR" %in% names(ac)) {
  setnames(ac, "#CHR", "chr")
}

required_ac <- c(
  "chr",
  "POS",
  "Count_A",
  "Count_C",
  "Count_G",
  "Count_T"
)

missing_ac <- setdiff(required_ac, names(ac))

if (length(missing_ac)) {
  stop(
    "Missing allele-count columns: ",
    paste(missing_ac, collapse = ", ")
  )
}

ac <- ac[
  POS >= start_bp &
    POS <= end_bp &
    sub("^chr", "", as.character(chr), ignore.case = TRUE) == chr_number
]

plot_phase_dir <- Sys.getenv(
  "ASCAT_PLOT_PHASE_DIR",
  unset = file.path(root, donor, "bam_to_beagle_output", "beagle_confident_v1")
)
vcf_path <- file.path(
  plot_phase_dir,
  paste0("chr", chr_number),
  paste0(donor, "_chr", chr_number, ".beagleoutput.vcf.gz")
)

if (!file.exists(vcf_path)) {
  stop("Missing phased VCF: ", vcf_path)
}

cat("Reading phased VCF:\n", vcf_path, "\n")

vcf <- fread(
  cmd = paste("gzip -dc", shQuote(vcf_path)),
  skip = "#CHROM",
  header = TRUE,
  select = c(1, 2, 4, 5, 10)
)

setnames(
  vcf,
  names(vcf),
  c("chr", "POS", "REF", "ALT", "GT")
)

vcf <- vcf[
  POS >= start_bp &
    POS <= end_bp &
    sub("^chr", "", as.character(chr), ignore.case = TRUE) == chr_number
]

vcf[, GT := sub(":.*$", "", GT)]

vcf <- vcf[
  GT %in% c("0|1", "1|0") &
    REF %in% c("A", "C", "G", "T") &
    ALT %in% c("A", "C", "G", "T")
]

vcf[
  ,
  `:=`(
    H1_allele = fifelse(GT == "0|1", REF, ALT),
    H2_allele = fifelse(GT == "0|1", ALT, REF)
  )
]

if (anyDuplicated(vcf$POS)) {
  warning("Duplicate phased VCF positions; retaining first.")
  vcf <- unique(vcf, by = "POS")
}

if (anyDuplicated(ac$POS)) {
  warning("Duplicate allele-count positions; aggregating counts.")
  
  ac <- ac[
    ,
    lapply(.SD, sum, na.rm = TRUE),
    by = .(POS),
    .SDcols = c(
      "Count_A",
      "Count_C",
      "Count_G",
      "Count_T"
    )
  ]
}

snp <- merge(
  ac,
  vcf[
    ,
    .(
      POS,
      REF,
      ALT,
      GT,
      H1_allele,
      H2_allele
    )
  ],
  by = "POS"
)

count_allele <- function(d, allele) {
  ans <- numeric(nrow(d))
  
  for (base in c("A", "C", "G", "T")) {
    ix <- which(allele == base)
    
    if (length(ix)) {
      ans[ix] <- d[[paste0("Count_", base)]][ix]
    }
  }
  
  ans
}

snp[, H1_count := count_allele(snp, H1_allele)]
snp[, H2_count := count_allele(snp, H2_allele)]
snp[, phased_depth := H1_count + H2_count]

snp <- snp[
  is.finite(phased_depth) &
    phased_depth > 0
]

if (!nrow(snp)) {
  stop("No phased SNPs with counts in selected region.")
}

snp[
  ,
  `:=`(
    position_mb = POS / 1e6,
    H1_fraction = H1_count / phased_depth
  )
]

ascat_intervals <- logr[
  ,
  .(
    ascat_bin_id,
    start,
    end
  )
]

setkey(ascat_intervals, start, end)

snp_positions <- snp[
  ,
  .(
    snp_id = .I,
    POS,
    snp_start = POS,
    snp_end = POS
  )
]

setkey(snp_positions, snp_start, snp_end)

overlap_hits <- foverlaps(
  snp_positions,
  ascat_intervals,
  by.x = c("snp_start", "snp_end"),
  by.y = c("start", "end"),
  type = "within",
  nomatch = 0L
)

hit_counts <- overlap_hits[
  ,
  .N,
  by = snp_id
]

ambiguous_ids <- hit_counts[N > 1, snp_id]

unassigned_ids <- setdiff(
  seq_len(nrow(snp)),
  unique(overlap_hits$snp_id)
)

cat("\nSNP-to-ASCAT-bin assignment:\n")
cat("Total phased SNPs:", nrow(snp), "\n")
cat("Assigned SNPs:", uniqueN(overlap_hits$snp_id), "\n")
cat("Unassigned SNPs:", length(unassigned_ids), "\n")
cat("Ambiguous SNPs:", length(ambiguous_ids), "\n")

if (length(ambiguous_ids)) {
  warning(
    length(ambiguous_ids),
    " SNPs overlap multiple ASCAT bins; excluding them."
  )
}

unique_hits <- overlap_hits[
  !snp_id %in% ambiguous_ids,
  .(
    snp_id,
    ascat_bin_id
  )
]

snp[, ascat_bin_id := NA_integer_]

snp[
  unique_hits$snp_id,
  ascat_bin_id := unique_hits$ascat_bin_id
]

assigned_snp <- snp[
  !is.na(ascat_bin_id)
]

if (!nrow(assigned_snp)) {
  stop("No phased SNPs uniquely assigned to ASCAT bins.")
}

binned_counts <- assigned_snp[
  ,
  .(
    H1_count = sum(H1_count),
    H2_count = sum(H2_count),
    total_depth = sum(phased_depth),
    n_snps = .N,
    median_snp_depth = median(phased_depth)
  ),
  by = ascat_bin_id
]

binned <- merge(
  logr[
    ,
    .(
      ascat_bin_id,
      start,
      end,
      position_mb,
      logR
    )
  ],
  binned_counts,
  by = "ascat_bin_id",
  all.x = TRUE,
  sort = FALSE
)

setorder(binned, ascat_bin_id)

binned[
  ,
  `:=`(
    H1_fraction = fifelse(
      !is.na(total_depth) & total_depth > 0,
      H1_count / total_depth,
      NA_real_
    ),
    H2_fraction = fifelse(
      !is.na(total_depth) & total_depth > 0,
      H2_count / total_depth,
      NA_real_
    )
  )
]

binned[
  ,
  `:=`(
    observed_folded_BAF = pmax(
      H1_fraction,
      H2_fraction
    ),
    bin_start_mb = start / 1e6,
    bin_end_mb = end / 1e6,
    bin_width_bp = end - start + 1
  )
]

observed_long <- melt(
  binned[
    is.finite(H1_fraction) &
      is.finite(H2_fraction)
  ],
  id.vars = c(
    "ascat_bin_id",
    "position_mb",
    "start",
    "end",
    "n_snps",
    "total_depth"
  ),
  measure.vars = c(
    "H1_fraction",
    "H2_fraction"
  ),
  variable.name = "haplotype",
  value.name = "fraction"
)

observed_long[
  ,
  cnv_state := "CN2"
]

if (nrow(non_cn2_plot)) {
  for (i in seq_len(nrow(non_cn2_plot))) {
    cnv_start_bp <- non_cn2_plot$xmin[i] * 1e6
    cnv_end_bp <- non_cn2_plot$xmax[i] * 1e6
    
    observed_long[
      end >= cnv_start_bp &
        start <= cnv_end_bp,
      cnv_state := non_cn2_plot$cn_state[i]
    ]
  }
}

observed_long[
  ,
  point_colour := fcase(
    cnv_state == "Loss", loss_point_colour,
    default = point_colour
  )
]

colour_for_position <- function(x) {
  result <- rep(point_colour, length(x))
  if (nrow(non_cn2_plot)) {
    for (i in seq_len(nrow(non_cn2_plot))) {
      inside <- x >= non_cn2_plot$xmin[i] & x <= non_cn2_plot$xmax[i]
      result[inside] <- if (non_cn2_plot$cn_state[i] == "Gain") {
        gain_point_colour
      } else {
        loss_point_colour
      }
    }
  }
  result
}

logr[, track_colour := colour_for_position(position_mb)]

colour_logr_segments <- function(d) {
  if (!nrow(d)) return(d[, track_colour := character()])
  boundaries <- sort(unique(c(non_cn2_plot$xmin, non_cn2_plot$xmax)))
  pieces <- lapply(seq_len(nrow(d)), function(i) {
    a <- d$display_start[i]
    b <- d$display_end[i]
    cuts <- sort(unique(c(a, boundaries[boundaries > a & boundaries < b], b)))
    if (length(cuts) < 2L) return(NULL)
    out <- d[rep(i, length(cuts) - 1L)]
    out[, `:=`(display_start = head(cuts, -1L), display_end = tail(cuts, -1L))]
    out[, track_colour := colour_for_position((display_start + display_end) / 2)]
    out
  })
  rbindlist(pieces, use.names = TRUE, fill = TRUE)
}

logr_segments_coloured <- colour_logr_segments(logr_segments_plot)
profile_plot[, track_colour := fifelse(
  cn_state == "Loss", loss_point_colour, line_colour
)]

make_shading_layer <- function() {
  if (!nrow(non_cn2_plot)) {
    return(NULL)
  }
  
  geom_rect(
    data = non_cn2_plot,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = cn_state
    ),
    inherit.aes = FALSE,
    alpha = shade_alpha,
    colour = NA,
    show.legend = FALSE
  )
}

base_theme <- theme_classic(base_size = 13) +
  theme(
    panel.background = element_rect(
      fill = "white",
      colour = NA
    ),
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),
    panel.border = element_rect(
      colour = "#444444",
      fill = NA,
      linewidth = panel_border_width
    ),
    axis.line = element_blank(),
    axis.ticks = element_line(
      colour = "#666666",
      linewidth = 0.4
    ),
    axis.text = element_text(
      colour = "#222222",
      size = axis_tick_size,
      face = "plain"
    ),
    axis.title = element_text(
      colour = "#111111",
      size = axis_title_size,
      face = "bold"
    ),
    axis.title.y = element_text(
      margin = margin(r = 9)
    ),
    axis.title.x = element_text(
      margin = margin(t = 10)
    ),
    plot.title = element_text(
      face = "bold",
      colour = "#111111",
      size = plot_title_size,
      margin = margin(b = 5)
    ),
    plot.margin = margin(2, 10, 2, 10),
    legend.position = "none"
  )

x_breaks <- pretty(c(start_mb, end_mb), n = 5)

x_breaks <- x_breaks[
  x_breaks >= start_mb &
    x_breaks <= end_mb
]

x_scale <- scale_x_continuous(
  limits = c(start_mb, end_mb),
  breaks = x_breaks,
  expand = c(0, 0)
)

hide_x <- theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.x = element_blank()
)

p_logr <- ggplot() +
  make_shading_layer() +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "#A0A0A0",
    linewidth = dashed_line_width
  ) +
  geom_point(
    data = logr,
    mapping = aes(
      x = position_mb,
      y = logR,
      colour = track_colour
    ),
    size = 2.1,
    alpha = 0.8
  ) +
  geom_segment(
    data = logr_segments_coloured,
    aes(
      x = display_start,
      xend = display_end,
      y = seg.mean,
      yend = seg.mean,
      colour = track_colour
    ),
    linewidth = 1.5,
    lineend = "butt"
  ) +
  scale_colour_identity() +
  scale_fill_manual(values = shade_colours) +
  x_scale +
  scale_y_continuous(
    breaks = c(-1.5, -1, -0.5, 0, 0.5, 1)
  ) +
  coord_cartesian(ylim = logr_ylim) +
  labs(
    title = paste0(
      cell, " | ",
      chromosome, ":",
      start_mb, "–", end_mb, " Mb"
    ),
    x = NULL,
    y = "logR"
  ) +
  base_theme +
  hide_x

p_cn <- ggplot() +
  make_shading_layer() +
  geom_hline(
    yintercept = 2,
    linetype = "dashed",
    colour = "#A0A0A0",
    linewidth = dashed_line_width
  ) +
  geom_segment(
    data = profile_plot[is.finite(total_copy_number)],
    aes(
      x = display_start,
      xend = display_end,
      y = total_copy_number,
      yend = total_copy_number,
      colour = track_colour
    ),
    linewidth = 1.8,
    lineend = "butt"
  ) +
  scale_colour_identity() +
  scale_fill_manual(values = shade_colours) +
  x_scale +
  scale_y_continuous(
    breaks = c(0, 1, 2, 3)
  ) +
  coord_cartesian(ylim = cn_ylim) +
  labs(
    x = NULL,
    y = "Copy number"
  ) +
  base_theme +
  hide_x

p_observed <- ggplot() +
  make_shading_layer() +
  geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    colour = "#A0A0A0",
    linewidth = dashed_line_width
  ) +
  geom_point(
    data = observed_long,
    aes(
      x = position_mb,
      y = fraction,
      colour = point_colour
    ),
    size = 2.3,
    alpha = 0.85
  ) +
  scale_colour_identity() +
  scale_fill_manual(values = shade_colours) +
  x_scale +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1)
  ) +
  coord_cartesian(ylim = allele_ylim) +
  labs(
    x = NULL,
    y = "Observed AF"
  ) +
  base_theme +
  hide_x

make_baf_plot <- function(model_name, value_column, show_x = FALSE) {
  d <- copy(profile_plot)
  
  d[, plot_baf := get(value_column)]
  d <- d[is.finite(plot_baf)]
  
  p <- ggplot() +
    make_shading_layer() +
    geom_hline(
      yintercept = 0.5,
      linetype = "dashed",
      colour = "#A0A0A0",
      linewidth = dashed_line_width
    ) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      colour = loss_point_colour,
      linewidth = dashed_line_width
    ) +
    geom_segment(
      data = d,
      aes(
        x = display_start,
        xend = display_end,
        y = plot_baf,
        yend = plot_baf,
        colour = track_colour
      ),
      linewidth = 1.8,
      lineend = "butt"
    ) +
    scale_colour_identity() +
    scale_fill_manual(values = shade_colours) +
    x_scale +
    scale_y_continuous(
      breaks = c(0.5, 0.75, 1.0),
      labels = function(x) sprintf("%.2f", x)
    ) +
    coord_cartesian(ylim = model_baf_ylim) +
    labs(
      x = if (show_x) {
        paste0(chromosome, " position (Mb)")
      } else {
        NULL
      },
      y = model_name
    ) +
    base_theme
  
  if (!show_x) {
    p <- p + hide_x
  }
  
  p
}

p_binomial <- make_baf_plot(
  model_name = "BAF\n(binomial)",
  value_column = "BAF"
)

p_beta <- make_baf_plot(
  model_name = "BAF\n(beta-binomial)",
  value_column = "BAF_beta_bin",
  show_x = TRUE
)

combined <- (
  p_logr /
    p_cn /
    p_observed /
    p_binomial /
    p_beta
) +
  plot_layout(
    heights = c(1.3, 0.65, 1.5, 1.0, 1.0)
  )

prefix <- paste0(
  donor, "_", cell, "_",
  chromosome, "_",
  start_mb, "-", end_mb, "Mb_",
  bin_size,
  "_regional_CN1_loss"
)

pdf_path <- file.path(
  outdir,
  paste0(prefix, ".pdf")
)

svg_path <- file.path(
  outdir,
  paste0(prefix, ".svg")
)

png_path <- file.path(
  outdir,
  paste0(prefix, ".png")
)

ggsave(
  pdf_path,
  combined,
  width = 12,
  height = 11.5
)

ggsave(
  png_path,
  combined,
  width = 12,
  height = 11.5,
  dpi = 300
)

ggsave(
  svg_path,
  combined,
  width = 12,
  height = 11.5,
  device = "svg"
)

print(combined)

fwrite(
  snp,
  file.path(
    outdir,
    paste0(prefix, "_phased_SNP_counts.tsv")
  ),
  sep = "\t"
)

fwrite(
  binned,
  file.path(
    outdir,
    paste0(prefix, "_matched_allele_counts.tsv")
  ),
  sep = "\t"
)

fwrite(
  logr,
  file.path(
    outdir,
    paste0(prefix, "_ASCAT_logR_bins.tsv")
  ),
  sep = "\t"
)

fwrite(
  logr_segments,
  file.path(
    outdir,
    paste0(prefix, "_logR_CBS_segments.tsv")
  ),
  sep = "\t"
)

fwrite(
  profile_plot,
  file.path(
    outdir,
    paste0(prefix, "_ASCAT_CN_BAF_segments.tsv")
  ),
  sep = "\t"
)

fwrite(
  observed_long,
  file.path(
    outdir,
    paste0(prefix, "_observed_allele_fractions.tsv")
  ),
  sep = "\t"
)

cat("\nSaved updated figure:\n")
cat(pdf_path, "\n")
cat(svg_path, "\n")
cat(png_path, "\n")

cat("\nASCAT fitted CN and BAF values:\n")

print(
  profile_plot[
    ,
    .(
      chr,
      startpos,
      endpos,
      total_copy_number,
      BAF,
      BAF_beta_bin,
      display_start,
      display_end
    )
  ]
)

cat("\nAllele-count summary:\n")
cat("Phased SNPs:", nrow(snp), "\n")
cat("ASCAT logR bins:", nrow(logr), "\n")
cat(
  "Bins with observed allele fractions:",
  sum(is.finite(binned$H1_fraction)),
  "\n"
)
cat("Unassigned SNPs:", length(unassigned_ids), "\n")
cat("Ambiguously assigned SNPs:", length(ambiguous_ids), "\n")