args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2L || length(args) > 3L) {
  stop(
    "Usage: Rscript export_simplified_cn_profile_image.R ",
    "<input.rds> <exact_cell_name> [output_directory]"
  )
}

suppressPackageStartupMessages(library(ASCAT.sc))

if (!requireNamespace("svglite", quietly = TRUE)) {
  stop("Package svglite is required. Run install.packages('svglite')")
}

rds_path <- normalizePath(args[1], mustWork = TRUE)
cell_name <- args[2]

outdir <- if (length(args) == 3L) {
  args[3]
} else {
  dirname(rds_path)
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

res <- readRDS(rds_path)

idx <- which(names(res$allTracks) == cell_name)

if (length(idx) != 1L) {
  cat("Matching cell names:\n")
  print(grep(
    cell_name,
    names(res$allTracks),
    value = TRUE,
    fixed = TRUE
  ))
  stop("Cell name must match exactly one ASCAT cell")
}

track <- res$allTracks.processed[[idx]]
solution <- res$allSolutions[[idx]]

if (is.null(track) || is.null(solution)) {
  stop("Missing ASCAT track or solution for selected cell")
}

cell_label <- sub("_.*$", "", cell_name)

gamma_used <- if (!is.null(res$gamma)) {
  res$gamma
} else {
  1
}

cn_lwd <- 6
cn_secondary_lwd <- 2.5
separator_lwd <- 1.2
border_lwd <- 2.2

y_bottom <- -0.8
y_top <- 6.8

width_px <- 5500
height_px <- 2496
dpi <- 300

expr_text <- function(x) {
  paste(deparse(x, width.cutoff = 500L), collapse = "")
}

changes <- new.env(parent = emptyenv())

changes$cn_primary <- 0L
changes$cn_secondary <- 0L
changes$header <- 0L
changes$footer <- 0L
changes$top_label <- 0L
changes$bottom_label <- 0L
changes$vertical <- 0L

modify_plot_code <- function(expr) {

  if (!is.call(expr)) {
    return(expr)
  }

  fname <- if (is.symbol(expr[[1L]])) {
    as.character(expr[[1L]])
  } else {
    ""
  }

  if (fname == "<-" && length(expr) == 3L) {

    lhs <- expr_text(expr[[2L]])

    if (lhs == "LWD_SEGS_INTEGER") {
      expr[[3L]] <- cn_lwd
      changes$cn_primary <- changes$cn_primary + 1L
    }

    if (lhs == "LWD_SEGS_INTEGER2") {
      expr[[3L]] <- cn_secondary_lwd
      changes$cn_secondary <- changes$cn_secondary + 1L
    }
  }

  if (fname == "segments") {

    lw <- expr[["lwd"]]

    if (
      is.numeric(lw) &&
      length(lw) == 1L &&
      lw == 25
    ) {

      y <- expr_text(expr[[3L]])

      if (y == "9") {
        changes$header <- changes$header + 1L
        return(quote(invisible(NULL)))
      }

      if (y == "-1") {
        changes$footer <- changes$footer + 1L
        return(quote(invisible(NULL)))
      }
    }

    if (
      is.numeric(lw) &&
      length(lw) == 1L &&
      lw == 4
    ) {

      y_start <- expr_text(expr[[3L]])

      if (y_start %in% c("-1.25", "-0.8")) {
        expr[[3L]] <- y_bottom
        expr[[5L]] <- y_top
        expr[["lwd"]] <- separator_lwd

        changes$vertical <- changes$vertical + 1L
      }
    }
  }

  if (
    fname == "text" &&
    !is.null(expr[["labels"]]) &&
    expr_text(expr[["labels"]]) == "labels[i]"
  ) {

    y <- expr_text(expr[["y"]])

    if (y == "-1") {
      changes$bottom_label <- changes$bottom_label + 1L
      return(quote(invisible(NULL)))
    }

    if (y == "ylim[2] + 1") {
      expr[["col"]] <- "black"
      changes$top_label <- changes$top_label + 1L
    }
  }

  if (length(expr) > 1L) {
    for (j in seq.int(2L, length(expr))) {
      expr[[j]] <- modify_plot_code(expr[[j]])
    }
  }

  expr
}

plotSolution_custom <- get(
  "plotSolution",
  envir = asNamespace("ASCAT.sc")
)

body(plotSolution_custom) <- modify_plot_code(
  body(plotSolution_custom)
)

observed <- c(
  cn_primary = changes$cn_primary,
  cn_secondary = changes$cn_secondary,
  header = changes$header,
  footer = changes$footer,
  top_label = changes$top_label,
  bottom_label = changes$bottom_label,
  vertical = changes$vertical
)

expected <- c(
  cn_primary = 1L,
  cn_secondary = 1L,
  header = 1L,
  footer = 1L,
  top_label = 1L,
  bottom_label = 1L,
  vertical = 3L
)

if (!identical(observed, expected)) {
  stop(
    "Unexpected ASCAT.sc plotting structure. No figures exported.\n",
    paste(names(observed), observed, collapse = ", ")
  )
}

chromosome_lengths <- vapply(
  track$lSegs,
  function(x) {
    max(x$output$loc.end) / 1e6
  },
  numeric(1)
)

genome_width <- sum(chromosome_lengths) + 4

draw_plot <- function() {

  plotSolution_custom(
    tracksSingle = track,
    purity = solution$purity,
    ploidy = solution$ploidy,
    sol = solution,
    gamma = gamma_used,
    ismale = isTRUE(
      tolower(as.character(res$sex[idx])) == "male"
    ),
    ismedian = FALSE,
    allchr = res$chr,
    rainbowChr = FALSE,
    is_pdf = FALSE,
    ylim = c(0, 6)
  )

  segments(
    x0 = 0,
    y0 = y_top,
    x1 = genome_width,
    y1 = y_top,
    col = "black",
    lwd = border_lwd,
    xpd = NA
  )

  segments(
    x0 = 0,
    y0 = y_bottom,
    x1 = genome_width,
    y1 = y_bottom,
    col = "black",
    lwd = border_lwd,
    xpd = NA
  )

  segments(
    x0 = 0,
    y0 = y_bottom,
    x1 = 0,
    y1 = y_top,
    col = "black",
    lwd = border_lwd,
    xpd = NA
  )

  segments(
    x0 = genome_width,
    y0 = y_bottom,
    x1 = genome_width,
    y1 = y_top,
    col = "black",
    lwd = border_lwd,
    xpd = NA
  )

  title(main = paste("Cell", cell_label), line = 0.5)
}

filename_base <- paste0(
  cell_label,
  "_ASCAT_black_clean"
)

output_svg <- file.path(
  outdir,
  paste0(filename_base, ".svg")
)

output_png <- file.path(
  outdir,
  paste0(filename_base, ".png")
)

svglite::svglite(
  file = output_svg,
  width = width_px / dpi,
  height = height_px / dpi,
  bg = "white"
)

draw_plot()
dev.off()

png(
  filename = output_png,
  width = width_px,
  height = height_px,
  units = "px",
  res = dpi,
  bg = "white"
)

draw_plot()
dev.off()

for (path in c(output_svg, output_png)) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    stop("Missing or empty output: ", path)
  }
}

cat("\nCell:", cell_name, "\n")
cat("PNG:", output_png, "\n")
cat("SVG:", output_svg, "\n")
cat("Both exports verified.\n")
