read_meth_file <- function(file, workers = 8, hdf5 = NULL, col_names = NULL) {
  bpparam <- BiocParallel::MulticoreParam(progressbar = TRUE, workers = workers)

  if (is.null(hdf5)) {
    hdf5 <- tempfile()
  } else if (file.exists(hdf5)) {
    return(HDF5Array::loadHDF5SummarizedExperiment(dir = hdf5))
  }

  if (is.null(col_names)) {
    col_names <- data.frame(row.names = seq_along(file))
  } else {
    col_names <- data.frame(row.names = col_names)
  }

  methfile <- read.bismark(files = file,
                           colData = col_names,
                           rmZeroCov = TRUE,
                           verbose = TRUE,
                           BPPARAM = bpparam,
                           BACKEND = "HDF5Array",
                           dir = hdf5,
                           nThread = 1)
  methfile
}

write_segmented_methylation <- function(segments, file, col_names = FALSE) {
  write.table(segments, file,
              append = TRUE,
              sep = "\t",
              quote = FALSE,
              row.names = FALSE, col.names = col_names)
}

write_segmented_bsseq <- function(cov, met, dir) {
  sample_names <- colnames(cov)[4:ncol(cov)]

  pos <- cov[, c(1, 2, 3)]
  colnames(pos) <- c("chr", "start", "end")
  cov_matrix <- as.matrix(cov[, 4:ncol(cov)])
  met_matrix <- as.matrix(met[, 4:ncol(met)])

  colnames(cov_matrix) <- sample_names
  colnames(met_matrix) <- sample_names

  bsseq_obj <- BSseq(gr = GRanges(pos), M = met_matrix, Cov = cov_matrix)
  saveHDF5SummarizedExperiment(bsseq_obj, dir = dir)
}

binom_seg_bp <- function(meth, cov, offset, lambda) {
  segments <- fusedsegmentation(meth,
                                C = cov,
                                lambda2 = lambda,
                                objective = "binomial",
                                l = 0)
  segments$end[1:(nrow(segments) - 1)] + offset
}

combine_two_bp_sets <- function(x, y){
  combine_two_bp_sets_(as.integer(x), as.integer(y))
}

segment_block_reduce <- function(data, rowp, colp, lambda){
  grid <- ArbitraryArrayGrid(tickmarks = list(rowp, colp))

  number_of_rows <- length(rowp)
  summary_list <- bplapply(seq_along(grid),
    function(b) {
      lambda_i <- (b-1) %/% number_of_rows + 1
      viewport <- grid[[b]]
      methylation <- rowSums(read_block(data@assays@data$M, viewport))
      coverage <- rowSums(read_block(data@assays@data$Cov, viewport))
      if (length(methylation) == 1)
        return(realize(matrix(end(ranges(viewport))[1], ncol = 1),
                       BACKEND = "HDF5Array"))
      coverage[coverage == 0] <- 1
      bps <- binom_seg_bp(methylation,
                          coverage,
                          start(ranges(viewport))[1] - 1,
                          lambda[lambda_i])
      realize(matrix(bps, ncol = 1), BACKEND = "HDF5Array")
    }
  )
  Reduce(combine_two_bp_sets, summary_list)
}

format_segments <- function(ranges, chrom_end_points, breakpoints) {
  segments <-
    data.frame(chrom = seqnames(ranges)[breakpoints],
               start = c(start(ranges)[1],
                   1 + start(ranges)[breakpoints[1:(length(breakpoints) - 1)]]),
               end = start(ranges)[breakpoints])

  segments_chrom_splits <-
    cumsum(
      Rle(as.factor(segments$chrom))@lengths
    )[1:(length(chrom_end_points) - 1)] + 1

  if (length(chrom_end_points) > 1) {
    segments[segments_chrom_splits, 2] <-
      start(ranges)[chrom_end_points[1:(length(chrom_end_points) - 1)] + 1]
  }
  segments
}

aggregate_meth_data <- function(data, segments, meth_colnames, type) {
  out <- getCoverage(data, regions = makeGRangesFromDataFrame(segments),
                     type = type,
                     what = "perRegionTotal")
  colnames(out) <- meth_colnames
  cbind(segments, out)
}

cluster_methylation <- function(infiles,
                                hdf5 = NULL,
                                regions = NULL,
                                chrom = NULL,
                                start = NULL,
                                end = NULL,
                                col_names = NULL,
                                binwidth = NULL) {
  if (!is.null(regions)) {
    if (class(regions) != "GRanges"){
      stop("regions should be GRanges.")
    }
    cat("Regions are given. Clustering will be done over the given regions!\n")
  } else if (is.null(start) || is.null(end) || is.null(chrom)){
    stop(paste("Either user-defined regions should be given through 'regions'",
               "argument or a genomic region should be deliminated through",
               "'chrom', 'start', and 'end' arguments."))
  } else {
    cat(paste0("Using genomic region: ", chrom, ": ", start, " - ", end, "\n"))
    breakpoints <- seq(start - 1, end, by = binwidth)
    if (breakpoints[length(breakpoints)] != end)
      breakpoints <- c(breakpoints, end)

    region_start <- breakpoints[1:(length(breakpoints) - 1)] + 1
    region_end <- breakpoints[2:length(breakpoints)]

    regions <- GRanges(data.frame(chr = chrom,
                                  start = region_start,
                                  end = region_end))
  }

  cat("Loading the data...\n")
  methdata <- read_meth_file(infiles, hdf5 = hdf5, col_names = col_names)

  met <- getCoverage(methdata,
                     regions = regions,
                     type = "M",
                     what = "perRegionTotal")

  cov <- getCoverage(methdata,
                     regions = regions,
                     type = "Cov",
                     what = "perRegionTotal")

  usable_sites <- which(apply(cov, 1, function(x) !any(is.na(x) | (x == 0))))
  if (length(usable_sites) == 0) {
    stop("Not enough sites with coverage. Please specify a bigger region.")
  }

  cat(paste0("The total number of sites with coverage: ",
             length(usable_sites), "\n"))

  beta <- met[usable_sites, ] / cov[usable_sites, ]

  beta <- t(beta)

  cat("Clustering the samples...\n")
  cat("Calculating BIC...\n")
  bic <- mclustBIC(beta, G = 1:round(nrow(beta) / 2))
  cat("Running EM algorithm...\n")
  gmm <- Mclust(beta, x = bic)
  gmm
}

compress_methylation_multi <- function(infiles, outfile, clusters,
                                       distance_threshold, lambda,
                                       col_names = NULL, region = NULL,
                                       hdf5 = NULL) {
  num_of_clusters <- length(unique(clusters))
  if (length(lambda) == 1) {
    lambda <- rep(lambda, num_of_clusters)
  } else if (length(lambda) != num_of_clusters) {
    stop("Length of lambda vector should be equal to the number of unique
     clusters!")
  }

  methdata <- read_meth_file(infiles, hdf5 = hdf5, col_names = col_names)

  if (!is.null(region)) {
    region_ii <-
      queryHits(findOverlaps(methdata@rowRanges,
                             makeGRangesFromDataFrame(region),
                             type = "within"))
    methdata <- methdata[region_ii]
  }

  gap_distances <- diff(start(methdata@rowRanges))
  gaps <- which(distance_threshold < gap_distances)

  chrom_splits <- cumsum(seqnames(methdata@rowRanges)@lengths)
  row_partitions <- combine_two_bp_sets(gaps, chrom_splits)

  ii <- order(clusters)
  clusters <- clusters[ii]
  col_partitions <-
    c(which(clusters[2:length(clusters)] != clusters[1:(length(clusters) - 1)]),
      length(clusters))

  breakpoints <- segment_block_reduce(methdata,
                                      row_partitions,
                                      col_partitions,
                                      lambda)
  breakpoints <- combine_two_bp_sets(row_partitions, breakpoints)
  segments <- format_segments(methdata@rowRanges, chrom_splits, breakpoints)

  if (!is.null(col_names)){
    infiles <- col_names
  }

  meth_out <- aggregate_meth_data(methdata, segments, infiles, "M")
  cov_out <- aggregate_meth_data(methdata, segments, infiles, "Cov")

  write_segmented_bsseq(cov_out, meth_out, paste0(outfile, ".bsseq"))
}
