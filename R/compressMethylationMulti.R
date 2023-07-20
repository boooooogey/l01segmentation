read_meth_file <- function(file, workers = 8, hdf5 = NULL, colnames = NULL) {
  bpparam <- BiocParallel::MulticoreParam(progressbar = TRUE, workers = workers)

  if (is.null(hdf5)) {
    hdf5 <- tempfile()
  } else if (file.exists(hdf5)) {
    return(HDF5Array::loadHDF5SummarizedExperiment(dir = hdf5))
  }

  if (is.null(colnames)) {
    colnames <- data.frame(row.names = seq_along(file))
  } else {
    colnames <- data.frame(row.names = colnames)
  }

  methfile <- read.bismark(files = file,
                           colData = colnames,
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

write_segmented_bsseq <- function(cov, meth, dir) {
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

binomialsegmentation_breakpoints <- function(meth, cov, offset, lambda){
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

  summary_list <- bplapply(seq_along(grid),
    function(b) {
      viewport <- grid[[b]]
      methylation <- rowSums(read_block(data@assays@data$M, viewport))
      coverage <- rowSums(read_block(data@assays@data$Cov, viewport))
      if (length(methylation) == 1)
        return(realize(matrix(end(ranges(viewport))[1], ncol = 1),
                       BACKEND = "HDF5Array"))
      coverage[coverage == 0] <- 1
      bps <- binomialsegmentation_breakpoints(methylation,
                                              coverage,
                                              start(ranges(viewport))[1] - 1,
                                              lambda)
      realize(matrix(bps, ncol = 1), BACKEND = "HDF5Array")
    }
  )
  Reduce(combine_two_bp_sets, summary_list)
}

format_segments <- function(ranges, chrom_end_points, breakpoints) {
  segments <- data.frame(chrom = seqnames(ranges)[breakpoints],
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

compressMethylationMulti <- function(infiles, outfile, clusters,
                                     distance_threshold, lambda,
                                     col_names = NULL, region = NULL,
                                     hdf5 = NULL) {
  meth_file <- paste0(outfile, ".M.bedGraph")
  cov_file <- paste0(outfile, ".cov.bedGraph")
  if (file.exists(meth_file)){
    unlink(meth_file)
  }
  if (file.exists(cov_file)){
    unlink(cov_file)
  }

  methdata <- read_meth_file(infiles, hdf5 = hdf5, colnames = col_names)

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

  #write_segmented_methylation(meth_out, meth_file, col.names = TRUE)
  #write_segmented_methylation(cov_out, cov_file, col.names = TRUE)
  write_segmented_bsseq(cov_out, meth_out, paste0(outfile, ".bsseq"))
}
