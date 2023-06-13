read_meth_file <- function(file, workers = 8, hdf5 = NULL){
  BPPARAM = BiocParallel::MulticoreParam(progressbar = TRUE, workers = workers)
  
  if(is.null(hdf5)){
    hdf5 = tempfile()
  }else if(file.exists(hdf5)){
    return(HDF5Array::loadHDF5SummarizedExperiment(dir=hdf5))
  }

  methfile = read.bismark(files = file,
                          colData = data.frame(row.names=format(1:length(file))),
                          rmZeroCov = T,
                          verbose = T,
                          BPPARAM = BPPARAM,
                          BACKEND = "HDF5Array",
                          dir = hdf5,
                          nThread=1)
  methfile
}

write_segmented_methylation <- function(segments, file, col.names=F){
  write.table(segments, file, append = T, sep = "\t", quote = F, row.names = F, col.names = col.names)
}

binomialsegmentation_breakpoints <- function(meth, cov, offset, lambda){
  segments = fusedsegmentation(meth, C=cov, lambda2=lambda, objective="binomial", l=0)#, format="full")
  segments$end[1:(nrow(segments)-1)] + offset
}

combine_two_bp_sets <- function(x, y){
  combine_two_bp_sets_(as.integer(x), as.integer(y))
}

segment_block_reduce <- function(data, rowp, colp, lambda){
  grid = ArbitraryArrayGrid(tickmarks = list(rowp, colp))
  
  summary_list = bplapply(seq_along(grid), 
                          function(b) {
                            viewport <- grid[[b]]
                            methylation <- rowSums(read_block(data@assays@data$M, viewport))
                            coverage <- rowSums(read_block(data@assays@data$Cov, viewport))
                            if(length(methylation)==1) return(realize(matrix(end(ranges(viewport))[1],ncol=1), BACKEND="HDF5Array"))
                            coverage[coverage==0] = 1
                            bps = binomialsegmentation_breakpoints(methylation, coverage, start(ranges(viewport))[1]-1, lambda)
                            realize(matrix(bps,ncol=1), BACKEND="HDF5Array")
                          }
  )
  Reduce(combine_two_bp_sets, summary_list)
} 

format_segments <- function(ranges, chrom_end_points, breakpoints){
  segments = data.frame(chrom = seqnames(ranges)[breakpoints], 
                        start = c(start(ranges)[1], 1 + start(ranges)[breakpoints[1:(length(breakpoints)-1)]]), 
                        end = start(ranges)[breakpoints])

  segments_chrom_splits = cumsum(Rle(as.factor(segments$chrom))@lengths)[1:(length(chrom_end_points) - 1)] +1
  
  if(length(chrom_end_points) > 1){
      segments[segments_chrom_splits, 2] = start(ranges)[chrom_end_points[1:(length(chrom_end_points) - 1)]+1] 
  }
  
  segments
}

aggregate_meth_data <- function(data, segments, meth_colnames, type){
  out = getCoverage(data, regions = makeGRangesFromDataFrame(segments), type=type, what="perRegionTotal")
  colnames(out) = meth_colnames
  cbind(segments, out)
}

compressMethylationMulti <- function(infiles, outfile, clusters, distance_threshold, lambda, col_names = NULL, region = NULL, hdf5 = NULL){
  meth_file = paste0(outfile, ".M.bedGraph")
  cov_file = paste0(outfile, ".cov.bedGraph")
  if(file.exists(meth_file)){
    unlink(meth_file)
  }
  if(file.exists(cov_file)){
    unlink(cov_file)
  }
  
  methdata = read_meth_file(infiles, hdf5 = hdf5)

  if(!is.null(region)){
    region_ii = queryHits(findOverlaps(methdata@rowRanges, makeGRangesFromDataFrame(region), type="within"))
    methdata = methdata[region_ii]
  }

  gap_distances = diff(start(methdata@rowRanges))
  gaps = which(distance_threshold < gap_distances)
  
  chrom_splits = cumsum(seqnames(methdata@rowRanges)@lengths)
  row_partitions = combine_two_bp_sets(gaps, chrom_splits)
  
  ii = order(clusters)
  clusters = clusters[ii]
  col_partitions = c(which(clusters[2:length(clusters)] != clusters[1:(length(clusters) - 1)]), length(clusters))

  breakpoints = segment_block_reduce(methdata, row_partitions, col_partitions, lambda)
  breakpoints = combine_two_bp_sets(row_partitions, breakpoints)
  segments = format_segments(methdata@rowRanges, chrom_splits, breakpoints)

  if(!is.null(col_names)){
      infiles = col_names
  }

  meth_out = aggregate_meth_data(methdata, segments, infiles, "M") 
  cov_out = aggregate_meth_data(methdata, segments, infiles, "Cov") 

  write_segmented_methylation(meth_out, meth_file, col.names=T)
  write_segmented_methylation(cov_out, cov_file, col.names=T)
}
