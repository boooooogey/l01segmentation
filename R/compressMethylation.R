read_meth_file <- function(file, workers = 8){
  BPPARAM = MulticoreParam(progressbar = TRUE, workers = workers)
  tmp = tempfile()
  
  methfile = read.bismark(files = file,
                          colData = data.frame(row.names="1"),
                          rmZeroCov = T,
                          verbose = T,
                          BPPARAM = BPPARAM,
                          BACKEND = "HDF5Array",
                          dir = tmp,
                          nThread=1)
}

write_segmented_methylation <- function(segments, file){
  write.table(segments, file, append = T, sep = "\t", quote = F, row.names = F, col.names = F)
}

get_cov_meth <- function(file, chrom, start, end){
  cov = getCoverage(file, 
                    regions = data.frame(chr=chrom, start=start, end=end), 
                    type = "Cov")[[1]]
  meth = getCoverage(file, 
                     regions = data.frame(chr=chrom, start=start, end=end), 
                     type = "M")[[1]]
  list("cov"=cov, "meth" = meth)
}

binomialsegmentation <- function(meth, cov, ranges, chrom, lambda){
  segments = fusedsegmentation(meth, C=cov, lambda2=lambda, objective="binomial", l=0)#, format="full")
  segments$start = ranges[segments$start]
  segments$end = ranges[segments$end]
  segments = add_column(segments, .before="start", chrom=chrom)
  segments
}

compressMethylation <- function(infile, outfile, distance_threshold, lambda){
  if(file.exists(outfile)){
    unlink(outfile)
  }
  
  methfile = read_meth_file(infile)
  
  rowranges = rowRanges(methfile)
  
  chrom_list = levels(seqnames(rowranges))
  
  for(chrom in chrom_list){
    rowranges_chrom = rowranges[seqnames(rowranges) == chrom]
    
    gap_distances = diff(start(rowranges_chrom))
    gaps = which(distance_threshold < gap_distances)
  
    for(i in 1:length(gaps)){
      if (i == 1){
        start = 1
      }else{
        start = end+1
      }
      end = gaps[i]
      tmp_rowranges = start(rowranges_chrom[start:end])
      n = length(tmp_rowranges)
      data = get_cov_meth(methfile, chrom, tmp_rowranges[1], tmp_rowranges[length(tmp_rowranges)])
      segments = binomialsegmentation(data$meth[1:n], data$cov[1:n], tmp_rowranges, chrom, lambda)
      write_segmented_methylation(segments, outfile)
    }
  }
}
