readFromBigwig <- function(path, sequence, start, end){
  which = GRanges(c(sequence), IRanges(c(start), c(end)))
  import(BigWigFile(path), which = which, as = "NumericList")[[1]]
}

readSequenceLenghtsFromBigwig <- function(path){
  seqlengths(BigWigFile(path))
}

wtab <- function(conn, d) {
  write.table(format(d, scientific = F), conn, col.names=F, row.names = F, sep = " ", quote = F, append = T)
  conn
}
  
#' @title compressBigwig
#' @description
#' Compresses the signal stored in the given BigWig file solving a L0 segmentation problem and stores the output into the given BedGraph file. 
#' @param in_bigwig_path Path to the BigWig file to be compressed.
#' @param out_bedgraph_path Path to the BedGraph file. The result of the partition will be saved in this file. 
#' @param partition_length The number of loci to be read into the memory from the BigWig file. The default is 1e6.
#' @param n_cores The number of cores to be used in the segmentation. The default is 4.
#' @param bin Size of the data bins. If it has a value other than NULL, the data is binned before segmentation.
#' @param ... To pass parameters to the L0 segmentation backend. For example, lambda2.
#'
compressBigwig <- function(in_bigwig_path, out_bedgraph_path, partition_length = 1e6, n_cores = 4, bin = NULL, ...){
  
  if(file.exists(out_bedgraph_path)){
    file.remove(out_bedgraph_path)
  }
  
  number_of_cores = detectCores()
  
  if (n_cores > number_of_cores){
    registerDoParallel(number_of_cores)
  } else{
    registerDoParallel(n_cores)
  }
  
  conn = file(out_bedgraph_path, "w")
  
  sequence_lengths = readSequenceLenghtsFromBigwig(in_bigwig_path)
  sequence_list = sort(names(sequence_lengths))
  
  for(partition_sequence in sequence_list){ 
    start = 1
    end = sequence_lengths[partition_sequence]
    number_of_iterations = ceiling((end - start + 1) / partition_length)

    foreach(i = 1:number_of_iterations, .init=conn, .combine='wtab') %dopar% {
      partition_start = (i-1) * partition_length + 1
      partition_end = min(i * partition_length, end)
      partition_size = partition_end - partition_start + 1
      data = readFromBigwig(in_bigwig_path, partition_sequence, partition_start, partition_end)

      if(!is.null(bin)){
        data = binVector(data, bin)
        bin_starts = data$start
        bin_ends = data$end
        data = data$value
      }

      bedg = data.frame(fusedsegmentation(data,...))

      if(!is.null(bin)){
        bedg$start = bin_starts[bedg$start] + partition_start - 1
        bedg$end = bin_ends[bedg$end] + partition_start - 1
      }else{
        bedg$start = bedg$start + partition_start - 2
        bedg$end = bedg$end + partition_start - 1
      }

      bedg$seq = partition_sequence
      bedg[, c("seq", "start", "end", "value")] #c(4, 2, 3, 1)] 
    }
  }

  on.exit(close(conn))
}

#' @title compressBigwigbyBinning
#' @description
#' Bins the signal stored in the given BigWig file and stores the output into the given BedGraph file. 
#' @param in_bigwig_path Path to the BigWig file to be compressed.
#' @param out_bedgraph_path Path to the BedGraph file. The binned data will be saved in this file. 
#' @param partition_length The number of loci to be read into the memory from the BigWig file. The default is 1e6.
#' @param n_cores The number of cores to be used in the segmentation. The default is 4.
#' @param bin Size of the data bins. 
#'
compressBigwigbyBinning <- function(in_bigwig_path, out_bedgraph_path, bin_size, partition_length = 1e6, n_cores = 4){
  
  if(file.exists(out_bedgraph_path)){
    file.remove(out_bedgraph_path)
  }
  
  number_of_cores = detectCores()
  
  if (n_cores > number_of_cores){
    registerDoParallel(number_of_cores)
  } else{
    registerDoParallel(n_cores)
  }
  
  conn = file(out_bedgraph_path, "w")
  
  sequence_lengths = readSequenceLenghtsFromBigwig(in_bigwig_path)
  sequence_list = sort(names(sequence_lengths))
  
  for(partition_sequence in sequence_list){ 
    start = 1
    end = sequence_lengths[partition_sequence]
    number_of_iterations = ceiling((end - start + 1) / partition_length)

    foreach(i = 1:number_of_iterations, .init=conn, .combine='wtab') %dopar% {
      partition_start = (i-1) * partition_length + 1
      partition_end = min(i * partition_length, end)
      partition_size = partition_end - partition_start + 1
      data = readFromBigwig(in_bigwig_path, partition_sequence, partition_start, partition_end)
      bedg = data.frame(binVector(data, 
                                  bin_size))
      bedg$seq = partition_sequence
      bedg$start = bedg$start + partition_start - 1
      bedg$end = bedg$end + partition_start - 1
      bedg[, c("seq", "start", "end", "value")] #c(4, 2, 3, 1)] 
    }

  }
  
  on.exit(close(conn))
}
