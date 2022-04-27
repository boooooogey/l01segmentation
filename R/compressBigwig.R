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
  
compressBigwig <- function(in_bigwig_path, out_bedgraph_path, lambda, partition_length = 1e6, n_cores = 4){
  
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
  
  lambda_vector = rep(lambda, partition_length-1)
  weight_vector = rep(1, partition_length)
  for(partition_sequence in sequence_list){ 
    start = 1
    end = sequence_lengths[partition_sequence]
    number_of_iterations = ceiling((end - start + 1) / partition_length)

    foreach(i = 1:number_of_iterations, .init=conn, .combine='wtab') %dopar% {
      partition_start = (i-1) * partition_length + 1
      partition_end = min(i * partition_length, end)
      partition_size = partition_end - partition_start + 1
      data = readFromBigwig(in_bigwig_path, partition_sequence, partition_start, partition_end)
      bedg = data.frame(L0PoissonApproximateCondensed(data, 
                                                      lambda_vector[1:(partition_size-1)], 
                                                      weight_vector[1:partition_size]))
      bedg$seq = partition_sequence
      bedg$start = bedg$start + partition_start - 1
      bedg$end = bedg$end + partition_start - 1
      bedg[dim(bedg)[1]:1, c(4, 2, 3, 1)] 
    }

  }
  
  on.exit(close(conn))
}
