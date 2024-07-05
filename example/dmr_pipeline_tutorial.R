library(l01segmentation)

# If the positions are apart by more than this threshold than put them in diffe-
# rent segments.
distance_threshold <- 10000

hdf5_file <- "P0_FB_vs_LV_segmentation_chrom2"

outfile <- "P0_FB_vs_LV_chrom2"

input_folder <- "./dmr_pipeline_data"
files <- file.path(input_folder, c("mc_P0_FB_1_2.txt.gz",
                                   "mc_P0_FB_2_2.txt.gz",
                                   "mc_P0_LV_1_2.txt.gz",
                                   "mc_P0_LV_2_2.txt.gz"))

names <- c("FB_1", "FB_2", "LV_1", "LV_2") 

chrom <- "2"

# Lambda for L0 segmentation, different values for different samples can be spe-
# cified.
lambda <- 75

#Cluster methylation samples based on binned data from a region. Uses GMM.
gmm <- cluster_methylation(infiles = files,
                           hdf5 = hdf5_file,
                           chrom = chrom,
                           start = 1000000,
                           end = 10000000,
                           binwidth = 10000,
                           col_names = names)

# Combined L0 segmentation of all labels from GMM clustering. Segmentation is
# stored as Bismark objects in hdf5 files.
compress_methylation_multi(infiles = files, 
                           outfile = paste0(outfile, "_", lambda), 
                           clusters = gmm$classification, 
                           distance_threshold = distance_threshold, 
                           lambda = lambda, 
                           col_names = names,
                           hdf5 = hdf5_file)
                         
# Load the L0 segmentation from the hdf5 file
bsseq_segmented <- HDF5Array::loadHDF5SummarizedExperiment(dir = paste0(outfile,
                                                                        "_",
                                                                        lambda,
                                                                        ".bsseq"
                                                                        ))

# Define groups for DML test
groups <- grepl("FB", colnames(bsseq_segmented)) * 1 +
           grepl("LV", colnames(bsseq_segmented)) * 2

number_of_cores <- 1L

# Carry out DML test
pvals <- l01segmentation::dml_test(bsseq_segmented, as.integer(groups), ncores = number_of_cores)

# Write p values in a file
write.table(x = pvals,
            file = paste0("LV_vs_FB_dmrs_chr", chrom, "_", lambda, ".tsv"),
            sep="\t", quote = F, row.names = F, col.names = T)