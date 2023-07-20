# DMLtest from DSS package modified to work with segmented data.
# github.com/haowulab/DSS/blob/7b6a3c4563827406b87dc60c3475bd919c9b48b2/R/DML.R#L11
DMLtest <- function(bsseq, labels, ncores) {

  x <- as.array(getCoverage(bsseq, type = "M"))
  n <- as.array(getCoverage(bsseq, type = "Cov"))

  x1 <- x[, labels == 1]
  n1 <- n[, labels == 1]
  x2 <- x[, labels == 2]
  n2 <- n[, labels == 2]

  allchr <- as.character(seqnames(bsseq))
  allstart <- start(bsseq)
  estprob1 <- compute.mean.noSmooth(x1, n1)
  estprob2 <- compute.mean.noSmooth(x2, n2)

  cat("Estimating dispersion for each CpG site, this will take a while ...\n")
  phi1 <- est.dispersion.BSseq(x1, n1, estprob1, ncores)
  phi2 <- est.dispersion.BSseq(x2, n2, estprob2, ncores)

  wt1 <- 1 / (1 + (n1 - 1) * phi1)
  wt1 <- wt1 / mean(wt1)
  wt2 <- 1 / (1 + (n2 - 1) * phi2)
  wt2 <- wt2 / mean(wt2)
  x1_wt <- x1 * wt1
  n1_wt <- n1 * wt1
  x2_wt <- x2 * wt2
  n2_wt <- n2 * wt2
  estprob1 <- compute.mean.noSmooth(x1_wt, n1_wt)
  estprob2 <- compute.mean.noSmooth(x2_wt, n2_wt)
  wald <- waldTest.DML(x1_wt, n1_wt, estprob1, phi1, x2_wt,
                       n2_wt, estprob2, phi2, smoothing = FALSE,
                       allchr = allchr,
                       allpos = allstart)
  segment_positions <- data.frame(bsseq@rowRanges)
  colnames(segment_positions)[1] <- "chr"
  wald$start <- wald$pos
  wald <- merge(wald, segment_positions, by = c("chr", "start"))
  wald <- wald[order(wald$chr, wald$start), ]
  return(wald[, c("chr", "start", "end", "mu1", "mu2", "diff", "diff.se",
                  "start", "phi1", "phi2", "pval", "fdr")])
}