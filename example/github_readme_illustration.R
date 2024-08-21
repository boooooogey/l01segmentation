library(ggplot2)

load("github_readme_example.rda")

chrom <- "chr6"
start <- 32934371
end <- 32951278

plot_track <- function(index, track, line = NULL) {
  df <- data.frame(index = index, track = track)
  p <- ggplot(df, aes(x = index, y = track, color = "Track")) + geom_point()
  if (!is.null(line)) {
    p <- p + geom_line(aes(y = line, color = "L0 approx"), linewidth = 2)
    p <- p + scale_color_manual(
      values = c("Track" = "black", "L0 approx" = "red")
    )
    p <- p + guides(colour=guide_legend(title = ""))
    p <- p + theme(legend.position = "bottom")
  }
  p + ylab("Count") + xlab("Index") + theme(text = element_text(size = 20))
}

p3 <- plot_track(
  start:end,
  example,
  l01segmentation::fusedsegmentation(
    example,
    lambda2 = 125,
    objective = "poisson",
    format = "full"
  )
)
