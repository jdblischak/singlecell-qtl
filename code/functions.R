# Utility function to force R to output a specific number of decimal places.
# Avoids Git thinking results have changed simply because of slight changes
# in insignificant digits.
# http://stackoverflow.com/a/12135122
specify_decimal <- function(x, k) format(round(x, k), nsmall = k)

# Show the upper left corner of a table
headl <- function(x, nr = 6, nc = 6) {
  stopifnot(class(x)[1] %in% c("data.frame", "data.table", "matrix", "tbl_df"),
            nr > 0, nc > 0)
  x[seq(nr), seq(nc)]
}

# Show the upper right corner of a table
headr <- function(x, nr = 6, nc = 6) {
  stopifnot(class(x)[1] %in% c("data.frame", "data.table", "matrix", "tbl_df"),
            nr > 0, nc > 0)
  xc <- ncol(x)
  x[seq(nr), seq(xc - nc + 1, xc)]
}

# Show the lower left corner of a table
taill <- function(x, nr = 6, nc = 6) {
  stopifnot(class(x)[1] %in% c("data.frame", "data.table", "matrix", "tbl_df"),
            nr > 0, nc > 0)
  xr <- nrow(x)
  x[seq(xr - nr + 1, xr), seq(nc)]
}

# Show the lower right corner of a table
tailr <- function(x, nr = 6, nc = 6) {
  stopifnot(class(x)[1] %in% c("data.frame", "data.table", "matrix", "tbl_df"),
            nr > 0, nc > 0)
  xr <- nrow(x)
  xc <- ncol(x)
  x[seq(xr - nr + 1, xr), seq(xc - nc + 1, xc)]
}

# Perform Principal Components Analysis (PCA).
#
# Args:
#   x: gene-by-sample matrix
#   retx, center, scale - see ?prcomp
#
# Returns a list with the following elements:
#   PCs - sample-by-PC matrix of principal components
#   explained - proportion of variance explained by each PC
#
# Reference: Zhang et al. 2009 (http://www.ncbi.nlm.nih.gov/pubmed/19763933)
run_pca <- function(x, retx = TRUE, center = TRUE, scale = TRUE) {
  library("testit")

  pca <- prcomp(t(x), retx = retx, center = center, scale. = scale)
  variances <- pca$sdev^2
  explained <- variances / sum(variances)
  assert("Variance explained is calculated correctly.",
         explained[1:2] - summary(pca)$importance[2, 1:2] < 0.0001)
  return(list(PCs = pca$x, explained = explained))
}

# Plot PCA results.
#
# Args:
#   x: numeric matrix of PCs
#   pcx: PC to plot on x-axis (default: 1)
#   pcy: PC to plot on y-axis (default: 2)
#   explained: numeric vector of fractions of variance explained by each PC
#   metadata: data frame or matrix that contains the metadata used to annotate
#             the plot
#   color, shape, size: column name of metadata used to pass column to ggplot
#                       aesthetic
#   factors: character vector which contains the column names of metadata that
#            need to be explicitly converted to a factor
#   ... : Additional arguments passed to geom_point
plot_pca <- function(x, pcx = 1, pcy = 2, explained = NULL, metadata = NULL,
                     color = NULL, shape = NULL, factors = NULL,
                     ...) {
  library("ggplot2")
  library("testit")

  # Prepare data frame to pass to ggplot
  if (!is.null(metadata)) {
    assert("PC and metadata have same number of rows.",
           nrow(x) == nrow(metadata))
    plot_data <- cbind(x, metadata)
    plot_data <- as.data.frame(plot_data)
    # Convert numeric factors to class "factor"
    for (f in factors) {
      plot_data[, f] <- as.factor(plot_data[, f])
    }
  } else {
    plot_data <- as.data.frame(x)
  }
  # Prepare axis labels
  if (!is.null(explained)) {
    assert("Number of PCs differs between x and explained.",
           length(explained) == ncol(x))
    xaxis <- sprintf("PC%d (%.2f%%)", pcx, round(explained[pcx] * 100, 2))
    yaxis <- sprintf("PC%d (%.2f%%)", pcy, round(explained[pcy] * 100, 2))
  } else {
    xaxis <- paste0("PC", pcx)
    yaxis <- paste0("PC", pcy)
  }
  # Plot
  p <- ggplot(plot_data, aes_string(x = paste0("PC", pcx),
                                    y = paste0("PC", pcy),
                                    color = color,
                                    shape = shape)) +
    geom_point(...) +
    labs(x = xaxis, y = yaxis)
  p
}
